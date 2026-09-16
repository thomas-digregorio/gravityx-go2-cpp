#include "gravityx/sparse_ac_economic.hpp"

#include "gravityx/fast_power_flow.hpp"
#include "gravityx/state_io.hpp"

#include <coin/IpIpoptApplication.hpp>
#include <coin/IpSolveStatistics.hpp>
#include <coin/IpTNLP.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace gravityx {
namespace {

constexpr double kInfinity = 1e19;

struct FlowTerm {
    double vf2{};
    double vt2{};
    double cross_cos{};
    double cross_sin{};
};

struct BranchCoefficients {
    std::array<FlowTerm, 4> term;
};

struct BranchEvaluation {
    std::array<double, 4> flow{};
    // Variable order: vm_from, vm_to, va_from, va_to.
    std::array<std::array<double, 4>, 4> derivative{};
};

using Complex = std::complex<double>;
using YRows = std::vector<std::map<int, Complex>>;

void add_admittance(YRows& rows, int row, int column, Complex value) {
    rows[row][column] += value;
}

std::pair<double, double> base_pg_bounds(
    const Generator& generator,
    int commitment,
    double delta_r) {
    if (commitment == 0) {
        return {0.0, 0.0};
    }
    const double previous = generator.status_prev == 0
        ? generator.pmin : generator.pg_prev;
    const double lower = std::max(
        generator.pmin, previous - delta_r * generator.prdmax);
    const double upper = std::min(
        generator.pmax, previous + delta_r * generator.prumax);
    if (lower > upper + 1e-12) {
        throw std::runtime_error(
            "empty sparse AC generator interval: " + generator.source_key);
    }
    return {lower, upper};
}

std::pair<double, double> base_load_bounds(
    const Load& load,
    double delta_r) {
    if (std::abs(load.pd_nominal) <= 1e-12) {
        return {load.tmin, load.tmax};
    }
    const double lower = std::max(
        load.tmin,
        (load.pd_prev - delta_r * load.prdmax) / load.pd_nominal);
    const double upper = std::min(
        load.tmax,
        (load.pd_prev + delta_r * load.prumax) / load.pd_nominal);
    if (lower > upper + 1e-12) {
        throw std::runtime_error(
            "empty sparse AC load interval: " + load.source_key);
    }
    return {lower, upper};
}

std::pair<double, double> corrective_pg_bounds(
    const Generator& generator, int commitment, double original_pg,
    double delta_r_ctg) {
    if (commitment == 0) {
        return {0.0, 0.0};
    }
    const double lower = std::max(
        generator.pmin, original_pg - delta_r_ctg * generator.prdmaxctg);
    const double upper = std::min(
        generator.pmax, original_pg + delta_r_ctg * generator.prumaxctg);
    if (lower > upper + 1e-12) {
        throw std::runtime_error("empty corrective generator interval: " +
                                 generator.source_key);
    }
    return {lower, upper};
}

std::pair<double, double> corrective_load_bounds(
    const Load& load, double original_factor, double delta_r_ctg) {
    if (std::abs(load.pd_nominal) <= 1e-12) {
        return {load.tmin, load.tmax};
    }
    const double previous = load.pd_nominal * original_factor;
    // Deliberately matches the independent source-semantics checker.
    const double lower = std::max(
        load.tmin, (previous - delta_r_ctg * load.prdmaxctg) / load.pd_nominal);
    const double upper = std::min(
        load.tmax, (previous + delta_r_ctg * load.prumaxctg) / load.pd_nominal);
    if (lower > upper + 1e-12) {
        throw std::runtime_error("empty corrective load interval: " + load.source_key);
    }
    return {lower, upper};
}

double rebuild_economic_candidate(
    const CaseData& data, const std::vector<int>& commitment,
    const AcState* original_base, AcState& state) {
    return original_base != nullptr
        ? rebuild_common_corrective_reference_state(data, *original_base, commitment, state)
        : rebuild_base_state_derived_fields(data, commitment, state);
}

ValidationReport validate_economic_candidate(
    const CaseData& data, const std::vector<int>& commitment,
    const AcState* original_base, const AcState& state) {
    if (original_base != nullptr) {
        ContingencyContext context;
        context.borrow_base_state(*original_base);
        return validate_state(data, ModelMode::ContingencySoft, state, commitment, context);
    }
    return validate_state(data, ModelMode::BaseSoft, state, commitment);
}

std::pair<double, double> pwl_value_slope(
    const std::vector<PwlPoint>& points,
    double power) {
    if (points.size() < 2) {
        throw std::runtime_error("sparse AC PWL curve has fewer than two points");
    }
    const auto slope = [&] (std::size_t segment) {
        const double width =
            points[segment + 1].mw - points[segment].mw;
        return std::abs(width) <= 1e-14
            ? 0.0
            : (points[segment + 1].cost - points[segment].cost) / width;
    };
    std::size_t segment = points.size() - 2;
    for (std::size_t i = 0; i + 1 < points.size(); ++i) {
        if (power <= points[i + 1].mw + 1e-12) {
            segment = i;
            break;
        }
    }
    const auto& left = points[segment];
    const double value = left.cost +
        slope(segment) * (power - left.mw);
    double derivative = slope(segment);
    if (segment + 1 < points.size() - 1 &&
        std::abs(power - points[segment + 1].mw) <= 1e-10) {
        derivative = 0.5 * (derivative + slope(segment + 1));
    }
    return {value, derivative};
}

struct EpigraphSegment {
    double slope{};
    double intercept{};
};

bool convex_epigraph_segments(
    const std::vector<PwlPoint>& points, double sign,
    std::vector<EpigraphSegment>& segments, double& scale) {
    segments.clear();
    scale = 1.0;
    if (points.size() < 2) {
        return false;
    }
    double previous_slope = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i + 1 < points.size(); ++i) {
        const auto& left = points[i];
        const auto& right = points[i + 1];
        const double width = right.mw - left.mw;
        const double slope = sign * (right.cost - left.cost) / width;
        const double intercept = sign * left.cost - slope * left.mw;
        // Never convexify a source curve or repair its slopes. Even a small
        // decreasing slope retains the original piecewise representation.
        if (!(width > 1e-14) || !std::isfinite(slope) ||
            !std::isfinite(intercept) || slope < previous_slope) {
            segments.clear();
            return false;
        }
        segments.push_back({slope, intercept});
        scale = std::max(scale, std::abs(slope));
        previous_slope = slope;
    }
    return true;
}

BranchCoefficients branch_coefficients(const Branch& branch) {
    const double denominator = branch.r * branch.r + branch.x * branch.x;
    const double g = denominator > 1e-20 ? branch.r / denominator : 0.0;
    const double b = denominator > 1e-20 ? -branch.x / denominator : 0.0;
    if (std::abs(branch.tap) <= 1e-12) {
        throw std::runtime_error(
            "zero tap ratio in sparse AC NLP: " + branch.source_key);
    }
    const double tm2 = branch.tap * branch.tap;
    const double tr = branch.tap * std::cos(branch.shift);
    const double ti = branch.tap * std::sin(branch.shift);
    const double from_g_self = branch.transformer
        ? g / tm2 + branch.g_fr
        : (g + branch.g_fr) / tm2;
    const double from_b_self = branch.transformer
        ? b / tm2 + branch.b_fr
        : (b + branch.b_fr) / tm2;
    BranchCoefficients result;
    result.term[0] = {
        from_g_self, 0.0,
        (-g * tr + b * ti) / tm2,
        (-b * tr - g * ti) / tm2};
    result.term[1] = {
        -from_b_self, 0.0,
        (b * tr + g * ti) / tm2,
        (-g * tr + b * ti) / tm2};
    result.term[2] = {
        0.0, g + branch.g_to,
        (-g * tr - b * ti) / tm2,
        (b * tr - g * ti) / tm2};
    result.term[3] = {
        0.0, -(b + branch.b_to),
        (b * tr - g * ti) / tm2,
        (g * tr + b * ti) / tm2};
    return result;
}

BranchEvaluation evaluate_branch(
    const BranchCoefficients& coefficients,
    double vm_from,
    double vm_to,
    double angle_difference) {
    const double cosine = std::cos(angle_difference);
    const double sine = std::sin(angle_difference);
    BranchEvaluation result;
    for (int component = 0; component < 4; ++component) {
        const auto& term = coefficients.term[component];
        const double cross =
            term.cross_cos * cosine + term.cross_sin * sine;
        const double angle =
            -term.cross_cos * sine + term.cross_sin * cosine;
        result.flow[component] =
            term.vf2 * vm_from * vm_from +
            term.vt2 * vm_to * vm_to +
            cross * vm_from * vm_to;
        result.derivative[component][0] =
            2.0 * term.vf2 * vm_from + cross * vm_to;
        result.derivative[component][1] =
            2.0 * term.vt2 * vm_to + cross * vm_from;
        result.derivative[component][2] =
            angle * vm_from * vm_to;
        result.derivative[component][3] =
            -result.derivative[component][2];
    }
    return result;
}

using FlowHessians = std::array<std::array<std::array<double, 4>, 4>, 4>;

FlowHessians branch_flow_hessians(
    const BranchCoefficients& coefficients, double vf, double vt, double angle) {
    FlowHessians result{};
    const double cosine = std::cos(angle), sine = std::sin(angle);
    for (int k = 0; k < 4; ++k) {
        const auto& term = coefficients.term[k];
        const double cross = term.cross_cos * cosine + term.cross_sin * sine;
        const double first = -term.cross_cos * sine + term.cross_sin * cosine;
        auto& h = result[k];
        h[0][0] = 2.0 * term.vf2;
        h[1][1] = 2.0 * term.vt2;
        h[1][0] = cross;
        h[2][0] = vt * first;
        h[3][0] = -vt * first;
        h[2][1] = vf * first;
        h[3][1] = -vf * first;
        h[2][2] = h[3][3] = -vf * vt * cross;
        h[3][2] = vf * vt * cross;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < i; ++j) h[j][i] = h[i][j];
        }
    }
    return result;
}

YRows build_ybus(const CaseData& data, const AcState& state) {
    YRows rows(data.buses.size());
    for (int i = 0; i < static_cast<int>(data.shunts.size()); ++i) {
        const auto& shunt = data.shunts[i];
        add_admittance(
            rows, shunt.bus, shunt.bus,
            {shunt.gs, effective_shunt_susceptance(data, state, i)});
    }
    for (const auto& branch : data.branches) {
        if (branch.status == 0) {
            continue;
        }
        const double denominator = branch.r * branch.r + branch.x * branch.x;
        const double g = denominator > 1e-20 ? branch.r / denominator : 0.0;
        const double b = denominator > 1e-20 ? -branch.x / denominator : 0.0;
        const double tm = branch.tap;
        if (std::abs(tm) <= 1e-12) {
            throw std::runtime_error(
                "zero tap ratio in sparse AC Ybus: " + branch.source_key);
        }
        const double tm2 = tm * tm;
        const Complex series(g, b);
        const Complex from_shunt(branch.g_fr, branch.b_fr);
        const Complex to_shunt(branch.g_to, branch.b_to);
        const Complex rotation = std::polar(1.0, branch.shift);
        const Complex yff = branch.transformer
            ? series / tm2 + from_shunt
            : (series + from_shunt) / tm2;
        const Complex yft = -series * rotation / tm;
        const Complex ytf = -series * std::conj(rotation) / tm;
        const Complex ytt = series + to_shunt;
        add_admittance(rows, branch.from, branch.from, yff);
        add_admittance(rows, branch.from, branch.to, yft);
        add_admittance(rows, branch.to, branch.from, ytf);
        add_admittance(rows, branch.to, branch.to, ytt);
    }
    return rows;
}

class SparseAcEconomicNlp final : public Ipopt::TNLP {
public:
    SparseAcEconomicNlp(
        const CaseData& data,
        std::vector<int> commitment,
        const AcState& start,
        double verification_tolerance,
        bool pwl_epigraph = false,
        const AcState* original_corrective_base = nullptr,
        bool exact_hessian = false)
        : data_(data),
          commitment_(std::move(commitment)),
          start_state_(start),
          original_corrective_base_(original_corrective_base),
          interval_duration_(original_corrective_base != nullptr
              ? data.delta_ctg : data.delta),
          verification_tolerance_(verification_tolerance),
          exact_hessian_(exact_hessian),
          nb_(static_cast<int>(data.buses.size())),
          ng_(static_cast<int>(data.generators.size())),
          nd_(static_cast<int>(data.loads.size())),
          nl_(static_cast<int>(data.branches.size())),
          vm_offset_(0),
          va_offset_(vm_offset_ + nb_),
          pg_offset_(va_offset_ + nb_),
          qg_offset_(pg_offset_ + ng_),
          demand_offset_(qg_offset_ + ng_),
          p_delta_offset_(demand_offset_ + nd_),
          q_delta_offset_(p_delta_offset_ + nb_),
          sm_offset_(q_delta_offset_ + nb_),
          variable_count_(sm_offset_ + nl_),
          balance_row_offset_(0),
          thermal_row_offset_(4 * nb_) {
        if (commitment_.size() != data.generators.size()) {
            throw std::runtime_error("sparse AC commitment dimension mismatch");
        }
        if (start.vm.size() != data.buses.size() ||
            start.va.size() != data.buses.size() ||
            start.pg.size() != data.generators.size() ||
            start.qg.size() != data.generators.size() ||
            start.demand_factor.size() != data.loads.size()) {
            throw std::runtime_error("sparse AC start dimension mismatch");
        }
        objective_scale_ = std::max({
            1.0,
            data.p_delta_cost_approx,
            data.q_delta_cost_approx,
            data.sm_cost_approx,
        });
        x_lower_.assign(variable_count_, -kInfinity);
        x_upper_.assign(variable_count_, kInfinity);
        start_x_.assign(variable_count_, 0.0);
        generator_points_.resize(ng_);
        load_points_.resize(nd_);
        generator_epigraph_column_.assign(ng_, -1);
        load_epigraph_column_.assign(nd_, -1);
        build_variables();
        if (pwl_epigraph) {
            build_pwl_epigraph();
        } else {
            pwl_original_curve_count_ = nd_ + static_cast<int>(
                std::count(commitment_.begin(), commitment_.end(), 1));
        }
        if (exact_hessian_ && (!pwl_epigraph || pwl_original_curve_count_ != 0)) {
            throw std::runtime_error(
                "exact Hessian requires all source curves to pass the smooth PWL epigraph check");
        }
        ybus_ = build_ybus(data_, start_state_);
        coefficients_.reserve(nl_);
        for (const auto& branch : data_.branches) {
            coefficients_.push_back(branch_coefficients(branch));
        }
        for (int branch = 0; branch < nl_; ++branch) {
            if (data_.branches[branch].status != 0) {
                active_branches_.push_back(branch);
            }
        }
        angle_row_offset_ = thermal_row_offset_ +
            2 * static_cast<int>(active_branches_.size());
        for (int branch : active_branches_) {
            const auto& item = data_.branches[branch];
            const double source_delta = original_corrective_base_ != nullptr
                ? original_corrective_base_->va[item.from] -
                    original_corrective_base_->va[item.to]
                : data_.buses[item.from].va_start - data_.buses[item.to].va_start;
            if (source_delta >= item.angmin &&
                source_delta <= item.angmax) {
                angle_branches_.push_back(branch);
            }
        }
        reference_row_offset_ = angle_row_offset_ +
            static_cast<int>(angle_branches_.size());
        for (int bus = 0; bus < nb_; ++bus) {
            if (data_.buses[bus].type == 3) {
                reference_buses_.push_back(bus);
            }
        }
        epigraph_row_offset_ = reference_row_offset_ +
            static_cast<int>(reference_buses_.size());
        constraint_count_ = epigraph_row_offset_ +
            static_cast<int>(epigraph_rows_.size());
        build_jacobian_structure();
        if (exact_hessian_) build_hessian_structure();
        std::vector<double> initial_constraints(constraint_count_);
        evaluate_constraints(start_x_.data(), initial_constraints.data());
        initial_constraint_violation_ = constraint_violation(
            start_x_.data(), initial_constraints.data());
    }

    bool get_nlp_info(
        Ipopt::Index& n,
        Ipopt::Index& m,
        Ipopt::Index& nnz_jac_g,
        Ipopt::Index& nnz_h_lag,
        IndexStyleEnum& index_style) override {
        n = variable_count_;
        m = constraint_count_;
        nnz_jac_g = static_cast<Ipopt::Index>(jacobian_rows_.size());
        nnz_h_lag = static_cast<Ipopt::Index>(hessian_rows_.size());
        index_style = C_STYLE;
        return true;
    }

    bool get_bounds_info(
        Ipopt::Index n,
        Ipopt::Number* x_l,
        Ipopt::Number* x_u,
        Ipopt::Index m,
        Ipopt::Number* g_l,
        Ipopt::Number* g_u) override {
        if (n != variable_count_ || m != constraint_count_) {
            return false;
        }
        std::copy(x_lower_.begin(), x_lower_.end(), x_l);
        std::copy(x_upper_.begin(), x_upper_.end(), x_u);
        std::fill(g_l, g_l + m, -kInfinity);
        std::fill(g_u, g_u + m, 0.0);
        for (int i = 0; i < static_cast<int>(angle_branches_.size()); ++i) {
            const auto& branch = data_.branches[angle_branches_[i]];
            g_l[angle_row_offset_ + i] = branch.angmin;
            g_u[angle_row_offset_ + i] = branch.angmax;
        }
        for (int i = 0; i < static_cast<int>(reference_buses_.size()); ++i) {
            g_l[reference_row_offset_ + i] = 0.0;
            g_u[reference_row_offset_ + i] = 0.0;
        }
        return true;
    }

    bool get_starting_point(
        Ipopt::Index n,
        bool init_x,
        Ipopt::Number* x,
        bool init_z,
        Ipopt::Number*,
        Ipopt::Number*,
        Ipopt::Index,
        bool init_lambda,
        Ipopt::Number*) override {
        if (n != variable_count_ || !init_x || init_z || init_lambda) {
            return false;
        }
        std::copy(start_x_.begin(), start_x_.end(), x);
        return true;
    }

    bool eval_f(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Number& objective) override {
        if (n != variable_count_) {
            return false;
        }
        objective = objective_value(x);
        return std::isfinite(objective);
    }

    bool eval_grad_f(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Number* gradient) override {
        if (n != variable_count_) {
            return false;
        }
        std::fill(gradient, gradient + n, 0.0);
        for (int i = 0; i < ng_; ++i) {
            if (commitment_[i] == 0 || generator_epigraph_column_[i] >= 0) {
                continue;
            }
            gradient[pg_offset_ + i] = interval_duration_ *
                pwl_value_slope(
                    generator_points_[i], x[pg_offset_ + i]).second /
                objective_scale_;
        }
        for (int i = 0; i < nd_; ++i) {
            if (load_epigraph_column_[i] >= 0) {
                continue;
            }
            const auto& load = data_.loads[i];
            gradient[demand_offset_ + i] = -interval_duration_ * load.pd_nominal *
                pwl_value_slope(
                    load_points_[i],
                    load.pd_nominal * x[demand_offset_ + i]).second /
                objective_scale_;
        }
        const double p_penalty =
            interval_duration_ * data_.p_delta_cost_approx / objective_scale_;
        const double q_penalty =
            interval_duration_ * data_.q_delta_cost_approx / objective_scale_;
        const double sm_penalty =
            interval_duration_ * data_.sm_cost_approx / objective_scale_;
        std::fill(
            gradient + p_delta_offset_,
            gradient + p_delta_offset_ + nb_, p_penalty);
        std::fill(
            gradient + q_delta_offset_,
            gradient + q_delta_offset_ + nb_, q_penalty);
        std::fill(
            gradient + sm_offset_,
            gradient + sm_offset_ + nl_, sm_penalty);
        for (const auto& [column, scale] : epigraph_weights_) {
            gradient[column] = interval_duration_ * scale / objective_scale_;
        }
        return true;
    }

    bool eval_g(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Index m,
        Ipopt::Number* constraints) override {
        if (n != variable_count_ || m != constraint_count_) {
            return false;
        }
        evaluate_constraints(x, constraints);
        return true;
    }

    bool eval_jac_g(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Index m,
        Ipopt::Index nele_jac,
        Ipopt::Index* i_row,
        Ipopt::Index* j_col,
        Ipopt::Number* values) override {
        if (n != variable_count_ || m != constraint_count_ ||
            nele_jac != static_cast<Ipopt::Index>(jacobian_rows_.size())) {
            return false;
        }
        if (values == nullptr) {
            std::copy(jacobian_rows_.begin(), jacobian_rows_.end(), i_row);
            std::copy(jacobian_columns_.begin(), jacobian_columns_.end(), j_col);
            return true;
        }
        fill_jacobian_values(x, values);
        return true;
    }

    bool eval_h(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool,
        Ipopt::Number, // The checked epigraph objective is linear: its Hessian is zero.
        Ipopt::Index m,
        const Ipopt::Number* lambda,
        bool,
        Ipopt::Index entries,
        Ipopt::Index* rows,
        Ipopt::Index* columns,
        Ipopt::Number* values) override {
        if (!exact_hessian_ || n != variable_count_ || m != constraint_count_ ||
            entries != static_cast<Ipopt::Index>(hessian_rows_.size())) return false;
        if (values == nullptr) {
            std::copy(hessian_rows_.begin(), hessian_rows_.end(), rows);
            std::copy(hessian_columns_.begin(), hessian_columns_.end(), columns);
            return true;
        }
        ++hessian_evaluations_;
        fill_hessian_values(x, lambda, values);
        return std::all_of(values, values + entries,
                           [](double value) { return std::isfinite(value); });
    }

    bool intermediate_callback(
        Ipopt::AlgorithmMode,
        Ipopt::Index iteration,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Number,
        Ipopt::Index,
        const Ipopt::IpoptData* ip_data,
        Ipopt::IpoptCalculatedQuantities* ip_cq) override {
        const auto callback_start = std::chrono::steady_clock::now();
        ++intermediate_callbacks_;
        try {
            std::vector<double> current_x(variable_count_);
            if (!get_curr_iterate(
                    ip_data, ip_cq, false,
                    variable_count_, current_x.data(), nullptr, nullptr,
                    constraint_count_, nullptr, nullptr)) {
                ++intermediate_capture_failures_;
                intermediate_capture_error_ =
                    "Ipopt did not expose the current unscaled primal iterate";
            } else {
                ++intermediate_iterates_retrieved_;
                SolveResult candidate;
                candidate.status = 0;
                candidate.iterations = static_cast<int>(iteration);
                candidate.state = state_from_x(current_x);
                candidate.objective = rebuild_economic_candidate(
                    data_, commitment_, original_corrective_base_, candidate.state);
                const auto validation = validate_economic_candidate(
                    data_, commitment_, original_corrective_base_, candidate.state);
                if (std::isfinite(candidate.objective) &&
                    validation.max_residual <= verification_tolerance_) {
                    ++intermediate_verified_candidates_;
                    if (!best_intermediate_found_ ||
                        candidate.objective >
                            best_intermediate_.objective + 1e-9) {
                        best_intermediate_found_ = true;
                        best_intermediate_iteration_ =
                            static_cast<int>(iteration);
                        best_intermediate_ = std::move(candidate);
                        best_intermediate_validation_ = validation;
                    }
                }
            }
        } catch (const std::exception& error) {
            ++intermediate_capture_failures_;
            intermediate_capture_error_ = error.what();
        } catch (...) {
            ++intermediate_capture_failures_;
            intermediate_capture_error_ =
                "unknown exception while checking an intermediate iterate";
        }
        intermediate_callback_seconds_ += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - callback_start).count();
        return true;
    }

    void finalize_solution(
        Ipopt::SolverReturn status,
        Ipopt::Index n,
        const Ipopt::Number* x,
        const Ipopt::Number*,
        const Ipopt::Number*,
        Ipopt::Index m,
        const Ipopt::Number* constraints,
        const Ipopt::Number*,
        Ipopt::Number objective,
        const Ipopt::IpoptData*,
        Ipopt::IpoptCalculatedQuantities*) override {
        solver_return_status_ = static_cast<int>(status);
        scaled_objective_ = objective;
        if (n != variable_count_ || m != constraint_count_ || x == nullptr) {
            return;
        }
        final_x_.assign(x, x + n);
        final_constraint_violation_ = constraints != nullptr
            ? constraint_violation(x, constraints)
            : kInfinity;
    }

    int variable_count() const { return variable_count_; }
    int constraint_count() const { return constraint_count_; }
    int pwl_epigraph_curve_count() const {
        return static_cast<int>(epigraph_weights_.size());
    }
    int pwl_epigraph_row_count() const {
        return static_cast<int>(epigraph_rows_.size());
    }
    int pwl_original_curve_count() const { return pwl_original_curve_count_; }
    int jacobian_nonzero_count() const {
        return static_cast<int>(jacobian_rows_.size());
    }
    int hessian_nonzero_count() const {
        return static_cast<int>(hessian_rows_.size());
    }
    int hessian_evaluations() const { return hessian_evaluations_; }
    int solver_return_status() const { return solver_return_status_; }
    double scaled_objective() const { return scaled_objective_; }
    double initial_constraint_violation() const {
        return initial_constraint_violation_;
    }
    double final_constraint_violation() const {
        return final_constraint_violation_;
    }
    bool has_final_x() const {
        return final_x_.size() == static_cast<std::size_t>(variable_count_);
    }

    AcState final_state() const {
        if (!has_final_x()) {
            throw std::runtime_error("sparse AC NLP did not return a primal point");
        }
        return state_from_x(final_x_);
    }

    int intermediate_callbacks() const { return intermediate_callbacks_; }
    int intermediate_iterates_retrieved() const {
        return intermediate_iterates_retrieved_;
    }
    int intermediate_verified_candidates() const {
        return intermediate_verified_candidates_;
    }
    int intermediate_capture_failures() const {
        return intermediate_capture_failures_;
    }
    double intermediate_callback_seconds() const {
        return intermediate_callback_seconds_;
    }
    bool best_intermediate_found() const {
        return best_intermediate_found_;
    }
    int best_intermediate_iteration() const {
        return best_intermediate_iteration_;
    }
    const SolveResult& best_intermediate() const {
        return best_intermediate_;
    }
    const ValidationReport& best_intermediate_validation() const {
        return best_intermediate_validation_;
    }
    const std::string& intermediate_capture_error() const {
        return intermediate_capture_error_;
    }

private:
    AcState state_from_x(const std::vector<double>& x) const {
        if (x.size() != static_cast<std::size_t>(variable_count_)) {
            throw std::runtime_error("sparse AC primal dimension mismatch");
        }
        AcState state = start_state_;
        state.vm.assign(x.begin() + vm_offset_, x.begin() + va_offset_);
        state.va.assign(x.begin() + va_offset_, x.begin() + pg_offset_);
        state.pg.assign(x.begin() + pg_offset_, x.begin() + qg_offset_);
        state.qg.assign(x.begin() + qg_offset_, x.begin() + demand_offset_);
        state.demand_factor.assign(
            x.begin() + demand_offset_, x.begin() + p_delta_offset_);
        state.p_delta.assign(
            x.begin() + p_delta_offset_, x.begin() + q_delta_offset_);
        state.q_delta.assign(
            x.begin() + q_delta_offset_, x.begin() + sm_offset_);
        state.sm_slack.assign(
            x.begin() + sm_offset_, x.begin() + sm_offset_ + nl_);
        return state;
    }

    void build_variables() {
        for (int bus = 0; bus < nb_; ++bus) {
            x_lower_[vm_offset_ + bus] = data_.buses[bus].vmin;
            x_upper_[vm_offset_ + bus] = data_.buses[bus].vmax;
            start_x_[vm_offset_ + bus] = std::clamp(
                start_state_.vm[bus], x_lower_[vm_offset_ + bus],
                x_upper_[vm_offset_ + bus]);
            start_x_[va_offset_ + bus] = start_state_.va[bus];
        }
        for (int i = 0; i < ng_; ++i) {
            if (commitment_[i] != 0 && commitment_[i] != 1) {
                throw std::runtime_error("sparse AC commitment is not binary");
            }
            const auto bounds = original_corrective_base_ != nullptr
                ? corrective_pg_bounds(data_.generators[i], commitment_[i],
                    original_corrective_base_->pg[i], data_.delta_r_ctg)
                : base_pg_bounds(data_.generators[i], commitment_[i], data_.delta_r);
            x_lower_[pg_offset_ + i] = bounds.first;
            x_upper_[pg_offset_ + i] = bounds.second;
            x_lower_[qg_offset_ + i] = commitment_[i]
                ? data_.generators[i].qmin : 0.0;
            x_upper_[qg_offset_ + i] = commitment_[i]
                ? data_.generators[i].qmax : 0.0;
            start_x_[pg_offset_ + i] = std::clamp(
                start_state_.pg[i], bounds.first, bounds.second);
            start_x_[qg_offset_ + i] = std::clamp(
                start_state_.qg[i], x_lower_[qg_offset_ + i],
                x_upper_[qg_offset_ + i]);
            if (commitment_[i]) {
                generator_points_[i] = active_pwl_points(
                    data_.generators[i].cost,
                    data_.generators[i].ncost,
                    bounds.first, bounds.second);
            }
        }
        for (int i = 0; i < nd_; ++i) {
            const auto bounds = original_corrective_base_ != nullptr
                ? corrective_load_bounds(data_.loads[i],
                    original_corrective_base_->demand_factor[i], data_.delta_r_ctg)
                : base_load_bounds(data_.loads[i], data_.delta_r);
            x_lower_[demand_offset_ + i] = bounds.first;
            x_upper_[demand_offset_ + i] = bounds.second;
            start_x_[demand_offset_ + i] = std::clamp(
                start_state_.demand_factor[i], bounds.first, bounds.second);
            load_points_[i] = active_pwl_points(
                data_.loads[i].cost, data_.loads[i].ncost,
                data_.loads[i].pd_min, data_.loads[i].pd_max);
        }
        for (int bus = 0; bus < nb_; ++bus) {
            x_lower_[p_delta_offset_ + bus] = 0.0;
            x_upper_[p_delta_offset_ + bus] = 0.5;
            x_lower_[q_delta_offset_ + bus] = 0.0;
            x_upper_[q_delta_offset_ + bus] = 0.5;
            start_x_[p_delta_offset_ + bus] = std::clamp(
                start_state_.p_delta[bus], 0.0, 0.5);
            start_x_[q_delta_offset_ + bus] = std::clamp(
                start_state_.q_delta[bus], 0.0, 0.5);
        }
        for (int branch = 0; branch < nl_; ++branch) {
            x_lower_[sm_offset_ + branch] = 0.0;
            x_upper_[sm_offset_ + branch] =
                data_.branches[branch].status == 0
                ? 0.0 : data_.sm_vio_limit;
            const double initial = start_state_.sm_slack.size() ==
                    data_.branches.size()
                ? start_state_.sm_slack[branch] : 0.0;
            start_x_[sm_offset_ + branch] = std::clamp(
                initial, x_lower_[sm_offset_ + branch],
                x_upper_[sm_offset_ + branch]);
        }
    }

    void build_pwl_epigraph() {
        const auto append_curve = [&](const std::vector<PwlPoint>& points,
                                      int physical_column, double power_factor,
                                      double sign, int& value_column) {
            std::vector<EpigraphSegment> segments;
            double scale = 1.0;
            if (!convex_epigraph_segments(points, sign, segments, scale)) {
                ++pwl_original_curve_count_;
                return;
            }
            value_column = variable_count_++;
            x_lower_.push_back(-kInfinity);
            x_upper_.push_back(kInfinity);
            start_x_.push_back(sign * pwl_value_slope(
                points, power_factor * start_x_[physical_column]).first / scale);
            epigraph_weights_.emplace_back(value_column, scale);
            for (const auto& segment : segments) {
                epigraph_rows_.push_back({
                    physical_column, value_column,
                    segment.slope * power_factor / scale,
                    segment.intercept / scale});
            }
        };
        for (int i = 0; i < ng_; ++i) {
            if (commitment_[i]) {
                append_curve(generator_points_[i], pg_offset_ + i, 1.0,
                             1.0, generator_epigraph_column_[i]);
            }
        }
        for (int i = 0; i < nd_; ++i) {
            // Minimize negative load benefit, a convex PWL function when
            // the supplied benefit curve is concave.
            append_curve(load_points_[i], demand_offset_ + i,
                         data_.loads[i].pd_nominal, -1.0,
                         load_epigraph_column_[i]);
        }
    }

    double objective_value(const double* x) const {
        double objective = 0.0;
        for (int i = 0; i < ng_; ++i) {
            if (commitment_[i] == 0) {
                continue;
            }
            if (generator_epigraph_column_[i] < 0) {
                objective += interval_duration_ * pwl_value_slope(
                    generator_points_[i], x[pg_offset_ + i]).first;
            }
            objective += interval_duration_ * data_.generators[i].oncost;
        }
        for (int i = 0; i < nd_; ++i) {
            if (load_epigraph_column_[i] < 0) {
                objective -= interval_duration_ * pwl_value_slope(
                    load_points_[i],
                    data_.loads[i].pd_nominal * x[demand_offset_ + i]).first;
            }
        }
        for (const auto& [column, scale] : epigraph_weights_) {
            objective += interval_duration_ * scale * x[column];
        }
        objective += interval_duration_ * data_.p_delta_cost_approx *
            std::accumulate(
                x + p_delta_offset_, x + p_delta_offset_ + nb_, 0.0);
        objective += interval_duration_ * data_.q_delta_cost_approx *
            std::accumulate(
                x + q_delta_offset_, x + q_delta_offset_ + nb_, 0.0);
        objective += interval_duration_ * data_.sm_cost_approx *
            std::accumulate(x + sm_offset_, x + sm_offset_ + nl_, 0.0);
        return objective / objective_scale_;
    }

    std::pair<double, double> network_injection(
        int bus,
        const double* x) const {
        const double vi = x[vm_offset_ + bus];
        const double ai = x[va_offset_ + bus];
        double p = 0.0;
        double q = 0.0;
        for (const auto& [other, admittance] : ybus_[bus]) {
            const double vj = x[vm_offset_ + other];
            const double angle = ai - x[va_offset_ + other];
            const double cosine = std::cos(angle);
            const double sine = std::sin(angle);
            const double g = admittance.real();
            const double b = admittance.imag();
            p += vi * vj * (g * cosine + b * sine);
            q += vi * vj * (g * sine - b * cosine);
        }
        return {p, q};
    }

    void evaluate_constraints(
        const double* x,
        double* constraints) const {
        for (int bus = 0; bus < nb_; ++bus) {
            auto [p, q] = network_injection(bus, x);
            for (int generator : data_.buses[bus].generators) {
                p -= x[pg_offset_ + generator];
                q -= x[qg_offset_ + generator];
            }
            for (int load : data_.buses[bus].loads) {
                p += data_.loads[load].pd_nominal *
                    x[demand_offset_ + load];
                q += data_.loads[load].qd_nominal *
                    x[demand_offset_ + load];
            }
            constraints[4 * bus] = p - x[p_delta_offset_ + bus];
            constraints[4 * bus + 1] = -p - x[p_delta_offset_ + bus];
            constraints[4 * bus + 2] = q - x[q_delta_offset_ + bus];
            constraints[4 * bus + 3] = -q - x[q_delta_offset_ + bus];
        }
        for (int position = 0;
             position < static_cast<int>(active_branches_.size());
             ++position) {
            const int branch_index = active_branches_[position];
            const auto& branch = data_.branches[branch_index];
            const auto flow = evaluate_branch(
                coefficients_[branch_index],
                x[vm_offset_ + branch.from],
                x[vm_offset_ + branch.to],
                x[va_offset_ + branch.from] -
                    x[va_offset_ + branch.to]);
            const double slack = x[sm_offset_ + branch_index];
            const double from_rating_voltage = branch.transformer
                ? 1.0 + slack
                : x[vm_offset_ + branch.from] + slack;
            const double to_rating_voltage = branch.transformer
                ? 1.0 + slack
                : x[vm_offset_ + branch.to] + slack;
            const double rating = original_corrective_base_ != nullptr
                ? branch.rate_c : branch.rate_a;
            const double rating_squared = rating * rating;
            constraints[thermal_row_offset_ + 2 * position] =
                flow.flow[0] * flow.flow[0] +
                flow.flow[1] * flow.flow[1] -
                rating_squared * from_rating_voltage * from_rating_voltage;
            constraints[thermal_row_offset_ + 2 * position + 1] =
                flow.flow[2] * flow.flow[2] +
                flow.flow[3] * flow.flow[3] -
                rating_squared * to_rating_voltage * to_rating_voltage;
        }
        for (int position = 0;
             position < static_cast<int>(angle_branches_.size());
             ++position) {
            const auto& branch = data_.branches[angle_branches_[position]];
            constraints[angle_row_offset_ + position] =
                x[va_offset_ + branch.from] -
                x[va_offset_ + branch.to];
        }
        for (int position = 0;
             position < static_cast<int>(reference_buses_.size());
             ++position) {
            constraints[reference_row_offset_ + position] =
                x[va_offset_ + reference_buses_[position]];
        }
        for (int i = 0; i < static_cast<int>(epigraph_rows_.size()); ++i) {
            const auto& row = epigraph_rows_[i];
            constraints[epigraph_row_offset_ + i] =
                row.slope * x[row.physical_column] + row.intercept -
                x[row.value_column];
        }
    }

    void append_structure(int row, int column) {
        jacobian_rows_.push_back(row);
        jacobian_columns_.push_back(column);
    }

    template <typename Emit>
    void visit_balance_jacobian(int bus, const double* x, Emit&& emit) const {
        const double vi = x == nullptr ? 1.0 : x[vm_offset_ + bus];
        const double ai = x == nullptr ? 0.0 : x[va_offset_ + bus];
        double dp_dvi = 0.0;
        double dq_dvi = 0.0;
        double dp_dai = 0.0;
        double dq_dai = 0.0;
        for (const auto& [other, admittance] : ybus_[bus]) {
            const double vj = x == nullptr ? 1.0 : x[vm_offset_ + other];
            const double aj = x == nullptr ? 0.0 : x[va_offset_ + other];
            const double angle = ai - aj;
            const double cosine = std::cos(angle);
            const double sine = std::sin(angle);
            const double g = admittance.real();
            const double b = admittance.imag();
            if (other == bus) {
                dp_dvi += 2.0 * vi * g;
                dq_dvi -= 2.0 * vi * b;
                continue;
            }
            const double active = g * cosine + b * sine;
            const double reactive = g * sine - b * cosine;
            dp_dvi += vj * active;
            dq_dvi += vj * reactive;
            dp_dai -= vi * vj * reactive;
            dq_dai += vi * vj * active;
        }
        const auto emit_network = [&] (int row, double sign, bool reactive) {
            emit(row, vm_offset_ + bus, sign *
                (reactive ? dq_dvi : dp_dvi));
            emit(row, va_offset_ + bus, sign *
                (reactive ? dq_dai : dp_dai));
            for (const auto& [other, admittance] : ybus_[bus]) {
                if (other == bus) {
                    continue;
                }
                const double vj = x == nullptr ? 1.0 : x[vm_offset_ + other];
                const double aj = x == nullptr ? 0.0 : x[va_offset_ + other];
                const double angle = ai - aj;
                const double cosine = std::cos(angle);
                const double sine = std::sin(angle);
                const double g = admittance.real();
                const double b = admittance.imag();
                const double active = g * cosine + b * sine;
                const double reactive_term = g * sine - b * cosine;
                const double d_voltage = vi *
                    (reactive ? reactive_term : active);
                const double d_angle = reactive
                    ? -vi * vj * active
                    : vi * vj * reactive_term;
                emit(row, vm_offset_ + other, sign * d_voltage);
                emit(row, va_offset_ + other, sign * d_angle);
            }
        };
        const int p_positive = 4 * bus;
        const int p_negative = p_positive + 1;
        const int q_positive = p_positive + 2;
        const int q_negative = p_positive + 3;
        emit_network(p_positive, 1.0, false);
        emit_network(p_negative, -1.0, false);
        emit_network(q_positive, 1.0, true);
        emit_network(q_negative, -1.0, true);
        for (int generator : data_.buses[bus].generators) {
            emit(p_positive, pg_offset_ + generator, -1.0);
            emit(p_negative, pg_offset_ + generator, 1.0);
            emit(q_positive, qg_offset_ + generator, -1.0);
            emit(q_negative, qg_offset_ + generator, 1.0);
        }
        for (int load : data_.buses[bus].loads) {
            emit(p_positive, demand_offset_ + load,
                 data_.loads[load].pd_nominal);
            emit(p_negative, demand_offset_ + load,
                 -data_.loads[load].pd_nominal);
            emit(q_positive, demand_offset_ + load,
                 data_.loads[load].qd_nominal);
            emit(q_negative, demand_offset_ + load,
                 -data_.loads[load].qd_nominal);
        }
        emit(p_positive, p_delta_offset_ + bus, -1.0);
        emit(p_negative, p_delta_offset_ + bus, -1.0);
        emit(q_positive, q_delta_offset_ + bus, -1.0);
        emit(q_negative, q_delta_offset_ + bus, -1.0);
    }

    template <typename Emit>
    void visit_thermal_jacobian(
        int position,
        const double* x,
        Emit&& emit) const {
        const int branch_index = active_branches_[position];
        const auto& branch = data_.branches[branch_index];
        const double vm_from = x == nullptr ? 1.0 : x[vm_offset_ + branch.from];
        const double vm_to = x == nullptr ? 1.0 : x[vm_offset_ + branch.to];
        const double va_from = x == nullptr ? 0.0 : x[va_offset_ + branch.from];
        const double va_to = x == nullptr ? 0.0 : x[va_offset_ + branch.to];
        const double slack = x == nullptr ? 0.0 : x[sm_offset_ + branch_index];
        const auto flow = evaluate_branch(
            coefficients_[branch_index], vm_from, vm_to, va_from - va_to);
        const double rating = original_corrective_base_ != nullptr
            ? branch.rate_c : branch.rate_a;
        const double rating_squared = rating * rating;
        for (int side = 0; side < 2; ++side) {
            const int active_component = side == 0 ? 0 : 2;
            const int reactive_component = active_component + 1;
            const int row = thermal_row_offset_ + 2 * position + side;
            for (int variable = 0; variable < 4; ++variable) {
                double derivative =
                    2.0 * flow.flow[active_component] *
                        flow.derivative[active_component][variable] +
                    2.0 * flow.flow[reactive_component] *
                        flow.derivative[reactive_component][variable];
                const bool terminal_voltage = !branch.transformer &&
                    ((side == 0 && variable == 0) ||
                     (side == 1 && variable == 1));
                if (terminal_voltage) {
                    const double voltage = side == 0 ? vm_from : vm_to;
                    derivative -= 2.0 * rating_squared *
                        (voltage + slack);
                }
                const int column = variable == 0
                    ? vm_offset_ + branch.from
                    : variable == 1
                        ? vm_offset_ + branch.to
                        : variable == 2
                            ? va_offset_ + branch.from
                            : va_offset_ + branch.to;
                emit(row, column, derivative);
            }
            const double rating_voltage = branch.transformer
                ? 1.0 + slack
                : (side == 0 ? vm_from : vm_to) + slack;
            emit(row, sm_offset_ + branch_index,
                 -2.0 * rating_squared * rating_voltage);
        }
    }

    void build_jacobian_structure() {
        const auto emit = [&] (int row, int column, double) {
            append_structure(row, column);
        };
        for (int bus = 0; bus < nb_; ++bus) {
            visit_balance_jacobian(bus, nullptr, emit);
        }
        for (int position = 0;
             position < static_cast<int>(active_branches_.size());
             ++position) {
            visit_thermal_jacobian(position, nullptr, emit);
        }
        for (int position = 0;
             position < static_cast<int>(angle_branches_.size());
             ++position) {
            const auto& branch = data_.branches[angle_branches_[position]];
            append_structure(angle_row_offset_ + position,
                             va_offset_ + branch.from);
            append_structure(angle_row_offset_ + position,
                             va_offset_ + branch.to);
        }
        for (int position = 0;
             position < static_cast<int>(reference_buses_.size());
             ++position) {
            append_structure(reference_row_offset_ + position,
                             va_offset_ + reference_buses_[position]);
        }
        for (int i = 0; i < static_cast<int>(epigraph_rows_.size()); ++i) {
            append_structure(epigraph_row_offset_ + i,
                             epigraph_rows_[i].physical_column);
            append_structure(epigraph_row_offset_ + i,
                             epigraph_rows_[i].value_column);
        }
    }

    void fill_jacobian_values(const double* x, double* values) const {
        std::size_t next = 0;
        const auto emit = [&] (int, int, double value) {
            values[next++] = value;
        };
        for (int bus = 0; bus < nb_; ++bus) {
            visit_balance_jacobian(bus, x, emit);
        }
        for (int position = 0;
             position < static_cast<int>(active_branches_.size());
             ++position) {
            visit_thermal_jacobian(position, x, emit);
        }
        for (int position = 0;
             position < static_cast<int>(angle_branches_.size());
             ++position) {
            emit(0, 0, 1.0);
            emit(0, 0, -1.0);
        }
        for (std::size_t position = 0;
             position < reference_buses_.size(); ++position) {
            emit(0, 0, 1.0);
        }
        for (const auto& row : epigraph_rows_) {
            emit(0, 0, row.slope);
            emit(0, 0, -1.0);
        }
        if (next != jacobian_rows_.size()) {
            throw std::runtime_error("sparse AC Jacobian fill count mismatch");
        }
    }

    void build_hessian_structure() {
        std::map<std::pair<int, int>, int> positions;
        const auto entry = [&](int first, int second) {
            const auto key = std::make_pair(std::max(first, second), std::min(first, second));
            const auto [it, inserted] = positions.emplace(key, static_cast<int>(positions.size()));
            if (inserted) {
                hessian_rows_.push_back(key.first);
                hessian_columns_.push_back(key.second);
            }
            return it->second;
        };
        for (int bus = 0; bus < nb_; ++bus) {
            hessian_vm_diagonal_.push_back(entry(vm_offset_ + bus, vm_offset_ + bus));
        }
        for (int branch_index : active_branches_) {
            const auto& branch = data_.branches[branch_index];
            HessianBranch pattern;
            pattern.columns = {vm_offset_ + branch.from, vm_offset_ + branch.to,
                               va_offset_ + branch.from, va_offset_ + branch.to,
                               sm_offset_ + branch_index};
            pattern.entries.fill(-1);
            int slot = 0;
            for (int i = 0; i < 5; ++i) {
                for (int j = 0; j <= i; ++j, ++slot) {
                    // Slack-angle curvature is identically zero. Transformer
                    // rating slack has no voltage cross term either.
                    if (i == 4 && j != 4 && (j >= 2 || branch.transformer)) continue;
                    pattern.entries[slot] = entry(pattern.columns[i], pattern.columns[j]);
                    // T' H T: a self-loop maps two local off-diagonal terms
                    // onto one global diagonal and must count both.
                    pattern.weights[slot] = i != j && pattern.columns[i] == pattern.columns[j]
                        ? 2.0 : 1.0;
                }
            }
            hessian_branch_patterns_.push_back(pattern);
        }
    }

    void fill_hessian_values(const double* x, const double* lambda, double* values) const {
        std::fill(values, values + hessian_rows_.size(), 0.0);
        for (int i = 0; i < static_cast<int>(data_.shunts.size()); ++i) {
            const auto& shunt = data_.shunts[i];
            const int row = 4 * shunt.bus;
            values[hessian_vm_diagonal_[shunt.bus]] += 2.0 * (
                (lambda[row] - lambda[row + 1]) * shunt.gs -
                (lambda[row + 2] - lambda[row + 3]) *
                    effective_shunt_susceptance(data_, start_state_, i));
        }
        for (int position = 0; position < static_cast<int>(active_branches_.size()); ++position) {
            const int branch_index = active_branches_[position];
            const auto& branch = data_.branches[branch_index];
            const auto& pattern = hessian_branch_patterns_[position];
            const double vf = x[pattern.columns[0]], vt = x[pattern.columns[1]];
            const double angle = x[pattern.columns[2]] - x[pattern.columns[3]];
            const auto flow = evaluate_branch(coefficients_[branch_index], vf, vt, angle);
            const auto second = branch_flow_hessians(coefficients_[branch_index], vf, vt, angle);
            const std::array<double, 4> balance_weight = {
                lambda[4 * branch.from] - lambda[4 * branch.from + 1],
                lambda[4 * branch.from + 2] - lambda[4 * branch.from + 3],
                lambda[4 * branch.to] - lambda[4 * branch.to + 1],
                lambda[4 * branch.to + 2] - lambda[4 * branch.to + 3]};
            std::array<std::array<double, 5>, 5> local{};
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j <= i; ++j) {
                    for (int k = 0; k < 4; ++k) local[i][j] += balance_weight[k] * second[k][i][j];
                }
            }
            const double rating = original_corrective_base_ != nullptr ? branch.rate_c : branch.rate_a;
            for (int side = 0; side < 2; ++side) {
                const double weight = 2.0 * lambda[thermal_row_offset_ + 2 * position + side];
                const int p = 2 * side, q = p + 1;
                for (int i = 0; i < 4; ++i) {
                    for (int j = 0; j <= i; ++j) {
                        local[i][j] += weight * (
                            flow.derivative[p][i] * flow.derivative[p][j] + flow.flow[p] * second[p][i][j] +
                            flow.derivative[q][i] * flow.derivative[q][j] + flow.flow[q] * second[q][i][j]);
                    }
                }
                local[4][4] -= weight * rating * rating;
                if (!branch.transformer) {
                    local[side][side] -= weight * rating * rating;
                    local[4][side] -= weight * rating * rating;
                }
            }
            int slot = 0;
            for (int i = 0; i < 5; ++i) {
                for (int j = 0; j <= i; ++j, ++slot) {
                    if (pattern.entries[slot] >= 0) {
                        values[pattern.entries[slot]] += pattern.weights[slot] * local[i][j];
                    }
                }
            }
        }
    }

    double constraint_violation(
        const double* x,
        const double* constraints) const {
        double violation = 0.0;
        for (int i = 0; i < variable_count_; ++i) {
            violation = std::max({
                violation,
                x_lower_[i] - x[i],
                x[i] - x_upper_[i],
            });
        }
        for (int row = 0; row < angle_row_offset_; ++row) {
            violation = std::max(violation, constraints[row]);
        }
        for (int i = 0; i < static_cast<int>(angle_branches_.size()); ++i) {
            const auto& branch = data_.branches[angle_branches_[i]];
            const double value = constraints[angle_row_offset_ + i];
            violation = std::max({
                violation, branch.angmin - value, value - branch.angmax});
        }
        for (int i = 0; i < static_cast<int>(reference_buses_.size()); ++i) {
            violation = std::max(
                violation,
                std::abs(constraints[reference_row_offset_ + i]));
        }
        for (int row = epigraph_row_offset_; row < constraint_count_; ++row) {
            violation = std::max(violation, constraints[row]);
        }
        return violation;
    }

    const CaseData& data_;
    std::vector<int> commitment_;
    AcState start_state_;
    const AcState* original_corrective_base_{};
    double interval_duration_{};
    double verification_tolerance_{};
    bool exact_hessian_{};
    int nb_{};
    int ng_{};
    int nd_{};
    int nl_{};
    int vm_offset_{};
    int va_offset_{};
    int pg_offset_{};
    int qg_offset_{};
    int demand_offset_{};
    int p_delta_offset_{};
    int q_delta_offset_{};
    int sm_offset_{};
    int variable_count_{};
    int balance_row_offset_{};
    int thermal_row_offset_{};
    int angle_row_offset_{};
    int reference_row_offset_{};
    int epigraph_row_offset_{};
    int constraint_count_{};
    int pwl_original_curve_count_{};
    double objective_scale_{1.0};
    std::vector<double> x_lower_;
    std::vector<double> x_upper_;
    std::vector<double> start_x_;
    std::vector<std::vector<PwlPoint>> generator_points_;
    std::vector<std::vector<PwlPoint>> load_points_;
    struct EpigraphRow {
        int physical_column{};
        int value_column{};
        double slope{};
        double intercept{};
    };
    std::vector<int> generator_epigraph_column_;
    std::vector<int> load_epigraph_column_;
    std::vector<EpigraphRow> epigraph_rows_;
    std::vector<std::pair<int, double>> epigraph_weights_;
    YRows ybus_;
    std::vector<BranchCoefficients> coefficients_;
    std::vector<int> active_branches_;
    std::vector<int> angle_branches_;
    std::vector<int> reference_buses_;
    std::vector<Ipopt::Index> jacobian_rows_;
    std::vector<Ipopt::Index> jacobian_columns_;
    struct HessianBranch {
        std::array<int, 5> columns{};
        std::array<int, 15> entries{};
        std::array<double, 15> weights{};
    };
    std::vector<Ipopt::Index> hessian_rows_;
    std::vector<Ipopt::Index> hessian_columns_;
    std::vector<int> hessian_vm_diagonal_;
    std::vector<HessianBranch> hessian_branch_patterns_;
    int hessian_evaluations_{};
    std::vector<double> final_x_;
    int solver_return_status_{-99};
    double scaled_objective_{};
    double initial_constraint_violation_{};
    double final_constraint_violation_{kInfinity};
    int intermediate_callbacks_{};
    int intermediate_iterates_retrieved_{};
    int intermediate_verified_candidates_{};
    int intermediate_capture_failures_{};
    int best_intermediate_iteration_{-1};
    double intermediate_callback_seconds_{};
    bool best_intermediate_found_{};
    SolveResult best_intermediate_;
    ValidationReport best_intermediate_validation_;
    std::string intermediate_capture_error_;
};

std::string application_status_string(Ipopt::ApplicationReturnStatus status) {
    switch (status) {
        case Ipopt::Solve_Succeeded: return "Solve_Succeeded";
        case Ipopt::Solved_To_Acceptable_Level:
            return "Solved_To_Acceptable_Level";
        case Ipopt::Infeasible_Problem_Detected:
            return "Infeasible_Problem_Detected";
        case Ipopt::Search_Direction_Becomes_Too_Small:
            return "Search_Direction_Becomes_Too_Small";
        case Ipopt::Maximum_Iterations_Exceeded:
            return "Maximum_Iterations_Exceeded";
        case Ipopt::Maximum_CpuTime_Exceeded:
            return "Maximum_CpuTime_Exceeded";
        case Ipopt::Maximum_WallTime_Exceeded:
            return "Maximum_WallTime_Exceeded";
        case Ipopt::Invalid_Number_Detected:
            return "Invalid_Number_Detected";
        default: return "Ipopt_status_" + std::to_string(static_cast<int>(status));
    }
}

}  // namespace

void run_sparse_ac_pwl_epigraph_regression(
    const CaseData& data, const std::vector<int>& commitment,
    const AcState& start) {
    const auto require = [](bool condition, const char* message) {
        if (!condition) {
            throw std::runtime_error(std::string("PWL epigraph regression: ") + message);
        }
    };
    for (double sign : {1.0, -1.0}) {
        const std::vector<PwlPoint> points{
            {2.0, sign * 3.0}, {3.0, sign * 5.0}, {5.0, sign * 13.0}};
        std::vector<EpigraphSegment> segments;
        double scale = 0.0;
        require(convex_epigraph_segments(points, sign, segments, scale),
                "valid convex signed cost rejected");
        for (double power : {2.0, 2.5, 3.0, 4.0, 5.0}) {
            double envelope = -kInfinity;
            for (const auto& segment : segments) {
                envelope = std::max(envelope,
                    segment.slope * power + segment.intercept);
            }
            require(std::abs(envelope - sign *
                pwl_value_slope(points, power).first) < 1e-12,
                "source cost changed at a breakpoint or segment interior");
        }
    }
    std::vector<EpigraphSegment> rejected;
    double scale = 0.0;
    require(!convex_epigraph_segments({{2.0, 0.0}, {3.0, 2.0}, {4.0, 3.0}},
                1.0, rejected, scale) && rejected.empty(),
            "nonconvex source curve was convexified");
    require(!convex_epigraph_segments({{2.0, 0.0}, {2.0, 1.0}},
                1.0, rejected, scale), "duplicate abscissa accepted");

    SparseAcEconomicNlp original(data, commitment, start, 1e-5, false);
    SparseAcEconomicNlp epigraph(data, commitment, start, 1e-5, true);
    Ipopt::Index n0, m0, j0, h0, n1, m1, j1, h1;
    Ipopt::TNLP::IndexStyleEnum style;
    original.get_nlp_info(n0, m0, j0, h0, style);
    epigraph.get_nlp_info(n1, m1, j1, h1, style);
    require(n1 > n0 && m1 > m0 && epigraph.pwl_epigraph_curve_count() > 0,
            "tiny fixture did not exercise the reformulation");
    std::vector<double> x0(n0), x1(n1), lo0(n0), hi0(n0), lo1(n1), hi1(n1);
    std::vector<double> gl0(m0), gu0(m0), gl1(m1), gu1(m1), g0(m0), g1(m1);
    original.get_bounds_info(n0, lo0.data(), hi0.data(), m0, gl0.data(), gu0.data());
    epigraph.get_bounds_info(n1, lo1.data(), hi1.data(), m1, gl1.data(), gu1.data());
    original.get_starting_point(n0, true, x0.data(), false, nullptr, nullptr,
                               m0, false, nullptr);
    epigraph.get_starting_point(n1, true, x1.data(), false, nullptr, nullptr,
                               m1, false, nullptr);
    require(std::equal(lo0.begin(), lo0.end(), lo1.begin()) &&
            std::equal(hi0.begin(), hi0.end(), hi1.begin()) &&
            std::equal(x0.begin(), x0.end(), x1.begin()) &&
            std::equal(gl0.begin(), gl0.end(), gl1.begin()) &&
            std::equal(gu0.begin(), gu0.end(), gu1.begin()),
            "source bounds, PMIN, physical start or row bounds changed");
    original.eval_g(n0, x0.data(), true, m0, g0.data());
    epigraph.eval_g(n1, x1.data(), true, m1, g1.data());
    require(std::equal(g0.begin(), g0.end(), g1.begin()),
            "physical constraints changed");
    for (int row = m0; row < m1; ++row) {
        require(g1[row] <= 1e-12, "mapped source cost is not epigraph feasible");
    }
    double f0 = 0.0, f1 = 0.0;
    original.eval_f(n0, x0.data(), true, f0);
    epigraph.eval_f(n1, x1.data(), true, f1);
    require(std::abs(f0 - f1) < 1e-12, "mapped objective changed");

    std::vector<Ipopt::Index> jr(j1), jc(j1);
    std::vector<double> values(j1), gradient(n1);
    epigraph.eval_jac_g(n1, x1.data(), true, m1, j1, jr.data(), jc.data(), nullptr);
    epigraph.eval_jac_g(n1, x1.data(), true, m1, j1, nullptr, nullptr, values.data());
    epigraph.eval_grad_f(n1, x1.data(), true, gradient.data());
    constexpr double step = 1e-6;
    for (int column = 0; column < n1; ++column) {
        auto lower = x1, upper = x1;
        lower[column] -= step;
        upper[column] += step;
        std::vector<double> lower_g(m1), upper_g(m1);
        epigraph.eval_g(n1, lower.data(), true, m1, lower_g.data());
        epigraph.eval_g(n1, upper.data(), true, m1, upper_g.data());
        for (int row = m0; row < m1; ++row) {
            double analytic = 0.0;
            for (int entry = j0; entry < j1; ++entry) {
                if (jr[entry] == row && jc[entry] == column) {
                    analytic += values[entry];
                }
            }
            require(std::abs(analytic - (upper_g[row] - lower_g[row]) /
                    (2.0 * step)) < 1e-7, "cost-row Jacobian mismatch");
        }
        double lower_f = 0.0, upper_f = 0.0;
        epigraph.eval_f(n1, lower.data(), true, lower_f);
        epigraph.eval_f(n1, upper.data(), true, upper_f);
        require(std::abs(gradient[column] - (upper_f - lower_f) /
                (2.0 * step)) < 1e-7, "objective gradient mismatch");
    }
}

void run_sparse_ac_corrective_reference_regression(
    const CaseData& data, const std::vector<int>& commitment,
    const AcState& original_base) {
    const auto require = [](bool condition, const char* message) {
        if (!condition) {
            throw std::runtime_error(std::string("corrective reference regression: ") + message);
        }
    };
    require(!data.generators.empty() && !data.loads.empty(), "fixture lacks controls");
    const auto frozen_base = ac_state_to_json(original_base);
    auto fixture = data;
    fixture.delta = 1.75;
    fixture.delta_ctg = 0.25;
    fixture.delta_r_ctg = 0.125;
    auto& generator = fixture.generators[0];
    generator.pmin = std::min(0.2, 0.5 * original_base.pg[0]);
    generator.pg_prev = original_base.pg[0] + 0.4;
    generator.prdmaxctg = 0.16;
    generator.prumaxctg = 0.24;
    auto& load = fixture.loads[0];
    load.tmin = 0.5;
    load.tmax = 1.5;
    load.pd_prev = load.pd_nominal * (original_base.demand_factor[0] + 0.3);
    load.prdmaxctg = 0.04;
    load.prumaxctg = 0.08;
    for (auto& branch : fixture.branches) {
        branch.rate_c = 1.5 * branch.rate_a;
    }
    AcState start = original_base;
    const double native_objective = rebuild_common_corrective_reference_state(
        fixture, original_base, commitment, start);
    SparseAcEconomicNlp model(fixture, commitment, start, 1e-5, true, &original_base);
    Ipopt::Index n, m, entries, h;
    Ipopt::TNLP::IndexStyleEnum style;
    model.get_nlp_info(n, m, entries, h, style);
    std::vector<double> x(n), lo(n), hi(n), gl(m), gu(m), g(m);
    model.get_bounds_info(n, lo.data(), hi.data(), m, gl.data(), gu.data());
    model.get_starting_point(n, true, x.data(), false, nullptr, nullptr,
                             m, false, nullptr);
    const int nb = static_cast<int>(fixture.buses.size());
    const int ng = static_cast<int>(fixture.generators.size());
    const int pg_column = 2 * nb;
    const int demand_column = 2 * nb + 2 * ng;
    require(generator.pmin > 0.0 &&
            std::abs(lo[pg_column] - std::max(generator.pmin,
                original_base.pg[0] - 0.125 * 0.16)) < 1e-12 &&
            std::abs(hi[pg_column] - std::min(generator.pmax,
                original_base.pg[0] + 0.125 * 0.24)) < 1e-12,
            "generator bounds/PMIN do not use original base and corrective ramps");
    require(std::abs(lo[demand_column] - std::max(load.tmin,
                original_base.demand_factor[0] - 0.125 * 0.04 / load.pd_nominal)) < 1e-12 &&
            std::abs(hi[demand_column] - std::min(load.tmax,
                original_base.demand_factor[0] + 0.125 * 0.08 / load.pd_nominal)) < 1e-12,
            "load bounds do not use original base and corrective ramps");
    require(corrective_pg_bounds(generator, 0, original_base.pg[0], 0.125) ==
                std::pair<double, double>{0.0, 0.0}, "uncommitted generator made available");
    double objective = 0.0;
    model.eval_f(n, x.data(), true, objective);
    const double objective_scale = std::max({1.0, fixture.p_delta_cost_approx,
        fixture.q_delta_cost_approx, fixture.sm_cost_approx});
    require(std::abs(-objective * objective_scale - native_objective) < 1e-7,
            "reference objective differs from independent corrective rebuild");
    model.eval_g(n, x.data(), true, m, g.data());
    int active_position = 0;
    for (std::size_t i = 0; i < fixture.branches.size(); ++i) {
        const auto& branch = fixture.branches[i];
        if (branch.status == 0) continue;
        const double limit = branch.rate_c * (branch.transformer
            ? 1.0 + start.sm_slack[i] : start.vm[branch.from] + start.sm_slack[i]);
        const double expected = start.pf[i] * start.pf[i] + start.qf[i] * start.qf[i]
            - limit * limit;
        require(std::abs(g[4 * nb + 2 * active_position] - expected) < 1e-10,
                "reference does not use source corrective thermal rating");
        ++active_position;
    }
    // Moving a candidate does not move any original-base-dependent bound.
    auto shifted = start;
    shifted.pg[0] += 0.01;
    shifted.demand_factor[0] += 0.005;
    SparseAcEconomicNlp shifted_model(
        fixture, commitment, shifted, 1e-5, true, &original_base);
    std::vector<double> shifted_lo(n), shifted_hi(n), shifted_gl(m), shifted_gu(m);
    shifted_model.get_bounds_info(n, shifted_lo.data(), shifted_hi.data(), m,
                                  shifted_gl.data(), shifted_gu.data());
    require(lo == shifted_lo && hi == shifted_hi && gl == shifted_gl && gu == shifted_gu,
            "candidate changed the corrective bound anchor");
    auto invalid = start;
    invalid.pg[0] = hi[pg_column] + 0.01;
    rebuild_common_corrective_reference_state(fixture, original_base, commitment, invalid);
    require(validate_economic_candidate(fixture, commitment, &original_base, invalid)
                .max_residual > 1e-5, "out-of-ramp candidate accepted");

    std::vector<Ipopt::Index> jr(entries), jc(entries);
    std::vector<double> values(entries), gradient(n);
    model.eval_jac_g(n, x.data(), true, m, entries, jr.data(), jc.data(), nullptr);
    model.eval_jac_g(n, x.data(), true, m, entries, nullptr, nullptr, values.data());
    model.eval_grad_f(n, x.data(), true, gradient.data());
    constexpr double step = 1e-6;
    for (int column = 0; column < n; ++column) {
        auto lower = x, upper = x;
        lower[column] -= step;
        upper[column] += step;
        std::vector<double> lower_g(m), upper_g(m), analytic(m, 0.0);
        model.eval_g(n, lower.data(), true, m, lower_g.data());
        model.eval_g(n, upper.data(), true, m, upper_g.data());
        for (int entry = 0; entry < entries; ++entry) {
            if (jc[entry] == column) analytic[jr[entry]] += values[entry];
        }
        for (int row = 0; row < m; ++row) {
            require(std::abs(analytic[row] - (upper_g[row] - lower_g[row]) /
                (2.0 * step)) < 1e-5, "corrective constraint Jacobian mismatch");
        }
        double lower_f = 0.0, upper_f = 0.0;
        model.eval_f(n, lower.data(), true, lower_f);
        model.eval_f(n, upper.data(), true, upper_f);
        require(std::abs(gradient[column] - (upper_f - lower_f) /
                (2.0 * step)) < 1e-7, "corrective objective gradient mismatch");
    }
    require(frozen_base == ac_state_to_json(original_base), "original base mutated");
}

void run_sparse_ac_exact_hessian_regression(
    const CaseData& data, const std::vector<int>& commitment, const AcState& start) {
    const auto require = [](bool condition, const char* message) {
        if (!condition) throw std::runtime_error(std::string("exact Hessian regression: ") + message);
    };
    require(data.buses.size() >= 2 && data.branches.size() >= 2, "fixture too small");
    auto fixture = data;
    fixture.delta = 1.75;
    fixture.delta_ctg = 0.25;
    fixture.branches[0].tap = 1.07;
    fixture.branches[0].shift = 0.13;
    fixture.branches[0].g_fr = 0.023;
    fixture.branches[0].b_fr = -0.017;
    fixture.branches[0].g_to = 0.009;
    fixture.branches[0].b_to = 0.021;
    fixture.branches[1].transformer = true;
    fixture.branches[1].tap = 0.93;
    fixture.branches[1].shift = -0.17;
    fixture.branches[1].g_fr = 0.019;
    fixture.branches[1].b_fr = 0.031;
    fixture.branches[1].g_to = -0.007;
    fixture.branches[1].b_to = -0.023;
    for (auto& branch : fixture.branches) branch.rate_c = 1.7 * branch.rate_a;
    Shunt shunt;
    shunt.bus = 0;
    shunt.gs = 0.037;
    shunt.bs = -0.051;
    fixture.shunts.push_back(shunt);
    for (const auto& branch : fixture.branches) {
        const Complex series = 1.0 / Complex(branch.r, branch.x);
        const Complex from_shunt(branch.g_fr, branch.b_fr), to_shunt(branch.g_to, branch.b_to);
        const double tap2 = branch.tap * branch.tap;
        const Complex rotation = std::polar(1.0, branch.shift);
        const Complex yff = branch.transformer ? series / tap2 + from_shunt
                                               : (series + from_shunt) / tap2;
        const Complex yft = -series * rotation / branch.tap;
        const Complex ytf = -series * std::conj(rotation) / branch.tap;
        const Complex ytt = series + to_shunt;
        const Complex vf = std::polar(0.97, 0.07), vt = std::polar(1.04, -0.11);
        const Complex sf = vf * std::conj(yff * vf + yft * vt);
        const Complex st = vt * std::conj(ytf * vf + ytt * vt);
        const auto analytic = evaluate_branch(branch_coefficients(branch), 0.97, 1.04, 0.18);
        const std::array<double, 4> complex_power = {sf.real(), sf.imag(), st.real(), st.imag()};
        for (int k = 0; k < 4; ++k) {
            require(std::abs(analytic.flow[k] - complex_power[k]) < 1e-12,
                    "terminal power differs from complex-admittance calculation");
        }
    }
    // Derivative probes need not be physically feasible, and never solve a
    // production case. Mutations here affect this copied tiny fixture only.
    for (int topology = 0; topology < 3; ++topology) {
        auto topology_data = fixture;
        if (topology == 1) topology_data.branches[1].to = topology_data.branches[1].from;
        if (topology == 2) topology_data.branches[1].status = 0;
        for (bool corrective : {false, true}) {
            SparseAcEconomicNlp model(topology_data, commitment, start, 1e-5,
                                       true, corrective ? &start : nullptr, true);
            SparseAcEconomicNlp approximate(topology_data, commitment, start, 1e-5,
                                             true, corrective ? &start : nullptr, false);
            Ipopt::Index n, m, nj, nh;
            Ipopt::TNLP::IndexStyleEnum style;
            model.get_nlp_info(n, m, nj, nh, style);
            require(nh > 0, "empty Hessian pattern");
            std::vector<double> x(n), lo(n), hi(n), gl(m), gu(m);
            model.get_starting_point(n, true, x.data(), false, nullptr, nullptr, m, false, nullptr);
            model.get_bounds_info(n, lo.data(), hi.data(), m, gl.data(), gu.data());
            std::vector<double> old_lo(n), old_hi(n), old_gl(m), old_gu(m);
            approximate.get_bounds_info(n, old_lo.data(), old_hi.data(), m, old_gl.data(), old_gu.data());
            require(lo == old_lo && hi == old_hi && gl == old_gl && gu == old_gu,
                    "derivative option changed bounds");
            x[0] = 0.97; x[1] = 1.04;
            x[topology_data.buses.size()] = 0.07;
            x[topology_data.buses.size() + 1] = -0.11;
            std::vector<double> g(m), old_g(m);
            model.eval_g(n, x.data(), true, m, g.data());
            approximate.eval_g(n, x.data(), true, m, old_g.data());
            require(g == old_g, "derivative option changed constraint values");
            std::vector<Ipopt::Index> jr(nj), jc(nj), hr(nh), hc(nh);
            model.eval_jac_g(n, nullptr, false, m, nj, jr.data(), jc.data(), nullptr);
            model.eval_h(n, nullptr, false, 1.3, m, nullptr, false, nh, hr.data(), hc.data(), nullptr);
            std::map<std::pair<int, int>, int> unique;
            for (int e = 0; e < nh; ++e) {
                require(hr[e] >= hc[e] && hr[e] < n && hc[e] >= 0, "invalid lower-triangle coordinate");
                require(unique.emplace(std::make_pair(hr[e], hc[e]), e).second, "duplicate Hessian coordinate");
            }
            for (int probe = 0; probe < 3; ++probe) {
                std::vector<double> lambda(m), h(nh), dense(n * n, 0.0);
                for (int row = 0; row < m; ++row) {
                    lambda[row] = probe == 2 ? 0.0 : 0.37 * std::sin(0.71 * (row + 1));
                    if (probe == 1 && row < 4 * static_cast<int>(topology_data.buses.size())) lambda[row] = 0.0;
                }
                require(model.eval_h(n, x.data(), true, 1.3, m, lambda.data(), true,
                        nh, nullptr, nullptr, h.data()), "Hessian evaluation failed");
                for (int e = 0; e < nh; ++e) {
                    dense[hr[e] * n + hc[e]] += h[e];
                    if (hr[e] != hc[e]) dense[hc[e] * n + hr[e]] += h[e];
                }
                const auto lagrangian_gradient = [&](const std::vector<double>& point) {
                    std::vector<double> gradient(n), jac(nj);
                    model.eval_grad_f(n, point.data(), true, gradient.data());
                    model.eval_jac_g(n, point.data(), true, m, nj, nullptr, nullptr, jac.data());
                    for (double& value : gradient) value *= 1.3;
                    for (int e = 0; e < nj; ++e) gradient[jc[e]] += lambda[jr[e]] * jac[e];
                    return gradient;
                };
                constexpr double step = 1e-6;
                for (int column = 0; column < n; ++column) {
                    auto lower = x, upper = x;
                    lower[column] -= step; upper[column] += step;
                    const auto left = lagrangian_gradient(lower), right = lagrangian_gradient(upper);
                    for (int row = 0; row < n; ++row) {
                        const double finite_difference = (right[row] - left[row]) / (2.0 * step);
                        const double exact = dense[row * n + column];
                        if (std::abs(exact - finite_difference) > 2e-5 + 2e-7 * std::abs(exact)) {
                            throw std::runtime_error("exact Hessian finite difference mismatch at " +
                                std::to_string(row) + "," + std::to_string(column) + ": " +
                                std::to_string(exact) + " vs " + std::to_string(finite_difference));
                        }
                    }
                }
            }
        }
    }
    for (bool epigraph : {false, true}) {
        auto nonsmooth = data;
        nonsmooth.generators[0].ncost = 3;
        nonsmooth.generators[0].cost = {0.0, 0.0, 1.0, 10.0, 2.0, 5.0};
        bool rejected = false;
        try {
            SparseAcEconomicNlp invalid(nonsmooth, commitment, start, 1e-5, epigraph, nullptr, true);
        } catch (const std::runtime_error&) { rejected = true; }
        require(rejected, "exact Hessian accepted unchecked nonsmooth source costs");
    }
}

nlohmann::json SparseAcEconomicResult::to_json(bool include_state) const {
    nlohmann::json result = {
        {"attempted", attempted},
        {"solver_initialized", solver_initialized},
        {"candidate_returned", candidate_returned},
        {"candidate_verified", candidate_verified},
        {"incumbent_verified", incumbent_verified},
        {"common_corrective_reference", common_corrective_reference},
        {"best_intermediate_found", best_intermediate_found},
        {"improved", improved},
        {"application_status", application_status},
        {"solver_return_status", solver_return_status},
        {"iterations", iterations},
        {"variable_count", variable_count},
        {"constraint_count", constraint_count},
        {"jacobian_nonzero_count", jacobian_nonzero_count},
        {"exact_hessian_enabled", exact_hessian_enabled},
        {"hessian_nonzero_count", hessian_nonzero_count},
        {"hessian_evaluations", hessian_evaluations},
        {"pwl_epigraph_enabled", pwl_epigraph_enabled},
        {"pwl_epigraph_curve_count", pwl_epigraph_curve_count},
        {"pwl_epigraph_row_count", pwl_epigraph_row_count},
        {"pwl_original_curve_count", pwl_original_curve_count},
        {"intermediate_callbacks", intermediate_callbacks},
        {"intermediate_iterates_retrieved", intermediate_iterates_retrieved},
        {"intermediate_verified_candidates", intermediate_verified_candidates},
        {"intermediate_capture_failures", intermediate_capture_failures},
        {"best_intermediate_iteration", best_intermediate_iteration},
        {"wall_seconds", wall_seconds},
        {"intermediate_callback_seconds", intermediate_callback_seconds},
        {"scaled_solver_objective", scaled_solver_objective},
        {"incumbent_objective", incumbent_objective},
        {"candidate_objective", candidate_objective},
        {"initial_constraint_violation", initial_constraint_violation},
        {"candidate_constraint_violation", candidate_constraint_violation},
        {"best_intermediate_objective", best_intermediate_objective},
        {"best_intermediate_max_residual", best_intermediate_max_residual},
        {"status", status},
        {"selected_source", selected_source},
        {"intermediate_capture_error", intermediate_capture_error},
        {"selected", solve_result_to_json(selected, include_state)},
        {"selected_validation", selected_validation.to_json()},
    };
    return result;
}

static SparseAcEconomicResult solve_sparse_ac_economic_impl(
    const CaseData& data,
    const std::vector<int>& commitment,
    const SolveResult& incumbent,
    const SparseAcEconomicOptions& options,
    const AcState* original_corrective_base) {
    if (!std::isfinite(options.time_limit_seconds) ||
        options.time_limit_seconds <= 0.0 ||
        !std::isfinite(options.tolerance) || options.tolerance <= 0.0 ||
        !std::isfinite(options.acceptable_tolerance) ||
        options.acceptable_tolerance <= 0.0) {
        throw std::runtime_error("invalid sparse AC economic options");
    }
    const auto wall_start = std::chrono::steady_clock::now();
    SparseAcEconomicResult output;
    output.attempted = true;
    output.common_corrective_reference = original_corrective_base != nullptr;
    output.selected = incumbent;
    output.selected.status = 0;
    output.selected.objective = rebuild_economic_candidate(
        data, commitment, original_corrective_base, output.selected.state);
    output.selected_validation = validate_economic_candidate(
        data, commitment, original_corrective_base, output.selected.state);
    output.incumbent_verified = std::isfinite(output.selected.objective) &&
        output.selected_validation.max_residual <= options.acceptable_tolerance;
    output.incumbent_objective = output.selected.objective;

    auto* raw_problem = new SparseAcEconomicNlp(
        data, commitment, output.selected.state,
        options.acceptable_tolerance, options.pwl_epigraph, original_corrective_base,
        options.exact_hessian);
    Ipopt::SmartPtr<Ipopt::TNLP> problem = raw_problem;
    output.variable_count = raw_problem->variable_count();
    output.constraint_count = raw_problem->constraint_count();
    output.jacobian_nonzero_count = raw_problem->jacobian_nonzero_count();
    output.exact_hessian_enabled = options.exact_hessian;
    output.hessian_nonzero_count = raw_problem->hessian_nonzero_count();
    output.pwl_epigraph_enabled = options.pwl_epigraph;
    output.pwl_epigraph_curve_count = raw_problem->pwl_epigraph_curve_count();
    output.pwl_epigraph_row_count = raw_problem->pwl_epigraph_row_count();
    output.pwl_original_curve_count = raw_problem->pwl_original_curve_count();
    output.initial_constraint_violation =
        raw_problem->initial_constraint_violation();

    auto application = IpoptApplicationFactory();
    application->Options()->SetIntegerValue(
        "print_level", options.print_level);
    application->Options()->SetStringValue(
        "linear_solver", "mumps");
    application->Options()->SetStringValue(
        "hessian_approximation", options.exact_hessian ? "exact" : "limited-memory");
    application->Options()->SetStringValue(
        "mu_strategy", "adaptive");
    application->Options()->SetStringValue(
        "nlp_scaling_method", "gradient-based");
    application->Options()->SetNumericValue("tol", options.tolerance);
    application->Options()->SetNumericValue(
        "acceptable_tol", options.acceptable_tolerance);
    application->Options()->SetIntegerValue("acceptable_iter", 3);
    application->Options()->SetNumericValue(
        "max_wall_time", options.time_limit_seconds);
    application->Options()->SetIntegerValue("max_iter", 500);
    application->Options()->SetStringValue("sb", "yes");
    const auto initialize_status = application->Initialize();
    output.solver_initialized = initialize_status == Ipopt::Solve_Succeeded;
    if (output.solver_initialized) {
        const auto status = application->OptimizeTNLP(problem);
        output.application_status = static_cast<int>(status);
        output.status = application_status_string(status);
        if (IsValid(application->Statistics())) {
            output.iterations = static_cast<int>(
                application->Statistics()->IterationCount());
        }
    } else {
        output.application_status = static_cast<int>(initialize_status);
        output.status = "Ipopt initialization failed";
    }
    output.solver_return_status = raw_problem->solver_return_status();
    output.hessian_evaluations = raw_problem->hessian_evaluations();
    output.scaled_solver_objective = raw_problem->scaled_objective();
    output.candidate_constraint_violation =
        raw_problem->final_constraint_violation();
    output.intermediate_callbacks = raw_problem->intermediate_callbacks();
    output.intermediate_iterates_retrieved =
        raw_problem->intermediate_iterates_retrieved();
    output.intermediate_verified_candidates =
        raw_problem->intermediate_verified_candidates();
    output.intermediate_capture_failures =
        raw_problem->intermediate_capture_failures();
    output.intermediate_callback_seconds =
        raw_problem->intermediate_callback_seconds();
    output.intermediate_capture_error =
        raw_problem->intermediate_capture_error();
    output.best_intermediate_found =
        raw_problem->best_intermediate_found();
    if (output.best_intermediate_found) {
        output.best_intermediate_iteration =
            raw_problem->best_intermediate_iteration();
        output.best_intermediate_objective =
            raw_problem->best_intermediate().objective;
        output.best_intermediate_max_residual =
            raw_problem->best_intermediate_validation().max_residual;
        if (output.selected_validation.max_residual > options.acceptable_tolerance ||
            raw_problem->best_intermediate().objective >
                output.selected.objective + 1e-9) {
            output.improved = true;
            output.selected = raw_problem->best_intermediate();
            output.selected_validation =
                raw_problem->best_intermediate_validation();
            output.selected_source = "verified_intermediate_iterate";
        }
    }
    output.candidate_returned = raw_problem->has_final_x();
    if (output.candidate_returned) {
        SolveResult candidate;
        candidate.status = 0;
        candidate.iterations = output.iterations;
        candidate.state = raw_problem->final_state();
        candidate.objective = rebuild_economic_candidate(
            data, commitment, original_corrective_base, candidate.state);
        const auto validation = validate_economic_candidate(
            data, commitment, original_corrective_base, candidate.state);
        output.candidate_objective = candidate.objective;
        output.candidate_verified = std::isfinite(candidate.objective) &&
            validation.max_residual <= options.acceptable_tolerance;
        if (output.candidate_verified &&
            (output.selected_validation.max_residual > options.acceptable_tolerance ||
             candidate.objective > output.selected.objective + 1e-9)) {
            output.improved = true;
            output.selected = std::move(candidate);
            output.selected_validation = validation;
            output.selected_source = "verified_final_iterate";
        }
    }
    output.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_start).count();
    output.selected.wall_seconds = output.wall_seconds;
    return output;
}

SparseAcEconomicResult solve_sparse_fixed_commitment_ac_economic(
    const CaseData& data, const std::vector<int>& commitment,
    const SolveResult& incumbent, const SparseAcEconomicOptions& options) {
    return solve_sparse_ac_economic_impl(data, commitment, incumbent, options, nullptr);
}

SparseAcEconomicResult solve_sparse_common_corrective_reference(
    const CaseData& data, const std::vector<int>& commitment,
    const SolveResult& original_base, const SparseAcEconomicOptions& options) {
    const auto validation = validate_state(
        data, ModelMode::BaseSoft, original_base.state, commitment);
    if (validation.max_residual > options.acceptable_tolerance) {
        throw std::runtime_error("common corrective reference requires a verified base");
    }
    return solve_sparse_ac_economic_impl(
        data, commitment, original_base, options, &original_base.state);
}

}  // namespace gravityx
