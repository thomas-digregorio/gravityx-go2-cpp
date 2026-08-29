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
        (b * tr - g * ti) / tm2,
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
        double verification_tolerance)
        : data_(data),
          commitment_(std::move(commitment)),
          start_state_(start),
          verification_tolerance_(verification_tolerance),
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
        build_variables();
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
            const double source_delta =
                data_.buses[item.from].va_start -
                data_.buses[item.to].va_start;
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
        constraint_count_ = reference_row_offset_ +
            static_cast<int>(reference_buses_.size());
        build_jacobian_structure();
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
        nnz_h_lag = 0;
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
            if (commitment_[i] == 0) {
                continue;
            }
            gradient[pg_offset_ + i] = data_.delta *
                pwl_value_slope(
                    generator_points_[i], x[pg_offset_ + i]).second /
                objective_scale_;
        }
        for (int i = 0; i < nd_; ++i) {
            const auto& load = data_.loads[i];
            gradient[demand_offset_ + i] = -data_.delta * load.pd_nominal *
                pwl_value_slope(
                    load_points_[i],
                    load.pd_nominal * x[demand_offset_ + i]).second /
                objective_scale_;
        }
        const double p_penalty =
            data_.delta * data_.p_delta_cost_approx / objective_scale_;
        const double q_penalty =
            data_.delta * data_.q_delta_cost_approx / objective_scale_;
        const double sm_penalty =
            data_.delta * data_.sm_cost_approx / objective_scale_;
        std::fill(
            gradient + p_delta_offset_,
            gradient + p_delta_offset_ + nb_, p_penalty);
        std::fill(
            gradient + q_delta_offset_,
            gradient + q_delta_offset_ + nb_, q_penalty);
        std::fill(
            gradient + sm_offset_,
            gradient + sm_offset_ + nl_, sm_penalty);
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
        Ipopt::Index,
        const Ipopt::Number*,
        bool,
        Ipopt::Number,
        Ipopt::Index,
        const Ipopt::Number*,
        bool,
        Ipopt::Index,
        Ipopt::Index*,
        Ipopt::Index*,
        Ipopt::Number*) override {
        return false;
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
                candidate.objective = rebuild_base_state_derived_fields(
                    data_, commitment_, candidate.state);
                const auto validation = validate_state(
                    data_, ModelMode::BaseSoft,
                    candidate.state, commitment_);
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
    int jacobian_nonzero_count() const {
        return static_cast<int>(jacobian_rows_.size());
    }
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
        state.sm_slack.assign(x.begin() + sm_offset_, x.end());
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
            const auto bounds = base_pg_bounds(
                data_.generators[i], commitment_[i], data_.delta_r);
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
            const auto bounds = base_load_bounds(data_.loads[i], data_.delta_r);
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

    double objective_value(const double* x) const {
        double objective = 0.0;
        for (int i = 0; i < ng_; ++i) {
            if (commitment_[i] == 0) {
                continue;
            }
            objective += data_.delta *
                pwl_value_slope(
                    generator_points_[i], x[pg_offset_ + i]).first;
            objective += data_.delta * data_.generators[i].oncost;
        }
        for (int i = 0; i < nd_; ++i) {
            objective -= data_.delta * pwl_value_slope(
                load_points_[i],
                data_.loads[i].pd_nominal * x[demand_offset_ + i]).first;
        }
        objective += data_.delta * data_.p_delta_cost_approx *
            std::accumulate(
                x + p_delta_offset_, x + p_delta_offset_ + nb_, 0.0);
        objective += data_.delta * data_.q_delta_cost_approx *
            std::accumulate(
                x + q_delta_offset_, x + q_delta_offset_ + nb_, 0.0);
        objective += data_.delta * data_.sm_cost_approx *
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
            const double rating_squared = branch.rate_a * branch.rate_a;
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
        const double rating_squared = branch.rate_a * branch.rate_a;
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
        if (next != jacobian_rows_.size()) {
            throw std::runtime_error("sparse AC Jacobian fill count mismatch");
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
        return violation;
    }

    const CaseData& data_;
    std::vector<int> commitment_;
    AcState start_state_;
    double verification_tolerance_{};
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
    int constraint_count_{};
    double objective_scale_{1.0};
    std::vector<double> x_lower_;
    std::vector<double> x_upper_;
    std::vector<double> start_x_;
    std::vector<std::vector<PwlPoint>> generator_points_;
    std::vector<std::vector<PwlPoint>> load_points_;
    YRows ybus_;
    std::vector<BranchCoefficients> coefficients_;
    std::vector<int> active_branches_;
    std::vector<int> angle_branches_;
    std::vector<int> reference_buses_;
    std::vector<Ipopt::Index> jacobian_rows_;
    std::vector<Ipopt::Index> jacobian_columns_;
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

nlohmann::json SparseAcEconomicResult::to_json(bool include_state) const {
    nlohmann::json result = {
        {"attempted", attempted},
        {"solver_initialized", solver_initialized},
        {"candidate_returned", candidate_returned},
        {"candidate_verified", candidate_verified},
        {"best_intermediate_found", best_intermediate_found},
        {"improved", improved},
        {"preserve_bound_active_start", preserve_bound_active_start},
        {"application_status", application_status},
        {"solver_return_status", solver_return_status},
        {"iterations", iterations},
        {"variable_count", variable_count},
        {"constraint_count", constraint_count},
        {"jacobian_nonzero_count", jacobian_nonzero_count},
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

SparseAcEconomicResult solve_sparse_fixed_commitment_ac_economic(
    const CaseData& data,
    const std::vector<int>& commitment,
    const SolveResult& incumbent,
    const SparseAcEconomicOptions& options) {
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
    output.preserve_bound_active_start =
        options.preserve_bound_active_start;
    output.selected = incumbent;
    output.selected.status = 0;
    output.selected.objective = rebuild_base_state_derived_fields(
        data, commitment, output.selected.state);
    output.selected_validation = validate_state(
        data, ModelMode::BaseSoft, output.selected.state, commitment);
    output.incumbent_objective = output.selected.objective;

    auto* raw_problem = new SparseAcEconomicNlp(
        data, commitment, output.selected.state,
        options.acceptable_tolerance);
    Ipopt::SmartPtr<Ipopt::TNLP> problem = raw_problem;
    output.variable_count = raw_problem->variable_count();
    output.constraint_count = raw_problem->constraint_count();
    output.jacobian_nonzero_count = raw_problem->jacobian_nonzero_count();
    output.initial_constraint_violation =
        raw_problem->initial_constraint_violation();

    auto application = IpoptApplicationFactory();
    application->Options()->SetIntegerValue(
        "print_level", options.print_level);
    application->Options()->SetStringValue(
        "linear_solver", "mumps");
    application->Options()->SetStringValue(
        "hessian_approximation", "limited-memory");
    application->Options()->SetStringValue(
        "mu_strategy", "adaptive");
    application->Options()->SetStringValue(
        "nlp_scaling_method", "gradient-based");
    if (options.preserve_bound_active_start) {
        // Keep the verified incumbent essentially where it was supplied.
        // These are initialization controls only: source bounds, nonlinear
        // equations, tolerances, and the independent acceptance gate remain
        // unchanged.
        constexpr double kInitialInterior = 1e-9;
        application->Options()->SetNumericValue(
            "bound_push", kInitialInterior);
        application->Options()->SetNumericValue(
            "bound_frac", kInitialInterior);
        application->Options()->SetNumericValue(
            "slack_bound_push", kInitialInterior);
        application->Options()->SetNumericValue(
            "slack_bound_frac", kInitialInterior);
        application->Options()->SetNumericValue(
            "bound_relax_factor", 0.0);
        application->Options()->SetStringValue(
            "honor_original_bounds", "yes");
        // Adaptive free-mode can increase mu by several orders of magnitude
        // after the first step.  On a high-quality bound-active incumbent
        // that recreates the same destructive central-path excursion the
        // small bound push is meant to avoid.  Monotone mode keeps the
        // deliberately small initial barrier and decreases it thereafter.
        application->Options()->SetStringValue(
            "mu_strategy", "monotone");
        application->Options()->SetNumericValue("mu_init", 1e-6);
    }
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
        if (raw_problem->best_intermediate().objective >
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
        candidate.objective = rebuild_base_state_derived_fields(
            data, commitment, candidate.state);
        const auto validation = validate_state(
            data, ModelMode::BaseSoft, candidate.state, commitment);
        output.candidate_objective = candidate.objective;
        output.candidate_verified = validation.max_residual <=
            options.acceptable_tolerance;
        if (output.candidate_verified &&
            candidate.objective > output.selected.objective + 1e-9) {
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

}  // namespace gravityx
