#include "gravityx/ac_model.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace gravityx {
namespace {

using gravity::Constraint;
using gravity::Model;
using gravity::R;
using gravity::func_;
using gravity::param;
using gravity::var;

constexpr double kBoundInfinity = 1e19;

void add_eq(Model& model, const std::string& name, func_ expression) {
    Constraint constraint(name);
    constraint = expression;
    model.add_constraint(constraint == 0.0);
}

void add_le(Model& model, const std::string& name, func_ expression) {
    Constraint constraint(name);
    constraint = expression;
    model.add_constraint(constraint <= 0.0);
}

void add_ge(Model& model, const std::string& name, func_ expression) {
    Constraint constraint(name);
    constraint = expression;
    model.add_constraint(constraint >= 0.0);
}

std::unique_ptr<var<double>> add_bounded_variable(
    Model& model,
    const std::string& name,
    const std::vector<double>& lower,
    const std::vector<double>& upper) {
    if (lower.size() != upper.size() || lower.empty()) {
        throw std::runtime_error("invalid bounds for variable " + name);
    }
    param<double> lb(name + "_lb");
    param<double> ub(name + "_ub");
    for (std::size_t i = 0; i < lower.size(); ++i) {
        lb.set_val(i, lower[i]);
        ub.set_val(i, upper[i]);
    }
    auto variable = std::make_unique<var<double>>(name, lb, ub);
    model.add_var(variable->in(R(lower.size())));
    // Gravity ACOPF revision f5af33e resizes the cached bound vectors in the
    // parameter-bounded var constructor but does not populate them.  Ipopt's
    // adapter reads those caches directly through get_lb/get_ub, so fill them
    // explicitly while preserving the source bounds exactly.
    for (std::size_t i = 0; i < lower.size(); ++i) {
        variable->_lb->_val->at(i) = lower[i];
        variable->_ub->_val->at(i) = upper[i];
    }
    variable->_lb->_evaluated = true;
    variable->_ub->_evaluated = true;
    return variable;
}

std::unique_ptr<var<double>> add_unbounded_variable(
    Model& model,
    const std::string& name,
    std::size_t count) {
    if (count == 0) {
        throw std::runtime_error("zero-length variable " + name);
    }
    auto variable = std::make_unique<var<double>>(name);
    model.add_var(variable->in(R(count)));
    return variable;
}

std::pair<double, double> base_pg_bounds(const Generator& gen, int status, double delta_r) {
    if (status == 0) {
        return {0.0, 0.0};
    }
    const double previous = gen.status_prev == 0 ? gen.pmin : gen.pg_prev;
    const double lower = std::max(gen.pmin, previous - delta_r * gen.prdmax);
    const double upper = std::min(gen.pmax, previous + delta_r * gen.prumax);
    if (lower > upper + 1e-12) {
        throw std::runtime_error("empty base ramp interval for generator " + gen.source_key);
    }
    return {lower, upper};
}

std::pair<double, double> base_load_bounds(const Load& load, double delta_r) {
    if (std::abs(load.pd_nominal) <= 1e-12) {
        return {load.tmin, load.tmax};
    }
    const double lower = std::max(load.tmin, (load.pd_prev - load.prdmax * delta_r) / load.pd_nominal);
    const double upper = std::min(load.tmax, (load.pd_prev + load.prumax * delta_r) / load.pd_nominal);
    if (lower > upper + 1e-12) {
        throw std::runtime_error("empty base ramp interval for load " + load.source_key);
    }
    return {lower, upper};
}

double clamp_to(double value, double lower, double upper) {
    return std::min(upper, std::max(lower, value));
}

std::vector<double> lambda_weights(const std::vector<PwlPoint>& points, double value) {
    std::vector<double> weights(points.size(), 0.0);
    if (value <= points.front().mw) {
        weights.front() = 1.0;
        return weights;
    }
    if (value >= points.back().mw) {
        weights.back() = 1.0;
        return weights;
    }
    for (std::size_t i = 0; i + 1 < points.size(); ++i) {
        if (value <= points[i + 1].mw) {
            const double width = points[i + 1].mw - points[i].mw;
            if (std::abs(width) <= 1e-14) {
                weights[i] = 1.0;
            } else {
                const double right = (value - points[i].mw) / width;
                weights[i] = 1.0 - right;
                weights[i + 1] = right;
            }
            return weights;
        }
    }
    weights.back() = 1.0;
    return weights;
}

}  // namespace

AcModel::AcModel(
    const CaseData& data,
    ModelMode mode,
    std::vector<int> fixed_status,
    std::optional<ContingencyContext> contingency)
    : data_(data),
      mode_(mode),
      fixed_status_(std::move(fixed_status)),
      contingency_(std::move(contingency)),
      model_(mode == ModelMode::BaseSoft ? "GO2 AC OPF soft" :
             mode == ModelMode::UnitCommitmentRelaxation ? "GO2 AC UC relaxation" :
             "GO2 corrective AC OPF soft") {
    if (fixed_status_.empty()) {
        fixed_status_.reserve(data_.generators.size());
        for (const auto& gen : data_.generators) {
            fixed_status_.push_back(gen.status_prev);
        }
    }
    if (fixed_status_.size() != data_.generators.size()) {
        throw std::runtime_error("fixed commitment vector has the wrong length");
    }
    if ((mode_ == ModelMode::ContingencySoft) != contingency_.has_value()) {
        throw std::runtime_error("contingency mode and context must be supplied together");
    }
    if (contingency_) {
        const auto& state = contingency_->base_state;
        if (state.vm.size() != data_.buses.size() || state.va.size() != data_.buses.size() ||
            state.pg.size() != data_.generators.size() || state.qg.size() != data_.generators.size() ||
            state.demand_factor.size() != data_.loads.size()) {
            throw std::runtime_error("contingency base state dimensions are invalid");
        }
    }
    build_variables();
    build_constraints_and_objective();
    initialize_source_point();
}

void AcModel::build_variables() {
    const std::size_t nb = data_.buses.size();
    const std::size_t ng = data_.generators.size();
    const std::size_t nd = data_.loads.size();
    const std::size_t nl = data_.branches.size();

    std::vector<double> lower(nb), upper(nb);
    for (std::size_t i = 0; i < nb; ++i) {
        lower[i] = data_.buses[i].vmin;
        upper[i] = data_.buses[i].vmax;
    }
    vm_ = add_bounded_variable(model_, "vm", lower, upper);
    va_ = add_unbounded_variable(model_, "va", nb);

    std::vector<double> pg_lower(ng), pg_upper(ng), qg_lower(ng), qg_upper(ng);
    gen_points_.resize(ng);
    gen_lambda_offset_.assign(ng, -1);
    int gen_lambda_count = 0;
    for (std::size_t i = 0; i < ng; ++i) {
        const auto& gen = data_.generators[i];
        const bool outaged = contingency_ && contingency_->outaged_generator == static_cast<int>(i);
        const bool active = !outaged &&
            (mode_ == ModelMode::UnitCommitmentRelaxation || fixed_status_[i] == 1);
        if (mode_ == ModelMode::UnitCommitmentRelaxation) {
            pg_lower[i] = std::min(0.0, gen.pmin);
            pg_upper[i] = std::max(0.0, gen.pmax);
            qg_lower[i] = std::min(0.0, gen.qmin);
            qg_upper[i] = std::max(0.0, gen.qmax);
        } else if (active) {
            const auto bounds = mode_ == ModelMode::ContingencySoft
                ? std::pair<double, double>{
                    std::max(gen.pmin, contingency_->base_state.pg[i] - data_.delta_r_ctg * gen.prdmaxctg),
                    std::min(gen.pmax, contingency_->base_state.pg[i] + data_.delta_r_ctg * gen.prumaxctg)}
                : base_pg_bounds(gen, 1, data_.delta_r);
            if (bounds.first > bounds.second + 1e-12) {
                throw std::runtime_error("empty generator ramp interval in " + gen.source_key);
            }
            pg_lower[i] = bounds.first;
            pg_upper[i] = bounds.second;
            qg_lower[i] = gen.qmin;
            qg_upper[i] = gen.qmax;
        } else {
            pg_lower[i] = pg_upper[i] = 0.0;
            qg_lower[i] = qg_upper[i] = 0.0;
        }
        if (active) {
            gen_points_[i] = active_pwl_points(gen.cost, gen.ncost, pg_lower[i], pg_upper[i]);
            gen_lambda_offset_[i] = gen_lambda_count;
            gen_lambda_count += static_cast<int>(gen_points_[i].size());
        }
    }
    pg_ = add_bounded_variable(model_, "pg", pg_lower, pg_upper);
    qg_ = add_bounded_variable(model_, "qg", qg_lower, qg_upper);

    std::vector<double> demand_lower(nd), demand_upper(nd);
    load_points_.resize(nd);
    load_lambda_offset_.resize(nd);
    int load_lambda_count = 0;
    for (std::size_t i = 0; i < nd; ++i) {
        const auto& load = data_.loads[i];
        if (mode_ == ModelMode::BaseSoft) {
            const auto bounds = base_load_bounds(load, data_.delta_r);
            demand_lower[i] = bounds.first;
            demand_upper[i] = bounds.second;
        } else if (mode_ == ModelMode::ContingencySoft) {
            const double previous = load.pd_nominal * contingency_->base_state.demand_factor[i];
            if (std::abs(load.pd_nominal) <= 1e-12) {
                demand_lower[i] = load.tmin;
                demand_upper[i] = load.tmax;
            } else {
                demand_lower[i] = std::max(
                    load.tmin, (previous - load.prdmaxctg * data_.delta_r_ctg) / load.pd_nominal);
                demand_upper[i] = std::min(
                    load.tmax, (previous + load.prumaxctg * data_.delta_r_ctg) / load.pd_nominal);
            }
        } else {
            demand_lower[i] = load.tmin;
            demand_upper[i] = load.tmax;
        }
        load_points_[i] = active_pwl_points(load.cost, load.ncost, load.pd_min, load.pd_max);
        load_lambda_offset_[i] = load_lambda_count;
        load_lambda_count += static_cast<int>(load_points_[i].size());
    }
    demand_ = add_bounded_variable(model_, "demand_factor", demand_lower, demand_upper);

    std::vector<double> flow_lower(nl), flow_upper(nl), slack_lower(nl, 0.0), slack_upper(nl, data_.sm_vio_limit);
    for (std::size_t i = 0; i < nl; ++i) {
        const bool outaged = contingency_ && contingency_->outaged_branch == static_cast<int>(i);
        flow_lower[i] = outaged ? 0.0 : -data_.branches[i].rate_a;
        flow_upper[i] = outaged ? 0.0 : data_.branches[i].rate_a;
        if (outaged) {
            slack_upper[i] = 0.0;
        }
    }
    pf_ = add_bounded_variable(model_, "pf", flow_lower, flow_upper);
    qf_ = add_bounded_variable(model_, "qf", flow_lower, flow_upper);
    pt_ = add_bounded_variable(model_, "pt", flow_lower, flow_upper);
    qt_ = add_bounded_variable(model_, "qt", flow_lower, flow_upper);
    sm_slack_ = add_bounded_variable(model_, "sm_slack", slack_lower, slack_upper);

    if (mode_ != ModelMode::UnitCommitmentRelaxation) {
        p_delta_ = add_bounded_variable(model_, "p_delta", std::vector<double>(nb, 0.0), std::vector<double>(nb, 0.5));
        q_delta_ = add_bounded_variable(model_, "q_delta", std::vector<double>(nb, 0.0), std::vector<double>(nb, 0.5));
    } else {
        std::vector<double> z_lower(ng, 0.0), z_upper(ng, 1.0);
        for (std::size_t i = 0; i < ng; ++i) {
            const auto& gen = data_.generators[i];
            if (gen.status_prev == 0 && gen.suqual == 0) {
                z_upper[i] = 0.0;
            }
            if (gen.status_prev == 1 && gen.sdqual == 0) {
                z_lower[i] = 1.0;
            }
        }
        commitment_ = add_bounded_variable(model_, "commitment", z_lower, z_upper);
        startup_ = add_bounded_variable(model_, "startup", std::vector<double>(ng, 0.0), std::vector<double>(ng, 1.0));
        shutdown_ = add_bounded_variable(model_, "shutdown", std::vector<double>(ng, 0.0), std::vector<double>(ng, 1.0));
    }

    if (gen_lambda_count <= 0 || load_lambda_count <= 0) {
        throw std::runtime_error("case did not produce piecewise-linear variables");
    }
    gen_lambda_ = add_bounded_variable(model_, "gen_cost_lambda",
        std::vector<double>(gen_lambda_count, 0.0), std::vector<double>(gen_lambda_count, 1.0));
    load_lambda_ = add_bounded_variable(model_, "load_value_lambda",
        std::vector<double>(load_lambda_count, 0.0), std::vector<double>(load_lambda_count, 1.0));
}

void AcModel::build_constraints_and_objective() {
    func_ objective;
    const double interval_delta = mode_ == ModelMode::ContingencySoft
        ? data_.delta_ctg
        : data_.delta;

    for (std::size_t i = 0; i < data_.generators.size(); ++i) {
        if (gen_lambda_offset_[i] < 0) {
            continue;
        }
        func_ sum_lambda;
        func_ power_expression;
        const int offset = gen_lambda_offset_[i];
        for (int j = 0; j < static_cast<int>(gen_points_[i].size()); ++j) {
            sum_lambda += (*gen_lambda_)(offset + j);
            power_expression += gen_points_[i][j].mw * (*gen_lambda_)(offset + j);
            objective -= interval_delta * gen_points_[i][j].cost * (*gen_lambda_)(offset + j);
        }
        add_eq(model_, "gen_lambda_sum_" + std::to_string(i), sum_lambda - 1.0);
        add_eq(model_, "gen_pwl_power_" + std::to_string(i), power_expression - (*pg_)(i));
    }

    for (std::size_t i = 0; i < data_.loads.size(); ++i) {
        func_ sum_lambda;
        func_ power_expression;
        const int offset = load_lambda_offset_[i];
        for (int j = 0; j < static_cast<int>(load_points_[i].size()); ++j) {
            sum_lambda += (*load_lambda_)(offset + j);
            power_expression += load_points_[i][j].mw * (*load_lambda_)(offset + j);
            objective += interval_delta * load_points_[i][j].cost * (*load_lambda_)(offset + j);
        }
        add_eq(model_, "load_lambda_sum_" + std::to_string(i), sum_lambda - 1.0);
        add_eq(model_, "load_pwl_power_" + std::to_string(i),
               power_expression - data_.loads[i].pd_nominal * (*demand_)(i));
    }

    if (mode_ == ModelMode::UnitCommitmentRelaxation) {
        for (std::size_t i = 0; i < data_.generators.size(); ++i) {
            const auto& gen = data_.generators[i];
            add_le(model_, "gen_pmax_" + std::to_string(i), (*pg_)(i) - gen.pmax * (*commitment_)(i));
            add_ge(model_, "gen_pmin_" + std::to_string(i), (*pg_)(i) - gen.pmin * (*commitment_)(i));
            add_le(model_, "gen_qmax_" + std::to_string(i), (*qg_)(i) - gen.qmax * (*commitment_)(i));
            add_ge(model_, "gen_qmin_" + std::to_string(i), (*qg_)(i) - gen.qmin * (*commitment_)(i));

            const double previous = gen.status_prev == 0 ? gen.pmin : gen.pg_prev;
            add_le(model_, "gen_ramp_up_" + std::to_string(i),
                   (*pg_)(i) - previous * (*commitment_)(i) - data_.delta_r * gen.prumax);
            add_ge(model_, "gen_ramp_down_" + std::to_string(i),
                   (*pg_)(i) - previous * (*commitment_)(i) + data_.delta_r * gen.prdmax);
            add_eq(model_, "gen_transition_" + std::to_string(i),
                   (*commitment_)(i) - static_cast<double>(gen.status_prev) - (*startup_)(i) + (*shutdown_)(i));
            add_le(model_, "gen_switch_exclusive_" + std::to_string(i), (*startup_)(i) + (*shutdown_)(i) - 1.0);

            objective -= data_.delta * gen.oncost * (*commitment_)(i);
            objective -= gen.sucost * (*startup_)(i);
            objective -= gen.sdcost * (*shutdown_)(i);
        }
        for (std::size_t i = 0; i < data_.loads.size(); ++i) {
            const auto& load = data_.loads[i];
            add_le(model_, "load_ramp_up_" + std::to_string(i),
                   load.pd_nominal * (*demand_)(i) - load.pd_prev - load.prumax * data_.delta_r);
            add_ge(model_, "load_ramp_down_" + std::to_string(i),
                   load.pd_nominal * (*demand_)(i) - load.pd_prev + load.prdmax * data_.delta_r);
        }
    } else {
        for (std::size_t i = 0; i < data_.generators.size(); ++i) {
            const bool outaged = contingency_ && contingency_->outaged_generator == static_cast<int>(i);
            if (!outaged && fixed_status_[i] == 1) {
                objective -= interval_delta * data_.generators[i].oncost;
            }
        }
    }

    for (std::size_t i = 0; i < data_.buses.size(); ++i) {
        if (data_.buses[i].type == 3) {
            add_eq(model_, "reference_angle_" + std::to_string(i), (*va_)(i));
        }
    }

    for (std::size_t i = 0; i < data_.buses.size(); ++i) {
        const auto& bus = data_.buses[i];
        func_ p_balance;
        func_ q_balance;
        for (int branch : bus.branches_from) {
            p_balance += (*pf_)(branch);
            q_balance += (*qf_)(branch);
        }
        for (int branch : bus.branches_to) {
            p_balance += (*pt_)(branch);
            q_balance += (*qt_)(branch);
        }
        for (int generator : bus.generators) {
            p_balance -= (*pg_)(generator);
            q_balance -= (*qg_)(generator);
        }
        for (int load : bus.loads) {
            p_balance += data_.loads[load].pd_nominal * (*demand_)(load);
            q_balance += data_.loads[load].qd_nominal * (*demand_)(load);
        }
        double gs = 0.0;
        double bs = 0.0;
        for (int shunt : bus.shunts) {
            gs += data_.shunts[shunt].gs;
            bs += data_.shunts[shunt].bs;
        }
        p_balance += gs * gravity::power((*vm_)(i), 2);
        q_balance += (-bs) * gravity::power((*vm_)(i), 2);

        if (mode_ != ModelMode::UnitCommitmentRelaxation) {
            add_le(model_, "p_balance_pos_" + std::to_string(i), p_balance - (*p_delta_)(i));
            add_le(model_, "p_balance_neg_" + std::to_string(i), -1.0 * p_balance - (*p_delta_)(i));
            add_le(model_, "q_balance_pos_" + std::to_string(i), q_balance - (*q_delta_)(i));
            add_le(model_, "q_balance_neg_" + std::to_string(i), -1.0 * q_balance - (*q_delta_)(i));
            objective -= interval_delta * data_.p_delta_cost_approx * (*p_delta_)(i);
            objective -= interval_delta * data_.q_delta_cost_approx * (*q_delta_)(i);
        } else {
            add_eq(model_, "p_balance_" + std::to_string(i), p_balance);
            add_eq(model_, "q_balance_" + std::to_string(i), q_balance);
        }
    }

    for (std::size_t i = 0; i < data_.branches.size(); ++i) {
        if (contingency_ && contingency_->outaged_branch == static_cast<int>(i)) {
            continue;
        }
        const auto& branch = data_.branches[i];
        const double denominator = branch.r * branch.r + branch.x * branch.x;
        const double g = denominator > 1e-20 ? branch.r / denominator : 0.0;
        const double b = denominator > 1e-20 ? -branch.x / denominator : 0.0;
        const double tm = branch.tap;
        if (std::abs(tm) <= 1e-12) {
            throw std::runtime_error("zero tap ratio on branch " + branch.source_key);
        }
        const double tr = tm * std::cos(branch.shift);
        const double ti = tm * std::sin(branch.shift);
        const double tm2 = tm * tm;
        const int f = branch.from;
        const int t = branch.to;

        const auto cross_cos_ft = (*vm_)(f) * (*vm_)(t) * gravity::cos((*va_)(f) - (*va_)(t));
        const auto cross_sin_ft = (*vm_)(f) * (*vm_)(t) * gravity::sin((*va_)(f) - (*va_)(t));
        const double from_g_self = branch.transformer ? g / tm2 + branch.g_fr : (g + branch.g_fr) / tm2;
        const double from_b_self = branch.transformer ? b / tm2 + branch.b_fr : (b + branch.b_fr) / tm2;

        func_ p_from = (*pf_)(i) + (-from_g_self) * gravity::power((*vm_)(f), 2);
        p_from += (-((-g * tr + b * ti) / tm2)) * cross_cos_ft;
        p_from += (-((-b * tr - g * ti) / tm2)) * cross_sin_ft;
        add_eq(model_, "ohms_pf_" + std::to_string(i), p_from);

        func_ q_from = (*qf_)(i) + from_b_self * gravity::power((*vm_)(f), 2);
        q_from += ((-b * tr - g * ti) / tm2) * cross_cos_ft;
        q_from += (-((-g * tr + b * ti) / tm2)) * cross_sin_ft;
        add_eq(model_, "ohms_qf_" + std::to_string(i), q_from);

        const auto cross_cos_tf = (*vm_)(t) * (*vm_)(f) * gravity::cos((*va_)(t) - (*va_)(f));
        const auto cross_sin_tf = (*vm_)(t) * (*vm_)(f) * gravity::sin((*va_)(t) - (*va_)(f));
        func_ p_to = (*pt_)(i) + (-(g + branch.g_to)) * gravity::power((*vm_)(t), 2);
        p_to += (-((-g * tr - b * ti) / tm2)) * cross_cos_tf;
        p_to += (-((-b * tr + g * ti) / tm2)) * cross_sin_tf;
        add_eq(model_, "ohms_pt_" + std::to_string(i), p_to);

        func_ q_to = (*qt_)(i) + (b + branch.b_to) * gravity::power((*vm_)(t), 2);
        q_to += ((-b * tr + g * ti) / tm2) * cross_cos_tf;
        q_to += (-((-g * tr - b * ti) / tm2)) * cross_sin_tf;
        add_eq(model_, "ohms_qt_" + std::to_string(i), q_to);

        const double start_delta = mode_ == ModelMode::ContingencySoft
            ? contingency_->base_state.va[f] - contingency_->base_state.va[t]
            : data_.buses[f].va_start - data_.buses[t].va_start;
        if (start_delta >= branch.angmin && start_delta <= branch.angmax) {
            add_le(model_, "angle_upper_" + std::to_string(i), (*va_)(f) - (*va_)(t) - branch.angmax);
            add_ge(model_, "angle_lower_" + std::to_string(i), (*va_)(f) - (*va_)(t) - branch.angmin);
        }

        func_ thermal_from = gravity::power((*pf_)(i), 2) + gravity::power((*qf_)(i), 2);
        func_ thermal_to = gravity::power((*pt_)(i), 2) + gravity::power((*qt_)(i), 2);
        if (branch.transformer) {
            thermal_from += (-branch.rate_a * branch.rate_a) * gravity::power(1.0 + (*sm_slack_)(i), 2);
            thermal_to += (-branch.rate_a * branch.rate_a) * gravity::power(1.0 + (*sm_slack_)(i), 2);
        } else {
            thermal_from += (-branch.rate_a * branch.rate_a) * gravity::power((*vm_)(f) + (*sm_slack_)(i), 2);
            thermal_to += (-branch.rate_a * branch.rate_a) * gravity::power((*vm_)(t) + (*sm_slack_)(i), 2);
        }
        add_le(model_, "thermal_from_" + std::to_string(i), thermal_from);
        add_le(model_, "thermal_to_" + std::to_string(i), thermal_to);
        objective -= interval_delta * data_.sm_cost_approx * (*sm_slack_)(i);
    }

    model_.max(objective);
}

void AcModel::initialize_source_point() {
    for (std::size_t i = 0; i < data_.buses.size(); ++i) {
        const double vm = mode_ == ModelMode::ContingencySoft
            ? contingency_->base_state.vm[i]
            : data_.buses[i].vm_start;
        const double va = mode_ == ModelMode::ContingencySoft
            ? contingency_->base_state.va[i]
            : data_.buses[i].va_start;
        vm_->initialize(i, clamp_to(vm, data_.buses[i].vmin, data_.buses[i].vmax));
        va_->initialize(i, va);
    }

    for (std::size_t i = 0; i < data_.generators.size(); ++i) {
        const auto& gen = data_.generators[i];
        double pg = mode_ == ModelMode::ContingencySoft
            ? contingency_->base_state.pg[i]
            : gen.pg_start;
        double qg = mode_ == ModelMode::ContingencySoft
            ? contingency_->base_state.qg[i]
            : gen.qg_start;
        if (mode_ == ModelMode::BaseSoft) {
            const auto bounds = base_pg_bounds(gen, fixed_status_[i], data_.delta_r);
            pg = clamp_to(pg, bounds.first, bounds.second);
            qg = fixed_status_[i] == 1 ? clamp_to(qg, gen.qmin, gen.qmax) : 0.0;
        } else if (mode_ == ModelMode::UnitCommitmentRelaxation) {
            const double z = (gen.status_prev == 0 && gen.suqual == 1) ? 1.0 : static_cast<double>(gen.status_prev);
            commitment_->initialize(i, z);
            startup_->initialize(i, std::max(0.0, z - gen.status_prev));
            shutdown_->initialize(i, std::max(0.0, gen.status_prev - z));
            if (z > 0.5 && gen.status_prev == 0) {
                pg = gen.pmin;
                qg = clamp_to(0.0, gen.qmin, gen.qmax);
            }
            pg = clamp_to(pg, std::min(0.0, gen.pmin), std::max(0.0, gen.pmax));
            qg = clamp_to(qg, std::min(0.0, gen.qmin), std::max(0.0, gen.qmax));
        } else {
            const bool outaged = contingency_->outaged_generator == static_cast<int>(i);
            if (outaged || fixed_status_[i] == 0) {
                pg = 0.0;
                qg = 0.0;
            } else {
                const double lower = std::max(
                    gen.pmin, contingency_->base_state.pg[i] - data_.delta_r_ctg * gen.prdmaxctg);
                const double upper = std::min(
                    gen.pmax, contingency_->base_state.pg[i] + data_.delta_r_ctg * gen.prumaxctg);
                pg = clamp_to(pg, lower, upper);
                qg = clamp_to(qg, gen.qmin, gen.qmax);
            }
        }
        pg_->initialize(i, pg);
        qg_->initialize(i, qg);

        if (gen_lambda_offset_[i] >= 0) {
            const auto weights = lambda_weights(gen_points_[i], pg);
            for (int j = 0; j < static_cast<int>(weights.size()); ++j) {
                gen_lambda_->initialize(gen_lambda_offset_[i] + j, weights[j]);
            }
        }
    }

    for (std::size_t i = 0; i < data_.loads.size(); ++i) {
        double z = mode_ == ModelMode::ContingencySoft
            ? contingency_->base_state.demand_factor[i]
            : data_.loads[i].z_start;
        if (mode_ == ModelMode::BaseSoft) {
            const auto bounds = base_load_bounds(data_.loads[i], data_.delta_r);
            z = clamp_to(z, bounds.first, bounds.second);
        } else if (mode_ == ModelMode::UnitCommitmentRelaxation) {
            z = clamp_to(z, data_.loads[i].tmin, data_.loads[i].tmax);
        } else {
            const auto& load = data_.loads[i];
            const double previous = load.pd_nominal * contingency_->base_state.demand_factor[i];
            const double lower = std::abs(load.pd_nominal) <= 1e-12
                ? load.tmin
                : std::max(load.tmin, (previous - load.prdmaxctg * data_.delta_r_ctg) / load.pd_nominal);
            const double upper = std::abs(load.pd_nominal) <= 1e-12
                ? load.tmax
                : std::min(load.tmax, (previous + load.prumaxctg * data_.delta_r_ctg) / load.pd_nominal);
            z = clamp_to(z, lower, upper);
        }
        demand_->initialize(i, z);
        const auto weights = lambda_weights(load_points_[i], data_.loads[i].pd_nominal * z);
        for (int j = 0; j < static_cast<int>(weights.size()); ++j) {
            load_lambda_->initialize(load_lambda_offset_[i] + j, weights[j]);
        }
    }

    for (std::size_t i = 0; i < data_.branches.size(); ++i) {
        const bool outaged = contingency_ && contingency_->outaged_branch == static_cast<int>(i);
        pf_->initialize(i, mode_ == ModelMode::ContingencySoft && !outaged ? contingency_->base_state.pf[i] : 0.0);
        qf_->initialize(i, mode_ == ModelMode::ContingencySoft && !outaged ? contingency_->base_state.qf[i] : 0.0);
        pt_->initialize(i, mode_ == ModelMode::ContingencySoft && !outaged ? contingency_->base_state.pt[i] : 0.0);
        qt_->initialize(i, mode_ == ModelMode::ContingencySoft && !outaged ? contingency_->base_state.qt[i] : 0.0);
        sm_slack_->initialize(i, 0.0);
    }
    if (mode_ != ModelMode::UnitCommitmentRelaxation) {
        p_delta_->initialize_all(0.0);
        q_delta_->initialize_all(0.0);
    }
}

SolveResult AcModel::solve(int print_level, double tolerance) {
    const auto start = std::chrono::steady_clock::now();
    gravity::solver nlp(model_, gravity::ipopt);
    const int status = nlp.run(print_level, false, tolerance, "mumps", "no");
    const auto finish = std::chrono::steady_clock::now();

    SolveResult result;
    result.status = status;
    result.objective = model_._obj_val;
    result.wall_seconds = std::chrono::duration<double>(finish - start).count();
    result.state = capture_state();
    return result;
}

AcState AcModel::capture_state() const {
    AcState state;
    const auto capture = [](const std::unique_ptr<var<double>>& variable, std::size_t count) {
        std::vector<double> values;
        if (!variable) {
            return values;
        }
        values.reserve(count);
        for (std::size_t i = 0; i < count; ++i) {
            values.push_back(variable->eval(i));
        }
        return values;
    };
    state.vm = capture(vm_, data_.buses.size());
    state.va = capture(va_, data_.buses.size());
    state.pg = capture(pg_, data_.generators.size());
    state.qg = capture(qg_, data_.generators.size());
    state.demand_factor = capture(demand_, data_.loads.size());
    state.pf = capture(pf_, data_.branches.size());
    state.qf = capture(qf_, data_.branches.size());
    state.pt = capture(pt_, data_.branches.size());
    state.qt = capture(qt_, data_.branches.size());
    state.sm_slack = capture(sm_slack_, data_.branches.size());
    state.p_delta = capture(p_delta_, data_.buses.size());
    state.q_delta = capture(q_delta_, data_.buses.size());
    state.commitment = capture(commitment_, data_.generators.size());
    state.startup = capture(startup_, data_.generators.size());
    state.shutdown = capture(shutdown_, data_.generators.size());
    state.gen_lambda = capture(gen_lambda_, gen_lambda_->get_nb_instances());
    state.load_lambda = capture(load_lambda_, load_lambda_->get_nb_instances());
    return state;
}

void AcModel::initialize_from(const AcState& state) {
    const auto initialize = [](std::unique_ptr<var<double>>& variable, const std::vector<double>& values) {
        if (!variable) {
            return;
        }
        for (std::size_t i = 0; i < values.size(); ++i) {
            variable->initialize(i, values[i]);
        }
    };
    initialize(vm_, state.vm);
    initialize(va_, state.va);
    initialize(demand_, state.demand_factor);
    initialize(pf_, state.pf);
    initialize(qf_, state.qf);
    initialize(pt_, state.pt);
    initialize(qt_, state.qt);
    initialize(sm_slack_, state.sm_slack);
    initialize(p_delta_, state.p_delta);
    initialize(q_delta_, state.q_delta);

    for (std::size_t i = 0; i < data_.generators.size() && i < state.pg.size(); ++i) {
        const auto& gen = data_.generators[i];
        double pg = state.pg[i];
        double qg = i < state.qg.size() ? state.qg[i] : 0.0;
        if (mode_ == ModelMode::BaseSoft) {
            const auto bounds = base_pg_bounds(gen, fixed_status_[i], data_.delta_r);
            pg = clamp_to(pg, bounds.first, bounds.second);
            qg = fixed_status_[i] == 1 ? clamp_to(qg, gen.qmin, gen.qmax) : 0.0;
        } else if (mode_ == ModelMode::UnitCommitmentRelaxation) {
            pg = clamp_to(pg, std::min(0.0, gen.pmin), std::max(0.0, gen.pmax));
            qg = clamp_to(qg, std::min(0.0, gen.qmin), std::max(0.0, gen.qmax));
        } else {
            const bool outaged = contingency_->outaged_generator == static_cast<int>(i);
            if (outaged || fixed_status_[i] == 0) {
                pg = 0.0;
                qg = 0.0;
            } else {
                const double lower = std::max(
                    gen.pmin, contingency_->base_state.pg[i] - data_.delta_r_ctg * gen.prdmaxctg);
                const double upper = std::min(
                    gen.pmax, contingency_->base_state.pg[i] + data_.delta_r_ctg * gen.prumaxctg);
                pg = clamp_to(pg, lower, upper);
                qg = clamp_to(qg, gen.qmin, gen.qmax);
            }
        }
        pg_->initialize(i, pg);
        qg_->initialize(i, qg);
        if (gen_lambda_offset_[i] >= 0) {
            const auto weights = lambda_weights(gen_points_[i], pg);
            for (int j = 0; j < static_cast<int>(weights.size()); ++j) {
                gen_lambda_->initialize(gen_lambda_offset_[i] + j, weights[j]);
            }
        }
    }
    for (std::size_t i = 0; i < data_.loads.size() && i < state.demand_factor.size(); ++i) {
        const double pd = data_.loads[i].pd_nominal * state.demand_factor[i];
        const auto weights = lambda_weights(load_points_[i], pd);
        for (int j = 0; j < static_cast<int>(weights.size()); ++j) {
            load_lambda_->initialize(load_lambda_offset_[i] + j, weights[j]);
        }
    }
    initialize(commitment_, state.commitment);
    initialize(startup_, state.startup);
    initialize(shutdown_, state.shutdown);
}

void AcModel::set_commitment_bound(int generator, int status) {
    if (!commitment_ || generator < 0 || generator >= static_cast<int>(data_.generators.size())) {
        throw std::runtime_error("invalid commitment bound update");
    }
    commitment_->_lb->_val->at(generator) = static_cast<double>(status);
    commitment_->_ub->_val->at(generator) = static_cast<double>(status);
    set_commitment_start(generator, status);
}

void AcModel::set_commitment_start(int generator, double value) {
    if (!commitment_ || generator < 0 || generator >= static_cast<int>(data_.generators.size())) {
        throw std::runtime_error("invalid commitment start update");
    }
    const auto& gen = data_.generators[generator];
    const double z = clamp_to(value, 0.0, 1.0);
    commitment_->initialize(generator, z);
    startup_->initialize(generator, std::max(0.0, z - gen.status_prev));
    shutdown_->initialize(generator, std::max(0.0, gen.status_prev - z));
}

}  // namespace gravityx
