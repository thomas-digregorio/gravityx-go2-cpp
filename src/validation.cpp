#include "gravityx/validation.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace gravityx {
namespace {

double positive_part(double value) {
    return std::max(0.0, value);
}

template <typename IdentityFactory>
void update_category(
    ValidationReport& report,
    double& category,
    double value,
    const char* name,
    bool capture_identity,
    IdentityFactory&& identity_factory) {
    value = std::abs(value);
    category = std::max(category, value);
    if (value > report.max_residual) {
        report.max_residual = value;
        report.worst_category = name;
        report.worst_identity = capture_identity
            ? identity_factory() : std::string();
    }
}

std::pair<double, double> base_pg_bounds(const Generator& gen, int status, double delta_r) {
    if (status == 0) {
        return {0.0, 0.0};
    }
    const double previous = gen.status_prev == 0 ? gen.pmin : gen.pg_prev;
    return {
        std::max(gen.pmin, previous - delta_r * gen.prdmax),
        std::min(gen.pmax, previous + delta_r * gen.prumax),
    };
}

std::pair<double, double> base_load_bounds(const Load& load, double delta_r) {
    if (std::abs(load.pd_nominal) <= 1e-12) {
        return {load.tmin, load.tmax};
    }
    return {
        std::max(load.tmin, (load.pd_prev - load.prdmax * delta_r) / load.pd_nominal),
        std::min(load.tmax, (load.pd_prev + load.prumax * delta_r) / load.pd_nominal),
    };
}

double bound_violation(double value, double lower, double upper) {
    return std::max(positive_part(lower - value), positive_part(value - upper));
}

void require_finite_vector(const std::vector<double>& values, const std::string& name) {
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (!std::isfinite(values[i])) {
            throw std::runtime_error(
                "state contains nonfinite " + name + " value at position " +
                std::to_string(i));
        }
    }
}

void require_finite_state(const AcState& state, const std::string& prefix) {
    require_finite_vector(state.vm, prefix + "vm");
    require_finite_vector(state.va, prefix + "va");
    require_finite_vector(state.pg, prefix + "pg");
    require_finite_vector(state.qg, prefix + "qg");
    require_finite_vector(state.demand_factor, prefix + "demand_factor");
    require_finite_vector(state.shunt_bs, prefix + "shunt_bs");
    require_finite_vector(state.pf, prefix + "pf");
    require_finite_vector(state.qf, prefix + "qf");
    require_finite_vector(state.pt, prefix + "pt");
    require_finite_vector(state.qt, prefix + "qt");
    require_finite_vector(state.sm_slack, prefix + "sm_slack");
    require_finite_vector(state.p_delta, prefix + "p_delta");
    require_finite_vector(state.q_delta, prefix + "q_delta");
    require_finite_vector(state.commitment, prefix + "commitment");
    require_finite_vector(state.startup, prefix + "startup");
    require_finite_vector(state.shutdown, prefix + "shutdown");
    require_finite_vector(state.gen_lambda, prefix + "gen_lambda");
    require_finite_vector(state.load_lambda, prefix + "load_lambda");
}

void validate_generator_pwl(
    ValidationReport& report,
    const Generator& gen,
    const AcState& state,
    std::size_t generator_index,
    double pg_lower,
    double pg_upper,
    int& lambda_offset,
    bool capture_identity) {
    const auto points = active_pwl_points(
        gen.cost, gen.ncost, pg_lower, pg_upper);
    double sum = 0.0;
    double power = 0.0;
    for (const auto& point : points) {
        if (lambda_offset >= static_cast<int>(state.gen_lambda.size())) {
            throw std::runtime_error("generator lambda vector is too short");
        }
        const double lambda = state.gen_lambda[lambda_offset++];
        if (!std::isfinite(lambda)) {
            throw std::runtime_error(
                "state contains nonfinite candidate.gen_lambda value at "
                "position " + std::to_string(lambda_offset - 1));
        }
        sum += lambda;
        power += point.mw * lambda;
        update_category(
            report, report.max_variable_bound_violation,
            bound_violation(lambda, 0.0, 1.0), "variable_bound",
            capture_identity,
            [&] { return "gen:" + gen.source_key + ":lambda"; });
    }
    update_category(
        report, report.max_pwl_sum_residual, sum - 1.0, "pwl_sum",
        capture_identity,
        [&] { return "gen:" + gen.source_key; });
    update_category(
        report, report.max_pwl_power_residual,
        power - state.pg[generator_index], "pwl_power", capture_identity,
        [&] { return "gen:" + gen.source_key; });
}

void validate_load_pwl(
    ValidationReport& report,
    const Load& load,
    const AcState& state,
    std::size_t load_index,
    int& lambda_offset,
    bool capture_identity) {
    const auto points = active_pwl_points(
        load.cost, load.ncost, load.pd_min, load.pd_max);
    double sum = 0.0;
    double power = 0.0;
    for (const auto& point : points) {
        if (lambda_offset >= static_cast<int>(state.load_lambda.size())) {
            throw std::runtime_error("load lambda vector is too short");
        }
        const double lambda = state.load_lambda[lambda_offset++];
        if (!std::isfinite(lambda)) {
            throw std::runtime_error(
                "state contains nonfinite candidate.load_lambda value at "
                "position " + std::to_string(lambda_offset - 1));
        }
        sum += lambda;
        power += point.mw * lambda;
        update_category(
            report, report.max_variable_bound_violation,
            bound_violation(lambda, 0.0, 1.0), "variable_bound",
            capture_identity,
            [&] { return "load:" + load.source_key + ":lambda"; });
    }
    update_category(
        report, report.max_pwl_sum_residual, sum - 1.0, "pwl_sum",
        capture_identity,
        [&] { return "load:" + load.source_key; });
    update_category(
        report, report.max_pwl_power_residual,
        power - load.pd_nominal * state.demand_factor[load_index],
        "pwl_power", capture_identity,
        [&] { return "load:" + load.source_key; });
}

void validate_branch_ohms(
    ValidationReport& report,
    const Branch& branch,
    const AcState& state,
    std::size_t branch_index,
    bool capture_identity) {
    const int f = branch.from;
    const int t = branch.to;
    double from_g_self = branch.flow_from_g_self;
    double from_b_self = branch.flow_from_b_self;
    double to_g_self = branch.flow_to_g_self;
    double to_b_self = branch.flow_to_b_self;
    double from_cross_cos = branch.flow_from_cross_cos;
    double from_cross_sin = branch.flow_from_cross_sin;
    double to_cross_cos = branch.flow_to_cross_cos;
    double to_cross_sin = branch.flow_to_cross_sin;
    if (!branch.flow_coefficients_valid) {
        const double denominator =
            branch.r * branch.r + branch.x * branch.x;
        const double g = denominator > 1e-20
            ? branch.r / denominator : 0.0;
        const double b = denominator > 1e-20
            ? -branch.x / denominator : 0.0;
        const double tap_squared = branch.tap * branch.tap;
        const double tap_real = branch.tap * std::cos(branch.shift);
        const double tap_imag = branch.tap * std::sin(branch.shift);
        from_g_self = branch.transformer
            ? g / tap_squared + branch.g_fr
            : (g + branch.g_fr) / tap_squared;
        from_b_self = branch.transformer
            ? b / tap_squared + branch.b_fr
            : (b + branch.b_fr) / tap_squared;
        to_g_self = g + branch.g_to;
        to_b_self = b + branch.b_to;
        from_cross_cos =
            (-g * tap_real + b * tap_imag) / tap_squared;
        from_cross_sin =
            (-b * tap_real - g * tap_imag) / tap_squared;
        to_cross_cos =
            (-g * tap_real - b * tap_imag) / tap_squared;
        to_cross_sin =
            (-b * tap_real + g * tap_imag) / tap_squared;
    }
    const double voltage_product = state.vm[f] * state.vm[t];
    const double angle_delta = state.va[f] - state.va[t];
    const double cross_cos_ft = voltage_product * std::cos(angle_delta);
    const double cross_sin_ft = voltage_product * std::sin(angle_delta);
    const double expected_pf =
        from_g_self * state.vm[f] * state.vm[f] +
        from_cross_cos * cross_cos_ft +
        from_cross_sin * cross_sin_ft;
    const double expected_qf =
        -from_b_self * state.vm[f] * state.vm[f] -
        from_cross_sin * cross_cos_ft +
        from_cross_cos * cross_sin_ft;
    const double expected_pt =
        to_g_self * state.vm[t] * state.vm[t] +
        to_cross_cos * cross_cos_ft -
        to_cross_sin * cross_sin_ft;
    const double expected_qt =
        -to_b_self * state.vm[t] * state.vm[t] -
        to_cross_sin * cross_cos_ft -
        to_cross_cos * cross_sin_ft;
    const double ohms = std::max({
        std::abs(state.pf[branch_index] - expected_pf),
        std::abs(state.qf[branch_index] - expected_qf),
        std::abs(state.pt[branch_index] - expected_pt),
        std::abs(state.qt[branch_index] - expected_qt)});
    if (!std::isfinite(ohms)) {
        throw std::runtime_error(
            "state contains nonfinite rebuilt branch Ohm-law residual at " +
            branch.source_key);
    }
    update_category(
        report, report.max_ohms_residual, ohms, "ohms", capture_identity,
        [&] { return "branch:" + branch.source_key; });
}

}  // namespace

nlohmann::json ValidationReport::to_json() const {
    return {
        {"max_variable_bound_violation", max_variable_bound_violation},
        {"max_pwl_sum_residual", max_pwl_sum_residual},
        {"max_pwl_power_residual", max_pwl_power_residual},
        {"max_reference_angle_residual", max_reference_angle_residual},
        {"max_generator_residual", max_generator_residual},
        {"max_load_ramp_violation", max_load_ramp_violation},
        {"max_active_balance_residual", max_active_balance_residual},
        {"max_reactive_balance_residual", max_reactive_balance_residual},
        {"max_ohms_residual", max_ohms_residual},
        {"max_angle_violation", max_angle_violation},
        {"max_flow_limit_violation", max_flow_limit_violation},
        {"max_residual", max_residual},
        {"worst_category", worst_category},
        {"worst_identity", worst_identity},
    };
}

static ValidationReport validate_state_impl(
    const CaseData& data,
    ModelMode mode,
    const AcState& state,
    const std::vector<int>& fixed_status_argument,
    const std::optional<ContingencyContext>& contingency,
    bool skip_rebuilt_economic_and_ohms_checks,
    bool capture_identity,
    double incumbent_max_residual =
        std::numeric_limits<double>::infinity(),
    int preferred_balance_bus = -1,
    int preferred_branch = -1) {
    const std::size_t nb = data.buses.size();
    const std::size_t ng = data.generators.size();
    const std::size_t nd = data.loads.size();
    const std::size_t nl = data.branches.size();
    if ((mode == ModelMode::ContingencySoft) != contingency.has_value()) {
        throw std::runtime_error("contingency validation requires exactly one contingency context");
    }
    if (state.vm.size() != nb || state.va.size() != nb || state.pg.size() != ng ||
        state.qg.size() != ng || state.demand_factor.size() != nd || state.pf.size() != nl ||
        state.qf.size() != nl || state.pt.size() != nl || state.qt.size() != nl ||
        state.sm_slack.size() != nl) {
        throw std::runtime_error("state dimensions do not match case dimensions");
    }
    if ((!state.shunt_bs.empty() &&
         state.shunt_bs.size() != data.shunts.size()) ||
        (!state.shunt_steps.empty() &&
         state.shunt_steps.size() != data.shunts.size())) {
        throw std::runtime_error(
            "switched-shunt state dimensions do not match case dimensions");
    }
    if (mode != ModelMode::UnitCommitmentRelaxation &&
        (state.p_delta.size() != nb || state.q_delta.size() != nb)) {
        throw std::runtime_error("soft-balance state dimensions do not match case dimensions");
    }
    if (contingency) {
        const auto& base = contingency->base_state;
        if (base.vm.size() != nb || base.va.size() != nb || base.pg.size() != ng ||
            base.qg.size() != ng || base.demand_factor.size() != nd ||
            base.pf.size() != nl || base.qf.size() != nl ||
            base.pt.size() != nl || base.qt.size() != nl) {
            throw std::runtime_error("contingency base state dimensions do not match case dimensions");
        }
    }
    require_finite_state(state, "candidate.");
    if (contingency && !skip_rebuilt_economic_and_ohms_checks) {
        // Rebuilt trial validation is called repeatedly inside one fast
        // contingency solve. Its immutable base state was fully validated
        // before screening and is checked again by the final independent
        // validator, so rescanning every base vector in every local-search
        // trial adds no new acceptance protection.
        require_finite_state(contingency->base_state, "base.");
    }

    std::vector<int> default_fixed_status;
    const std::vector<int>* fixed_status = &fixed_status_argument;
    if (fixed_status_argument.empty()) {
        default_fixed_status.reserve(ng);
        for (const auto& gen : data.generators) {
            default_fixed_status.push_back(gen.status_prev);
        }
        fixed_status = &default_fixed_status;
    }
    if (fixed_status->size() != ng) {
        throw std::runtime_error("fixed commitment vector has the wrong length");
    }

    ValidationReport report;
    const auto rejection_proven = [&]() {
        return std::isfinite(incumbent_max_residual) &&
            report.max_residual + 1e-10 >= incumbent_max_residual;
    };

    for (std::size_t i = 0; i < nb; ++i) {
        update_category(report, report.max_variable_bound_violation,
            bound_violation(state.vm[i], data.buses[i].vmin, data.buses[i].vmax),
            "variable_bound", capture_identity,
            [&] { return "bus:" + data.buses[i].source_key + ":vm"; });
        if (data.buses[i].type == 3) {
            update_category(report, report.max_reference_angle_residual,
                state.va[i], "reference_angle", capture_identity,
                [&] { return "bus:" + data.buses[i].source_key; });
        }
        if (rejection_proven()) {
            return report;
        }
    }

    for (std::size_t i = 0; i < data.shunts.size(); ++i) {
        const auto& shunt = data.shunts[i];
        if (state.shunt_steps.empty()) {
            continue;
        }
        if (state.shunt_steps[i].size() <
            shunt.block_maximum_steps.size()) {
            throw std::runtime_error(
                "switched-shunt block vector is too short: " +
                shunt.source_key);
        }
        double expected_bs = 0.0;
        for (std::size_t block = 0;
             block < shunt.block_maximum_steps.size(); ++block) {
            const int step = state.shunt_steps[i][block];
            update_category(
                report, report.max_variable_bound_violation,
                bound_violation(
                    static_cast<double>(step), 0.0,
                    static_cast<double>(
                        shunt.block_maximum_steps[block])),
                "variable_bound", capture_identity,
                [&] {
                    return "shunt:" + shunt.source_key + ":block:" +
                        std::to_string(block);
                });
            expected_bs += static_cast<double>(step) *
                shunt.block_susceptance[block];
        }
        if (!shunt.dispatchable) {
            expected_bs = shunt.bs;
        }
        update_category(
            report, report.max_variable_bound_violation,
            effective_shunt_susceptance(data, state, static_cast<int>(i)) -
                expected_bs,
            "variable_bound", capture_identity,
            [&] { return "shunt:" + shunt.source_key + ":bs"; });
    }

    int gen_lambda_offset = 0;
    for (std::size_t i = 0; i < ng; ++i) {
        const auto& gen = data.generators[i];
        const bool outaged = contingency &&
            contingency->outaged_generator == static_cast<int>(i);
        const bool active = !outaged &&
            (mode == ModelMode::UnitCommitmentRelaxation || (*fixed_status)[i] == 1);
        double pg_lower = 0.0;
        double pg_upper = 0.0;
        double qg_lower = 0.0;
        double qg_upper = 0.0;
        if (mode == ModelMode::UnitCommitmentRelaxation) {
            pg_lower = std::min(0.0, gen.pmin);
            pg_upper = std::max(0.0, gen.pmax);
            qg_lower = std::min(0.0, gen.qmin);
            qg_upper = std::max(0.0, gen.qmax);
        } else if (active) {
            if (mode == ModelMode::ContingencySoft) {
                pg_lower = std::max(
                    gen.pmin,
                    contingency->base_state.pg[i] - data.delta_r_ctg * gen.prdmaxctg);
                pg_upper = std::min(
                    gen.pmax,
                    contingency->base_state.pg[i] + data.delta_r_ctg * gen.prumaxctg);
            } else {
                std::tie(pg_lower, pg_upper) = base_pg_bounds(gen, 1, data.delta_r);
            }
            qg_lower = gen.qmin;
            qg_upper = gen.qmax;
        }
        update_category(report, report.max_variable_bound_violation,
            bound_violation(state.pg[i], pg_lower, pg_upper), "variable_bound",
            capture_identity,
            [&] { return "gen:" + gen.source_key + ":pg"; });
        update_category(report, report.max_variable_bound_violation,
            bound_violation(state.qg[i], qg_lower, qg_upper), "variable_bound",
            capture_identity,
            [&] { return "gen:" + gen.source_key + ":qg"; });

        if (active && !skip_rebuilt_economic_and_ohms_checks) {
            validate_generator_pwl(
                report, gen, state, i, pg_lower, pg_upper,
                gen_lambda_offset, capture_identity);
        }

        if (mode == ModelMode::UnitCommitmentRelaxation) {
            if (state.commitment.size() != ng || state.startup.size() != ng || state.shutdown.size() != ng) {
                throw std::runtime_error("unit-commitment state vectors have the wrong length");
            }
            const double z = state.commitment[i];
            const double su = state.startup[i];
            const double sd = state.shutdown[i];
            update_category(report, report.max_generator_residual,
                positive_part(state.pg[i] - gen.pmax * z), "generator",
                capture_identity,
                [&] { return "gen:" + gen.source_key + ":pmax"; });
            update_category(report, report.max_generator_residual,
                positive_part(gen.pmin * z - state.pg[i]), "generator",
                capture_identity,
                [&] { return "gen:" + gen.source_key + ":pmin"; });
            update_category(report, report.max_generator_residual,
                positive_part(state.qg[i] - gen.qmax * z), "generator",
                capture_identity,
                [&] { return "gen:" + gen.source_key + ":qmax"; });
            update_category(report, report.max_generator_residual,
                positive_part(gen.qmin * z - state.qg[i]), "generator",
                capture_identity,
                [&] { return "gen:" + gen.source_key + ":qmin"; });
            const double previous = gen.status_prev == 0 ? gen.pmin : gen.pg_prev;
            update_category(report, report.max_generator_residual,
                positive_part(state.pg[i] - previous * z - data.delta_r * gen.prumax),
                "generator", capture_identity,
                [&] { return "gen:" + gen.source_key + ":ramp_up"; });
            update_category(report, report.max_generator_residual,
                positive_part(previous * z - data.delta_r * gen.prdmax - state.pg[i]),
                "generator", capture_identity,
                [&] { return "gen:" + gen.source_key + ":ramp_down"; });
            update_category(report, report.max_generator_residual,
                z - gen.status_prev - su + sd, "generator", capture_identity,
                [&] { return "gen:" + gen.source_key + ":transition"; });
            update_category(report, report.max_generator_residual,
                positive_part(su + sd - 1.0), "generator", capture_identity,
                [&] { return "gen:" + gen.source_key + ":exclusive"; });
        }
        if (rejection_proven()) {
            return report;
        }
    }

    int load_lambda_offset = 0;
    for (std::size_t i = 0; i < nd; ++i) {
        const auto& load = data.loads[i];
        std::pair<double, double> bounds;
        if (mode == ModelMode::BaseSoft) {
            bounds = base_load_bounds(load, data.delta_r);
        } else if (mode == ModelMode::ContingencySoft) {
            const double previous = load.pd_nominal * contingency->base_state.demand_factor[i];
            bounds = std::abs(load.pd_nominal) <= 1e-12
                ? std::pair<double, double>{load.tmin, load.tmax}
                : std::pair<double, double>{
                    std::max(load.tmin,
                        (previous - load.prdmaxctg * data.delta_r_ctg) / load.pd_nominal),
                    std::min(load.tmax,
                        (previous + load.prumaxctg * data.delta_r_ctg) / load.pd_nominal)};
        } else {
            bounds = {load.tmin, load.tmax};
        }
        update_category(report, report.max_variable_bound_violation,
            bound_violation(state.demand_factor[i], bounds.first, bounds.second),
            "variable_bound", capture_identity,
            [&] { return "load:" + load.source_key + ":factor"; });
        if (!skip_rebuilt_economic_and_ohms_checks) {
            validate_load_pwl(
                report, load, state, i, load_lambda_offset,
                capture_identity);
        }
        if (mode == ModelMode::UnitCommitmentRelaxation) {
            update_category(report, report.max_load_ramp_violation,
                positive_part(load.pd_nominal * state.demand_factor[i] - load.pd_prev - load.prumax * data.delta_r),
                "load_ramp", capture_identity,
                [&] { return "load:" + load.source_key + ":up"; });
            update_category(report, report.max_load_ramp_violation,
                positive_part(load.pd_prev - load.prdmax * data.delta_r - load.pd_nominal * state.demand_factor[i]),
                "load_ramp", capture_identity,
                [&] { return "load:" + load.source_key + ":down"; });
        } else if (mode == ModelMode::ContingencySoft) {
            const double previous = load.pd_nominal * contingency->base_state.demand_factor[i];
            const double current = load.pd_nominal * state.demand_factor[i];
            update_category(report, report.max_load_ramp_violation,
                positive_part(current - previous - load.prumaxctg * data.delta_r_ctg),
                "load_ramp", capture_identity,
                [&] {
                    return "load:" + load.source_key + ":contingency_up";
                });
            update_category(report, report.max_load_ramp_violation,
                positive_part(previous - load.prdmaxctg * data.delta_r_ctg - current),
                "load_ramp", capture_identity,
                [&] {
                    return "load:" + load.source_key + ":contingency_down";
                });
        }
        if (rejection_proven()) {
            return report;
        }
    }

    const auto validate_balance_bus = [&](std::size_t i) {
        const auto& bus = data.buses[i];
        double p = 0.0;
        double q = 0.0;
        for (int branch : bus.branches_from) {
            p += state.pf[branch];
            q += state.qf[branch];
        }
        for (int branch : bus.branches_to) {
            p += state.pt[branch];
            q += state.qt[branch];
        }
        for (int generator : bus.generators) {
            p -= state.pg[generator];
            q -= state.qg[generator];
        }
        for (int load : bus.loads) {
            p += data.loads[load].pd_nominal * state.demand_factor[load];
            q += data.loads[load].qd_nominal * state.demand_factor[load];
        }
        for (int shunt : bus.shunts) {
            p += data.shunts[shunt].gs * state.vm[i] * state.vm[i];
            q -= effective_shunt_susceptance(data, state, shunt) *
                state.vm[i] * state.vm[i];
        }
        const double p_residual = mode != ModelMode::UnitCommitmentRelaxation
            ? positive_part(std::abs(p) - state.p_delta[i])
            : std::abs(p);
        const double q_residual = mode != ModelMode::UnitCommitmentRelaxation
            ? positive_part(std::abs(q) - state.q_delta[i])
            : std::abs(q);
        if (mode != ModelMode::UnitCommitmentRelaxation) {
            update_category(report, report.max_variable_bound_violation,
                bound_violation(state.p_delta[i], 0.0, 0.5),
                "variable_bound", capture_identity,
                [&] { return "bus:" + bus.source_key + ":p_delta"; });
            update_category(report, report.max_variable_bound_violation,
                bound_violation(state.q_delta[i], 0.0, 0.5),
                "variable_bound", capture_identity,
                [&] { return "bus:" + bus.source_key + ":q_delta"; });
        }
        update_category(report, report.max_active_balance_residual,
            p_residual, "active_balance", capture_identity,
            [&] { return "bus:" + bus.source_key; });
        update_category(report, report.max_reactive_balance_residual,
            q_residual, "reactive_balance", capture_identity,
            [&] { return "bus:" + bus.source_key; });
    };
    if (preferred_balance_bus >= 0 &&
        preferred_balance_bus < static_cast<int>(nb)) {
        validate_balance_bus(static_cast<std::size_t>(preferred_balance_bus));
        if (rejection_proven()) {
            return report;
        }
    }
    for (std::size_t i = 0; i < nb; ++i) {
        if (static_cast<int>(i) == preferred_balance_bus) {
            continue;
        }
        validate_balance_bus(i);
        if (rejection_proven()) {
            return report;
        }
    }

    const auto validate_branch = [&](std::size_t i) {
        const auto& branch = data.branches[i];
        const bool outaged = contingency &&
            contingency->outaged_branch == static_cast<int>(i);
        const bool unavailable = branch.status == 0 || outaged;
        const double rating = mode == ModelMode::ContingencySoft
            ? branch.rate_c : branch.rate_a;
        const double flow_lower = unavailable ? 0.0 : -rating;
        const double flow_upper = unavailable ? 0.0 : rating;
        for (const auto& item : std::array<std::pair<const char*, double>, 4>{{
                 {"pf", state.pf[i]}, {"qf", state.qf[i]},
                 {"pt", state.pt[i]}, {"qt", state.qt[i]}}}) {
            update_category(report, report.max_variable_bound_violation,
                bound_violation(item.second, flow_lower, flow_upper),
                "variable_bound", capture_identity,
                [&] {
                    return "branch:" + branch.source_key + ":" + item.first;
                });
        }
        update_category(report, report.max_variable_bound_violation,
            bound_violation(
                state.sm_slack[i], 0.0,
                unavailable ? 0.0 : data.sm_vio_limit),
            "variable_bound", capture_identity,
            [&] { return "branch:" + branch.source_key + ":sm_slack"; });
        if (unavailable) {
            return;
        }
        const int f = branch.from;
        const int t = branch.to;
        if (!skip_rebuilt_economic_and_ohms_checks) {
            validate_branch_ohms(
                report, branch, state, i, capture_identity);
        }

        const double start_delta = mode == ModelMode::ContingencySoft
            ? contingency->base_state.va[f] - contingency->base_state.va[t]
            : data.buses[f].va_start - data.buses[t].va_start;
        if (start_delta >= branch.angmin && start_delta <= branch.angmax) {
            const double angle = state.va[f] - state.va[t];
            update_category(report, report.max_angle_violation,
                std::max(positive_part(angle - branch.angmax), positive_part(branch.angmin - angle)),
                "angle", capture_identity,
                [&] { return "branch:" + branch.source_key; });
        }
        const double from_squared = state.pf[i] * state.pf[i] + state.qf[i] * state.qf[i];
        const double to_squared = state.pt[i] * state.pt[i] + state.qt[i] * state.qt[i];
        const double from_limit = branch.transformer
            ? rating * rating * std::pow(1.0 + state.sm_slack[i], 2)
            : rating * rating * std::pow(state.vm[f] + state.sm_slack[i], 2);
        const double to_limit = branch.transformer
            ? rating * rating * std::pow(1.0 + state.sm_slack[i], 2)
            : rating * rating * std::pow(state.vm[t] + state.sm_slack[i], 2);
        update_category(report, report.max_flow_limit_violation,
            std::max(positive_part(from_squared - from_limit), positive_part(to_squared - to_limit)),
            "flow_limit", capture_identity,
            [&] { return "branch:" + branch.source_key; });
    };
    if (preferred_branch >= 0 && preferred_branch < static_cast<int>(nl)) {
        validate_branch(static_cast<std::size_t>(preferred_branch));
        if (rejection_proven()) {
            return report;
        }
    }
    for (std::size_t i = 0; i < nl; ++i) {
        if (static_cast<int>(i) == preferred_branch) {
            continue;
        }
        validate_branch(i);
        if (rejection_proven()) {
            return report;
        }
    }

    return report;
}

ValidationReport validate_state(
    const CaseData& data,
    ModelMode mode,
    const AcState& state,
    const std::vector<int>& fixed_status,
    const std::optional<ContingencyContext>& contingency) {
    return validate_state_impl(
        data, mode, state, fixed_status, contingency, false, true);
}

ValidationReport validate_rebuilt_contingency_trial(
    const CaseData& data,
    const AcState& state,
    const std::vector<int>& fixed_status,
    const ContingencyContext& contingency) {
    return validate_state_impl(
        data, ModelMode::ContingencySoft, state, fixed_status,
        contingency, true, false);
}

ValidationReport validate_rebuilt_contingency_trial_until_rejected(
    const CaseData& data,
    const AcState& state,
    const std::vector<int>& fixed_status,
    const ContingencyContext& contingency,
    double incumbent_max_residual,
    int preferred_balance_bus,
    int preferred_branch) {
    return validate_state_impl(
        data, ModelMode::ContingencySoft, state, fixed_status,
        contingency, true, false, incumbent_max_residual,
        preferred_balance_bus, preferred_branch);
}

ValidationReport validate_rebuilt_contingency_predictor(
    const CaseData& data,
    const AcState& state,
    const std::vector<int>& fixed_status,
    const ContingencyContext& contingency) {
    return validate_state_impl(
        data, ModelMode::ContingencySoft, state, fixed_status,
        contingency, true, true);
}

ValidationReport validate_rebuilt_contingency_economic_and_ohms(
    const CaseData& data,
    const AcState& state,
    const std::vector<int>& fixed_status_argument,
    const ContingencyContext& contingency,
    const ValidationReport& physical_report) {
    const std::size_t nb = data.buses.size();
    const std::size_t ng = data.generators.size();
    const std::size_t nd = data.loads.size();
    const std::size_t nl = data.branches.size();
    if (state.vm.size() != nb || state.va.size() != nb || state.pg.size() != ng ||
        state.demand_factor.size() != nd || state.pf.size() != nl ||
        state.qf.size() != nl || state.pt.size() != nl || state.qt.size() != nl ||
        contingency.base_state.pg.size() != ng) {
        throw std::runtime_error(
            "rebuilt economic/Ohm validation dimensions do not match case");
    }

    std::vector<int> default_fixed_status;
    const std::vector<int>* fixed_status = &fixed_status_argument;
    if (fixed_status_argument.empty()) {
        default_fixed_status.reserve(ng);
        for (const auto& gen : data.generators) {
            default_fixed_status.push_back(gen.status_prev);
        }
        fixed_status = &default_fixed_status;
    }
    if (fixed_status->size() != ng) {
        throw std::runtime_error("fixed commitment vector has the wrong length");
    }

    ValidationReport report = physical_report;
    int gen_lambda_offset = 0;
    for (std::size_t i = 0; i < ng; ++i) {
        const auto& gen = data.generators[i];
        const bool outaged =
            contingency.outaged_generator == static_cast<int>(i);
        const bool active = !outaged && (*fixed_status)[i] == 1;
        if (!active) {
            continue;
        }
        const double pg_lower = std::max(
            gen.pmin,
            contingency.base_state.pg[i] -
                data.delta_r_ctg * gen.prdmaxctg);
        const double pg_upper = std::min(
            gen.pmax,
            contingency.base_state.pg[i] +
                data.delta_r_ctg * gen.prumaxctg);
        validate_generator_pwl(
            report, gen, state, i, pg_lower, pg_upper,
            gen_lambda_offset, true);
    }

    int load_lambda_offset = 0;
    for (std::size_t i = 0; i < nd; ++i) {
        validate_load_pwl(
            report, data.loads[i], state, i, load_lambda_offset, true);
    }

    for (std::size_t i = 0; i < nl; ++i) {
        const auto& branch = data.branches[i];
        const bool outaged =
            contingency.outaged_branch == static_cast<int>(i);
        if (branch.status == 0 || outaged) {
            continue;
        }
        validate_branch_ohms(report, branch, state, i, true);
    }
    return report;
}

}  // namespace gravityx
