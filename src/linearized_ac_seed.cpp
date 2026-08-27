#include <highs/Highs.h>

#include "gravityx/linearized_ac_seed.hpp"

#include "gravityx/state_io.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace gravityx {
namespace {

struct FlowValues {
    double pf{};
    double qf{};
    double pt{};
    double qt{};
};

struct AffineFlow {
    double constant{};
    double vm_from{};
    double vm_to{};
    double va_from{};
    double va_to{};
};

struct LinearizedBranch {
    AffineFlow pf;
    AffineFlow qf;
    AffineFlow pt;
    AffineFlow qt;
};

struct SparseRow {
    double lower{};
    double upper{};
    std::vector<std::pair<HighsInt, double>> entries;
};

struct LinearizedModelAudit {
    bool passed{true};
    std::string failure;
    int tiny_matrix_entries_removed{};
    double maximum_tiny_matrix_entry_removed{};
};

double piecewise_linear_value(
    const std::vector<PwlPoint>& points,
    double power) {
    if (points.size() < 2) {
        throw std::runtime_error(
            "piecewise-linear value requires at least two points");
    }
    std::size_t segment = points.size() - 2;
    for (std::size_t i = 0; i + 1 < points.size(); ++i) {
        if (power <= points[i + 1].mw + 1e-12) {
            segment = i;
            break;
        }
    }
    const auto& left = points[segment];
    const auto& right = points[segment + 1];
    const double width = right.mw - left.mw;
    if (std::abs(width) <= 1e-14) {
        return left.cost;
    }
    return left.cost +
        (right.cost - left.cost) * (power - left.mw) / width;
}

double piecewise_linear_subgradient(
    const std::vector<PwlPoint>& points,
    double power) {
    if (points.size() < 2) {
        throw std::runtime_error(
            "piecewise-linear slope requires at least two points");
    }
    const auto slope = [&] (std::size_t segment) {
        const double width =
            points[segment + 1].mw - points[segment].mw;
        if (std::abs(width) <= 1e-14) {
            return 0.0;
        }
        return (points[segment + 1].cost - points[segment].cost) / width;
    };
    for (std::size_t point = 1; point + 1 < points.size(); ++point) {
        if (std::abs(power - points[point].mw) <= 1e-10) {
            return 0.5 * (slope(point - 1) + slope(point));
        }
    }
    for (std::size_t segment = 0;
         segment + 1 < points.size(); ++segment) {
        if (power <= points[segment + 1].mw + 1e-12) {
            return slope(segment);
        }
    }
    return slope(points.size() - 2);
}

FlowValues branch_flows(
    const Branch& branch,
    double vm_from,
    double vm_to,
    double va_from,
    double va_to) {
    const double denominator = branch.r * branch.r + branch.x * branch.x;
    const double g = denominator > 1e-20 ? branch.r / denominator : 0.0;
    const double b = denominator > 1e-20 ? -branch.x / denominator : 0.0;
    const double tm = branch.tap;
    if (std::abs(tm) <= 1e-12) {
        throw std::runtime_error(
            "zero tap ratio in linearized AC seed: " + branch.source_key);
    }
    const double tm2 = tm * tm;
    const double tr = tm * std::cos(branch.shift);
    const double ti = tm * std::sin(branch.shift);
    const double delta_ft = va_from - va_to;
    const double delta_tf = -delta_ft;
    const double cross_ft = vm_from * vm_to;
    const double from_g_self = branch.transformer
        ? g / tm2 + branch.g_fr : (g + branch.g_fr) / tm2;
    const double from_b_self = branch.transformer
        ? b / tm2 + branch.b_fr : (b + branch.b_fr) / tm2;

    FlowValues value;
    value.pf = from_g_self * vm_from * vm_from
        + ((-g * tr + b * ti) / tm2) * cross_ft * std::cos(delta_ft)
        + ((-b * tr - g * ti) / tm2) * cross_ft * std::sin(delta_ft);
    value.qf = -from_b_self * vm_from * vm_from
        - ((-b * tr + g * ti) / tm2) * cross_ft * std::cos(delta_ft)
        + ((-g * tr + b * ti) / tm2) * cross_ft * std::sin(delta_ft);
    value.pt = (g + branch.g_to) * vm_to * vm_to
        + ((-g * tr - b * ti) / tm2) * cross_ft * std::cos(delta_tf)
        + ((-b * tr + g * ti) / tm2) * cross_ft * std::sin(delta_tf);
    value.qt = -(b + branch.b_to) * vm_to * vm_to
        - ((-b * tr + g * ti) / tm2) * cross_ft * std::cos(delta_tf)
        + ((-g * tr - b * ti) / tm2) * cross_ft * std::sin(delta_tf);
    return value;
}

double flow_component(const FlowValues& value, int component) {
    switch (component) {
        case 0: return value.pf;
        case 1: return value.qf;
        case 2: return value.pt;
        case 3: return value.qt;
        default: throw std::runtime_error("invalid branch-flow component");
    }
}

LinearizedBranch linearize_branch(
    const Branch& branch,
    double vm_from,
    double vm_to,
    double va_from,
    double va_to) {
    constexpr double h = 1e-6;
    const auto base = branch_flows(
        branch, vm_from, vm_to, va_from, va_to);
    std::array<std::array<double, 4>, 4> derivative{};
    const std::array<double, 4> point{vm_from, vm_to, va_from, va_to};
    for (int variable = 0; variable < 4; ++variable) {
        auto plus = point;
        auto minus = point;
        plus[variable] += h;
        minus[variable] -= h;
        const auto plus_flow = branch_flows(
            branch, plus[0], plus[1], plus[2], plus[3]);
        const auto minus_flow = branch_flows(
            branch, minus[0], minus[1], minus[2], minus[3]);
        for (int component = 0; component < 4; ++component) {
            derivative[component][variable] =
                (flow_component(plus_flow, component)
                    - flow_component(minus_flow, component)) / (2.0 * h);
        }
    }
    LinearizedBranch result;
    std::array<AffineFlow*, 4> output{
        &result.pf, &result.qf, &result.pt, &result.qt};
    for (int component = 0; component < 4; ++component) {
        auto& affine = *output[component];
        affine.vm_from = derivative[component][0];
        affine.vm_to = derivative[component][1];
        affine.va_from = derivative[component][2];
        affine.va_to = derivative[component][3];
        affine.constant = flow_component(base, component)
            - affine.vm_from * vm_from
            - affine.vm_to * vm_to
            - affine.va_from * va_from
            - affine.va_to * va_to;
    }
    return result;
}

void append(SparseRow& row, HighsInt column, double value) {
    // Preserve non-finite values until the explicit model preflight so they
    // are reported rather than silently discarded by a floating-point
    // comparison with NaN.
    if (!std::isfinite(value) || std::abs(value) > 1e-14) {
        row.entries.emplace_back(column, value);
    }
}

void normalize(SparseRow& row) {
    std::sort(row.entries.begin(), row.entries.end(),
        [](const auto& left, const auto& right) {
            return left.first < right.first;
        });
    std::vector<std::pair<HighsInt, double>> combined;
    for (const auto& entry : row.entries) {
        if (!combined.empty() && combined.back().first == entry.first) {
            combined.back().second += entry.second;
        } else {
            combined.push_back(entry);
        }
    }
    combined.erase(
        std::remove_if(combined.begin(), combined.end(),
            [](const auto& entry) { return std::abs(entry.second) <= 1e-14; }),
        combined.end());
    row.entries = std::move(combined);
}

LinearizedModelAudit audit_and_prune_model(
    const std::vector<double>& lower,
    const std::vector<double>& upper,
    const std::vector<double>& cost,
    std::vector<SparseRow>& rows,
    double small_matrix_value) {
    LinearizedModelAudit audit;
    if (lower.size() != upper.size() || lower.size() != cost.size()) {
        audit.passed = false;
        audit.failure = "column vector sizes disagree";
        return audit;
    }
    if (!std::isfinite(small_matrix_value) || small_matrix_value <= 0.0) {
        audit.passed = false;
        audit.failure = "small matrix value is not positive and finite";
        return audit;
    }
    for (std::size_t column = 0; column < lower.size(); ++column) {
        if (std::isnan(lower[column]) ||
            std::isnan(upper[column]) ||
            !std::isfinite(cost[column])) {
            audit.passed = false;
            audit.failure = "non-finite column data at column " +
                std::to_string(column);
            return audit;
        }
        if (lower[column] > upper[column]) {
            audit.passed = false;
            audit.failure = "empty column interval at column " +
                std::to_string(column);
            return audit;
        }
    }
    for (std::size_t row_index = 0; row_index < rows.size(); ++row_index) {
        auto& row = rows[row_index];
        normalize(row);
        if (std::isnan(row.lower) || std::isnan(row.upper)) {
            audit.passed = false;
            audit.failure = "non-finite row bound at row " +
                std::to_string(row_index);
            return audit;
        }
        if (row.lower > row.upper) {
            audit.passed = false;
            audit.failure = "empty row interval at row " +
                std::to_string(row_index);
            return audit;
        }
        for (const auto& [column, value] : row.entries) {
            if (column < 0 ||
                column >= static_cast<HighsInt>(lower.size())) {
                audit.passed = false;
                audit.failure = "column index out of range at row " +
                    std::to_string(row_index);
                return audit;
            }
            if (!std::isfinite(value)) {
                audit.passed = false;
                audit.failure = "non-finite matrix coefficient at row " +
                    std::to_string(row_index) + ", column " +
                    std::to_string(column);
                return audit;
            }
        }
        row.entries.erase(
            std::remove_if(
                row.entries.begin(), row.entries.end(),
                [&](const auto& entry) {
                    const double magnitude = std::abs(entry.second);
                    if (magnitude <= small_matrix_value) {
                        ++audit.tiny_matrix_entries_removed;
                        audit.maximum_tiny_matrix_entry_removed = std::max(
                            audit.maximum_tiny_matrix_entry_removed,
                            magnitude);
                        return true;
                    }
                    return false;
                }),
            row.entries.end());
    }
    return audit;
}

}  // namespace

nlohmann::json LinearizedAcSeedResult::to_json(bool include_state) const {
    nlohmann::json value = {
        {"success", success},
        {"wall_seconds", wall_seconds},
        {"economic_objective", economic_objective},
        {"compact_economic_objective", compact_economic_objective},
        {"projected_balance_slack", projected_balance_slack},
        {"branch_security_rows_omitted", branch_security_rows_omitted},
        {"branch_security_subset_count", branch_security_subset_count},
        {"feasibility_only", feasibility_only},
        {"elastic_balance_phase_one", elastic_balance_phase_one},
        {"primal_start_attempted", primal_start_attempted},
        {"primal_start_accepted", primal_start_accepted},
        {"primal_start_status", primal_start_status},
        {"primal_basis_attempted", primal_basis_attempted},
        {"primal_basis_accepted", primal_basis_accepted},
        {"primal_basis_status", primal_basis_status},
        {"presolve_enabled", presolve_enabled},
        {"primal_simplex_bound_perturbation_multiplier",
         primal_simplex_bound_perturbation_multiplier},
        {"simplex_strategy", simplex_strategy},
        {"maximum_column_scale", maximum_column_scale},
        {"maximum_row_scale", maximum_row_scale},
        {"objective_scale", objective_scale},
        {"voltage_trust_radius", voltage_trust_radius},
        {"angle_trust_radius", angle_trust_radius},
        {"projected_reference_voltage_count",
         projected_reference_voltage_count},
        {"maximum_reference_voltage_projection",
         maximum_reference_voltage_projection},
        {"maximum_balance_slack", maximum_balance_slack},
        {"total_balance_slack", total_balance_slack},
        {"solution_value_valid", solution_value_valid},
        {"info_valid", info_valid},
        {"accepted_feasible_nonoptimal_phase_one",
         accepted_feasible_nonoptimal_phase_one},
        {"accepted_feasible_nonoptimal_economic",
         accepted_feasible_nonoptimal_economic},
        {"accepted_approximate_economic_direction",
         accepted_approximate_economic_direction},
        {"resident_segment_count", resident_segment_count},
        {"feasible_segment_snapshot_count",
         feasible_segment_snapshot_count},
        {"recovered_feasible_segment_snapshot",
         recovered_feasible_segment_snapshot},
        {"canonicalized_segment_snapshot",
         canonicalized_segment_snapshot},
        {"canonicalized_snapshot_max_primal_infeasibility",
         canonicalized_snapshot_max_primal_infeasibility},
        {"canonicalized_snapshot_objective",
         canonicalized_snapshot_objective},
        {"terminal_run_status", terminal_run_status},
        {"terminal_model_status", terminal_model_status},
        {"resident_segments", resident_segments},
        {"time_limit_seconds", time_limit_seconds},
        {"ipm_optimality_tolerance", ipm_optimality_tolerance},
        {"row_count", row_count},
        {"column_count", column_count},
        {"nonzero_count", nonzero_count},
        {"model_preflight_passed", model_preflight_passed},
        {"model_preflight_failure", model_preflight_failure},
        {"tiny_matrix_entries_removed", tiny_matrix_entries_removed},
        {"maximum_tiny_matrix_entry_removed",
         maximum_tiny_matrix_entry_removed},
        {"small_matrix_value", small_matrix_value},
        {"add_vars_status", add_vars_status},
        {"change_cols_cost_status", change_cols_cost_status},
        {"add_rows_status", add_rows_status},
        {"model_load_warning", model_load_warning},
        {"model_construction_success", model_construction_success},
        {"model_load_failure_call", model_load_failure_call},
        {"run_status", run_status},
        {"model_status", model_status},
        {"primal_solution_status", primal_solution_status},
        {"num_primal_infeasibilities", num_primal_infeasibilities},
        {"max_primal_infeasibility", max_primal_infeasibility},
        {"status", status},
        {"iterations", iterations},
        {"objective", objective},
    };
    if (include_state) {
        value["state"] = ac_state_to_json(state);
    }
    return value;
}

LinearizedAcSeedResult solve_linearized_ac_seed(
    const CaseData& data,
    const AcState& reference,
    const std::vector<int>& commitment,
    double balance_slack_limit,
    const std::optional<ContingencyContext>& contingency,
    bool project_balance_slack,
    bool request_lightweight_large_seed,
    double time_limit_seconds,
    bool omit_branch_security_rows,
    bool feasibility_only,
    const std::vector<int>& branch_security_subset,
    bool force_projected_balance_phase_one,
    double voltage_trust_radius_override,
    double angle_trust_radius_override,
    bool economic_objective) {
    const auto wall_start = std::chrono::steady_clock::now();
    const int nb = static_cast<int>(data.buses.size());
    const int ng = static_cast<int>(data.generators.size());
    const int nd = static_cast<int>(data.loads.size());
    const int nl = static_cast<int>(data.branches.size());
    if (commitment.size() != data.generators.size() ||
        reference.vm.size() != data.buses.size() ||
        reference.va.size() != data.buses.size() ||
        reference.pg.size() != data.generators.size() ||
        reference.qg.size() != data.generators.size() ||
        reference.demand_factor.size() != data.loads.size()) {
        throw std::runtime_error("linearized AC seed received invalid dimensions");
    }
    if (economic_objective &&
        (reference.pf.size() != data.branches.size() ||
         reference.qf.size() != data.branches.size() ||
         reference.pt.size() != data.branches.size() ||
         reference.qt.size() != data.branches.size())) {
        throw std::runtime_error(
            "linearized economic seed requires rebuilt branch-flow state");
    }
    if (!(balance_slack_limit > 0.0 && balance_slack_limit < 0.5)) {
        throw std::runtime_error("linearized AC balance slack must be in (0, 0.5)");
    }
    if (!std::isfinite(time_limit_seconds) || time_limit_seconds <= 0.0) {
        throw std::runtime_error("linearized AC time limit must be positive");
    }
    if (economic_objective && feasibility_only) {
        throw std::runtime_error(
            "linearized economic objective is incompatible with feasibility-only mode");
    }
    // A caller may omit all but a selected set of branch-security rows for
    // either a base or contingency Phase-I model.  Such a reduced model is
    // candidate generation only; callers must rebuild nonlinear branch flows
    // and apply the complete independent validator before acceptance.
    std::vector<unsigned char> selected_branch_security(
        static_cast<std::size_t>(nl), 0);
    int branch_security_subset_count = 0;
    for (const int branch : branch_security_subset) {
        if (branch < 0 || branch >= nl) {
            throw std::runtime_error(
                "linearized AC branch-security subset index is out of range");
        }
        if (!selected_branch_security[static_cast<std::size_t>(branch)]) {
            selected_branch_security[static_cast<std::size_t>(branch)] = 1;
            ++branch_security_subset_count;
        }
    }
    if (contingency) {
        const auto& base = contingency->effective_base_state();
        if (base.pg.size() != data.generators.size() ||
            base.demand_factor.size() != data.loads.size() ||
            base.va.size() != data.buses.size()) {
            throw std::runtime_error(
                "linearized AC contingency base state has invalid dimensions");
        }
    }
    const int outaged_generator = contingency
        ? contingency->outaged_generator : -1;
    const int outaged_branch = contingency
        ? contingency->outaged_branch : -1;
    // At 8k-bus scale this LP is only a nonlinear-feasibility seed.  Including
    // every linearized branch box and angle row more than doubles its size and
    // dominated the five-minute contingency budget.  The returned candidate is
    // never accepted on this reduced model: sparse Newton follows immediately,
    // then the independent nonlinear validator checks every source branch flow
    // and angle constraint.  Smaller cases retain the original full seed LP.
    const char* full_seed_override = std::getenv("GRAVITYX_FULL_LINEAR_SEED");
    const bool force_full_seed = full_seed_override != nullptr &&
        std::string(full_seed_override) != "0";
    const bool lightweight_large_seed = !force_full_seed &&
        ((contingency.has_value() && nb >= 8000) ||
         request_lightweight_large_seed);
    const bool compact_economic_objective =
        economic_objective && omit_branch_security_rows;
    // A zero-objective feasibility LP gave large-case simplex no useful
    // direction and no readily available feasible basis.  Give the 19k-bus
    // Phase-I model explicit positive/negative balance violations instead.
    // Their L1 objective measures progress, and their identity columns let us
    // construct a valid primal start analytically.  The elastic solution is
    // only a nonlinear-repair seed; the unchanged complete validator remains
    // the sole acceptance gate.
    // Use an elastic objective only when balance itself is the unresolved
    // large-case defect.  Once exact violated security rows are supplied, the
    // reference already satisfies the soft balance band; a much smaller local
    // zero-objective feasibility model is sufficient and avoids the large
    // elastic simplex system.
    const bool elastic_balance_phase_one = feasibility_only && nb >= 16000 &&
        branch_security_subset.empty() &&
        !force_projected_balance_phase_one;
    const bool targeted_security_repair = feasibility_only && nb >= 16000 &&
        !branch_security_subset.empty();
    const double default_voltage_trust_radius =
        targeted_security_repair ? 0.01 : 0.05;
    const double default_angle_trust_radius =
        targeted_security_repair ? 0.006 : 0.15;
    const double voltage_trust_radius = voltage_trust_radius_override > 0.0
        ? voltage_trust_radius_override : default_voltage_trust_radius;
    const double angle_trust_radius = angle_trust_radius_override > 0.0
        ? angle_trust_radius_override : default_angle_trust_radius;
    // Sparse Newton is allowed to return an intermediate point outside a
    // source voltage bound.  Centering a local LP trust box on that point can
    // create an empty interval before HiGHS sees the model.  Project only the
    // LP linearization voltage into [VMIN, VMAX]; every returned state is still
    // checked against the unchanged complete nonlinear validator.
    std::vector<double> linearization_vm = reference.vm;
    int projected_reference_voltage_count = 0;
    double maximum_reference_voltage_projection = 0.0;
    for (int i = 0; i < nb; ++i) {
        const double projected = std::clamp(
            reference.vm[i], data.buses[i].vmin, data.buses[i].vmax);
        const double movement = std::abs(projected - reference.vm[i]);
        if (movement > 1e-12) {
            ++projected_reference_voltage_count;
            maximum_reference_voltage_projection = std::max(
                maximum_reference_voltage_projection, movement);
        }
        linearization_vm[i] = projected;
    }
    const bool projected_balance_slack =
        project_balance_slack &&
        (lightweight_large_seed || feasibility_only) &&
        !elastic_balance_phase_one;
    const int balance_slack_column_count =
        (!feasibility_only || elastic_balance_phase_one) ? nb : 0;

    const int vm_offset = 0;
    const int va_offset = vm_offset + nb;
    const int pg_offset = va_offset + nb;
    const int qg_offset = pg_offset + ng;
    const int demand_offset = qg_offset + ng;
    const int p_pos_offset = demand_offset + nd;
    const int p_neg_offset = p_pos_offset + balance_slack_column_count;
    const int q_pos_offset = p_neg_offset + balance_slack_column_count;
    const int q_neg_offset = q_pos_offset + balance_slack_column_count;
    const bool movement_objective = !feasibility_only && !economic_objective;
    const int dpg_offset = q_neg_offset + balance_slack_column_count;
    const int dqg_offset = dpg_offset + (movement_objective ? ng : 0);
    const int dload_offset = dqg_offset + (movement_objective ? ng : 0);
    const int gen_cost_offset = dload_offset + (movement_objective ? nd : 0);
    const int load_value_offset =
        gen_cost_offset +
        (economic_objective && !compact_economic_objective ? ng : 0);
    const int sm_offset =
        load_value_offset +
        (economic_objective && !compact_economic_objective ? nd : 0);
    const int compact_vm_down_offset = sm_offset +
        (economic_objective && !compact_economic_objective ? nl : 0);
    const int compact_va_down_offset = compact_vm_down_offset +
        (compact_economic_objective ? nb : 0);
    const int compact_pg_down_offset = compact_va_down_offset +
        (compact_economic_objective ? nb : 0);
    const int compact_qg_down_offset = compact_pg_down_offset +
        (compact_economic_objective ? ng : 0);
    const int compact_demand_down_offset = compact_qg_down_offset +
        (compact_economic_objective ? ng : 0);
    const int column_count = compact_demand_down_offset +
        (compact_economic_objective ? nd : 0);

    std::vector<double> lower(column_count, 0.0);
    std::vector<double> upper(column_count, kHighsInf);
    std::vector<double> cost(column_count, 0.0);
    for (int i = 0; i < nb; ++i) {
        lower[vm_offset + i] = lightweight_large_seed
            ? std::max(
                data.buses[i].vmin,
                linearization_vm[i] - voltage_trust_radius)
            : data.buses[i].vmin;
        upper[vm_offset + i] = lightweight_large_seed
            ? std::min(
                data.buses[i].vmax,
                linearization_vm[i] + voltage_trust_radius)
            : data.buses[i].vmax;
        lower[va_offset + i] = lightweight_large_seed
            ? std::max(
                -10.0,
                reference.va[i] - angle_trust_radius)
            : -10.0;
        upper[va_offset + i] = lightweight_large_seed
            ? std::min(
                10.0,
                reference.va[i] + angle_trust_radius)
            : 10.0;
        if (data.buses[i].type == 3) {
            lower[va_offset + i] = 0.0;
            upper[va_offset + i] = 0.0;
        }
        if (!feasibility_only || elastic_balance_phase_one) {
            const double explicit_slack_upper = elastic_balance_phase_one
                ? kHighsInf
                : (lightweight_large_seed && !economic_objective
                    ? 0.0 : balance_slack_limit);
            upper[p_pos_offset + i] = explicit_slack_upper;
            upper[p_neg_offset + i] = explicit_slack_upper;
            upper[q_pos_offset + i] = explicit_slack_upper;
            upper[q_neg_offset + i] = explicit_slack_upper;
            if (elastic_balance_phase_one) {
                cost[p_pos_offset + i] = 1.0;
                cost[p_neg_offset + i] = 1.0;
                cost[q_pos_offset + i] = 1.0;
                cost[q_neg_offset + i] = 1.0;
            } else if (economic_objective) {
                cost[p_pos_offset + i] =
                    data.delta * data.p_delta_cost_approx;
                cost[p_neg_offset + i] =
                    data.delta * data.p_delta_cost_approx;
                cost[q_pos_offset + i] =
                    data.delta * data.q_delta_cost_approx;
                cost[q_neg_offset + i] =
                    data.delta * data.q_delta_cost_approx;
            } else if (!lightweight_large_seed) {
                cost[p_pos_offset + i] = 1e6;
                cost[p_neg_offset + i] = 1e6;
                cost[q_pos_offset + i] = 1e6;
                cost[q_neg_offset + i] = 1e6;
            }
        }
    }
    for (int i = 0; i < ng; ++i) {
        const auto& gen = data.generators[i];
        if (commitment[i] == 0 || i == outaged_generator) {
            lower[pg_offset + i] = upper[pg_offset + i] = 0.0;
            lower[qg_offset + i] = upper[qg_offset + i] = 0.0;
        } else {
            if (contingency) {
                lower[pg_offset + i] = std::max(
                    gen.pmin,
                    contingency->effective_base_state().pg[i]
                        - data.delta_r_ctg * gen.prdmaxctg);
                upper[pg_offset + i] = std::min(
                    gen.pmax,
                    contingency->effective_base_state().pg[i]
                        + data.delta_r_ctg * gen.prumaxctg);
            } else {
                const double previous =
                    gen.status_prev == 0 ? gen.pmin : gen.pg_prev;
                lower[pg_offset + i] = std::max(
                    gen.pmin, previous - data.delta_r * gen.prdmax);
                upper[pg_offset + i] = std::min(
                    gen.pmax, previous + data.delta_r * gen.prumax);
            }
            if (lower[pg_offset + i] > upper[pg_offset + i] + 1e-12) {
                throw std::runtime_error(
                    "empty generator interval in linearized AC seed: "
                    + gen.source_key);
            }
            lower[qg_offset + i] = gen.qmin;
            upper[qg_offset + i] = gen.qmax;
        }
        if (movement_objective) {
            cost[dpg_offset + i] = 1.0;
            cost[dqg_offset + i] = 0.1;
        }
    }
    for (int i = 0; i < nd; ++i) {
        const auto& load = data.loads[i];
        lower[demand_offset + i] = load.tmin;
        upper[demand_offset + i] = load.tmax;
        if (std::abs(load.pd_nominal) > 1e-12) {
            const double previous = contingency
                ? load.pd_nominal
                    * contingency->effective_base_state().demand_factor[i]
                : load.pd_prev;
            const double down = contingency
                ? load.prdmaxctg * data.delta_r_ctg
                : load.prdmax * data.delta_r;
            const double up = contingency
                ? load.prumaxctg * data.delta_r_ctg
                : load.prumax * data.delta_r;
            lower[demand_offset + i] = std::max(
                lower[demand_offset + i],
                (previous - down) / load.pd_nominal);
            upper[demand_offset + i] = std::min(
                upper[demand_offset + i],
                (previous + up) / load.pd_nominal);
        }
        if (lower[demand_offset + i] > upper[demand_offset + i] + 1e-12) {
            throw std::runtime_error(
                "empty load interval in linearized AC seed: "
                + load.source_key);
        }
        if (movement_objective) {
            cost[dload_offset + i] = 1.0;
        }
    }
    if (economic_objective && !compact_economic_objective) {
        for (int i = 0; i < ng; ++i) {
            if (commitment[i] == 0) {
                lower[gen_cost_offset + i] = 0.0;
                upper[gen_cost_offset + i] = 0.0;
            } else {
                lower[gen_cost_offset + i] = -kHighsInf;
                upper[gen_cost_offset + i] = kHighsInf;
            }
            cost[gen_cost_offset + i] = data.delta;
        }
        for (int i = 0; i < nd; ++i) {
            lower[load_value_offset + i] = -kHighsInf;
            upper[load_value_offset + i] = kHighsInf;
            cost[load_value_offset + i] = -data.delta;
        }
        for (int i = 0; i < nl; ++i) {
            lower[sm_offset + i] = 0.0;
            upper[sm_offset + i] = data.branches[i].status == 0
                ? 0.0 : data.sm_vio_limit;
            cost[sm_offset + i] =
                data.delta * data.sm_cost_approx;
        }
    } else if (compact_economic_objective) {
        // A subgradient of a convex generator-cost curve is a global affine
        // under-estimator; a subgradient of a concave load-value curve is a
        // global affine over-estimator.  These coefficients therefore give a
        // compact first-order market-surplus direction without PWL auxiliary
        // columns.  The canonical nonlinear rebuild computes the exact source
        // PWL objective before any candidate can be accepted.
        for (int i = 0; i < ng; ++i) {
            if (commitment[i] == 0) {
                continue;
            }
            const auto& gen = data.generators[i];
            const auto points = active_pwl_points(
                gen.cost, gen.ncost,
                lower[pg_offset + i], upper[pg_offset + i]);
            cost[pg_offset + i] = data.delta *
                piecewise_linear_subgradient(points, reference.pg[i]);
        }
        for (int i = 0; i < nd; ++i) {
            const auto& load = data.loads[i];
            const auto points = active_pwl_points(
                load.cost, load.ncost, load.pd_min, load.pd_max);
            cost[demand_offset + i] = -data.delta * load.pd_nominal *
                piecewise_linear_subgradient(
                    points,
                    load.pd_nominal * reference.demand_factor[i]);
        }
    }
    if (compact_economic_objective) {
        // Express each movement as nonnegative up/down columns.  Zero
        // movement is then at a bound for every physical column; one of the
        // two balance-slack columns can be basic in each nodal equation,
        // yielding an explicit feasible simplex basis at the incumbent.
        const auto split_movement_bounds = [&] (
            int up_column,
            int down_column,
            double reference_value) {
            const double absolute_lower = lower[up_column];
            const double absolute_upper = upper[up_column];
            lower[up_column] = 0.0;
            upper[up_column] = std::max(
                0.0, absolute_upper - reference_value);
            lower[down_column] = 0.0;
            upper[down_column] = std::max(
                0.0, reference_value - absolute_lower);
            cost[down_column] = -cost[up_column];
        };
        for (int i = 0; i < nb; ++i) {
            split_movement_bounds(
                vm_offset + i, compact_vm_down_offset + i,
                linearization_vm[i]);
            split_movement_bounds(
                va_offset + i, compact_va_down_offset + i,
                reference.va[i]);
        }
        for (int i = 0; i < ng; ++i) {
            split_movement_bounds(
                pg_offset + i, compact_pg_down_offset + i,
                reference.pg[i]);
            split_movement_bounds(
                qg_offset + i, compact_qg_down_offset + i,
                reference.qg[i]);
        }
        for (int i = 0; i < nd; ++i) {
            split_movement_bounds(
                demand_offset + i, compact_demand_down_offset + i,
                reference.demand_factor[i]);
        }
    }

    std::vector<LinearizedBranch> linearized(nl);
    for (int i = 0; i < nl; ++i) {
        if (i == outaged_branch || data.branches[i].status == 0) {
            continue;
        }
        const auto& branch = data.branches[i];
        linearized[i] = linearize_branch(
            branch,
            linearization_vm[branch.from], linearization_vm[branch.to],
            reference.va[branch.from], reference.va[branch.to]);
    }

    std::vector<SparseRow> rows;
    rows.reserve(2 * nb + 5 * nl + 4 * ng + 2 * nd);
    for (int bus = 0; bus < nb; ++bus) {
        SparseRow active;
        SparseRow reactive;
        double active_constant = 0.0;
        double reactive_constant = 0.0;
        for (int branch_index : data.buses[bus].branches_from) {
            if (branch_index == outaged_branch ||
                data.branches[branch_index].status == 0) {
                continue;
            }
            const auto& branch = data.branches[branch_index];
            const auto& p = linearized[branch_index].pf;
            const auto& q = linearized[branch_index].qf;
            active_constant += p.constant;
            reactive_constant += q.constant;
            append(active, vm_offset + branch.from, p.vm_from);
            append(active, vm_offset + branch.to, p.vm_to);
            append(active, va_offset + branch.from, p.va_from);
            append(active, va_offset + branch.to, p.va_to);
            append(reactive, vm_offset + branch.from, q.vm_from);
            append(reactive, vm_offset + branch.to, q.vm_to);
            append(reactive, va_offset + branch.from, q.va_from);
            append(reactive, va_offset + branch.to, q.va_to);
        }
        for (int branch_index : data.buses[bus].branches_to) {
            if (branch_index == outaged_branch ||
                data.branches[branch_index].status == 0) {
                continue;
            }
            const auto& branch = data.branches[branch_index];
            const auto& p = linearized[branch_index].pt;
            const auto& q = linearized[branch_index].qt;
            active_constant += p.constant;
            reactive_constant += q.constant;
            append(active, vm_offset + branch.from, p.vm_from);
            append(active, vm_offset + branch.to, p.vm_to);
            append(active, va_offset + branch.from, p.va_from);
            append(active, va_offset + branch.to, p.va_to);
            append(reactive, vm_offset + branch.from, q.vm_from);
            append(reactive, vm_offset + branch.to, q.vm_to);
            append(reactive, va_offset + branch.from, q.va_from);
            append(reactive, va_offset + branch.to, q.va_to);
        }
        for (int generator : data.buses[bus].generators) {
            append(active, pg_offset + generator, -1.0);
            append(reactive, qg_offset + generator, -1.0);
        }
        for (int load : data.buses[bus].loads) {
            append(active, demand_offset + load, data.loads[load].pd_nominal);
            append(reactive, demand_offset + load, data.loads[load].qd_nominal);
        }
        double gs = 0.0;
        double bs = 0.0;
        for (int shunt : data.buses[bus].shunts) {
            gs += data.shunts[shunt].gs;
            // The nonlinear model and independent validator use the
            // step-dependent shunt state, not the immutable source starting
            // value.  Linearizing around data.shunts[shunt].bs after a
            // switched-shunt move makes the supposedly feasible incumbent
            // carry an artificial reactive-balance violation in this LP.
            bs += effective_shunt_susceptance(data, reference, shunt);
        }
        active_constant -=
            gs * linearization_vm[bus] * linearization_vm[bus];
        reactive_constant +=
            bs * linearization_vm[bus] * linearization_vm[bus];
        append(active, vm_offset + bus, 2.0 * gs * linearization_vm[bus]);
        append(reactive, vm_offset + bus, -2.0 * bs * linearization_vm[bus]);
        if (compact_economic_objective) {
            // Compute the affine row's reference residual from the canonical
            // rebuilt physical flows.  Reconstructing it as
            // constant + J*x suffers catastrophic cancellation on the 19k
            // near-zero-impedance circuits: terms around 1e5 previously left
            // a false residual above the exact 0.5 source slack bound even
            // though the independently verified residual was about 1e-5.
            double active_reference_residual = 0.0;
            double reactive_reference_residual = 0.0;
            for (int branch_index : data.buses[bus].branches_from) {
                if (branch_index == outaged_branch ||
                    data.branches[branch_index].status == 0) {
                    continue;
                }
                active_reference_residual += reference.pf[branch_index];
                reactive_reference_residual += reference.qf[branch_index];
            }
            for (int branch_index : data.buses[bus].branches_to) {
                if (branch_index == outaged_branch ||
                    data.branches[branch_index].status == 0) {
                    continue;
                }
                active_reference_residual += reference.pt[branch_index];
                reactive_reference_residual += reference.qt[branch_index];
            }
            for (int generator : data.buses[bus].generators) {
                active_reference_residual -= reference.pg[generator];
                reactive_reference_residual -= reference.qg[generator];
            }
            for (int load : data.buses[bus].loads) {
                active_reference_residual +=
                    data.loads[load].pd_nominal *
                    reference.demand_factor[load];
                reactive_reference_residual +=
                    data.loads[load].qd_nominal *
                    reference.demand_factor[load];
            }
            for (int shunt : data.buses[bus].shunts) {
                active_reference_residual +=
                    data.shunts[shunt].gs * linearization_vm[bus] *
                    linearization_vm[bus];
                reactive_reference_residual -=
                    effective_shunt_susceptance(data, reference, shunt) *
                    linearization_vm[bus] * linearization_vm[bus];
            }
            const auto down_column = [&] (HighsInt column) -> HighsInt {
                if (column >= vm_offset && column < va_offset) {
                    return compact_vm_down_offset + column - vm_offset;
                }
                if (column >= va_offset && column < pg_offset) {
                    return compact_va_down_offset + column - va_offset;
                }
                if (column >= pg_offset && column < qg_offset) {
                    return compact_pg_down_offset + column - pg_offset;
                }
                if (column >= qg_offset && column < demand_offset) {
                    return compact_qg_down_offset + column - qg_offset;
                }
                if (column >= demand_offset && column < p_pos_offset) {
                    return compact_demand_down_offset +
                        column - demand_offset;
                }
                throw std::runtime_error(
                    "unexpected compact economic movement column");
            };
            const auto active_up_entries = active.entries;
            const auto reactive_up_entries = reactive.entries;
            for (const auto& [column, coefficient] : active_up_entries) {
                append(active, down_column(column), -coefficient);
            }
            for (const auto& [column, coefficient] : reactive_up_entries) {
                append(reactive, down_column(column), -coefficient);
            }
            append(active, p_pos_offset + bus, -1.0);
            append(active, p_neg_offset + bus, 1.0);
            append(reactive, q_pos_offset + bus, -1.0);
            append(reactive, q_neg_offset + bus, 1.0);
            active.lower = active.upper = -active_reference_residual;
            reactive.lower = reactive.upper = -reactive_reference_residual;
        } else if (projected_balance_slack) {
            // Project the two bounded positive/negative balance-slack columns
            // onto the balance row.  This has the identical feasible set in
            // the physical seed variables but removes four active columns per
            // bus from the large contingency IPM system.
            active.lower = -active_constant - balance_slack_limit;
            active.upper = -active_constant + balance_slack_limit;
            reactive.lower = -reactive_constant - balance_slack_limit;
            reactive.upper = -reactive_constant + balance_slack_limit;
        } else if (!feasibility_only || elastic_balance_phase_one) {
            append(active, p_pos_offset + bus, -1.0);
            append(active, p_neg_offset + bus, 1.0);
            append(reactive, q_pos_offset + bus, -1.0);
            append(reactive, q_neg_offset + bus, 1.0);
            active.lower = active.upper = -active_constant;
            reactive.lower = reactive.upper = -reactive_constant;
        } else {
            // Feasibility-only mode has no explicit balance-slack columns.
            // Without projection, this is the same zero-slack equality that
            // the previous fixed-zero columns represented.
            active.lower = active.upper = -active_constant;
            reactive.lower = reactive.upper = -reactive_constant;
        }
        normalize(active);
        normalize(reactive);
        rows.push_back(std::move(active));
        rows.push_back(std::move(reactive));
    }

    if (economic_objective && !compact_economic_objective) {
        // The source generator curves are convex and the source load-value
        // curves are concave.  Epigraph/hypograph rows therefore represent
        // their complete PWL functions exactly without lambda or SOS2
        // columns.  The accepted nonlinear candidate is still rebuilt with
        // the canonical lambda semantics below this direction-generator LP.
        for (int i = 0; i < ng; ++i) {
            if (commitment[i] == 0) {
                continue;
            }
            const auto& gen = data.generators[i];
            const auto points = active_pwl_points(
                gen.cost, gen.ncost,
                lower[pg_offset + i], upper[pg_offset + i]);
            double previous_slope = -kHighsInf;
            for (std::size_t segment = 0;
                 segment + 1 < points.size(); ++segment) {
                const double width =
                    points[segment + 1].mw - points[segment].mw;
                if (std::abs(width) <= 1e-14) {
                    continue;
                }
                const double slope =
                    (points[segment + 1].cost - points[segment].cost) /
                    width;
                if (slope + 1e-9 < previous_slope) {
                    throw std::runtime_error(
                        "linearized economic generator curve is not convex: " +
                        gen.source_key);
                }
                previous_slope = slope;
                const double intercept =
                    points[segment].cost - slope * points[segment].mw;
                SparseRow epigraph;
                epigraph.lower = intercept;
                epigraph.upper = kHighsInf;
                append(epigraph, gen_cost_offset + i, 1.0);
                append(epigraph, pg_offset + i, -slope);
                rows.push_back(std::move(epigraph));
            }
        }
        for (int i = 0; i < nd; ++i) {
            const auto& load = data.loads[i];
            const auto points = active_pwl_points(
                load.cost, load.ncost, load.pd_min, load.pd_max);
            double previous_slope = kHighsInf;
            for (std::size_t segment = 0;
                 segment + 1 < points.size(); ++segment) {
                const double width =
                    points[segment + 1].mw - points[segment].mw;
                if (std::abs(width) <= 1e-14) {
                    continue;
                }
                const double slope =
                    (points[segment + 1].cost - points[segment].cost) /
                    width;
                if (slope > previous_slope + 1e-9) {
                    throw std::runtime_error(
                        "linearized economic load curve is not concave: " +
                        load.source_key);
                }
                previous_slope = slope;
                const double intercept =
                    points[segment].cost - slope * points[segment].mw;
                SparseRow hypograph;
                hypograph.lower = -kHighsInf;
                hypograph.upper = intercept;
                append(hypograph, load_value_offset + i, 1.0);
                append(
                    hypograph, demand_offset + i,
                    -slope * load.pd_nominal);
                rows.push_back(std::move(hypograph));
            }
        }
    }

    const auto affine_range = [&](const AffineFlow& flow, const Branch& branch) {
        double minimum = flow.constant;
        double maximum = flow.constant;
        const std::array<std::pair<int, double>, 4> terms{{
            {vm_offset + branch.from, flow.vm_from},
            {vm_offset + branch.to, flow.vm_to},
            {va_offset + branch.from, flow.va_from},
            {va_offset + branch.to, flow.va_to},
        }};
        for (const auto& [column, coefficient] : terms) {
            if (coefficient >= 0.0) {
                minimum += coefficient * lower[column];
                maximum += coefficient * upper[column];
            } else {
                minimum += coefficient * upper[column];
                maximum += coefficient * lower[column];
            }
        }
        return std::pair<double, double>{minimum, maximum};
    };

    for (int i = 0; i < nl; ++i) {
        if (omit_branch_security_rows &&
            !selected_branch_security[static_cast<std::size_t>(i)]) {
            continue;
        }
        if (i == outaged_branch || data.branches[i].status == 0) {
            continue;
        }
        const auto& branch = data.branches[i];
        const double rating = contingency ? branch.rate_c : branch.rate_a;
        const std::array<const AffineFlow*, 4> flows{
            &linearized[i].pf, &linearized[i].qf,
            &linearized[i].pt, &linearized[i].qt};
        const double from_component_bound = branch_terminal_component_bound(
            data, branch, rating, true);
        const double to_component_bound = branch_terminal_component_bound(
            data, branch, rating, false);
        const std::array<double, 4> component_bounds{
            from_component_bound, from_component_bound,
            to_component_bound, to_component_bound};
        for (std::size_t side = 0; side < flows.size(); ++side) {
            const auto* flow = flows[side];
            const double component_bound = component_bounds[side];
            if (lightweight_large_seed) {
                const auto [minimum, maximum] = affine_range(*flow, branch);
                if (minimum >= -component_bound &&
                    maximum <= component_bound) {
                    continue;
                }
            }
            SparseRow row;
            row.lower = -component_bound - flow->constant;
            row.upper = component_bound - flow->constant;
            append(row, vm_offset + branch.from, flow->vm_from);
            append(row, vm_offset + branch.to, flow->vm_to);
            append(row, va_offset + branch.from, flow->va_from);
            append(row, va_offset + branch.to, flow->va_to);
            normalize(row);
            rows.push_back(std::move(row));
        }
        if (economic_objective && !compact_economic_objective) {
            const auto base_flow = branch_flows(
                branch,
                linearization_vm[branch.from],
                linearization_vm[branch.to],
                reference.va[branch.from],
                reference.va[branch.to]);
            const double reference_sm =
                reference.sm_slack.size() == data.branches.size()
                ? std::clamp(
                    reference.sm_slack[i], 0.0, data.sm_vio_limit)
                : 0.0;
            const auto append_thermal_linearization = [&] (
                const AffineFlow& active_flow,
                const AffineFlow& reactive_flow,
                double active_base,
                double reactive_base,
                int terminal_bus) {
                const double terminal_voltage = branch.transformer
                    ? 1.0 : linearization_vm[terminal_bus];
                const double rating_term = terminal_voltage + reference_sm;
                const double rating_squared = rating * rating;
                const double rating_gradient =
                    -2.0 * rating_squared * rating_term;
                const double vm_from_gradient =
                    2.0 * active_base * active_flow.vm_from +
                    2.0 * reactive_base * reactive_flow.vm_from +
                    (!branch.transformer && terminal_bus == branch.from
                        ? rating_gradient : 0.0);
                const double vm_to_gradient =
                    2.0 * active_base * active_flow.vm_to +
                    2.0 * reactive_base * reactive_flow.vm_to +
                    (!branch.transformer && terminal_bus == branch.to
                        ? rating_gradient : 0.0);
                const double va_from_gradient =
                    2.0 * active_base * active_flow.va_from +
                    2.0 * reactive_base * reactive_flow.va_from;
                const double va_to_gradient =
                    2.0 * active_base * active_flow.va_to +
                    2.0 * reactive_base * reactive_flow.va_to;
                const double sm_gradient = rating_gradient;
                const double base_value =
                    active_base * active_base +
                    reactive_base * reactive_base -
                    rating_squared * rating_term * rating_term;
                const double constant = base_value
                    - vm_from_gradient * linearization_vm[branch.from]
                    - vm_to_gradient * linearization_vm[branch.to]
                    - va_from_gradient * reference.va[branch.from]
                    - va_to_gradient * reference.va[branch.to]
                    - sm_gradient * reference_sm;
                SparseRow thermal;
                thermal.lower = -kHighsInf;
                thermal.upper = -constant;
                append(
                    thermal, vm_offset + branch.from,
                    vm_from_gradient);
                append(
                    thermal, vm_offset + branch.to,
                    vm_to_gradient);
                append(
                    thermal, va_offset + branch.from,
                    va_from_gradient);
                append(
                    thermal, va_offset + branch.to,
                    va_to_gradient);
                append(thermal, sm_offset + i, sm_gradient);
                rows.push_back(std::move(thermal));
            };
            append_thermal_linearization(
                linearized[i].pf, linearized[i].qf,
                base_flow.pf, base_flow.qf, branch.from);
            append_thermal_linearization(
                linearized[i].pt, linearized[i].qt,
                base_flow.pt, base_flow.qt, branch.to);
        }
        const double source_delta = contingency
            ? contingency->effective_base_state().va[branch.from]
                - contingency->effective_base_state().va[branch.to]
            : data.buses[branch.from].va_start
                - data.buses[branch.to].va_start;
        const double angle_minimum =
            lower[va_offset + branch.from] - upper[va_offset + branch.to];
        const double angle_maximum =
            upper[va_offset + branch.from] - lower[va_offset + branch.to];
        const bool angle_row_is_redundant =
            lightweight_large_seed &&
            angle_minimum >= branch.angmin && angle_maximum <= branch.angmax;
        if (source_delta >= branch.angmin && source_delta <= branch.angmax &&
            !angle_row_is_redundant) {
            SparseRow angle;
            angle.lower = branch.angmin;
            angle.upper = branch.angmax;
            append(angle, va_offset + branch.from, 1.0);
            append(angle, va_offset + branch.to, -1.0);
            rows.push_back(std::move(angle));
        }
    }

    // Absolute deviation variables and their epigraph rows define the
    // ordinary seed objective.  They have no mathematical role in the
    // zero-objective Phase-I problem and substantially enlarge its normal
    // equations at 19k buses, so omit the rows entirely in feasibility mode.
    // Their now-isolated zero-cost columns are removed by HiGHS presolve.
    if (movement_objective) {
        for (int i = 0; i < ng; ++i) {
            SparseRow pg_positive;
            pg_positive.lower = -kHighsInf;
            pg_positive.upper = reference.pg[i];
            append(pg_positive, pg_offset + i, 1.0);
            append(pg_positive, dpg_offset + i, -1.0);
            rows.push_back(std::move(pg_positive));
            SparseRow pg_negative;
            pg_negative.lower = -kHighsInf;
            pg_negative.upper = -reference.pg[i];
            append(pg_negative, pg_offset + i, -1.0);
            append(pg_negative, dpg_offset + i, -1.0);
            rows.push_back(std::move(pg_negative));
            SparseRow qg_positive;
            qg_positive.lower = -kHighsInf;
            qg_positive.upper = reference.qg[i];
            append(qg_positive, qg_offset + i, 1.0);
            append(qg_positive, dqg_offset + i, -1.0);
            rows.push_back(std::move(qg_positive));
            SparseRow qg_negative;
            qg_negative.lower = -kHighsInf;
            qg_negative.upper = -reference.qg[i];
            append(qg_negative, qg_offset + i, -1.0);
            append(qg_negative, dqg_offset + i, -1.0);
            rows.push_back(std::move(qg_negative));
        }
        for (int i = 0; i < nd; ++i) {
            const double coefficient = data.loads[i].pd_nominal;
            const double reference_power =
                coefficient * reference.demand_factor[i];
            SparseRow positive;
            positive.lower = -kHighsInf;
            positive.upper = reference_power;
            append(positive, demand_offset + i, coefficient);
            append(positive, dload_offset + i, -1.0);
            rows.push_back(std::move(positive));
            SparseRow negative;
            negative.lower = -kHighsInf;
            negative.upper = -reference_power;
            append(negative, demand_offset + i, -coefficient);
            append(negative, dload_offset + i, -1.0);
            rows.push_back(std::move(negative));
        }
    }

    std::vector<double> column_scale(
        static_cast<std::size_t>(column_count), 1.0);
    if (elastic_balance_phase_one) {
        // The 19k cases contain near-zero-impedance circuits whose angle
        // derivatives reach roughly 2e5.  Long primal-simplex pivot sequences
        // on the unscaled columns eventually violated nonnegative elastic
        // bounds.  Scale only the zero-cost physical Phase-I columns, capped
        // at 1e4 so small but meaningful coefficients remain above HiGHS'
        // minimum accepted matrix value.  This is an exact variable change:
        // y_j = scale_j * x_j.
        constexpr double kMaximumPhysicalColumnScale = 1e4;
        const auto scalable_column = [&] (HighsInt column) {
            return column < p_pos_offset;
        };
        for (const auto& row : rows) {
            for (const auto& [column, coefficient] : row.entries) {
                if (!scalable_column(column)) {
                    continue;
                }
                column_scale[column] = std::min(
                    kMaximumPhysicalColumnScale,
                    std::max(column_scale[column], std::abs(coefficient)));
            }
        }
        for (int column = 0; column < column_count; ++column) {
            if (!scalable_column(column)) {
                continue;
            }
            lower[column] *= column_scale[column];
            upper[column] *= column_scale[column];
            cost[column] /= column_scale[column];
        }
        for (auto& row : rows) {
            for (auto& [column, coefficient] : row.entries) {
                if (scalable_column(column)) {
                    coefficient /= column_scale[column];
                }
            }
        }
    }

    double maximum_row_scale = 1.0;
    double objective_scale = 1.0;
    if (compact_economic_objective) {
        // This is an exact positive scaling of each equality, not a model
        // relaxation.  Near-zero-impedance circuits create derivatives near
        // 2e5 in the raw nodal rows; limiting the external row scale to 1e4
        // keeps the explicit balance-slack pivots numerically visible while
        // reducing the large unscaled residuals observed at a time limit.
        constexpr double kMaximumCompactRowScale = 1e4;
        for (auto& row : rows) {
            double row_scale = 1.0;
            for (const auto& [column, coefficient] : row.entries) {
                static_cast<void>(column);
                row_scale = std::max(row_scale, std::abs(coefficient));
            }
            row_scale = std::min(row_scale, kMaximumCompactRowScale);
            maximum_row_scale = std::max(maximum_row_scale, row_scale);
            if (row_scale <= 1.0) {
                continue;
            }
            for (auto& [column, coefficient] : row.entries) {
                static_cast<void>(column);
                coefficient /= row_scale;
            }
            if (std::isfinite(row.lower)) {
                row.lower /= row_scale;
            }
            if (std::isfinite(row.upper)) {
                row.upper /= row_scale;
            }
        }
        for (double coefficient : cost) {
            objective_scale = std::max(
                objective_scale, std::abs(coefficient));
        }
        if (objective_scale > 1.0) {
            for (double& coefficient : cost) {
                coefficient /= objective_scale;
            }
        }
    }

    const char* highs_log = std::getenv("GRAVITYX_HIGHS_LOG");
    // The 8k-bus trust-region seed has a much smaller, better-scaled presolved
    // system than the original full contingency LP.  IPM solves this form in
    // less than half the measured dual-simplex time and returns a candidate
    // that the nonlinear correction can validate in one round.
    const char* solver_override = std::getenv("GRAVITYX_LINEAR_SEED_SOLVER");
    const std::string solver = solver_override != nullptr
        ? std::string(solver_override)
        : ((feasibility_only && nb >= 16000) || economic_objective
            ? "simplex" : "ipm");
    int simplex_strategy = elastic_balance_phase_one ? 1 : 4;
    const char* simplex_strategy_override =
        std::getenv("GRAVITYX_LINEAR_SEED_SIMPLEX_STRATEGY");
    if (simplex_strategy_override != nullptr) {
        simplex_strategy = std::stoi(simplex_strategy_override);
    }
    // The feasibility-only LP has an identically zero objective, so every
    // primal-feasible point is already an exact optimum for Phase I.  A looser
    // IPM dual-gap stopping test avoids spending tens of seconds proving an
    // economically meaningless zero-objective dual certificate.  The returned
    // primal must still satisfy the unchanged 1e-8 HiGHS feasibility test and
    // the exact nonlinear validator below this layer.
    constexpr double kOrdinaryIpmOptimalityTolerance = 1e-8;
    constexpr double kPhaseOneIpmOptimalityTolerance = 1e-4;
    const double ipm_optimality_tolerance = feasibility_only
        ? kPhaseOneIpmOptimalityTolerance
        : kOrdinaryIpmOptimalityTolerance;
    // Match sparse preprocessing to the exact value HiGHS will use.  HiGHS
    // returns kWarning after dropping smaller entries.  The previous wrapper
    // treated that warning as a construction failure even though the model
    // was valid and had been loaded.  Pruning the same entries ourselves is
    // mathematically identical to the solver's effective matrix and makes the
    // transformation explicit and auditable.
    constexpr double kOrdinarySmallMatrixValue = 1e-9;
    constexpr double kElasticSmallMatrixValue = 1e-12;
    const double small_matrix_value =
        (elastic_balance_phase_one || compact_economic_objective)
        ? kElasticSmallMatrixValue : kOrdinarySmallMatrixValue;
    const auto model_audit = audit_and_prune_model(
        lower, upper, cost, rows, small_matrix_value);

    std::vector<double> row_lower;
    std::vector<double> row_upper;
    std::vector<HighsInt> starts;
    std::vector<HighsInt> indices;
    std::vector<double> values;
    row_lower.reserve(rows.size());
    row_upper.reserve(rows.size());
    starts.reserve(rows.size() + 1);
    starts.push_back(0);
    if (model_audit.passed) {
        for (const auto& row : rows) {
            row_lower.push_back(row.lower);
            row_upper.push_back(row.upper);
            for (const auto& [column, value] : row.entries) {
                indices.push_back(column);
                values.push_back(value);
            }
            starts.push_back(static_cast<HighsInt>(indices.size()));
        }
    }

    LinearizedAcSeedResult output;
    output.economic_objective = economic_objective;
    output.compact_economic_objective = compact_economic_objective;
    output.projected_balance_slack = projected_balance_slack;
    output.branch_security_rows_omitted = omit_branch_security_rows;
    output.branch_security_subset_count = branch_security_subset_count;
    output.feasibility_only = feasibility_only;
    output.elastic_balance_phase_one = elastic_balance_phase_one;
    output.simplex_strategy = simplex_strategy;
    output.maximum_column_scale = *std::max_element(
        column_scale.begin(), column_scale.end());
    output.maximum_row_scale = maximum_row_scale;
    output.objective_scale = objective_scale;
    output.voltage_trust_radius = lightweight_large_seed
        ? voltage_trust_radius : 0.0;
    output.angle_trust_radius = lightweight_large_seed
        ? angle_trust_radius : 0.0;
    output.projected_reference_voltage_count =
        projected_reference_voltage_count;
    output.maximum_reference_voltage_projection =
        maximum_reference_voltage_projection;
    output.time_limit_seconds = time_limit_seconds;
    output.ipm_optimality_tolerance = ipm_optimality_tolerance;
    output.row_count = static_cast<int>(rows.size());
    output.column_count = column_count;
    output.nonzero_count = static_cast<int>(indices.size());
    output.model_preflight_passed = model_audit.passed;
    output.model_preflight_failure = model_audit.failure;
    output.tiny_matrix_entries_removed =
        model_audit.tiny_matrix_entries_removed;
    output.maximum_tiny_matrix_entry_removed =
        model_audit.maximum_tiny_matrix_entry_removed;
    output.small_matrix_value = small_matrix_value;
    output.run_status = -99;
    output.model_status = -99;
    output.primal_solution_status = -99;
    const auto finish_model_failure = [&](const std::string& call,
                                          const std::string& status) {
        output.success = false;
        output.model_construction_success = false;
        output.model_load_failure_call = call;
        output.status = status;
        output.state = reference;
        output.wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        if (highs_log != nullptr && std::string(highs_log) != "0") {
            std::cerr << "GRAVITYX_LINEAR_SEED_RESULT "
                      << output.to_json(false).dump() << '\n';
        }
        return output;
    };
    if (!model_audit.passed) {
        return finish_model_failure(
            "preflight", "Model preflight failed: " + model_audit.failure);
    }

    Highs highs;
    highs.setOptionValue(
        "output_flag", highs_log != nullptr && std::string(highs_log) != "0");
    highs.setOptionValue("threads", 1);
    // The compact economic model has an exact incumbent primal and an
    // explicit diagonal balance-slack basis.  Preserve that basis in the
    // original column space: presolve otherwise remaps or discards its basic
    // elastic columns before the deliberately short resident simplex solve.
    output.presolve_enabled = !compact_economic_objective;
    highs.setOptionValue(
        "presolve", output.presolve_enabled ? "on" : "off");
    highs.setOptionValue("small_matrix_value", small_matrix_value);
    highs.setOptionValue("solver", solver);
    if (solver == "simplex" && (feasibility_only || economic_objective)) {
        highs.setOptionValue("simplex_strategy", simplex_strategy);
        if (compact_economic_objective) {
            // A time-limited primal-simplex return must remain within the
            // exact source bounds.  Bound perturbation had produced terminal
            // movement columns outside those bounds even though the live log
            // reported a feasible iterate.  Disabling it changes no model
            // bound or acceptance tolerance.
            output.primal_simplex_bound_perturbation_multiplier = 0.0;
            highs.setOptionValue(
                "primal_simplex_bound_perturbation_multiplier", 0.0);
        }
    }
    highs.setOptionValue("run_crossover", "off");
    highs.setOptionValue("time_limit", time_limit_seconds);
    highs.setOptionValue("primal_feasibility_tolerance", 1e-8);
    highs.setOptionValue("dual_feasibility_tolerance", 1e-8);
    highs.setOptionValue(
        "ipm_optimality_tolerance", ipm_optimality_tolerance);

    const HighsStatus add_vars_status =
        highs.addVars(column_count, lower.data(), upper.data());
    output.add_vars_status = static_cast<int>(add_vars_status);
    output.model_load_warning = add_vars_status == HighsStatus::kWarning;
    if (add_vars_status == HighsStatus::kError) {
        return finish_model_failure(
            "addVars", "HiGHS addVars returned kError");
    }
    const HighsStatus change_cols_cost_status =
        highs.changeColsCost(0, column_count - 1, cost.data());
    output.change_cols_cost_status =
        static_cast<int>(change_cols_cost_status);
    output.model_load_warning = output.model_load_warning ||
        change_cols_cost_status == HighsStatus::kWarning;
    if (change_cols_cost_status == HighsStatus::kError) {
        return finish_model_failure(
            "changeColsCost", "HiGHS changeColsCost returned kError");
    }
    const HighsStatus add_rows_status = highs.addRows(
        static_cast<HighsInt>(rows.size()),
        row_lower.data(), row_upper.data(),
        static_cast<HighsInt>(indices.size()), starts.data(),
        indices.data(), values.data());
    output.add_rows_status = static_cast<int>(add_rows_status);
    output.model_load_warning = output.model_load_warning ||
        add_rows_status == HighsStatus::kWarning;
    if (add_rows_status == HighsStatus::kError) {
        return finish_model_failure(
            "addRows", "HiGHS addRows returned kError");
    }
    output.model_construction_success = true;
    bool primal_start_attempted = false;
    HighsStatus primal_start_status = HighsStatus::kOk;
    bool primal_basis_attempted = false;
    HighsStatus primal_basis_status = HighsStatus::kOk;
    const char* elastic_start_override =
        std::getenv("GRAVITYX_ELASTIC_PHASE_ONE_START");
    const bool elastic_start_enabled =
        elastic_start_override != nullptr &&
        std::string(elastic_start_override) != "0";
    if (elastic_balance_phase_one && elastic_start_enabled) {
        primal_start_attempted = true;
        std::vector<double> primal_start(
            static_cast<std::size_t>(column_count), 0.0);
        HighsBasis primal_basis;
        primal_basis.alien = true;
        primal_basis.useful = true;
        primal_basis.col_status.assign(
            static_cast<std::size_t>(column_count),
            HighsBasisStatus::kLower);
        primal_basis.row_status.assign(
            rows.size(), HighsBasisStatus::kBasic);
        const auto initialize_nonbasic_column = [&](int column) {
            if (std::isfinite(lower[column])) {
                primal_start[column] = lower[column];
                primal_basis.col_status[column] = HighsBasisStatus::kLower;
            } else if (std::isfinite(upper[column])) {
                primal_start[column] = upper[column];
                primal_basis.col_status[column] = HighsBasisStatus::kUpper;
            } else {
                primal_start[column] = 0.0;
                primal_basis.col_status[column] = HighsBasisStatus::kZero;
            }
        };
        for (int column = 0; column < column_count; ++column) {
            initialize_nonbasic_column(column);
        }
        const auto set_balance_slack = [&](int row_index,
                                           int positive_column,
                                           int negative_column) {
            double activity = 0.0;
            for (const auto& [column, coefficient] : rows[row_index].entries) {
                activity += coefficient * primal_start[column];
            }
            const double residual = activity - rows[row_index].lower;
            if (residual >= 0.0) {
                // The positive column has coefficient -1.
                primal_start[positive_column] = residual;
                primal_basis.col_status[positive_column] =
                    HighsBasisStatus::kBasic;
            } else {
                // The negative column has coefficient +1.
                primal_start[negative_column] = -residual;
                primal_basis.col_status[negative_column] =
                    HighsBasisStatus::kBasic;
            }
            // A basic elastic column replaces the implicit row slack.
            primal_basis.row_status[row_index] = HighsBasisStatus::kLower;
        };
        for (int bus = 0; bus < nb; ++bus) {
            set_balance_slack(
                2 * bus, p_pos_offset + bus, p_neg_offset + bus);
            set_balance_slack(
                2 * bus + 1, q_pos_offset + bus, q_neg_offset + bus);
        }
        std::vector<HighsInt> start_indices(
            static_cast<std::size_t>(column_count));
        std::iota(start_indices.begin(), start_indices.end(), HighsInt{0});
        primal_start_status = highs.setSolution(
            static_cast<HighsInt>(column_count),
            start_indices.data(), primal_start.data());
        primal_basis_attempted = true;
        primal_basis_status = highs.setBasis(
            primal_basis, "elastic_balance_diagonal");
    } else if (economic_objective) {
        // The verified nonlinear incumbent is also feasible for its own
        // first-order model.  Supply it explicitly so primal simplex starts
        // from a useful market dispatch rather than rebuilding one cold.
        primal_start_attempted = true;
        std::vector<double> primal_start(
            static_cast<std::size_t>(column_count), 0.0);
        for (int bus = 0; bus < nb; ++bus) {
            primal_start[vm_offset + bus] = compact_economic_objective
                ? 0.0 : linearization_vm[bus];
            primal_start[va_offset + bus] = compact_economic_objective
                ? 0.0 : reference.va[bus];
        }
        for (int i = 0; i < ng; ++i) {
            primal_start[pg_offset + i] = compact_economic_objective
                ? 0.0 : reference.pg[i];
            primal_start[qg_offset + i] = compact_economic_objective
                ? 0.0 : reference.qg[i];
            if (!compact_economic_objective && commitment[i] != 0) {
                const auto& gen = data.generators[i];
                const auto points = active_pwl_points(
                    gen.cost, gen.ncost,
                    lower[pg_offset + i], upper[pg_offset + i]);
                primal_start[gen_cost_offset + i] =
                    piecewise_linear_value(points, reference.pg[i]);
            }
        }
        for (int i = 0; i < nd; ++i) {
            primal_start[demand_offset + i] = compact_economic_objective
                ? 0.0 : reference.demand_factor[i];
            if (!compact_economic_objective) {
                const auto& load = data.loads[i];
                const auto points = active_pwl_points(
                    load.cost, load.ncost, load.pd_min, load.pd_max);
                primal_start[load_value_offset + i] = piecewise_linear_value(
                    points, load.pd_nominal * reference.demand_factor[i]);
            }
        }
        if (!compact_economic_objective) {
            for (int i = 0; i < nl; ++i) {
                primal_start[sm_offset + i] =
                    reference.sm_slack.size() == data.branches.size()
                    ? std::clamp(
                        reference.sm_slack[i], 0.0, data.sm_vio_limit)
                    : 0.0;
            }
        }
        std::vector<int> balance_basic_columns(
            static_cast<std::size_t>(2 * nb), -1);
        const auto set_balance_slack = [&] (
            int row_index, int positive_column, int negative_column) {
            double activity = 0.0;
            double positive_coefficient = 0.0;
            double negative_coefficient = 0.0;
            for (const auto& [column, coefficient] : rows[row_index].entries) {
                activity += coefficient * primal_start[column];
                if (column == positive_column) {
                    positive_coefficient = coefficient;
                } else if (column == negative_column) {
                    negative_coefficient = coefficient;
                }
            }
            const double residual = activity - rows[row_index].lower;
            if (residual >= 0.0) {
                if (positive_coefficient >= 0.0) {
                    throw std::runtime_error(
                        "economic positive balance slack has invalid coefficient");
                }
                primal_start[positive_column] =
                    residual / -positive_coefficient;
                balance_basic_columns[row_index] = positive_column;
            } else {
                if (negative_coefficient <= 0.0) {
                    throw std::runtime_error(
                        "economic negative balance slack has invalid coefficient");
                }
                primal_start[negative_column] =
                    -residual / negative_coefficient;
                balance_basic_columns[row_index] = negative_column;
            }
        };
        for (int bus = 0; bus < nb; ++bus) {
            set_balance_slack(
                2 * bus, p_pos_offset + bus, p_neg_offset + bus);
            set_balance_slack(
                2 * bus + 1, q_pos_offset + bus, q_neg_offset + bus);
        }
        std::vector<HighsInt> start_indices(
            static_cast<std::size_t>(column_count));
        std::iota(start_indices.begin(), start_indices.end(), HighsInt{0});
        primal_start_status = highs.setSolution(
            static_cast<HighsInt>(column_count),
            start_indices.data(), primal_start.data());
        if (compact_economic_objective) {
            primal_basis_attempted = true;
            HighsBasis primal_basis;
            primal_basis.alien = true;
            primal_basis.useful = true;
            primal_basis.col_status.assign(
                static_cast<std::size_t>(column_count),
                HighsBasisStatus::kZero);
            primal_basis.row_status.assign(
                rows.size(), HighsBasisStatus::kLower);
            for (int column = 0; column < column_count; ++column) {
                if (std::abs(lower[column]) <= 1e-12) {
                    primal_basis.col_status[column] =
                        HighsBasisStatus::kLower;
                } else if (std::abs(upper[column]) <= 1e-12) {
                    primal_basis.col_status[column] =
                        HighsBasisStatus::kUpper;
                } else if (!(lower[column] < 0.0 && upper[column] > 0.0)) {
                    throw std::runtime_error(
                        "compact economic zero movement is outside a column bound");
                }
            }
            for (int row = 0; row < 2 * nb; ++row) {
                if (balance_basic_columns[row] < 0) {
                    throw std::runtime_error(
                        "compact economic balance basis is incomplete");
                }
                primal_basis.col_status[balance_basic_columns[row]] =
                    HighsBasisStatus::kBasic;
            }
            primal_basis_status = highs.setBasis(
                primal_basis, "compact_economic_balance_diagonal");
        }
    }
    struct ResidentSolveSnapshot {
        HighsStatus run_status{HighsStatus::kError};
        HighsModelStatus model_status{HighsModelStatus::kNotset};
        HighsSolution solution;
        HighsInfo info;
        bool canonicalized{};
        double canonical_max_primal_infeasibility{kHighsInf};
        double canonical_objective{kHighsInf};
    };
    ResidentSolveSnapshot terminal;
    ResidentSolveSnapshot selected_snapshot;
    bool selected_snapshot_valid = false;
    double selected_snapshot_objective = kHighsInf;
    constexpr double kEconomicDirectionSnapshotTolerance = 1e-4;
    const auto canonicalize_compact_direction = [&] (
        HighsSolution& candidate,
        double& maximum_primal_infeasibility,
        double& canonical_objective) {
        if (!candidate.value_valid ||
            candidate.col_value.size() !=
                static_cast<std::size_t>(column_count)) {
            return false;
        }
        for (int column = 0; column < column_count; ++column) {
            double& value = candidate.col_value[column];
            if (!std::isfinite(value)) {
                return false;
            }
            value = std::clamp(value, lower[column], upper[column]);
        }
        // These compact models contain a dedicated positive/negative slack
        // pair for every nodal equality.  Recompute those columns exactly
        // from the returned physical direction.  This corrects only the
        // unscaled simplex-return artifact; it neither changes a physical
        // control nor relaxes any model or independent validation tolerance.
        for (int bus = 0; bus < nb; ++bus) {
            for (int balance = 0; balance < 2; ++balance) {
                const int row_index = 2 * bus + balance;
                if (row_index >= static_cast<int>(rows.size())) {
                    return false;
                }
                const int positive_column = balance == 0
                    ? p_pos_offset + bus : q_pos_offset + bus;
                const int negative_column = balance == 0
                    ? p_neg_offset + bus : q_neg_offset + bus;
                candidate.col_value[positive_column] = 0.0;
                candidate.col_value[negative_column] = 0.0;
                double activity = 0.0;
                for (const auto& [column, coefficient] :
                     rows[row_index].entries) {
                    activity += coefficient * candidate.col_value[column];
                }
                const double required_change =
                    rows[row_index].lower - activity;
                double positive_coefficient = 0.0;
                double negative_coefficient = 0.0;
                for (const auto& [column, coefficient] :
                     rows[row_index].entries) {
                    if (column == positive_column) {
                        positive_coefficient = coefficient;
                    } else if (column == negative_column) {
                        negative_coefficient = coefficient;
                    }
                }
                if (required_change <= 0.0) {
                    if (positive_coefficient >= 0.0) {
                        return false;
                    }
                    const double value =
                        required_change / positive_coefficient;
                    if (value > upper[positive_column] + 1e-12) {
                        return false;
                    }
                    candidate.col_value[positive_column] = value;
                } else {
                    if (negative_coefficient <= 0.0) {
                        return false;
                    }
                    const double value =
                        required_change / negative_coefficient;
                    if (value > upper[negative_column] + 1e-12) {
                        return false;
                    }
                    candidate.col_value[negative_column] = value;
                }
            }
        }
        maximum_primal_infeasibility = 0.0;
        for (int column = 0; column < column_count; ++column) {
            maximum_primal_infeasibility = std::max({
                maximum_primal_infeasibility,
                lower[column] - candidate.col_value[column],
                candidate.col_value[column] - upper[column],
            });
        }
        for (std::size_t row_index = 0;
             row_index < rows.size(); ++row_index) {
            double activity = 0.0;
            for (const auto& [column, coefficient] :
                 rows[row_index].entries) {
                activity += coefficient * candidate.col_value[column];
            }
            maximum_primal_infeasibility = std::max({
                maximum_primal_infeasibility,
                rows[row_index].lower - activity,
                activity - rows[row_index].upper,
            });
        }
        canonical_objective = 0.0;
        for (int column = 0; column < column_count; ++column) {
            canonical_objective +=
                cost[column] * candidate.col_value[column];
        }
        return std::isfinite(maximum_primal_infeasibility) &&
            maximum_primal_infeasibility <= 1e-8 &&
            std::isfinite(canonical_objective);
    };
    const bool segmented_resident_simplex =
        compact_economic_objective && solver == "simplex";
    const int maximum_segments = segmented_resident_simplex
        ? std::max(1, static_cast<int>(
            std::ceil(time_limit_seconds / 10.0)) + 1)
        : 1;
    int prior_simplex_iterations = -1;
    for (int segment = 1; segment <= maximum_segments; ++segment) {
        if (segmented_resident_simplex) {
            const double total_elapsed = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - wall_start).count();
            const double remaining = time_limit_seconds - total_elapsed;
            if (remaining <= 0.05) {
                break;
            }
            // HiGHS' run timer is resident and cumulative.  Advancing the
            // deadline in short increments makes it return a valid feasible
            // simplex iterate before a later long pivot sequence can fail
            // numerically.  The model and basis remain resident between
            // calls, so this is continuation, not a cold restart.
            const double segment_seconds = std::min(10.0, remaining);
            highs.setOptionValue(
                "time_limit", highs.getRunTime() + segment_seconds);
        }
        const auto segment_wall_start = std::chrono::steady_clock::now();
        terminal.run_status = highs.run();
        terminal.model_status = highs.getModelStatus();
        terminal.solution = highs.getSolution();
        terminal.info = highs.getInfo();
        const double segment_wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - segment_wall_start).count();
        ++output.resident_segment_count;

        const bool shape_valid = terminal.solution.value_valid &&
            terminal.solution.col_value.size() ==
                static_cast<std::size_t>(column_count);
        const bool solver_direction_snapshot_valid =
            terminal.run_status != HighsStatus::kError &&
            shape_valid && terminal.info.valid &&
            terminal.info.num_primal_infeasibilities >= 0 &&
            std::isfinite(terminal.info.max_primal_infeasibility) &&
            terminal.info.max_primal_infeasibility <=
                kEconomicDirectionSnapshotTolerance &&
            std::isfinite(terminal.info.objective_function_value);
        if (compact_economic_objective && shape_valid) {
            terminal.canonicalized = canonicalize_compact_direction(
                terminal.solution,
                terminal.canonical_max_primal_infeasibility,
                terminal.canonical_objective);
        }
        const bool direction_snapshot_valid =
            solver_direction_snapshot_valid || terminal.canonicalized;
        const double snapshot_objective = terminal.canonicalized
            ? terminal.canonical_objective
            : terminal.info.objective_function_value;
        if (direction_snapshot_valid) {
            ++output.feasible_segment_snapshot_count;
            if (!selected_snapshot_valid ||
                snapshot_objective < selected_snapshot_objective) {
                selected_snapshot = terminal;
                selected_snapshot_valid = true;
                selected_snapshot_objective = snapshot_objective;
            }
        }
        output.resident_segments.push_back({
            {"segment", segment},
            {"wall_seconds", segment_wall_seconds},
            {"run_status", static_cast<int>(terminal.run_status)},
            {"model_status", static_cast<int>(terminal.model_status)},
            {"solution_value_valid", terminal.solution.value_valid},
            {"info_valid", terminal.info.valid},
            {"simplex_iterations",
             static_cast<int>(terminal.info.simplex_iteration_count)},
            {"objective", terminal.info.objective_function_value},
            {"num_primal_infeasibilities",
             static_cast<int>(terminal.info.num_primal_infeasibilities)},
            {"max_primal_infeasibility",
             terminal.info.max_primal_infeasibility},
            {"canonicalized", terminal.canonicalized},
            {"canonical_max_primal_infeasibility",
             terminal.canonical_max_primal_infeasibility},
            {"canonical_objective", terminal.canonical_objective},
            {"snapshot_accepted", direction_snapshot_valid},
        });

        if (!segmented_resident_simplex ||
            terminal.model_status == HighsModelStatus::kOptimal ||
            terminal.run_status == HighsStatus::kError) {
            break;
        }
        const int current_simplex_iterations = static_cast<int>(
            terminal.info.simplex_iteration_count);
        if (segment > 1 && segment_wall_seconds < 0.05 &&
            current_simplex_iterations <= prior_simplex_iterations) {
            break;
        }
        prior_simplex_iterations = current_simplex_iterations;
    }
    output.terminal_run_status = static_cast<int>(terminal.run_status);
    output.terminal_model_status = static_cast<int>(terminal.model_status);
    if (segmented_resident_simplex && selected_snapshot_valid) {
        output.recovered_feasible_segment_snapshot =
            terminal.run_status == HighsStatus::kError ||
            !terminal.solution.value_valid || !terminal.info.valid ||
            (terminal.canonicalized
                ? terminal.canonical_objective
                : terminal.info.objective_function_value) !=
                selected_snapshot_objective;
        terminal = std::move(selected_snapshot);
    }
    output.canonicalized_segment_snapshot = terminal.canonicalized;
    output.canonicalized_snapshot_max_primal_infeasibility =
        terminal.canonical_max_primal_infeasibility;
    output.canonicalized_snapshot_objective =
        terminal.canonical_objective;
    const auto run_status = terminal.run_status;
    const auto model_status = terminal.model_status;
    const auto& solution = terminal.solution;
    const auto& info = terminal.info;

    output.projected_balance_slack = projected_balance_slack;
    output.branch_security_rows_omitted = omit_branch_security_rows;
    output.branch_security_subset_count = branch_security_subset_count;
    output.feasibility_only = feasibility_only;
    output.elastic_balance_phase_one = elastic_balance_phase_one;
    output.primal_start_attempted = primal_start_attempted;
    output.primal_start_status = static_cast<int>(primal_start_status);
    output.primal_start_accepted =
        primal_start_attempted && primal_start_status == HighsStatus::kOk;
    output.primal_basis_attempted = primal_basis_attempted;
    output.primal_basis_status = static_cast<int>(primal_basis_status);
    output.primal_basis_accepted =
        primal_basis_attempted && primal_basis_status == HighsStatus::kOk;
    output.simplex_strategy = simplex_strategy;
    output.maximum_column_scale = *std::max_element(
        column_scale.begin(), column_scale.end());
    output.maximum_row_scale = maximum_row_scale;
    output.objective_scale = objective_scale;
    output.voltage_trust_radius = lightweight_large_seed
        ? voltage_trust_radius : 0.0;
    output.angle_trust_radius = lightweight_large_seed
        ? angle_trust_radius : 0.0;
    output.projected_reference_voltage_count =
        projected_reference_voltage_count;
    output.maximum_reference_voltage_projection =
        maximum_reference_voltage_projection;
    output.time_limit_seconds = time_limit_seconds;
    output.ipm_optimality_tolerance = ipm_optimality_tolerance;
    output.row_count = static_cast<int>(rows.size());
    output.column_count = column_count;
    output.nonzero_count = static_cast<int>(indices.size());
    output.run_status = static_cast<int>(run_status);
    output.model_status = static_cast<int>(model_status);
    output.solution_value_valid = solution.value_valid;
    output.info_valid = info.valid;
    output.primal_solution_status =
        static_cast<int>(info.primal_solution_status);
    output.num_primal_infeasibilities =
        static_cast<int>(info.num_primal_infeasibilities);
    output.max_primal_infeasibility = info.max_primal_infeasibility;
    output.status = highs.modelStatusToString(model_status);
    output.iterations = static_cast<int>(
        info.ipm_iteration_count > 0
            ? info.ipm_iteration_count
            : info.simplex_iteration_count);
    output.objective = terminal.canonicalized
        ? terminal.canonical_objective
        : info.objective_function_value;
    constexpr double kRequiredPrimalFeasibilityTolerance = 1e-8;
    const bool solution_shape_valid = solution.value_valid &&
        solution.col_value.size() == static_cast<std::size_t>(column_count);
    if ((elastic_balance_phase_one || economic_objective) &&
        solution_shape_valid) {
        for (int bus = 0; bus < nb; ++bus) {
            const double active_slack =
                solution.col_value[p_pos_offset + bus] +
                solution.col_value[p_neg_offset + bus];
            const double reactive_slack =
                solution.col_value[q_pos_offset + bus] +
                solution.col_value[q_neg_offset + bus];
            output.maximum_balance_slack = std::max({
                output.maximum_balance_slack,
                active_slack,
                reactive_slack,
            });
            output.total_balance_slack += active_slack + reactive_slack;
        }
    }
    const bool optimal_solution = run_status != HighsStatus::kError &&
        model_status == HighsModelStatus::kOptimal && solution_shape_valid;
    // Phase I is used only to generate a candidate for the exact nonlinear
    // repair and complete validator.  A zero-objective IPM may hit its time
    // limit after finding a fully feasible primal but before certifying dual
    // optimality.  Preserve that useful point only when HiGHS itself reports
    // no primal violations at the unchanged 1e-8 LP tolerance.  Ordinary cost
    // LPs continue to require an optimal model status.
    const bool feasible_nonoptimal_phase_one = feasibility_only &&
        run_status != HighsStatus::kError && solution_shape_valid && info.valid &&
        info.primal_solution_status == kSolutionStatusFeasible &&
        info.num_primal_infeasibilities == 0 &&
        std::isfinite(info.max_primal_infeasibility) &&
        info.max_primal_infeasibility <= kRequiredPrimalFeasibilityTolerance;
    output.accepted_feasible_nonoptimal_phase_one =
        feasible_nonoptimal_phase_one && !optimal_solution;
    const bool feasible_nonoptimal_economic = economic_objective &&
        run_status != HighsStatus::kError && solution_shape_valid && info.valid &&
        info.primal_solution_status == kSolutionStatusFeasible &&
        info.num_primal_infeasibilities == 0 &&
        std::isfinite(info.max_primal_infeasibility) &&
        info.max_primal_infeasibility <= kRequiredPrimalFeasibilityTolerance;
    output.accepted_feasible_nonoptimal_economic =
        feasible_nonoptimal_economic && !optimal_solution;
    // A time-limited primal-simplex iterate can be slightly infeasible after
    // HiGHS unscales it even when the live simplex path reports Pr: 0.  This
    // point is useful only as a search direction.  The caller must rebuild the
    // nonlinear state and pass the unchanged independent 1e-5 validation gate
    // before it can affect the incumbent.
    constexpr double kApproximateEconomicDirectionTolerance = 1e-4;
    const bool approximate_economic_direction = economic_objective &&
        solution_shape_valid &&
        (terminal.canonicalized ||
         (run_status != HighsStatus::kError && info.valid &&
          info.num_primal_infeasibilities >= 0 &&
          std::isfinite(info.max_primal_infeasibility) &&
          info.max_primal_infeasibility <=
              kApproximateEconomicDirectionTolerance));
    output.accepted_approximate_economic_direction =
        approximate_economic_direction && !optimal_solution &&
        !feasible_nonoptimal_economic;
    output.success = optimal_solution || feasible_nonoptimal_phase_one ||
        feasible_nonoptimal_economic || approximate_economic_direction;
    output.state = reference;
    if (output.success) {
        for (int i = 0; i < nb; ++i) {
            output.state.vm[i] = solution.col_value[vm_offset + i] /
                column_scale[vm_offset + i];
            output.state.va[i] = solution.col_value[va_offset + i] /
                column_scale[va_offset + i];
            if (compact_economic_objective) {
                output.state.vm[i] += linearization_vm[i] -
                    solution.col_value[compact_vm_down_offset + i] /
                        column_scale[compact_vm_down_offset + i];
                output.state.va[i] += reference.va[i] -
                    solution.col_value[compact_va_down_offset + i] /
                        column_scale[compact_va_down_offset + i];
            }
        }
        for (int i = 0; i < ng; ++i) {
            output.state.pg[i] = solution.col_value[pg_offset + i] /
                column_scale[pg_offset + i];
            output.state.qg[i] = solution.col_value[qg_offset + i] /
                column_scale[qg_offset + i];
            if (compact_economic_objective) {
                output.state.pg[i] += reference.pg[i] -
                    solution.col_value[compact_pg_down_offset + i] /
                        column_scale[compact_pg_down_offset + i];
                output.state.qg[i] += reference.qg[i] -
                    solution.col_value[compact_qg_down_offset + i] /
                        column_scale[compact_qg_down_offset + i];
            }
        }
        if (outaged_generator >= 0) {
            output.state.pg[outaged_generator] = 0.0;
            output.state.qg[outaged_generator] = 0.0;
        }
        for (int i = 0; i < nd; ++i) {
            output.state.demand_factor[i] =
                solution.col_value[demand_offset + i] /
                column_scale[demand_offset + i];
            if (compact_economic_objective) {
                output.state.demand_factor[i] +=
                    reference.demand_factor[i] -
                    solution.col_value[compact_demand_down_offset + i] /
                        column_scale[compact_demand_down_offset + i];
            }
        }
        if (economic_objective && !compact_economic_objective) {
            output.state.sm_slack.resize(data.branches.size());
            for (int i = 0; i < nl; ++i) {
                output.state.sm_slack[i] =
                    solution.col_value[sm_offset + i] /
                    column_scale[sm_offset + i];
            }
        }
    }
    output.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_start).count();
    if (highs_log != nullptr && std::string(highs_log) != "0") {
        std::cerr << "GRAVITYX_LINEAR_SEED_RESULT "
                  << output.to_json(false).dump() << '\n';
    }
    return output;
}

}  // namespace gravityx
