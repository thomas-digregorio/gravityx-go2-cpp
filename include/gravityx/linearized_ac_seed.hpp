#pragma once

#include "gravityx/ac_model.hpp"
#include "gravityx/case_data.hpp"

#include <nlohmann/json.hpp>

#include <optional>
#include <string>
#include <vector>

namespace gravityx {

struct LinearizedAcSeedResult {
    bool success{};
    double wall_seconds{};
    bool economic_objective{};
    bool compact_economic_objective{};
    bool projected_balance_slack{};
    bool branch_security_rows_omitted{};
    int branch_security_subset_count{};
    bool feasibility_only{};
    bool elastic_balance_phase_one{};
    bool primal_start_attempted{};
    bool primal_start_accepted{};
    int primal_start_status{};
    bool primal_start_finite{true};
    int primal_start_nonfinite_count{};
    double primal_start_maximum_column_violation{};
    int primal_start_worst_column{-1};
    double primal_start_maximum_row_violation{};
    int primal_start_worst_row{-1};
    bool primal_basis_attempted{};
    bool primal_basis_accepted{};
    int primal_basis_status{};
    bool presolve_enabled{true};
    double primal_simplex_bound_perturbation_multiplier{-1.0};
    int simplex_strategy{};
    double maximum_column_scale{1.0};
    double maximum_row_scale{1.0};
    double objective_scale{1.0};
    double voltage_trust_radius{};
    double angle_trust_radius{};
    int projected_reference_voltage_count{};
    double maximum_reference_voltage_projection{};
    double maximum_balance_slack{};
    double total_balance_slack{};
    bool solution_value_valid{};
    bool info_valid{};
    bool accepted_feasible_nonoptimal_phase_one{};
    bool accepted_feasible_nonoptimal_economic{};
    bool accepted_approximate_economic_direction{};
    int resident_segment_count{};
    int feasible_segment_snapshot_count{};
    bool recovered_feasible_segment_snapshot{};
    bool canonicalized_segment_snapshot{};
    double canonicalized_snapshot_max_primal_infeasibility{};
    double canonicalized_snapshot_objective{};
    int terminal_run_status{-99};
    int terminal_model_status{-99};
    nlohmann::json resident_segments = nlohmann::json::array();
    double time_limit_seconds{};
    double ipm_optimality_tolerance{};
    int row_count{};
    int column_count{};
    int nonzero_count{};
    bool model_preflight_passed{};
    std::string model_preflight_failure;
    int tiny_matrix_entries_removed{};
    double maximum_tiny_matrix_entry_removed{};
    double small_matrix_value{};
    int add_vars_status{-99};
    int change_cols_cost_status{-99};
    int add_rows_status{-99};
    bool model_load_warning{};
    bool model_construction_success{};
    std::string model_load_failure_call;
    int run_status{};
    int model_status{};
    int primal_solution_status{};
    int num_primal_infeasibilities{};
    double max_primal_infeasibility{};
    int iterations{};
    double objective{};
    std::string status;
    AcState state;

    nlohmann::json to_json(bool include_state = false) const;
};

LinearizedAcSeedResult solve_linearized_ac_seed(
    const CaseData& data,
    const AcState& reference,
    const std::vector<int>& commitment,
    double balance_slack_limit = 0.49,
    const std::optional<ContingencyContext>& contingency = std::nullopt,
    bool project_balance_slack = false,
    bool request_lightweight_large_seed = false,
    double time_limit_seconds = 60.0,
    bool omit_branch_security_rows = false,
    bool feasibility_only = false,
    const std::vector<int>& branch_security_subset = {},
    bool force_projected_balance_phase_one = false,
    double voltage_trust_radius_override = -1.0,
    double angle_trust_radius_override = -1.0,
    bool economic_objective = false);

}  // namespace gravityx
