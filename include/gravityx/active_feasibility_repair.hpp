#pragma once

#include "gravityx/ac_model.hpp"
#include "gravityx/case_data.hpp"

#include <nlohmann/json.hpp>

#include <string>
#include <vector>

namespace gravityx {

// Tiny deterministic regression of numerical reconstruction, without a case solve.
void run_active_repair_canonicalization_regression();

struct ActiveFeasibilityRepairResult {
    bool success{};
    bool accepted_feasible_nonoptimal{};
    double wall_seconds{};
    double time_limit_seconds{};
    double balance_slack_limit{};
    double angle_trust_radius{};
    double voltage_trust_radius{};
    bool include_reactive{};
    bool current_security_rows_only{};
    bool include_component_box_rows{true};
    bool minimize_balance_slack{};
    bool local_search{};
    int local_control_bus_count{};
    int solver_row_count{};
    int solver_column_count{};
    int evaluated_branch_derivatives{};
    int fixed_branch_derivatives_skipped{};
    double branch_derivative_seconds{};
    int row_count{};
    int column_count{};
    int nonzero_count{};
    int branch_security_row_count{};
    int simplex_strategy{-1};
    int simplex_iteration_limit{-1};
    bool primal_start_attempted{};
    int primal_start_status{};
    double primal_start_maximum_row_violation{};
    double primal_start_maximum_column_violation{};
    bool primal_basis_attempted{};
    int primal_basis_status{};
    bool presolve_enabled{true};
    int run_status{};
    int model_status{};
    int primal_solution_status{};
    int num_primal_infeasibilities{};
    int iterations{};
    double max_primal_infeasibility{};
    double maximum_linearized_violation{};
    double maximum_column_violation{};
    bool finite_solution_values{};
    bool canonicalized_terminal_solution{};
    double maximum_angle_change{};
    double maximum_voltage_change{};
    double maximum_generation_change{};
    double maximum_reactive_generation_change{};
    double maximum_load_change{};
    double objective{};
    std::string solver;
    std::string status;
    AcState state;

    nlohmann::json to_json(bool include_state = false) const;
};

// Construct a compact source-bounded P/Q feasibility LP around an AC
// contingency state.  It changes voltage angles and magnitudes, active and
// reactive generator dispatch, and source-authorized corrective load.  The
// returned point is a candidate: callers must rebuild the nonlinear AC fields
// and run the complete independent validator before accepting it.
ActiveFeasibilityRepairResult solve_linearized_active_feasibility_repair(
    const CaseData& data,
    const AcState& reference,
    const std::vector<int>& commitment,
    const ContingencyContext& contingency,
    double balance_slack_limit = 0.25,
    double angle_trust_radius = 0.15,
    double time_limit_seconds = 5.0,
    double voltage_trust_radius = 0.02,
    bool include_reactive = true,
    bool current_security_rows_only = false,
    bool include_component_box_rows = true,
    bool minimize_balance_slack = false,
    // Optional search restriction, NOT a reduced acceptance model. 1 permits
    // local controls and imbalance, 2 permits only boundary imbalance, and
    // 0 fixes all controls and signed imbalance to the reference. Returned
    // candidates still require the complete original nonlinear validator.
    const std::vector<unsigned char>* local_bus_mask = nullptr,
    // False is the full-derivative oracle for tiny component comparisons.
    bool omit_fixed_derivatives = true);

// Source topology only, with a one-hop fixed-control boundary. Empty means
// the requested neighborhood exceeds the bounded work limit.
std::vector<unsigned char> contingency_repair_neighborhood(
    const CaseData& data, const Contingency& contingency,
    int depth = 3, std::size_t maximum_control_buses = 512);

}  // namespace gravityx
