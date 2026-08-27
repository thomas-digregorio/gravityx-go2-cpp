#pragma once

#include "gravityx/ac_model.hpp"
#include "gravityx/case_data.hpp"

#include <nlohmann/json.hpp>

#include <string>
#include <vector>

namespace gravityx {

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
    int row_count{};
    int column_count{};
    int nonzero_count{};
    int branch_security_row_count{};
    int simplex_strategy{-1};
    int run_status{};
    int model_status{};
    int primal_solution_status{};
    int num_primal_infeasibilities{};
    int iterations{};
    double max_primal_infeasibility{};
    double maximum_linearized_violation{};
    double maximum_column_violation{};
    bool finite_solution_values{};
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
    bool current_security_rows_only = false);

}  // namespace gravityx
