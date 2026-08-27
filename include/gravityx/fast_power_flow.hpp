#pragma once

#include "gravityx/ac_model.hpp"
#include "gravityx/validation.hpp"

#include <nlohmann/json.hpp>

#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace gravityx {

struct FastPowerFlowOptions {
    bool distributed_balance_polish{true};
    bool fixed_jacobian_screen_only{false};
    bool economic_balance_polish{false};
    bool minimize_active_balance_slack{false};
    bool minimize_reactive_balance_slack{false};
    double balance_cleanup_fraction{1.0};
    bool capture_diagnostics{false};
    int max_economic_balance_polish_iterations{4};
    int max_economic_linearized_polish_rounds{1};
    double economic_linearized_polish_seconds{0.75};
    double economic_linearized_trigger_slack{0.05};
    double economic_balance_polish_objective_threshold{
        std::numeric_limits<double>::infinity()};
    int max_economic_linearized_phase_two_rounds{1};
    double economic_linearized_phase_two_seconds{0.75};
    int max_newton_iterations{50};
    int max_active_redispatch_passes{20};
    int max_reactive_limit_passes{8};
    double newton_tolerance{1e-10};
    double validation_tolerance{1e-5};
};

struct FastPowerFlowResult {
    bool converged{};
    bool feasible{};
    bool direct_candidate_attempted{};
    bool direct_candidate_selected{};
    ValidationReport direct_candidate_validation;
    bool local_balance_candidate_attempted{};
    bool local_balance_candidate_selected{};
    int local_balance_backtracking_attempts{};
    double local_balance_selected_step{};
    ValidationReport local_balance_candidate_validation;
    bool fixed_jacobian_predictor_attempted{};
    bool fixed_jacobian_predictor_selected{};
    int fixed_jacobian_predictor_iterations{};
    double fixed_jacobian_predictor_preparation_seconds{};
    ValidationReport fixed_jacobian_predictor_validation;
    nlohmann::json fixed_jacobian_predictor_trace = nlohmann::json::array();
    bool economic_balance_polish_attempted{};
    bool economic_balance_polish_threshold_passed{};
    double economic_balance_polish_objective_threshold{};
    bool economic_balance_polish_selected{};
    int economic_balance_polish_iterations{};
    int economic_balance_polish_backtracking_attempts{};
    double economic_balance_polish_objective_before{};
    double economic_balance_polish_objective_after{};
    double economic_balance_polish_active_slack_before{};
    double economic_balance_polish_active_slack_after{};
    double economic_balance_polish_reactive_slack_before{};
    double economic_balance_polish_reactive_slack_after{};
    ValidationReport economic_balance_polish_validation;
    nlohmann::json economic_balance_polish_trace = nlohmann::json::array();
    bool newton_candidate_selected{};
    ValidationReport newton_candidate_validation;
    bool active_only_newton_attempted{};
    bool active_only_newton_selected{};
    bool active_only_newton_converged{};
    int active_only_newton_iterations{};
    ValidationReport active_only_newton_validation;
    int active_only_backtracking_attempts{};
    double active_only_selected_step{};
    ValidationReport active_only_backtracking_validation;
    bool reactive_only_newton_attempted{};
    bool reactive_only_newton_selected{};
    bool reactive_only_newton_converged{};
    int reactive_only_newton_iterations{};
    ValidationReport reactive_only_newton_validation;
    int reactive_only_backtracking_attempts{};
    double reactive_only_selected_step{};
    nlohmann::json reactive_only_trace = nlohmann::json::array();
    bool distributed_balance_polish_attempted{};
    bool distributed_balance_polish_selected{};
    int distributed_balance_polish_iterations{};
    int distributed_balance_voltage_projections{};
    std::string distributed_balance_polish_failure_reason;
    ValidationReport distributed_balance_polish_validation;
    bool best_intermediate_candidate_selected{};
    std::string best_intermediate_candidate_source;
    ValidationReport best_intermediate_candidate_validation;
    int newton_iterations{};
    double initial_newton_residual{};
    int active_redispatch_passes{};
    int reactive_limit_passes{};
    double wall_seconds{};
    std::string failure_reason;
    SolveResult solve;
    ValidationReport validation;

    nlohmann::json to_json() const;
};

struct ValidatedSourceBaseResult {
    bool feasible{};
    double wall_seconds{};
    SolveResult solve;
    ValidationReport validation;

    nlohmann::json to_json() const;
};

ValidatedSourceBaseResult build_validated_source_base(
    const CaseData& data,
    std::vector<int> commitment,
    double validation_tolerance = 1e-5);

void run_fast_power_flow_topology_cache_regression();

double rebuild_base_state_derived_fields(
    const CaseData& data,
    const std::vector<int>& commitment,
    AcState& state,
    double balance_slack_upper = 0.5);

double rebuild_contingency_state_derived_fields(
    const CaseData& data,
    const AcState& base_state,
    const std::vector<int>& commitment,
    const Contingency& contingency,
    AcState& state,
    double balance_slack_upper = 0.5);

class FastContingencyPowerFlow {
public:
    FastContingencyPowerFlow(
        const CaseData& data,
        const AcState& base_state,
        std::vector<int> commitment,
        FastPowerFlowOptions options = {});
    ~FastContingencyPowerFlow();

    FastContingencyPowerFlow(const FastContingencyPowerFlow&) = delete;
    FastContingencyPowerFlow& operator=(const FastContingencyPowerFlow&) = delete;

    FastPowerFlowResult solve_base() const;
    FastPowerFlowResult solve(const Contingency& contingency) const;
    FastPowerFlowResult solve(
        const Contingency& contingency,
        const AcState& initial_state) const;
    FastPowerFlowResult screen_candidate(
        const Contingency& contingency,
        const AcState& candidate_state) const;

private:
    const CaseData& data_;
    const AcState& base_state_;
    std::vector<int> commitment_;
    FastPowerFlowOptions options_;
    std::vector<std::vector<int>> base_components_;
    std::vector<unsigned char> bridge_branch_;
    struct FixedJacobianPredictorCache;
    mutable std::unique_ptr<FixedJacobianPredictorCache> predictor_cache_;

    FastPowerFlowResult solve_impl(
        const Contingency* contingency,
        const AcState* initial_state = nullptr,
        bool supplied_candidate_direct_only = false) const;
};

}  // namespace gravityx
