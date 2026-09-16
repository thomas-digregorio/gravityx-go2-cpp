#pragma once

#include "gravityx/ac_model.hpp"
#include "gravityx/validation.hpp"

#include <nlohmann/json.hpp>

#include <cstddef>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace gravityx {

double default_fixed_jacobian_screen_seconds(
    std::size_t bus_count, ContingencyType contingency_type);

struct FastPowerFlowOptions {
    bool distributed_balance_polish{true};
    bool enable_fixed_jacobian_predictor{true};
    // Only factorization-derived unit-vector solves; never prior solutions.
    bool cache_outage_inverse_rows{true};
    bool adaptive_jacobian_refresh{false};
    int adaptive_jacobian_refresh_window{4};
    int max_adaptive_jacobian_refreshes{2};
    bool bounded_voltage_extrapolation{false};
    // Candidate controls only; flows, slacks and PWL fields are rebuilt before
    // their respective physical/final checks. False retains the copy oracle.
    bool controls_only_trial_copy{false};
    // Sum freshly validated terminal flows for polish injections. False keeps
    // the complex-admittance oracle; neither path is an acceptance check.
    bool reuse_polish_branch_flows{false};
    // Disabling this uses the existing decoupled correction directly, allowing
    // tiny tests to cover the production fallback with both injection paths.
    bool coupled_polish_correction{true};
    // Size routing policy. Tiny tests can exercise the identical predictor
    // and economic-incumbent path without manufacturing a large network.
    std::size_t fixed_jacobian_minimum_bus_count{16000};
    // Linearization only. Source corrective bounds and every validator keep
    // using the constructor's immutable base_state, never this candidate.
    const AcState* fixed_jacobian_linearization_state{};
    bool fixed_jacobian_screen_only{false};
    // Explicit cooperative override. Infinity selects the class policy for
    // first-stage screens; ordinary repair calls remain unbounded here.
    // The caller's process deadline remains the hard limit.
    double fixed_jacobian_time_limit_seconds{
        std::numeric_limits<double>::infinity()};
    bool economic_balance_polish{false};
    bool minimize_active_balance_slack{false};
    bool minimize_reactive_balance_slack{false};
    double balance_cleanup_fraction{1.0};
    bool capture_diagnostics{false};
    int max_economic_balance_polish_iterations{4};
    // Work-allocation target on summed P/Q slack, not a feasibility tolerance.
    double economic_balance_polish_stop_slack{1e-7};
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

// Opt-in cached economic cleanup plus bounded corrective-feasibility probes.
// This never invokes either optional per-contingency LP polish stage.
inline void enable_cached_economic_polish(FastPowerFlowOptions& options) {
    options.economic_balance_polish = true;
    options.adaptive_jacobian_refresh = true;
    options.bounded_voltage_extrapolation = true;
    options.controls_only_trial_copy = true;
    options.reuse_polish_branch_flows = true;
    options.max_economic_balance_polish_iterations = 3;
    options.economic_balance_polish_stop_slack = 0.025;
    options.economic_balance_polish_objective_threshold =
        std::numeric_limits<double>::infinity();
    options.max_economic_linearized_polish_rounds = 0;
    options.max_economic_linearized_phase_two_rounds = 0;
}

inline bool needs_contingency_economic_cleanup(
    const FastPowerFlowOptions& options, std::size_t bus_count,
    bool base_mode, bool direct_only) {
    // A seed-bank probe is a validator, not another optimization attempt.
    return !base_mode && !direct_only &&
        bus_count >= options.fixed_jacobian_minimum_bus_count &&
        options.economic_balance_polish;
}

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
    bool fixed_jacobian_budget_exhausted{};
    double effective_predictor_time_limit_seconds{
        std::numeric_limits<double>::infinity()};
    int fixed_jacobian_predictor_iterations{};
    double fixed_jacobian_predictor_preparation_seconds{};
    int outage_update_requests{};
    int outage_update_basis_cache_hits{};
    int outage_update_basis_cache_misses{};
    int outage_update_rhs_columns{};
    std::size_t outage_update_basis_cache_bytes{};
    double outage_update_seconds{};
    double outage_update_rhs_seconds{};
    double economic_balance_polish_seconds{};
    double economic_balance_polish_correction_seconds{};
    int economic_balance_polish_flow_reuses{};
    int economic_balance_polish_ybus_builds{};
    double economic_balance_polish_injection_seconds{};
    double economic_balance_polish_ybus_seconds{};
    std::size_t corrective_trial_copy_count{};
    std::size_t corrective_trial_control_copy_count{};
    std::size_t corrective_trial_copy_avoided_bytes{};
    double corrective_trial_copy_seconds{};
    int adaptive_jacobian_refresh_attempts{};
    int adaptive_jacobian_refresh_selected{};
    double adaptive_jacobian_refresh_seconds{};
    double adaptive_jacobian_refresh_best_before{};
    double adaptive_jacobian_refresh_best_after{};
    int voltage_extrapolation_searches{};
    int voltage_extrapolation_trials{};
    int voltage_extrapolation_selected{};
    double voltage_extrapolation_seconds{};
    double voltage_extrapolation_largest_scale{};
    double voltage_extrapolation_best_before{};
    double voltage_extrapolation_best_after{};
    ValidationReport fixed_jacobian_predictor_validation;
    nlohmann::json fixed_jacobian_predictor_trace = nlohmann::json::array();
    bool economic_balance_polish_attempted{};
    bool economic_direct_candidate_verified{};
    bool economic_direct_incumbent_selected{};
    double economic_direct_candidate_objective{};
    bool common_reference_base_candidate_verified{};
    bool common_reference_base_candidate_selected{};
    double common_reference_base_candidate_objective{};
    bool economic_balance_polish_threshold_passed{};
    double economic_balance_polish_objective_threshold{};
    bool economic_balance_polish_selected{};
    int economic_balance_polish_iterations{};
    int economic_balance_polish_backtracking_attempts{};
    int economic_balance_polish_trial_count{};
    int economic_balance_polish_physical_rejections{};
    int economic_balance_polish_economic_checks{};
    double economic_balance_polish_physical_check_seconds{};
    double economic_balance_polish_economic_check_seconds{};
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
    nlohmann::json economic_summary_json() const;
    nlohmann::json runtime_profile_json() const;
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
void run_outage_inverse_row_cache_regression();
void run_adaptive_jacobian_policy_regression();
void run_voltage_extrapolation_policy_regression();
void run_corrective_trial_copy_regression();
void run_polish_flow_injection_regression(
    const CaseData& data, const AcState& original_base);
void run_voltage_extrapolation_physics_regression(
    const CaseData& data, const std::vector<int>& commitment,
    const AcState& original_base, const Contingency& contingency,
    const AcState& known_feasible_state);

void run_economic_polish_trial_regression(
    const CaseData& data, const std::vector<int>& commitment,
    const AcState& original_base);

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

// Candidate generation only: original corrective bounds/cost interval and
// ratings, with no component removed. This is never a source contingency.
double rebuild_common_corrective_reference_state(
    const CaseData& data, const AcState& original_base,
    const std::vector<int>& commitment, AcState& state);

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
