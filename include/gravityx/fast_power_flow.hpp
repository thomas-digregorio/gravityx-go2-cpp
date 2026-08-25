#pragma once

#include "gravityx/ac_model.hpp"
#include "gravityx/validation.hpp"

#include <nlohmann/json.hpp>

#include <memory>
#include <string>
#include <vector>

namespace gravityx {

struct FastPowerFlowOptions {
    bool distributed_balance_polish{true};
    bool fixed_jacobian_screen_only{false};
    bool capture_diagnostics{false};
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
    bool fixed_jacobian_predictor_attempted{};
    bool fixed_jacobian_predictor_selected{};
    int fixed_jacobian_predictor_iterations{};
    double fixed_jacobian_predictor_preparation_seconds{};
    ValidationReport fixed_jacobian_predictor_validation;
    nlohmann::json fixed_jacobian_predictor_trace = nlohmann::json::array();
    bool newton_candidate_selected{};
    ValidationReport newton_candidate_validation;
    bool active_only_newton_attempted{};
    bool active_only_newton_selected{};
    bool active_only_newton_converged{};
    int active_only_newton_iterations{};
    ValidationReport active_only_newton_validation;
    bool reactive_only_newton_attempted{};
    bool reactive_only_newton_selected{};
    bool reactive_only_newton_converged{};
    int reactive_only_newton_iterations{};
    ValidationReport reactive_only_newton_validation;
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
    struct FixedJacobianPredictorCache;
    mutable std::unique_ptr<FixedJacobianPredictorCache> predictor_cache_;

    FastPowerFlowResult solve_impl(
        const Contingency* contingency,
        const AcState* initial_state = nullptr,
        bool supplied_candidate_direct_only = false) const;
};

}  // namespace gravityx
