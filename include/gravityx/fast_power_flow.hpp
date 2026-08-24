#pragma once

#include "gravityx/ac_model.hpp"
#include "gravityx/validation.hpp"

#include <nlohmann/json.hpp>

#include <string>
#include <vector>

namespace gravityx {

struct FastPowerFlowOptions {
    bool distributed_balance_polish{true};
    int max_newton_iterations{50};
    int max_active_redispatch_passes{20};
    int max_reactive_limit_passes{8};
    double newton_tolerance{1e-10};
    double validation_tolerance{1e-5};
};

struct FastPowerFlowResult {
    bool converged{};
    bool feasible{};
    bool distributed_balance_polish_attempted{};
    bool distributed_balance_polish_selected{};
    int distributed_balance_polish_iterations{};
    int distributed_balance_voltage_projections{};
    std::string distributed_balance_polish_failure_reason;
    ValidationReport distributed_balance_polish_validation;
    int newton_iterations{};
    int active_redispatch_passes{};
    int reactive_limit_passes{};
    double wall_seconds{};
    std::string failure_reason;
    SolveResult solve;
    ValidationReport validation;

    nlohmann::json to_json() const;
};

class FastContingencyPowerFlow {
public:
    FastContingencyPowerFlow(
        const CaseData& data,
        const AcState& base_state,
        std::vector<int> commitment,
        FastPowerFlowOptions options = {});

    FastPowerFlowResult solve(const Contingency& contingency) const;

private:
    const CaseData& data_;
    const AcState& base_state_;
    std::vector<int> commitment_;
    FastPowerFlowOptions options_;
};

}  // namespace gravityx
