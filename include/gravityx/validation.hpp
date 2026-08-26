#pragma once

#include "gravityx/ac_model.hpp"
#include "gravityx/case_data.hpp"

#include <nlohmann/json.hpp>

#include <optional>
#include <string>

namespace gravityx {

struct ValidationReport {
    double max_variable_bound_violation{};
    double max_pwl_sum_residual{};
    double max_pwl_power_residual{};
    double max_reference_angle_residual{};
    double max_generator_residual{};
    double max_load_ramp_violation{};
    double max_active_balance_residual{};
    double max_reactive_balance_residual{};
    double max_ohms_residual{};
    double max_angle_violation{};
    double max_flow_limit_violation{};
    double max_residual{};
    std::string worst_category;
    std::string worst_identity;

    nlohmann::json to_json() const;
};

ValidationReport validate_state(
    const CaseData& data,
    ModelMode mode,
    const AcState& state,
    const std::vector<int>& fixed_status = {},
    const std::optional<ContingencyContext>& contingency = std::nullopt);

// Trial corrections have just had their branch flows and balance slacks
// rebuilt.  Rank them without rechecking PWL interpolation or Ohm's law;
// every selected state is rebuilt economically and passed to validate_state
// before it can be accepted.
ValidationReport validate_rebuilt_contingency_trial(
    const CaseData& data,
    const AcState& state,
    const std::vector<int>& fixed_status,
    const ContingencyContext& contingency);

}  // namespace gravityx
