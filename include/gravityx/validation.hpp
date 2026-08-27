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
// rebuilt. Rank them without rechecking PWL interpolation or Ohm's law.
ValidationReport validate_rebuilt_contingency_trial(
    const CaseData& data,
    const AcState& state,
    const std::vector<int>& fixed_status,
    const ContingencyContext& contingency);

// Rank an internal rebuilt trial against a known incumbent residual. The
// returned report is complete whenever the trial can improve the incumbent.
// Once a partial maximum proves that the caller's strict 1e-10 improvement
// test cannot pass, the remaining invariant scans are skipped.
ValidationReport validate_rebuilt_contingency_trial_until_rejected(
    const CaseData& data,
    const AcState& state,
    const std::vector<int>& fixed_status,
    const ContingencyContext& contingency,
    double incumbent_max_residual,
    int preferred_balance_bus = -1,
    int preferred_branch = -1);

// Predictor iteration selection also needs the exact identity of the
// dominant physical residual (for example, whether a variable-bound maximum
// is a branch apparent-power slack).  Preserve that routing information while
// retaining the same rebuilt-field shortcuts as trial validation.
ValidationReport validate_rebuilt_contingency_predictor(
    const CaseData& data,
    const AcState& state,
    const std::vector<int>& fixed_status,
    const ContingencyContext& contingency);

// Complete a predictor's already-issued physical validation certificate after
// its deterministic economic fields have been rebuilt. This checks exactly
// the omitted PWL interpolation and branch Ohm-law constraints and merges them
// into the physical report without rescanning the unchanged physical state.
ValidationReport validate_rebuilt_contingency_economic_and_ohms(
    const CaseData& data,
    const AcState& state,
    const std::vector<int>& fixed_status,
    const ContingencyContext& contingency,
    const ValidationReport& physical_report);

}  // namespace gravityx
