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
    bool projected_balance_slack{};
    double time_limit_seconds{};
    int model_status{};
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
    double time_limit_seconds = 60.0);

}  // namespace gravityx
