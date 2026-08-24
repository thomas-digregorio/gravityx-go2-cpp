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
    const std::optional<ContingencyContext>& contingency = std::nullopt);

}  // namespace gravityx
