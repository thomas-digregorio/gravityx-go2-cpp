#pragma once

#include "gravityx/ac_model.hpp"

#include <nlohmann/json.hpp>

namespace gravityx {

nlohmann::json ac_state_to_json(const AcState& state);
AcState ac_state_from_json(const nlohmann::json& value);
nlohmann::json solve_result_to_json(const SolveResult& result, bool include_state = false);

}  // namespace gravityx
