#include "gravityx/state_io.hpp"

#include <stdexcept>
#include <string>
#include <vector>

namespace gravityx {
namespace {

std::vector<double> vector_or_empty(const nlohmann::json& value, const char* key) {
    const auto item = value.find(key);
    if (item == value.end() || item->is_null()) {
        return {};
    }
    if (!item->is_array()) {
        throw std::runtime_error(std::string("state field is not an array: ") + key);
    }
    return item->get<std::vector<double>>();
}

}  // namespace

nlohmann::json ac_state_to_json(const AcState& state) {
    return {
        {"vm", state.vm},
        {"va", state.va},
        {"pg", state.pg},
        {"qg", state.qg},
        {"demand_factor", state.demand_factor},
        {"pf", state.pf},
        {"qf", state.qf},
        {"pt", state.pt},
        {"qt", state.qt},
        {"sm_slack", state.sm_slack},
        {"p_delta", state.p_delta},
        {"q_delta", state.q_delta},
        {"commitment", state.commitment},
        {"startup", state.startup},
        {"shutdown", state.shutdown},
        {"gen_lambda", state.gen_lambda},
        {"load_lambda", state.load_lambda},
    };
}

nlohmann::json ac_submission_state_to_json(const AcState& state) {
    return {
        {"vm", state.vm},
        {"va", state.va},
        {"pg", state.pg},
        {"qg", state.qg},
        {"demand_factor", state.demand_factor},
    };
}

AcState ac_state_from_json(const nlohmann::json& value) {
    if (!value.is_object()) {
        throw std::runtime_error("AC state must be a JSON object");
    }
    AcState state;
    state.vm = vector_or_empty(value, "vm");
    state.va = vector_or_empty(value, "va");
    state.pg = vector_or_empty(value, "pg");
    state.qg = vector_or_empty(value, "qg");
    state.demand_factor = vector_or_empty(value, "demand_factor");
    state.pf = vector_or_empty(value, "pf");
    state.qf = vector_or_empty(value, "qf");
    state.pt = vector_or_empty(value, "pt");
    state.qt = vector_or_empty(value, "qt");
    state.sm_slack = vector_or_empty(value, "sm_slack");
    state.p_delta = vector_or_empty(value, "p_delta");
    state.q_delta = vector_or_empty(value, "q_delta");
    state.commitment = vector_or_empty(value, "commitment");
    state.startup = vector_or_empty(value, "startup");
    state.shutdown = vector_or_empty(value, "shutdown");
    state.gen_lambda = vector_or_empty(value, "gen_lambda");
    state.load_lambda = vector_or_empty(value, "load_lambda");
    return state;
}

nlohmann::json solve_result_to_json(const SolveResult& result, bool include_state) {
    nlohmann::json value = {
        {"status", result.status},
        {"objective", result.objective},
        {"wall_seconds", result.wall_seconds},
        {"iterations", result.iterations},
        {"resident_reoptimization", result.resident_reoptimization},
        {"acceptable_termination_enabled", result.acceptable_termination_enabled},
    };
    if (include_state) {
        value["state"] = ac_state_to_json(result.state);
    }
    return value;
}

nlohmann::json solve_result_to_submission_json(const SolveResult& result) {
    auto value = solve_result_to_json(result, false);
    value["state"] = ac_submission_state_to_json(result.state);
    return value;
}

}  // namespace gravityx
