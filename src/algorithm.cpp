#include "gravityx/algorithm.hpp"
#include "gravityx/state_io.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <tuple>

namespace gravityx {
namespace {

bool solver_succeeded(int status) {
    return status == 0 || status == 1;
}

bool state_is_acceptable(const SolveResult& result, const ValidationReport& validation, double tolerance) {
    return solver_succeeded(result.status) && std::isfinite(result.objective) &&
        validation.max_residual <= tolerance;
}

int rounded_status(double value, int prior, double threshold) {
    if (std::abs(value - threshold) <= 1e-10) {
        return prior;
    }
    return value > threshold ? 1 : 0;
}

}  // namespace

nlohmann::json IbrRound::to_json() const {
    return {
        {"round", round},
        {"batch", batch},
        {"proposed_status", proposed_status},
        {"fallback_to_prior", fallback_to_prior},
        {"solve", solve_result_to_json(solve)},
        {"validation", validation.to_json()},
    };
}

nlohmann::json IbrResult::to_json(bool include_state) const {
    nlohmann::json rounds_json = nlohmann::json::array();
    for (const auto& round : rounds) {
        rounds_json.push_back(round.to_json());
    }
    nlohmann::json result = {
        {"success", success},
        {"candidate_accepted", candidate_accepted},
        {"wall_seconds", wall_seconds},
        {"switching_cost", switching_cost},
        {"candidate_proxy", candidate_proxy},
        {"base", solve_result_to_json(base)},
        {"base_validation", base_validation.to_json()},
        {"rounds", rounds_json},
        {"fixed_repair", solve_result_to_json(fixed_repair)},
        {"fixed_validation", fixed_validation.to_json()},
        {"commitment", commitment},
    };
    if (include_state) {
        result["selected_state"] = ac_state_to_json(selected_state);
    }
    return result;
}

IbrResult run_iterative_batch_rounding(const CaseData& data, const IbrOptions& options) {
    if (options.batch_count <= 0 || !(options.threshold > 0.0 && options.threshold < 1.0)) {
        throw std::runtime_error("invalid iterative batch-rounding options");
    }
    const auto wall_start = std::chrono::steady_clock::now();
    IbrResult output;

    std::vector<int> prior;
    std::vector<int> eligible;
    prior.reserve(data.generators.size());
    for (int i = 0; i < static_cast<int>(data.generators.size()); ++i) {
        const auto& gen = data.generators[i];
        prior.push_back(gen.status_prev);
        if ((gen.status_prev == 0 && gen.suqual == 1) ||
            (gen.status_prev == 1 && gen.sdqual == 1)) {
            eligible.push_back(i);
        }
    }

    AcModel base_model(data, ModelMode::BaseSoft, prior);
    output.base = base_model.solve(options.print_level, options.tolerance);
    output.base_validation = validate_state(data, ModelMode::BaseSoft, output.base.state, prior);
    if (!state_is_acceptable(output.base, output.base_validation, options.validation_tolerance)) {
        output.commitment = prior;
        output.selected_state = output.base.state;
        output.wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        return output;
    }

    AcModel relaxation(data, ModelMode::UnitCommitmentRelaxation);
    relaxation.initialize_from(output.base.state);
    IbrRound initial;
    initial.round = 0;
    initial.solve = relaxation.solve(options.print_level, options.tolerance);
    initial.validation = validate_state(data, ModelMode::UnitCommitmentRelaxation, initial.solve.state);
    output.rounds.push_back(initial);
    if (!state_is_acceptable(initial.solve, initial.validation, options.validation_tolerance)) {
        output.commitment = prior;
        output.selected_state = output.base.state;
        output.wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        return output;
    }

    std::sort(eligible.begin(), eligible.end(), [&](int left, int right) {
        const double left_distance = std::abs(initial.solve.state.commitment[left] - options.threshold);
        const double right_distance = std::abs(initial.solve.state.commitment[right] - options.threshold);
        if (std::abs(left_distance - right_distance) > 1e-15) {
            return left_distance > right_distance;
        }
        return data.generators[left].index < data.generators[right].index;
    });

    const int batch_count = std::min(options.batch_count, std::max(1, static_cast<int>(eligible.size())));
    std::vector<std::vector<int>> batches(batch_count);
    for (int position = 0; position < static_cast<int>(eligible.size()); ++position) {
        const int batch = std::min(batch_count - 1,
            position * batch_count / static_cast<int>(eligible.size()));
        batches[batch].push_back(eligible[position]);
    }

    std::vector<int> fixed(data.generators.size(), -1);
    AcState last_state = initial.solve.state;
    for (int batch_index = 0; batch_index < batch_count; ++batch_index) {
        if (batches[batch_index].empty()) {
            continue;
        }
        IbrRound round;
        round.round = batch_index + 1;
        round.batch = batches[batch_index];
        relaxation.initialize_from(last_state);
        for (int generator : round.batch) {
            const int proposal = rounded_status(
                last_state.commitment[generator], prior[generator], options.threshold);
            round.proposed_status.push_back(proposal);
            fixed[generator] = proposal;
            relaxation.set_commitment_bound(generator, proposal);
        }
        round.solve = relaxation.solve(options.print_level, options.tolerance);
        round.validation = validate_state(data, ModelMode::UnitCommitmentRelaxation, round.solve.state);
        bool acceptable = state_is_acceptable(round.solve, round.validation, options.validation_tolerance);
        if (acceptable) {
            for (int generator : round.batch) {
                acceptable = acceptable &&
                    std::abs(round.solve.state.commitment[generator] - fixed[generator]) <= options.validation_tolerance;
            }
        }
        if (!acceptable) {
            round.fallback_to_prior = true;
            relaxation.initialize_from(last_state);
            for (int generator : round.batch) {
                fixed[generator] = prior[generator];
                relaxation.set_commitment_bound(generator, prior[generator]);
            }
            round.solve = relaxation.solve(options.print_level, options.tolerance);
            round.validation = validate_state(data, ModelMode::UnitCommitmentRelaxation, round.solve.state);
            acceptable = state_is_acceptable(round.solve, round.validation, options.validation_tolerance);
        }
        output.rounds.push_back(round);
        if (!acceptable) {
            output.commitment = prior;
            output.selected_state = output.base.state;
            output.wall_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - wall_start).count();
            return output;
        }
        last_state = round.solve.state;
    }

    output.commitment = prior;
    for (int generator : eligible) {
        if (fixed[generator] >= 0) {
            output.commitment[generator] = fixed[generator];
        }
    }
    for (std::size_t i = 0; i < data.generators.size(); ++i) {
        if (output.commitment[i] != prior[i]) {
            output.switching_cost += output.commitment[i] == 1
                ? data.generators[i].sucost
                : data.generators[i].sdcost;
        }
    }

    AcModel fixed_model(data, ModelMode::BaseSoft, output.commitment);
    fixed_model.initialize_from(last_state);
    output.fixed_repair = fixed_model.solve(options.print_level, options.tolerance);
    output.fixed_validation = validate_state(
        data, ModelMode::BaseSoft, output.fixed_repair.state, output.commitment);
    output.candidate_proxy = output.fixed_repair.objective - output.switching_cost;
    output.candidate_accepted = state_is_acceptable(
        output.fixed_repair, output.fixed_validation, options.validation_tolerance) &&
        output.candidate_proxy > output.base.objective + options.validation_tolerance;

    if (output.candidate_accepted) {
        output.selected_state = output.fixed_repair.state;
    } else {
        output.commitment = prior;
        output.selected_state = output.base.state;
    }
    output.success = true;
    output.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_start).count();
    return output;
}

}  // namespace gravityx
