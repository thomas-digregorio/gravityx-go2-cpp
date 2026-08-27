#include "gravityx/algorithm.hpp"
#include "gravityx/fast_power_flow.hpp"
#include "gravityx/linearized_ac_seed.hpp"
#include "gravityx/state_io.hpp"

#include <algorithm>
#include <array>
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
    return solver_succeeded(result.status) &&
        validated_candidate_is_feasible(result, validation, tolerance);
}

int rounded_status(double value, int prior, double threshold) {
    if (std::abs(value - threshold) <= 1e-10) {
        return prior;
    }
    return value > threshold ? 1 : 0;
}

}  // namespace

bool validated_candidate_is_feasible(
    const SolveResult& result,
    const ValidationReport& validation,
    double tolerance) {
    return tolerance >= 0.0 && std::isfinite(result.objective) &&
        std::isfinite(validation.max_residual) &&
        validation.max_residual <= tolerance;
}

bool verified_economic_candidate_improves_incumbent(
    const SolveResult& incumbent,
    const ValidationReport& incumbent_validation,
    const SolveResult& candidate,
    const ValidationReport& candidate_validation,
    double validation_tolerance,
    double objective_tolerance) {
    return objective_tolerance >= 0.0 &&
        validated_candidate_is_feasible(
            incumbent, incumbent_validation, validation_tolerance) &&
        validated_candidate_is_feasible(
            candidate, candidate_validation, validation_tolerance) &&
        candidate.objective > incumbent.objective + objective_tolerance;
}

nlohmann::json SparseEconomicRefinementResult::to_json(
    bool include_state) const {
    nlohmann::json result = {
        {"incumbent_verified", incumbent_verified},
        {"attempted", attempted},
        {"improved", improved},
        {"time_limit_reached", time_limit_reached},
        {"wall_seconds", wall_seconds},
        {"incumbent_objective", incumbent_objective},
        {"selected_objective", selected_objective},
        {"selected", solve_result_to_json(selected, include_state)},
        {"selected_validation", selected_validation.to_json()},
        {"rounds", rounds},
    };
    return result;
}

SparseEconomicRefinementResult refine_fixed_commitment_sparse(
    const CaseData& data,
    const std::vector<int>& commitment,
    const SolveResult& incumbent,
    const SparseEconomicRefinementOptions& options) {
    if (commitment.size() != data.generators.size()) {
        throw std::runtime_error(
            "sparse economic refinement commitment size mismatch");
    }
    if (!std::isfinite(options.time_limit_seconds) ||
        options.time_limit_seconds <= 0.0 ||
        !std::isfinite(options.validation_tolerance) ||
        options.validation_tolerance < 0.0 ||
        !std::isfinite(options.objective_tolerance) ||
        options.objective_tolerance < 0.0 ||
        options.maximum_rounds <= 0 ||
        options.maximum_linear_economic_rounds < 0 ||
        !std::isfinite(options.linear_economic_time_limit_seconds) ||
        options.linear_economic_time_limit_seconds <= 0.0 ||
        !std::isfinite(options.linear_economic_voltage_trust_radius) ||
        options.linear_economic_voltage_trust_radius <= 0.0 ||
        !std::isfinite(options.linear_economic_angle_trust_radius) ||
        options.linear_economic_angle_trust_radius <= 0.0 ||
        options.maximum_voltage_coordinate_passes < 0 ||
        options.voltage_coordinate_bus_count <= 0) {
        throw std::runtime_error(
            "invalid sparse economic refinement options");
    }

    const auto wall_start = std::chrono::steady_clock::now();
    SparseEconomicRefinementResult output;
    output.selected = incumbent;
    output.selected.status = 0;
    output.selected.objective = rebuild_base_state_derived_fields(
        data, commitment, output.selected.state);
    output.selected_validation = validate_state(
        data, ModelMode::BaseSoft, output.selected.state, commitment);
    output.incumbent_objective = output.selected.objective;
    output.selected_objective = output.selected.objective;
    output.incumbent_verified = validated_candidate_is_feasible(
        output.selected, output.selected_validation,
        options.validation_tolerance);
    if (!output.incumbent_verified) {
        output.wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        return output;
    }

    auto reference = output.selected.state;
    for (int round = 1; round <= options.maximum_rounds; ++round) {
        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        const double remaining = options.time_limit_seconds - elapsed;
        if (remaining <= 0.05) {
            output.time_limit_reached = true;
            break;
        }
        output.attempted = true;
        FastPowerFlowOptions fast_options;
        fast_options.minimize_active_balance_slack = true;
        fast_options.minimize_reactive_balance_slack = true;
        const double incumbent_active_balance_slack = std::accumulate(
            output.selected.state.p_delta.begin(),
            output.selected.state.p_delta.end(), 0.0);
        const double incumbent_reactive_balance_slack = std::accumulate(
            output.selected.state.q_delta.begin(),
            output.selected.state.q_delta.end(), 0.0);
        fast_options.max_newton_iterations = 40;
        fast_options.max_active_redispatch_passes = 12;
        std::vector<double> cleanup_fractions;
        if (incumbent_reactive_balance_slack > 30.0) {
            cleanup_fractions = {1.0};
        } else if (incumbent_active_balance_slack <= 0.05) {
            cleanup_fractions = {0.2, 0.1, 0.05, 0.02};
        } else {
            cleanup_fractions = {0.2};
        }
        FastPowerFlowResult fast_candidate;
        bool have_fast_candidate = false;
        bool selected_fast_candidate_feasible = false;
        double selected_cleanup_fraction = cleanup_fractions.front();
        nlohmann::json fraction_trials = nlohmann::json::array();
        for (const double cleanup_fraction : cleanup_fractions) {
            fast_options.balance_cleanup_fraction = cleanup_fraction;
            FastContingencyPowerFlow fast_cleanup(
                data, reference, commitment, fast_options);
            auto trial = fast_cleanup.solve_base();
            trial.solve.objective = rebuild_base_state_derived_fields(
                data, commitment, trial.solve.state);
            trial.validation = validate_state(
                data, ModelMode::BaseSoft,
                trial.solve.state, commitment);
            const bool trial_feasible = validated_candidate_is_feasible(
                trial.solve, trial.validation,
                options.validation_tolerance);
            fraction_trials.push_back({
                {"fraction", cleanup_fraction},
                {"objective", trial.solve.objective},
                {"feasible", trial_feasible},
                {"active_balance_slack", std::accumulate(
                    trial.solve.state.p_delta.begin(),
                    trial.solve.state.p_delta.end(), 0.0)},
                {"reactive_balance_slack", std::accumulate(
                    trial.solve.state.q_delta.begin(),
                    trial.solve.state.q_delta.end(), 0.0)},
                {"validation", trial.validation.to_json()},
                {"wall_seconds", trial.wall_seconds},
            });
            const bool select_trial = !have_fast_candidate ||
                (trial_feasible && !selected_fast_candidate_feasible) ||
                (trial_feasible == selected_fast_candidate_feasible &&
                 trial.solve.objective >
                    fast_candidate.solve.objective +
                        options.objective_tolerance);
            if (select_trial) {
                fast_candidate = std::move(trial);
                selected_cleanup_fraction = cleanup_fraction;
                selected_fast_candidate_feasible = trial_feasible;
                have_fast_candidate = true;
            }
        }
        nlohmann::json fast_round = {
            {"round", round},
            {"remaining_seconds_at_start", remaining},
            {"incumbent_reactive_balance_slack",
             incumbent_reactive_balance_slack},
            {"balance_cleanup_fraction", selected_cleanup_fraction},
            {"fraction_trials", std::move(fraction_trials)},
            {"fast_balance_cleanup", fast_candidate.to_json()},
        };
        const bool fast_improved =
            verified_economic_candidate_improves_incumbent(
                output.selected, output.selected_validation,
                fast_candidate.solve, fast_candidate.validation,
                options.validation_tolerance,
                options.objective_tolerance);
        fast_round["accepted"] = fast_improved;
        if (fast_improved) {
            output.selected = std::move(fast_candidate.solve);
            output.selected_validation = fast_candidate.validation;
            output.selected_objective = output.selected.objective;
            output.improved = true;
            reference = output.selected.state;
            const double active_balance_slack = std::accumulate(
                output.selected.state.p_delta.begin(),
                output.selected.state.p_delta.end(), 0.0);
            fast_round["selected_objective"] = output.selected.objective;
            fast_round["selected_active_balance_slack"] =
                active_balance_slack;
            output.rounds.push_back(std::move(fast_round));
            if (active_balance_slack <= 1e-9) {
                break;
            }
            continue;
        }
        // The function admitted the incumbent only after a complete nonlinear
        // validation.  Once the fast cleanup no longer improves that verified
        // state, another feasibility Phase I is redundant and can consume the
        // entire economic budget.  Continue directly to the actual market-
        // surplus LP below; it already carries bounded balance-slack columns
        // and every proposed state must pass the same independent validator.
        if (output.incumbent_verified) {
            fast_round["status"] =
                "verified_balance_cleanup_stagnated_phase_one_skipped";
            fast_round["phase_one_skipped"] = true;
            output.rounds.push_back(std::move(fast_round));
            break;
        }
        output.rounds.push_back(std::move(fast_round));

        const double elapsed_after_fast = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        const double phase_one_remaining =
            options.time_limit_seconds - elapsed_after_fast;
        if (phase_one_remaining <= 0.05) {
            output.time_limit_reached = true;
            break;
        }
        const auto linear = solve_linearized_ac_seed(
            data, reference, commitment, 0.49, std::nullopt,
            true, true, phase_one_remaining, true, true);
        nlohmann::json round_json = {
            {"round", round},
            {"remaining_seconds_at_start", phase_one_remaining},
            {"phase_one", linear.to_json(false)},
            {"candidates", nlohmann::json::array()},
        };
        if (!linear.success) {
            round_json["status"] = "phase_one_failed";
            output.rounds.push_back(std::move(round_json));
            break;
        }

        SolveResult best_round = output.selected;
        ValidationReport best_round_validation = output.selected_validation;
        double best_round_fraction = 0.0;
        const std::array<double, 5> fractions{1.0, 0.75, 0.5, 0.25, 0.1};
        for (const double fraction : fractions) {
            AcState candidate_state = reference;
            const auto interpolate = [fraction](
                const std::vector<double>& from,
                const std::vector<double>& to,
                std::vector<double>& target) {
                target.resize(from.size());
                for (std::size_t i = 0; i < from.size(); ++i) {
                    target[i] = from[i] + fraction * (to[i] - from[i]);
                }
            };
            interpolate(reference.vm, linear.state.vm, candidate_state.vm);
            interpolate(reference.va, linear.state.va, candidate_state.va);
            interpolate(reference.pg, linear.state.pg, candidate_state.pg);
            interpolate(reference.qg, linear.state.qg, candidate_state.qg);
            interpolate(
                reference.demand_factor, linear.state.demand_factor,
                candidate_state.demand_factor);

            SolveResult candidate;
            candidate.status = 0;
            candidate.iterations = linear.iterations;
            candidate.wall_seconds = linear.wall_seconds;
            candidate.state = std::move(candidate_state);
            candidate.objective = rebuild_base_state_derived_fields(
                data, commitment, candidate.state);
            const auto candidate_validation = validate_state(
                data, ModelMode::BaseSoft, candidate.state, commitment);
            const bool accepted =
                verified_economic_candidate_improves_incumbent(
                    best_round, best_round_validation,
                    candidate, candidate_validation,
                    options.validation_tolerance,
                    options.objective_tolerance);
            round_json["candidates"].push_back({
                {"source", "phase_one_blend"},
                {"fraction", fraction},
                {"objective", candidate.objective},
                {"accepted", accepted},
                {"validation", candidate_validation.to_json()},
            });
            if (accepted) {
                best_round = std::move(candidate);
                best_round_validation = candidate_validation;
                best_round_fraction = fraction;
            }
        }

        const double elapsed_before_repair = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        if (elapsed_before_repair + 0.05 < options.time_limit_seconds) {
            FastContingencyPowerFlow nonlinear_repair(
                data, linear.state, commitment);
            auto repaired = nonlinear_repair.solve_base();
            repaired.solve.objective = rebuild_base_state_derived_fields(
                data, commitment, repaired.solve.state);
            repaired.validation = validate_state(
                data, ModelMode::BaseSoft,
                repaired.solve.state, commitment);
            const bool accepted =
                verified_economic_candidate_improves_incumbent(
                    best_round, best_round_validation,
                    repaired.solve, repaired.validation,
                    options.validation_tolerance,
                    options.objective_tolerance);
            round_json["candidates"].push_back({
                {"source", "phase_one_sparse_newton"},
                {"objective", repaired.solve.objective},
                {"accepted", accepted},
                {"validation", repaired.validation.to_json()},
                {"repair", repaired.to_json()},
            });
            if (accepted) {
                best_round = std::move(repaired.solve);
                best_round_validation = repaired.validation;
                best_round_fraction = -1.0;
            }
        }

        const bool round_improved =
            verified_economic_candidate_improves_incumbent(
                output.selected, output.selected_validation,
                best_round, best_round_validation,
                options.validation_tolerance,
                options.objective_tolerance);
        round_json["accepted"] = round_improved;
        round_json["accepted_fraction"] = best_round_fraction;
        round_json["selected_objective"] = best_round.objective;
        round_json["selected_validation"] =
            best_round_validation.to_json();
        output.rounds.push_back(std::move(round_json));
        if (!round_improved) {
            break;
        }
        output.selected = std::move(best_round);
        output.selected_validation = best_round_validation;
        output.selected_objective = output.selected.objective;
        output.improved = true;
        reference = output.selected.state;

        const double active_balance_slack = std::accumulate(
            output.selected.state.p_delta.begin(),
            output.selected.state.p_delta.end(), 0.0);
        if (active_balance_slack <= 1e-9) {
            break;
        }
    }

    // Once the sparse feasibility cleanup has produced a verified incumbent,
    // use a bounded sequential LP to optimize the actual fixed-commitment GO2
    // market-surplus terms.  The LP contains exact PWL generator/load
    // economics and first-order AC balance/thermal rows.  It is only a
    // direction generator: every blend and nonlinear repair is rebuilt and
    // independently validated before it may replace the incumbent.
    for (int economic_round = 1;
         economic_round <= options.maximum_linear_economic_rounds;
         ++economic_round) {
        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        const double remaining = options.time_limit_seconds - elapsed;
        // Reserve a small tail for the exact voltage-coordinate search.  On
        // the 19k pilot those independently verified nonlinear moves produced
        // most of the base-score gain in well under two seconds; letting an
        // interrupted LP consume the final fraction of the budget discarded
        // that high-value deterministic cleanup.
        constexpr double kVerifiedCoordinateReserveSeconds = 2.0;
        if (remaining <= kVerifiedCoordinateReserveSeconds + 0.25) {
            break;
        }
        const double lp_time_limit = std::min(
            options.linear_economic_time_limit_seconds,
            remaining - kVerifiedCoordinateReserveSeconds);
        output.attempted = true;
        const auto linear = solve_linearized_ac_seed(
            data, reference, commitment, 0.499999, std::nullopt,
            false, true, lp_time_limit, true, false, {}, false,
            options.linear_economic_voltage_trust_radius,
            options.linear_economic_angle_trust_radius, true);
        nlohmann::json round_json = {
            {"source", "linearized_exact_economic"},
            {"round", economic_round},
            {"remaining_seconds_at_start", remaining},
            {"linear_model", linear.to_json(false)},
            {"candidates", nlohmann::json::array()},
        };
        if (!linear.success) {
            round_json["status"] = "linear_economic_solve_failed";
            output.rounds.push_back(std::move(round_json));
            break;
        }

        SolveResult best = output.selected;
        ValidationReport best_validation = output.selected_validation;
        double best_fraction = 0.0;
        const std::array<double, 8> fractions{
            1.0, 0.75, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625};
        for (const double fraction : fractions) {
            AcState candidate_state = reference;
            const auto interpolate = [fraction] (
                const std::vector<double>& from,
                const std::vector<double>& to,
                std::vector<double>& target) {
                target.resize(from.size());
                for (std::size_t i = 0; i < from.size(); ++i) {
                    target[i] = from[i] + fraction * (to[i] - from[i]);
                }
            };
            interpolate(reference.vm, linear.state.vm, candidate_state.vm);
            interpolate(reference.va, linear.state.va, candidate_state.va);
            interpolate(reference.pg, linear.state.pg, candidate_state.pg);
            interpolate(reference.qg, linear.state.qg, candidate_state.qg);
            interpolate(
                reference.demand_factor, linear.state.demand_factor,
                candidate_state.demand_factor);

            SolveResult candidate;
            candidate.status = 0;
            candidate.iterations = linear.iterations;
            candidate.wall_seconds = linear.wall_seconds;
            candidate.state = std::move(candidate_state);
            candidate.objective = rebuild_base_state_derived_fields(
                data, commitment, candidate.state);
            const auto candidate_validation = validate_state(
                data, ModelMode::BaseSoft, candidate.state, commitment);
            const bool accepted =
                verified_economic_candidate_improves_incumbent(
                    best, best_validation,
                    candidate, candidate_validation,
                    options.validation_tolerance,
                    options.objective_tolerance);
            round_json["candidates"].push_back({
                {"source", "linear_economic_blend"},
                {"fraction", fraction},
                {"objective", candidate.objective},
                {"accepted", accepted},
                {"validation", candidate_validation.to_json()},
            });
            if (accepted) {
                best = std::move(candidate);
                best_validation = candidate_validation;
                best_fraction = fraction;
            }
        }

        const double elapsed_before_repair = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        if (options.time_limit_seconds - elapsed_before_repair > 0.5) {
            FastPowerFlowOptions repair_options;
            repair_options.minimize_active_balance_slack = true;
            repair_options.minimize_reactive_balance_slack = true;
            repair_options.balance_cleanup_fraction = 0.5;
            repair_options.max_newton_iterations = 30;
            repair_options.max_active_redispatch_passes = 8;
            FastContingencyPowerFlow repair(
                data, linear.state, commitment, repair_options);
            auto repaired = repair.solve_base();
            repaired.solve.objective = rebuild_base_state_derived_fields(
                data, commitment, repaired.solve.state);
            repaired.validation = validate_state(
                data, ModelMode::BaseSoft,
                repaired.solve.state, commitment);
            const bool accepted =
                verified_economic_candidate_improves_incumbent(
                    best, best_validation,
                    repaired.solve, repaired.validation,
                    options.validation_tolerance,
                    options.objective_tolerance);
            round_json["candidates"].push_back({
                {"source", "linear_economic_sparse_repair"},
                {"objective", repaired.solve.objective},
                {"accepted", accepted},
                {"validation", repaired.validation.to_json()},
                {"repair", repaired.to_json()},
            });
            if (accepted) {
                best = std::move(repaired.solve);
                best_validation = repaired.validation;
                best_fraction = -1.0;
            }
        }

        const bool improved = verified_economic_candidate_improves_incumbent(
            output.selected, output.selected_validation,
            best, best_validation,
            options.validation_tolerance,
            options.objective_tolerance);
        round_json["accepted"] = improved;
        round_json["accepted_fraction"] = best_fraction;
        round_json["selected_objective"] = best.objective;
        round_json["selected_validation"] = best_validation.to_json();
        output.rounds.push_back(std::move(round_json));
        if (!improved) {
            break;
        }
        output.selected = std::move(best);
        output.selected_validation = best_validation;
        output.selected_objective = output.selected.objective;
        output.improved = true;
        reference = output.selected.state;
    }

    const std::array<double, 6> voltage_changes{
        -0.01, -0.005, -0.001, 0.001, 0.005, 0.01};
    for (int pass = 1;
         pass <= options.maximum_voltage_coordinate_passes;
         ++pass) {
        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        if (options.time_limit_seconds - elapsed <= 0.25) {
            output.time_limit_reached = true;
            break;
        }
        std::vector<int> candidate_buses(data.buses.size());
        std::iota(candidate_buses.begin(), candidate_buses.end(), 0);
        const int retained_bus_count = std::min(
            options.voltage_coordinate_bus_count,
            static_cast<int>(candidate_buses.size()));
        std::partial_sort(
            candidate_buses.begin(),
            candidate_buses.begin() + retained_bus_count,
            candidate_buses.end(),
            [&](int left, int right) {
                const double left_slack =
                    output.selected.state.q_delta[left];
                const double right_slack =
                    output.selected.state.q_delta[right];
                if (left_slack != right_slack) {
                    return left_slack > right_slack;
                }
                return left < right;
            });
        candidate_buses.resize(retained_bus_count);

        SolveResult best = output.selected;
        ValidationReport best_validation = output.selected_validation;
        int best_bus = -1;
        double best_change = 0.0;
        int evaluated = 0;
        int feasible = 0;
        struct VoltageProposal {
            int bus{};
            double change{};
            double improvement{};
        };
        std::vector<VoltageProposal> proposals;
        for (int bus : candidate_buses) {
            double bus_best_objective = output.selected.objective;
            double bus_best_change = 0.0;
            for (double change : voltage_changes) {
                const double proposed =
                    output.selected.state.vm[bus] + change;
                if (proposed < data.buses[bus].vmin - 1e-12 ||
                    proposed > data.buses[bus].vmax + 1e-12) {
                    continue;
                }
                const double trial_elapsed = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - wall_start).count();
                if (options.time_limit_seconds - trial_elapsed <= 0.1) {
                    output.time_limit_reached = true;
                    break;
                }
                ++evaluated;
                SolveResult trial;
                trial.status = 0;
                trial.state = output.selected.state;
                trial.state.vm[bus] = proposed;
                trial.objective = rebuild_base_state_derived_fields(
                    data, commitment, trial.state, 0.5);
                const auto trial_validation = validate_state(
                    data, ModelMode::BaseSoft,
                    trial.state, commitment);
                if (trial_validation.max_residual <=
                    options.validation_tolerance) {
                    ++feasible;
                }
                if (trial_validation.max_residual <=
                        options.validation_tolerance &&
                    trial.objective > bus_best_objective +
                        options.objective_tolerance) {
                    bus_best_objective = trial.objective;
                    bus_best_change = change;
                }
                if (verified_economic_candidate_improves_incumbent(
                        best, best_validation,
                        trial, trial_validation,
                        options.validation_tolerance,
                        options.objective_tolerance)) {
                    best = std::move(trial);
                    best_validation = trial_validation;
                    best_bus = bus;
                    best_change = change;
                }
            }
            if (bus_best_change != 0.0) {
                proposals.push_back({
                    bus,
                    bus_best_change,
                    bus_best_objective - output.selected.objective,
                });
            }
            if (output.time_limit_reached) {
                break;
            }
        }
        std::sort(
            proposals.begin(), proposals.end(),
            [](const VoltageProposal& left,
               const VoltageProposal& right) {
                if (left.improvement != right.improvement) {
                    return left.improvement > right.improvement;
                }
                return left.bus < right.bus;
            });
        std::vector<unsigned char> unavailable(data.buses.size(), 0);
        std::vector<VoltageProposal> batch_proposals;
        for (const auto& proposal : proposals) {
            if (unavailable[proposal.bus] != 0) {
                continue;
            }
            batch_proposals.push_back(proposal);
            unavailable[proposal.bus] = 1;
            for (int branch : data.buses[proposal.bus].branches_from) {
                unavailable[data.branches[branch].to] = 1;
            }
            for (int branch : data.buses[proposal.bus].branches_to) {
                unavailable[data.branches[branch].from] = 1;
            }
        }
        double selected_batch_scale = 0.0;
        nlohmann::json selected_batch_buses = nlohmann::json::array();
        for (const double scale :
             std::array<double, 5>{1.0, 0.75, 0.5, 0.25, 0.125}) {
            if (batch_proposals.size() <= 1) {
                break;
            }
            SolveResult batch;
            batch.status = 0;
            batch.state = output.selected.state;
            for (const auto& proposal : batch_proposals) {
                batch.state.vm[proposal.bus] +=
                    scale * proposal.change;
            }
            batch.objective = rebuild_base_state_derived_fields(
                data, commitment, batch.state, 0.5);
            const auto batch_validation = validate_state(
                data, ModelMode::BaseSoft,
                batch.state, commitment);
            if (!verified_economic_candidate_improves_incumbent(
                    best, best_validation,
                    batch, batch_validation,
                    options.validation_tolerance,
                    options.objective_tolerance)) {
                continue;
            }
            best = std::move(batch);
            best_validation = batch_validation;
            best_bus = -2;
            best_change = scale;
            selected_batch_scale = scale;
            selected_batch_buses = nlohmann::json::array();
            for (const auto& proposal : batch_proposals) {
                selected_batch_buses.push_back({
                    {"bus", data.buses[proposal.bus].source_key},
                    {"unscaled_change", proposal.change},
                    {"individual_improvement", proposal.improvement},
                });
            }
            break;
        }
        nlohmann::json coordinate_round = {
            {"source", "exact_voltage_coordinate"},
            {"pass", pass},
            {"candidate_bus_count", retained_bus_count},
            {"evaluated_candidates", evaluated},
            {"feasible_candidates", feasible},
            {"individual_improving_coordinates", proposals.size()},
            {"nonadjacent_batch_size", batch_proposals.size()},
            {"accepted", best_bus != -1},
        };
        if (best_bus == -1) {
            coordinate_round["status"] = "no_improving_coordinate";
            output.rounds.push_back(std::move(coordinate_round));
            break;
        }
        coordinate_round["selected_bus"] = best_bus >= 0
            ? nlohmann::json(data.buses[best_bus].source_key)
            : nlohmann::json("nonadjacent_batch");
        coordinate_round["selected_change"] = best_change;
        coordinate_round["selected_batch_scale"] = selected_batch_scale;
        coordinate_round["selected_batch_buses"] =
            std::move(selected_batch_buses);
        coordinate_round["selected_objective_before_polish"] =
            best.objective;
        coordinate_round["selected_validation_before_polish"] =
            best_validation.to_json();

        FastPowerFlowOptions polish_options;
        polish_options.minimize_active_balance_slack = true;
        polish_options.minimize_reactive_balance_slack = false;
        polish_options.max_newton_iterations = 12;
        polish_options.max_active_redispatch_passes = 6;
        FastContingencyPowerFlow polish(
            data, best.state, commitment, polish_options);
        auto polished = polish.solve_base();
        polished.solve.objective = rebuild_base_state_derived_fields(
            data, commitment, polished.solve.state, 0.5);
        polished.validation = validate_state(
            data, ModelMode::BaseSoft,
            polished.solve.state, commitment);
        const bool polish_accepted =
            verified_economic_candidate_improves_incumbent(
                best, best_validation,
                polished.solve, polished.validation,
                options.validation_tolerance,
                options.objective_tolerance);
        coordinate_round["polish_accepted"] = polish_accepted;
        coordinate_round["polish_wall_seconds"] = polished.wall_seconds;
        coordinate_round["polish_validation"] =
            polished.validation.to_json();
        if (polish_accepted) {
            best = std::move(polished.solve);
            best_validation = polished.validation;
        }
        coordinate_round["selected_objective"] = best.objective;
        coordinate_round["selected_validation"] =
            best_validation.to_json();
        output.rounds.push_back(std::move(coordinate_round));
        output.selected = std::move(best);
        output.selected_validation = best_validation;
        output.selected_objective = output.selected.objective;
        output.improved = true;
        reference = output.selected.state;
    }
    output.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_start).count();
    output.time_limit_reached = output.time_limit_reached ||
        output.wall_seconds >= options.time_limit_seconds;
    output.selected.wall_seconds = output.wall_seconds;
    output.selected_objective = output.selected.objective;
    return output;
}

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

    if (options.source_status_only) {
        output.commitment = prior;
        output.selected_state = output.base.state;
        output.success = true;
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
