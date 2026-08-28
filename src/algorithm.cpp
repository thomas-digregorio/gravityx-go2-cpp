#include <highs/Highs.h>
#include <Eigen/SparseLU>

#include "gravityx/algorithm.hpp"
#include "gravityx/fast_power_flow.hpp"
#include "gravityx/linearized_ac_seed.hpp"
#include "gravityx/state_io.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <utility>

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

double base_commitment_transition_cost(
    const CaseData& data,
    const std::vector<int>& commitment) {
    if (commitment.size() != data.generators.size()) {
        throw std::runtime_error(
            "base commitment transition-cost size mismatch");
    }
    double cost = 0.0;
    for (int generator = 0;
         generator < static_cast<int>(data.generators.size());
         ++generator) {
        if (commitment[generator] != 0 && commitment[generator] != 1) {
            throw std::runtime_error(
                "base commitment transition-cost status is not binary");
        }
        const auto& source = data.generators[generator];
        if (commitment[generator] > source.status_prev) {
            cost += source.sucost;
        } else if (commitment[generator] < source.status_prev) {
            cost += source.sdcost;
        }
    }
    return cost;
}

nlohmann::json BusInjectionCommitmentResult::to_json(
    bool include_state) const {
    return {
        {"attempted", attempted},
        {"solver_feasible", solver_feasible},
        {"candidate_verified", candidate_verified},
        {"improved", improved},
        {"wall_seconds", wall_seconds},
        {"incumbent_objective", incumbent_objective},
        {"candidate_objective", candidate_objective},
        {"incumbent_transition_cost", incumbent_transition_cost},
        {"candidate_transition_cost", candidate_transition_cost},
        {"incumbent_official_proxy", incumbent_official_proxy},
        {"candidate_official_proxy", candidate_official_proxy},
        {"maximum_bus_active_injection_change",
         maximum_bus_active_injection_change},
        {"maximum_bus_reactive_injection_change",
         maximum_bus_reactive_injection_change},
        {"online_before", online_before},
        {"online_after", online_after},
        {"shutdown_count", shutdown_count},
        {"row_count", row_count},
        {"column_count", column_count},
        {"nonzero_count", nonzero_count},
        {"run_status", run_status},
        {"model_status", model_status},
        {"primal_solution_status", primal_solution_status},
        {"mip_node_count", mip_node_count},
        {"mip_gap", mip_gap},
        {"mip_dual_bound", mip_dual_bound},
        {"status", status},
        {"commitment", commitment},
        {"selected", solve_result_to_json(selected, include_state)},
        {"selected_validation", selected_validation.to_json()},
    };
}

BusInjectionCommitmentResult refine_commitment_preserving_bus_injections(
    const CaseData& data,
    const std::vector<int>& incumbent_commitment,
    const SolveResult& incumbent,
    const BusInjectionCommitmentOptions& options) {
    const auto wall_start = std::chrono::steady_clock::now();
    BusInjectionCommitmentResult output;
    output.attempted = true;
    output.commitment = incumbent_commitment;
    output.selected = incumbent;
    if (incumbent_commitment.size() != data.generators.size() ||
        incumbent.state.pg.size() != data.generators.size() ||
        incumbent.state.qg.size() != data.generators.size() ||
        !std::isfinite(options.time_limit_seconds) ||
        options.time_limit_seconds <= 0.0 ||
        !std::isfinite(options.mip_relative_gap) ||
        options.mip_relative_gap < 0.0 ||
        !std::isfinite(options.validation_tolerance) ||
        options.validation_tolerance < 0.0 ||
        !std::isfinite(options.objective_tolerance) ||
        options.objective_tolerance < 0.0) {
        output.status = "invalid_input";
        output.wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        return output;
    }

    output.selected.objective = rebuild_base_state_derived_fields(
        data, incumbent_commitment, output.selected.state);
    output.selected_validation = validate_state(
        data, ModelMode::BaseSoft, output.selected.state,
        incumbent_commitment);
    output.incumbent_objective = output.selected.objective;
    output.candidate_objective = output.incumbent_objective;
    output.incumbent_transition_cost = base_commitment_transition_cost(
        data, incumbent_commitment);
    output.candidate_transition_cost = output.incumbent_transition_cost;
    output.incumbent_official_proxy = output.incumbent_objective -
        output.incumbent_transition_cost;
    output.candidate_official_proxy = output.incumbent_official_proxy;
    output.online_before = static_cast<int>(std::count(
        incumbent_commitment.begin(), incumbent_commitment.end(), 1));
    output.online_after = output.online_before;
    if (!validated_candidate_is_feasible(
            output.selected, output.selected_validation,
            options.validation_tolerance)) {
        output.status = "incumbent_not_verified";
        output.wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        return output;
    }

    struct SparseRow {
        double lower{-kHighsInf};
        double upper{kHighsInf};
        std::vector<std::pair<HighsInt, double>> entries;
    };
    const int ng = static_cast<int>(data.generators.size());
    const int nb = static_cast<int>(data.buses.size());
    const int u_offset = 0;
    const int q_offset = ng;
    const int headroom_offset = 2 * ng;
    int next_column = 3 * ng;
    std::vector<int> lambda_offset(static_cast<std::size_t>(ng), 0);
    std::vector<std::vector<PwlPoint>> points(static_cast<std::size_t>(ng));
    std::vector<double> p_lower(static_cast<std::size_t>(ng), 0.0);
    std::vector<double> p_upper(static_cast<std::size_t>(ng), 0.0);
    for (int generator = 0; generator < ng; ++generator) {
        const auto& source = data.generators[generator];
        const double previous = source.status_prev == 0
            ? source.pmin : source.pg_prev;
        p_lower[generator] = std::max(
            source.pmin, previous - data.delta_r * source.prdmax);
        p_upper[generator] = std::min(
            source.pmax, previous + data.delta_r * source.prumax);
        if (p_lower[generator] > p_upper[generator] + 1e-12) {
            output.status = "empty_generator_interval";
            output.wall_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - wall_start).count();
            return output;
        }
        points[generator] = active_pwl_points(
            source.cost, source.ncost,
            p_lower[generator], p_upper[generator]);
        lambda_offset[generator] = next_column;
        next_column += static_cast<int>(points[generator].size());
    }
    const int column_count = next_column;
    output.column_count = column_count;
    std::vector<double> lower(
        static_cast<std::size_t>(column_count), 0.0);
    std::vector<double> upper(
        static_cast<std::size_t>(column_count), 1.0);
    std::vector<double> cost(
        static_cast<std::size_t>(column_count), 0.0);
    for (int generator = 0; generator < ng; ++generator) {
        const auto& source = data.generators[generator];
        const bool can_start = source.status_prev == 0 && source.suqual == 1;
        const bool can_stop = source.status_prev == 1 && source.sdqual == 1;
        if ((source.status_prev == 0 && !can_start) ||
            (source.status_prev == 1 && !can_stop)) {
            lower[u_offset + generator] = source.status_prev;
            upper[u_offset + generator] = source.status_prev;
        }
        lower[q_offset + generator] = std::min(0.0, source.qmin);
        upper[q_offset + generator] = std::max(0.0, source.qmax);
        lower[headroom_offset + generator] = 0.0;
        upper[headroom_offset + generator] = std::max(
            0.0, data.delta_r_ctg * source.prumaxctg);
        cost[u_offset + generator] = data.delta * source.oncost +
            (source.status_prev == 0 ? source.sucost : -source.sdcost);
        for (int point = 0;
             point < static_cast<int>(points[generator].size());
             ++point) {
            cost[lambda_offset[generator] + point] =
                data.delta * points[generator][point].cost;
        }
    }

    std::vector<SparseRow> rows;
    rows.reserve(static_cast<std::size_t>(5 * ng + 3 * nb));
    const auto append = [](SparseRow& row, int column, double coefficient) {
        if (std::abs(coefficient) > 1e-14) {
            row.entries.emplace_back(column, coefficient);
        }
    };
    const auto append_pg = [&](SparseRow& row, int generator,
                               double multiplier = 1.0) {
        for (int point = 0;
             point < static_cast<int>(points[generator].size());
             ++point) {
            append(
                row, lambda_offset[generator] + point,
                multiplier * points[generator][point].mw);
        }
    };
    for (int generator = 0; generator < ng; ++generator) {
        const auto& source = data.generators[generator];
        SparseRow lambda_sum;
        lambda_sum.lower = 0.0;
        lambda_sum.upper = 0.0;
        append(lambda_sum, u_offset + generator, -1.0);
        for (int point = 0;
             point < static_cast<int>(points[generator].size());
             ++point) {
            append(lambda_sum, lambda_offset[generator] + point, 1.0);
        }
        rows.push_back(std::move(lambda_sum));

        SparseRow pg_min;
        pg_min.lower = 0.0;
        append_pg(pg_min, generator);
        append(pg_min, u_offset + generator, -p_lower[generator]);
        rows.push_back(std::move(pg_min));
        SparseRow pg_max;
        pg_max.upper = 0.0;
        append_pg(pg_max, generator);
        append(pg_max, u_offset + generator, -p_upper[generator]);
        rows.push_back(std::move(pg_max));

        SparseRow q_min;
        q_min.lower = 0.0;
        append(q_min, q_offset + generator, 1.0);
        append(q_min, u_offset + generator, -source.qmin);
        rows.push_back(std::move(q_min));
        SparseRow q_max;
        q_max.upper = 0.0;
        append(q_max, q_offset + generator, 1.0);
        append(q_max, u_offset + generator, -source.qmax);
        rows.push_back(std::move(q_max));

        SparseRow physical_headroom;
        physical_headroom.upper = 0.0;
        append(
            physical_headroom, headroom_offset + generator, 1.0);
        append_pg(physical_headroom, generator);
        append(
            physical_headroom, u_offset + generator, -source.pmax);
        rows.push_back(std::move(physical_headroom));
        SparseRow ramp_headroom;
        ramp_headroom.upper = 0.0;
        append(ramp_headroom, headroom_offset + generator, 1.0);
        append(
            ramp_headroom, u_offset + generator,
            -std::max(0.0, data.delta_r_ctg * source.prumaxctg));
        rows.push_back(std::move(ramp_headroom));
    }

    std::vector<double> target_p(static_cast<std::size_t>(nb), 0.0);
    std::vector<double> target_q(static_cast<std::size_t>(nb), 0.0);
    for (int generator = 0; generator < ng; ++generator) {
        const int bus = data.generators[generator].bus;
        target_p[bus] += incumbent.state.pg[generator];
        target_q[bus] += incumbent.state.qg[generator];
    }
    for (int bus = 0; bus < nb; ++bus) {
        if (data.buses[bus].generators.empty()) {
            continue;
        }
        SparseRow active_total;
        active_total.lower = target_p[bus];
        active_total.upper = target_p[bus];
        SparseRow reactive_total;
        reactive_total.lower = target_q[bus];
        reactive_total.upper = target_q[bus];
        SparseRow online_count;
        online_count.lower = 1.0;
        bool has_available_generator = false;
        for (int generator : data.buses[bus].generators) {
            append_pg(active_total, generator);
            append(reactive_total, q_offset + generator, 1.0);
            append(online_count, u_offset + generator, 1.0);
            has_available_generator = has_available_generator ||
                upper[u_offset + generator] >= 1.0 - 1e-12;
        }
        rows.push_back(std::move(active_total));
        rows.push_back(std::move(reactive_total));
        if (has_available_generator) {
            rows.push_back(std::move(online_count));
        }
    }

    if (options.enforce_generator_contingency_headroom) {
        std::vector<unsigned char> added(
            static_cast<std::size_t>(ng), 0);
        for (const auto& contingency : data.contingencies) {
            if (contingency.type != ContingencyType::Generator ||
                contingency.component < 0 || contingency.component >= ng ||
                added[contingency.component]) {
                continue;
            }
            const int outaged = contingency.component;
            added[outaged] = 1;
            SparseRow reserve;
            reserve.lower = 0.0;
            for (int generator = 0; generator < ng; ++generator) {
                if (generator != outaged) {
                    append(
                        reserve, headroom_offset + generator, 1.0);
                }
            }
            append_pg(reserve, outaged, -1.0);
            rows.push_back(std::move(reserve));

            double original_maximum_active = 0.0;
            double original_maximum_reactive = 0.0;
            const int bus = data.generators[outaged].bus;
            for (int generator : data.buses[bus].generators) {
                if (incumbent_commitment[generator] == 0) {
                    continue;
                }
                original_maximum_active = std::max(
                    original_maximum_active,
                    std::abs(incumbent.state.pg[generator]));
                original_maximum_reactive = std::max(
                    original_maximum_reactive,
                    std::abs(incumbent.state.qg[generator]));
            }
            SparseRow active_loss_cap;
            active_loss_cap.lower = -original_maximum_active;
            active_loss_cap.upper = original_maximum_active;
            append_pg(active_loss_cap, outaged);
            rows.push_back(std::move(active_loss_cap));
            SparseRow reactive_loss_cap;
            reactive_loss_cap.lower = -original_maximum_reactive;
            reactive_loss_cap.upper = original_maximum_reactive;
            append(reactive_loss_cap, q_offset + outaged, 1.0);
            rows.push_back(std::move(reactive_loss_cap));
        }
    }

    std::vector<double> row_lower;
    std::vector<double> row_upper;
    std::vector<HighsInt> starts;
    std::vector<HighsInt> indices;
    std::vector<double> values;
    row_lower.reserve(rows.size());
    row_upper.reserve(rows.size());
    starts.reserve(rows.size() + 1);
    starts.push_back(0);
    for (const auto& row : rows) {
        row_lower.push_back(row.lower);
        row_upper.push_back(row.upper);
        for (const auto& [column, coefficient] : row.entries) {
            indices.push_back(column);
            values.push_back(coefficient);
        }
        starts.push_back(static_cast<HighsInt>(indices.size()));
    }
    output.row_count = static_cast<int>(rows.size());
    output.nonzero_count = static_cast<int>(indices.size());

    Highs highs;
    const char* highs_log = std::getenv("GRAVITYX_HIGHS_LOG");
    highs.setOptionValue(
        "output_flag", highs_log != nullptr && std::string(highs_log) != "0");
    highs.setOptionValue("threads", 1);
    highs.setOptionValue("presolve", "on");
    highs.setOptionValue("time_limit", options.time_limit_seconds);
    highs.setOptionValue("mip_rel_gap", options.mip_relative_gap);
    highs.setOptionValue("primal_feasibility_tolerance", 1e-8);
    highs.setOptionValue("dual_feasibility_tolerance", 1e-8);
    const bool model_loaded =
        highs.addVars(column_count, lower.data(), upper.data()) ==
            HighsStatus::kOk &&
        highs.changeColsCost(0, column_count - 1, cost.data()) ==
            HighsStatus::kOk &&
        highs.addRows(
            static_cast<HighsInt>(rows.size()), row_lower.data(),
            row_upper.data(), static_cast<HighsInt>(indices.size()),
            starts.data(), indices.data(), values.data()) ==
            HighsStatus::kOk;
    if (!model_loaded) {
        output.status = "model_construction_failed";
        output.wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        return output;
    }
    std::vector<HighsInt> integer_columns(static_cast<std::size_t>(ng));
    std::iota(integer_columns.begin(), integer_columns.end(), HighsInt{0});
    std::vector<HighsVarType> integer_types(
        static_cast<std::size_t>(ng), HighsVarType::kInteger);
    if (highs.changeColsIntegrality(
            ng, integer_columns.data(), integer_types.data()) !=
        HighsStatus::kOk) {
        output.status = "integrality_construction_failed";
        output.wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        return output;
    }

    const HighsStatus run_status = highs.run();
    const HighsModelStatus model_status = highs.getModelStatus();
    const auto& solution = highs.getSolution();
    const auto& info = highs.getInfo();
    output.run_status = static_cast<int>(run_status);
    output.model_status = static_cast<int>(model_status);
    output.primal_solution_status =
        static_cast<int>(info.primal_solution_status);
    output.mip_node_count = static_cast<int>(info.mip_node_count);
    output.mip_gap = info.mip_gap;
    output.mip_dual_bound = info.mip_dual_bound;
    output.status = highs.modelStatusToString(model_status);
    output.solver_feasible = solution.value_valid &&
        solution.col_value.size() == static_cast<std::size_t>(column_count);
    if (!output.solver_feasible) {
        output.wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        return output;
    }

    std::vector<int> candidate_commitment(static_cast<std::size_t>(ng), 0);
    AcState candidate_state = output.selected.state;
    candidate_state.commitment.assign(static_cast<std::size_t>(ng), 0.0);
    candidate_state.startup.assign(static_cast<std::size_t>(ng), 0.0);
    candidate_state.shutdown.assign(static_cast<std::size_t>(ng), 0.0);
    for (int generator = 0; generator < ng; ++generator) {
        const double relaxed_status = solution.col_value[u_offset + generator];
        if (std::abs(relaxed_status - std::round(relaxed_status)) > 1e-6) {
            output.solver_feasible = false;
            output.status = "nonintegral_incumbent";
            output.wall_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - wall_start).count();
            return output;
        }
        candidate_commitment[generator] =
            static_cast<int>(std::round(relaxed_status));
        candidate_state.commitment[generator] =
            static_cast<double>(candidate_commitment[generator]);
        candidate_state.startup[generator] = std::max(
            0, candidate_commitment[generator] -
                data.generators[generator].status_prev);
        candidate_state.shutdown[generator] = std::max(
            0, data.generators[generator].status_prev -
                candidate_commitment[generator]);
        double pg = 0.0;
        for (int point = 0;
             point < static_cast<int>(points[generator].size());
             ++point) {
            pg += points[generator][point].mw *
                solution.col_value[lambda_offset[generator] + point];
        }
        candidate_state.pg[generator] = pg;
        candidate_state.qg[generator] =
            solution.col_value[q_offset + generator];
    }
    SolveResult candidate;
    candidate.status = 0;
    candidate.state = std::move(candidate_state);
    candidate.objective = rebuild_base_state_derived_fields(
        data, candidate_commitment, candidate.state);
    const auto candidate_validation = validate_state(
        data, ModelMode::BaseSoft, candidate.state,
        candidate_commitment);
    output.candidate_verified = validated_candidate_is_feasible(
        candidate, candidate_validation, options.validation_tolerance);
    output.candidate_objective = candidate.objective;
    output.candidate_transition_cost = base_commitment_transition_cost(
        data, candidate_commitment);
    output.candidate_official_proxy = output.candidate_objective -
        output.candidate_transition_cost;
    output.online_after = static_cast<int>(std::count(
        candidate_commitment.begin(), candidate_commitment.end(), 1));
    output.shutdown_count = 0;
    for (int generator = 0; generator < ng; ++generator) {
        output.shutdown_count +=
            incumbent_commitment[generator] == 1 &&
            candidate_commitment[generator] == 0;
    }
    for (int bus = 0; bus < nb; ++bus) {
        double candidate_p = 0.0;
        double candidate_q = 0.0;
        for (int generator : data.buses[bus].generators) {
            candidate_p += candidate.state.pg[generator];
            candidate_q += candidate.state.qg[generator];
        }
        output.maximum_bus_active_injection_change = std::max(
            output.maximum_bus_active_injection_change,
            std::abs(candidate_p - target_p[bus]));
        output.maximum_bus_reactive_injection_change = std::max(
            output.maximum_bus_reactive_injection_change,
            std::abs(candidate_q - target_q[bus]));
    }
    output.improved = output.candidate_verified &&
        output.candidate_official_proxy >
            output.incumbent_official_proxy + options.objective_tolerance;
    if (output.improved) {
        output.commitment = std::move(candidate_commitment);
        output.selected = std::move(candidate);
        output.selected_validation = candidate_validation;
    }
    output.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_start).count();
    output.selected.wall_seconds = output.wall_seconds;
    return output;
}

nlohmann::json ComponentEconomicDispatchResult::to_json(
    bool include_state) const {
    return {
        {"incumbent_verified", incumbent_verified},
        {"attempted", attempted},
        {"solver_feasible", solver_feasible},
        {"solver_optimal", solver_optimal},
        {"improved", improved},
        {"time_limit_reached", time_limit_reached},
        {"primal_start_attempted", primal_start_attempted},
        {"primal_start_accepted", primal_start_accepted},
        {"wall_seconds", wall_seconds},
        {"solver_wall_seconds", solver_wall_seconds},
        {"incumbent_objective", incumbent_objective},
        {"relaxed_market_surplus", relaxed_market_surplus},
        {"selected_objective", selected_objective},
        {"selected_fraction", selected_fraction},
        {"component_count", component_count},
        {"rounds_completed", rounds_completed},
        {"row_count", row_count},
        {"column_count", column_count},
        {"nonzero_count", nonzero_count},
        {"run_status", run_status},
        {"model_status", model_status},
        {"primal_solution_status", primal_solution_status},
        {"primal_start_status", primal_start_status},
        {"simplex_iterations", simplex_iterations},
        {"ipm_iterations", ipm_iterations},
        {"status", status},
        {"trials", trials},
        {"selected", solve_result_to_json(selected, include_state)},
        {"selected_validation", selected_validation.to_json()},
    };
}

ComponentEconomicDispatchResult refine_fixed_commitment_component_economic(
    const CaseData& data,
    const std::vector<int>& commitment,
    const SolveResult& incumbent,
    const ComponentEconomicDispatchOptions& options) {
    const auto wall_start = std::chrono::steady_clock::now();
    ComponentEconomicDispatchResult output;
    output.selected = incumbent;
    if (commitment.size() != data.generators.size() ||
        incumbent.state.pg.size() != data.generators.size() ||
        incumbent.state.qg.size() != data.generators.size() ||
        incumbent.state.demand_factor.size() != data.loads.size() ||
        incumbent.state.vm.size() != data.buses.size() ||
        incumbent.state.va.size() != data.buses.size() ||
        !std::isfinite(options.time_limit_seconds) ||
        options.time_limit_seconds <= 0.0 ||
        !std::isfinite(options.validation_tolerance) ||
        options.validation_tolerance < 0.0 ||
        !std::isfinite(options.objective_tolerance) ||
        options.objective_tolerance < 0.0 ||
        options.maximum_rounds <= 0 ||
        options.maximum_candidate_trials <= 0) {
        output.status = "invalid_input";
        output.wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        return output;
    }
    for (int value : commitment) {
        if (value != 0 && value != 1) {
            output.status = "nonbinary_commitment";
            output.wall_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - wall_start).count();
            return output;
        }
    }

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
        output.status = "incumbent_not_verified";
        output.wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        return output;
    }
    output.attempted = true;

    const int nb = static_cast<int>(data.buses.size());
    const int ng = static_cast<int>(data.generators.size());
    const int nd = static_cast<int>(data.loads.size());

    // Deterministic disjoint-set construction of the in-service connected
    // components.  Each component retains the incumbent's net real-power
    // injection, which in turn retains its aggregate AC loss/shunt burden.
    std::vector<int> parent(static_cast<std::size_t>(nb));
    std::iota(parent.begin(), parent.end(), 0);
    const auto find_root = [&](int start) {
        int root = start;
        while (parent[root] != root) {
            root = parent[root];
        }
        int node = start;
        while (parent[node] != node) {
            const int next = parent[node];
            parent[node] = root;
            node = next;
        }
        return root;
    };
    for (const auto& branch : data.branches) {
        if (!branch.present || branch.status == 0 ||
            branch.from < 0 || branch.from >= nb ||
            branch.to < 0 || branch.to >= nb) {
            continue;
        }
        const int left = find_root(branch.from);
        const int right = find_root(branch.to);
        if (left != right) {
            parent[std::max(left, right)] = std::min(left, right);
        }
    }
    std::vector<int> component_of_bus(static_cast<std::size_t>(nb), -1);
    std::vector<int> root_to_component(static_cast<std::size_t>(nb), -1);
    std::vector<int> component_reference;
    for (int bus = 0; bus < nb; ++bus) {
        const int root = find_root(bus);
        if (root_to_component[root] < 0) {
            root_to_component[root] = output.component_count++;
            component_reference.push_back(bus);
        }
        component_of_bus[bus] = root_to_component[root];
    }

    std::vector<std::vector<PwlPoint>> generator_points(
        static_cast<std::size_t>(ng));
    std::vector<std::vector<PwlPoint>> load_points(
        static_cast<std::size_t>(nd));
    std::vector<int> generator_offset(static_cast<std::size_t>(ng), -1);
    std::vector<int> load_offset(static_cast<std::size_t>(nd), -1);
    std::vector<double> generator_lower(static_cast<std::size_t>(ng), 0.0);
    std::vector<double> generator_upper(static_cast<std::size_t>(ng), 0.0);
    std::vector<double> load_mw_lower(static_cast<std::size_t>(nd), 0.0);
    std::vector<double> load_mw_upper(static_cast<std::size_t>(nd), 0.0);
    int next_column = 0;
    for (int generator = 0; generator < ng; ++generator) {
        if (commitment[generator] == 0) {
            continue;
        }
        const auto& source = data.generators[generator];
        const double previous = source.status_prev == 0
            ? source.pmin : source.pg_prev;
        generator_lower[generator] = std::max(
            source.pmin, previous - data.delta_r * source.prdmax);
        generator_upper[generator] = std::min(
            source.pmax, previous + data.delta_r * source.prumax);
        if (generator_lower[generator] >
                generator_upper[generator] + 1e-12) {
            output.status = "empty_generator_interval";
            output.wall_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - wall_start).count();
            return output;
        }
        generator_points[generator] = active_pwl_points(
            source.cost, source.ncost,
            generator_lower[generator], generator_upper[generator]);
        generator_offset[generator] = next_column;
        next_column += static_cast<int>(generator_points[generator].size());
    }
    for (int load = 0; load < nd; ++load) {
        const auto& source = data.loads[load];
        double factor_lower = source.tmin;
        double factor_upper = source.tmax;
        if (std::abs(source.pd_nominal) > 1e-12) {
            factor_lower = std::max(
                factor_lower,
                (source.pd_prev - source.prdmax * data.delta_r) /
                    source.pd_nominal);
            factor_upper = std::min(
                factor_upper,
                (source.pd_prev + source.prumax * data.delta_r) /
                    source.pd_nominal);
        }
        if (factor_lower > factor_upper + 1e-12) {
            output.status = "empty_load_interval";
            output.wall_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - wall_start).count();
            return output;
        }
        load_mw_lower[load] = source.pd_nominal * factor_lower;
        load_mw_upper[load] = source.pd_nominal * factor_upper;
        if (load_mw_lower[load] > load_mw_upper[load]) {
            std::swap(load_mw_lower[load], load_mw_upper[load]);
        }
        load_points[load] = active_pwl_points(
            source.cost, source.ncost, source.pd_min, source.pd_max);
        load_offset[load] = next_column;
        next_column += static_cast<int>(load_points[load].size());
    }
    output.column_count = next_column;
    if (next_column <= 0) {
        output.status = "empty_relaxation";
        output.wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        return output;
    }

    struct SparseRow {
        double lower{-kHighsInf};
        double upper{kHighsInf};
        std::vector<std::pair<HighsInt, double>> entries;
    };
    std::vector<SparseRow> rows;
    rows.reserve(static_cast<std::size_t>(2 * (ng + nd) +
        output.component_count));
    const auto append = [](SparseRow& row, int column, double coefficient) {
        if (std::abs(coefficient) > 1e-14) {
            row.entries.emplace_back(column, coefficient);
        }
    };
    const auto append_points = [&](SparseRow& row, int offset,
                                   const std::vector<PwlPoint>& points,
                                   double multiplier = 1.0) {
        for (int point = 0; point < static_cast<int>(points.size()); ++point) {
            append(row, offset + point, multiplier * points[point].mw);
        }
    };

    std::vector<double> lower(
        static_cast<std::size_t>(next_column), 0.0);
    std::vector<double> upper(
        static_cast<std::size_t>(next_column), 1.0);
    std::vector<double> cost(
        static_cast<std::size_t>(next_column), 0.0);
    for (int generator = 0; generator < ng; ++generator) {
        if (generator_offset[generator] < 0) {
            continue;
        }
        SparseRow lambda_sum;
        lambda_sum.lower = 1.0;
        lambda_sum.upper = 1.0;
        for (int point = 0;
             point < static_cast<int>(generator_points[generator].size());
             ++point) {
            const int column = generator_offset[generator] + point;
            append(lambda_sum, column, 1.0);
            cost[column] = data.delta *
                generator_points[generator][point].cost;
        }
        rows.push_back(std::move(lambda_sum));
        SparseRow power_bounds;
        power_bounds.lower = generator_lower[generator];
        power_bounds.upper = generator_upper[generator];
        append_points(
            power_bounds, generator_offset[generator],
            generator_points[generator]);
        rows.push_back(std::move(power_bounds));
    }
    for (int load = 0; load < nd; ++load) {
        SparseRow lambda_sum;
        lambda_sum.lower = 1.0;
        lambda_sum.upper = 1.0;
        for (int point = 0;
             point < static_cast<int>(load_points[load].size());
             ++point) {
            const int column = load_offset[load] + point;
            append(lambda_sum, column, 1.0);
            cost[column] = -data.delta * load_points[load][point].cost;
        }
        rows.push_back(std::move(lambda_sum));
        SparseRow power_bounds;
        power_bounds.lower = load_mw_lower[load];
        power_bounds.upper = load_mw_upper[load];
        append_points(power_bounds, load_offset[load], load_points[load]);
        rows.push_back(std::move(power_bounds));
    }

    std::vector<SparseRow> component_rows(
        static_cast<std::size_t>(output.component_count));
    std::vector<double> component_target(
        static_cast<std::size_t>(output.component_count), 0.0);
    for (int generator = 0; generator < ng; ++generator) {
        const int component = component_of_bus[data.generators[generator].bus];
        component_target[component] += incumbent.state.pg[generator];
        if (generator_offset[generator] >= 0) {
            append_points(
                component_rows[component], generator_offset[generator],
                generator_points[generator], 1.0);
        }
    }
    for (int load = 0; load < nd; ++load) {
        const int component = component_of_bus[data.loads[load].bus];
        component_target[component] -= data.loads[load].pd_nominal *
            incumbent.state.demand_factor[load];
        append_points(
            component_rows[component], load_offset[load],
            load_points[load], -1.0);
    }
    for (int component = 0;
         component < output.component_count; ++component) {
        component_rows[component].lower = component_target[component];
        component_rows[component].upper = component_target[component];
        rows.push_back(std::move(component_rows[component]));
    }

    std::vector<double> row_lower;
    std::vector<double> row_upper;
    std::vector<HighsInt> starts;
    std::vector<HighsInt> indices;
    std::vector<double> values;
    row_lower.reserve(rows.size());
    row_upper.reserve(rows.size());
    starts.reserve(rows.size() + 1);
    starts.push_back(0);
    for (const auto& row : rows) {
        row_lower.push_back(row.lower);
        row_upper.push_back(row.upper);
        for (const auto& [column, coefficient] : row.entries) {
            indices.push_back(column);
            values.push_back(coefficient);
        }
        starts.push_back(static_cast<HighsInt>(indices.size()));
    }
    output.row_count = static_cast<int>(rows.size());
    output.nonzero_count = static_cast<int>(indices.size());

    Highs highs;
    const char* highs_log = std::getenv("GRAVITYX_HIGHS_LOG");
    highs.setOptionValue(
        "output_flag", highs_log != nullptr && std::string(highs_log) != "0");
    highs.setOptionValue("threads", 1);
    highs.setOptionValue("presolve", "on");
    highs.setOptionValue("solver", "simplex");
    highs.setOptionValue("primal_feasibility_tolerance", 1e-8);
    highs.setOptionValue("dual_feasibility_tolerance", 1e-8);
    const double elapsed_before_solver = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_start).count();
    const double solver_limit = std::max(
        0.05, std::min(
            0.6 * options.time_limit_seconds,
            options.time_limit_seconds - elapsed_before_solver - 0.05));
    highs.setOptionValue("time_limit", solver_limit);
    const bool model_loaded =
        highs.addVars(next_column, lower.data(), upper.data()) ==
            HighsStatus::kOk &&
        highs.changeColsCost(0, next_column - 1, cost.data()) ==
            HighsStatus::kOk &&
        highs.addRows(
            static_cast<HighsInt>(rows.size()), row_lower.data(),
            row_upper.data(), static_cast<HighsInt>(indices.size()),
            starts.data(), indices.data(), values.data()) ==
            HighsStatus::kOk;
    if (!model_loaded) {
        output.status = "model_construction_failed";
        output.wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        return output;
    }

    const auto interpolation = [](const std::vector<PwlPoint>& points,
                                  double target) {
        std::vector<double> weights(points.size(), 0.0);
        if (target <= points.front().mw) {
            weights.front() = 1.0;
            return weights;
        }
        if (target >= points.back().mw) {
            weights.back() = 1.0;
            return weights;
        }
        for (int point = 0;
             point + 1 < static_cast<int>(points.size()); ++point) {
            if (target > points[point + 1].mw) {
                continue;
            }
            const double width = points[point + 1].mw - points[point].mw;
            if (std::abs(width) <= 1e-14) {
                weights[point] = 1.0;
            } else {
                const double right = (target - points[point].mw) / width;
                weights[point] = 1.0 - right;
                weights[point + 1] = right;
            }
            return weights;
        }
        weights.back() = 1.0;
        return weights;
    };
    std::vector<double> primal_start(
        static_cast<std::size_t>(next_column), 0.0);
    for (int generator = 0; generator < ng; ++generator) {
        if (generator_offset[generator] < 0) {
            continue;
        }
        const auto weights = interpolation(
            generator_points[generator], incumbent.state.pg[generator]);
        std::copy(
            weights.begin(), weights.end(),
            primal_start.begin() + generator_offset[generator]);
    }
    for (int load = 0; load < nd; ++load) {
        const auto weights = interpolation(
            load_points[load], data.loads[load].pd_nominal *
                incumbent.state.demand_factor[load]);
        std::copy(
            weights.begin(), weights.end(),
            primal_start.begin() + load_offset[load]);
    }
    std::vector<HighsInt> start_indices(
        static_cast<std::size_t>(next_column));
    std::iota(start_indices.begin(), start_indices.end(), HighsInt{0});
    output.primal_start_attempted = true;
    const HighsStatus start_status = highs.setSolution(
        next_column, start_indices.data(), primal_start.data());
    output.primal_start_status = static_cast<int>(start_status);
    output.primal_start_accepted = start_status == HighsStatus::kOk;

    double on_cost = 0.0;
    for (int generator = 0; generator < ng; ++generator) {
        if (commitment[generator] != 0) {
            on_cost += data.delta * data.generators[generator].oncost;
        }
    }
    const int component_row_offset =
        static_cast<int>(rows.size()) - output.component_count;
    bool every_solve_optimal = true;
    for (int round = 1; round <= options.maximum_rounds; ++round) {
        const double elapsed_before_round = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        const double remaining_before_round =
            options.time_limit_seconds - elapsed_before_round;
        if (remaining_before_round <= 0.02) {
            output.time_limit_reached = true;
            break;
        }

        const SolveResult round_incumbent = output.selected;
        std::fill(component_target.begin(), component_target.end(), 0.0);
        for (int generator = 0; generator < ng; ++generator) {
            const int component =
                component_of_bus[data.generators[generator].bus];
            component_target[component] +=
                round_incumbent.state.pg[generator];
        }
        for (int load = 0; load < nd; ++load) {
            const int component = component_of_bus[data.loads[load].bus];
            component_target[component] -= data.loads[load].pd_nominal *
                round_incumbent.state.demand_factor[load];
        }
        bool bounds_updated = true;
        for (int component = 0;
             component < output.component_count; ++component) {
            bounds_updated = bounds_updated &&
                highs.changeRowBounds(
                    component_row_offset + component,
                    component_target[component],
                    component_target[component]) == HighsStatus::kOk;
        }
        if (!bounds_updated) {
            output.status = "component_balance_update_failed";
            break;
        }
        highs.setOptionValue(
            "time_limit", std::max(
                0.01, std::min(
                    0.5 * remaining_before_round,
                    remaining_before_round - 0.01)));

        const auto solver_start = std::chrono::steady_clock::now();
        const HighsStatus run_status = highs.run();
        output.solver_wall_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - solver_start).count();
        const HighsModelStatus model_status = highs.getModelStatus();
        const auto& solution = highs.getSolution();
        const auto& info = highs.getInfo();
        output.run_status = static_cast<int>(run_status);
        output.model_status = static_cast<int>(model_status);
        output.primal_solution_status =
            static_cast<int>(info.primal_solution_status);
        output.simplex_iterations +=
            static_cast<int>(info.simplex_iteration_count);
        output.ipm_iterations += static_cast<int>(info.ipm_iteration_count);
        output.status = highs.modelStatusToString(model_status);
        every_solve_optimal = every_solve_optimal &&
            model_status == HighsModelStatus::kOptimal;
        output.solver_optimal = every_solve_optimal;
        output.time_limit_reached = output.time_limit_reached ||
            model_status == HighsModelStatus::kTimeLimit;
        const bool round_solver_feasible = solution.value_valid &&
            solution.col_value.size() ==
                static_cast<std::size_t>(next_column);
        output.solver_feasible = output.solver_feasible ||
            round_solver_feasible;
        if (std::isfinite(info.objective_function_value)) {
            output.relaxed_market_surplus =
                -info.objective_function_value - on_cost;
        }
        if (!round_solver_feasible) {
            break;
        }
        ++output.rounds_completed;

        std::vector<double> target_pg(static_cast<std::size_t>(ng), 0.0);
        std::vector<double> target_demand(static_cast<std::size_t>(nd), 0.0);
        for (int generator = 0; generator < ng; ++generator) {
            if (generator_offset[generator] < 0) {
                continue;
            }
            for (int point = 0;
                 point < static_cast<int>(generator_points[generator].size());
                 ++point) {
                target_pg[generator] +=
                    generator_points[generator][point].mw *
                    solution.col_value[
                        generator_offset[generator] + point];
            }
        }
        for (int load = 0; load < nd; ++load) {
            double mw = 0.0;
            for (int point = 0;
                 point < static_cast<int>(load_points[load].size());
                 ++point) {
                mw += load_points[load][point].mw *
                    solution.col_value[load_offset[load] + point];
            }
            target_demand[load] =
                std::abs(data.loads[load].pd_nominal) > 1e-12
                ? mw / data.loads[load].pd_nominal
                : round_incumbent.state.demand_factor[load];
        }

        // The component relaxation changes injections at thousands of buses.
        // Starting nonlinear repair from the incumbent angles makes those
        // changes look like enormous local imbalances, which previously
        // triggered defensive load/generator projection and discarded most
        // of the economic gain.  Solve the reduced active-angle Newton system
        // once, with one reference angle removed per connected component, so
        // each line-search proposal begins near its own power-flow manifold.
        const auto angle_predictor_start = std::chrono::steady_clock::now();
        std::vector<unsigned char> is_reference(
            static_cast<std::size_t>(nb), 0);
        for (int bus : component_reference) {
            is_reference[bus] = 1;
        }
        std::vector<int> reduced_index(static_cast<std::size_t>(nb), -1);
        int reduced_count = 0;
        for (int bus = 0; bus < nb; ++bus) {
            if (is_reference[bus] == 0) {
                reduced_index[bus] = reduced_count++;
            }
        }
        using SparseMatrix = Eigen::SparseMatrix<
            double, Eigen::ColMajor, int>;
        std::vector<Eigen::Triplet<double, int>> triplets;
        triplets.reserve(4 * data.branches.size());
        const auto add_angle_entry = [&](int row_bus, int column_bus,
                                         double coefficient) {
            if (reduced_index[row_bus] >= 0 &&
                reduced_index[column_bus] >= 0 &&
                std::isfinite(coefficient) &&
                std::abs(coefficient) > 1e-14) {
                triplets.emplace_back(
                    reduced_index[row_bus],
                    reduced_index[column_bus], coefficient);
            }
        };
        for (const auto& branch : data.branches) {
            if (!branch.present || branch.status == 0) {
                continue;
            }
            double from_cross_cos = branch.flow_from_cross_cos;
            double from_cross_sin = branch.flow_from_cross_sin;
            double to_cross_cos = branch.flow_to_cross_cos;
            double to_cross_sin = branch.flow_to_cross_sin;
            if (!branch.flow_coefficients_valid) {
                const double denominator =
                    branch.r * branch.r + branch.x * branch.x;
                const double g = denominator > 1e-20
                    ? branch.r / denominator : 0.0;
                const double b = denominator > 1e-20
                    ? -branch.x / denominator : 0.0;
                // MATPOWER encodes a zero tap ratio as the nominal 1.0 ratio.
                // Parsed production cases normally carry precomputed flow
                // coefficients, but keep the fallback numerically safe and
                // semantically identical for tiny fixtures.
                const double effective_tap =
                    std::abs(branch.tap) > 1e-12 ? branch.tap : 1.0;
                const double tap_squared = effective_tap * effective_tap;
                const double tap_real =
                    effective_tap * std::cos(branch.shift);
                const double tap_imag =
                    effective_tap * std::sin(branch.shift);
                from_cross_cos =
                    (-g * tap_real + b * tap_imag) / tap_squared;
                from_cross_sin =
                    (-b * tap_real - g * tap_imag) / tap_squared;
                to_cross_cos =
                    (-g * tap_real - b * tap_imag) / tap_squared;
                to_cross_sin =
                    (-b * tap_real + g * tap_imag) / tap_squared;
            }
            const int from = branch.from;
            const int to = branch.to;
            const double angle =
                round_incumbent.state.va[from] -
                round_incumbent.state.va[to];
            const double voltage_product =
                round_incumbent.state.vm[from] *
                round_incumbent.state.vm[to];
            const double derivative_from = voltage_product *
                (-from_cross_cos * std::sin(angle) +
                 from_cross_sin * std::cos(angle));
            const double derivative_to = voltage_product *
                (-to_cross_cos * std::sin(angle) -
                 to_cross_sin * std::cos(angle));
            add_angle_entry(from, from, derivative_from);
            add_angle_entry(from, to, -derivative_from);
            add_angle_entry(to, from, derivative_to);
            add_angle_entry(to, to, -derivative_to);
        }
        SparseMatrix active_angle_jacobian(reduced_count, reduced_count);
        active_angle_jacobian.setFromTriplets(
            triplets.begin(), triplets.end());
        active_angle_jacobian.makeCompressed();
        Eigen::VectorXd injection_change = Eigen::VectorXd::Zero(
            reduced_count);
        for (int generator = 0; generator < ng; ++generator) {
            const int bus = data.generators[generator].bus;
            if (reduced_index[bus] >= 0) {
                injection_change[reduced_index[bus]] +=
                    target_pg[generator] -
                    round_incumbent.state.pg[generator];
            }
        }
        for (int load = 0; load < nd; ++load) {
            const int bus = data.loads[load].bus;
            if (reduced_index[bus] >= 0) {
                injection_change[reduced_index[bus]] -=
                    data.loads[load].pd_nominal *
                    (target_demand[load] -
                     round_incumbent.state.demand_factor[load]);
            }
        }
        std::vector<double> angle_correction(
            static_cast<std::size_t>(nb), 0.0);
        // A network made entirely of one-bus components has no free angle
        // variables.  Its zero correction is already the exact predictor.
        bool angle_predictor_success = reduced_count == 0;
        double maximum_angle_correction = 0.0;
        if (reduced_count > 0) {
            Eigen::SparseLU<
                SparseMatrix, Eigen::COLAMDOrdering<int>> angle_factorization;
            angle_factorization.analyzePattern(active_angle_jacobian);
            angle_factorization.factorize(active_angle_jacobian);
            angle_predictor_success =
                angle_factorization.info() == Eigen::Success;
            if (!angle_predictor_success) {
                // Retain the incumbent-angle fallback.  The candidate still
                // passes through nonlinear repair and the independent
                // validator before it can be accepted.
            } else {
            const Eigen::VectorXd reduced_angle =
                angle_factorization.solve(injection_change);
            angle_predictor_success =
                angle_factorization.info() == Eigen::Success &&
                reduced_angle.allFinite();
            if (angle_predictor_success) {
                for (int bus = 0; bus < nb; ++bus) {
                    if (reduced_index[bus] < 0) {
                        continue;
                    }
                    angle_correction[bus] =
                        reduced_angle[reduced_index[bus]];
                    maximum_angle_correction = std::max(
                        maximum_angle_correction,
                        std::abs(angle_correction[bus]));
                }
            }
            }
        }
        const double angle_predictor_seconds =
            std::chrono::duration<double>(
                std::chrono::steady_clock::now() -
                angle_predictor_start).count();

        // The component LP omits network limits, so its full injection move
        // can be far outside the AC-feasible neighborhood.  Start the
        // deterministic backtracking at a small-angle trust boundary rather
        // than spending most of the budget on predictably overlarge trials.
        // Powers of two preserve deterministic trial identities.
        constexpr double kMaximumPredictedAngleStep = 0.01;
        double initial_fraction = 1.0;
        if (angle_predictor_success &&
            maximum_angle_correction > kMaximumPredictedAngleStep) {
            while (initial_fraction * maximum_angle_correction >
                       kMaximumPredictedAngleStep &&
                   initial_fraction > std::ldexp(
                       1.0, 1 - options.maximum_candidate_trials)) {
                initial_fraction *= 0.5;
            }
        }

        bool round_improved = false;
        for (int trial_index = 0;
             trial_index < options.maximum_candidate_trials; ++trial_index) {
            const double elapsed = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - wall_start).count();
            if (elapsed >= options.time_limit_seconds - 0.01) {
                output.time_limit_reached = true;
                break;
            }
            const double fraction =
                initial_fraction * std::ldexp(1.0, -trial_index);
            SolveResult proposal;
            proposal.status = 0;
            proposal.state = round_incumbent.state;
            for (int generator = 0; generator < ng; ++generator) {
                proposal.state.pg[generator] =
                    round_incumbent.state.pg[generator] + fraction *
                    (target_pg[generator] -
                     round_incumbent.state.pg[generator]);
            }
            for (int load = 0; load < nd; ++load) {
                proposal.state.demand_factor[load] =
                    round_incumbent.state.demand_factor[load] + fraction *
                    (target_demand[load] -
                     round_incumbent.state.demand_factor[load]);
            }
            if (angle_predictor_success) {
                for (int bus = 0; bus < nb; ++bus) {
                    proposal.state.va[bus] +=
                        fraction * angle_correction[bus];
                }
            }
            proposal.objective = rebuild_base_state_derived_fields(
                data, commitment, proposal.state);
            auto proposal_validation = validate_state(
                data, ModelMode::BaseSoft, proposal.state, commitment);
            const double raw_proposal_objective = proposal.objective;
            const auto raw_proposal_validation = proposal_validation;
            const double raw_active_balance_slack = std::accumulate(
                proposal.state.p_delta.begin(),
                proposal.state.p_delta.end(), 0.0);
            const double raw_reactive_balance_slack = std::accumulate(
                proposal.state.q_delta.begin(),
                proposal.state.q_delta.end(), 0.0);
            const double raw_thermal_slack = std::accumulate(
                proposal.state.sm_slack.begin(),
                proposal.state.sm_slack.end(), 0.0);
            const double raw_penalty = data.delta *
                (data.p_delta_cost_approx * raw_active_balance_slack +
                 data.q_delta_cost_approx * raw_reactive_balance_slack +
                 data.sm_cost_approx * raw_thermal_slack);
            const double raw_market_surplus_without_penalty =
                raw_proposal_objective + raw_penalty;
            const auto proposed_pg = proposal.state.pg;
            const auto proposed_demand = proposal.state.demand_factor;
            bool nonlinear_repair_attempted = false;
            double nonlinear_repair_seconds = 0.0;
            bool nonlinear_repair_feasible = false;
            bool nonlinear_repair_converged = false;
            bool nonlinear_repair_newton_selected = false;
            int nonlinear_repair_newton_iterations = 0;
            int nonlinear_repair_active_passes = 0;
            int nonlinear_repair_reactive_passes = 0;
            std::string nonlinear_repair_failure_reason;
            if (raw_active_balance_slack + raw_reactive_balance_slack > 1e-8 ||
                !validated_candidate_is_feasible(
                    proposal, proposal_validation,
                    options.validation_tolerance)) {
                nonlinear_repair_attempted = true;
                FastPowerFlowOptions fast_options;
                // The component LP has already selected the economic P/load
                // target.  Ordinary AC Newton should route those injections
                // first; the balance-minimization path locally rewrites the
                // target before attempting the power flow and was observed to
                // discard nearly all of the relaxation gain.
                fast_options.minimize_active_balance_slack = true;
                fast_options.minimize_reactive_balance_slack = true;
                fast_options.skip_balance_cleanup_prepasses = true;
                fast_options.max_newton_iterations = 30;
                fast_options.max_active_redispatch_passes = 12;
                fast_options.max_reactive_limit_passes = 8;
                FastContingencyPowerFlow repair(
                    data, proposal.state, commitment, fast_options);
                auto repaired = repair.solve_base();
                nonlinear_repair_seconds = repaired.wall_seconds;
                nonlinear_repair_feasible = repaired.feasible;
                nonlinear_repair_converged = repaired.converged;
                nonlinear_repair_newton_selected =
                    repaired.newton_candidate_selected;
                nonlinear_repair_newton_iterations =
                    repaired.newton_iterations;
                nonlinear_repair_active_passes =
                    repaired.active_redispatch_passes;
                nonlinear_repair_reactive_passes =
                    repaired.reactive_limit_passes;
                nonlinear_repair_failure_reason = repaired.failure_reason;
                repaired.solve.objective = rebuild_base_state_derived_fields(
                    data, commitment, repaired.solve.state);
                repaired.validation = validate_state(
                    data, ModelMode::BaseSoft,
                    repaired.solve.state, commitment);
                proposal = std::move(repaired.solve);
                proposal_validation = repaired.validation;
            }
            const double candidate_active_balance_slack = std::accumulate(
                proposal.state.p_delta.begin(),
                proposal.state.p_delta.end(), 0.0);
            const double candidate_reactive_balance_slack = std::accumulate(
                proposal.state.q_delta.begin(),
                proposal.state.q_delta.end(), 0.0);
            const double candidate_thermal_slack = std::accumulate(
                proposal.state.sm_slack.begin(),
                proposal.state.sm_slack.end(), 0.0);
            const double candidate_penalty = data.delta *
                (data.p_delta_cost_approx * candidate_active_balance_slack +
                 data.q_delta_cost_approx * candidate_reactive_balance_slack +
                 data.sm_cost_approx * candidate_thermal_slack);
            double generator_movement_from_proposal = 0.0;
            for (int generator = 0; generator < ng; ++generator) {
                generator_movement_from_proposal += std::abs(
                    proposal.state.pg[generator] - proposed_pg[generator]);
            }
            double load_movement_from_proposal = 0.0;
            for (int load = 0; load < nd; ++load) {
                load_movement_from_proposal += std::abs(
                    data.loads[load].pd_nominal *
                    (proposal.state.demand_factor[load] -
                     proposed_demand[load]));
            }
            const bool accepted =
                verified_economic_candidate_improves_incumbent(
                    output.selected, output.selected_validation,
                    proposal, proposal_validation,
                    options.validation_tolerance,
                    options.objective_tolerance);
            output.trials.push_back({
                {"round", round},
                {"trial", trial_index + 1},
                {"fraction", fraction},
                {"angle_predictor_success", angle_predictor_success},
                {"angle_predictor_seconds", angle_predictor_seconds},
                {"maximum_full_angle_correction",
                 maximum_angle_correction},
                {"initial_candidate_fraction", initial_fraction},
                {"raw_objective", raw_proposal_objective},
                {"raw_validation", raw_proposal_validation.to_json()},
                {"raw_active_balance_slack", raw_active_balance_slack},
                {"raw_reactive_balance_slack", raw_reactive_balance_slack},
                {"raw_thermal_slack", raw_thermal_slack},
                {"raw_penalty", raw_penalty},
                {"raw_market_surplus_without_penalty",
                 raw_market_surplus_without_penalty},
                {"nonlinear_repair_attempted", nonlinear_repair_attempted},
                {"nonlinear_repair_seconds", nonlinear_repair_seconds},
                {"nonlinear_repair_feasible", nonlinear_repair_feasible},
                {"nonlinear_repair_converged", nonlinear_repair_converged},
                {"nonlinear_repair_newton_selected",
                 nonlinear_repair_newton_selected},
                {"nonlinear_repair_newton_iterations",
                 nonlinear_repair_newton_iterations},
                {"nonlinear_repair_active_passes",
                 nonlinear_repair_active_passes},
                {"nonlinear_repair_reactive_passes",
                 nonlinear_repair_reactive_passes},
                {"nonlinear_repair_failure_reason",
                 nonlinear_repair_failure_reason},
                {"candidate_objective", proposal.objective},
                {"candidate_penalty", candidate_penalty},
                {"candidate_market_surplus_without_penalty",
                 proposal.objective + candidate_penalty},
                {"generator_movement_from_proposal",
                 generator_movement_from_proposal},
                {"load_movement_from_proposal",
                 load_movement_from_proposal},
                {"candidate_validation", proposal_validation.to_json()},
                {"accepted", accepted},
            });
            if (accepted) {
                output.selected = std::move(proposal);
                output.selected_validation = proposal_validation;
                output.selected_objective = output.selected.objective;
                output.selected_fraction = fraction;
                output.improved = true;
                round_improved = true;
                // This is the largest trusted fraction that passed every
                // acceptance gate.  Recompute the relaxation from the new
                // verified point rather than testing smaller points along the
                // now-stale direction.
                break;
            }
        }
        if (output.time_limit_reached || !round_improved) {
            break;
        }
    }

    output.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_start).count();
    output.time_limit_reached = output.time_limit_reached ||
        output.wall_seconds >= options.time_limit_seconds;
    output.selected.wall_seconds = output.wall_seconds;
    output.selected_objective = output.selected.objective;
    return output;
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
