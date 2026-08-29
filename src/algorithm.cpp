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

struct CommitmentFlowCoefficients {
    double from_g_self{};
    double from_b_self{};
    double to_g_self{};
    double to_b_self{};
    double from_cross_cos{};
    double from_cross_sin{};
    double to_cross_cos{};
    double to_cross_sin{};
};

struct CommitmentBranchAngleDerivative {
    double pf{};
    double qf{};
    double pt{};
    double qt{};
    double pf_vf{};
    double pf_vt{};
    double qf_vf{};
    double qf_vt{};
    double pt_vf{};
    double pt_vt{};
    double qt_vf{};
    double qt_vt{};
};

CommitmentFlowCoefficients commitment_flow_coefficients(
    const Branch& branch) {
    if (branch.flow_coefficients_valid) {
        return {
            branch.flow_from_g_self,
            branch.flow_from_b_self,
            branch.flow_to_g_self,
            branch.flow_to_b_self,
            branch.flow_from_cross_cos,
            branch.flow_from_cross_sin,
            branch.flow_to_cross_cos,
            branch.flow_to_cross_sin,
        };
    }
    const double denominator = branch.r * branch.r + branch.x * branch.x;
    const double g = denominator > 1e-20 ? branch.r / denominator : 0.0;
    const double b = denominator > 1e-20 ? -branch.x / denominator : 0.0;
    if (std::abs(branch.tap) <= 1e-12) {
        throw std::runtime_error(
            "zero tap in component commitment network model: " +
            branch.source_key);
    }
    const double tap_squared = branch.tap * branch.tap;
    const double tap_real = branch.tap * std::cos(branch.shift);
    const double tap_imag = branch.tap * std::sin(branch.shift);
    return {
        branch.transformer
            ? g / tap_squared + branch.g_fr
            : (g + branch.g_fr) / tap_squared,
        branch.transformer
            ? b / tap_squared + branch.b_fr
            : (b + branch.b_fr) / tap_squared,
        g + branch.g_to,
        b + branch.b_to,
        (-g * tap_real + b * tap_imag) / tap_squared,
        (-b * tap_real - g * tap_imag) / tap_squared,
        (-g * tap_real - b * tap_imag) / tap_squared,
        (-b * tap_real + g * tap_imag) / tap_squared,
    };
}

CommitmentBranchAngleDerivative commitment_branch_angle_derivative(
    const Branch& branch,
    const AcState& reference) {
    const auto coefficient = commitment_flow_coefficients(branch);
    const double vf = reference.vm[branch.from];
    const double vt = reference.vm[branch.to];
    const double angle =
        reference.va[branch.from] - reference.va[branch.to];
    const double sine = std::sin(angle);
    const double cosine = std::cos(angle);
    const double voltage_product = vf * vt;
    const double pf_cross =
        coefficient.from_cross_cos * cosine +
        coefficient.from_cross_sin * sine;
    const double qf_cross =
        -coefficient.from_cross_sin * cosine +
        coefficient.from_cross_cos * sine;
    const double pt_cross =
        coefficient.to_cross_cos * cosine -
        coefficient.to_cross_sin * sine;
    const double qt_cross =
        -coefficient.to_cross_sin * cosine -
        coefficient.to_cross_cos * sine;
    return {
        voltage_product *
            (-coefficient.from_cross_cos * sine +
             coefficient.from_cross_sin * cosine),
        voltage_product *
            (coefficient.from_cross_sin * sine +
             coefficient.from_cross_cos * cosine),
        voltage_product *
            (-coefficient.to_cross_cos * sine -
             coefficient.to_cross_sin * cosine),
        voltage_product *
            (coefficient.to_cross_sin * sine -
             coefficient.to_cross_cos * cosine),
        2.0 * coefficient.from_g_self * vf + pf_cross * vt,
        pf_cross * vf,
        -2.0 * coefficient.from_b_self * vf + qf_cross * vt,
        qf_cross * vf,
        pt_cross * vt,
        2.0 * coefficient.to_g_self * vt + pt_cross * vf,
        qt_cross * vt,
        -2.0 * coefficient.to_b_self * vt + qt_cross * vf,
    };
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

nlohmann::json ComponentCommitmentResult::to_json(
    bool include_state) const {
    return {
        {"incumbent_verified", incumbent_verified},
        {"attempted", attempted},
        {"solver_feasible", solver_feasible},
        {"solver_optimal", solver_optimal},
        {"mip_start_attempted", mip_start_attempted},
        {"mip_start_accepted", mip_start_accepted},
        {"compact_pwl_formulation", compact_pwl_formulation},
        {"candidate_repair_attempted", candidate_repair_attempted},
        {"candidate_repair_preserved_dispatch",
         candidate_repair_preserved_dispatch},
        {"candidate_repair_feasible", candidate_repair_feasible},
        {"candidate_repair_converged", candidate_repair_converged},
        {"candidate_verified", candidate_verified},
        {"improved", improved},
        {"time_limit_reached", time_limit_reached},
        {"wall_seconds", wall_seconds},
        {"solver_wall_seconds", solver_wall_seconds},
        {"candidate_repair_wall_seconds", candidate_repair_wall_seconds},
        {"incumbent_objective", incumbent_objective},
        {"raw_candidate_objective", raw_candidate_objective},
        {"candidate_objective", candidate_objective},
        {"incumbent_penalty_slack", incumbent_penalty_slack},
        {"raw_candidate_penalty_slack", raw_candidate_penalty_slack},
        {"candidate_penalty_slack", candidate_penalty_slack},
        {"incumbent_transition_cost", incumbent_transition_cost},
        {"candidate_transition_cost", candidate_transition_cost},
        {"incumbent_official_proxy", incumbent_official_proxy},
        {"candidate_official_proxy", candidate_official_proxy},
        {"solver_objective", solver_objective},
        {"maximum_milp_residual", maximum_milp_residual},
        {"mip_start_maximum_column_violation",
         mip_start_maximum_column_violation},
        {"mip_start_maximum_row_violation",
         mip_start_maximum_row_violation},
        {"candidate_headroom_residual", candidate_headroom_residual},
        {"maximum_milp_residual_identity",
         maximum_milp_residual_identity},
        {"used_near_incumbent_dispatch", used_near_incumbent_dispatch},
        {"component_count", component_count},
        {"bus_active_injection_trust_rows",
         bus_active_injection_trust_rows},
        {"linearized_angle_columns", linearized_angle_columns},
        {"linearized_voltage_columns", linearized_voltage_columns},
        {"linearized_reactive_generation_columns",
         linearized_reactive_generation_columns},
        {"linearized_active_balance_rows",
         linearized_active_balance_rows},
        {"linearized_reactive_balance_rows",
         linearized_reactive_balance_rows},
        {"linearized_reactive_capability_rows",
         linearized_reactive_capability_rows},
        {"linearized_angle_limit_rows", linearized_angle_limit_rows},
        {"linearized_thermal_rows", linearized_thermal_rows},
        {"generator_contingency_headroom_rows",
         generator_contingency_headroom_rows},
        {"candidate_generator_count", candidate_generator_count},
        {"fixed_incumbent_generator_count",
         fixed_incumbent_generator_count},
        {"online_before", online_before},
        {"online_after", online_after},
        {"startup_count", startup_count},
        {"shutdown_count", shutdown_count},
        {"commitment_change_count", commitment_change_count},
        {"row_count", row_count},
        {"column_count", column_count},
        {"nonzero_count", nonzero_count},
        {"run_status", run_status},
        {"model_status", model_status},
        {"primal_solution_status", primal_solution_status},
        {"mip_start_status", mip_start_status},
        {"mip_start_worst_column", mip_start_worst_column},
        {"mip_start_worst_row", mip_start_worst_row},
        {"mip_node_count", mip_node_count},
        {"mip_gap", mip_gap},
        {"mip_dual_bound", mip_dual_bound},
        {"status", status},
        {"candidate_repair_failure_reason",
         candidate_repair_failure_reason},
        {"candidate_commitment", candidate_commitment},
        {"selected_commitment", selected_commitment},
        {"candidate", solve_result_to_json(candidate, include_state)},
        {"raw_candidate_validation", raw_candidate_validation.to_json()},
        {"candidate_validation", candidate_validation.to_json()},
        {"selected", solve_result_to_json(selected, include_state)},
        {"selected_validation", selected_validation.to_json()},
    };
}

ComponentCommitmentResult refine_component_economic_commitment(
    const CaseData& data,
    const std::vector<int>& incumbent_commitment,
    const SolveResult& incumbent,
    const ComponentCommitmentOptions& options) {
    const auto wall_start = std::chrono::steady_clock::now();
    ComponentCommitmentResult output;
    output.selected = incumbent;
    output.selected.status = 0;
    output.selected_commitment = incumbent_commitment;
    output.candidate_commitment = incumbent_commitment;
    const int nb = static_cast<int>(data.buses.size());
    const int ng = static_cast<int>(data.generators.size());
    const int nd = static_cast<int>(data.loads.size());
    if (incumbent_commitment.size() != data.generators.size() ||
        incumbent.state.pg.size() != data.generators.size() ||
        incumbent.state.qg.size() != data.generators.size() ||
        incumbent.state.demand_factor.size() != data.loads.size() ||
        incumbent.state.vm.size() != data.buses.size() ||
        incumbent.state.va.size() != data.buses.size() ||
        !std::isfinite(options.time_limit_seconds) ||
        options.time_limit_seconds <= 0.0 ||
        !std::isfinite(options.mip_relative_gap) ||
        options.mip_relative_gap < 0.0 ||
        !std::isfinite(options.validation_tolerance) ||
        options.validation_tolerance < 0.0 ||
        !std::isfinite(options.objective_tolerance) ||
        options.objective_tolerance < 0.0 ||
        options.maximum_commitment_changes < 0 ||
        options.maximum_candidate_generators < 0 ||
        !std::isfinite(options.bus_active_injection_trust_radius) ||
        !std::isfinite(options.angle_trust_radius) ||
        !std::isfinite(options.voltage_trust_radius) ||
        !std::isfinite(options.thermal_row_utilization_threshold) ||
        !std::isfinite(options.network_movement_penalty) ||
        (options.linearized_active_network &&
         options.angle_trust_radius <= 0.0) ||
        (options.linearized_reactive_network &&
         (!options.linearized_active_network ||
          options.voltage_trust_radius <= 0.0)) ||
        options.thermal_row_utilization_threshold < 0.0 ||
        options.network_movement_penalty < 0.0) {
        output.status = "invalid_input";
        output.wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        return output;
    }
    for (int value : incumbent_commitment) {
        if (value != 0 && value != 1) {
            output.status = "nonbinary_incumbent_commitment";
            output.wall_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - wall_start).count();
            return output;
        }
    }

    output.selected.objective = rebuild_base_state_derived_fields(
        data, incumbent_commitment, output.selected.state);
    output.selected_validation = validate_state(
        data, ModelMode::BaseSoft, output.selected.state,
        incumbent_commitment);
    output.incumbent_verified = validated_candidate_is_feasible(
        output.selected, output.selected_validation,
        options.validation_tolerance);
    output.incumbent_objective = output.selected.objective;
    const auto total_penalty_slack = [](const AcState& state) {
        return std::accumulate(
                   state.p_delta.begin(), state.p_delta.end(), 0.0) +
            std::accumulate(
                   state.q_delta.begin(), state.q_delta.end(), 0.0) +
            std::accumulate(
                   state.sm_slack.begin(), state.sm_slack.end(), 0.0);
    };
    output.incumbent_penalty_slack =
        total_penalty_slack(output.selected.state);
    output.incumbent_transition_cost = base_commitment_transition_cost(
        data, incumbent_commitment);
    output.incumbent_official_proxy = output.incumbent_objective -
        output.incumbent_transition_cost;
    output.candidate_objective = output.incumbent_objective;
    output.candidate_transition_cost = output.incumbent_transition_cost;
    output.candidate_official_proxy = output.incumbent_official_proxy;
    output.online_before = static_cast<int>(std::count(
        incumbent_commitment.begin(), incumbent_commitment.end(), 1));
    output.online_after = output.online_before;
    if (!output.incumbent_verified) {
        output.status = "incumbent_not_verified";
        output.wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        return output;
    }
    output.attempted = true;

    // Deterministic disjoint-set construction of the in-service components.
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
    for (int bus = 0; bus < nb; ++bus) {
        const int root = find_root(bus);
        if (root_to_component[root] < 0) {
            root_to_component[root] = output.component_count++;
        }
        component_of_bus[bus] = root_to_component[root];
    }

    struct SparseRow {
        double lower{-kHighsInf};
        double upper{kHighsInf};
        std::vector<std::pair<HighsInt, double>> entries;
    };
    const int u_offset = 0;
    const int headroom_offset = ng;
    int next_column = 2 * ng;
    std::vector<int> generator_offset(static_cast<std::size_t>(ng), -1);
    std::vector<int> generator_cost_offset(
        static_cast<std::size_t>(ng), -1);
    std::vector<int> load_offset(static_cast<std::size_t>(nd), -1);
    std::vector<int> load_value_offset(
        static_cast<std::size_t>(nd), -1);
    std::vector<std::vector<PwlPoint>> generator_points(
        static_cast<std::size_t>(ng));
    std::vector<std::vector<PwlPoint>> load_points(
        static_cast<std::size_t>(nd));
    std::vector<double> generator_lower(static_cast<std::size_t>(ng), 0.0);
    std::vector<double> generator_upper(static_cast<std::size_t>(ng), 0.0);
    std::vector<double> load_lower(static_cast<std::size_t>(nd), 0.0);
    std::vector<double> load_upper(static_cast<std::size_t>(nd), 0.0);
    for (int generator = 0; generator < ng; ++generator) {
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
            source.cost, source.ncost, generator_lower[generator],
            generator_upper[generator]);
        generator_offset[generator] = next_column;
        if (options.compact_pwl_formulation) {
            ++next_column;
            generator_cost_offset[generator] = next_column++;
        } else {
            next_column += static_cast<int>(
                generator_points[generator].size());
        }
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
        load_lower[load] = source.pd_nominal * factor_lower;
        load_upper[load] = source.pd_nominal * factor_upper;
        if (load_lower[load] > load_upper[load]) {
            std::swap(load_lower[load], load_upper[load]);
        }
        load_points[load] = active_pwl_points(
            source.cost, source.ncost, source.pd_min, source.pd_max);
        load_offset[load] = next_column;
        if (options.compact_pwl_formulation) {
            ++next_column;
            load_value_offset[load] = next_column++;
        } else {
            next_column += static_cast<int>(load_points[load].size());
        }
    }
    if (options.compact_pwl_formulation) {
        const auto slopes_ordered = [](
            const std::vector<PwlPoint>& points, bool convex) {
            double previous = convex ? -kHighsInf : kHighsInf;
            for (int point = 0;
                 point + 1 < static_cast<int>(points.size()); ++point) {
                const double width =
                    points[point + 1].mw - points[point].mw;
                if (width <= 1e-14) {
                    return false;
                }
                const double slope =
                    (points[point + 1].cost - points[point].cost) /
                    width;
                const double tolerance = 1e-9 * std::max({
                    1.0, std::abs(previous), std::abs(slope)});
                if ((convex && slope + tolerance < previous) ||
                    (!convex && slope - tolerance > previous)) {
                    return false;
                }
                previous = slope;
            }
            return true;
        };
        for (int generator = 0; generator < ng; ++generator) {
            if (!slopes_ordered(generator_points[generator], true)) {
                output.status = "compact_generator_curve_not_convex";
                output.wall_seconds = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - wall_start).count();
                return output;
            }
        }
        for (int load = 0; load < nd; ++load) {
            if (!slopes_ordered(load_points[load], false)) {
                output.status = "compact_load_curve_not_concave";
                output.wall_seconds = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - wall_start).count();
                return output;
            }
        }
        output.compact_pwl_formulation = true;
    }
    int angle_offset = -1;
    int voltage_offset = -1;
    int reactive_generation_offset = -1;
    int angle_absolute_offset = -1;
    int voltage_absolute_offset = -1;
    if (options.linearized_active_network) {
        angle_offset = next_column;
        next_column += nb;
        if (options.linearized_reactive_network) {
            voltage_offset = next_column;
            next_column += nb;
            reactive_generation_offset = next_column;
            next_column += ng;
        }
        if (options.network_movement_penalty > 0.0) {
            angle_absolute_offset = next_column;
            next_column += nb;
            if (options.linearized_reactive_network) {
                voltage_absolute_offset = next_column;
                next_column += nb;
            }
        }
        output.linearized_angle_columns = nb;
        if (options.linearized_reactive_network) {
            output.linearized_voltage_columns = nb;
            output.linearized_reactive_generation_columns = ng;
        }
    }
    output.column_count = next_column;

    std::vector<double> lower(
        static_cast<std::size_t>(next_column), 0.0);
    std::vector<double> upper(
        static_cast<std::size_t>(next_column), 1.0);
    std::vector<double> cost(
        static_cast<std::size_t>(next_column), 0.0);
    if (options.linearized_active_network) {
        std::vector<int> component_anchor(
            static_cast<std::size_t>(output.component_count), -1);
        for (int bus = 0; bus < nb; ++bus) {
            const int component = component_of_bus[bus];
            if (component_anchor[component] < 0 ||
                data.buses[bus].type == 3) {
                component_anchor[component] = bus;
            }
        }
        for (int bus = 0; bus < nb; ++bus) {
            lower[angle_offset + bus] = -options.angle_trust_radius;
            upper[angle_offset + bus] = options.angle_trust_radius;
            if (data.buses[bus].type == 3 ||
                component_anchor[component_of_bus[bus]] == bus) {
                lower[angle_offset + bus] = 0.0;
                upper[angle_offset + bus] = 0.0;
            }
            if (options.linearized_reactive_network) {
                lower[voltage_offset + bus] = std::max(
                    -options.voltage_trust_radius,
                    data.buses[bus].vmin - output.selected.state.vm[bus]);
                upper[voltage_offset + bus] = std::min(
                    options.voltage_trust_radius,
                    data.buses[bus].vmax - output.selected.state.vm[bus]);
            }
            if (options.network_movement_penalty > 0.0) {
                upper[angle_absolute_offset + bus] =
                    options.angle_trust_radius;
                cost[angle_absolute_offset + bus] =
                    options.network_movement_penalty;
                if (options.linearized_reactive_network) {
                    upper[voltage_absolute_offset + bus] =
                        options.voltage_trust_radius;
                    cost[voltage_absolute_offset + bus] =
                        options.network_movement_penalty;
                }
            }
        }
        for (int generator = 0;
             options.linearized_reactive_network && generator < ng;
             ++generator) {
            const auto& source = data.generators[generator];
            if (source.qmin > source.qmax + 1e-12) {
                output.status = "empty_generator_reactive_interval";
                output.wall_seconds = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - wall_start).count();
                return output;
            }
            lower[reactive_generation_offset + generator] =
                std::min(0.0, source.qmin);
            upper[reactive_generation_offset + generator] =
                std::max(0.0, source.qmax);
        }
    }
    double transition_objective_constant = 0.0;
    for (int generator = 0; generator < ng; ++generator) {
        const auto& source = data.generators[generator];
        if (source.status_prev == 0 && source.suqual == 0) {
            upper[u_offset + generator] = 0.0;
        }
        if (source.status_prev == 1 && source.sdqual == 0) {
            lower[u_offset + generator] = 1.0;
        }
        upper[headroom_offset + generator] = std::max(
            0.0, data.delta_r_ctg * source.prumaxctg);
        cost[u_offset + generator] = data.delta * source.oncost;
        if (source.status_prev == 0) {
            cost[u_offset + generator] += source.sucost;
        } else {
            cost[u_offset + generator] -= source.sdcost;
            transition_objective_constant += source.sdcost;
        }
        if (options.compact_pwl_formulation) {
            lower[generator_offset[generator]] = std::min(
                0.0, generator_lower[generator]);
            upper[generator_offset[generator]] = std::max(
                0.0, generator_upper[generator]);
            lower[generator_cost_offset[generator]] = -kHighsInf;
            upper[generator_cost_offset[generator]] = kHighsInf;
            cost[generator_cost_offset[generator]] = data.delta;
        } else {
            for (int point = 0;
                 point < static_cast<int>(
                     generator_points[generator].size()); ++point) {
                cost[generator_offset[generator] + point] =
                    data.delta * generator_points[generator][point].cost;
            }
        }
    }
    for (int load = 0; load < nd; ++load) {
        if (options.compact_pwl_formulation) {
            lower[load_offset[load]] = load_lower[load];
            upper[load_offset[load]] = load_upper[load];
            lower[load_value_offset[load]] = -kHighsInf;
            upper[load_value_offset[load]] = kHighsInf;
            cost[load_value_offset[load]] = -data.delta;
        } else {
            for (int point = 0;
                 point < static_cast<int>(load_points[load].size());
                 ++point) {
                cost[load_offset[load] + point] =
                    -data.delta * load_points[load][point].cost;
            }
        }
    }

    if (options.maximum_candidate_generators > 0) {
        struct RankedGenerator {
            int generator{};
            double score{};
        };
        const auto pwl_value = [](const std::vector<PwlPoint>& points,
                                  double power) {
            if (power <= points.front().mw) {
                return points.front().cost;
            }
            for (int point = 0;
                 point + 1 < static_cast<int>(points.size()); ++point) {
                if (power > points[point + 1].mw) {
                    continue;
                }
                const double width =
                    points[point + 1].mw - points[point].mw;
                if (std::abs(width) <= 1e-14) {
                    return points[point].cost;
                }
                const double fraction =
                    (power - points[point].mw) / width;
                return points[point].cost + fraction *
                    (points[point + 1].cost - points[point].cost);
            }
            return points.back().cost;
        };
        std::vector<RankedGenerator> shutdown_candidates;
        std::vector<RankedGenerator> startup_candidates;
        for (int generator = 0; generator < ng; ++generator) {
            const auto& source = data.generators[generator];
            if (incumbent_commitment[generator] == 1 &&
                source.sdqual != 0) {
                // Approximate the source cost released by a shutdown.  The
                // network MILP remains authoritative about whether the lost
                // P/Q can actually be replaced.
                shutdown_candidates.push_back({
                    generator,
                    data.delta * source.oncost +
                        data.delta * pwl_value(
                            generator_points[generator],
                            output.selected.state.pg[generator]) -
                        source.sdcost,
                });
            } else if (incumbent_commitment[generator] == 0 &&
                       source.suqual != 0) {
                const double usable_power = std::max(
                    1e-8, generator_upper[generator] -
                        generator_lower[generator]);
                const double incremental_cost =
                    data.delta * source.oncost + source.sucost +
                    data.delta * (
                        pwl_value(
                            generator_points[generator],
                            generator_upper[generator]) -
                        pwl_value(
                            generator_points[generator],
                            generator_lower[generator]));
                // Larger score means a cheaper startup candidate.
                startup_candidates.push_back({
                    generator, -incremental_cost / usable_power});
            }
        }
        const auto rank_order = [](const RankedGenerator& left,
                                   const RankedGenerator& right) {
            if (left.score != right.score) {
                return left.score > right.score;
            }
            return left.generator < right.generator;
        };
        std::sort(
            shutdown_candidates.begin(), shutdown_candidates.end(),
            rank_order);
        std::sort(
            startup_candidates.begin(), startup_candidates.end(),
            rank_order);
        const int limit = std::min(
            options.maximum_candidate_generators,
            static_cast<int>(shutdown_candidates.size() +
                startup_candidates.size()));
        int startup_limit = startup_candidates.empty()
            ? 0
            : std::min(
                static_cast<int>(startup_candidates.size()),
                std::max(1, limit / 4));
        int shutdown_limit = std::min(
            static_cast<int>(shutdown_candidates.size()),
            limit - startup_limit);
        int remaining = limit - startup_limit - shutdown_limit;
        if (remaining > 0) {
            const int extra_startup = std::min(
                remaining,
                static_cast<int>(startup_candidates.size()) -
                    startup_limit);
            startup_limit += extra_startup;
            remaining -= extra_startup;
        }
        if (remaining > 0) {
            shutdown_limit += std::min(
                remaining,
                static_cast<int>(shutdown_candidates.size()) -
                    shutdown_limit);
        }
        std::vector<unsigned char> candidate(
            static_cast<std::size_t>(ng), 0);
        for (int rank = 0; rank < shutdown_limit; ++rank) {
            candidate[shutdown_candidates[rank].generator] = 1;
        }
        for (int rank = 0; rank < startup_limit; ++rank) {
            candidate[startup_candidates[rank].generator] = 1;
        }
        for (int generator = 0; generator < ng; ++generator) {
            if (candidate[generator]) {
                ++output.candidate_generator_count;
                continue;
            }
            lower[u_offset + generator] = incumbent_commitment[generator];
            upper[u_offset + generator] = incumbent_commitment[generator];
            ++output.fixed_incumbent_generator_count;
        }
    } else {
        output.candidate_generator_count = ng;
    }

    std::vector<SparseRow> rows;
    rows.reserve(static_cast<std::size_t>(
        (options.linearized_active_network ? 7 : 5) * ng + 2 * nd +
        output.component_count + 6 * nb + 4 * data.branches.size() +
        data.contingencies.size()));
    const auto append = [](SparseRow& row, int column, double coefficient) {
        if (std::abs(coefficient) > 1e-14) {
            row.entries.emplace_back(column, coefficient);
        }
    };
    const auto append_pg = [&](SparseRow& row, int generator,
                               double multiplier = 1.0) {
        if (options.compact_pwl_formulation) {
            append(row, generator_offset[generator], multiplier);
        } else {
            for (int point = 0;
                 point < static_cast<int>(
                     generator_points[generator].size()); ++point) {
                append(
                    row, generator_offset[generator] + point,
                    multiplier * generator_points[generator][point].mw);
            }
        }
    };
    const auto append_load = [&](SparseRow& row, int load,
                                 double multiplier = 1.0) {
        if (options.compact_pwl_formulation) {
            append(row, load_offset[load], multiplier);
        } else {
            for (int point = 0;
                 point < static_cast<int>(load_points[load].size());
                 ++point) {
                append(
                    row, load_offset[load] + point,
                    multiplier * load_points[load][point].mw);
            }
        }
    };
    for (int generator = 0; generator < ng; ++generator) {
        const auto& source = data.generators[generator];
        if (!options.compact_pwl_formulation) {
            SparseRow lambda_sum;
            lambda_sum.lower = 0.0;
            lambda_sum.upper = 0.0;
            append(lambda_sum, u_offset + generator, -1.0);
            for (int point = 0;
                 point < static_cast<int>(
                     generator_points[generator].size()); ++point) {
                append(
                    lambda_sum,
                    generator_offset[generator] + point, 1.0);
            }
            rows.push_back(std::move(lambda_sum));
        } else {
            for (int point = 0;
                 point + 1 < static_cast<int>(
                     generator_points[generator].size()); ++point) {
                const auto& left = generator_points[generator][point];
                const auto& right =
                    generator_points[generator][point + 1];
                const double slope =
                    (right.cost - left.cost) /
                    (right.mw - left.mw);
                const double intercept = left.cost - slope * left.mw;
                SparseRow epigraph;
                epigraph.lower = 0.0;
                append(
                    epigraph, generator_cost_offset[generator], 1.0);
                append_pg(epigraph, generator, -slope);
                append(
                    epigraph, u_offset + generator, -intercept);
                rows.push_back(std::move(epigraph));
            }
        }

        // active_pwl_points() deliberately retains a bracketing point just
        // outside each active interval.  These rows are therefore required:
        // lambda_sum == u alone would otherwise permit a committed unit below
        // its exact source PMIN (or above its exact source/ramp PMAX).
        SparseRow pg_min;
        pg_min.lower = 0.0;
        append_pg(pg_min, generator, 1.0);
        append(
            pg_min, u_offset + generator,
            -generator_lower[generator]);
        rows.push_back(std::move(pg_min));
        SparseRow pg_max;
        pg_max.upper = 0.0;
        append_pg(pg_max, generator, 1.0);
        append(
            pg_max, u_offset + generator,
            -generator_upper[generator]);
        rows.push_back(std::move(pg_max));

        SparseRow physical_headroom;
        physical_headroom.upper = 0.0;
        append(physical_headroom, headroom_offset + generator, 1.0);
        append_pg(physical_headroom, generator, 1.0);
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

        if (options.linearized_reactive_network) {
            SparseRow reactive_minimum;
            reactive_minimum.lower = 0.0;
            append(
                reactive_minimum,
                reactive_generation_offset + generator, 1.0);
            append(
                reactive_minimum, u_offset + generator, -source.qmin);
            rows.push_back(std::move(reactive_minimum));
            SparseRow reactive_maximum;
            reactive_maximum.upper = 0.0;
            append(
                reactive_maximum,
                reactive_generation_offset + generator, 1.0);
            append(
                reactive_maximum, u_offset + generator, -source.qmax);
            rows.push_back(std::move(reactive_maximum));
        }
    }
    for (int load = 0; load < nd; ++load) {
        if (!options.compact_pwl_formulation) {
            SparseRow lambda_sum;
            lambda_sum.lower = 1.0;
            lambda_sum.upper = 1.0;
            for (int point = 0;
                 point < static_cast<int>(load_points[load].size());
                 ++point) {
                append(lambda_sum, load_offset[load] + point, 1.0);
            }
            rows.push_back(std::move(lambda_sum));
        } else {
            for (int point = 0;
                 point + 1 < static_cast<int>(load_points[load].size());
                 ++point) {
                const auto& left = load_points[load][point];
                const auto& right = load_points[load][point + 1];
                const double slope =
                    (right.cost - left.cost) /
                    (right.mw - left.mw);
                const double intercept = left.cost - slope * left.mw;
                SparseRow hypograph;
                hypograph.upper = intercept;
                append(hypograph, load_value_offset[load], 1.0);
                append_load(hypograph, load, -slope);
                rows.push_back(std::move(hypograph));
            }
        }
        SparseRow load_bounds;
        load_bounds.lower = load_lower[load];
        load_bounds.upper = load_upper[load];
        append_load(load_bounds, load);
        rows.push_back(std::move(load_bounds));
    }
    if (options.linearized_active_network &&
        options.network_movement_penalty > 0.0) {
        for (int bus = 0; bus < nb; ++bus) {
            SparseRow angle_positive;
            angle_positive.lower = 0.0;
            append(
                angle_positive, angle_absolute_offset + bus, 1.0);
            append(angle_positive, angle_offset + bus, -1.0);
            rows.push_back(std::move(angle_positive));
            SparseRow angle_negative;
            angle_negative.lower = 0.0;
            append(
                angle_negative, angle_absolute_offset + bus, 1.0);
            append(angle_negative, angle_offset + bus, 1.0);
            rows.push_back(std::move(angle_negative));
            if (options.linearized_reactive_network) {
                SparseRow voltage_positive;
                voltage_positive.lower = 0.0;
                append(
                    voltage_positive,
                    voltage_absolute_offset + bus, 1.0);
                append(voltage_positive, voltage_offset + bus, -1.0);
                rows.push_back(std::move(voltage_positive));
                SparseRow voltage_negative;
                voltage_negative.lower = 0.0;
                append(
                    voltage_negative,
                    voltage_absolute_offset + bus, 1.0);
                append(voltage_negative, voltage_offset + bus, 1.0);
                rows.push_back(std::move(voltage_negative));
            }
        }
    }

    std::vector<SparseRow> component_rows(
        static_cast<std::size_t>(output.component_count));
    std::vector<double> component_target(
        static_cast<std::size_t>(output.component_count), 0.0);
    for (int generator = 0; generator < ng; ++generator) {
        const int component =
            component_of_bus[data.generators[generator].bus];
        component_target[component] += incumbent.state.pg[generator];
        append_pg(component_rows[component], generator, 1.0);
    }
    for (int load = 0; load < nd; ++load) {
        const int component = component_of_bus[data.loads[load].bus];
        component_target[component] -= data.loads[load].pd_nominal *
            incumbent.state.demand_factor[load];
        append_load(component_rows[component], load, -1.0);
    }
    if (!options.linearized_active_network) {
        for (int component = 0;
             component < output.component_count; ++component) {
            component_rows[component].lower = component_target[component];
            component_rows[component].upper = component_target[component];
            rows.push_back(std::move(component_rows[component]));
        }
    } else {
        std::vector<CommitmentBranchAngleDerivative> branch_derivative(
            data.branches.size());
        for (int index = 0;
             index < static_cast<int>(data.branches.size()); ++index) {
            const auto& branch = data.branches[index];
            if (!branch.present || branch.status == 0) {
                continue;
            }
            branch_derivative[index] =
                commitment_branch_angle_derivative(
                    branch, output.selected.state);
        }

        // Preserve the independently verified incumbent balance residual at
        // every bus.  Thus the proposal cannot manufacture additional P/Q
        // balance slack: dispatch and responsive load movements must be
        // supported by the first-order AC angle response.
        for (int bus = 0; bus < nb; ++bus) {
            SparseRow active;
            SparseRow reactive;
            double active_target = 0.0;
            double reactive_target = 0.0;
            for (int branch_index : data.buses[bus].branches_from) {
                const auto& branch = data.branches[branch_index];
                if (!branch.present || branch.status == 0) {
                    continue;
                }
                const auto& derivative = branch_derivative[branch_index];
                append(
                    active, angle_offset + branch.from, derivative.pf);
                append(
                    active, angle_offset + branch.to, -derivative.pf);
                if (options.linearized_reactive_network) {
                    append(
                        active, voltage_offset + branch.from,
                        derivative.pf_vf);
                    append(
                        active, voltage_offset + branch.to,
                        derivative.pf_vt);
                    append(
                        reactive, angle_offset + branch.from,
                        derivative.qf);
                    append(
                        reactive, angle_offset + branch.to,
                        -derivative.qf);
                    append(
                        reactive, voltage_offset + branch.from,
                        derivative.qf_vf);
                    append(
                        reactive, voltage_offset + branch.to,
                        derivative.qf_vt);
                }
            }
            for (int branch_index : data.buses[bus].branches_to) {
                const auto& branch = data.branches[branch_index];
                if (!branch.present || branch.status == 0) {
                    continue;
                }
                const auto& derivative = branch_derivative[branch_index];
                append(
                    active, angle_offset + branch.from, derivative.pt);
                append(
                    active, angle_offset + branch.to, -derivative.pt);
                if (options.linearized_reactive_network) {
                    append(
                        active, voltage_offset + branch.from,
                        derivative.pt_vf);
                    append(
                        active, voltage_offset + branch.to,
                        derivative.pt_vt);
                    append(
                        reactive, angle_offset + branch.from,
                        derivative.qt);
                    append(
                        reactive, angle_offset + branch.to,
                        -derivative.qt);
                    append(
                        reactive, voltage_offset + branch.from,
                        derivative.qt_vf);
                    append(
                        reactive, voltage_offset + branch.to,
                        derivative.qt_vt);
                }
            }
            for (int generator : data.buses[bus].generators) {
                append_pg(active, generator, -1.0);
                active_target -= output.selected.state.pg[generator];
                if (options.linearized_reactive_network) {
                    append(
                        reactive,
                        reactive_generation_offset + generator, -1.0);
                    reactive_target -=
                        output.selected.state.qg[generator];
                }
            }
            for (int load : data.buses[bus].loads) {
                append_load(active, load, 1.0);
                active_target += data.loads[load].pd_nominal *
                    output.selected.state.demand_factor[load];
                if (options.linearized_reactive_network &&
                    std::abs(data.loads[load].pd_nominal) > 1e-12) {
                    const double ratio = data.loads[load].qd_nominal /
                        data.loads[load].pd_nominal;
                    append_load(reactive, load, ratio);
                    reactive_target += data.loads[load].qd_nominal *
                        output.selected.state.demand_factor[load];
                }
            }
            if (options.linearized_reactive_network) {
                double shunt_conductance = 0.0;
                double shunt_susceptance = 0.0;
                for (int shunt : data.buses[bus].shunts) {
                    shunt_conductance += data.shunts[shunt].gs;
                    shunt_susceptance += effective_shunt_susceptance(
                        data, output.selected.state, shunt);
                }
                append(
                    active, voltage_offset + bus,
                    2.0 * shunt_conductance *
                        output.selected.state.vm[bus]);
                append(
                    reactive, voltage_offset + bus,
                    -2.0 * shunt_susceptance *
                        output.selected.state.vm[bus]);
            }
            active.lower = active.upper = active_target;
            rows.push_back(std::move(active));
            ++output.linearized_active_balance_rows;
            if (options.linearized_reactive_network) {
                reactive.lower = reactive.upper = reactive_target;
                rows.push_back(std::move(reactive));
                ++output.linearized_reactive_balance_rows;
            }
        }
        if (!options.linearized_reactive_network) {
            for (int bus = 0; bus < nb; ++bus) {
                if (data.buses[bus].generators.empty()) {
                    continue;
                }
                double target_q = 0.0;
                SparseRow lower_capability;
                SparseRow upper_capability;
                for (int generator : data.buses[bus].generators) {
                    target_q += output.selected.state.qg[generator];
                    append(
                        lower_capability, u_offset + generator,
                        data.generators[generator].qmax);
                    append(
                        upper_capability, u_offset + generator,
                        data.generators[generator].qmin);
                }
                lower_capability.lower = target_q;
                upper_capability.upper = target_q;
                rows.push_back(std::move(lower_capability));
                rows.push_back(std::move(upper_capability));
                output.linearized_reactive_capability_rows += 2;
            }
        }

        for (int index = 0;
             index < static_cast<int>(data.branches.size()); ++index) {
            const auto& branch = data.branches[index];
            if (!branch.present || branch.status == 0) {
                continue;
            }
            const double source_delta =
                data.buses[branch.from].va_start -
                data.buses[branch.to].va_start;
            const double reference_delta =
                output.selected.state.va[branch.from] -
                output.selected.state.va[branch.to];
            const bool trust_box_can_reach_angle_limit =
                reference_delta - 2.0 * options.angle_trust_radius <
                    branch.angmin ||
                reference_delta + 2.0 * options.angle_trust_radius >
                    branch.angmax;
            if (source_delta >= branch.angmin &&
                source_delta <= branch.angmax &&
                trust_box_can_reach_angle_limit) {
                SparseRow angle_limit;
                angle_limit.lower = branch.angmin - reference_delta;
                angle_limit.upper = branch.angmax - reference_delta;
                append(angle_limit, angle_offset + branch.from, 1.0);
                append(angle_limit, angle_offset + branch.to, -1.0);
                rows.push_back(std::move(angle_limit));
                ++output.linearized_angle_limit_rows;
            }
            if (branch.rate_a <= 1e-12) {
                continue;
            }
            const double slack = output.selected.state.sm_slack.size() ==
                    data.branches.size()
                ? output.selected.state.sm_slack[index]
                : 0.0;
            const double from_limit = branch.rate_a *
                (branch.transformer
                    ? 1.0 + slack
                    : output.selected.state.vm[branch.from] + slack);
            const double to_limit = branch.rate_a *
                (branch.transformer
                    ? 1.0 + slack
                    : output.selected.state.vm[branch.to] + slack);
            const double from_magnitude = std::hypot(
                output.selected.state.pf[index],
                output.selected.state.qf[index]);
            const double to_magnitude = std::hypot(
                output.selected.state.pt[index],
                output.selected.state.qt[index]);
            const auto add_thermal_row = [&](double magnitude,
                                             double limit,
                                             double p,
                                             double q,
                                             double dp,
                                             double dq,
                                             double dp_vf,
                                             double dq_vf,
                                             double dp_vt,
                                             double dq_vt,
                                             bool from_terminal) {
                if (limit <= 1e-12 ||
                    magnitude / limit <
                        options.thermal_row_utilization_threshold) {
                    return;
                }
                const double derivative = magnitude > 1e-12
                    ? (p * dp + q * dq) / magnitude
                    : 0.0;
                SparseRow thermal;
                thermal.upper = std::max(limit, magnitude) - magnitude;
                append(
                    thermal, angle_offset + branch.from, derivative);
                append(
                    thermal, angle_offset + branch.to, -derivative);
                if (options.linearized_reactive_network) {
                    double voltage_from_derivative = magnitude > 1e-12
                        ? (p * dp_vf + q * dq_vf) / magnitude
                        : 0.0;
                    double voltage_to_derivative = magnitude > 1e-12
                        ? (p * dp_vt + q * dq_vt) / magnitude
                        : 0.0;
                    if (!branch.transformer) {
                        if (from_terminal) {
                            voltage_from_derivative -= branch.rate_a;
                        } else {
                            voltage_to_derivative -= branch.rate_a;
                        }
                    }
                    append(
                        thermal, voltage_offset + branch.from,
                        voltage_from_derivative);
                    append(
                        thermal, voltage_offset + branch.to,
                        voltage_to_derivative);
                }
                rows.push_back(std::move(thermal));
                ++output.linearized_thermal_rows;
            };
            const auto& derivative = branch_derivative[index];
            add_thermal_row(
                from_magnitude, from_limit,
                output.selected.state.pf[index],
                output.selected.state.qf[index],
                derivative.pf, derivative.qf,
                derivative.pf_vf, derivative.qf_vf,
                derivative.pf_vt, derivative.qf_vt, true);
            add_thermal_row(
                to_magnitude, to_limit,
                output.selected.state.pt[index],
                output.selected.state.qt[index],
                derivative.pt, derivative.qt,
                derivative.pt_vf, derivative.qt_vf,
                derivative.pt_vt, derivative.qt_vt, false);
        }
    }

    if (options.bus_active_injection_trust_radius >= 0.0) {
        for (int bus = 0; bus < nb; ++bus) {
            if (data.buses[bus].generators.empty() &&
                data.buses[bus].loads.empty()) {
                continue;
            }
            SparseRow injection;
            double incumbent_injection = 0.0;
            for (int generator : data.buses[bus].generators) {
                incumbent_injection += output.selected.state.pg[generator];
                append_pg(injection, generator, 1.0);
            }
            for (int load : data.buses[bus].loads) {
                incumbent_injection -= data.loads[load].pd_nominal *
                    output.selected.state.demand_factor[load];
                append_load(injection, load, -1.0);
            }
            injection.lower = incumbent_injection -
                options.bus_active_injection_trust_radius;
            injection.upper = incumbent_injection +
                options.bus_active_injection_trust_radius;
            rows.push_back(std::move(injection));
            ++output.bus_active_injection_trust_rows;
        }
    }

    if (options.maximum_commitment_changes > 0) {
        SparseRow hamming_trust;
        int incumbent_online = 0;
        for (int generator = 0; generator < ng; ++generator) {
            if (incumbent_commitment[generator] == 1) {
                ++incumbent_online;
                append(hamming_trust, u_offset + generator, -1.0);
            } else {
                append(hamming_trust, u_offset + generator, 1.0);
            }
        }
        // sum_{u0=1}(1-u) + sum_{u0=0}u <= K
        hamming_trust.upper =
            options.maximum_commitment_changes - incumbent_online;
        rows.push_back(std::move(hamming_trust));
    }

    std::vector<int> headroom_outages;
    if (options.enforce_generator_contingency_headroom) {
        std::vector<unsigned char> added(static_cast<std::size_t>(ng), 0);
        for (const auto& contingency : data.contingencies) {
            if (contingency.type != ContingencyType::Generator ||
                contingency.component < 0 ||
                contingency.component >= ng ||
                added[contingency.component]) {
                continue;
            }
            const int outaged = contingency.component;
            added[outaged] = 1;
            headroom_outages.push_back(outaged);
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
        }
    }
    output.generator_contingency_headroom_rows =
        static_cast<int>(headroom_outages.size());

    // Parallel circuits and co-located devices naturally create duplicate
    // coefficients in the nodal rows.  Canonicalize every sparse row before
    // handing it to HiGHS; this also removes cancellation noise from the
    // first-order network equations.
    for (auto& row : rows) {
        std::sort(
            row.entries.begin(), row.entries.end(),
            [](const auto& left, const auto& right) {
                return left.first < right.first;
            });
        std::vector<std::pair<HighsInt, double>> compressed;
        compressed.reserve(row.entries.size());
        for (const auto& [column, coefficient] : row.entries) {
            if (!compressed.empty() &&
                compressed.back().first == column) {
                compressed.back().second += coefficient;
            } else {
                compressed.emplace_back(column, coefficient);
            }
        }
        compressed.erase(
            std::remove_if(
                compressed.begin(), compressed.end(),
                [](const auto& entry) {
                    // HiGHS deliberately drops matrix entries at or below
                    // its 1e-9 small-matrix threshold and reports a warning.
                    // Remove the same cancellation remnants ourselves so a
                    // harmless warning is not mistaken for construction
                    // failure by the strict status gate below.
                    return std::abs(entry.second) <= 1e-9;
                }),
            compressed.end());
        row.entries = std::move(compressed);
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
    highs.setOptionValue("time_limit", std::max(
        0.05, options.time_limit_seconds - 0.05));
    highs.setOptionValue("mip_rel_gap", options.mip_relative_gap);
    highs.setOptionValue("primal_feasibility_tolerance", 1e-8);
    highs.setOptionValue("dual_feasibility_tolerance", 1e-8);
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

    const auto interpolation = [](const std::vector<PwlPoint>& points,
                                  double target) {
        std::vector<double> weights(points.size(), 0.0);
        if (points.empty()) {
            return weights;
        }
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
            const double width =
                points[point + 1].mw - points[point].mw;
            if (std::abs(width) <= 1e-14) {
                weights[point] = 1.0;
            } else {
                const double right =
                    (target - points[point].mw) / width;
                weights[point] = 1.0 - right;
                weights[point + 1] = right;
            }
            return weights;
        }
        weights.back() = 1.0;
        return weights;
    };
    std::vector<double> mip_start(
        static_cast<std::size_t>(next_column), 0.0);
    for (int generator = 0; generator < ng; ++generator) {
        const int committed = incumbent_commitment[generator];
        mip_start[u_offset + generator] = committed;
        if (committed == 0) {
            continue;
        }
        const double maximum_headroom = std::max(
            0.0, std::min(
                upper[headroom_offset + generator],
                data.generators[generator].pmax -
                    output.selected.state.pg[generator]));
        mip_start[headroom_offset + generator] = maximum_headroom;
        const auto weights = interpolation(
            generator_points[generator],
            output.selected.state.pg[generator]);
        if (options.compact_pwl_formulation) {
            mip_start[generator_offset[generator]] =
                output.selected.state.pg[generator];
            for (int point = 0;
                 point < static_cast<int>(weights.size()); ++point) {
                mip_start[generator_cost_offset[generator]] +=
                    weights[point] *
                    generator_points[generator][point].cost;
            }
        } else {
            std::copy(
                weights.begin(), weights.end(),
                mip_start.begin() + generator_offset[generator]);
        }
        if (options.linearized_reactive_network) {
            mip_start[reactive_generation_offset + generator] =
                output.selected.state.qg[generator];
        }
    }
    for (int load = 0; load < nd; ++load) {
        const auto weights = interpolation(
            load_points[load], data.loads[load].pd_nominal *
                output.selected.state.demand_factor[load]);
        if (options.compact_pwl_formulation) {
            mip_start[load_offset[load]] =
                data.loads[load].pd_nominal *
                output.selected.state.demand_factor[load];
            for (int point = 0;
                 point < static_cast<int>(weights.size()); ++point) {
                mip_start[load_value_offset[load]] +=
                    weights[point] * load_points[load][point].cost;
            }
        } else {
            std::copy(
                weights.begin(), weights.end(),
                mip_start.begin() + load_offset[load]);
        }
    }
    for (int column = 0; column < next_column; ++column) {
        const double violation = std::max({
            0.0,
            lower[column] - mip_start[column],
            mip_start[column] - upper[column],
        });
        if (violation > output.mip_start_maximum_column_violation) {
            output.mip_start_maximum_column_violation = violation;
            output.mip_start_worst_column = column;
        }
    }
    for (int row_index = 0;
         row_index < static_cast<int>(rows.size()); ++row_index) {
        double activity = 0.0;
        for (const auto& [column, coefficient] : rows[row_index].entries) {
            activity += coefficient * mip_start[column];
        }
        const double violation = std::max({
            0.0,
            rows[row_index].lower - activity,
            activity - rows[row_index].upper,
        });
        if (violation > output.mip_start_maximum_row_violation) {
            output.mip_start_maximum_row_violation = violation;
            output.mip_start_worst_row = row_index;
        }
    }
    std::vector<HighsInt> mip_start_indices(
        static_cast<std::size_t>(next_column));
    std::iota(
        mip_start_indices.begin(), mip_start_indices.end(), HighsInt{0});
    output.mip_start_attempted = true;
    const HighsStatus mip_start_status = highs.setSolution(
        next_column, mip_start_indices.data(), mip_start.data());
    output.mip_start_status = static_cast<int>(mip_start_status);
    output.mip_start_accepted =
        mip_start_status == HighsStatus::kOk;

    const auto solver_start = std::chrono::steady_clock::now();
    const HighsStatus run_status = highs.run();
    output.solver_wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - solver_start).count();
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
    output.solver_objective = std::isfinite(info.objective_function_value)
        ? info.objective_function_value + transition_objective_constant
        : info.objective_function_value;
    output.status = highs.modelStatusToString(model_status);
    output.solver_optimal = model_status == HighsModelStatus::kOptimal;
    output.time_limit_reached = model_status == HighsModelStatus::kTimeLimit;
    output.solver_feasible = solution.value_valid &&
        solution.col_value.size() == static_cast<std::size_t>(next_column);
    if (!output.solver_feasible) {
        output.wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        return output;
    }

    output.candidate_commitment.assign(static_cast<std::size_t>(ng), 0);
    std::vector<double> candidate_pg(static_cast<std::size_t>(ng), 0.0);
    std::vector<double> candidate_qg(static_cast<std::size_t>(ng), 0.0);
    std::vector<double> candidate_angle_delta(
        static_cast<std::size_t>(nb), 0.0);
    std::vector<double> candidate_voltage_delta(
        static_cast<std::size_t>(nb), 0.0);
    std::vector<double> candidate_load_mw(static_cast<std::size_t>(nd), 0.0);
    std::vector<double> candidate_headroom(static_cast<std::size_t>(ng), 0.0);
    const auto update_milp_residual =
        [&](double residual, const std::string& identity) {
            if (std::isfinite(residual) &&
                residual > output.maximum_milp_residual) {
                output.maximum_milp_residual = residual;
                output.maximum_milp_residual_identity = identity;
            }
        };
    for (int generator = 0; generator < ng; ++generator) {
        const double u = solution.col_value[u_offset + generator];
        const double rounded = std::round(u);
        update_milp_residual(
            std::abs(u - rounded),
            "generator:" + data.generators[generator].source_key +
                ":integrality");
        output.candidate_commitment[generator] =
            static_cast<int>(rounded);
        candidate_headroom[generator] =
            solution.col_value[headroom_offset + generator];
        if (options.linearized_reactive_network) {
            candidate_qg[generator] = solution.col_value[
                reactive_generation_offset + generator];
        }
        if (options.compact_pwl_formulation) {
            candidate_pg[generator] =
                solution.col_value[generator_offset[generator]];
            const double candidate_cost =
                solution.col_value[generator_cost_offset[generator]];
            for (int point = 0;
                 point + 1 < static_cast<int>(
                     generator_points[generator].size()); ++point) {
                const auto& left = generator_points[generator][point];
                const auto& right =
                    generator_points[generator][point + 1];
                const double slope =
                    (right.cost - left.cost) /
                    (right.mw - left.mw);
                const double intercept = left.cost - slope * left.mw;
                update_milp_residual(
                    std::max(
                        0.0, slope * candidate_pg[generator] +
                            intercept * rounded - candidate_cost),
                    "generator:" +
                        data.generators[generator].source_key +
                        ":cost_epigraph");
            }
        } else {
            double lambda_sum = 0.0;
            for (int point = 0;
                 point < static_cast<int>(
                     generator_points[generator].size()); ++point) {
                const double weight = solution.col_value[
                    generator_offset[generator] + point];
                lambda_sum += weight;
                candidate_pg[generator] +=
                    generator_points[generator][point].mw * weight;
            }
            update_milp_residual(
                std::abs(lambda_sum - rounded),
                "generator:" + data.generators[generator].source_key +
                    ":lambda_sum");
        }
        const auto& source = data.generators[generator];
        if ((source.status_prev == 0 && source.suqual == 0 && rounded != 0.0) ||
            (source.status_prev == 1 && source.sdqual == 0 && rounded != 1.0)) {
            update_milp_residual(
                1.0, "generator:" + source.source_key + ":eligibility");
        }
        update_milp_residual(
            std::max(0.0, generator_lower[generator] * rounded -
                candidate_pg[generator]),
            "generator:" + source.source_key + ":pmin");
        update_milp_residual(
            std::max(0.0, candidate_pg[generator] -
                generator_upper[generator] * rounded),
            "generator:" + source.source_key + ":pmax");
        update_milp_residual(
            std::max(0.0, candidate_headroom[generator] +
                candidate_pg[generator] - source.pmax * rounded),
            "generator:" + source.source_key + ":physical_headroom");
        update_milp_residual(
            std::max(0.0, candidate_headroom[generator] -
                std::max(0.0, data.delta_r_ctg * source.prumaxctg) *
                    rounded),
            "generator:" + source.source_key + ":ramp_headroom");
        update_milp_residual(
            std::max(0.0, -candidate_headroom[generator]),
            "generator:" + source.source_key + ":negative_headroom");
        if (options.linearized_reactive_network) {
            update_milp_residual(
                std::max(
                    0.0, source.qmin * rounded -
                        candidate_qg[generator]),
                "generator:" + source.source_key + ":qmin");
            update_milp_residual(
                std::max(
                    0.0, candidate_qg[generator] -
                        source.qmax * rounded),
                "generator:" + source.source_key + ":qmax");
        }
    }
    if (options.linearized_active_network) {
        for (int bus = 0; bus < nb; ++bus) {
            candidate_angle_delta[bus] =
                solution.col_value[angle_offset + bus];
            update_milp_residual(
                std::max(
                    0.0, std::abs(candidate_angle_delta[bus]) -
                        options.angle_trust_radius),
                "bus:" + data.buses[bus].source_key +
                    ":angle_trust");
            if (data.buses[bus].type == 3) {
                update_milp_residual(
                    std::abs(candidate_angle_delta[bus]),
                    "bus:" + data.buses[bus].source_key +
                        ":reference_angle");
            }
            if (options.linearized_reactive_network) {
                candidate_voltage_delta[bus] =
                    solution.col_value[voltage_offset + bus];
                update_milp_residual(
                    std::max({
                        0.0,
                        data.buses[bus].vmin -
                            output.selected.state.vm[bus] -
                            candidate_voltage_delta[bus],
                        output.selected.state.vm[bus] +
                            candidate_voltage_delta[bus] -
                            data.buses[bus].vmax,
                        std::abs(candidate_voltage_delta[bus]) -
                            options.voltage_trust_radius,
                    }),
                    "bus:" + data.buses[bus].source_key +
                        ":voltage_trust");
            }
        }
    }
    for (int load = 0; load < nd; ++load) {
        if (options.compact_pwl_formulation) {
            candidate_load_mw[load] =
                solution.col_value[load_offset[load]];
            const double candidate_value =
                solution.col_value[load_value_offset[load]];
            for (int point = 0;
                 point + 1 < static_cast<int>(load_points[load].size());
                 ++point) {
                const auto& left = load_points[load][point];
                const auto& right = load_points[load][point + 1];
                const double slope =
                    (right.cost - left.cost) /
                    (right.mw - left.mw);
                const double intercept = left.cost - slope * left.mw;
                update_milp_residual(
                    std::max(
                        0.0, candidate_value -
                            slope * candidate_load_mw[load] - intercept),
                    "load:" + data.loads[load].source_key +
                        ":value_hypograph");
            }
        } else {
            double lambda_sum = 0.0;
            for (int point = 0;
                 point < static_cast<int>(load_points[load].size());
                 ++point) {
                const double weight =
                    solution.col_value[load_offset[load] + point];
                lambda_sum += weight;
                candidate_load_mw[load] +=
                    load_points[load][point].mw * weight;
            }
            update_milp_residual(
                std::abs(lambda_sum - 1.0),
                "load:" + data.loads[load].source_key + ":lambda_sum");
        }
        update_milp_residual(
            std::max(0.0, load_lower[load] - candidate_load_mw[load]),
            "load:" + data.loads[load].source_key + ":lower");
        update_milp_residual(
            std::max(0.0, candidate_load_mw[load] - load_upper[load]),
            "load:" + data.loads[load].source_key + ":upper");
    }
    if (options.linearized_active_network &&
        !options.linearized_reactive_network) {
        for (int bus = 0; bus < nb; ++bus) {
            if (data.buses[bus].generators.empty()) {
                continue;
            }
            double target_q = 0.0;
            double lower_capability = 0.0;
            double upper_capability = 0.0;
            for (int generator : data.buses[bus].generators) {
                target_q += output.selected.state.qg[generator];
                lower_capability += data.generators[generator].qmin *
                    output.candidate_commitment[generator];
                upper_capability += data.generators[generator].qmax *
                    output.candidate_commitment[generator];
            }
            update_milp_residual(
                std::max({
                    0.0,
                    lower_capability - target_q,
                    target_q - upper_capability,
                }),
                "bus:" + data.buses[bus].source_key +
                    ":reactive_capability");
        }
    }
    std::vector<double> decoded_component_balance(
        static_cast<std::size_t>(output.component_count), 0.0);
    for (int generator = 0; generator < ng; ++generator) {
        decoded_component_balance[
            component_of_bus[data.generators[generator].bus]] +=
                candidate_pg[generator];
    }
    for (int load = 0; load < nd; ++load) {
        decoded_component_balance[
            component_of_bus[data.loads[load].bus]] -= candidate_load_mw[load];
    }
    if (!options.linearized_active_network) {
        for (int component = 0;
             component < output.component_count; ++component) {
            update_milp_residual(
                std::abs(decoded_component_balance[component] -
                    component_target[component]),
                "component:" + std::to_string(component) + ":balance");
        }
    } else {
        std::vector<CommitmentBranchAngleDerivative> branch_derivative(
            data.branches.size());
        for (int index = 0;
             index < static_cast<int>(data.branches.size()); ++index) {
            const auto& branch = data.branches[index];
            if (!branch.present || branch.status == 0) {
                continue;
            }
            branch_derivative[index] =
                commitment_branch_angle_derivative(
                    branch, output.selected.state);
        }
        for (int bus = 0; bus < nb; ++bus) {
            double active_residual = 0.0;
            double reactive_residual = 0.0;
            for (int branch_index : data.buses[bus].branches_from) {
                const auto& branch = data.branches[branch_index];
                if (!branch.present || branch.status == 0) {
                    continue;
                }
                const double delta =
                    candidate_angle_delta[branch.from] -
                    candidate_angle_delta[branch.to];
                active_residual +=
                    branch_derivative[branch_index].pf * delta;
                if (options.linearized_reactive_network) {
                    reactive_residual +=
                        branch_derivative[branch_index].qf * delta;
                    active_residual +=
                        branch_derivative[branch_index].pf_vf *
                            candidate_voltage_delta[branch.from] +
                        branch_derivative[branch_index].pf_vt *
                            candidate_voltage_delta[branch.to];
                    reactive_residual +=
                        branch_derivative[branch_index].qf_vf *
                            candidate_voltage_delta[branch.from] +
                        branch_derivative[branch_index].qf_vt *
                            candidate_voltage_delta[branch.to];
                }
            }
            for (int branch_index : data.buses[bus].branches_to) {
                const auto& branch = data.branches[branch_index];
                if (!branch.present || branch.status == 0) {
                    continue;
                }
                const double delta =
                    candidate_angle_delta[branch.from] -
                    candidate_angle_delta[branch.to];
                active_residual +=
                    branch_derivative[branch_index].pt * delta;
                if (options.linearized_reactive_network) {
                    reactive_residual +=
                        branch_derivative[branch_index].qt * delta;
                    active_residual +=
                        branch_derivative[branch_index].pt_vf *
                            candidate_voltage_delta[branch.from] +
                        branch_derivative[branch_index].pt_vt *
                            candidate_voltage_delta[branch.to];
                    reactive_residual +=
                        branch_derivative[branch_index].qt_vf *
                            candidate_voltage_delta[branch.from] +
                        branch_derivative[branch_index].qt_vt *
                            candidate_voltage_delta[branch.to];
                }
            }
            for (int generator : data.buses[bus].generators) {
                active_residual -= candidate_pg[generator] -
                    output.selected.state.pg[generator];
                if (options.linearized_reactive_network) {
                    reactive_residual -= candidate_qg[generator] -
                        output.selected.state.qg[generator];
                }
            }
            for (int load : data.buses[bus].loads) {
                const double incumbent_load =
                    data.loads[load].pd_nominal *
                    output.selected.state.demand_factor[load];
                const double active_change =
                    candidate_load_mw[load] - incumbent_load;
                active_residual += active_change;
                if (options.linearized_reactive_network &&
                    std::abs(data.loads[load].pd_nominal) > 1e-12) {
                    reactive_residual += active_change *
                        data.loads[load].qd_nominal /
                        data.loads[load].pd_nominal;
                }
            }
            if (options.linearized_reactive_network) {
                double shunt_conductance = 0.0;
                double shunt_susceptance = 0.0;
                for (int shunt : data.buses[bus].shunts) {
                    shunt_conductance += data.shunts[shunt].gs;
                    shunt_susceptance += effective_shunt_susceptance(
                        data, output.selected.state, shunt);
                }
                active_residual += 2.0 * shunt_conductance *
                    output.selected.state.vm[bus] *
                    candidate_voltage_delta[bus];
                reactive_residual -= 2.0 * shunt_susceptance *
                    output.selected.state.vm[bus] *
                    candidate_voltage_delta[bus];
            }
            update_milp_residual(
                std::abs(active_residual),
                "bus:" + data.buses[bus].source_key +
                    ":linearized_active_balance");
            if (options.linearized_reactive_network) {
                update_milp_residual(
                    std::abs(reactive_residual),
                    "bus:" + data.buses[bus].source_key +
                        ":linearized_reactive_balance");
            }
        }
        for (int index = 0;
             index < static_cast<int>(data.branches.size()); ++index) {
            const auto& branch = data.branches[index];
            if (!branch.present || branch.status == 0) {
                continue;
            }
            const double source_delta =
                data.buses[branch.from].va_start -
                data.buses[branch.to].va_start;
            const double candidate_delta =
                output.selected.state.va[branch.from] -
                output.selected.state.va[branch.to] +
                candidate_angle_delta[branch.from] -
                candidate_angle_delta[branch.to];
            if (source_delta >= branch.angmin &&
                source_delta <= branch.angmax) {
                update_milp_residual(
                    std::max({
                        0.0,
                        branch.angmin - candidate_delta,
                        candidate_delta - branch.angmax,
                    }),
                    "branch:" + branch.source_key + ":angle_limit");
            }
            if (branch.rate_a <= 1e-12) {
                continue;
            }
            const double slack = output.selected.state.sm_slack.size() ==
                    data.branches.size()
                ? output.selected.state.sm_slack[index]
                : 0.0;
            const double from_limit = branch.rate_a *
                (branch.transformer
                    ? 1.0 + slack
                    : output.selected.state.vm[branch.from] + slack);
            const double to_limit = branch.rate_a *
                (branch.transformer
                    ? 1.0 + slack
                    : output.selected.state.vm[branch.to] + slack);
            const double from_magnitude = std::hypot(
                output.selected.state.pf[index],
                output.selected.state.qf[index]);
            const double to_magnitude = std::hypot(
                output.selected.state.pt[index],
                output.selected.state.qt[index]);
            const double angle_change =
                candidate_angle_delta[branch.from] -
                candidate_angle_delta[branch.to];
            const auto check_terminal = [&](double magnitude,
                                            double limit,
                                            double p,
                                            double q,
                                            double dp,
                                            double dq,
                                            double dp_vf,
                                            double dq_vf,
                                            double dp_vt,
                                            double dq_vt,
                                            bool from_terminal,
                                            const char* terminal) {
                if (limit <= 1e-12 ||
                    magnitude / limit <
                        options.thermal_row_utilization_threshold) {
                    return;
                }
                const double slope = magnitude > 1e-12
                    ? (p * dp + q * dq) / magnitude
                    : 0.0;
                double voltage_response = 0.0;
                if (options.linearized_reactive_network) {
                    double voltage_from_slope = magnitude > 1e-12
                        ? (p * dp_vf + q * dq_vf) / magnitude
                        : 0.0;
                    double voltage_to_slope = magnitude > 1e-12
                        ? (p * dp_vt + q * dq_vt) / magnitude
                        : 0.0;
                    if (!branch.transformer) {
                        if (from_terminal) {
                            voltage_from_slope -= branch.rate_a;
                        } else {
                            voltage_to_slope -= branch.rate_a;
                        }
                    }
                    voltage_response =
                        voltage_from_slope *
                            candidate_voltage_delta[branch.from] +
                        voltage_to_slope *
                            candidate_voltage_delta[branch.to];
                }
                update_milp_residual(
                    std::max(
                        0.0, magnitude + slope * angle_change -
                            std::max(limit, magnitude) +
                            voltage_response),
                    "branch:" + branch.source_key + ":thermal_" +
                        terminal);
            };
            check_terminal(
                from_magnitude, from_limit,
                output.selected.state.pf[index],
                output.selected.state.qf[index],
                branch_derivative[index].pf,
                branch_derivative[index].qf,
                branch_derivative[index].pf_vf,
                branch_derivative[index].qf_vf,
                branch_derivative[index].pf_vt,
                branch_derivative[index].qf_vt,
                true, "from");
            check_terminal(
                to_magnitude, to_limit,
                output.selected.state.pt[index],
                output.selected.state.qt[index],
                branch_derivative[index].pt,
                branch_derivative[index].qt,
                branch_derivative[index].pt_vf,
                branch_derivative[index].qt_vf,
                branch_derivative[index].pt_vt,
                branch_derivative[index].qt_vt,
                false, "to");
        }
    }
    if (options.bus_active_injection_trust_radius >= 0.0) {
        for (int bus = 0; bus < nb; ++bus) {
            double incumbent_injection = 0.0;
            double candidate_injection = 0.0;
            for (int generator : data.buses[bus].generators) {
                incumbent_injection += output.selected.state.pg[generator];
                candidate_injection += candidate_pg[generator];
            }
            for (int load : data.buses[bus].loads) {
                incumbent_injection -= data.loads[load].pd_nominal *
                    output.selected.state.demand_factor[load];
                candidate_injection -= candidate_load_mw[load];
            }
            update_milp_residual(
                std::max(
                    0.0, std::abs(candidate_injection -
                        incumbent_injection) -
                        options.bus_active_injection_trust_radius),
                "bus:" + data.buses[bus].source_key +
                    ":active_injection_trust");
        }
    }
    for (int outaged : headroom_outages) {
        double available = 0.0;
        for (int generator = 0; generator < ng; ++generator) {
            if (generator != outaged) {
                available += candidate_headroom[generator];
            }
        }
        update_milp_residual(
            std::max(0.0, candidate_pg[outaged] - available),
            "generator:" + data.generators[outaged].source_key +
                ":contingency_headroom");
    }
    if (output.maximum_milp_residual > 1e-6) {
        output.solver_feasible = false;
        output.status = "decoded_milp_residual";
        output.wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        return output;
    }

    output.online_after = static_cast<int>(std::count(
        output.candidate_commitment.begin(),
        output.candidate_commitment.end(), 1));
    for (int generator = 0; generator < ng; ++generator) {
        output.startup_count +=
            data.generators[generator].status_prev == 0 &&
            output.candidate_commitment[generator] == 1;
        output.shutdown_count +=
            data.generators[generator].status_prev == 1 &&
            output.candidate_commitment[generator] == 0;
        output.commitment_change_count +=
            output.candidate_commitment[generator] !=
            incumbent_commitment[generator];
    }

    std::vector<double> proposal_pg = candidate_pg;
    std::vector<double> proposal_load_mw = candidate_load_mw;
    if (options.initialize_near_incumbent_dispatch &&
        !options.linearized_active_network) {
        output.used_near_incumbent_dispatch = true;
        for (int generator = 0; generator < ng; ++generator) {
            proposal_pg[generator] = output.candidate_commitment[generator] == 0
                ? 0.0
                : std::clamp(
                    output.selected.state.pg[generator],
                    generator_lower[generator],
                    generator_upper[generator]);
        }
        for (int load = 0; load < nd; ++load) {
            proposal_load_mw[load] = data.loads[load].pd_nominal *
                output.selected.state.demand_factor[load];
        }
        const auto adjust_generation_total =
            [&](const std::vector<int>& generators, double target) {
                double total_lower = 0.0;
                double total_upper = 0.0;
                double current = 0.0;
                for (int generator : generators) {
                    total_lower += generator_lower[generator];
                    total_upper += generator_upper[generator];
                    current += proposal_pg[generator];
                }
                if (target < total_lower - 1e-8 ||
                    target > total_upper + 1e-8) {
                    return false;
                }
                for (int pass = 0;
                     pass < 8 && std::abs(target - current) > 1e-9;
                     ++pass) {
                    const double difference = target - current;
                    double room = 0.0;
                    for (int generator : generators) {
                        room += difference > 0.0
                            ? generator_upper[generator] -
                                proposal_pg[generator]
                            : proposal_pg[generator] -
                                generator_lower[generator];
                    }
                    if (room <= 1e-12) {
                        break;
                    }
                    for (int generator : generators) {
                        const double individual = difference > 0.0
                            ? generator_upper[generator] -
                                proposal_pg[generator]
                            : proposal_pg[generator] -
                                generator_lower[generator];
                        proposal_pg[generator] +=
                            difference * individual / room;
                        proposal_pg[generator] = std::clamp(
                            proposal_pg[generator],
                            generator_lower[generator],
                            generator_upper[generator]);
                    }
                    current = 0.0;
                    for (int generator : generators) {
                        current += proposal_pg[generator];
                    }
                }
                return std::abs(target - current) <= 1e-7;
            };

        // Preserve a bus's aggregate real injection whenever its remaining
        // online units can do so.  Any unavoidable local deficit is then
        // spread only within the same connected component.
        for (int bus = 0; bus < nb; ++bus) {
            std::vector<int> online;
            double target = 0.0;
            for (int generator : data.buses[bus].generators) {
                target += output.selected.state.pg[generator];
                if (output.candidate_commitment[generator] == 1) {
                    online.push_back(generator);
                }
            }
            if (!online.empty()) {
                adjust_generation_total(online, target);
            }
        }
        bool proximity_dispatch_feasible = true;
        for (int component = 0;
             component < output.component_count; ++component) {
            std::vector<int> online;
            double target = 0.0;
            for (int generator = 0; generator < ng; ++generator) {
                if (component_of_bus[data.generators[generator].bus] !=
                    component) {
                    continue;
                }
                target += output.selected.state.pg[generator];
                if (output.candidate_commitment[generator] == 1) {
                    online.push_back(generator);
                }
            }
            proximity_dispatch_feasible = proximity_dispatch_feasible &&
                !online.empty() && adjust_generation_total(online, target);
        }
        if (!proximity_dispatch_feasible) {
            proposal_pg = candidate_pg;
            proposal_load_mw = candidate_load_mw;
            output.used_near_incumbent_dispatch = false;
        }
    }

    output.candidate.status = 0;
    output.candidate.state = output.selected.state;
    output.candidate.state.pg = proposal_pg;
    if (options.linearized_active_network) {
        for (int bus = 0; bus < nb; ++bus) {
            output.candidate.state.va[bus] +=
                candidate_angle_delta[bus];
        }
        if (options.linearized_reactive_network) {
            output.candidate.state.qg = candidate_qg;
            for (int bus = 0; bus < nb; ++bus) {
                output.candidate.state.vm[bus] +=
                    candidate_voltage_delta[bus];
            }
        }
    }
    output.candidate.state.commitment.assign(
        static_cast<std::size_t>(ng), 0.0);
    output.candidate.state.startup.assign(
        static_cast<std::size_t>(ng), 0.0);
    output.candidate.state.shutdown.assign(
        static_cast<std::size_t>(ng), 0.0);
    for (int generator = 0; generator < ng; ++generator) {
        const int status = output.candidate_commitment[generator];
        output.candidate.state.commitment[generator] = status;
        output.candidate.state.startup[generator] = std::max(
            0, status - data.generators[generator].status_prev);
        output.candidate.state.shutdown[generator] = std::max(
            0, data.generators[generator].status_prev - status);
        if (status == 0) {
            output.candidate.state.qg[generator] = 0.0;
        } else {
            output.candidate.state.qg[generator] = std::clamp(
                output.candidate.state.qg[generator],
                data.generators[generator].qmin,
                data.generators[generator].qmax);
        }
    }
    for (int load = 0; load < nd; ++load) {
        output.candidate.state.demand_factor[load] =
            std::abs(data.loads[load].pd_nominal) > 1e-12
            ? proposal_load_mw[load] / data.loads[load].pd_nominal
            : output.selected.state.demand_factor[load];
    }

    // Preserve each incumbent bus's total Q whenever the candidate's online
    // units at that bus can supply it.  This is only a deterministic AC start;
    // the nonlinear repair and validator remain authoritative.
    if (!options.linearized_reactive_network) {
        for (int bus = 0; bus < nb; ++bus) {
            std::vector<int> online;
            double target_q = 0.0;
            double lower_q = 0.0;
            double upper_q = 0.0;
            for (int generator : data.buses[bus].generators) {
                target_q += output.selected.state.qg[generator];
                if (output.candidate_commitment[generator] == 0) {
                    continue;
                }
                online.push_back(generator);
                lower_q += data.generators[generator].qmin;
                upper_q += data.generators[generator].qmax;
            }
            if (online.empty() || target_q < lower_q - 1e-9 ||
                target_q > upper_q + 1e-9) {
                continue;
            }
            double current_q = 0.0;
            for (int generator : online) {
                current_q += output.candidate.state.qg[generator];
            }
            for (int pass = 0; pass < 4 &&
                 std::abs(target_q - current_q) > 1e-10; ++pass) {
                const double difference = target_q - current_q;
                double room = 0.0;
                for (int generator : online) {
                    room += difference > 0.0
                        ? data.generators[generator].qmax -
                            output.candidate.state.qg[generator]
                        : output.candidate.state.qg[generator] -
                            data.generators[generator].qmin;
                }
                if (room <= 1e-12) {
                    break;
                }
                for (int generator : online) {
                    const double individual = difference > 0.0
                        ? data.generators[generator].qmax -
                            output.candidate.state.qg[generator]
                        : output.candidate.state.qg[generator] -
                            data.generators[generator].qmin;
                    output.candidate.state.qg[generator] +=
                        difference * individual / room;
                    output.candidate.state.qg[generator] = std::clamp(
                        output.candidate.state.qg[generator],
                        data.generators[generator].qmin,
                        data.generators[generator].qmax);
                }
                current_q = 0.0;
                for (int generator : online) {
                    current_q += output.candidate.state.qg[generator];
                }
            }
        }
    }

    output.candidate.objective = rebuild_base_state_derived_fields(
        data, output.candidate_commitment, output.candidate.state);
    output.candidate_validation = validate_state(
        data, ModelMode::BaseSoft, output.candidate.state,
        output.candidate_commitment);
    output.raw_candidate_objective = output.candidate.objective;
    output.raw_candidate_validation = output.candidate_validation;
    output.raw_candidate_penalty_slack =
        total_penalty_slack(output.candidate.state);
    output.candidate_verified = validated_candidate_is_feasible(
        output.candidate, output.candidate_validation,
        options.validation_tolerance) &&
        output.raw_candidate_penalty_slack <=
            output.incumbent_penalty_slack + 1e-8;
    const double elapsed_before_repair = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_start).count();
    if (options.repair_candidate && !output.candidate_verified &&
        elapsed_before_repair < options.time_limit_seconds - 0.02) {
        output.candidate_repair_attempted = true;
        output.candidate_repair_preserved_dispatch =
            options.preserve_candidate_dispatch_during_repair;
        FastPowerFlowOptions fast_options;
        fast_options.minimize_active_balance_slack = true;
        fast_options.minimize_reactive_balance_slack = true;
        fast_options.economic_balance_polish = true;
        fast_options.max_economic_balance_polish_iterations = 8;
        fast_options.economic_balance_polish_wall_seconds = std::max(
            0.05, std::min(
                1.0, options.time_limit_seconds -
                    elapsed_before_repair));
        // The default remains the legacy feasibility-first repair.  Economic
        // proposal consumers can instead preserve the component-balanced MILP
        // dispatch; exact Newton/Q-limit recovery and the complete validator
        // still decide whether that candidate may replace the incumbent.
        fast_options.skip_balance_cleanup_prepasses =
            options.preserve_candidate_dispatch_during_repair;
        fast_options.max_newton_iterations = 30;
        fast_options.max_active_redispatch_passes = 12;
        fast_options.max_reactive_limit_passes = 8;
        FastContingencyPowerFlow repair(
            data, output.candidate.state,
            output.candidate_commitment, fast_options);
        auto repaired = repair.solve_base();
        output.candidate_repair_wall_seconds = repaired.wall_seconds;
        output.candidate_repair_feasible = repaired.feasible;
        output.candidate_repair_converged = repaired.converged;
        output.candidate_repair_failure_reason = repaired.failure_reason;
        repaired.solve.status = 0;
        repaired.solve.objective = rebuild_base_state_derived_fields(
            data, output.candidate_commitment, repaired.solve.state);
        repaired.validation = validate_state(
            data, ModelMode::BaseSoft, repaired.solve.state,
            output.candidate_commitment);
        output.candidate = std::move(repaired.solve);
        output.candidate_validation = repaired.validation;
        output.candidate_penalty_slack =
            total_penalty_slack(output.candidate.state);
        output.candidate_verified = validated_candidate_is_feasible(
            output.candidate, output.candidate_validation,
            options.validation_tolerance) &&
            output.candidate_penalty_slack <=
                output.incumbent_penalty_slack + 1e-8;
    }

    if (!output.candidate_repair_attempted) {
        output.candidate_penalty_slack =
            output.raw_candidate_penalty_slack;
    }

    output.candidate_objective = output.candidate.objective;
    output.candidate_transition_cost = base_commitment_transition_cost(
        data, output.candidate_commitment);
    output.candidate_official_proxy = output.candidate_objective -
        output.candidate_transition_cost;
    if (options.enforce_generator_contingency_headroom) {
        std::vector<double> final_headroom(static_cast<std::size_t>(ng), 0.0);
        for (int generator = 0; generator < ng; ++generator) {
            if (output.candidate_commitment[generator] == 0) {
                continue;
            }
            const auto& source = data.generators[generator];
            final_headroom[generator] = std::max(
                0.0, std::min(
                    source.pmax - output.candidate.state.pg[generator],
                    std::max(
                        0.0, data.delta_r_ctg * source.prumaxctg)));
        }
        for (int outaged : headroom_outages) {
            double available = 0.0;
            for (int generator = 0; generator < ng; ++generator) {
                if (generator != outaged) {
                    available += final_headroom[generator];
                }
            }
            output.candidate_headroom_residual = std::max(
                output.candidate_headroom_residual,
                std::max(
                    0.0, output.candidate.state.pg[outaged] - available));
        }
    }
    output.improved = output.candidate_verified &&
        output.candidate_headroom_residual <= 1e-6 &&
        output.candidate_official_proxy >
            output.incumbent_official_proxy + options.objective_tolerance;
    if (output.improved) {
        output.selected = output.candidate;
        output.selected_validation = output.candidate_validation;
        output.selected_commitment = output.candidate_commitment;
    }
    output.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_start).count();
    output.time_limit_reached = output.time_limit_reached ||
        output.wall_seconds >= options.time_limit_seconds;
    output.selected.wall_seconds = output.wall_seconds;
    return output;
}

nlohmann::json GreedyCommitmentSearchResult::to_json(
    bool include_state) const {
    return {
        {"incumbent_verified", incumbent_verified},
        {"proposal_attempted", proposal_attempted},
        {"proposal_feasible", proposal_feasible},
        {"improved", improved},
        {"time_limit_reached", time_limit_reached},
        {"wall_seconds", wall_seconds},
        {"proposal_wall_seconds", proposal_wall_seconds},
        {"candidate_repair_wall_seconds", candidate_repair_wall_seconds},
        {"incumbent_objective", incumbent_objective},
        {"selected_objective", selected_objective},
        {"incumbent_official_proxy", incumbent_official_proxy},
        {"selected_official_proxy", selected_official_proxy},
        {"incumbent_penalty_slack", incumbent_penalty_slack},
        {"selected_penalty_slack", selected_penalty_slack},
        {"selected_headroom_residual", selected_headroom_residual},
        {"proposal_change_count", proposal_change_count},
        {"candidate_pool_size", candidate_pool_size},
        {"rounds_completed", rounds_completed},
        {"candidates_attempted", candidates_attempted},
        {"candidates_precheck_rejected", candidates_precheck_rejected},
        {"candidates_repaired", candidates_repaired},
        {"candidates_verified", candidates_verified},
        {"accepted_moves", accepted_moves},
        {"candidate_order_refreshes", candidate_order_refreshes},
        {"first_improvement_selections", first_improvement_selections},
        {"status", status},
        {"selected_commitment", selected_commitment},
        {"selected", solve_result_to_json(selected, include_state)},
        {"selected_validation", selected_validation.to_json()},
        {"trials", trials},
    };
}

GreedyCommitmentSearchResult refine_greedy_economic_commitment(
    const CaseData& data,
    const std::vector<int>& incumbent_commitment,
    const SolveResult& incumbent,
    const GreedyCommitmentSearchOptions& options) {
    const auto wall_start = std::chrono::steady_clock::now();
    const auto elapsed = [&]() {
        return std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
    };
    GreedyCommitmentSearchResult output;
    output.selected = incumbent;
    output.selected.status = 0;
    output.selected_commitment = incumbent_commitment;
    const int nb = static_cast<int>(data.buses.size());
    const int ng = static_cast<int>(data.generators.size());
    const int nd = static_cast<int>(data.loads.size());
    if (incumbent_commitment.size() != data.generators.size() ||
        incumbent.state.pg.size() != data.generators.size() ||
        incumbent.state.qg.size() != data.generators.size() ||
        incumbent.state.demand_factor.size() != data.loads.size() ||
        incumbent.state.vm.size() != data.buses.size() ||
        incumbent.state.va.size() != data.buses.size() ||
        !std::isfinite(options.time_limit_seconds) ||
        options.time_limit_seconds <= 0.0 ||
        !std::isfinite(options.proposal_time_limit_seconds) ||
        options.proposal_time_limit_seconds <= 0.0 ||
        !std::isfinite(options.validation_tolerance) ||
        options.validation_tolerance < 0.0 ||
        !std::isfinite(options.objective_tolerance) ||
        options.objective_tolerance < 0.0 ||
        options.maximum_rounds <= 0 ||
        options.maximum_candidates_per_round <= 0) {
        output.status = "invalid_input";
        output.wall_seconds = elapsed();
        return output;
    }
    for (int value : incumbent_commitment) {
        if (value != 0 && value != 1) {
            output.status = "nonbinary_incumbent_commitment";
            output.wall_seconds = elapsed();
            return output;
        }
    }
    const auto total_penalty_slack = [](const AcState& state) {
        return std::accumulate(
                   state.p_delta.begin(), state.p_delta.end(), 0.0) +
            std::accumulate(
                   state.q_delta.begin(), state.q_delta.end(), 0.0) +
            std::accumulate(
                   state.sm_slack.begin(), state.sm_slack.end(), 0.0);
    };
    output.selected.objective = rebuild_base_state_derived_fields(
        data, output.selected_commitment, output.selected.state);
    output.selected_validation = validate_state(
        data, ModelMode::BaseSoft, output.selected.state,
        output.selected_commitment);
    output.incumbent_verified = validated_candidate_is_feasible(
        output.selected, output.selected_validation,
        options.validation_tolerance);
    output.incumbent_objective = output.selected.objective;
    output.selected_objective = output.selected.objective;
    output.incumbent_penalty_slack =
        total_penalty_slack(output.selected.state);
    output.selected_penalty_slack = output.incumbent_penalty_slack;
    output.incumbent_official_proxy = output.incumbent_objective -
        base_commitment_transition_cost(data, output.selected_commitment);
    output.selected_official_proxy = output.incumbent_official_proxy;
    if (!output.incumbent_verified) {
        output.status = "incumbent_not_verified";
        output.wall_seconds = elapsed();
        return output;
    }

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
    std::vector<int> root_to_component(static_cast<std::size_t>(nb), -1);
    std::vector<int> component_of_bus(static_cast<std::size_t>(nb), -1);
    int component_count = 0;
    for (int bus = 0; bus < nb; ++bus) {
        const int root = find_root(bus);
        if (root_to_component[root] < 0) {
            root_to_component[root] = component_count++;
        }
        component_of_bus[bus] = root_to_component[root];
    }
    std::vector<std::vector<int>> component_generators(
        static_cast<std::size_t>(component_count));
    for (int generator = 0; generator < ng; ++generator) {
        component_generators[
            component_of_bus[data.generators[generator].bus]].push_back(
                generator);
    }

    std::vector<double> generator_lower(static_cast<std::size_t>(ng), 0.0);
    std::vector<double> generator_upper(static_cast<std::size_t>(ng), 0.0);
    std::vector<std::vector<PwlPoint>> generator_points(
        static_cast<std::size_t>(ng));
    for (int generator = 0; generator < ng; ++generator) {
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
            output.wall_seconds = elapsed();
            return output;
        }
        generator_points[generator] = active_pwl_points(
            source.cost, source.ncost, generator_lower[generator],
            generator_upper[generator]);
    }
    const auto pwl_value = [](const std::vector<PwlPoint>& points,
                              double power) {
        if (power <= points.front().mw) {
            return points.front().cost;
        }
        for (int point = 0;
             point + 1 < static_cast<int>(points.size()); ++point) {
            if (power > points[point + 1].mw) {
                continue;
            }
            const double width =
                points[point + 1].mw - points[point].mw;
            if (std::abs(width) <= 1e-14) {
                return points[point].cost;
            }
            const double fraction =
                (power - points[point].mw) / width;
            return points[point].cost + fraction *
                (points[point + 1].cost - points[point].cost);
        }
        return points.back().cost;
    };
    const auto marginal_cost = [&](int generator, double power) {
        const auto& points = generator_points[generator];
        for (int point = 0;
             point + 1 < static_cast<int>(points.size()); ++point) {
            if (power > points[point + 1].mw + 1e-12) {
                continue;
            }
            const double width = points[point + 1].mw - points[point].mw;
            return std::abs(width) <= 1e-14
                ? 0.0
                : (points[point + 1].cost - points[point].cost) /
                    width;
        }
        const int last = static_cast<int>(points.size()) - 1;
        const double width = points[last].mw - points[last - 1].mw;
        return std::abs(width) <= 1e-14
            ? 0.0
            : (points[last].cost - points[last - 1].cost) / width;
    };

    output.proposal_attempted = true;
    ComponentCommitmentOptions proposal_options;
    proposal_options.time_limit_seconds = std::min(
        options.proposal_time_limit_seconds,
        std::max(0.05, options.time_limit_seconds - elapsed() - 0.05));
    proposal_options.mip_relative_gap = 1e-3;
    proposal_options.validation_tolerance = options.validation_tolerance;
    proposal_options.objective_tolerance = options.objective_tolerance;
    proposal_options.enforce_generator_contingency_headroom =
        options.enforce_generator_contingency_headroom;
    proposal_options.initialize_near_incumbent_dispatch = false;
    proposal_options.repair_candidate = false;
    const auto proposal_start = std::chrono::steady_clock::now();
    const auto proposal = refine_component_economic_commitment(
        data, incumbent_commitment, output.selected, proposal_options);
    output.proposal_wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - proposal_start).count();
    output.proposal_feasible = proposal.solver_feasible;

    struct ToggleCandidate {
        int generator{};
        int target{};
        double score{};
        bool from_proposal{};
    };
    std::vector<ToggleCandidate> pool;
    std::vector<unsigned char> pooled(static_cast<std::size_t>(ng), 0);
    const auto eligible = [&](int generator, int target) {
        const auto& source = data.generators[generator];
        return target == 1
            ? source.suqual != 0
            : source.sdqual != 0;
    };
    std::vector<ToggleCandidate> ranked;
    for (int generator = 0; generator < ng; ++generator) {
        const auto& source = data.generators[generator];
        if (incumbent_commitment[generator] == 1 &&
            source.sdqual != 0) {
            ranked.push_back({
                generator, 0,
                data.delta * source.oncost +
                    data.delta * pwl_value(
                        generator_points[generator],
                        output.selected.state.pg[generator]) -
                    source.sdcost,
                false,
            });
        } else if (incumbent_commitment[generator] == 0 &&
                   source.suqual != 0) {
            const double usable_power = std::max(
                1e-8, generator_upper[generator] -
                    generator_lower[generator]);
            const double incremental_cost =
                data.delta * source.oncost + source.sucost +
                data.delta * (
                    pwl_value(
                        generator_points[generator],
                        generator_upper[generator]) -
                    pwl_value(
                        generator_points[generator],
                        generator_lower[generator]));
            ranked.push_back({
                generator, 1, -incremental_cost / usable_power, false});
        }
    }
    std::sort(
        ranked.begin(), ranked.end(),
        [](const ToggleCandidate& left, const ToggleCandidate& right) {
            if (left.score != right.score) {
                return left.score > right.score;
            }
            return left.generator < right.generator;
        });
    std::vector<double> economic_score_by_generator(
        static_cast<std::size_t>(ng), 0.0);
    for (const auto& candidate : ranked) {
        economic_score_by_generator[candidate.generator] = candidate.score;
    }
    if (proposal.solver_feasible &&
        proposal.candidate_commitment.size() == data.generators.size()) {
        for (int generator = 0; generator < ng; ++generator) {
            const int target = proposal.candidate_commitment[generator];
            if (target == incumbent_commitment[generator] ||
                !eligible(generator, target)) {
                continue;
            }
            ++output.proposal_change_count;
            pool.push_back({
                generator, target,
                economic_score_by_generator[generator], true});
            pooled[generator] = 1;
        }
    }
    const int pool_limit = std::max(
        options.maximum_candidates_per_round,
        options.maximum_candidates_per_round * options.maximum_rounds);
    for (const auto& candidate : ranked) {
        if (static_cast<int>(pool.size()) >= pool_limit) {
            break;
        }
        if (pooled[candidate.generator]) {
            continue;
        }
        pool.push_back(candidate);
        pooled[candidate.generator] = 1;
    }
    output.candidate_pool_size = static_cast<int>(pool.size());

    std::vector<int> current_commitment = output.selected_commitment;
    SolveResult current = output.selected;
    ValidationReport current_validation = output.selected_validation;
    double current_proxy = output.selected_official_proxy;
    double current_penalty = output.selected_penalty_slack;
    double current_headroom_residual = 0.0;

    const auto predicted_reactive_initialization_residual = [&]
        (const ToggleCandidate& toggle,
         const SolveResult& point,
         const std::vector<int>& commitment) {
        const int bus = data.generators[toggle.generator].bus;
        double target_q = 0.0;
        double lower_q = 0.0;
        double upper_q = 0.0;
        for (int generator : data.buses[bus].generators) {
            target_q += point.state.qg[generator];
            const int status = generator == toggle.generator
                ? toggle.target
                : commitment[generator];
            if (status == 0) {
                continue;
            }
            lower_q += data.generators[generator].qmin;
            upper_q += data.generators[generator].qmax;
        }
        return std::max({0.0, lower_q - target_q, target_q - upper_q});
    };

    const auto headroom_residual = [&](const SolveResult& point,
                                       const std::vector<int>& commitment) {
        if (!options.enforce_generator_contingency_headroom) {
            return 0.0;
        }
        std::vector<double> available(static_cast<std::size_t>(ng), 0.0);
        for (int generator = 0; generator < ng; ++generator) {
            if (commitment[generator] == 0) {
                continue;
            }
            const auto& source = data.generators[generator];
            available[generator] = std::max(
                0.0, std::min(
                    source.pmax - point.state.pg[generator],
                    std::max(
                        0.0, data.delta_r_ctg * source.prumaxctg)));
        }
        std::vector<unsigned char> checked(static_cast<std::size_t>(ng), 0);
        double residual = 0.0;
        for (const auto& contingency : data.contingencies) {
            if (contingency.type != ContingencyType::Generator ||
                contingency.component < 0 ||
                contingency.component >= ng ||
                checked[contingency.component]) {
                continue;
            }
            const int outaged = contingency.component;
            checked[outaged] = 1;
            double total = 0.0;
            for (int generator = 0; generator < ng; ++generator) {
                if (generator != outaged) {
                    total += available[generator];
                }
            }
            residual = std::max(
                residual,
                std::max(0.0, point.state.pg[outaged] - total));
        }
        return residual;
    };

    for (int round = 1; round <= options.maximum_rounds; ++round) {
        if (elapsed() >= options.time_limit_seconds - 0.02) {
            output.time_limit_reached = true;
            break;
        }
        bool found = false;
        double best_proxy = current_proxy;
        double best_penalty = current_penalty;
        double best_headroom = current_headroom_residual;
        std::vector<int> best_commitment;
        SolveResult best;
        ValidationReport best_validation;
        int evaluated_this_round = 0;
        struct OrderedToggle {
            const ToggleCandidate* toggle{};
            double reactive_residual{};
        };
        std::vector<OrderedToggle> ordered;
        ordered.reserve(pool.size());
        for (const auto& toggle : pool) {
            if (current_commitment[toggle.generator] == toggle.target) {
                continue;
            }
            ordered.push_back({
                &toggle,
                predicted_reactive_initialization_residual(
                    toggle, current, current_commitment)});
        }
        constexpr double kReactiveEasyTolerance = 1e-9;
        std::sort(
            ordered.begin(), ordered.end(),
            [&](const OrderedToggle& left, const OrderedToggle& right) {
                const bool left_easy =
                    left.reactive_residual <= kReactiveEasyTolerance;
                const bool right_easy =
                    right.reactive_residual <= kReactiveEasyTolerance;
                if (left_easy != right_easy) {
                    return left_easy;
                }
                if (left.toggle->from_proposal !=
                    right.toggle->from_proposal) {
                    return left.toggle->from_proposal;
                }
                if (left_easy && left.toggle->score !=
                        right.toggle->score) {
                    return left.toggle->score > right.toggle->score;
                }
                if (!left_easy && left.reactive_residual !=
                        right.reactive_residual) {
                    return left.reactive_residual <
                        right.reactive_residual;
                }
                if (left.toggle->score != right.toggle->score) {
                    return left.toggle->score > right.toggle->score;
                }
                return left.toggle->generator <
                    right.toggle->generator;
            });
        ++output.candidate_order_refreshes;
        for (const auto& ordered_toggle : ordered) {
            const auto& toggle = *ordered_toggle.toggle;
            if (evaluated_this_round >=
                    options.maximum_candidates_per_round ||
                elapsed() >= options.time_limit_seconds - 0.02) {
                break;
            }
            if (current_commitment[toggle.generator] == toggle.target) {
                continue;
            }
            ++evaluated_this_round;
            ++output.candidates_attempted;
            const auto trial_start = std::chrono::steady_clock::now();
            nlohmann::json trial = {
                {"round", round},
                {"generator", toggle.generator},
                {"source_key",
                 data.generators[toggle.generator].source_key},
                {"from", current_commitment[toggle.generator]},
                {"to", toggle.target},
                {"from_proposal", toggle.from_proposal},
                {"rank_score", toggle.score},
                {"predicted_q_initialization_residual",
                 ordered_toggle.reactive_residual},
            };
            std::vector<int> candidate_commitment = current_commitment;
            candidate_commitment[toggle.generator] = toggle.target;
            AcState candidate_state = current.state;
            candidate_state.commitment.assign(
                static_cast<std::size_t>(ng), 0.0);
            candidate_state.startup.assign(
                static_cast<std::size_t>(ng), 0.0);
            candidate_state.shutdown.assign(
                static_cast<std::size_t>(ng), 0.0);
            for (int generator = 0; generator < ng; ++generator) {
                const int status = candidate_commitment[generator];
                candidate_state.commitment[generator] = status;
                candidate_state.startup[generator] = std::max(
                    0, status - data.generators[generator].status_prev);
                candidate_state.shutdown[generator] = std::max(
                    0, data.generators[generator].status_prev - status);
                if (status == 0) {
                    candidate_state.pg[generator] = 0.0;
                    candidate_state.qg[generator] = 0.0;
                } else if (current_commitment[generator] == 0) {
                    candidate_state.pg[generator] =
                        generator_lower[generator];
                    candidate_state.qg[generator] = std::clamp(
                        data.generators[generator].qg_start,
                        data.generators[generator].qmin,
                        data.generators[generator].qmax);
                } else {
                    candidate_state.pg[generator] = std::clamp(
                        candidate_state.pg[generator],
                        generator_lower[generator],
                        generator_upper[generator]);
                    candidate_state.qg[generator] = std::clamp(
                        candidate_state.qg[generator],
                        data.generators[generator].qmin,
                        data.generators[generator].qmax);
                }
            }
            const int component = component_of_bus[
                data.generators[toggle.generator].bus];
            double target_generation = 0.0;
            double candidate_generation = 0.0;
            std::vector<int> adjustable;
            for (int generator : component_generators[component]) {
                target_generation += current.state.pg[generator];
                candidate_generation += candidate_state.pg[generator];
                if (candidate_commitment[generator] == 1) {
                    adjustable.push_back(generator);
                }
            }
            double generation_difference =
                target_generation - candidate_generation;
            std::sort(
                adjustable.begin(), adjustable.end(),
                [&](int left, int right) {
                    const double left_cost = marginal_cost(
                        left, candidate_state.pg[left]);
                    const double right_cost = marginal_cost(
                        right, candidate_state.pg[right]);
                    if (left_cost != right_cost) {
                        return generation_difference > 0.0
                            ? left_cost < right_cost
                            : left_cost > right_cost;
                    }
                    return left < right;
                });
            for (int generator : adjustable) {
                if (std::abs(generation_difference) <= 1e-9) {
                    break;
                }
                const double room = generation_difference > 0.0
                    ? generator_upper[generator] -
                        candidate_state.pg[generator]
                    : candidate_state.pg[generator] -
                        generator_lower[generator];
                const double movement = std::min(
                    std::abs(generation_difference), std::max(0.0, room));
                candidate_state.pg[generator] +=
                    std::copysign(movement, generation_difference);
                generation_difference -=
                    std::copysign(movement, generation_difference);
            }
            bool precheck_passed =
                std::abs(generation_difference) <= 1e-7;
            double maximum_q_initialization_residual = 0.0;
            int q_initialization_changed_bus_count = 0;
            int q_initialization_empty_bus_count = 0;
            for (int bus = 0; bus < nb; ++bus) {
                double target_q = 0.0;
                double candidate_q = 0.0;
                double lower_q = 0.0;
                double upper_q = 0.0;
                std::vector<int> online;
                for (int generator : data.buses[bus].generators) {
                    target_q += current.state.qg[generator];
                    if (candidate_commitment[generator] == 0) {
                        continue;
                    }
                    online.push_back(generator);
                    candidate_q += candidate_state.qg[generator];
                    lower_q += data.generators[generator].qmin;
                    upper_q += data.generators[generator].qmax;
                }
                if (online.empty()) {
                    maximum_q_initialization_residual = std::max(
                        maximum_q_initialization_residual,
                        std::abs(target_q));
                    if (std::abs(target_q) > 1e-9) {
                        ++q_initialization_changed_bus_count;
                        ++q_initialization_empty_bus_count;
                    }
                    continue;
                }
                const double attainable_target = std::clamp(
                    target_q, lower_q, upper_q);
                double difference = attainable_target - candidate_q;
                for (int pass = 0; pass < 4 &&
                     std::abs(difference) > 1e-10; ++pass) {
                    double room = 0.0;
                    for (int generator : online) {
                        room += difference > 0.0
                            ? data.generators[generator].qmax -
                                candidate_state.qg[generator]
                            : candidate_state.qg[generator] -
                                data.generators[generator].qmin;
                    }
                    if (room <= 1e-12) {
                        break;
                    }
                    for (int generator : online) {
                        const double individual = difference > 0.0
                            ? data.generators[generator].qmax -
                                candidate_state.qg[generator]
                            : candidate_state.qg[generator] -
                                data.generators[generator].qmin;
                        candidate_state.qg[generator] +=
                            difference * individual / room;
                        candidate_state.qg[generator] = std::clamp(
                            candidate_state.qg[generator],
                            data.generators[generator].qmin,
                            data.generators[generator].qmax);
                    }
                    candidate_q = 0.0;
                    for (int generator : online) {
                        candidate_q += candidate_state.qg[generator];
                    }
                    difference = attainable_target - candidate_q;
                }
                const double residual = std::abs(target_q - candidate_q);
                maximum_q_initialization_residual = std::max(
                    maximum_q_initialization_residual, residual);
                if (residual > 1e-9) {
                    ++q_initialization_changed_bus_count;
                }
            }
            trial["maximum_q_initialization_residual"] =
                maximum_q_initialization_residual;
            trial["q_initialization_changed_bus_count"] =
                q_initialization_changed_bus_count;
            trial["q_initialization_empty_bus_count"] =
                q_initialization_empty_bus_count;
            if (!precheck_passed) {
                ++output.candidates_precheck_rejected;
                trial["status"] = "precheck_rejected";
                trial["wall_seconds"] = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - trial_start).count();
                output.trials.push_back(std::move(trial));
                continue;
            }

            rebuild_base_state_derived_fields(
                data, candidate_commitment, candidate_state);
            FastPowerFlowOptions fast_options;
            fast_options.minimize_active_balance_slack = true;
            fast_options.minimize_reactive_balance_slack = true;
            fast_options.skip_balance_cleanup_prepasses = false;
            fast_options.max_newton_iterations = 40;
            fast_options.max_active_redispatch_passes = 16;
            fast_options.max_reactive_limit_passes = 10;
            FastContingencyPowerFlow repair(
                data, candidate_state, candidate_commitment, fast_options);
            const auto repair_start = std::chrono::steady_clock::now();
            auto repaired = repair.solve_base();
            const double repair_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - repair_start).count();
            output.candidate_repair_wall_seconds += repair_seconds;
            ++output.candidates_repaired;
            repaired.solve.status = 0;
            repaired.solve.objective = rebuild_base_state_derived_fields(
                data, candidate_commitment, repaired.solve.state);
            repaired.validation = validate_state(
                data, ModelMode::BaseSoft, repaired.solve.state,
                candidate_commitment);
            const double penalty =
                total_penalty_slack(repaired.solve.state);
            const double headroom =
                headroom_residual(repaired.solve, candidate_commitment);
            const double proxy = repaired.solve.objective -
                base_commitment_transition_cost(
                    data, candidate_commitment);
            const bool verified = validated_candidate_is_feasible(
                    repaired.solve, repaired.validation,
                    options.validation_tolerance) &&
                penalty <= current_penalty + 1e-8 &&
                headroom <= 1e-6;
            if (verified) {
                ++output.candidates_verified;
            }
            const bool improves = verified &&
                proxy > best_proxy + options.objective_tolerance;
            trial["status"] = verified ? "verified" : "rejected";
            trial["repair_feasible"] = repaired.feasible;
            trial["repair_converged"] = repaired.converged;
            trial["repair_failure_reason"] = repaired.failure_reason;
            trial["objective"] = repaired.solve.objective;
            trial["official_proxy"] = proxy;
            trial["penalty_slack"] = penalty;
            trial["headroom_residual"] = headroom;
            trial["max_residual"] = repaired.validation.max_residual;
            trial["improves_round_incumbent"] = improves;
            trial["wall_seconds"] = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - trial_start).count();
            output.trials.push_back(std::move(trial));
            if (improves) {
                found = true;
                best_proxy = proxy;
                best_penalty = penalty;
                best_headroom = headroom;
                best_commitment = std::move(candidate_commitment);
                best = std::move(repaired.solve);
                best_validation = repaired.validation;
                if (options.accept_first_improving_toggle) {
                    ++output.first_improvement_selections;
                    break;
                }
            }
        }
        if (!found) {
            output.status = "no_improving_verified_toggle";
            break;
        }
        current_commitment = std::move(best_commitment);
        current = std::move(best);
        current_validation = best_validation;
        current_proxy = best_proxy;
        current_penalty = best_penalty;
        current_headroom_residual = best_headroom;
        ++output.rounds_completed;
        ++output.accepted_moves;
    }

    output.selected_commitment = std::move(current_commitment);
    output.selected = std::move(current);
    output.selected_validation = current_validation;
    output.selected_objective = output.selected.objective;
    output.selected_official_proxy = current_proxy;
    output.selected_penalty_slack = current_penalty;
    output.selected_headroom_residual = current_headroom_residual;
    output.improved = output.selected_official_proxy >
        output.incumbent_official_proxy + options.objective_tolerance;
    if (output.status.empty()) {
        output.status = output.time_limit_reached
            ? "time_limit_reached"
            : (output.improved ? "completed" : "no_improvement");
    }
    output.wall_seconds = elapsed();
    output.time_limit_reached = output.time_limit_reached ||
        output.wall_seconds >= options.time_limit_seconds;
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
