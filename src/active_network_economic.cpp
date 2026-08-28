#include <highs/Highs.h>

#include "gravityx/algorithm.hpp"
#include "gravityx/fast_power_flow.hpp"
#include "gravityx/state_io.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace gravityx {
namespace {

struct SparseRow {
    double lower{-kHighsInf};
    double upper{kHighsInf};
    std::vector<std::pair<HighsInt, double>> entries;
};

struct FlowCoefficients {
    double from_g_self{};
    double from_b_self{};
    double to_g_self{};
    double to_b_self{};
    double from_cross_cos{};
    double from_cross_sin{};
    double to_cross_cos{};
    double to_cross_sin{};
};

FlowCoefficients flow_coefficients(const Branch& branch) {
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
            "zero tap in active-network economic model: " +
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

double pwl_slope(const std::vector<PwlPoint>& points, double power) {
    if (points.size() < 2) {
        throw std::runtime_error(
            "active-network economic slope requires two PWL points");
    }
    const auto segment_slope = [&](int segment) {
        const double width =
            points[segment + 1].mw - points[segment].mw;
        return std::abs(width) <= 1e-14
            ? 0.0
            : (points[segment + 1].cost - points[segment].cost) / width;
    };
    for (int point = 1;
         point + 1 < static_cast<int>(points.size()); ++point) {
        if (std::abs(power - points[point].mw) <= 1e-10) {
            return 0.5 *
                (segment_slope(point - 1) + segment_slope(point));
        }
    }
    for (int segment = 0;
         segment + 1 < static_cast<int>(points.size()); ++segment) {
        if (power <= points[segment + 1].mw + 1e-12) {
            return segment_slope(segment);
        }
    }
    return segment_slope(static_cast<int>(points.size()) - 2);
}

void append(SparseRow& row, int column, double coefficient) {
    if (!std::isfinite(coefficient) || std::abs(coefficient) > 1e-14) {
        row.entries.emplace_back(column, coefficient);
    }
}

void normalize_row(SparseRow& row) {
    std::sort(
        row.entries.begin(), row.entries.end(),
        [](const auto& left, const auto& right) {
            return left.first < right.first;
        });
    std::vector<std::pair<HighsInt, double>> combined;
    combined.reserve(row.entries.size());
    for (const auto& entry : row.entries) {
        if (!combined.empty() && combined.back().first == entry.first) {
            combined.back().second += entry.second;
        } else {
            combined.push_back(entry);
        }
    }
    combined.erase(
        std::remove_if(
            combined.begin(), combined.end(),
            [](const auto& entry) {
                return std::abs(entry.second) <= 1e-12;
            }),
        combined.end());
    double scale = 1.0;
    for (const auto& [column, coefficient] : combined) {
        static_cast<void>(column);
        scale = std::max(scale, std::abs(coefficient));
    }
    if (scale > 1.0) {
        for (auto& [column, coefficient] : combined) {
            static_cast<void>(column);
            coefficient /= scale;
        }
        if (std::isfinite(row.lower)) {
            row.lower /= scale;
        }
        if (std::isfinite(row.upper)) {
            row.upper /= scale;
        }
    }
    row.entries = std::move(combined);
}

}  // namespace

nlohmann::json ActiveNetworkEconomicDispatchResult::to_json(
    bool include_state) const {
    return {
        {"incumbent_verified", incumbent_verified},
        {"attempted", attempted},
        {"solver_feasible", solver_feasible},
        {"all_solver_rounds_optimal", all_solver_rounds_optimal},
        {"improved", improved},
        {"time_limit_reached", time_limit_reached},
        {"wall_seconds", wall_seconds},
        {"solver_wall_seconds", solver_wall_seconds},
        {"incumbent_objective", incumbent_objective},
        {"selected_objective", selected_objective},
        {"selected_fraction", selected_fraction},
        {"maximum_selected_angle_change", maximum_selected_angle_change},
        {"component_count", component_count},
        {"rounds_completed", rounds_completed},
        {"thermal_row_count", thermal_row_count},
        {"row_count", row_count},
        {"column_count", column_count},
        {"nonzero_count", nonzero_count},
        {"simplex_iterations", simplex_iterations},
        {"ipm_iterations", ipm_iterations},
        {"status", status},
        {"rounds", rounds},
        {"trials", trials},
        {"selected", solve_result_to_json(selected, include_state)},
        {"selected_validation", selected_validation.to_json()},
    };
}

ActiveNetworkEconomicDispatchResult
refine_fixed_commitment_active_network_economic(
    const CaseData& data,
    const std::vector<int>& commitment,
    const SolveResult& incumbent,
    const ActiveNetworkEconomicDispatchOptions& options) {
    const auto wall_start = std::chrono::steady_clock::now();
    ActiveNetworkEconomicDispatchResult output;
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
        !std::isfinite(options.angle_trust_radius) ||
        options.angle_trust_radius <= 0.0 ||
        !std::isfinite(options.thermal_row_utilization_threshold) ||
        options.thermal_row_utilization_threshold < 0.0 ||
        options.thermal_row_utilization_threshold > 1.0 ||
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
    // All physical controls are represented as nonnegative up/down movements
    // from the verified reference.  The zero vector is therefore an exact
    // feasible primal point and supports a deterministic sparse crash basis.
    const int pg_up_offset = 0;
    const int pg_down_offset = pg_up_offset + ng;
    const int demand_up_offset = pg_down_offset + ng;
    const int demand_down_offset = demand_up_offset + nd;
    const int angle_up_offset = demand_down_offset + nd;
    const int angle_down_offset = angle_up_offset + nb;
    const int column_count = angle_down_offset + nb;
    output.column_count = column_count;

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
        if (!branch.present || branch.status == 0) {
            continue;
        }
        const int left = find_root(branch.from);
        const int right = find_root(branch.to);
        if (left != right) {
            parent[std::max(left, right)] = std::min(left, right);
        }
    }
    std::vector<int> component_reference;
    std::vector<int> root_to_component(static_cast<std::size_t>(nb), -1);
    std::vector<int> component_of_bus(static_cast<std::size_t>(nb), -1);
    for (int bus = 0; bus < nb; ++bus) {
        const int root = find_root(bus);
        if (root_to_component[root] < 0) {
            root_to_component[root] = output.component_count++;
            component_reference.push_back(bus);
        }
        component_of_bus[bus] = root_to_component[root];
    }

    for (int round = 1; round <= options.maximum_rounds; ++round) {
        const double elapsed_before_round = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        const double remaining =
            options.time_limit_seconds - elapsed_before_round;
        if (remaining <= 0.05) {
            output.time_limit_reached = true;
            break;
        }
        const SolveResult reference = output.selected;

        std::vector<double> lower(
            static_cast<std::size_t>(column_count), 0.0);
        std::vector<double> upper(
            static_cast<std::size_t>(column_count), 0.0);
        std::vector<double> cost(
            static_cast<std::size_t>(column_count), 0.0);
        for (int generator = 0; generator < ng; ++generator) {
            const auto& source = data.generators[generator];
            if (commitment[generator] == 0) {
                continue;
            }
            const double previous = source.status_prev == 0
                ? source.pmin : source.pg_prev;
            const double physical_lower = std::max(
                source.pmin, previous - data.delta_r * source.prdmax);
            const double physical_upper = std::min(
                source.pmax, previous + data.delta_r * source.prumax);
            if (reference.state.pg[generator] < physical_lower - 1e-8 ||
                reference.state.pg[generator] > physical_upper + 1e-8) {
                output.status = "reference_generator_outside_source_bounds";
                break;
            }
            upper[pg_up_offset + generator] = std::max(
                0.0, physical_upper - reference.state.pg[generator]);
            upper[pg_down_offset + generator] = std::max(
                0.0, reference.state.pg[generator] - physical_lower);
            const auto points = active_pwl_points(
                source.cost, source.ncost,
                physical_lower, physical_upper);
            const double slope = data.delta *
                pwl_slope(points, reference.state.pg[generator]);
            cost[pg_up_offset + generator] = slope;
            cost[pg_down_offset + generator] = -slope;
        }
        if (output.status ==
            "reference_generator_outside_source_bounds") {
            break;
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
            if (reference.state.demand_factor[load] < factor_lower - 1e-8 ||
                reference.state.demand_factor[load] > factor_upper + 1e-8) {
                output.status = "reference_load_outside_source_bounds";
                break;
            }
            upper[demand_up_offset + load] = std::max(
                0.0,
                factor_upper - reference.state.demand_factor[load]);
            upper[demand_down_offset + load] = std::max(
                0.0,
                reference.state.demand_factor[load] - factor_lower);
            const auto points = active_pwl_points(
                source.cost, source.ncost,
                source.pd_min, source.pd_max);
            const double slope = -data.delta *
                source.pd_nominal * pwl_slope(
                    points, source.pd_nominal *
                        reference.state.demand_factor[load]);
            cost[demand_up_offset + load] = slope;
            cost[demand_down_offset + load] = -slope;
        }
        if (output.status == "reference_load_outside_source_bounds") {
            break;
        }
        for (int bus = 0; bus < nb; ++bus) {
            upper[angle_up_offset + bus] = options.angle_trust_radius;
            upper[angle_down_offset + bus] = options.angle_trust_radius;
        }
        for (int bus : component_reference) {
            upper[angle_up_offset + bus] = 0.0;
            upper[angle_down_offset + bus] = 0.0;
        }
        double objective_scale = 1.0;
        for (double coefficient : cost) {
            objective_scale = std::max(
                objective_scale, std::abs(coefficient));
        }
        for (double& coefficient : cost) {
            coefficient /= objective_scale;
        }

        std::vector<double> dpf(data.branches.size(), 0.0);
        std::vector<double> dqf(data.branches.size(), 0.0);
        std::vector<double> dpt(data.branches.size(), 0.0);
        std::vector<double> dqt(data.branches.size(), 0.0);
        for (int index = 0;
             index < static_cast<int>(data.branches.size()); ++index) {
            const auto& branch = data.branches[index];
            if (!branch.present || branch.status == 0) {
                continue;
            }
            const auto coefficient = flow_coefficients(branch);
            const double angle =
                reference.state.va[branch.from] -
                reference.state.va[branch.to];
            const double voltage_product =
                reference.state.vm[branch.from] *
                reference.state.vm[branch.to];
            const double sine = std::sin(angle);
            const double cosine = std::cos(angle);
            dpf[index] = voltage_product *
                (-coefficient.from_cross_cos * sine +
                 coefficient.from_cross_sin * cosine);
            dqf[index] = voltage_product *
                (coefficient.from_cross_sin * sine +
                 coefficient.from_cross_cos * cosine);
            dpt[index] = voltage_product *
                (-coefficient.to_cross_cos * sine -
                 coefficient.to_cross_sin * cosine);
            dqt[index] = voltage_product *
                (coefficient.to_cross_sin * sine -
                 coefficient.to_cross_cos * cosine);
        }

        std::vector<SparseRow> rows;
        rows.reserve(static_cast<std::size_t>(
            nb + data.branches.size()));
        for (int bus = 0; bus < nb; ++bus) {
            SparseRow balance;
            for (int generator : data.buses[bus].generators) {
                append(balance, pg_up_offset + generator, 1.0);
                append(balance, pg_down_offset + generator, -1.0);
            }
            for (int load : data.buses[bus].loads) {
                append(
                    balance, demand_up_offset + load,
                    -data.loads[load].pd_nominal);
                append(
                    balance, demand_down_offset + load,
                    data.loads[load].pd_nominal);
            }
            for (int branch_index : data.buses[bus].branches_from) {
                if (!data.branches[branch_index].present ||
                    data.branches[branch_index].status == 0) {
                    continue;
                }
                const auto& branch = data.branches[branch_index];
                append(
                    balance, angle_up_offset + branch.from,
                    -dpf[branch_index]);
                append(
                    balance, angle_down_offset + branch.from,
                    dpf[branch_index]);
                append(
                    balance, angle_up_offset + branch.to,
                    dpf[branch_index]);
                append(
                    balance, angle_down_offset + branch.to,
                    -dpf[branch_index]);
            }
            for (int branch_index : data.buses[bus].branches_to) {
                if (!data.branches[branch_index].present ||
                    data.branches[branch_index].status == 0) {
                    continue;
                }
                const auto& branch = data.branches[branch_index];
                append(
                    balance, angle_up_offset + branch.from,
                    -dpt[branch_index]);
                append(
                    balance, angle_down_offset + branch.from,
                    dpt[branch_index]);
                append(
                    balance, angle_up_offset + branch.to,
                    dpt[branch_index]);
                append(
                    balance, angle_down_offset + branch.to,
                    -dpt[branch_index]);
            }
            balance.lower = 0.0;
            balance.upper = 0.0;
            normalize_row(balance);
            rows.push_back(std::move(balance));
        }

        int thermal_rows = 0;
        for (int index = 0;
             index < static_cast<int>(data.branches.size()); ++index) {
            const auto& branch = data.branches[index];
            if (!branch.present || branch.status == 0) {
                continue;
            }
            const double start_delta =
                data.buses[branch.from].va_start -
                data.buses[branch.to].va_start;
            const double reference_delta =
                reference.state.va[branch.from] -
                reference.state.va[branch.to];
            // Each endpoint can move by at most the angle trust radius.  If
            // that whole box remains inside the source angle interval, the
            // explicit branch row is exactly redundant and is omitted.
            const bool trust_box_can_reach_angle_limit =
                reference_delta - 2.0 * options.angle_trust_radius <
                    branch.angmin ||
                reference_delta + 2.0 * options.angle_trust_radius >
                    branch.angmax;
            if (start_delta >= branch.angmin &&
                start_delta <= branch.angmax &&
                trust_box_can_reach_angle_limit) {
                SparseRow angle_limit;
                angle_limit.lower = branch.angmin -
                    reference_delta;
                angle_limit.upper = branch.angmax -
                    reference_delta;
                append(angle_limit, angle_up_offset + branch.from, 1.0);
                append(angle_limit, angle_down_offset + branch.from, -1.0);
                append(angle_limit, angle_up_offset + branch.to, -1.0);
                append(angle_limit, angle_down_offset + branch.to, 1.0);
                rows.push_back(std::move(angle_limit));
            }
            if (branch.rate_a <= 1e-12) {
                continue;
            }
            const double slack = reference.state.sm_slack.size() ==
                    data.branches.size()
                ? reference.state.sm_slack[index]
                : 0.0;
            const double from_limit = branch.rate_a *
                (branch.transformer
                    ? 1.0 + slack
                    : reference.state.vm[branch.from] + slack);
            const double to_limit = branch.rate_a *
                (branch.transformer
                    ? 1.0 + slack
                    : reference.state.vm[branch.to] + slack);
            const double from_magnitude = std::hypot(
                reference.state.pf[index], reference.state.qf[index]);
            const double to_magnitude = std::hypot(
                reference.state.pt[index], reference.state.qt[index]);
            const auto add_thermal_row = [&](double magnitude,
                                             double limit,
                                             double p,
                                             double q,
                                             double dp,
                                             double dq) {
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
                    thermal, angle_up_offset + branch.from, derivative);
                append(
                    thermal, angle_down_offset + branch.from, -derivative);
                append(
                    thermal, angle_up_offset + branch.to, -derivative);
                append(
                    thermal, angle_down_offset + branch.to, derivative);
                normalize_row(thermal);
                rows.push_back(std::move(thermal));
                ++thermal_rows;
            };
            add_thermal_row(
                from_magnitude, from_limit,
                reference.state.pf[index], reference.state.qf[index],
                dpf[index], dqf[index]);
            add_thermal_row(
                to_magnitude, to_limit,
                reference.state.pt[index], reference.state.qt[index],
                dpt[index], dqt[index]);
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
        bool finite_model = true;
        for (const auto& row : rows) {
            finite_model = finite_model && !std::isnan(row.lower) &&
                !std::isnan(row.upper) && row.lower <= row.upper;
            row_lower.push_back(row.lower);
            row_upper.push_back(row.upper);
            for (const auto& [column, coefficient] : row.entries) {
                finite_model = finite_model && column >= 0 &&
                    column < column_count && std::isfinite(coefficient);
                indices.push_back(column);
                values.push_back(coefficient);
            }
            starts.push_back(static_cast<HighsInt>(indices.size()));
        }
        if (!finite_model) {
            output.status = "nonfinite_linearized_model";
            break;
        }
        output.row_count = static_cast<int>(rows.size());
        output.nonzero_count = static_cast<int>(indices.size());
        output.thermal_row_count = thermal_rows;

        Highs highs;
        const char* highs_log = std::getenv("GRAVITYX_HIGHS_LOG");
        highs.setOptionValue(
            "output_flag",
            highs_log != nullptr && std::string(highs_log) != "0");
        highs.setOptionValue("threads", 1);
        const char* solver_override =
            std::getenv("GRAVITYX_ACTIVE_ECONOMIC_SOLVER");
        const std::string solver = solver_override != nullptr
            ? std::string(solver_override) : "simplex";
        highs.setOptionValue("presolve", solver == "simplex" ? "off" : "on");
        highs.setOptionValue("solver", solver);
        if (solver == "simplex") {
            highs.setOptionValue("simplex_strategy", 4);
            highs.setOptionValue(
                "primal_simplex_bound_perturbation_multiplier", 0.0);
        }
        highs.setOptionValue("primal_feasibility_tolerance", 1e-8);
        highs.setOptionValue("dual_feasibility_tolerance", 1e-8);
        highs.setOptionValue(
            "time_limit", std::max(
                0.05, std::min(0.5 * remaining, remaining - 0.02)));
        const bool model_loaded =
            highs.addVars(column_count, lower.data(), upper.data()) ==
                HighsStatus::kOk &&
            highs.changeColsCost(
                0, column_count - 1, cost.data()) == HighsStatus::kOk &&
            highs.addRows(
                static_cast<HighsInt>(rows.size()), row_lower.data(),
                row_upper.data(), static_cast<HighsInt>(indices.size()),
                starts.data(), indices.data(), values.data()) ==
                HighsStatus::kOk;
        if (!model_loaded) {
            output.status = "model_construction_failed";
            break;
        }
        std::vector<double> primal_start(
            static_cast<std::size_t>(column_count), 0.0);
        std::vector<HighsInt> start_indices(
            static_cast<std::size_t>(column_count));
        std::iota(start_indices.begin(), start_indices.end(), HighsInt{0});
        const HighsStatus start_status = highs.setSolution(
            column_count, start_indices.data(), primal_start.data());
        HighsStatus basis_status = HighsStatus::kOk;
        bool basis_attempted = false;
        if (solver == "simplex") {
            basis_attempted = true;
            HighsBasis basis;
            basis.alien = true;
            basis.useful = true;
            basis.col_status.assign(
                static_cast<std::size_t>(column_count),
                HighsBasisStatus::kLower);
            basis.row_status.assign(
                rows.size(), HighsBasisStatus::kBasic);
            for (int bus = 0; bus < nb; ++bus) {
                basis.row_status[bus] = HighsBasisStatus::kLower;
            }
            std::vector<unsigned char> is_reference(
                static_cast<std::size_t>(nb), 0);
            for (int bus : component_reference) {
                is_reference[bus] = 1;
            }
            for (int bus = 0; bus < nb; ++bus) {
                if (is_reference[bus] == 0) {
                    basis.col_status[angle_up_offset + bus] =
                        HighsBasisStatus::kBasic;
                }
            }
            std::vector<unsigned char> injection_chosen(
                static_cast<std::size_t>(output.component_count), 0);
            const auto choose_injection = [&](int component, int column) {
                if (injection_chosen[component] == 0 &&
                    upper[column] > 1e-12) {
                    basis.col_status[column] = HighsBasisStatus::kBasic;
                    injection_chosen[component] = 1;
                }
            };
            for (int generator = 0; generator < ng; ++generator) {
                const int component = component_of_bus[
                    data.generators[generator].bus];
                choose_injection(
                    component, pg_up_offset + generator);
                choose_injection(
                    component, pg_down_offset + generator);
            }
            for (int load = 0; load < nd; ++load) {
                const int component = component_of_bus[data.loads[load].bus];
                choose_injection(component, demand_up_offset + load);
                choose_injection(component, demand_down_offset + load);
            }
            const bool complete_basis = std::all_of(
                injection_chosen.begin(), injection_chosen.end(),
                [](unsigned char value) { return value != 0; });
            basis_status = complete_basis
                ? highs.setBasis(basis, "active_network_tree_crash")
                : HighsStatus::kError;
        }

        const auto solver_start = std::chrono::steady_clock::now();
        const HighsStatus run_status = highs.run();
        const double solver_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - solver_start).count();
        output.solver_wall_seconds += solver_seconds;
        const auto model_status = highs.getModelStatus();
        const auto& solution = highs.getSolution();
        const auto& info = highs.getInfo();
        const bool round_optimal =
            model_status == HighsModelStatus::kOptimal;
        output.all_solver_rounds_optimal =
            output.all_solver_rounds_optimal && round_optimal;
        output.time_limit_reached = output.time_limit_reached ||
            model_status == HighsModelStatus::kTimeLimit;
        constexpr double kRequiredPrimalFeasibility = 1e-8;
        const bool round_feasible = info.valid && solution.value_valid &&
            solution.col_value.size() ==
                static_cast<std::size_t>(column_count) &&
            info.primal_solution_status == kSolutionStatusFeasible &&
            std::isfinite(info.max_primal_infeasibility) &&
            info.max_primal_infeasibility <=
                kRequiredPrimalFeasibility;
        output.solver_feasible = output.solver_feasible || round_feasible;
        output.simplex_iterations +=
            static_cast<int>(info.simplex_iteration_count);
        output.ipm_iterations += static_cast<int>(info.ipm_iteration_count);
        output.status = highs.modelStatusToString(model_status);
        output.rounds.push_back({
            {"round", round},
            {"row_count", rows.size()},
            {"column_count", column_count},
            {"nonzero_count", indices.size()},
            {"thermal_row_count", thermal_rows},
            {"objective_scale", objective_scale},
            {"solver", solver},
            {"primal_start_status", static_cast<int>(start_status)},
            {"basis_attempted", basis_attempted},
            {"basis_status", static_cast<int>(basis_status)},
            {"run_status", static_cast<int>(run_status)},
            {"model_status", static_cast<int>(model_status)},
            {"status", output.status},
            {"optimal", round_optimal},
            {"solution_value_valid", solution.value_valid},
            {"info_valid", info.valid},
            {"primal_solution_status", info.primal_solution_status},
            {"max_primal_infeasibility", info.max_primal_infeasibility},
            {"simplex_iterations", info.simplex_iteration_count},
            {"ipm_iterations", info.ipm_iteration_count},
            {"solver_wall_seconds", solver_seconds},
        });
        if (!round_feasible) {
            break;
        }
        ++output.rounds_completed;

        AcState target = reference.state;
        for (int generator = 0; generator < ng; ++generator) {
            target.pg[generator] +=
                solution.col_value[pg_up_offset + generator] -
                solution.col_value[pg_down_offset + generator];
        }
        for (int load = 0; load < nd; ++load) {
            target.demand_factor[load] +=
                solution.col_value[demand_up_offset + load] -
                solution.col_value[demand_down_offset + load];
        }
        for (int bus = 0; bus < nb; ++bus) {
            target.va[bus] +=
                solution.col_value[angle_up_offset + bus] -
                solution.col_value[angle_down_offset + bus];
        }

        bool round_improved = false;
        for (int trial = 0;
             trial < options.maximum_candidate_trials; ++trial) {
            const double trial_elapsed = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - wall_start).count();
            if (trial_elapsed >= options.time_limit_seconds - 0.01) {
                output.time_limit_reached = true;
                break;
            }
            const double fraction = std::ldexp(1.0, -trial);
            SolveResult proposal;
            proposal.status = 0;
            proposal.state = reference.state;
            for (int generator = 0; generator < ng; ++generator) {
                proposal.state.pg[generator] += fraction *
                    (target.pg[generator] - reference.state.pg[generator]);
            }
            for (int load = 0; load < nd; ++load) {
                proposal.state.demand_factor[load] += fraction *
                    (target.demand_factor[load] -
                     reference.state.demand_factor[load]);
            }
            for (int bus = 0; bus < nb; ++bus) {
                proposal.state.va[bus] += fraction *
                    (target.va[bus] - reference.state.va[bus]);
            }
            proposal.objective = rebuild_base_state_derived_fields(
                data, commitment, proposal.state);
            auto proposal_validation = validate_state(
                data, ModelMode::BaseSoft, proposal.state, commitment);
            const double raw_objective = proposal.objective;
            const auto raw_validation = proposal_validation;
            const double active_slack = std::accumulate(
                proposal.state.p_delta.begin(),
                proposal.state.p_delta.end(), 0.0);
            const double reactive_slack = std::accumulate(
                proposal.state.q_delta.begin(),
                proposal.state.q_delta.end(), 0.0);
            bool repair_attempted = false;
            double repair_seconds = 0.0;
            if (active_slack + reactive_slack > 1e-8 ||
                !validated_candidate_is_feasible(
                    proposal, proposal_validation,
                    options.validation_tolerance)) {
                repair_attempted = true;
                FastPowerFlowOptions fast_options;
                fast_options.minimize_active_balance_slack = true;
                fast_options.minimize_reactive_balance_slack = true;
                fast_options.max_newton_iterations = 40;
                fast_options.max_active_redispatch_passes = 12;
                fast_options.max_reactive_limit_passes = 8;
                FastContingencyPowerFlow repair(
                    data, proposal.state, commitment, fast_options);
                auto repaired = repair.solve_base();
                repair_seconds = repaired.wall_seconds;
                repaired.solve.objective = rebuild_base_state_derived_fields(
                    data, commitment, repaired.solve.state);
                repaired.validation = validate_state(
                    data, ModelMode::BaseSoft,
                    repaired.solve.state, commitment);
                proposal = std::move(repaired.solve);
                proposal_validation = repaired.validation;
            }
            const bool accepted =
                verified_economic_candidate_improves_incumbent(
                    output.selected, output.selected_validation,
                    proposal, proposal_validation,
                    options.validation_tolerance,
                    options.objective_tolerance);
            output.trials.push_back({
                {"round", round},
                {"trial", trial + 1},
                {"fraction", fraction},
                {"raw_objective", raw_objective},
                {"raw_active_balance_slack", active_slack},
                {"raw_reactive_balance_slack", reactive_slack},
                {"raw_validation", raw_validation.to_json()},
                {"repair_attempted", repair_attempted},
                {"repair_seconds", repair_seconds},
                {"candidate_objective", proposal.objective},
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
            }
        }
        if (output.time_limit_reached || !round_improved) {
            break;
        }
    }

    for (int bus = 0; bus < nb; ++bus) {
        output.maximum_selected_angle_change = std::max(
            output.maximum_selected_angle_change,
            std::abs(
                output.selected.state.va[bus] - incumbent.state.va[bus]));
    }
    output.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_start).count();
    output.time_limit_reached = output.time_limit_reached ||
        output.wall_seconds >= options.time_limit_seconds;
    output.selected.wall_seconds = output.wall_seconds;
    output.selected_objective = output.selected.objective;
    return output;
}

}  // namespace gravityx
