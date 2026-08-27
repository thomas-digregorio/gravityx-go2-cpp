#include <highs/Highs.h>

#include "gravityx/active_feasibility_repair.hpp"

#include "gravityx/state_io.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

namespace gravityx {
namespace {

struct SparseRow {
    double lower{};
    double upper{};
    std::vector<std::pair<HighsInt, double>> entries;
};

struct AngleFlowDerivative {
    double pf{};
    double qf{};
    double pt{};
    double qt{};
};

struct FlowValues {
    double pf{};
    double qf{};
    double pt{};
    double qt{};
};

struct VoltageFlowDerivative {
    FlowValues from;
    FlowValues to;
};

void append(SparseRow& row, HighsInt column, double value) {
    if (column >= 0 && std::abs(value) > 1e-14) {
        row.entries.emplace_back(column, value);
    }
}

void normalize(SparseRow& row) {
    std::sort(row.entries.begin(), row.entries.end(),
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
        std::remove_if(combined.begin(), combined.end(),
            [](const auto& entry) {
                return std::abs(entry.second) <= 1e-14;
            }),
        combined.end());
    row.entries = std::move(combined);
}

std::vector<std::vector<int>> connected_components(
    const CaseData& data,
    int outaged_branch) {
    std::vector<std::vector<int>> adjacency(data.buses.size());
    for (int i = 0; i < static_cast<int>(data.branches.size()); ++i) {
        if (i == outaged_branch || data.branches[i].status == 0) {
            continue;
        }
        const auto& branch = data.branches[i];
        adjacency[branch.from].push_back(branch.to);
        adjacency[branch.to].push_back(branch.from);
    }
    std::vector<unsigned char> visited(data.buses.size(), 0);
    std::vector<std::vector<int>> result;
    for (int start = 0; start < static_cast<int>(data.buses.size()); ++start) {
        if (visited[start]) {
            continue;
        }
        std::queue<int> frontier;
        std::vector<int> component;
        frontier.push(start);
        visited[start] = 1;
        while (!frontier.empty()) {
            const int bus = frontier.front();
            frontier.pop();
            component.push_back(bus);
            for (int neighbor : adjacency[bus]) {
                if (!visited[neighbor]) {
                    visited[neighbor] = 1;
                    frontier.push(neighbor);
                }
            }
        }
        result.push_back(std::move(component));
    }
    return result;
}

AngleFlowDerivative angle_flow_derivative(
    const Branch& branch,
    const AcState& state) {
    const double denominator = branch.r * branch.r + branch.x * branch.x;
    if (denominator <= 1e-20 || std::abs(branch.tap) <= 1e-12) {
        return {};
    }
    const double g = branch.r / denominator;
    const double b = -branch.x / denominator;
    const double tm2 = branch.tap * branch.tap;
    const double tr = branch.tap * std::cos(branch.shift);
    const double ti = branch.tap * std::sin(branch.shift);
    const double angle = state.va[branch.from] - state.va[branch.to];
    const double sine = std::sin(angle);
    const double cosine = std::cos(angle);
    const double cross = state.vm[branch.from] * state.vm[branch.to];

    const double pf_cos = ((-g * tr + b * ti) / tm2) * cross;
    const double pf_sin = ((-b * tr - g * ti) / tm2) * cross;
    const double qf_cos = -((-b * tr + g * ti) / tm2) * cross;
    const double qf_sin = ((-g * tr + b * ti) / tm2) * cross;
    const double pt_cos = ((-g * tr - b * ti) / tm2) * cross;
    const double pt_sin = ((-b * tr + g * ti) / tm2) * cross;
    const double qt_cos = -((-b * tr - g * ti) / tm2) * cross;
    const double qt_sin = ((-g * tr - b * ti) / tm2) * cross;
    return {
        -pf_cos * sine + pf_sin * cosine,
        -qf_cos * sine + qf_sin * cosine,
        -pt_cos * sine - pt_sin * cosine,
        -qt_cos * sine - qt_sin * cosine,
    };
}

FlowValues branch_flow_values(
    const Branch& branch,
    double vm_from,
    double vm_to,
    double va_from,
    double va_to) {
    const double denominator = branch.r * branch.r + branch.x * branch.x;
    const double g = denominator > 1e-20 ? branch.r / denominator : 0.0;
    const double b = denominator > 1e-20 ? -branch.x / denominator : 0.0;
    if (std::abs(branch.tap) <= 1e-12) {
        return {};
    }
    const double tm2 = branch.tap * branch.tap;
    const double tr = branch.tap * std::cos(branch.shift);
    const double ti = branch.tap * std::sin(branch.shift);
    const double delta = va_from - va_to;
    const double cross = vm_from * vm_to;
    const double from_g_self = branch.transformer
        ? g / tm2 + branch.g_fr : (g + branch.g_fr) / tm2;
    const double from_b_self = branch.transformer
        ? b / tm2 + branch.b_fr : (b + branch.b_fr) / tm2;
    FlowValues value;
    value.pf = from_g_self * vm_from * vm_from
        + ((-g * tr + b * ti) / tm2) * cross * std::cos(delta)
        + ((-b * tr - g * ti) / tm2) * cross * std::sin(delta);
    value.qf = -from_b_self * vm_from * vm_from
        - ((-b * tr + g * ti) / tm2) * cross * std::cos(delta)
        + ((-g * tr + b * ti) / tm2) * cross * std::sin(delta);
    value.pt = (g + branch.g_to) * vm_to * vm_to
        + ((-g * tr - b * ti) / tm2) * cross * std::cos(-delta)
        + ((-b * tr + g * ti) / tm2) * cross * std::sin(-delta);
    value.qt = -(b + branch.b_to) * vm_to * vm_to
        - ((-b * tr - g * ti) / tm2) * cross * std::cos(-delta)
        + ((-g * tr - b * ti) / tm2) * cross * std::sin(-delta);
    return value;
}

VoltageFlowDerivative voltage_flow_derivative(
    const Branch& branch,
    const AcState& state) {
    constexpr double h = 1e-6;
    const int from = branch.from;
    const int to = branch.to;
    const auto derivative = [h](const FlowValues& plus,
                               const FlowValues& minus) {
        return FlowValues{
            (plus.pf - minus.pf) / (2.0 * h),
            (plus.qf - minus.qf) / (2.0 * h),
            (plus.pt - minus.pt) / (2.0 * h),
            (plus.qt - minus.qt) / (2.0 * h),
        };
    };
    const auto from_plus = branch_flow_values(
        branch, state.vm[from] + h, state.vm[to],
        state.va[from], state.va[to]);
    const auto from_minus = branch_flow_values(
        branch, state.vm[from] - h, state.vm[to],
        state.va[from], state.va[to]);
    const auto to_plus = branch_flow_values(
        branch, state.vm[from], state.vm[to] + h,
        state.va[from], state.va[to]);
    const auto to_minus = branch_flow_values(
        branch, state.vm[from], state.vm[to] - h,
        state.va[from], state.va[to]);
    return {
        derivative(from_plus, from_minus),
        derivative(to_plus, to_minus),
    };
}

double row_violation(const SparseRow& row, const std::vector<double>& value) {
    double activity = 0.0;
    for (const auto& [column, coefficient] : row.entries) {
        activity += coefficient * value[column];
    }
    return std::max({0.0, row.lower - activity, activity - row.upper});
}

}  // namespace

nlohmann::json ActiveFeasibilityRepairResult::to_json(
    bool include_state) const {
    nlohmann::json value = {
        {"success", success},
        {"accepted_feasible_nonoptimal", accepted_feasible_nonoptimal},
        {"wall_seconds", wall_seconds},
        {"time_limit_seconds", time_limit_seconds},
        {"balance_slack_limit", balance_slack_limit},
        {"angle_trust_radius", angle_trust_radius},
        {"voltage_trust_radius", voltage_trust_radius},
        {"include_reactive", include_reactive},
        {"current_security_rows_only", current_security_rows_only},
        {"row_count", row_count},
        {"column_count", column_count},
        {"nonzero_count", nonzero_count},
        {"branch_security_row_count", branch_security_row_count},
        {"simplex_strategy", simplex_strategy},
        {"run_status", run_status},
        {"model_status", model_status},
        {"primal_solution_status", primal_solution_status},
        {"num_primal_infeasibilities", num_primal_infeasibilities},
        {"iterations", iterations},
        {"max_primal_infeasibility", max_primal_infeasibility},
        {"maximum_linearized_violation", maximum_linearized_violation},
        {"maximum_column_violation", maximum_column_violation},
        {"finite_solution_values", finite_solution_values},
        {"maximum_angle_change", maximum_angle_change},
        {"maximum_voltage_change", maximum_voltage_change},
        {"maximum_generation_change", maximum_generation_change},
        {"maximum_reactive_generation_change",
         maximum_reactive_generation_change},
        {"maximum_load_change", maximum_load_change},
        {"objective", objective},
        {"solver", solver},
        {"status", status},
    };
    if (include_state) {
        value["state"] = ac_state_to_json(state);
    }
    return value;
}

ActiveFeasibilityRepairResult solve_linearized_active_feasibility_repair(
    const CaseData& data,
    const AcState& reference,
    const std::vector<int>& commitment,
    const ContingencyContext& contingency,
    double balance_slack_limit,
    double angle_trust_radius,
    double time_limit_seconds,
    double voltage_trust_radius,
    bool include_reactive,
    bool current_security_rows_only) {
    const auto wall_start = std::chrono::steady_clock::now();
    ActiveFeasibilityRepairResult output;
    output.balance_slack_limit = balance_slack_limit;
    output.angle_trust_radius = angle_trust_radius;
    output.voltage_trust_radius = voltage_trust_radius;
    output.include_reactive = include_reactive;
    output.current_security_rows_only = current_security_rows_only;
    output.time_limit_seconds = time_limit_seconds;
    output.state = reference;
    const int nb = static_cast<int>(data.buses.size());
    const int ng = static_cast<int>(data.generators.size());
    const int nd = static_cast<int>(data.loads.size());
    const int nl = static_cast<int>(data.branches.size());
    if (commitment.size() != data.generators.size() ||
        contingency.effective_base_state().pg.size() !=
            data.generators.size() ||
        contingency.effective_base_state().demand_factor.size() !=
            data.loads.size() ||
        reference.vm.size() != data.buses.size() ||
        reference.va.size() != data.buses.size() ||
        reference.pg.size() != data.generators.size() ||
        reference.qg.size() != data.generators.size() ||
        reference.demand_factor.size() != data.loads.size() ||
        reference.pf.size() != data.branches.size() ||
        reference.qf.size() != data.branches.size() ||
        reference.pt.size() != data.branches.size() ||
        reference.qt.size() != data.branches.size()) {
        output.status = "invalid_dimensions";
        output.wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        return output;
    }
    if (!std::isfinite(balance_slack_limit) || balance_slack_limit < 0.0 ||
        balance_slack_limit > 0.5 || !std::isfinite(angle_trust_radius) ||
        angle_trust_radius <= 0.0 ||
        !std::isfinite(voltage_trust_radius) ||
        voltage_trust_radius <= 0.0 ||
        !std::isfinite(time_limit_seconds) ||
        time_limit_seconds <= 0.0) {
        output.status = "invalid_options";
        output.wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        return output;
    }

    const int outaged_generator = contingency.outaged_generator;
    const int outaged_branch = contingency.outaged_branch;
    std::vector<unsigned char> active(static_cast<std::size_t>(ng), 0);
    std::vector<double> p_lower(static_cast<std::size_t>(ng), 0.0);
    std::vector<double> p_upper(static_cast<std::size_t>(ng), 0.0);
    std::vector<double> q_lower(static_cast<std::size_t>(ng), 0.0);
    std::vector<double> q_upper(static_cast<std::size_t>(ng), 0.0);
    for (int i = 0; i < ng; ++i) {
        active[i] = commitment[i] == 1 && i != outaged_generator;
        if (!active[i]) {
            output.state.pg[i] = 0.0;
            output.state.qg[i] = 0.0;
            continue;
        }
        const auto& generator = data.generators[i];
        p_lower[i] = std::max(
            generator.pmin,
            contingency.effective_base_state().pg[i] -
                data.delta_r_ctg * generator.prdmaxctg);
        p_upper[i] = std::min(
            generator.pmax,
            contingency.effective_base_state().pg[i] +
                data.delta_r_ctg * generator.prumaxctg);
        q_lower[i] = generator.qmin;
        q_upper[i] = generator.qmax;
        if (p_lower[i] > p_upper[i] + 1e-12 ||
            q_lower[i] > q_upper[i] + 1e-12) {
            output.status = "empty_generator_interval";
            output.wall_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - wall_start).count();
            return output;
        }
        output.state.pg[i] = std::clamp(
            reference.pg[i], p_lower[i], p_upper[i]);
        output.state.qg[i] = std::clamp(
            reference.qg[i], q_lower[i], q_upper[i]);
    }

    std::vector<double> load_factor_lower(static_cast<std::size_t>(nd), 0.0);
    std::vector<double> load_factor_upper(static_cast<std::size_t>(nd), 0.0);
    std::vector<double> load_factor(static_cast<std::size_t>(nd), 0.0);
    for (int i = 0; i < nd; ++i) {
        const auto& load = data.loads[i];
        const double previous = load.pd_nominal *
            contingency.effective_base_state().demand_factor[i];
        if (std::abs(load.pd_nominal) <= 1e-12) {
            load_factor_lower[i] = load.tmin;
            load_factor_upper[i] = load.tmax;
            load_factor[i] = std::clamp(
                reference.demand_factor[i], load.tmin, load.tmax);
            output.state.demand_factor[i] = load_factor[i];
            continue;
        }
        const double factor_lower = std::max(
            load.tmin,
            (previous - data.delta_r_ctg * load.prdmaxctg) /
                load.pd_nominal);
        const double factor_upper = std::min(
            load.tmax,
            (previous + data.delta_r_ctg * load.prumaxctg) /
                load.pd_nominal);
        if (factor_lower > factor_upper + 1e-12) {
            output.status = "empty_load_interval";
            output.wall_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - wall_start).count();
            return output;
        }
        load_factor_lower[i] = factor_lower;
        load_factor_upper[i] = factor_upper;
        load_factor[i] = std::clamp(
            reference.demand_factor[i], factor_lower, factor_upper);
        output.state.demand_factor[i] = load_factor[i];
    }

    std::vector<int> angle_index(static_cast<std::size_t>(nb), -1);
    int angle_count = 0;
    for (const auto& component : connected_components(data, outaged_branch)) {
        int reference_bus = -1;
        for (int bus : component) {
            bool available_generation = false;
            for (int generator : data.buses[bus].generators) {
                available_generation = available_generation || active[generator];
            }
            if (available_generation && data.buses[bus].type == 3) {
                reference_bus = bus;
                break;
            }
        }
        if (reference_bus < 0) {
            for (int bus : component) {
                for (int generator : data.buses[bus].generators) {
                    if (active[generator]) {
                        reference_bus = bus;
                        break;
                    }
                }
                if (reference_bus >= 0) {
                    break;
                }
            }
        }
        if (reference_bus < 0) {
            output.status = "component_without_available_generator";
            output.wall_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - wall_start).count();
            return output;
        }
        for (int bus : component) {
            if (bus != reference_bus) {
                angle_index[bus] = angle_count++;
            }
        }
    }

    std::vector<AngleFlowDerivative> derivative(
        static_cast<std::size_t>(nl));
    std::vector<VoltageFlowDerivative> voltage_derivative(
        static_cast<std::size_t>(nl));
    std::vector<double> angle_scale(
        static_cast<std::size_t>(angle_count), 1.0);
    std::vector<double> voltage_scale(static_cast<std::size_t>(nb), 1.0);
    for (int i = 0; i < nl; ++i) {
        if (i == outaged_branch || data.branches[i].status == 0) {
            continue;
        }
        derivative[i] = angle_flow_derivative(data.branches[i], output.state);
        voltage_derivative[i] = voltage_flow_derivative(
            data.branches[i], output.state);
        const double maximum = std::max({
            std::abs(derivative[i].pf), std::abs(derivative[i].qf),
            std::abs(derivative[i].pt), std::abs(derivative[i].qt),
        });
        for (int bus : {data.branches[i].from, data.branches[i].to}) {
            if (angle_index[bus] >= 0) {
                angle_scale[angle_index[bus]] = std::min(
                    1e4, std::max(angle_scale[angle_index[bus]], maximum));
            }
        }
        const auto maximum_voltage_derivative = [](const FlowValues& value) {
            return std::max({
                std::abs(value.pf), std::abs(value.qf),
                std::abs(value.pt), std::abs(value.qt),
            });
        };
        voltage_scale[data.branches[i].from] = std::min(
            1e4, std::max(
                voltage_scale[data.branches[i].from],
                maximum_voltage_derivative(voltage_derivative[i].from)));
        voltage_scale[data.branches[i].to] = std::min(
            1e4, std::max(
                voltage_scale[data.branches[i].to],
                maximum_voltage_derivative(voltage_derivative[i].to)));
    }
    for (int i = 0; i < static_cast<int>(data.shunts.size()); ++i) {
        const auto& shunt = data.shunts[i];
        const double derivative_magnitude = 2.0 * output.state.vm[shunt.bus] *
            std::max(
                std::abs(shunt.gs),
                std::abs(effective_shunt_susceptance(
                    data, output.state, i)));
        voltage_scale[shunt.bus] = std::min(
            1e4, std::max(
                voltage_scale[shunt.bus], derivative_magnitude));
    }

    const int angle_up_offset = 0;
    const int angle_down_offset = angle_up_offset + angle_count;
    const int voltage_variable_count = include_reactive ? nb : 0;
    const int reactive_generator_variable_count = include_reactive ? ng : 0;
    const int voltage_up_offset = angle_down_offset + angle_count;
    const int voltage_down_offset =
        voltage_up_offset + voltage_variable_count;
    const int generator_up_offset =
        voltage_down_offset + voltage_variable_count;
    const int generator_down_offset = generator_up_offset + ng;
    const int reactive_generator_up_offset =
        generator_down_offset + ng;
    const int reactive_generator_down_offset =
        reactive_generator_up_offset + reactive_generator_variable_count;
    const int load_up_offset =
        reactive_generator_down_offset + reactive_generator_variable_count;
    const int load_down_offset = load_up_offset + nd;
    const HighsInt column_count = load_down_offset + nd;
    output.column_count = static_cast<int>(column_count);
    std::vector<double> lower(static_cast<std::size_t>(column_count), 0.0);
    std::vector<double> upper(static_cast<std::size_t>(column_count), 0.0);
    std::vector<double> cost(static_cast<std::size_t>(column_count), 0.0);
    for (int i = 0; i < angle_count; ++i) {
        const double scaled_trust = angle_trust_radius * angle_scale[i];
        upper[angle_up_offset + i] = scaled_trust;
        upper[angle_down_offset + i] = scaled_trust;
        cost[angle_up_offset + i] = 1e-3 / angle_scale[i];
        cost[angle_down_offset + i] = 1e-3 / angle_scale[i];
    }
    for (int i = 0; include_reactive && i < nb; ++i) {
        const double required_upward_move = std::max(
            0.0, data.buses[i].vmin - output.state.vm[i]);
        const double required_downward_move = std::max(
            0.0, output.state.vm[i] - data.buses[i].vmax);
        if (required_upward_move > voltage_trust_radius + 1e-12 ||
            required_downward_move > voltage_trust_radius + 1e-12) {
            output.status = "reference_voltage_outside_trust_region";
            output.wall_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - wall_start).count();
            return output;
        }
        const double upward_room = std::min(
            voltage_trust_radius,
            std::max(0.0, data.buses[i].vmax - output.state.vm[i]));
        const double downward_room = std::min(
            voltage_trust_radius,
            std::max(0.0, output.state.vm[i] - data.buses[i].vmin));
        upper[voltage_up_offset + i] =
            upward_room * voltage_scale[i];
        upper[voltage_down_offset + i] =
            downward_room * voltage_scale[i];
        // The split movement variables describe the LP point itself.  If the
        // nonlinear reference is outside a source voltage bound, require the
        // LP to return inside that bound rather than clamping only after the
        // balance rows have been solved.  Post-solve clamping alone changes
        // branch injections and invalidates the claimed linearized point.
        lower[voltage_up_offset + i] =
            required_upward_move * voltage_scale[i];
        lower[voltage_down_offset + i] =
            required_downward_move * voltage_scale[i];
        cost[voltage_up_offset + i] = 1e-3 / voltage_scale[i];
        cost[voltage_down_offset + i] = 1e-3 / voltage_scale[i];
    }
    for (int i = 0; i < ng; ++i) {
        if (!active[i]) {
            continue;
        }
        upper[generator_up_offset + i] = std::max(
            0.0, p_upper[i] - output.state.pg[i]);
        upper[generator_down_offset + i] = std::max(
            0.0, output.state.pg[i] - p_lower[i]);
        cost[generator_up_offset + i] = 1.0;
        cost[generator_down_offset + i] = 1.0;
        if (include_reactive) {
            upper[reactive_generator_up_offset + i] = std::max(
                0.0, q_upper[i] - output.state.qg[i]);
            upper[reactive_generator_down_offset + i] = std::max(
                0.0, output.state.qg[i] - q_lower[i]);
            cost[reactive_generator_up_offset + i] = 0.1;
            cost[reactive_generator_down_offset + i] = 0.1;
        }
    }
    for (int i = 0; i < nd; ++i) {
        upper[load_up_offset + i] = std::max(
            0.0, load_factor_upper[i] - load_factor[i]);
        upper[load_down_offset + i] = std::max(
            0.0, load_factor[i] - load_factor_lower[i]);
        const double load_movement_cost = 10.0 * std::max({
            std::abs(data.loads[i].pd_nominal),
            std::abs(data.loads[i].qd_nominal), 1e-3});
        cost[load_up_offset + i] = load_movement_cost;
        cost[load_down_offset + i] = load_movement_cost;
    }

    const auto append_angle = [&](SparseRow& row, int bus, double value) {
        const int index = angle_index[bus];
        if (index < 0 || std::abs(value) <= 1e-14) {
            return;
        }
        const double scaled = value / angle_scale[index];
        append(row, angle_up_offset + index, scaled);
        append(row, angle_down_offset + index, -scaled);
    };
    const auto append_voltage = [&](SparseRow& row, int bus, double value) {
        if (!include_reactive || std::abs(value) <= 1e-14) {
            return;
        }
        const double scaled = value / voltage_scale[bus];
        append(row, voltage_up_offset + bus, scaled);
        append(row, voltage_down_offset + bus, -scaled);
    };

    std::vector<SparseRow> active_balance_rows(static_cast<std::size_t>(nb));
    std::vector<SparseRow> reactive_balance_rows(static_cast<std::size_t>(nb));
    std::vector<double> current_active_balance(
        static_cast<std::size_t>(nb), 0.0);
    std::vector<double> current_reactive_balance(
        static_cast<std::size_t>(nb), 0.0);
    for (int bus = 0; bus < nb; ++bus) {
        for (int generator : data.buses[bus].generators) {
            current_active_balance[bus] -= output.state.pg[generator];
            current_reactive_balance[bus] -= output.state.qg[generator];
            if (active[generator]) {
                append(
                    active_balance_rows[bus],
                    generator_up_offset + generator, -1.0);
                append(
                    active_balance_rows[bus],
                    generator_down_offset + generator, 1.0);
                if (include_reactive) {
                    append(
                        reactive_balance_rows[bus],
                        reactive_generator_up_offset + generator, -1.0);
                    append(
                        reactive_balance_rows[bus],
                        reactive_generator_down_offset + generator, 1.0);
                }
            }
        }
        for (int load_index : data.buses[bus].loads) {
            const auto& load = data.loads[load_index];
            current_active_balance[bus] +=
                load.pd_nominal * load_factor[load_index];
            current_reactive_balance[bus] +=
                load.qd_nominal * load_factor[load_index];
            append(
                active_balance_rows[bus], load_up_offset + load_index,
                load.pd_nominal);
            append(
                active_balance_rows[bus], load_down_offset + load_index,
                -load.pd_nominal);
            if (include_reactive) {
                append(
                    reactive_balance_rows[bus], load_up_offset + load_index,
                    load.qd_nominal);
                append(
                    reactive_balance_rows[bus], load_down_offset + load_index,
                    -load.qd_nominal);
            }
        }
        for (int shunt_index : data.buses[bus].shunts) {
            const auto& shunt = data.shunts[shunt_index];
            const double vm = output.state.vm[bus];
            const double bs = effective_shunt_susceptance(
                data, output.state, shunt_index);
            current_active_balance[bus] += shunt.gs * vm * vm;
            current_reactive_balance[bus] -= bs * vm * vm;
            append_voltage(
                active_balance_rows[bus], bus, 2.0 * shunt.gs * vm);
            append_voltage(
                reactive_balance_rows[bus], bus, -2.0 * bs * vm);
        }
    }
    for (int i = 0; i < nl; ++i) {
        if (i == outaged_branch || data.branches[i].status == 0) {
            continue;
        }
        const auto& branch = data.branches[i];
        current_active_balance[branch.from] += output.state.pf[i];
        current_reactive_balance[branch.from] += output.state.qf[i];
        current_active_balance[branch.to] += output.state.pt[i];
        current_reactive_balance[branch.to] += output.state.qt[i];
        append_angle(
            active_balance_rows[branch.from], branch.from,
            derivative[i].pf);
        append_angle(
            active_balance_rows[branch.from], branch.to,
            -derivative[i].pf);
        append_voltage(
            active_balance_rows[branch.from], branch.from,
            voltage_derivative[i].from.pf);
        append_voltage(
            active_balance_rows[branch.from], branch.to,
            voltage_derivative[i].to.pf);
        append_angle(
            reactive_balance_rows[branch.from], branch.from,
            derivative[i].qf);
        append_angle(
            reactive_balance_rows[branch.from], branch.to,
            -derivative[i].qf);
        append_voltage(
            reactive_balance_rows[branch.from], branch.from,
            voltage_derivative[i].from.qf);
        append_voltage(
            reactive_balance_rows[branch.from], branch.to,
            voltage_derivative[i].to.qf);
        append_angle(
            active_balance_rows[branch.to], branch.from,
            derivative[i].pt);
        append_angle(
            active_balance_rows[branch.to], branch.to,
            -derivative[i].pt);
        append_voltage(
            active_balance_rows[branch.to], branch.from,
            voltage_derivative[i].from.pt);
        append_voltage(
            active_balance_rows[branch.to], branch.to,
            voltage_derivative[i].to.pt);
        append_angle(
            reactive_balance_rows[branch.to], branch.from,
            derivative[i].qt);
        append_angle(
            reactive_balance_rows[branch.to], branch.to,
            -derivative[i].qt);
        append_voltage(
            reactive_balance_rows[branch.to], branch.from,
            voltage_derivative[i].from.qt);
        append_voltage(
            reactive_balance_rows[branch.to], branch.to,
            voltage_derivative[i].to.qt);
    }

    std::vector<SparseRow> rows;
    rows.reserve(static_cast<std::size_t>(2 * nb) + 4096);
    for (int bus = 0; bus < nb; ++bus) {
        auto active_row = std::move(active_balance_rows[bus]);
        active_row.lower =
            -balance_slack_limit - current_active_balance[bus];
        active_row.upper =
            balance_slack_limit - current_active_balance[bus];
        normalize(active_row);
        rows.push_back(std::move(active_row));
        if (include_reactive) {
            auto reactive_row = std::move(reactive_balance_rows[bus]);
            reactive_row.lower =
                -balance_slack_limit - current_reactive_balance[bus];
            reactive_row.upper =
                balance_slack_limit - current_reactive_balance[bus];
            normalize(reactive_row);
            rows.push_back(std::move(reactive_row));
        }
    }

    const auto append_security_row = [&](
        double current,
        double angle_derivative,
        double voltage_from_derivative,
        double voltage_to_derivative,
        const Branch& branch,
        double target) {
        if (current_security_rows_only &&
            current >= -target && current <= target) {
            return;
        }
        const double angle_range = angle_trust_radius *
            ((angle_index[branch.from] >= 0 ? 1.0 : 0.0) +
             (angle_index[branch.to] >= 0 ? 1.0 : 0.0));
        const double movement =
            std::abs(angle_derivative) * angle_range +
            voltage_trust_radius *
                (std::abs(voltage_from_derivative) +
                 std::abs(voltage_to_derivative));
        if (current - movement >= -target &&
            current + movement <= target) {
            return;
        }
        SparseRow row;
        row.lower = -target - current;
        row.upper = target - current;
        append_angle(row, branch.from, angle_derivative);
        append_angle(row, branch.to, -angle_derivative);
        append_voltage(
            row, branch.from, voltage_from_derivative);
        append_voltage(
            row, branch.to, voltage_to_derivative);
        normalize(row);
        rows.push_back(std::move(row));
        ++output.branch_security_row_count;
    };
    const auto append_apparent_row = [&](
        double p,
        double q,
        double angle_tangent,
        double voltage_from_tangent,
        double voltage_to_tangent,
        const Branch& branch,
        double limit_squared) {
        const double angle_range = angle_trust_radius *
            ((angle_index[branch.from] >= 0 ? 1.0 : 0.0) +
             (angle_index[branch.to] >= 0 ? 1.0 : 0.0));
        const double current_squared = p * p + q * q;
        if (current_security_rows_only &&
            current_squared <= limit_squared) {
            return;
        }
        const double movement =
            std::abs(angle_tangent) * angle_range +
            voltage_trust_radius *
                (std::abs(voltage_from_tangent) +
                 std::abs(voltage_to_tangent));
        if (current_squared + movement <= limit_squared) {
            return;
        }
        SparseRow row;
        row.lower = -kHighsInf;
        row.upper = limit_squared - current_squared;
        append_angle(row, branch.from, angle_tangent);
        append_angle(row, branch.to, -angle_tangent);
        append_voltage(
            row, branch.from, voltage_from_tangent);
        append_voltage(
            row, branch.to, voltage_to_tangent);
        normalize(row);
        rows.push_back(std::move(row));
        ++output.branch_security_row_count;
    };

    for (int i = 0; i < nl; ++i) {
        if (i == outaged_branch || data.branches[i].status == 0) {
            continue;
        }
        const auto& branch = data.branches[i];
        const double box_margin = std::min(
            1e-4, std::max(0.0, 0.01 * branch.rate_c));
        const double box_target = std::max(0.0, branch.rate_c - box_margin);
        append_security_row(
            output.state.pf[i], derivative[i].pf,
            voltage_derivative[i].from.pf,
            voltage_derivative[i].to.pf, branch, box_target);
        append_security_row(
            output.state.qf[i], derivative[i].qf,
            voltage_derivative[i].from.qf,
            voltage_derivative[i].to.qf, branch, box_target);
        append_security_row(
            output.state.pt[i], derivative[i].pt,
            voltage_derivative[i].from.pt,
            voltage_derivative[i].to.pt, branch, box_target);
        append_security_row(
            output.state.qt[i], derivative[i].qt,
            voltage_derivative[i].from.qt,
            voltage_derivative[i].to.qt, branch, box_target);
        const double from_scale = branch.transformer
            ? 1.0 + data.sm_vio_limit
            : output.state.vm[branch.from] + data.sm_vio_limit;
        const double to_scale = branch.transformer
            ? 1.0 + data.sm_vio_limit
            : output.state.vm[branch.to] + data.sm_vio_limit;
        append_apparent_row(
            output.state.pf[i], output.state.qf[i],
            2.0 * (output.state.pf[i] * derivative[i].pf +
                   output.state.qf[i] * derivative[i].qf),
            2.0 * (output.state.pf[i] *
                       voltage_derivative[i].from.pf +
                   output.state.qf[i] *
                       voltage_derivative[i].from.qf) -
                (branch.transformer ? 0.0 :
                    2.0 * branch.rate_c * branch.rate_c * from_scale),
            2.0 * (output.state.pf[i] *
                       voltage_derivative[i].to.pf +
                   output.state.qf[i] *
                       voltage_derivative[i].to.qf),
            branch,
            branch.rate_c * branch.rate_c * from_scale * from_scale);
        append_apparent_row(
            output.state.pt[i], output.state.qt[i],
            2.0 * (output.state.pt[i] * derivative[i].pt +
                   output.state.qt[i] * derivative[i].qt),
            2.0 * (output.state.pt[i] *
                       voltage_derivative[i].from.pt +
                   output.state.qt[i] *
                       voltage_derivative[i].from.qt),
            2.0 * (output.state.pt[i] *
                       voltage_derivative[i].to.pt +
                   output.state.qt[i] *
                       voltage_derivative[i].to.qt) -
                (branch.transformer ? 0.0 :
                    2.0 * branch.rate_c * branch.rate_c * to_scale),
            branch,
            branch.rate_c * branch.rate_c * to_scale * to_scale);

        const double source_delta =
            contingency.effective_base_state().va[branch.from] -
            contingency.effective_base_state().va[branch.to];
        if (source_delta < branch.angmin || source_delta > branch.angmax) {
            continue;
        }
        const double current_delta =
            output.state.va[branch.from] - output.state.va[branch.to];
        if (current_security_rows_only &&
            current_delta >= branch.angmin &&
            current_delta <= branch.angmax) {
            continue;
        }
        const double delta_range = angle_trust_radius *
            ((angle_index[branch.from] >= 0 ? 1.0 : 0.0) +
             (angle_index[branch.to] >= 0 ? 1.0 : 0.0));
        if (current_delta - delta_range >= branch.angmin &&
            current_delta + delta_range <= branch.angmax) {
            continue;
        }
        SparseRow row;
        row.lower = branch.angmin - current_delta;
        row.upper = branch.angmax - current_delta;
        append_angle(row, branch.from, 1.0);
        append_angle(row, branch.to, -1.0);
        normalize(row);
        rows.push_back(std::move(row));
        ++output.branch_security_row_count;
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
    for (auto& row : rows) {
        normalize(row);
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
    highs.setOptionValue("run_crossover", "off");
    highs.setOptionValue("small_matrix_value", 1e-12);
    highs.setOptionValue("time_limit", time_limit_seconds);
    highs.setOptionValue("primal_feasibility_tolerance", 1e-8);
    highs.setOptionValue("dual_feasibility_tolerance", 1e-8);
    highs.setOptionValue("ipm_optimality_tolerance", 1e-6);
    const char* solver_override =
        std::getenv("GRAVITYX_ACTIVE_REPAIR_SOLVER");
    output.solver = solver_override != nullptr
        ? std::string(solver_override) : "simplex";
    highs.setOptionValue("solver", output.solver);
    if (current_security_rows_only && output.solver == "simplex") {
        // This repair is a feasibility oracle.  Dual simplex can spend the
        // complete short deadline proving movement optimality without ever
        // exposing a primal-feasible basis; primal simplex prioritizes the
        // candidate that the nonlinear validator actually needs.
        output.simplex_strategy = 4;
        highs.setOptionValue("simplex_strategy", output.simplex_strategy);
    }
    if (highs.addVars(column_count, lower.data(), upper.data()) !=
            HighsStatus::kOk ||
        highs.changeColsCost(
            0, column_count - 1, cost.data()) != HighsStatus::kOk ||
        highs.addRows(
            static_cast<HighsInt>(rows.size()), row_lower.data(),
            row_upper.data(), static_cast<HighsInt>(indices.size()),
            starts.data(), indices.data(), values.data()) != HighsStatus::kOk) {
        output.status = "model_construction_failed";
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
    output.num_primal_infeasibilities =
        static_cast<int>(info.num_primal_infeasibilities);
    output.max_primal_infeasibility = info.max_primal_infeasibility;
    output.iterations = static_cast<int>(
        info.ipm_iteration_count > 0
            ? info.ipm_iteration_count : info.simplex_iteration_count);
    output.objective = info.objective_function_value;
    output.status = highs.modelStatusToString(model_status);
    const bool shape_valid = solution.value_valid &&
        solution.col_value.size() == static_cast<std::size_t>(column_count);
    if (shape_valid) {
        output.finite_solution_values = true;
        for (HighsInt column = 0; column < column_count; ++column) {
            const double value = solution.col_value[column];
            if (!std::isfinite(value)) {
                output.finite_solution_values = false;
                break;
            }
            output.maximum_column_violation = std::max(
                output.maximum_column_violation,
                std::max(
                    std::max(0.0, lower[column] - value),
                    std::max(0.0, value - upper[column])));
        }
        if (output.finite_solution_values) {
            for (const auto& row : rows) {
                output.maximum_linearized_violation = std::max(
                    output.maximum_linearized_violation,
                    row_violation(row, solution.col_value));
            }
        }
    }
    const bool optimal = run_status != HighsStatus::kError &&
        model_status == HighsModelStatus::kOptimal && shape_valid;
    const bool feasible_nonoptimal = run_status != HighsStatus::kError &&
        shape_valid && info.valid &&
        info.primal_solution_status == kSolutionStatusFeasible &&
        info.num_primal_infeasibilities == 0 &&
        std::isfinite(info.max_primal_infeasibility) &&
        info.max_primal_infeasibility <= 1e-8;
    const bool independently_feasible = shape_valid &&
        output.finite_solution_values &&
        output.maximum_column_violation <= 1e-7 &&
        output.maximum_linearized_violation <= 1e-7;
    output.accepted_feasible_nonoptimal =
        (feasible_nonoptimal || independently_feasible) && !optimal;
    output.success = optimal || feasible_nonoptimal || independently_feasible;
    if (output.success) {
        for (int bus = 0; bus < nb; ++bus) {
            const int index = angle_index[bus];
            if (index >= 0) {
                const double change =
                    (solution.col_value[angle_up_offset + index] -
                     solution.col_value[angle_down_offset + index]) /
                    angle_scale[index];
                output.state.va[bus] += change;
                output.maximum_angle_change = std::max(
                    output.maximum_angle_change, std::abs(change));
            }
            if (include_reactive) {
                const double voltage_change =
                    (solution.col_value[voltage_up_offset + bus] -
                     solution.col_value[voltage_down_offset + bus]) /
                    voltage_scale[bus];
                output.state.vm[bus] = std::clamp(
                    output.state.vm[bus] + voltage_change,
                    data.buses[bus].vmin, data.buses[bus].vmax);
                output.maximum_voltage_change = std::max(
                    output.maximum_voltage_change,
                    std::abs(voltage_change));
            }
        }
        for (int i = 0; i < ng; ++i) {
            if (!active[i]) {
                output.state.pg[i] = 0.0;
                output.state.qg[i] = 0.0;
                continue;
            }
            const double prior = output.state.pg[i];
            output.state.pg[i] +=
                solution.col_value[generator_up_offset + i] -
                solution.col_value[generator_down_offset + i];
            output.state.pg[i] = std::clamp(
                output.state.pg[i], p_lower[i], p_upper[i]);
            output.maximum_generation_change = std::max(
                output.maximum_generation_change,
                std::abs(output.state.pg[i] - prior));
            if (include_reactive) {
                const double prior_reactive = output.state.qg[i];
                output.state.qg[i] +=
                    solution.col_value[reactive_generator_up_offset + i] -
                    solution.col_value[reactive_generator_down_offset + i];
                output.state.qg[i] = std::clamp(
                    output.state.qg[i], q_lower[i], q_upper[i]);
                output.maximum_reactive_generation_change = std::max(
                    output.maximum_reactive_generation_change,
                    std::abs(output.state.qg[i] - prior_reactive));
            }
        }
        for (int i = 0; i < nd; ++i) {
            const double prior = load_factor[i];
            const double current = std::clamp(
                prior + solution.col_value[load_up_offset + i] -
                    solution.col_value[load_down_offset + i],
                load_factor_lower[i], load_factor_upper[i]);
            output.state.demand_factor[i] = current;
            output.maximum_load_change = std::max(
                output.maximum_load_change,
                std::max(
                    std::abs(data.loads[i].pd_nominal),
                    std::abs(data.loads[i].qd_nominal)) *
                    std::abs(current - prior));
        }
    }
    output.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_start).count();
    return output;
}

}  // namespace gravityx
