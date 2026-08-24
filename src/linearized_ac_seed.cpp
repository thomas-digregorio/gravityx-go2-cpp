#include <highs/Highs.h>

#include "gravityx/linearized_ac_seed.hpp"

#include "gravityx/state_io.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace gravityx {
namespace {

struct FlowValues {
    double pf{};
    double qf{};
    double pt{};
    double qt{};
};

struct AffineFlow {
    double constant{};
    double vm_from{};
    double vm_to{};
    double va_from{};
    double va_to{};
};

struct LinearizedBranch {
    AffineFlow pf;
    AffineFlow qf;
    AffineFlow pt;
    AffineFlow qt;
};

struct SparseRow {
    double lower{};
    double upper{};
    std::vector<std::pair<HighsInt, double>> entries;
};

FlowValues branch_flows(
    const Branch& branch,
    double vm_from,
    double vm_to,
    double va_from,
    double va_to) {
    const double denominator = branch.r * branch.r + branch.x * branch.x;
    const double g = denominator > 1e-20 ? branch.r / denominator : 0.0;
    const double b = denominator > 1e-20 ? -branch.x / denominator : 0.0;
    const double tm = branch.tap;
    if (std::abs(tm) <= 1e-12) {
        throw std::runtime_error(
            "zero tap ratio in linearized AC seed: " + branch.source_key);
    }
    const double tm2 = tm * tm;
    const double tr = tm * std::cos(branch.shift);
    const double ti = tm * std::sin(branch.shift);
    const double delta_ft = va_from - va_to;
    const double delta_tf = -delta_ft;
    const double cross_ft = vm_from * vm_to;
    const double from_g_self = branch.transformer
        ? g / tm2 + branch.g_fr : (g + branch.g_fr) / tm2;
    const double from_b_self = branch.transformer
        ? b / tm2 + branch.b_fr : (b + branch.b_fr) / tm2;

    FlowValues value;
    value.pf = from_g_self * vm_from * vm_from
        + ((-g * tr + b * ti) / tm2) * cross_ft * std::cos(delta_ft)
        + ((-b * tr - g * ti) / tm2) * cross_ft * std::sin(delta_ft);
    value.qf = -from_b_self * vm_from * vm_from
        - ((-b * tr + g * ti) / tm2) * cross_ft * std::cos(delta_ft)
        + ((-g * tr + b * ti) / tm2) * cross_ft * std::sin(delta_ft);
    value.pt = (g + branch.g_to) * vm_to * vm_to
        + ((-g * tr - b * ti) / tm2) * cross_ft * std::cos(delta_tf)
        + ((-b * tr + g * ti) / tm2) * cross_ft * std::sin(delta_tf);
    value.qt = -(b + branch.b_to) * vm_to * vm_to
        - ((-b * tr + g * ti) / tm2) * cross_ft * std::cos(delta_tf)
        + ((-g * tr - b * ti) / tm2) * cross_ft * std::sin(delta_tf);
    return value;
}

double flow_component(const FlowValues& value, int component) {
    switch (component) {
        case 0: return value.pf;
        case 1: return value.qf;
        case 2: return value.pt;
        case 3: return value.qt;
        default: throw std::runtime_error("invalid branch-flow component");
    }
}

LinearizedBranch linearize_branch(
    const Branch& branch,
    double vm_from,
    double vm_to,
    double va_from,
    double va_to) {
    constexpr double h = 1e-6;
    const auto base = branch_flows(
        branch, vm_from, vm_to, va_from, va_to);
    std::array<std::array<double, 4>, 4> derivative{};
    const std::array<double, 4> point{vm_from, vm_to, va_from, va_to};
    for (int variable = 0; variable < 4; ++variable) {
        auto plus = point;
        auto minus = point;
        plus[variable] += h;
        minus[variable] -= h;
        const auto plus_flow = branch_flows(
            branch, plus[0], plus[1], plus[2], plus[3]);
        const auto minus_flow = branch_flows(
            branch, minus[0], minus[1], minus[2], minus[3]);
        for (int component = 0; component < 4; ++component) {
            derivative[component][variable] =
                (flow_component(plus_flow, component)
                    - flow_component(minus_flow, component)) / (2.0 * h);
        }
    }
    LinearizedBranch result;
    std::array<AffineFlow*, 4> output{
        &result.pf, &result.qf, &result.pt, &result.qt};
    for (int component = 0; component < 4; ++component) {
        auto& affine = *output[component];
        affine.vm_from = derivative[component][0];
        affine.vm_to = derivative[component][1];
        affine.va_from = derivative[component][2];
        affine.va_to = derivative[component][3];
        affine.constant = flow_component(base, component)
            - affine.vm_from * vm_from
            - affine.vm_to * vm_to
            - affine.va_from * va_from
            - affine.va_to * va_to;
    }
    return result;
}

void append(SparseRow& row, HighsInt column, double value) {
    if (std::abs(value) > 1e-14) {
        row.entries.emplace_back(column, value);
    }
}

void normalize(SparseRow& row) {
    std::sort(row.entries.begin(), row.entries.end(),
        [](const auto& left, const auto& right) {
            return left.first < right.first;
        });
    std::vector<std::pair<HighsInt, double>> combined;
    for (const auto& entry : row.entries) {
        if (!combined.empty() && combined.back().first == entry.first) {
            combined.back().second += entry.second;
        } else {
            combined.push_back(entry);
        }
    }
    combined.erase(
        std::remove_if(combined.begin(), combined.end(),
            [](const auto& entry) { return std::abs(entry.second) <= 1e-14; }),
        combined.end());
    row.entries = std::move(combined);
}

}  // namespace

nlohmann::json LinearizedAcSeedResult::to_json(bool include_state) const {
    nlohmann::json value = {
        {"success", success},
        {"wall_seconds", wall_seconds},
        {"model_status", model_status},
        {"status", status},
        {"iterations", iterations},
        {"objective", objective},
    };
    if (include_state) {
        value["state"] = ac_state_to_json(state);
    }
    return value;
}

LinearizedAcSeedResult solve_linearized_ac_seed(
    const CaseData& data,
    const AcState& reference,
    const std::vector<int>& commitment,
    double balance_slack_limit,
    const std::optional<ContingencyContext>& contingency) {
    const auto wall_start = std::chrono::steady_clock::now();
    const int nb = static_cast<int>(data.buses.size());
    const int ng = static_cast<int>(data.generators.size());
    const int nd = static_cast<int>(data.loads.size());
    const int nl = static_cast<int>(data.branches.size());
    if (commitment.size() != data.generators.size() ||
        reference.vm.size() != data.buses.size() ||
        reference.va.size() != data.buses.size() ||
        reference.pg.size() != data.generators.size() ||
        reference.qg.size() != data.generators.size() ||
        reference.demand_factor.size() != data.loads.size()) {
        throw std::runtime_error("linearized AC seed received invalid dimensions");
    }
    if (!(balance_slack_limit > 0.0 && balance_slack_limit < 0.5)) {
        throw std::runtime_error("linearized AC balance slack must be in (0, 0.5)");
    }
    if (contingency) {
        const auto& base = contingency->base_state;
        if (base.pg.size() != data.generators.size() ||
            base.demand_factor.size() != data.loads.size() ||
            base.va.size() != data.buses.size()) {
            throw std::runtime_error(
                "linearized AC contingency base state has invalid dimensions");
        }
    }
    const int outaged_generator = contingency
        ? contingency->outaged_generator : -1;
    const int outaged_branch = contingency
        ? contingency->outaged_branch : -1;

    const int vm_offset = 0;
    const int va_offset = vm_offset + nb;
    const int pg_offset = va_offset + nb;
    const int qg_offset = pg_offset + ng;
    const int demand_offset = qg_offset + ng;
    const int p_pos_offset = demand_offset + nd;
    const int p_neg_offset = p_pos_offset + nb;
    const int q_pos_offset = p_neg_offset + nb;
    const int q_neg_offset = q_pos_offset + nb;
    const int dpg_offset = q_neg_offset + nb;
    const int dqg_offset = dpg_offset + ng;
    const int dload_offset = dqg_offset + ng;
    const int column_count = dload_offset + nd;

    std::vector<double> lower(column_count, 0.0);
    std::vector<double> upper(column_count, kHighsInf);
    std::vector<double> cost(column_count, 0.0);
    for (int i = 0; i < nb; ++i) {
        lower[vm_offset + i] = data.buses[i].vmin;
        upper[vm_offset + i] = data.buses[i].vmax;
        lower[va_offset + i] = -10.0;
        upper[va_offset + i] = 10.0;
        if (data.buses[i].type == 3) {
            lower[va_offset + i] = 0.0;
            upper[va_offset + i] = 0.0;
        }
        upper[p_pos_offset + i] = balance_slack_limit;
        upper[p_neg_offset + i] = balance_slack_limit;
        upper[q_pos_offset + i] = balance_slack_limit;
        upper[q_neg_offset + i] = balance_slack_limit;
        cost[p_pos_offset + i] = 1e6;
        cost[p_neg_offset + i] = 1e6;
        cost[q_pos_offset + i] = 1e6;
        cost[q_neg_offset + i] = 1e6;
    }
    for (int i = 0; i < ng; ++i) {
        const auto& gen = data.generators[i];
        if (commitment[i] == 0 || i == outaged_generator) {
            lower[pg_offset + i] = upper[pg_offset + i] = 0.0;
            lower[qg_offset + i] = upper[qg_offset + i] = 0.0;
        } else {
            if (contingency) {
                lower[pg_offset + i] = std::max(
                    gen.pmin,
                    contingency->base_state.pg[i]
                        - data.delta_r_ctg * gen.prdmaxctg);
                upper[pg_offset + i] = std::min(
                    gen.pmax,
                    contingency->base_state.pg[i]
                        + data.delta_r_ctg * gen.prumaxctg);
            } else {
                const double previous =
                    gen.status_prev == 0 ? gen.pmin : gen.pg_prev;
                lower[pg_offset + i] = std::max(
                    gen.pmin, previous - data.delta_r * gen.prdmax);
                upper[pg_offset + i] = std::min(
                    gen.pmax, previous + data.delta_r * gen.prumax);
            }
            if (lower[pg_offset + i] > upper[pg_offset + i] + 1e-12) {
                throw std::runtime_error(
                    "empty generator interval in linearized AC seed: "
                    + gen.source_key);
            }
            lower[qg_offset + i] = gen.qmin;
            upper[qg_offset + i] = gen.qmax;
        }
        cost[dpg_offset + i] = 1.0;
        cost[dqg_offset + i] = 0.1;
    }
    for (int i = 0; i < nd; ++i) {
        const auto& load = data.loads[i];
        lower[demand_offset + i] = load.tmin;
        upper[demand_offset + i] = load.tmax;
        if (std::abs(load.pd_nominal) > 1e-12) {
            const double previous = contingency
                ? load.pd_nominal
                    * contingency->base_state.demand_factor[i]
                : load.pd_prev;
            const double down = contingency
                ? load.prdmaxctg * data.delta_r_ctg
                : load.prdmax * data.delta_r;
            const double up = contingency
                ? load.prumaxctg * data.delta_r_ctg
                : load.prumax * data.delta_r;
            lower[demand_offset + i] = std::max(
                lower[demand_offset + i],
                (previous - down) / load.pd_nominal);
            upper[demand_offset + i] = std::min(
                upper[demand_offset + i],
                (previous + up) / load.pd_nominal);
        }
        if (lower[demand_offset + i] > upper[demand_offset + i] + 1e-12) {
            throw std::runtime_error(
                "empty load interval in linearized AC seed: "
                + load.source_key);
        }
        cost[dload_offset + i] = 1.0;
    }

    std::vector<LinearizedBranch> linearized(nl);
    for (int i = 0; i < nl; ++i) {
        if (i == outaged_branch) {
            continue;
        }
        const auto& branch = data.branches[i];
        linearized[i] = linearize_branch(
            branch,
            reference.vm[branch.from], reference.vm[branch.to],
            reference.va[branch.from], reference.va[branch.to]);
    }

    std::vector<SparseRow> rows;
    rows.reserve(2 * nb + 5 * nl + 4 * ng + 2 * nd);
    for (int bus = 0; bus < nb; ++bus) {
        SparseRow active;
        SparseRow reactive;
        double active_constant = 0.0;
        double reactive_constant = 0.0;
        for (int branch_index : data.buses[bus].branches_from) {
            if (branch_index == outaged_branch) {
                continue;
            }
            const auto& branch = data.branches[branch_index];
            const auto& p = linearized[branch_index].pf;
            const auto& q = linearized[branch_index].qf;
            active_constant += p.constant;
            reactive_constant += q.constant;
            append(active, vm_offset + branch.from, p.vm_from);
            append(active, vm_offset + branch.to, p.vm_to);
            append(active, va_offset + branch.from, p.va_from);
            append(active, va_offset + branch.to, p.va_to);
            append(reactive, vm_offset + branch.from, q.vm_from);
            append(reactive, vm_offset + branch.to, q.vm_to);
            append(reactive, va_offset + branch.from, q.va_from);
            append(reactive, va_offset + branch.to, q.va_to);
        }
        for (int branch_index : data.buses[bus].branches_to) {
            if (branch_index == outaged_branch) {
                continue;
            }
            const auto& branch = data.branches[branch_index];
            const auto& p = linearized[branch_index].pt;
            const auto& q = linearized[branch_index].qt;
            active_constant += p.constant;
            reactive_constant += q.constant;
            append(active, vm_offset + branch.from, p.vm_from);
            append(active, vm_offset + branch.to, p.vm_to);
            append(active, va_offset + branch.from, p.va_from);
            append(active, va_offset + branch.to, p.va_to);
            append(reactive, vm_offset + branch.from, q.vm_from);
            append(reactive, vm_offset + branch.to, q.vm_to);
            append(reactive, va_offset + branch.from, q.va_from);
            append(reactive, va_offset + branch.to, q.va_to);
        }
        for (int generator : data.buses[bus].generators) {
            append(active, pg_offset + generator, -1.0);
            append(reactive, qg_offset + generator, -1.0);
        }
        for (int load : data.buses[bus].loads) {
            append(active, demand_offset + load, data.loads[load].pd_nominal);
            append(reactive, demand_offset + load, data.loads[load].qd_nominal);
        }
        double gs = 0.0;
        double bs = 0.0;
        for (int shunt : data.buses[bus].shunts) {
            gs += data.shunts[shunt].gs;
            bs += data.shunts[shunt].bs;
        }
        active_constant -= gs * reference.vm[bus] * reference.vm[bus];
        reactive_constant += bs * reference.vm[bus] * reference.vm[bus];
        append(active, vm_offset + bus, 2.0 * gs * reference.vm[bus]);
        append(reactive, vm_offset + bus, -2.0 * bs * reference.vm[bus]);
        append(active, p_pos_offset + bus, -1.0);
        append(active, p_neg_offset + bus, 1.0);
        append(reactive, q_pos_offset + bus, -1.0);
        append(reactive, q_neg_offset + bus, 1.0);
        active.lower = active.upper = -active_constant;
        reactive.lower = reactive.upper = -reactive_constant;
        normalize(active);
        normalize(reactive);
        rows.push_back(std::move(active));
        rows.push_back(std::move(reactive));
    }

    for (int i = 0; i < nl; ++i) {
        if (i == outaged_branch) {
            continue;
        }
        const auto& branch = data.branches[i];
        const std::array<const AffineFlow*, 4> flows{
            &linearized[i].pf, &linearized[i].qf,
            &linearized[i].pt, &linearized[i].qt};
        for (const auto* flow : flows) {
            SparseRow row;
            row.lower = -branch.rate_a - flow->constant;
            row.upper = branch.rate_a - flow->constant;
            append(row, vm_offset + branch.from, flow->vm_from);
            append(row, vm_offset + branch.to, flow->vm_to);
            append(row, va_offset + branch.from, flow->va_from);
            append(row, va_offset + branch.to, flow->va_to);
            normalize(row);
            rows.push_back(std::move(row));
        }
        const double source_delta = contingency
            ? contingency->base_state.va[branch.from]
                - contingency->base_state.va[branch.to]
            : data.buses[branch.from].va_start
                - data.buses[branch.to].va_start;
        if (source_delta >= branch.angmin && source_delta <= branch.angmax) {
            SparseRow angle;
            angle.lower = branch.angmin;
            angle.upper = branch.angmax;
            append(angle, va_offset + branch.from, 1.0);
            append(angle, va_offset + branch.to, -1.0);
            rows.push_back(std::move(angle));
        }
    }

    for (int i = 0; i < ng; ++i) {
        SparseRow pg_positive;
        pg_positive.lower = -kHighsInf;
        pg_positive.upper = reference.pg[i];
        append(pg_positive, pg_offset + i, 1.0);
        append(pg_positive, dpg_offset + i, -1.0);
        rows.push_back(std::move(pg_positive));
        SparseRow pg_negative;
        pg_negative.lower = -kHighsInf;
        pg_negative.upper = -reference.pg[i];
        append(pg_negative, pg_offset + i, -1.0);
        append(pg_negative, dpg_offset + i, -1.0);
        rows.push_back(std::move(pg_negative));
        SparseRow qg_positive;
        qg_positive.lower = -kHighsInf;
        qg_positive.upper = reference.qg[i];
        append(qg_positive, qg_offset + i, 1.0);
        append(qg_positive, dqg_offset + i, -1.0);
        rows.push_back(std::move(qg_positive));
        SparseRow qg_negative;
        qg_negative.lower = -kHighsInf;
        qg_negative.upper = -reference.qg[i];
        append(qg_negative, qg_offset + i, -1.0);
        append(qg_negative, dqg_offset + i, -1.0);
        rows.push_back(std::move(qg_negative));
    }
    for (int i = 0; i < nd; ++i) {
        const double coefficient = data.loads[i].pd_nominal;
        const double reference_power = coefficient * reference.demand_factor[i];
        SparseRow positive;
        positive.lower = -kHighsInf;
        positive.upper = reference_power;
        append(positive, demand_offset + i, coefficient);
        append(positive, dload_offset + i, -1.0);
        rows.push_back(std::move(positive));
        SparseRow negative;
        negative.lower = -kHighsInf;
        negative.upper = -reference_power;
        append(negative, demand_offset + i, -coefficient);
        append(negative, dload_offset + i, -1.0);
        rows.push_back(std::move(negative));
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
        for (const auto& [column, value] : row.entries) {
            indices.push_back(column);
            values.push_back(value);
        }
        starts.push_back(static_cast<HighsInt>(indices.size()));
    }

    Highs highs;
    highs.setOptionValue("output_flag", false);
    highs.setOptionValue("threads", 1);
    highs.setOptionValue("presolve", "on");
    highs.setOptionValue("solver", "ipm");
    highs.setOptionValue("run_crossover", "off");
    highs.setOptionValue("time_limit", 30.0);
    highs.setOptionValue("primal_feasibility_tolerance", 1e-8);
    highs.setOptionValue("dual_feasibility_tolerance", 1e-8);
    if (highs.addVars(column_count, lower.data(), upper.data()) != HighsStatus::kOk ||
        highs.changeColsCost(0, column_count - 1, cost.data()) != HighsStatus::kOk ||
        highs.addRows(
            static_cast<HighsInt>(rows.size()),
            row_lower.data(), row_upper.data(),
            static_cast<HighsInt>(indices.size()), starts.data(),
            indices.data(), values.data()) != HighsStatus::kOk) {
        throw std::runtime_error("failed to construct the linearized AC HiGHS model");
    }
    const auto run_status = highs.run();
    const auto model_status = highs.getModelStatus();
    const auto& solution = highs.getSolution();
    const auto& info = highs.getInfo();

    LinearizedAcSeedResult output;
    output.model_status = static_cast<int>(model_status);
    output.status = highs.modelStatusToString(model_status);
    output.iterations = static_cast<int>(
        info.ipm_iteration_count > 0
            ? info.ipm_iteration_count
            : info.simplex_iteration_count);
    output.objective = info.objective_function_value;
    output.success = run_status == HighsStatus::kOk &&
        model_status == HighsModelStatus::kOptimal && solution.value_valid &&
        solution.col_value.size() == static_cast<std::size_t>(column_count);
    output.state = reference;
    if (output.success) {
        for (int i = 0; i < nb; ++i) {
            output.state.vm[i] = solution.col_value[vm_offset + i];
            output.state.va[i] = solution.col_value[va_offset + i];
        }
        for (int i = 0; i < ng; ++i) {
            output.state.pg[i] = solution.col_value[pg_offset + i];
            output.state.qg[i] = solution.col_value[qg_offset + i];
        }
        if (outaged_generator >= 0) {
            output.state.pg[outaged_generator] = 0.0;
            output.state.qg[outaged_generator] = 0.0;
        }
        for (int i = 0; i < nd; ++i) {
            output.state.demand_factor[i] =
                solution.col_value[demand_offset + i];
        }
    }
    output.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_start).count();
    return output;
}

}  // namespace gravityx
