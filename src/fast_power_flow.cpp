#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include "gravityx/fast_power_flow.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <utility>

namespace gravityx {
namespace {

using Complex = std::complex<double>;
using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;
using YRows = std::vector<std::map<int, Complex>>;

constexpr double kAllocationTolerance = 1e-7;

void add_admittance(YRows& rows, int row, int column, Complex value) {
    rows[row][column] += value;
}

YRows build_ybus(const CaseData& data, int outaged_branch) {
    YRows rows(data.buses.size());
    for (const auto& shunt : data.shunts) {
        add_admittance(rows, shunt.bus, shunt.bus, {shunt.gs, shunt.bs});
    }
    for (int i = 0; i < static_cast<int>(data.branches.size()); ++i) {
        if (i == outaged_branch) {
            continue;
        }
        const auto& branch = data.branches[i];
        const double denominator = branch.r * branch.r + branch.x * branch.x;
        const double g = denominator > 1e-20 ? branch.r / denominator : 0.0;
        const double b = denominator > 1e-20 ? -branch.x / denominator : 0.0;
        const Complex series(g, b);
        const double tm = branch.tap;
        if (std::abs(tm) <= 1e-12) {
            throw std::runtime_error("zero tap ratio in fast power flow: " + branch.source_key);
        }
        const double tm2 = tm * tm;
        const Complex from_shunt(branch.g_fr, branch.b_fr);
        const Complex to_shunt(branch.g_to, branch.b_to);
        const Complex rotation = std::polar(1.0, branch.shift);
        const Complex yff = branch.transformer
            ? series / tm2 + from_shunt
            : (series + from_shunt) / tm2;
        const Complex yft = -series * rotation / tm;
        const Complex ytf = -series * std::conj(rotation) / tm;
        const Complex ytt = series + to_shunt;
        add_admittance(rows, branch.from, branch.from, yff);
        add_admittance(rows, branch.from, branch.to, yft);
        add_admittance(rows, branch.to, branch.from, ytf);
        add_admittance(rows, branch.to, branch.to, ytt);
    }
    return rows;
}

void network_injections(
    const YRows& ybus,
    const std::vector<double>& vm,
    const std::vector<double>& va,
    std::vector<double>& p,
    std::vector<double>& q) {
    p.assign(vm.size(), 0.0);
    q.assign(vm.size(), 0.0);
    for (int i = 0; i < static_cast<int>(vm.size()); ++i) {
        Complex current(0.0, 0.0);
        const Complex vi = std::polar(vm[i], va[i]);
        for (const auto& [j, admittance] : ybus[i]) {
            current += admittance * std::polar(vm[j], va[j]);
        }
        const Complex injection = vi * std::conj(current);
        p[i] = injection.real();
        q[i] = injection.imag();
    }
}

double mismatch_norm(
    const std::vector<double>& p_spec,
    const std::vector<double>& q_spec,
    const std::vector<double>& p,
    const std::vector<double>& q,
    const std::vector<int>& angle_index,
    const std::vector<int>& voltage_index) {
    double result = 0.0;
    for (int i = 0; i < static_cast<int>(p.size()); ++i) {
        if (angle_index[i] >= 0) {
            result = std::max(result, std::abs(p_spec[i] - p[i]));
        }
        if (voltage_index[i] >= 0) {
            result = std::max(result, std::abs(q_spec[i] - q[i]));
        }
    }
    return result;
}

struct NewtonResult {
    bool converged{};
    int iterations{};
    std::string failure_reason;
};

NewtonResult run_newton(
    const CaseData& data,
    const YRows& ybus,
    const std::vector<bool>& slack,
    const std::vector<bool>& pq,
    const std::vector<double>& p_spec,
    const std::vector<double>& q_spec,
    int max_iterations,
    double tolerance,
    std::vector<double>& vm,
    std::vector<double>& va) {
    const int nb = static_cast<int>(data.buses.size());
    std::vector<int> angle_index(nb, -1);
    std::vector<int> voltage_index(nb, -1);
    int angle_count = 0;
    int voltage_count = 0;
    for (int i = 0; i < nb; ++i) {
        if (!slack[i]) {
            angle_index[i] = angle_count++;
        }
        if (pq[i]) {
            voltage_index[i] = voltage_count++;
        }
    }
    const int dimension = angle_count + voltage_count;
    if (dimension <= 0) {
        return {true, 0, {}};
    }

    std::vector<double> p;
    std::vector<double> q;
    for (int iteration = 0; iteration <= max_iterations; ++iteration) {
        network_injections(ybus, vm, va, p, q);
        const double norm = mismatch_norm(
            p_spec, q_spec, p, q, angle_index, voltage_index);
        if (norm <= tolerance) {
            return {true, iteration, {}};
        }
        if (iteration == max_iterations) {
            return {false, iteration,
                "Newton iteration limit; residual=" + std::to_string(norm)};
        }

        Eigen::VectorXd mismatch(dimension);
        mismatch.setZero();
        for (int i = 0; i < nb; ++i) {
            if (angle_index[i] >= 0) {
                mismatch[angle_index[i]] = p_spec[i] - p[i];
            }
            if (voltage_index[i] >= 0) {
                mismatch[angle_count + voltage_index[i]] = q_spec[i] - q[i];
            }
        }

        std::vector<Triplet> entries;
        entries.reserve(static_cast<std::size_t>(dimension) * 8);
        for (int i = 0; i < nb; ++i) {
            const double vi = vm[i];
            if (vi <= 1e-8) {
                return {false, iteration, "nonpositive voltage magnitude"};
            }
            const auto diagonal = ybus[i].find(i);
            const double gii = diagonal == ybus[i].end() ? 0.0 : diagonal->second.real();
            const double bii = diagonal == ybus[i].end() ? 0.0 : diagonal->second.imag();

            if (angle_index[i] >= 0) {
                const int row = angle_index[i];
                entries.emplace_back(row, angle_index[i], -q[i] - bii * vi * vi);
                if (voltage_index[i] >= 0) {
                    entries.emplace_back(
                        row, angle_count + voltage_index[i], p[i] / vi + gii * vi);
                }
                for (const auto& [j, admittance] : ybus[i]) {
                    if (j == i) {
                        continue;
                    }
                    const double delta = va[i] - va[j];
                    const double gij = admittance.real();
                    const double bij = admittance.imag();
                    if (angle_index[j] >= 0) {
                        entries.emplace_back(
                            row, angle_index[j],
                            vi * vm[j] * (gij * std::sin(delta) - bij * std::cos(delta)));
                    }
                    if (voltage_index[j] >= 0) {
                        entries.emplace_back(
                            row, angle_count + voltage_index[j],
                            vi * (gij * std::cos(delta) + bij * std::sin(delta)));
                    }
                }
            }

            if (voltage_index[i] >= 0) {
                const int row = angle_count + voltage_index[i];
                if (angle_index[i] >= 0) {
                    entries.emplace_back(row, angle_index[i], p[i] - gii * vi * vi);
                }
                entries.emplace_back(
                    row, angle_count + voltage_index[i], q[i] / vi - bii * vi);
                for (const auto& [j, admittance] : ybus[i]) {
                    if (j == i) {
                        continue;
                    }
                    const double delta = va[i] - va[j];
                    const double gij = admittance.real();
                    const double bij = admittance.imag();
                    if (angle_index[j] >= 0) {
                        entries.emplace_back(
                            row, angle_index[j],
                            -vi * vm[j] * (gij * std::cos(delta) + bij * std::sin(delta)));
                    }
                    if (voltage_index[j] >= 0) {
                        entries.emplace_back(
                            row, angle_count + voltage_index[j],
                            vi * (gij * std::sin(delta) - bij * std::cos(delta)));
                    }
                }
            }
        }

        SparseMatrix jacobian(dimension, dimension);
        jacobian.setFromTriplets(entries.begin(), entries.end());
        Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>> factorization;
        factorization.analyzePattern(jacobian);
        factorization.factorize(jacobian);
        if (factorization.info() != Eigen::Success) {
            return {false, iteration, "sparse Jacobian factorization failed"};
        }
        const Eigen::VectorXd step = factorization.solve(mismatch);
        if (factorization.info() != Eigen::Success || !step.allFinite()) {
            return {false, iteration, "sparse Newton solve failed"};
        }

        double scale = 1.0;
        double max_angle_step = 0.0;
        double max_voltage_step = 0.0;
        for (int i = 0; i < nb; ++i) {
            if (angle_index[i] >= 0) {
                max_angle_step = std::max(max_angle_step, std::abs(step[angle_index[i]]));
            }
            if (voltage_index[i] >= 0) {
                const double dv = step[angle_count + voltage_index[i]];
                max_voltage_step = std::max(max_voltage_step, std::abs(dv));
            }
        }
        if (max_angle_step > 0.5) {
            scale = std::min(scale, 0.5 / max_angle_step);
        }
        if (max_voltage_step > 0.1) {
            scale = std::min(scale, 0.1 / max_voltage_step);
        }
        scale = std::clamp(scale, 1e-4, 1.0);

        const auto prior_vm = vm;
        const auto prior_va = va;
        const double prior_norm = norm;
        bool accepted = false;
        for (int line_search = 0; line_search < 12; ++line_search) {
            vm = prior_vm;
            va = prior_va;
            for (int i = 0; i < nb; ++i) {
                if (angle_index[i] >= 0) {
                    va[i] += scale * step[angle_index[i]];
                }
                if (voltage_index[i] >= 0) {
                    vm[i] += scale * step[angle_count + voltage_index[i]];
                }
            }
            network_injections(ybus, vm, va, p, q);
            const double trial_norm = mismatch_norm(
                p_spec, q_spec, p, q, angle_index, voltage_index);
            if (std::isfinite(trial_norm) && trial_norm < prior_norm) {
                accepted = true;
                break;
            }
            scale *= 0.5;
        }
        if (!accepted) {
            vm = prior_vm;
            va = prior_va;
            return {false, iteration, "Newton line search stagnated"};
        }
    }
    return {false, max_iterations, "Newton failed"};
}

bool allocate_total(
    const std::vector<int>& indices,
    const std::vector<double>& lower,
    const std::vector<double>& upper,
    const std::vector<double>& preferred,
    double target,
    std::vector<double>& values) {
    double total_lower = 0.0;
    double total_upper = 0.0;
    for (int index : indices) {
        total_lower += lower[index];
        total_upper += upper[index];
    }
    if (target < total_lower - kAllocationTolerance ||
        target > total_upper + kAllocationTolerance) {
        return false;
    }
    for (int index : indices) {
        values[index] = std::clamp(preferred[index], lower[index], upper[index]);
    }
    for (int pass = 0; pass < 8; ++pass) {
        double current = 0.0;
        for (int index : indices) {
            current += values[index];
        }
        const double difference = target - current;
        if (std::abs(difference) <= kAllocationTolerance) {
            return true;
        }
        double room = 0.0;
        for (int index : indices) {
            room += difference > 0.0
                ? std::max(0.0, upper[index] - values[index])
                : std::max(0.0, values[index] - lower[index]);
        }
        if (room <= kAllocationTolerance) {
            break;
        }
        for (int index : indices) {
            const double individual_room = difference > 0.0
                ? std::max(0.0, upper[index] - values[index])
                : std::max(0.0, values[index] - lower[index]);
            values[index] += difference * individual_room / room;
            values[index] = std::clamp(values[index], lower[index], upper[index]);
        }
    }
    double current = 0.0;
    for (int index : indices) {
        current += values[index];
    }
    return std::abs(target - current) <= 1e-7;
}

std::vector<double> interpolation_weights(
    const std::vector<PwlPoint>& points,
    double target) {
    if (points.empty()) {
        throw std::runtime_error("cannot interpolate an empty PWL curve");
    }
    std::vector<double> weights(points.size(), 0.0);
    if (target <= points.front().mw) {
        weights.front() = 1.0;
        return weights;
    }
    if (target >= points.back().mw) {
        weights.back() = 1.0;
        return weights;
    }
    for (std::size_t i = 1; i < points.size(); ++i) {
        if (target <= points[i].mw) {
            const double width = points[i].mw - points[i - 1].mw;
            if (std::abs(width) <= 1e-14) {
                weights[i] = 1.0;
            } else {
                const double right = (target - points[i - 1].mw) / width;
                weights[i - 1] = 1.0 - right;
                weights[i] = right;
            }
            return weights;
        }
    }
    weights.back() = 1.0;
    return weights;
}

void append_weights(
    std::vector<double>& destination,
    const std::vector<PwlPoint>& points,
    double target) {
    const auto weights = interpolation_weights(points, target);
    destination.insert(destination.end(), weights.begin(), weights.end());
}

void compute_branch_flows(
    const CaseData& data,
    int outaged_branch,
    AcState& state) {
    const int nl = static_cast<int>(data.branches.size());
    state.pf.assign(nl, 0.0);
    state.qf.assign(nl, 0.0);
    state.pt.assign(nl, 0.0);
    state.qt.assign(nl, 0.0);
    state.sm_slack.assign(nl, 0.0);
    for (int i = 0; i < nl; ++i) {
        if (i == outaged_branch) {
            continue;
        }
        const auto& branch = data.branches[i];
        const double denominator = branch.r * branch.r + branch.x * branch.x;
        const double g = denominator > 1e-20 ? branch.r / denominator : 0.0;
        const double b = denominator > 1e-20 ? -branch.x / denominator : 0.0;
        const double tm = branch.tap;
        const double tm2 = tm * tm;
        const double tr = tm * std::cos(branch.shift);
        const double ti = tm * std::sin(branch.shift);
        const int f = branch.from;
        const int t = branch.to;
        const double cross_cos_ft = state.vm[f] * state.vm[t]
            * std::cos(state.va[f] - state.va[t]);
        const double cross_sin_ft = state.vm[f] * state.vm[t]
            * std::sin(state.va[f] - state.va[t]);
        const double from_g_self = branch.transformer
            ? g / tm2 + branch.g_fr : (g + branch.g_fr) / tm2;
        const double from_b_self = branch.transformer
            ? b / tm2 + branch.b_fr : (b + branch.b_fr) / tm2;
        state.pf[i] = from_g_self * state.vm[f] * state.vm[f]
            + ((-g * tr + b * ti) / tm2) * cross_cos_ft
            + ((-b * tr - g * ti) / tm2) * cross_sin_ft;
        state.qf[i] = -from_b_self * state.vm[f] * state.vm[f]
            - ((-b * tr - g * ti) / tm2) * cross_cos_ft
            + ((-g * tr + b * ti) / tm2) * cross_sin_ft;
        const double cross_cos_tf = state.vm[t] * state.vm[f]
            * std::cos(state.va[t] - state.va[f]);
        const double cross_sin_tf = state.vm[t] * state.vm[f]
            * std::sin(state.va[t] - state.va[f]);
        state.pt[i] = (g + branch.g_to) * state.vm[t] * state.vm[t]
            + ((-g * tr - b * ti) / tm2) * cross_cos_tf
            + ((-b * tr + g * ti) / tm2) * cross_sin_tf;
        state.qt[i] = -(b + branch.b_to) * state.vm[t] * state.vm[t]
            - ((-b * tr + g * ti) / tm2) * cross_cos_tf
            + ((-g * tr - b * ti) / tm2) * cross_sin_tf;

        const double from_magnitude = std::hypot(state.pf[i], state.qf[i]);
        const double to_magnitude = std::hypot(state.pt[i], state.qt[i]);
        if (branch.rate_a > 1e-12) {
            const double from_required = branch.transformer
                ? from_magnitude / branch.rate_a - 1.0
                : from_magnitude / branch.rate_a - state.vm[f];
            const double to_required = branch.transformer
                ? to_magnitude / branch.rate_a - 1.0
                : to_magnitude / branch.rate_a - state.vm[t];
            state.sm_slack[i] = std::max({0.0, from_required, to_required});
        }
    }
}

std::vector<std::vector<int>> connected_components(
    const CaseData& data,
    int outaged_branch) {
    const int nb = static_cast<int>(data.buses.size());
    std::vector<std::vector<int>> adjacency(nb);
    for (int i = 0; i < static_cast<int>(data.branches.size()); ++i) {
        if (i == outaged_branch) {
            continue;
        }
        const auto& branch = data.branches[i];
        adjacency[branch.from].push_back(branch.to);
        adjacency[branch.to].push_back(branch.from);
    }
    std::vector<bool> visited(nb, false);
    std::vector<std::vector<int>> result;
    for (int start = 0; start < nb; ++start) {
        if (visited[start]) {
            continue;
        }
        std::vector<int> component;
        std::queue<int> frontier;
        frontier.push(start);
        visited[start] = true;
        while (!frontier.empty()) {
            const int bus = frontier.front();
            frontier.pop();
            component.push_back(bus);
            for (int neighbor : adjacency[bus]) {
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    frontier.push(neighbor);
                }
            }
        }
        result.push_back(std::move(component));
    }
    return result;
}

}  // namespace

nlohmann::json FastPowerFlowResult::to_json() const {
    return {
        {"converged", converged},
        {"feasible", feasible},
        {"newton_iterations", newton_iterations},
        {"active_redispatch_passes", active_redispatch_passes},
        {"reactive_limit_passes", reactive_limit_passes},
        {"wall_seconds", wall_seconds},
        {"failure_reason", failure_reason},
        {"validation", validation.to_json()},
    };
}

FastContingencyPowerFlow::FastContingencyPowerFlow(
    const CaseData& data,
    const AcState& base_state,
    std::vector<int> commitment,
    FastPowerFlowOptions options)
    : data_(data),
      base_state_(base_state),
      commitment_(std::move(commitment)),
      options_(options) {
    if (commitment_.size() != data_.generators.size()) {
        throw std::runtime_error("fast power flow commitment has wrong length");
    }
    if (base_state_.vm.size() != data_.buses.size() ||
        base_state_.pg.size() != data_.generators.size()) {
        throw std::runtime_error("fast power flow base state has wrong dimensions");
    }
}

FastPowerFlowResult FastContingencyPowerFlow::solve(
    const Contingency& contingency) const {
    const auto wall_start = std::chrono::steady_clock::now();
    FastPowerFlowResult output;
    const int nb = static_cast<int>(data_.buses.size());
    const int ng = static_cast<int>(data_.generators.size());
    const int outaged_generator = contingency.type == ContingencyType::Generator
        ? contingency.component : -1;
    const int outaged_branch = contingency.type == ContingencyType::Branch
        ? contingency.component : -1;

    std::vector<bool> active(ng, false);
    std::vector<double> p_lower(ng, 0.0), p_upper(ng, 0.0);
    std::vector<double> q_lower(ng, 0.0), q_upper(ng, 0.0);
    std::vector<double> pg(ng, 0.0), qg(ng, 0.0);
    for (int i = 0; i < ng; ++i) {
        active[i] = commitment_[i] == 1 && i != outaged_generator;
        if (!active[i]) {
            continue;
        }
        const auto& gen = data_.generators[i];
        p_lower[i] = std::max(
            gen.pmin, base_state_.pg[i] - data_.delta_r_ctg * gen.prdmaxctg);
        p_upper[i] = std::min(
            gen.pmax, base_state_.pg[i] + data_.delta_r_ctg * gen.prumaxctg);
        q_lower[i] = gen.qmin;
        q_upper[i] = gen.qmax;
        if (p_lower[i] > p_upper[i] + 1e-12 || q_lower[i] > q_upper[i] + 1e-12) {
            output.failure_reason = "empty generator contingency interval";
            output.wall_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - wall_start).count();
            return output;
        }
        pg[i] = std::clamp(base_state_.pg[i], p_lower[i], p_upper[i]);
        qg[i] = std::clamp(base_state_.qg[i], q_lower[i], q_upper[i]);
    }

    const auto components = connected_components(data_, outaged_branch);
    std::vector<bool> slack(nb, false);
    std::vector<bool> pq(nb, true);
    std::vector<int> component_of(nb, -1);
    std::vector<std::vector<int>> active_at_bus(nb);
    for (int i = 0; i < ng; ++i) {
        if (active[i]) {
            active_at_bus[data_.generators[i].bus].push_back(i);
        }
    }
    for (int c = 0; c < static_cast<int>(components.size()); ++c) {
        int reference = -1;
        for (int bus : components[c]) {
            component_of[bus] = c;
            if (!active_at_bus[bus].empty()) {
                pq[bus] = false;
                if (data_.buses[bus].type == 3) {
                    reference = bus;
                }
            }
        }
        if (reference < 0) {
            for (int bus : components[c]) {
                if (!active_at_bus[bus].empty()) {
                    reference = bus;
                    break;
                }
            }
        }
        if (reference < 0) {
            output.failure_reason = "network island has no committed generator";
            output.wall_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - wall_start).count();
            return output;
        }
        slack[reference] = true;
    }

    for (const auto& component : components) {
        std::vector<int> participants;
        double target = 0.0;
        double lower = 0.0;
        double upper = 0.0;
        for (int bus : component) {
            for (int gen : data_.buses[bus].generators) {
                if (commitment_[gen] == 1) {
                    target += base_state_.pg[gen];
                }
                if (active[gen]) {
                    participants.push_back(gen);
                    lower += p_lower[gen];
                    upper += p_upper[gen];
                }
            }
        }
        if (!participants.empty()) {
            target = std::clamp(target, lower, upper);
            allocate_total(
                participants, p_lower, p_upper,
                base_state_.pg, target, pg);
        }
    }

    const YRows ybus = build_ybus(data_, outaged_branch);
    std::vector<double> vm = base_state_.vm;
    std::vector<double> va = base_state_.va;
    for (int i = 0; i < nb; ++i) {
        vm[i] = std::clamp(vm[i], data_.buses[i].vmin, data_.buses[i].vmax);
        if (data_.buses[i].type == 3) {
            va[i] = 0.0;
        }
    }

    std::vector<double> load_p(nb, 0.0), load_q(nb, 0.0);
    for (int i = 0; i < static_cast<int>(data_.loads.size()); ++i) {
        const auto& load = data_.loads[i];
        load_p[load.bus] += load.pd_nominal * base_state_.demand_factor[i];
        load_q[load.bus] += load.qd_nominal * base_state_.demand_factor[i];
    }
    std::vector<double> fixed_q_bus(nb, 0.0);
    for (int bus = 0; bus < nb; ++bus) {
        for (int gen : active_at_bus[bus]) {
            fixed_q_bus[bus] += qg[gen];
        }
    }

    bool fully_solved = false;
    bool approximate_candidate = false;
    std::vector<double> p_network, q_network;
    for (int active_pass = 0;
         active_pass <= options_.max_active_redispatch_passes && !fully_solved;
         ++active_pass) {
        output.active_redispatch_passes = active_pass;
        bool reactive_solved = false;
        for (int q_pass = 0; q_pass <= options_.max_reactive_limit_passes; ++q_pass) {
            output.reactive_limit_passes = std::max(output.reactive_limit_passes, q_pass);
            std::vector<double> p_spec(nb, 0.0), q_spec(nb, 0.0);
            for (int bus = 0; bus < nb; ++bus) {
                p_spec[bus] = -load_p[bus];
                q_spec[bus] = fixed_q_bus[bus] - load_q[bus];
                for (int gen : active_at_bus[bus]) {
                    p_spec[bus] += pg[gen];
                }
            }
            const NewtonResult newton = run_newton(
                data_, ybus, slack, pq, p_spec, q_spec,
                options_.max_newton_iterations, options_.newton_tolerance,
                vm, va);
            output.newton_iterations += newton.iterations;
            if (!newton.converged) {
                output.failure_reason = newton.failure_reason;
                approximate_candidate = true;
                reactive_solved = true;
                break;
            }
            bool projected_voltage = false;
            for (int bus = 0; bus < nb; ++bus) {
                if (!pq[bus]) {
                    continue;
                }
                const double projected = std::clamp(
                    vm[bus], data_.buses[bus].vmin, data_.buses[bus].vmax);
                if (std::abs(projected - vm[bus]) > 1e-9) {
                    vm[bus] = projected;
                    pq[bus] = false;
                    projected_voltage = true;
                }
            }
            if (projected_voltage) {
                continue;
            }
            network_injections(ybus, vm, va, p_network, q_network);

            bool added_q_limit = false;
            for (int bus = 0; bus < nb; ++bus) {
                if (active_at_bus[bus].empty() || pq[bus]) {
                    continue;
                }
                double lower = 0.0;
                double upper = 0.0;
                for (int gen : active_at_bus[bus]) {
                    lower += q_lower[gen];
                    upper += q_upper[gen];
                }
                const double required = q_network[bus] + load_q[bus];
                if (required < lower - 1e-7 || required > upper + 1e-7) {
                    fixed_q_bus[bus] = std::clamp(required, lower, upper);
                    pq[bus] = true;
                    added_q_limit = true;
                }
            }
            if (!added_q_limit) {
                reactive_solved = true;
                break;
            }
        }
        if (!reactive_solved) {
            output.failure_reason = "reactive-limit pass limit";
            break;
        }
        if (approximate_candidate) {
            fully_solved = true;
            break;
        }

        network_injections(ybus, vm, va, p_network, q_network);
        bool adjusted_active = false;
        for (int c = 0; c < static_cast<int>(components.size()); ++c) {
            int slack_bus = -1;
            for (int bus : components[c]) {
                if (slack[bus]) {
                    slack_bus = bus;
                    break;
                }
            }
            const double required_slack = p_network[slack_bus] + load_p[slack_bus];
            double slack_lower = 0.0;
            double slack_upper = 0.0;
            for (int gen : active_at_bus[slack_bus]) {
                slack_lower += p_lower[gen];
                slack_upper += p_upper[gen];
            }
            if (required_slack >= slack_lower - 1e-8 &&
                required_slack <= slack_upper + 1e-8) {
                continue;
            }
            const double slack_target = std::clamp(
                required_slack, slack_lower, slack_upper);
            const double change_needed = required_slack - slack_target;
            std::vector<int> participants;
            for (int bus : components[c]) {
                if (bus == slack_bus) {
                    continue;
                }
                participants.insert(
                    participants.end(), active_at_bus[bus].begin(), active_at_bus[bus].end());
            }
            double current = 0.0;
            for (int gen : participants) {
                current += pg[gen];
            }
            if (!allocate_total(
                    participants, p_lower, p_upper, pg,
                    current + change_needed, pg)) {
                output.failure_reason = "insufficient corrective active-power range; validating bounded-residual candidate";
                approximate_candidate = true;
                adjusted_active = false;
                break;
            }
            adjusted_active = true;
        }
        if (!adjusted_active) {
            fully_solved = true;
        }
    }
    if (!fully_solved) {
        approximate_candidate = true;
        fully_solved = true;
        output.failure_reason = "active-redispatch pass limit; validating bounded-residual candidate";
    }

    network_injections(ybus, vm, va, p_network, q_network);
    if (!approximate_candidate) {
        for (int bus = 0; bus < nb; ++bus) {
            if (active_at_bus[bus].empty()) {
                continue;
            }
            double required_p = 0.0;
            for (int gen : active_at_bus[bus]) {
                required_p += pg[gen];
            }
            if (slack[bus]) {
                required_p = p_network[bus] + load_p[bus];
            }
            const double required_q = pq[bus]
                ? fixed_q_bus[bus] : q_network[bus] + load_q[bus];
            auto proposed_pg = pg;
            auto proposed_qg = qg;
            if (!allocate_total(
                    active_at_bus[bus], p_lower, p_upper,
                    base_state_.pg, required_p, proposed_pg) ||
                !allocate_total(
                    active_at_bus[bus], q_lower, q_upper,
                    base_state_.qg, required_q, proposed_qg)) {
                approximate_candidate = true;
                output.failure_reason = "generator allocation failed; validating bounded-residual candidate";
                break;
            }
            pg = std::move(proposed_pg);
            qg = std::move(proposed_qg);
        }
    }

    for (const auto& component : components) {
        int angle_reference = -1;
        for (int bus : component) {
            if (data_.buses[bus].type == 3) {
                angle_reference = bus;
                break;
            }
        }
        if (angle_reference >= 0) {
            const double offset = va[angle_reference];
            for (int bus : component) {
                va[bus] -= offset;
            }
        }
    }
    for (int bus = 0; bus < nb; ++bus) {
        vm[bus] = std::clamp(
            vm[bus], data_.buses[bus].vmin, data_.buses[bus].vmax);
    }

    AcState state = base_state_;
    state.vm = std::move(vm);
    state.va = std::move(va);
    state.pg = std::move(pg);
    state.qg = std::move(qg);
    compute_branch_flows(data_, outaged_branch, state);
    state.p_delta.assign(nb, 0.0);
    state.q_delta.assign(nb, 0.0);
    for (int bus = 0; bus < nb; ++bus) {
        double p_balance = 0.0;
        double q_balance = 0.0;
        for (int branch : data_.buses[bus].branches_from) {
            p_balance += state.pf[branch];
            q_balance += state.qf[branch];
        }
        for (int branch : data_.buses[bus].branches_to) {
            p_balance += state.pt[branch];
            q_balance += state.qt[branch];
        }
        for (int gen : data_.buses[bus].generators) {
            p_balance -= state.pg[gen];
            q_balance -= state.qg[gen];
        }
        for (int load : data_.buses[bus].loads) {
            p_balance += data_.loads[load].pd_nominal * state.demand_factor[load];
            q_balance += data_.loads[load].qd_nominal * state.demand_factor[load];
        }
        for (int shunt : data_.buses[bus].shunts) {
            p_balance += data_.shunts[shunt].gs * state.vm[bus] * state.vm[bus];
            q_balance -= data_.shunts[shunt].bs * state.vm[bus] * state.vm[bus];
        }
        state.p_delta[bus] = std::min(0.5, std::abs(p_balance));
        state.q_delta[bus] = std::min(0.5, std::abs(q_balance));
    }

    const int nd = static_cast<int>(data_.loads.size());
    std::vector<double> load_power_lower(nd, 0.0);
    std::vector<double> load_power_upper(nd, 0.0);
    std::vector<double> load_power_preferred(nd, 0.0);
    std::vector<double> load_reactive_lower(nd, 0.0);
    std::vector<double> load_reactive_upper(nd, 0.0);
    std::vector<double> load_reactive_preferred(nd, 0.0);
    for (int i = 0; i < nd; ++i) {
        const auto& load = data_.loads[i];
        const double previous = load.pd_nominal * base_state_.demand_factor[i];
        const double factor_lower = std::abs(load.pd_nominal) <= 1e-12
            ? load.tmin
            : std::max(load.tmin,
                (previous - load.prdmaxctg * data_.delta_r_ctg) / load.pd_nominal);
        const double factor_upper = std::abs(load.pd_nominal) <= 1e-12
            ? load.tmax
            : std::min(load.tmax,
                (previous + load.prumaxctg * data_.delta_r_ctg) / load.pd_nominal);
        load_power_lower[i] = std::min(
            load.pd_nominal * factor_lower, load.pd_nominal * factor_upper);
        load_power_upper[i] = std::max(
            load.pd_nominal * factor_lower, load.pd_nominal * factor_upper);
        load_power_preferred[i] = previous;
        load_reactive_lower[i] = std::min(
            load.qd_nominal * factor_lower, load.qd_nominal * factor_upper);
        load_reactive_upper[i] = std::max(
            load.qd_nominal * factor_lower, load.qd_nominal * factor_upper);
        load_reactive_preferred[i] =
            load.qd_nominal * base_state_.demand_factor[i];
    }

    const auto refresh_pwl_weights = [&]() {
        state.gen_lambda.clear();
        for (int i = 0; i < ng; ++i) {
            if (!active[i]) {
                continue;
            }
            const auto points = active_pwl_points(
                data_.generators[i].cost, data_.generators[i].ncost,
                p_lower[i], p_upper[i]);
            append_weights(state.gen_lambda, points, state.pg[i]);
        }
        state.load_lambda.clear();
        for (int i = 0; i < nd; ++i) {
            const auto& load = data_.loads[i];
            const auto points = active_pwl_points(
                load.cost, load.ncost, load.pd_min, load.pd_max);
            append_weights(
                state.load_lambda, points,
                load.pd_nominal * state.demand_factor[i]);
        }
    };

    const auto refresh_network_fields = [&]() {
        compute_branch_flows(data_, outaged_branch, state);
        const auto compute_balances = [&]() {
            std::pair<std::vector<double>, std::vector<double>> balances{
                std::vector<double>(nb, 0.0), std::vector<double>(nb, 0.0)};
            auto& p_balance = balances.first;
            auto& q_balance = balances.second;
            for (int bus = 0; bus < nb; ++bus) {
                for (int branch : data_.buses[bus].branches_from) {
                    p_balance[bus] += state.pf[branch];
                    q_balance[bus] += state.qf[branch];
                }
                for (int branch : data_.buses[bus].branches_to) {
                    p_balance[bus] += state.pt[branch];
                    q_balance[bus] += state.qt[branch];
                }
                for (int gen : data_.buses[bus].generators) {
                    p_balance[bus] -= state.pg[gen];
                    q_balance[bus] -= state.qg[gen];
                }
                for (int load : data_.buses[bus].loads) {
                    p_balance[bus] += data_.loads[load].pd_nominal
                        * state.demand_factor[load];
                    q_balance[bus] += data_.loads[load].qd_nominal
                        * state.demand_factor[load];
                }
                for (int shunt : data_.buses[bus].shunts) {
                    p_balance[bus] += data_.shunts[shunt].gs
                        * state.vm[bus] * state.vm[bus];
                    q_balance[bus] -= data_.shunts[shunt].bs
                        * state.vm[bus] * state.vm[bus];
                }
            }
            return balances;
        };

        auto balances = compute_balances();
        for (int bus = 0; bus < nb; ++bus) {
            const double q_excess = std::copysign(
                std::max(0.0, std::abs(balances.second[bus]) - 0.45),
                balances.second[bus]);
            if (std::abs(q_excess) > 0.0 && !active_at_bus[bus].empty()) {
                double current = 0.0;
                for (int gen : active_at_bus[bus]) {
                    current += state.qg[gen];
                }
                auto proposed = state.qg;
                if (allocate_total(
                        active_at_bus[bus], q_lower, q_upper,
                        base_state_.qg, current + q_excess, proposed)) {
                    state.qg = std::move(proposed);
                }
            }
        }

        balances = compute_balances();
        std::vector<double> load_reactive(nd, 0.0);
        for (int i = 0; i < nd; ++i) {
            load_reactive[i] = data_.loads[i].qd_nominal * state.demand_factor[i];
        }
        for (int bus = 0; bus < nb; ++bus) {
            const double q_excess = std::copysign(
                std::max(0.0, std::abs(balances.second[bus]) - 0.45),
                balances.second[bus]);
            if (std::abs(q_excess) <= 0.0 || data_.buses[bus].loads.empty()) {
                continue;
            }
            double current = 0.0;
            double total_lower = 0.0;
            double total_upper = 0.0;
            for (int load : data_.buses[bus].loads) {
                current += load_reactive[load];
                total_lower += load_reactive_lower[load];
                total_upper += load_reactive_upper[load];
            }
            const double target = std::clamp(
                current - q_excess, total_lower, total_upper);
            if (allocate_total(
                    data_.buses[bus].loads,
                    load_reactive_lower, load_reactive_upper,
                    load_reactive_preferred, target, load_reactive)) {
                for (int load : data_.buses[bus].loads) {
                    if (std::abs(data_.loads[load].qd_nominal) > 1e-12) {
                        state.demand_factor[load] =
                            load_reactive[load] / data_.loads[load].qd_nominal;
                    }
                }
            }
        }

        balances = compute_balances();
        for (int bus = 0; bus < nb; ++bus) {
            const double p_excess = std::copysign(
                std::max(0.0, std::abs(balances.first[bus]) - 0.45),
                balances.first[bus]);
            if (std::abs(p_excess) <= 0.0 || active_at_bus[bus].empty()) {
                continue;
            }
            double current = 0.0;
            for (int gen : active_at_bus[bus]) {
                current += state.pg[gen];
            }
            auto proposed = state.pg;
            if (allocate_total(
                    active_at_bus[bus], p_lower, p_upper,
                    base_state_.pg, current + p_excess, proposed)) {
                state.pg = std::move(proposed);
            }
        }

        balances = compute_balances();
        std::vector<double> load_power(nd, 0.0);
        for (int i = 0; i < nd; ++i) {
            load_power[i] = data_.loads[i].pd_nominal * state.demand_factor[i];
        }
        for (int bus = 0; bus < nb; ++bus) {
            const double p_excess = std::copysign(
                std::max(0.0, std::abs(balances.first[bus]) - 0.45),
                balances.first[bus]);
            if (std::abs(p_excess) <= 0.0 || data_.buses[bus].loads.empty()) {
                continue;
            }
            double current = 0.0;
            double total_lower = 0.0;
            double total_upper = 0.0;
            for (int load : data_.buses[bus].loads) {
                current += load_power[load];
                total_lower += load_power_lower[load];
                total_upper += load_power_upper[load];
            }
            const double target = std::clamp(
                current - p_excess, total_lower, total_upper);
            allocate_total(
                data_.buses[bus].loads,
                load_power_lower, load_power_upper,
                load_power_preferred, target, load_power);
            for (int load : data_.buses[bus].loads) {
                if (std::abs(data_.loads[load].pd_nominal) > 1e-12) {
                    state.demand_factor[load] =
                        load_power[load] / data_.loads[load].pd_nominal;
                }
            }
        }

        balances = compute_balances();
        state.p_delta.assign(nb, 0.0);
        state.q_delta.assign(nb, 0.0);
        for (int bus = 0; bus < nb; ++bus) {
            state.p_delta[bus] = std::min(0.5, std::abs(balances.first[bus]));
            state.q_delta[bus] = std::min(0.5, std::abs(balances.second[bus]));
        }
        refresh_pwl_weights();
    };

    refresh_network_fields();

    ContingencyContext context;
    context.base_state = base_state_;
    context.outaged_generator = outaged_generator;
    context.outaged_branch = outaged_branch;
    output.solve.status = 0;
    output.solve.objective = 0.0;
    output.solve.iterations = output.newton_iterations;
    output.solve.state = std::move(state);
    output.converged = !approximate_candidate;
    output.validation = validate_state(
        data_, ModelMode::ContingencySoft,
        output.solve.state, commitment_, context);
    for (int scaling_pass = 0;
         scaling_pass < 6 && output.validation.max_residual > options_.validation_tolerance;
         ++scaling_pass) {
        if (output.validation.worst_category != "variable_bound" &&
            output.validation.worst_category != "flow_limit") {
            break;
        }
        double maximum_component_ratio = 1.0;
        for (int i = 0; i < static_cast<int>(data_.branches.size()); ++i) {
            if (i == outaged_branch || data_.branches[i].rate_a <= 1e-12) {
                continue;
            }
            maximum_component_ratio = std::max(maximum_component_ratio,
                std::max({std::abs(output.solve.state.pf[i]),
                          std::abs(output.solve.state.qf[i]),
                          std::abs(output.solve.state.pt[i]),
                          std::abs(output.solve.state.qt[i])})
                    / data_.branches[i].rate_a);
        }
        const double scale = std::min(0.9995,
            0.9999 / std::sqrt(maximum_component_ratio));
        bool changed = false;
        for (int bus = 0; bus < nb; ++bus) {
            const double reduced = std::max(
                data_.buses[bus].vmin, output.solve.state.vm[bus] * scale);
            changed = changed || reduced < output.solve.state.vm[bus] - 1e-12;
            output.solve.state.vm[bus] = reduced;
        }
        if (!changed) {
            break;
        }
        state = output.solve.state;
        refresh_network_fields();
        output.solve.state = state;
        output.validation = validate_state(
            data_, ModelMode::ContingencySoft,
            output.solve.state, commitment_, context);
    }
    output.feasible = output.validation.max_residual <= options_.validation_tolerance;
    if (!output.feasible) {
        output.failure_reason = "independent validation failed: "
            + output.validation.worst_category + " at "
            + output.validation.worst_identity;
    } else if (approximate_candidate) {
        output.failure_reason.clear();
    }
    output.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_start).count();
    output.solve.wall_seconds = output.wall_seconds;
    return output;
}

}  // namespace gravityx
