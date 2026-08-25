#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include "gravityx/active_feasibility_repair.hpp"
#include "gravityx/fast_power_flow.hpp"
#include "gravityx/state_io.hpp"

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

void ensure_shunt_control_state(const CaseData& data, AcState& state) {
    if (state.shunt_steps.size() != data.shunts.size()) {
        state.shunt_steps.clear();
        state.shunt_steps.reserve(data.shunts.size());
        for (const auto& shunt : data.shunts) {
            state.shunt_steps.push_back(shunt.steps);
        }
    }
    if (state.shunt_bs.size() != data.shunts.size()) {
        state.shunt_bs.clear();
        state.shunt_bs.reserve(data.shunts.size());
        for (const auto& shunt : data.shunts) {
            state.shunt_bs.push_back(shunt.bs);
        }
    }
}

YRows build_ybus(
    const CaseData& data,
    int outaged_branch,
    const AcState* state = nullptr) {
    YRows rows(data.buses.size());
    for (int i = 0; i < static_cast<int>(data.shunts.size()); ++i) {
        const auto& shunt = data.shunts[i];
        const double bs = state != nullptr
            ? effective_shunt_susceptance(data, *state, i)
            : shunt.bs;
        add_admittance(rows, shunt.bus, shunt.bus, {shunt.gs, bs});
    }
    for (int i = 0; i < static_cast<int>(data.branches.size()); ++i) {
        if (i == outaged_branch || data.branches[i].status == 0) {
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

NewtonResult run_distributed_active_newton(
    const CaseData& data,
    const YRows& ybus,
    const std::vector<bool>& slack,
    const std::vector<bool>& pq,
    const std::vector<int>& component_of,
    const std::vector<double>& active_slack_weights,
    const std::vector<double>& p_spec,
    const std::vector<double>& q_spec,
    int max_iterations,
    double tolerance,
    std::vector<double>& vm,
    std::vector<double>& va) {
    const int nb = static_cast<int>(data.buses.size());
    if (component_of.size() != static_cast<std::size_t>(nb) ||
        active_slack_weights.size() != static_cast<std::size_t>(nb)) {
        return {false, 0, "distributed-slack dimensions are inconsistent"};
    }

    int component_count = 0;
    for (int component : component_of) {
        component_count = std::max(component_count, component + 1);
    }
    if (component_count <= 0) {
        return {false, 0, "distributed-slack component set is empty"};
    }

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
    const int voltage_offset = angle_count;
    const int delta_offset = angle_count + voltage_count;
    const int dimension = delta_offset + component_count;
    if (dimension != nb + voltage_count) {
        return {false, 0, "distributed-slack Newton system is not square"};
    }

    std::vector<double> distributed_delta(component_count, 0.0);
    std::vector<double> distributed_weight_sum(component_count, 0.0);
    std::vector<double> p;
    std::vector<double> q;
    network_injections(ybus, vm, va, p, q);
    for (int i = 0; i < nb; ++i) {
        distributed_delta[component_of[i]] += p[i] - p_spec[i];
        distributed_weight_sum[component_of[i]] += active_slack_weights[i];
    }
    for (int c = 0; c < component_count; ++c) {
        if (std::abs(distributed_weight_sum[c]) <= 1e-12) {
            return {false, 0,
                "distributed-slack weights have zero component sum"};
        }
        distributed_delta[c] /= distributed_weight_sum[c];
    }
    const auto residual_norm = [&]() {
        double result = 0.0;
        for (int i = 0; i < nb; ++i) {
            const double distributed_target = p_spec[i]
                + active_slack_weights[i] * distributed_delta[component_of[i]];
            result = std::max(result, std::abs(distributed_target - p[i]));
            if (voltage_index[i] >= 0) {
                result = std::max(result, std::abs(q_spec[i] - q[i]));
            }
        }
        return result;
    };

    for (int iteration = 0; iteration <= max_iterations; ++iteration) {
        network_injections(ybus, vm, va, p, q);
        const double norm = residual_norm();
        if (norm <= tolerance) {
            return {true, iteration, {}};
        }
        if (iteration == max_iterations) {
            return {false, iteration,
                "distributed-slack Newton iteration limit; residual="
                    + std::to_string(norm)};
        }

        Eigen::VectorXd mismatch(dimension);
        mismatch.setZero();
        for (int i = 0; i < nb; ++i) {
            mismatch[i] = p_spec[i]
                + active_slack_weights[i] * distributed_delta[component_of[i]]
                - p[i];
            if (voltage_index[i] >= 0) {
                mismatch[nb + voltage_index[i]] = q_spec[i] - q[i];
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
            const double gii = diagonal == ybus[i].end()
                ? 0.0 : diagonal->second.real();
            const double bii = diagonal == ybus[i].end()
                ? 0.0 : diagonal->second.imag();

            const int p_row = i;
            if (angle_index[i] >= 0) {
                entries.emplace_back(
                    p_row, angle_index[i], -q[i] - bii * vi * vi);
            }
            if (voltage_index[i] >= 0) {
                entries.emplace_back(
                    p_row, voltage_offset + voltage_index[i],
                    p[i] / vi + gii * vi);
            }
            if (std::abs(active_slack_weights[i]) > 0.0) {
                entries.emplace_back(
                    p_row, delta_offset + component_of[i],
                    -active_slack_weights[i]);
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
                        p_row, angle_index[j],
                        vi * vm[j]
                            * (gij * std::sin(delta) - bij * std::cos(delta)));
                }
                if (voltage_index[j] >= 0) {
                    entries.emplace_back(
                        p_row, voltage_offset + voltage_index[j],
                        vi * (gij * std::cos(delta) + bij * std::sin(delta)));
                }
            }

            if (voltage_index[i] < 0) {
                continue;
            }
            const int q_row = nb + voltage_index[i];
            if (angle_index[i] >= 0) {
                entries.emplace_back(
                    q_row, angle_index[i], p[i] - gii * vi * vi);
            }
            entries.emplace_back(
                q_row, voltage_offset + voltage_index[i],
                q[i] / vi - bii * vi);
            for (const auto& [j, admittance] : ybus[i]) {
                if (j == i) {
                    continue;
                }
                const double delta = va[i] - va[j];
                const double gij = admittance.real();
                const double bij = admittance.imag();
                if (angle_index[j] >= 0) {
                    entries.emplace_back(
                        q_row, angle_index[j],
                        -vi * vm[j]
                            * (gij * std::cos(delta) + bij * std::sin(delta)));
                }
                if (voltage_index[j] >= 0) {
                    entries.emplace_back(
                        q_row, voltage_offset + voltage_index[j],
                        vi * (gij * std::sin(delta) - bij * std::cos(delta)));
                }
            }
        }

        SparseMatrix jacobian(dimension, dimension);
        jacobian.setFromTriplets(entries.begin(), entries.end());
        Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>> factorization;
        factorization.analyzePattern(jacobian);
        factorization.factorize(jacobian);
        if (factorization.info() != Eigen::Success) {
            return {false, iteration,
                "distributed-slack Jacobian factorization failed"};
        }
        const Eigen::VectorXd step = factorization.solve(mismatch);
        if (factorization.info() != Eigen::Success || !step.allFinite()) {
            return {false, iteration, "distributed-slack Newton solve failed"};
        }

        double scale = 1.0;
        double max_angle_step = 0.0;
        double max_voltage_step = 0.0;
        for (int i = 0; i < nb; ++i) {
            if (angle_index[i] >= 0) {
                max_angle_step = std::max(
                    max_angle_step, std::abs(step[angle_index[i]]));
            }
            if (voltage_index[i] >= 0) {
                max_voltage_step = std::max(
                    max_voltage_step,
                    std::abs(step[voltage_offset + voltage_index[i]]));
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
        const auto prior_delta = distributed_delta;
        const double prior_norm = norm;
        bool accepted = false;
        for (int line_search = 0; line_search < 12; ++line_search) {
            vm = prior_vm;
            va = prior_va;
            distributed_delta = prior_delta;
            for (int i = 0; i < nb; ++i) {
                if (angle_index[i] >= 0) {
                    va[i] += scale * step[angle_index[i]];
                }
                if (voltage_index[i] >= 0) {
                    vm[i] += scale * step[voltage_offset + voltage_index[i]];
                }
            }
            for (int c = 0; c < component_count; ++c) {
                distributed_delta[c] += scale * step[delta_offset + c];
            }
            network_injections(ybus, vm, va, p, q);
            const double trial_norm = residual_norm();
            if (std::isfinite(trial_norm) && trial_norm < prior_norm) {
                accepted = true;
                break;
            }
            scale *= 0.5;
        }
        if (!accepted) {
            vm = prior_vm;
            va = prior_va;
            return {false, iteration,
                "distributed-slack Newton line search stagnated"};
        }
    }
    return {false, max_iterations, "distributed-slack Newton failed"};
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
    bool contingency_mode,
    AcState& state) {
    const int nl = static_cast<int>(data.branches.size());
    state.pf.assign(nl, 0.0);
    state.qf.assign(nl, 0.0);
    state.pt.assign(nl, 0.0);
    state.qt.assign(nl, 0.0);
    state.sm_slack.assign(nl, 0.0);
    for (int i = 0; i < nl; ++i) {
        if (i == outaged_branch || data.branches[i].status == 0) {
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
        const double rating = contingency_mode ? branch.rate_c : branch.rate_a;
        if (rating > 1e-12) {
            const double from_required = branch.transformer
                ? from_magnitude / rating - 1.0
                : from_magnitude / rating - state.vm[f];
            const double to_required = branch.transformer
                ? to_magnitude / rating - 1.0
                : to_magnitude / rating - state.vm[t];
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
        if (i == outaged_branch || data.branches[i].status == 0) {
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

struct FastContingencyPowerFlow::FixedJacobianPredictorCache {
    struct LowRankUpdate {
        bool enabled{};
        Eigen::MatrixXd difference_rows;
        Eigen::MatrixXd correction;
    };

    bool valid{};
    double preparation_seconds{};
    std::vector<int> angle_index;
    std::vector<int> voltage_index;
    int angle_count{};
    int voltage_count{};
    Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>> factorization;
    bool active_valid{};
    SparseMatrix active_jacobian;
    Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>>
        active_factorization;
    bool reactive_valid{};
    SparseMatrix reactive_jacobian;
    Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>>
        reactive_factorization;
    LowRankUpdate active_outage_update;
    LowRankUpdate reactive_outage_update;

    FixedJacobianPredictorCache(
        const CaseData& data,
        const AcState& reference_state,
        const std::vector<int>& commitment,
        int outaged_branch = -1) {
        const auto started = std::chrono::steady_clock::now();
        const int nb = static_cast<int>(data.buses.size());
        std::vector<bool> slack(static_cast<std::size_t>(nb), false);
        for (const auto& component :
             connected_components(data, outaged_branch)) {
            int reference = -1;
            for (int bus : component) {
                bool has_committed_generator = false;
                for (int generator : data.buses[bus].generators) {
                    has_committed_generator = has_committed_generator ||
                        commitment[generator] == 1;
                }
                if (has_committed_generator && data.buses[bus].type == 3) {
                    reference = bus;
                    break;
                }
            }
            if (reference < 0) {
                for (int bus : component) {
                    for (int generator : data.buses[bus].generators) {
                        if (commitment[generator] == 1) {
                            reference = bus;
                            break;
                        }
                    }
                    if (reference >= 0) {
                        break;
                    }
                }
            }
            if (reference < 0) {
                preparation_seconds = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - started).count();
                return;
            }
            slack[reference] = true;
        }

        angle_index.assign(static_cast<std::size_t>(nb), -1);
        voltage_index.assign(static_cast<std::size_t>(nb), -1);
        for (int bus = 0; bus < nb; ++bus) {
            if (!slack[bus]) {
                angle_index[bus] = angle_count++;
            }
            // The reference angle removes one active-balance equation and
            // one angle variable per connected component.  Voltage magnitude
            // remains corrective at that bus, so its reactive-balance
            // equation and voltage variable must stay in the square system.
            voltage_index[bus] = voltage_count++;
        }
        const int dimension = angle_count + voltage_count;
        if (dimension <= 0) {
            preparation_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - started).count();
            return;
        }

        const auto ybus = build_ybus(
            data, outaged_branch, &reference_state);
        std::vector<double> p;
        std::vector<double> q;
        network_injections(
            ybus, reference_state.vm, reference_state.va, p, q);
        std::vector<Triplet> entries;
        entries.reserve(static_cast<std::size_t>(dimension) * 8);
        std::vector<Triplet> active_entries;
        active_entries.reserve(static_cast<std::size_t>(angle_count) * 4);
        std::vector<Triplet> reactive_entries;
        reactive_entries.reserve(static_cast<std::size_t>(voltage_count) * 4);
        for (int i = 0; i < nb; ++i) {
            const double vi = reference_state.vm[i];
            if (vi <= 1e-8) {
                preparation_seconds = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - started).count();
                return;
            }
            const auto diagonal = ybus[i].find(i);
            const double gii = diagonal == ybus[i].end()
                ? 0.0 : diagonal->second.real();
            const double bii = diagonal == ybus[i].end()
                ? 0.0 : diagonal->second.imag();
            if (angle_index[i] >= 0) {
                const int row = angle_index[i];
                const double active_angle_diagonal =
                    -q[i] - bii * vi * vi;
                entries.emplace_back(
                    row, angle_index[i], active_angle_diagonal);
                active_entries.emplace_back(
                    row, angle_index[i], active_angle_diagonal);
                entries.emplace_back(
                    row, angle_count + voltage_index[i],
                    p[i] / vi + gii * vi);
                for (const auto& [j, admittance] : ybus[i]) {
                    if (j == i) {
                        continue;
                    }
                    const double delta =
                        reference_state.va[i] - reference_state.va[j];
                    const double gij = admittance.real();
                    const double bij = admittance.imag();
                    if (angle_index[j] >= 0) {
                        const double active_angle_off_diagonal =
                            vi * reference_state.vm[j] *
                            (gij * std::sin(delta) -
                             bij * std::cos(delta));
                        entries.emplace_back(
                            row, angle_index[j],
                            active_angle_off_diagonal);
                        active_entries.emplace_back(
                            row, angle_index[j],
                            active_angle_off_diagonal);
                    }
                    if (voltage_index[j] >= 0) {
                        entries.emplace_back(
                            row, angle_count + voltage_index[j],
                            vi * (gij * std::cos(delta) +
                                  bij * std::sin(delta)));
                    }
                }
            }
            if (voltage_index[i] >= 0) {
                const int row = angle_count + voltage_index[i];
                if (angle_index[i] >= 0) {
                    entries.emplace_back(
                        row, angle_index[i], p[i] - gii * vi * vi);
                }
                const double reactive_voltage_diagonal =
                    q[i] / vi - bii * vi;
                entries.emplace_back(
                    row, angle_count + voltage_index[i],
                    reactive_voltage_diagonal);
                reactive_entries.emplace_back(
                    voltage_index[i], voltage_index[i],
                    reactive_voltage_diagonal);
                for (const auto& [j, admittance] : ybus[i]) {
                    if (j == i) {
                        continue;
                    }
                    const double delta =
                        reference_state.va[i] - reference_state.va[j];
                    const double gij = admittance.real();
                    const double bij = admittance.imag();
                    if (angle_index[j] >= 0) {
                        entries.emplace_back(
                            row, angle_index[j],
                            -vi * reference_state.vm[j] *
                                (gij * std::cos(delta) +
                                 bij * std::sin(delta)));
                    }
                    if (voltage_index[j] >= 0) {
                        const double reactive_voltage_off_diagonal =
                            vi * (gij * std::sin(delta) -
                                  bij * std::cos(delta));
                        entries.emplace_back(
                            row, angle_count + voltage_index[j],
                            reactive_voltage_off_diagonal);
                        reactive_entries.emplace_back(
                            voltage_index[i], voltage_index[j],
                            reactive_voltage_off_diagonal);
                    }
                }
            }
        }
        SparseMatrix jacobian(dimension, dimension);
        jacobian.setFromTriplets(entries.begin(), entries.end());
        active_jacobian.resize(angle_count, angle_count);
        active_jacobian.setFromTriplets(
            active_entries.begin(), active_entries.end());
        reactive_jacobian.resize(voltage_count, voltage_count);
        reactive_jacobian.setFromTriplets(
            reactive_entries.begin(), reactive_entries.end());
        factorization.analyzePattern(jacobian);
        factorization.factorize(jacobian);
        valid = factorization.info() == Eigen::Success;
        active_factorization.analyzePattern(active_jacobian);
        active_factorization.factorize(active_jacobian);
        active_valid = active_factorization.info() == Eigen::Success;
        reactive_factorization.analyzePattern(reactive_jacobian);
        reactive_factorization.factorize(reactive_jacobian);
        reactive_valid = reactive_factorization.info() == Eigen::Success;
        preparation_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();
    }

    template <typename Factorization>
    static bool prepare_low_rank_update(
        const SparseMatrix& difference,
        Factorization& base_factorization,
        LowRankUpdate& update) {
        update = {};
        const int dimension = difference.rows();
        if (dimension <= 0 || difference.cols() != dimension) {
            return false;
        }
        using RowSparseMatrix =
            Eigen::SparseMatrix<double, Eigen::RowMajor>;
        const RowSparseMatrix row_difference = difference;
        std::vector<int> changed_rows;
        for (int row = 0; row < row_difference.outerSize(); ++row) {
            bool changed = false;
            for (RowSparseMatrix::InnerIterator entry(row_difference, row);
                 entry; ++entry) {
                if (std::abs(entry.value()) > 1e-12) {
                    changed = true;
                    break;
                }
            }
            if (changed) {
                changed_rows.push_back(row);
            }
        }
        // Removing one AC branch changes only the P/Q equations at its two
        // terminal buses. More rows indicate an island/reference change, for
        // which the base factorization is not a valid Woodbury parent.
        if (changed_rows.size() > 4) {
            return false;
        }
        const int rank = static_cast<int>(changed_rows.size());
        update.enabled = true;
        if (rank == 0) {
            return true;
        }
        update.difference_rows = Eigen::MatrixXd::Zero(rank, dimension);
        Eigen::MatrixXd row_selector = Eigen::MatrixXd::Zero(dimension, rank);
        for (int position = 0; position < rank; ++position) {
            const int row = changed_rows[position];
            row_selector(row, position) = 1.0;
            for (RowSparseMatrix::InnerIterator entry(row_difference, row);
                 entry; ++entry) {
                update.difference_rows(position, entry.col()) = entry.value();
            }
        }
        const Eigen::MatrixXd inverse_rows =
            base_factorization.solve(row_selector);
        if (base_factorization.info() != Eigen::Success ||
            !inverse_rows.allFinite()) {
            update = {};
            return false;
        }
        const Eigen::MatrixXd small_system =
            Eigen::MatrixXd::Identity(rank, rank) +
            update.difference_rows * inverse_rows;
        Eigen::FullPivLU<Eigen::MatrixXd> small_factorization(small_system);
        small_factorization.setThreshold(1e-10);
        if (!small_factorization.isInvertible()) {
            update = {};
            return false;
        }
        update.correction = inverse_rows * small_factorization.inverse();
        if (!update.correction.allFinite()) {
            update = {};
            return false;
        }
        return true;
    }

    bool configure_branch_outage_update(
        const CaseData& data,
        const AcState& reference_state,
        int outaged_branch) {
        active_outage_update = {};
        reactive_outage_update = {};
        if (outaged_branch < 0 || !valid || !active_valid ||
            !reactive_valid) {
            return false;
        }
        if (outaged_branch >= static_cast<int>(data.branches.size())) {
            return false;
        }
        const auto& branch = data.branches[outaged_branch];
        if (branch.status == 0 || branch.from == branch.to ||
            std::abs(branch.tap) <= 1e-12) {
            return false;
        }
        const double denominator =
            branch.r * branch.r + branch.x * branch.x;
        if (denominator <= 1e-20) {
            return false;
        }
        const double g = branch.r / denominator;
        const double b = -branch.x / denominator;
        const double tm2 = branch.tap * branch.tap;
        const double tr = branch.tap * std::cos(branch.shift);
        const double ti = branch.tap * std::sin(branch.shift);
        const double from_g_self = branch.transformer
            ? g / tm2 + branch.g_fr : (g + branch.g_fr) / tm2;
        const double from_b_self = branch.transformer
            ? b / tm2 + branch.b_fr : (b + branch.b_fr) / tm2;
        const auto terminal_flows = [&](const std::array<double, 4>& point) {
            const double vf = point[0];
            const double af = point[1];
            const double vt = point[2];
            const double at = point[3];
            const double cross_cos_ft =
                vf * vt * std::cos(af - at);
            const double cross_sin_ft =
                vf * vt * std::sin(af - at);
            const double cross_cos_tf =
                vt * vf * std::cos(at - af);
            const double cross_sin_tf =
                vt * vf * std::sin(at - af);
            return std::array<double, 4>{
                from_g_self * vf * vf +
                    ((-g * tr + b * ti) / tm2) * cross_cos_ft +
                    ((-b * tr - g * ti) / tm2) * cross_sin_ft,
                -from_b_self * vf * vf -
                    ((-b * tr - g * ti) / tm2) * cross_cos_ft +
                    ((-g * tr + b * ti) / tm2) * cross_sin_ft,
                (g + branch.g_to) * vt * vt +
                    ((-g * tr - b * ti) / tm2) * cross_cos_tf +
                    ((-b * tr + g * ti) / tm2) * cross_sin_tf,
                -(b + branch.b_to) * vt * vt -
                    ((-b * tr + g * ti) / tm2) * cross_cos_tf +
                    ((-g * tr - b * ti) / tm2) * cross_sin_tf,
            };
        };
        const int from = branch.from;
        const int to = branch.to;
        const std::array<double, 4> reference{
            reference_state.vm[from], reference_state.va[from],
            reference_state.vm[to], reference_state.va[to]};
        Eigen::Matrix<double, 4, 4> derivative;
        constexpr double kDifferenceStep = 1e-6;
        for (int variable = 0; variable < 4; ++variable) {
            auto lower = reference;
            auto upper = reference;
            lower[variable] -= kDifferenceStep;
            upper[variable] += kDifferenceStep;
            const auto lower_flow = terminal_flows(lower);
            const auto upper_flow = terminal_flows(upper);
            for (int flow = 0; flow < 4; ++flow) {
                derivative(flow, variable) =
                    (upper_flow[flow] - lower_flow[flow]) /
                    (2.0 * kDifferenceStep);
            }
        }

        std::vector<Triplet> active_difference_entries;
        const std::array<int, 2> terminal_buses{from, to};
        const std::array<int, 2> active_flow_rows{0, 2};
        for (int terminal = 0; terminal < 2; ++terminal) {
            const int row = angle_index[terminal_buses[terminal]];
            if (row < 0) {
                continue;
            }
            for (int variable_terminal = 0; variable_terminal < 2;
                 ++variable_terminal) {
                const int column = angle_index[
                    terminal_buses[variable_terminal]];
                if (column >= 0) {
                    active_difference_entries.emplace_back(
                        row, column,
                        -derivative(
                            active_flow_rows[terminal],
                            1 + 2 * variable_terminal));
                }
            }
        }
        SparseMatrix active_difference(angle_count, angle_count);
        active_difference.setFromTriplets(
            active_difference_entries.begin(),
            active_difference_entries.end());

        std::vector<Triplet> reactive_difference_entries;
        const std::array<int, 2> reactive_flow_rows{1, 3};
        for (int terminal = 0; terminal < 2; ++terminal) {
            const int row = voltage_index[terminal_buses[terminal]];
            for (int variable_terminal = 0; variable_terminal < 2;
                 ++variable_terminal) {
                const int column = voltage_index[
                    terminal_buses[variable_terminal]];
                reactive_difference_entries.emplace_back(
                    row, column,
                    -derivative(
                        reactive_flow_rows[terminal],
                        2 * variable_terminal));
            }
        }
        SparseMatrix reactive_difference(voltage_count, voltage_count);
        reactive_difference.setFromTriplets(
            reactive_difference_entries.begin(),
            reactive_difference_entries.end());

        const bool active_ready = prepare_low_rank_update(
            active_difference, active_factorization, active_outage_update);
        const bool reactive_ready = prepare_low_rank_update(
            reactive_difference, reactive_factorization,
            reactive_outage_update);
        if (!active_ready || !reactive_ready) {
            active_outage_update = {};
            reactive_outage_update = {};
            return false;
        }
        return true;
    }

    template <typename Factorization>
    static Eigen::VectorXd solve_with_update(
        Factorization& base_factorization,
        const Eigen::VectorXd& right_hand_side,
        const LowRankUpdate& update) {
        Eigen::VectorXd solution =
            base_factorization.solve(right_hand_side);
        if (base_factorization.info() != Eigen::Success ||
            !solution.allFinite() || !update.enabled ||
            update.difference_rows.rows() == 0) {
            return solution;
        }
        solution -= update.correction *
            (update.difference_rows * solution);
        return solution;
    }

    bool apply_correction(
        const CaseData& data,
        const std::vector<double>& p_spec,
        const std::vector<double>& q_spec,
        const std::vector<double>& p_network,
        const std::vector<double>& q_network,
        std::vector<double>& vm,
        std::vector<double>& va,
        double damping = 1.0) {
        if (!valid) {
            return false;
        }
        Eigen::VectorXd mismatch(angle_count + voltage_count);
        mismatch.setZero();
        for (int bus = 0; bus < static_cast<int>(data.buses.size()); ++bus) {
            if (angle_index[bus] >= 0) {
                mismatch[angle_index[bus]] =
                    p_spec[bus] - p_network[bus];
            }
            if (voltage_index[bus] >= 0) {
                mismatch[angle_count + voltage_index[bus]] =
                    q_spec[bus] - q_network[bus];
            }
        }
        const Eigen::VectorXd step = factorization.solve(mismatch);
        if (factorization.info() != Eigen::Success || !step.allFinite()) {
            return false;
        }
        double scale = 1.0;
        double maximum_angle_step = 0.0;
        double maximum_voltage_step = 0.0;
        for (int bus = 0; bus < static_cast<int>(data.buses.size()); ++bus) {
            if (angle_index[bus] >= 0) {
                maximum_angle_step = std::max(
                    maximum_angle_step,
                    std::abs(step[angle_index[bus]]));
            }
            if (voltage_index[bus] >= 0) {
                maximum_voltage_step = std::max(
                    maximum_voltage_step,
                    std::abs(step[angle_count + voltage_index[bus]]));
            }
        }
        if (maximum_angle_step > 0.35) {
            scale = std::min(scale, 0.35 / maximum_angle_step);
        }
        if (maximum_voltage_step > 0.08) {
            scale = std::min(scale, 0.08 / maximum_voltage_step);
        }
        scale = std::clamp(scale * damping, 1e-4, 1.0);
        for (int bus = 0; bus < static_cast<int>(data.buses.size()); ++bus) {
            if (angle_index[bus] >= 0) {
                va[bus] += scale * step[angle_index[bus]];
            }
            if (voltage_index[bus] >= 0) {
                vm[bus] = std::clamp(
                    vm[bus] +
                        scale * step[angle_count + voltage_index[bus]],
                    data.buses[bus].vmin, data.buses[bus].vmax);
            }
        }
        return true;
    }

    bool apply_active_correction(
        const CaseData& data,
        const std::vector<double>& p_spec,
        const std::vector<double>& p_network,
        std::vector<double>& va,
        double damping = 1.0) {
        if (!active_valid) {
            return false;
        }
        Eigen::VectorXd mismatch(angle_count);
        mismatch.setZero();
        for (int bus = 0; bus < static_cast<int>(data.buses.size()); ++bus) {
            if (angle_index[bus] >= 0) {
                mismatch[angle_index[bus]] =
                    p_spec[bus] - p_network[bus];
            }
        }
        const Eigen::VectorXd step = solve_with_update(
            active_factorization, mismatch, active_outage_update);
        if (active_factorization.info() != Eigen::Success ||
            !step.allFinite()) {
            return false;
        }
        double scale = 1.0;
        double maximum_angle_step = 0.0;
        for (int bus = 0; bus < static_cast<int>(data.buses.size()); ++bus) {
            if (angle_index[bus] >= 0) {
                maximum_angle_step = std::max(
                    maximum_angle_step,
                    std::abs(step[angle_index[bus]]));
            }
        }
        if (maximum_angle_step > 0.35) {
            scale = 0.35 / maximum_angle_step;
        }
        scale = std::clamp(scale * damping, 1e-4, 1.0);
        for (int bus = 0; bus < static_cast<int>(data.buses.size()); ++bus) {
            if (angle_index[bus] >= 0) {
                va[bus] += scale * step[angle_index[bus]];
            }
        }
        return true;
    }

    bool apply_active_flow_redispatch(
        const CaseData& data,
        int branch_index,
        bool from_side,
        double desired_flow_change,
        const std::vector<double>& p_lower,
        const std::vector<double>& p_upper,
        const AcState& reference_state,
        std::vector<double>& pg,
        std::vector<double>& p_spec) {
        if (!active_valid || branch_index < 0 ||
            branch_index >= static_cast<int>(data.branches.size()) ||
            pg.size() != data.generators.size() ||
            p_lower.size() != data.generators.size() ||
            p_upper.size() != data.generators.size() ||
            p_spec.size() != data.buses.size()) {
            return false;
        }
        const auto& branch = data.branches[branch_index];
        const double denominator =
            branch.r * branch.r + branch.x * branch.x;
        if (denominator <= 1e-20 || branch.tap <= 1e-12) {
            return false;
        }
        const int from = branch.from;
        const int to = branch.to;
        const double g = branch.r / denominator;
        const double b = -branch.x / denominator;
        const double tm2 = branch.tap * branch.tap;
        const double tr = branch.tap * std::cos(branch.shift);
        const double ti = branch.tap * std::sin(branch.shift);
        const double angle =
            reference_state.va[from] - reference_state.va[to];
        const double voltage_product =
            reference_state.vm[from] * reference_state.vm[to];
        double derivative = 0.0;
        if (from_side) {
            const double cosine_coefficient =
                ((-g * tr + b * ti) / tm2) * voltage_product;
            const double sine_coefficient =
                ((-b * tr - g * ti) / tm2) * voltage_product;
            derivative =
                -cosine_coefficient * std::sin(angle) +
                sine_coefficient * std::cos(angle);
        } else {
            const double cosine_coefficient =
                ((-g * tr - b * ti) / tm2) * voltage_product;
            const double sine_coefficient =
                ((-b * tr + g * ti) / tm2) * voltage_product;
            derivative =
                -cosine_coefficient * std::sin(angle) -
                sine_coefficient * std::cos(angle);
        }
        if (std::abs(derivative) <= 1e-10) {
            return false;
        }
        Eigen::VectorXd flow_angle_gradient(angle_count);
        flow_angle_gradient.setZero();
        if (angle_index[from] >= 0) {
            flow_angle_gradient[angle_index[from]] += derivative;
        }
        if (angle_index[to] >= 0) {
            flow_angle_gradient[angle_index[to]] -= derivative;
        }
        const Eigen::VectorXd injection_sensitivity =
            active_factorization.transpose().solve(flow_angle_gradient);
        if (active_factorization.info() != Eigen::Success ||
            !injection_sensitivity.allFinite()) {
            return false;
        }
        const auto generator_sensitivity = [&](int generator) {
            const int bus = data.generators[generator].bus;
            return angle_index[bus] >= 0
                ? injection_sensitivity[angle_index[bus]] : 0.0;
        };
        int selected_up = -1;
        int selected_down = -1;
        double selected_coefficient = 0.0;
        for (int up = 0;
             up < static_cast<int>(data.generators.size()); ++up) {
            if (p_upper[up] - pg[up] <= 1e-8) {
                continue;
            }
            const double up_sensitivity = generator_sensitivity(up);
            for (int down = 0;
                 down < static_cast<int>(data.generators.size()); ++down) {
                if (down == up || pg[down] - p_lower[down] <= 1e-8) {
                    continue;
                }
                const double coefficient =
                    up_sensitivity - generator_sensitivity(down);
                if (std::abs(coefficient) <= 1e-10 ||
                    desired_flow_change / coefficient <= 0.0) {
                    continue;
                }
                if (selected_up < 0 ||
                    std::abs(coefficient) >
                        std::abs(selected_coefficient)) {
                    selected_up = up;
                    selected_down = down;
                    selected_coefficient = coefficient;
                }
            }
        }
        if (selected_up < 0 || selected_down < 0) {
            return false;
        }
        const double requested_change =
            desired_flow_change / selected_coefficient;
        const double redispatch = std::min({
            requested_change,
            p_upper[selected_up] - pg[selected_up],
            pg[selected_down] - p_lower[selected_down],
            1.0,
        });
        if (redispatch <= 1e-8) {
            return false;
        }
        pg[selected_up] += redispatch;
        pg[selected_down] -= redispatch;
        p_spec[data.generators[selected_up].bus] += redispatch;
        p_spec[data.generators[selected_down].bus] -= redispatch;
        return true;
    }

    bool apply_reactive_correction(
        const CaseData& data,
        const std::vector<double>& q_spec,
        const std::vector<double>& q_network,
        std::vector<double>& vm,
        double damping = 1.0) {
        if (!reactive_valid) {
            return false;
        }
        Eigen::VectorXd mismatch(voltage_count);
        mismatch.setZero();
        for (int bus = 0; bus < static_cast<int>(data.buses.size()); ++bus) {
            if (voltage_index[bus] >= 0) {
                mismatch[voltage_index[bus]] =
                    q_spec[bus] - q_network[bus];
            }
        }
        const Eigen::VectorXd step = solve_with_update(
            reactive_factorization, mismatch, reactive_outage_update);
        if (reactive_factorization.info() != Eigen::Success ||
            !step.allFinite()) {
            return false;
        }
        double scale = 1.0;
        double maximum_voltage_step = 0.0;
        for (int bus = 0; bus < static_cast<int>(data.buses.size()); ++bus) {
            if (voltage_index[bus] >= 0) {
                maximum_voltage_step = std::max(
                    maximum_voltage_step,
                    std::abs(step[voltage_index[bus]]));
            }
        }
        if (maximum_voltage_step > 0.08) {
            scale = 0.08 / maximum_voltage_step;
        }
        scale = std::clamp(scale * damping, 1e-4, 1.0);
        for (int bus = 0; bus < static_cast<int>(data.buses.size()); ++bus) {
            if (voltage_index[bus] >= 0) {
                vm[bus] = std::clamp(
                    vm[bus] + scale * step[voltage_index[bus]],
                    data.buses[bus].vmin, data.buses[bus].vmax);
            }
        }
        return true;
    }

    bool apply_local_reactive_correction(
        const CaseData& data,
        int target_bus,
        double target_network_q_change,
        std::vector<double>& vm,
        double damping = 1.0) {
        if (!reactive_valid || target_bus < 0 ||
            target_bus >= static_cast<int>(data.buses.size()) ||
            voltage_index[target_bus] < 0) {
            return false;
        }
        Eigen::VectorXd mismatch(voltage_count);
        mismatch.setZero();
        mismatch[voltage_index[target_bus]] = target_network_q_change;
        const Eigen::VectorXd step = solve_with_update(
            reactive_factorization, mismatch, reactive_outage_update);
        if (reactive_factorization.info() != Eigen::Success ||
            !step.allFinite()) {
            return false;
        }
        double scale = 1.0;
        double maximum_voltage_step = 0.0;
        for (int bus = 0; bus < static_cast<int>(data.buses.size()); ++bus) {
            if (voltage_index[bus] >= 0) {
                maximum_voltage_step = std::max(
                    maximum_voltage_step,
                    std::abs(step[voltage_index[bus]]));
            }
        }
        if (maximum_voltage_step > 0.08) {
            scale = 0.08 / maximum_voltage_step;
        }
        scale = std::clamp(scale * damping, 1e-4, 1.0);
        for (int bus = 0; bus < static_cast<int>(data.buses.size()); ++bus) {
            if (voltage_index[bus] >= 0) {
                vm[bus] = std::clamp(
                    vm[bus] + scale * step[voltage_index[bus]],
                    data.buses[bus].vmin, data.buses[bus].vmax);
            }
        }
        return true;
    }

    bool apply_reactive_band_correction(
        const CaseData& data,
        const std::vector<double>& q_balance,
        std::vector<double>& vm,
        double damping = 1.0) {
        if (!reactive_valid ||
            q_balance.size() != data.buses.size()) {
            return false;
        }
        Eigen::VectorXd mismatch(voltage_count);
        mismatch.setZero();
        for (int bus = 0; bus < static_cast<int>(data.buses.size()); ++bus) {
            if (voltage_index[bus] >= 0) {
                mismatch[voltage_index[bus]] =
                    std::clamp(q_balance[bus], -0.49, 0.49) -
                    q_balance[bus];
            }
        }
        if (mismatch.lpNorm<Eigen::Infinity>() <= 1e-12) {
            return false;
        }
        const Eigen::VectorXd step = solve_with_update(
            reactive_factorization, mismatch, reactive_outage_update);
        if (reactive_factorization.info() != Eigen::Success ||
            !step.allFinite()) {
            return false;
        }
        double scale = 1.0;
        double maximum_voltage_step = 0.0;
        for (int bus = 0; bus < static_cast<int>(data.buses.size()); ++bus) {
            if (voltage_index[bus] >= 0) {
                maximum_voltage_step = std::max(
                    maximum_voltage_step,
                    std::abs(step[voltage_index[bus]]));
            }
        }
        if (maximum_voltage_step > 0.08) {
            scale = 0.08 / maximum_voltage_step;
        }
        scale = std::clamp(scale * damping, 1e-4, 1.0);
        for (int bus = 0; bus < static_cast<int>(data.buses.size()); ++bus) {
            if (voltage_index[bus] >= 0) {
                vm[bus] = std::clamp(
                    vm[bus] + scale * step[voltage_index[bus]],
                    data.buses[bus].vmin, data.buses[bus].vmax);
            }
        }
        return true;
    }

    bool apply_local_reactive_least_squares(
        const CaseData& data,
        const std::vector<double>& q_balance,
        std::vector<double>& vm,
        double damping = 1.0) {
        if (!reactive_valid ||
            q_balance.size() != data.buses.size()) {
            return false;
        }
        std::vector<std::pair<double, int>> ranked_violations;
        for (int bus = 0; bus < static_cast<int>(data.buses.size()); ++bus) {
            const double excess = std::abs(q_balance[bus]) - 0.49;
            if (excess > 1e-8) {
                ranked_violations.emplace_back(excess, bus);
            }
        }
        if (ranked_violations.empty()) {
            return false;
        }
        std::sort(
            ranked_violations.begin(), ranked_violations.end(),
            [](const auto& left, const auto& right) {
                if (left.first != right.first) {
                    return left.first > right.first;
                }
                return left.second < right.second;
            });
        constexpr std::size_t kMaximumViolationRows = 32;
        if (ranked_violations.size() > kMaximumViolationRows) {
            ranked_violations.resize(kMaximumViolationRows);
        }

        std::vector<int> control_buses;
        for (const auto& [unused_excess, bus] : ranked_violations) {
            static_cast<void>(unused_excess);
            control_buses.push_back(bus);
            for (int branch : data.buses[bus].branches_from) {
                if (data.branches[branch].status != 0) {
                    control_buses.push_back(data.branches[branch].to);
                }
            }
            for (int branch : data.buses[bus].branches_to) {
                if (data.branches[branch].status != 0) {
                    control_buses.push_back(data.branches[branch].from);
                }
            }
        }
        std::sort(control_buses.begin(), control_buses.end());
        control_buses.erase(
            std::unique(control_buses.begin(), control_buses.end()),
            control_buses.end());
        constexpr std::size_t kMaximumControlBuses = 256;
        if (control_buses.size() > kMaximumControlBuses) {
            control_buses.resize(kMaximumControlBuses);
        }
        if (control_buses.empty()) {
            return false;
        }

        Eigen::MatrixXd sensitivity(
            static_cast<Eigen::Index>(ranked_violations.size()),
            static_cast<Eigen::Index>(control_buses.size()));
        Eigen::VectorXd target(
            static_cast<Eigen::Index>(ranked_violations.size()));
        for (std::size_t row = 0; row < ranked_violations.size(); ++row) {
            const int bus = ranked_violations[row].second;
            target[static_cast<Eigen::Index>(row)] =
                std::clamp(q_balance[bus], -0.49, 0.49) -
                q_balance[bus];
            for (std::size_t column = 0;
                 column < control_buses.size(); ++column) {
                sensitivity(
                    static_cast<Eigen::Index>(row),
                    static_cast<Eigen::Index>(column)) =
                    reactive_jacobian.coeff(
                        voltage_index[bus],
                        voltage_index[control_buses[column]]);
            }
        }
        Eigen::MatrixXd normal =
            sensitivity.transpose() * sensitivity;
        const double maximum_diagonal =
            std::max(1e-12, normal.diagonal().maxCoeff());
        normal.diagonal().array() += 1e-6 * maximum_diagonal;
        const Eigen::VectorXd right_hand_side =
            sensitivity.transpose() * target;
        const Eigen::VectorXd step =
            normal.ldlt().solve(right_hand_side);
        if (!step.allFinite()) {
            return false;
        }
        double scale = 1.0;
        const double maximum_step = step.lpNorm<Eigen::Infinity>();
        if (maximum_step > 0.02) {
            scale = 0.02 / maximum_step;
        }
        scale = std::clamp(scale * damping, 1e-4, 1.0);
        for (std::size_t column = 0;
             column < control_buses.size(); ++column) {
            const int bus = control_buses[column];
            vm[bus] = std::clamp(
                vm[bus] +
                    scale * step[static_cast<Eigen::Index>(column)],
                data.buses[bus].vmin, data.buses[bus].vmax);
        }
        return true;
    }
};

nlohmann::json FastPowerFlowResult::to_json() const {
    return {
        {"converged", converged},
        {"feasible", feasible},
        {"direct_candidate_attempted", direct_candidate_attempted},
        {"direct_candidate_selected", direct_candidate_selected},
        {"direct_candidate_validation", direct_candidate_validation.to_json()},
        {"fixed_jacobian_predictor_attempted",
         fixed_jacobian_predictor_attempted},
        {"fixed_jacobian_predictor_selected",
         fixed_jacobian_predictor_selected},
        {"fixed_jacobian_predictor_iterations",
         fixed_jacobian_predictor_iterations},
        {"fixed_jacobian_predictor_preparation_seconds",
         fixed_jacobian_predictor_preparation_seconds},
        {"fixed_jacobian_predictor_validation",
         fixed_jacobian_predictor_validation.to_json()},
        {"fixed_jacobian_predictor_trace",
         fixed_jacobian_predictor_trace},
        {"newton_candidate_selected", newton_candidate_selected},
        {"newton_candidate_validation", newton_candidate_validation.to_json()},
        {"active_only_newton_attempted", active_only_newton_attempted},
        {"active_only_newton_selected", active_only_newton_selected},
        {"active_only_newton_converged", active_only_newton_converged},
        {"active_only_newton_iterations", active_only_newton_iterations},
        {"active_only_newton_validation", active_only_newton_validation.to_json()},
        {"reactive_only_newton_attempted", reactive_only_newton_attempted},
        {"reactive_only_newton_selected", reactive_only_newton_selected},
        {"reactive_only_newton_converged", reactive_only_newton_converged},
        {"reactive_only_newton_iterations", reactive_only_newton_iterations},
        {"reactive_only_newton_validation", reactive_only_newton_validation.to_json()},
        {"distributed_balance_polish_attempted", distributed_balance_polish_attempted},
        {"distributed_balance_polish_selected", distributed_balance_polish_selected},
        {"distributed_balance_polish_iterations", distributed_balance_polish_iterations},
        {"distributed_balance_voltage_projections", distributed_balance_voltage_projections},
        {"distributed_balance_polish_failure_reason", distributed_balance_polish_failure_reason},
        {"distributed_balance_polish_validation", distributed_balance_polish_validation.to_json()},
        {"best_intermediate_candidate_selected",
         best_intermediate_candidate_selected},
        {"best_intermediate_candidate_source",
         best_intermediate_candidate_source},
        {"best_intermediate_candidate_validation",
         best_intermediate_candidate_validation.to_json()},
        {"newton_iterations", newton_iterations},
        {"initial_newton_residual", initial_newton_residual},
        {"active_redispatch_passes", active_redispatch_passes},
        {"reactive_limit_passes", reactive_limit_passes},
        {"wall_seconds", wall_seconds},
        {"failure_reason", failure_reason},
        {"validation", validation.to_json()},
    };
}

nlohmann::json ValidatedSourceBaseResult::to_json() const {
    return {
        {"feasible", feasible},
        {"wall_seconds", wall_seconds},
        {"solve", solve_result_to_json(solve, true)},
        {"validation", validation.to_json()},
    };
}

double rebuild_base_state_derived_fields(
    const CaseData& data,
    const std::vector<int>& commitment,
    AcState& state,
    double balance_slack_upper) {
    if (commitment.size() != data.generators.size() ||
        state.vm.size() != data.buses.size() ||
        state.va.size() != data.buses.size() ||
        state.pg.size() != data.generators.size() ||
        state.qg.size() != data.generators.size() ||
        state.demand_factor.size() != data.loads.size()) {
        throw std::runtime_error(
            "cannot rebuild a dimensionally invalid base state");
    }
    if (!std::isfinite(balance_slack_upper) ||
        balance_slack_upper < 0.0 || balance_slack_upper > 0.5) {
        throw std::runtime_error("invalid base-state balance slack upper bound");
    }
    ensure_shunt_control_state(data, state);

    double objective = 0.0;
    state.gen_lambda.clear();
    for (int i = 0; i < static_cast<int>(data.generators.size()); ++i) {
        const auto& gen = data.generators[i];
        if (commitment[i] != 0 && commitment[i] != 1) {
            throw std::runtime_error("base-state commitment is not binary");
        }
        if (commitment[i] == 0) {
            continue;
        }
        const double previous = gen.status_prev == 0 ? gen.pmin : gen.pg_prev;
        const double lower = std::max(
            gen.pmin, previous - data.delta_r * gen.prdmax);
        const double upper = std::min(
            gen.pmax, previous + data.delta_r * gen.prumax);
        const auto points = active_pwl_points(gen.cost, gen.ncost, lower, upper);
        const auto weights = interpolation_weights(points, state.pg[i]);
        state.gen_lambda.insert(
            state.gen_lambda.end(), weights.begin(), weights.end());
        for (int j = 0; j < static_cast<int>(points.size()); ++j) {
            objective -= data.delta * points[j].cost * weights[j];
        }
        objective -= data.delta * gen.oncost;
    }

    state.load_lambda.clear();
    for (int i = 0; i < static_cast<int>(data.loads.size()); ++i) {
        const auto& load = data.loads[i];
        const auto points = active_pwl_points(
            load.cost, load.ncost, load.pd_min, load.pd_max);
        const auto weights = interpolation_weights(
            points, load.pd_nominal * state.demand_factor[i]);
        state.load_lambda.insert(
            state.load_lambda.end(), weights.begin(), weights.end());
        for (int j = 0; j < static_cast<int>(points.size()); ++j) {
            objective += data.delta * points[j].cost * weights[j];
        }
    }

    compute_branch_flows(data, -1, false, state);
    const auto balance = nodal_balance_slack_seed(
        data, state, balance_slack_upper, 1e-7);
    state.p_delta = balance.active;
    state.q_delta = balance.reactive;
    for (double value : state.sm_slack) {
        objective -= data.delta * data.sm_cost_approx * value;
    }
    for (double value : state.p_delta) {
        objective -= data.delta * data.p_delta_cost_approx * value;
    }
    for (double value : state.q_delta) {
        objective -= data.delta * data.q_delta_cost_approx * value;
    }
    return objective;
}

double rebuild_contingency_state_derived_fields(
    const CaseData& data,
    const AcState& base_state,
    const std::vector<int>& commitment,
    const Contingency& contingency,
    AcState& state,
    double balance_slack_upper) {
    if (commitment.size() != data.generators.size() ||
        base_state.pg.size() != data.generators.size() ||
        base_state.demand_factor.size() != data.loads.size() ||
        state.vm.size() != data.buses.size() ||
        state.va.size() != data.buses.size() ||
        state.pg.size() != data.generators.size() ||
        state.qg.size() != data.generators.size() ||
        state.demand_factor.size() != data.loads.size()) {
        throw std::runtime_error(
            "cannot rebuild a dimensionally invalid contingency state");
    }
    if (!std::isfinite(balance_slack_upper) ||
        balance_slack_upper < 0.0 || balance_slack_upper > 0.5) {
        throw std::runtime_error(
            "invalid contingency-state balance slack upper bound");
    }
    ensure_shunt_control_state(data, state);
    const int outaged_generator =
        contingency.type == ContingencyType::Generator
        ? contingency.component : -1;
    const int outaged_branch =
        contingency.type == ContingencyType::Branch
        ? contingency.component : -1;

    double objective = 0.0;
    state.gen_lambda.clear();
    for (int i = 0; i < static_cast<int>(data.generators.size()); ++i) {
        const auto& gen = data.generators[i];
        const bool active = commitment[i] == 1 && i != outaged_generator;
        if (!active) {
            state.pg[i] = 0.0;
            state.qg[i] = 0.0;
            continue;
        }
        const double lower = std::max(
            gen.pmin,
            base_state.pg[i] - data.delta_r_ctg * gen.prdmaxctg);
        const double upper = std::min(
            gen.pmax,
            base_state.pg[i] + data.delta_r_ctg * gen.prumaxctg);
        const auto points = active_pwl_points(gen.cost, gen.ncost, lower, upper);
        const auto weights = interpolation_weights(points, state.pg[i]);
        state.gen_lambda.insert(
            state.gen_lambda.end(), weights.begin(), weights.end());
        for (int j = 0; j < static_cast<int>(points.size()); ++j) {
            objective -= data.delta_ctg * points[j].cost * weights[j];
        }
        objective -= data.delta_ctg * gen.oncost;
    }

    state.load_lambda.clear();
    for (int i = 0; i < static_cast<int>(data.loads.size()); ++i) {
        const auto& load = data.loads[i];
        const auto points = active_pwl_points(
            load.cost, load.ncost, load.pd_min, load.pd_max);
        const auto weights = interpolation_weights(
            points, load.pd_nominal * state.demand_factor[i]);
        state.load_lambda.insert(
            state.load_lambda.end(), weights.begin(), weights.end());
        for (int j = 0; j < static_cast<int>(points.size()); ++j) {
            objective += data.delta_ctg * points[j].cost * weights[j];
        }
    }

    compute_branch_flows(data, outaged_branch, true, state);
    const auto balance = nodal_balance_slack_seed(
        data, state, balance_slack_upper, 1e-7);
    state.p_delta = balance.active;
    state.q_delta = balance.reactive;
    for (double value : state.sm_slack) {
        objective -= data.delta_ctg * data.sm_cost_approx * value;
    }
    for (double value : state.p_delta) {
        objective -= data.delta_ctg * data.p_delta_cost_approx * value;
    }
    for (double value : state.q_delta) {
        objective -= data.delta_ctg * data.q_delta_cost_approx * value;
    }
    return objective;
}

ValidatedSourceBaseResult build_validated_source_base(
    const CaseData& data,
    std::vector<int> commitment,
    double validation_tolerance) {
    const auto wall_start = std::chrono::steady_clock::now();
    if (commitment.size() != data.generators.size()) {
        throw std::runtime_error("validated source base commitment has wrong length");
    }
    if (!std::isfinite(validation_tolerance) || validation_tolerance < 0.0) {
        throw std::runtime_error("invalid source-base validation tolerance");
    }

    AcState state;
    state.vm.resize(data.buses.size());
    state.va.resize(data.buses.size());
    for (int i = 0; i < static_cast<int>(data.buses.size()); ++i) {
        state.vm[i] = std::clamp(
            data.buses[i].vm_start, data.buses[i].vmin, data.buses[i].vmax);
        state.va[i] = data.buses[i].va_start;
    }
    for (const auto& component : connected_components(data, -1)) {
        int reference = -1;
        for (int bus : component) {
            if (data.buses[bus].type == 3) {
                reference = bus;
                break;
            }
        }
        if (reference >= 0) {
            const double offset = state.va[reference];
            for (int bus : component) {
                state.va[bus] -= offset;
            }
        }
    }

    state.pg.assign(data.generators.size(), 0.0);
    state.qg.assign(data.generators.size(), 0.0);
    double objective = 0.0;
    for (int i = 0; i < static_cast<int>(data.generators.size()); ++i) {
        const auto& gen = data.generators[i];
        if (commitment[i] != 0 && commitment[i] != 1) {
            throw std::runtime_error("source-base commitment is not binary");
        }
        if (commitment[i] == 0) {
            continue;
        }
        const double previous = gen.status_prev == 0 ? gen.pmin : gen.pg_prev;
        const double lower = std::max(
            gen.pmin, previous - data.delta_r * gen.prdmax);
        const double upper = std::min(
            gen.pmax, previous + data.delta_r * gen.prumax);
        if (lower > upper + 1e-12) {
            throw std::runtime_error(
                "empty source-base generator interval: " + gen.source_key);
        }
        state.pg[i] = std::clamp(gen.pg_start, lower, upper);
        state.qg[i] = std::clamp(gen.qg_start, gen.qmin, gen.qmax);
        const auto points = active_pwl_points(gen.cost, gen.ncost, lower, upper);
        const auto weights = interpolation_weights(points, state.pg[i]);
        state.gen_lambda.insert(
            state.gen_lambda.end(), weights.begin(), weights.end());
        for (int j = 0; j < static_cast<int>(points.size()); ++j) {
            objective -= data.delta * points[j].cost * weights[j];
        }
        objective -= data.delta * gen.oncost;
    }

    state.demand_factor.resize(data.loads.size());
    for (int i = 0; i < static_cast<int>(data.loads.size()); ++i) {
        const auto& load = data.loads[i];
        double lower = load.tmin;
        double upper = load.tmax;
        if (std::abs(load.pd_nominal) > 1e-12) {
            lower = std::max(
                lower,
                (load.pd_prev - load.prdmax * data.delta_r) / load.pd_nominal);
            upper = std::min(
                upper,
                (load.pd_prev + load.prumax * data.delta_r) / load.pd_nominal);
        }
        if (lower > upper + 1e-12) {
            throw std::runtime_error(
                "empty source-base load interval: " + load.source_key);
        }
        state.demand_factor[i] = std::clamp(load.z_start, lower, upper);
        const auto points = active_pwl_points(
            load.cost, load.ncost, load.pd_min, load.pd_max);
        const auto weights = interpolation_weights(
            points, load.pd_nominal * state.demand_factor[i]);
        state.load_lambda.insert(
            state.load_lambda.end(), weights.begin(), weights.end());
        for (int j = 0; j < static_cast<int>(points.size()); ++j) {
            objective += data.delta * points[j].cost * weights[j];
        }
    }

    compute_branch_flows(data, -1, false, state);
    const auto balance = nodal_balance_slack_seed(data, state, 0.5);
    state.p_delta = balance.active;
    state.q_delta = balance.reactive;
    for (double value : state.sm_slack) {
        objective -= data.delta * data.sm_cost_approx * value;
    }
    for (double value : state.p_delta) {
        objective -= data.delta * data.p_delta_cost_approx * value;
    }
    for (double value : state.q_delta) {
        objective -= data.delta * data.q_delta_cost_approx * value;
    }

    ValidatedSourceBaseResult output;
    output.solve.objective = objective;
    output.solve.iterations = 0;
    output.solve.state = std::move(state);
    output.validation = validate_state(
        data, ModelMode::BaseSoft, output.solve.state, commitment);
    output.feasible = std::isfinite(objective) &&
        std::isfinite(output.validation.max_residual) &&
        output.validation.max_residual <= validation_tolerance;
    output.solve.status = output.feasible ? 0 : 2;
    output.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_start).count();
    output.solve.wall_seconds = output.wall_seconds;
    return output;
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

FastContingencyPowerFlow::~FastContingencyPowerFlow() = default;

FastPowerFlowResult FastContingencyPowerFlow::solve(
    const Contingency& contingency) const {
    return solve_impl(&contingency);
}

FastPowerFlowResult FastContingencyPowerFlow::solve(
    const Contingency& contingency,
    const AcState& initial_state) const {
    return solve_impl(&contingency, &initial_state);
}

FastPowerFlowResult FastContingencyPowerFlow::screen_candidate(
    const Contingency& contingency,
    const AcState& candidate_state) const {
    return solve_impl(&contingency, &candidate_state, true);
}

FastPowerFlowResult FastContingencyPowerFlow::solve_base() const {
    return solve_impl(nullptr);
}

FastPowerFlowResult FastContingencyPowerFlow::solve_impl(
    const Contingency* contingency,
    const AcState* supplied_initial_state,
    bool supplied_candidate_direct_only) const {
    const auto wall_start = std::chrono::steady_clock::now();
    FastPowerFlowResult output;
    const int nb = static_cast<int>(data_.buses.size());
    const int ng = static_cast<int>(data_.generators.size());
    const bool base_mode = contingency == nullptr;
    const auto branch_rating = [&](int index) {
        return base_mode
            ? data_.branches[index].rate_a
            : data_.branches[index].rate_c;
    };
    const int outaged_generator = !base_mode &&
        contingency->type == ContingencyType::Generator
        ? contingency->component : -1;
    const int outaged_branch = !base_mode &&
        contingency->type == ContingencyType::Branch
        ? contingency->component : -1;
    const AcState& initial_state = supplied_initial_state
        ? *supplied_initial_state : base_state_;
    if (initial_state.vm.size() != data_.buses.size() ||
        initial_state.va.size() != data_.buses.size() ||
        initial_state.pg.size() != data_.generators.size() ||
        initial_state.qg.size() != data_.generators.size() ||
        initial_state.demand_factor.size() != data_.loads.size()) {
        output.failure_reason = "fast power flow initial state has wrong dimensions";
        output.wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        return output;
    }

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
        if (base_mode) {
            const double previous = gen.status_prev == 0 ? gen.pmin : gen.pg_prev;
            p_lower[i] = std::max(
                gen.pmin, previous - data_.delta_r * gen.prdmax);
            p_upper[i] = std::min(
                gen.pmax, previous + data_.delta_r * gen.prumax);
        } else {
            p_lower[i] = std::max(
                gen.pmin,
                base_state_.pg[i] - data_.delta_r_ctg * gen.prdmaxctg);
            p_upper[i] = std::min(
                gen.pmax,
                base_state_.pg[i] + data_.delta_r_ctg * gen.prumaxctg);
        }
        q_lower[i] = gen.qmin;
        q_upper[i] = gen.qmax;
        if (p_lower[i] > p_upper[i] + 1e-12 || q_lower[i] > q_upper[i] + 1e-12) {
            output.failure_reason = "empty generator contingency interval";
            output.wall_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - wall_start).count();
            return output;
        }
        pg[i] = std::clamp(initial_state.pg[i], p_lower[i], p_upper[i]);
        qg[i] = std::clamp(initial_state.qg[i], q_lower[i], q_upper[i]);
    }

    std::optional<ContingencyContext> direct_context;
    if (!base_mode) {
        direct_context = ContingencyContext{};
        direct_context->base_state = base_state_;
        direct_context->outaged_generator = outaged_generator;
        direct_context->outaged_branch = outaged_branch;
    }
    AcState direct_state = initial_state;
    direct_state.pg = pg;
    direct_state.qg = qg;
    compute_branch_flows(data_, outaged_branch, !base_mode, direct_state);
    const auto direct_balance = nodal_balance_slack_seed(
        data_, direct_state, 0.5, 1e-7);
    direct_state.p_delta = direct_balance.active;
    direct_state.q_delta = direct_balance.reactive;
    direct_state.gen_lambda.clear();
    for (int i = 0; i < ng; ++i) {
        if (!active[i]) {
            continue;
        }
        const auto points = active_pwl_points(
            data_.generators[i].cost, data_.generators[i].ncost,
            p_lower[i], p_upper[i]);
        append_weights(direct_state.gen_lambda, points, direct_state.pg[i]);
    }
    direct_state.load_lambda.clear();
    for (int i = 0; i < static_cast<int>(data_.loads.size()); ++i) {
        const auto& load = data_.loads[i];
        const auto points = active_pwl_points(
            load.cost, load.ncost, load.pd_min, load.pd_max);
        append_weights(
            direct_state.load_lambda, points,
            load.pd_nominal * direct_state.demand_factor[i]);
    }
    output.direct_candidate_attempted = true;
    output.direct_candidate_validation = validate_state(
        data_,
        base_mode ? ModelMode::BaseSoft : ModelMode::ContingencySoft,
        direct_state, commitment_, direct_context);
    AcState best_intermediate_state = direct_state;
    ValidationReport best_intermediate_validation =
        output.direct_candidate_validation;
    std::string best_intermediate_source = "direct_candidate";
    const auto retain_best_intermediate = [&best_intermediate_state,
                                           &best_intermediate_validation,
                                           &best_intermediate_source](
                                              const AcState& candidate,
                                              const ValidationReport& validation,
                                              const std::string& source) {
        if (validation.max_residual + 1e-12 <
            best_intermediate_validation.max_residual) {
            best_intermediate_state = candidate;
            best_intermediate_validation = validation;
            best_intermediate_source = source;
        }
    };
    if (output.direct_candidate_validation.max_residual <=
        options_.validation_tolerance) {
        output.converged = true;
        output.feasible = true;
        output.direct_candidate_selected = true;
        output.solve.status = 0;
        output.solve.objective = base_mode
            ? rebuild_base_state_derived_fields(
                data_, commitment_, direct_state, 0.5)
            : rebuild_contingency_state_derived_fields(
                data_, base_state_, commitment_, *contingency,
                direct_state);
        output.solve.iterations = 0;
        output.solve.state = std::move(direct_state);
        output.validation = output.direct_candidate_validation;
        output.wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        output.solve.wall_seconds = output.wall_seconds;
        return output;
    }
    if (supplied_candidate_direct_only) {
        output.solve.status = 2;
        output.solve.iterations = 0;
        output.solve.state = direct_state;
        output.solve.objective = rebuild_contingency_state_derived_fields(
            data_, base_state_, commitment_, *contingency,
            output.solve.state);
        output.validation = validate_state(
            data_, ModelMode::ContingencySoft, output.solve.state,
            commitment_, direct_context);
        output.failure_reason =
            "supplied corrective candidate failed independent validation: "
            + output.validation.worst_category + " at "
            + output.validation.worst_identity;
        output.wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        output.solve.wall_seconds = output.wall_seconds;
        return output;
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

    std::vector<double> predictor_demand_factor =
        initial_state.demand_factor;
    std::vector<double> predictor_load_power_lower(data_.loads.size(), 0.0);
    std::vector<double> predictor_load_power_upper(data_.loads.size(), 0.0);
    std::vector<double> predictor_load_power_preferred(data_.loads.size(), 0.0);
    std::vector<double> predictor_load_power(data_.loads.size(), 0.0);
    std::vector<double> predictor_load_reactive_lower(
        data_.loads.size(), 0.0);
    std::vector<double> predictor_load_reactive_upper(
        data_.loads.size(), 0.0);
    for (int load_index = 0;
         load_index < static_cast<int>(data_.loads.size()); ++load_index) {
        const auto& load = data_.loads[load_index];
        predictor_load_power_preferred[load_index] =
            load.pd_nominal * initial_state.demand_factor[load_index];
        predictor_load_power[load_index] =
            predictor_load_power_preferred[load_index];
        if (std::abs(load.pd_nominal) <= 1e-12) {
            predictor_load_power_lower[load_index] =
                predictor_load_power_preferred[load_index];
            predictor_load_power_upper[load_index] =
                predictor_load_power_preferred[load_index];
            const double reactive = load.qd_nominal *
                initial_state.demand_factor[load_index];
            predictor_load_reactive_lower[load_index] = reactive;
            predictor_load_reactive_upper[load_index] = reactive;
            continue;
        }
        const double prior = predictor_load_power_preferred[load_index];
        const double factor_lower = std::max(
            load.tmin,
            (prior - data_.delta_r_ctg * load.prdmaxctg) /
                load.pd_nominal);
        const double factor_upper = std::min(
            load.tmax,
            (prior + data_.delta_r_ctg * load.prumaxctg) /
                load.pd_nominal);
        const double lower_power = load.pd_nominal * factor_lower;
        const double upper_power = load.pd_nominal * factor_upper;
        predictor_load_power_lower[load_index] =
            std::min(lower_power, upper_power);
        predictor_load_power_upper[load_index] =
            std::max(lower_power, upper_power);
        predictor_load_reactive_lower[load_index] = std::min(
            load.qd_nominal * factor_lower,
            load.qd_nominal * factor_upper);
        predictor_load_reactive_upper[load_index] = std::max(
            load.qd_nominal * factor_lower,
            load.qd_nominal * factor_upper);
    }

    for (const auto& component : components) {
        std::vector<int> participants;
        std::vector<int> participating_loads;
        double active_target = 0.0;
        double active_lower = 0.0;
        double active_upper = 0.0;
        double reactive_target = 0.0;
        double reactive_lower = 0.0;
        double reactive_upper = 0.0;
        double prior_load_total = 0.0;
        double load_lower_total = 0.0;
        double load_upper_total = 0.0;
        for (int bus : component) {
            for (int gen : data_.buses[bus].generators) {
                if (commitment_[gen] == 1) {
                    active_target += initial_state.pg[gen];
                    reactive_target += initial_state.qg[gen];
                }
                if (active[gen]) {
                    participants.push_back(gen);
                    active_lower += p_lower[gen];
                    active_upper += p_upper[gen];
                    reactive_lower += q_lower[gen];
                    reactive_upper += q_upper[gen];
                }
            }
            for (int load_index : data_.buses[bus].loads) {
                prior_load_total +=
                    predictor_load_power_preferred[load_index];
                if (std::abs(data_.loads[load_index].pd_nominal) > 1e-12) {
                    participating_loads.push_back(load_index);
                    load_lower_total +=
                        predictor_load_power_lower[load_index];
                    load_upper_total +=
                        predictor_load_power_upper[load_index];
                }
            }
        }
        if (!participants.empty()) {
            const double prior_generation_total = active_target;
            active_target = std::clamp(
                active_target, active_lower, active_upper);
            allocate_total(
                participants, p_lower, p_upper,
                initial_state.pg, active_target, pg);
            reactive_target = std::clamp(
                reactive_target, reactive_lower, reactive_upper);
            allocate_total(
                participants, q_lower, q_upper,
                initial_state.qg, reactive_target, qg);
            double adjusted_generation_total = 0.0;
            for (int generator : participants) {
                adjusted_generation_total += pg[generator];
            }
            if (!participating_loads.empty()) {
                const double desired_load_total = std::clamp(
                    prior_load_total + adjusted_generation_total -
                        prior_generation_total,
                    load_lower_total, load_upper_total);
                allocate_total(
                    participating_loads,
                    predictor_load_power_lower,
                    predictor_load_power_upper,
                    predictor_load_power_preferred, desired_load_total,
                    predictor_load_power);
            }
        }
    }
    for (int load_index = 0;
         load_index < static_cast<int>(data_.loads.size()); ++load_index) {
        if (std::abs(data_.loads[load_index].pd_nominal) > 1e-12) {
            predictor_demand_factor[load_index] =
                predictor_load_power[load_index] /
                    data_.loads[load_index].pd_nominal;
        }
    }

    if (!base_mode && nb >= 16000) {
        output.fixed_jacobian_predictor_attempted = true;
        if (!predictor_cache_) {
            predictor_cache_ = std::make_unique<FixedJacobianPredictorCache>(
                data_, base_state_, commitment_);
            output.fixed_jacobian_predictor_preparation_seconds =
                predictor_cache_->preparation_seconds;
        }
        bool branch_outage_low_rank_update = false;
        if (predictor_cache_ && outaged_branch < 0) {
            predictor_cache_->configure_branch_outage_update(
                data_, base_state_, -1);
        }
        if (predictor_cache_ && predictor_cache_->valid) {
            AcState predictor_state = initial_state;
            predictor_state.pg = pg;
            predictor_state.qg = qg;
            predictor_state.demand_factor = predictor_demand_factor;
            ValidationReport predictor_validation =
                output.direct_candidate_validation;
            std::unique_ptr<FixedJacobianPredictorCache>
                contingency_predictor_cache;
            int active_feasibility_repair_attempts = 0;
            // Difficult generator outages can enter a late, monotonically
            // shrinking cycle between active- and reactive-flow limits after
            // the global feasibility repair.  Keep the inexpensive local
            // corrections alive long enough to cross the exact validator's
            // tolerance instead of stopping at an arbitrary 48-step edge.
            constexpr int kMaximumFixedJacobianIterations = 96;
            for (int predictor_iteration = 0;
                 predictor_iteration <= kMaximumFixedJacobianIterations;
                 ++predictor_iteration) {
                rebuild_contingency_state_derived_fields(
                    data_, base_state_, commitment_, *contingency,
                    predictor_state);

                std::vector<double> p_network(
                    static_cast<std::size_t>(nb), 0.0);
                std::vector<double> q_network(
                    static_cast<std::size_t>(nb), 0.0);
                for (int branch_index = 0;
                     branch_index < static_cast<int>(data_.branches.size());
                     ++branch_index) {
                    if (branch_index == outaged_branch ||
                        data_.branches[branch_index].status == 0) {
                        continue;
                    }
                    const auto& branch = data_.branches[branch_index];
                    p_network[branch.from] += predictor_state.pf[branch_index];
                    q_network[branch.from] += predictor_state.qf[branch_index];
                    p_network[branch.to] += predictor_state.pt[branch_index];
                    q_network[branch.to] += predictor_state.qt[branch_index];
                }
                for (int shunt_index = 0;
                     shunt_index < static_cast<int>(data_.shunts.size());
                     ++shunt_index) {
                    const auto& shunt = data_.shunts[shunt_index];
                    const double vm2 =
                        predictor_state.vm[shunt.bus] *
                        predictor_state.vm[shunt.bus];
                    p_network[shunt.bus] += shunt.gs * vm2;
                    q_network[shunt.bus] -= effective_shunt_susceptance(
                        data_, predictor_state, shunt_index) * vm2;
                }

                // Keep the component-wide active redispatch fixed while the
                // chord step routes it through the network.  Re-projecting P
                // independently at every generator bus changes the component
                // total and leaves that drift at the omitted reference-bus
                // equation.  Reactive generation can still be projected
                // locally because each Q equation remains in the predictor.
                for (int bus = 0; bus < nb; ++bus) {
                    if (active_at_bus[bus].empty()) {
                        continue;
                    }
                    double generated_q = 0.0;
                    double load_at_bus_q = 0.0;
                    for (int generator : active_at_bus[bus]) {
                        generated_q += predictor_state.qg[generator];
                    }
                    for (int load : data_.buses[bus].loads) {
                        load_at_bus_q += data_.loads[load].qd_nominal *
                            predictor_state.demand_factor[load];
                    }
                    const double q_balance =
                        q_network[bus] - generated_q + load_at_bus_q;
                    auto proposed_qg = predictor_state.qg;
                    if (allocate_total(
                            active_at_bus[bus], q_lower, q_upper,
                            initial_state.qg,
                            generated_q + q_balance -
                                std::clamp(q_balance, -0.49, 0.49),
                            proposed_qg)) {
                        predictor_state.qg = std::move(proposed_qg);
                    }
                }

                // A PQ load is also a source-authorized corrective control.
                // If local Q generation cannot put a bus inside the soft
                // balance band, move only the loads at that bus and respect
                // the exact contingency ramp and tmin/tmax bounds computed
                // above. This is especially important for voltage-limited
                // load buses where a voltage correction has no remaining
                // feasible direction.
                std::vector<double> predictor_load_reactive(
                    data_.loads.size(), 0.0);
                for (int load_index = 0;
                     load_index < static_cast<int>(data_.loads.size());
                     ++load_index) {
                    predictor_load_reactive[load_index] =
                        data_.loads[load_index].qd_nominal *
                        predictor_state.demand_factor[load_index];
                }
                const auto predictor_load_reactive_preferred =
                    predictor_load_reactive;
                for (int bus = 0; bus < nb; ++bus) {
                    if (data_.buses[bus].loads.empty()) {
                        continue;
                    }
                    double generated_q = 0.0;
                    double loaded_q = 0.0;
                    double total_lower = 0.0;
                    double total_upper = 0.0;
                    for (int generator : active_at_bus[bus]) {
                        generated_q += predictor_state.qg[generator];
                    }
                    for (int load : data_.buses[bus].loads) {
                        loaded_q += predictor_load_reactive[load];
                        total_lower +=
                            predictor_load_reactive_lower[load];
                        total_upper +=
                            predictor_load_reactive_upper[load];
                    }
                    const double q_balance =
                        q_network[bus] - generated_q + loaded_q;
                    const double q_excess = std::copysign(
                        std::max(0.0, std::abs(q_balance) - 0.49),
                        q_balance);
                    if (std::abs(q_excess) <= 1e-12) {
                        continue;
                    }
                    const double target = std::clamp(
                        loaded_q - q_excess, total_lower, total_upper);
                    if (!allocate_total(
                            data_.buses[bus].loads,
                            predictor_load_reactive_lower,
                            predictor_load_reactive_upper,
                            predictor_load_reactive_preferred,
                            target, predictor_load_reactive)) {
                        continue;
                    }
                    for (int load : data_.buses[bus].loads) {
                        if (std::abs(data_.loads[load].qd_nominal) > 1e-12) {
                            predictor_state.demand_factor[load] =
                                predictor_load_reactive[load] /
                                data_.loads[load].qd_nominal;
                        }
                    }
                }

                // Project an excessive local active-power mismatch onto the
                // source-authorized contingency controls before asking the
                // angle step to route power through the network.  Generation
                // is used first within its exact contingency PMIN/PMAX and
                // ramp bounds; any remainder is assigned to flexible load
                // within its exact tmin/tmax and corrective-ramp bounds.  The
                // 0.49 target stays strictly inside the source 0.5 balance
                // slack limit and leaves a small numerical margin.
                std::vector<double> predictor_load_active(
                    data_.loads.size(), 0.0);
                for (int load_index = 0;
                     load_index < static_cast<int>(data_.loads.size());
                     ++load_index) {
                    predictor_load_active[load_index] =
                        data_.loads[load_index].pd_nominal *
                        predictor_state.demand_factor[load_index];
                }
                for (int bus = 0; bus < nb; ++bus) {
                    double generated_p = 0.0;
                    double loaded_p = 0.0;
                    for (int generator : active_at_bus[bus]) {
                        generated_p += predictor_state.pg[generator];
                    }
                    for (int load : data_.buses[bus].loads) {
                        loaded_p += predictor_load_active[load];
                    }

                    double p_balance =
                        p_network[bus] - generated_p + loaded_p;
                    double p_excess = std::copysign(
                        std::max(0.0, std::abs(p_balance) - 0.49),
                        p_balance);
                    if (std::abs(p_excess) <= 1e-12) {
                        continue;
                    }

                    if (!active_at_bus[bus].empty()) {
                        double total_lower = 0.0;
                        double total_upper = 0.0;
                        for (int generator : active_at_bus[bus]) {
                            total_lower += p_lower[generator];
                            total_upper += p_upper[generator];
                        }
                        const double target_generation = std::clamp(
                            generated_p + p_excess,
                            total_lower, total_upper);
                        const auto preferred_pg = predictor_state.pg;
                        auto proposed_pg = predictor_state.pg;
                        if (allocate_total(
                                active_at_bus[bus], p_lower, p_upper,
                                preferred_pg, target_generation,
                                proposed_pg)) {
                            predictor_state.pg = std::move(proposed_pg);
                            generated_p = 0.0;
                            for (int generator : active_at_bus[bus]) {
                                generated_p += predictor_state.pg[generator];
                            }
                            p_balance =
                                p_network[bus] - generated_p + loaded_p;
                            p_excess = std::copysign(
                                std::max(
                                    0.0, std::abs(p_balance) - 0.49),
                                p_balance);
                        }
                    }

                    if (std::abs(p_excess) <= 1e-12 ||
                        data_.buses[bus].loads.empty()) {
                        continue;
                    }
                    double total_lower = 0.0;
                    double total_upper = 0.0;
                    for (int load : data_.buses[bus].loads) {
                        total_lower += predictor_load_power_lower[load];
                        total_upper += predictor_load_power_upper[load];
                    }
                    const double target_load = std::clamp(
                        loaded_p - p_excess, total_lower, total_upper);
                    const auto preferred_load = predictor_load_active;
                    if (!allocate_total(
                            data_.buses[bus].loads,
                            predictor_load_power_lower,
                            predictor_load_power_upper,
                            preferred_load, target_load,
                            predictor_load_active)) {
                        continue;
                    }
                    for (int load : data_.buses[bus].loads) {
                        if (std::abs(data_.loads[load].pd_nominal) > 1e-12) {
                            predictor_state.demand_factor[load] =
                                predictor_load_active[load] /
                                data_.loads[load].pd_nominal;
                        }
                    }
                }

                const double predictor_objective =
                    rebuild_contingency_state_derived_fields(
                        data_, base_state_, commitment_, *contingency,
                        predictor_state);
                predictor_validation = validate_state(
                    data_, ModelMode::ContingencySoft,
                    predictor_state, commitment_, direct_context);
                output.fixed_jacobian_predictor_iterations =
                    predictor_iteration;
                output.fixed_jacobian_predictor_validation =
                    predictor_validation;
                double predictor_generation_total = 0.0;
                for (double value : predictor_state.pg) {
                    predictor_generation_total += value;
                }
                double predictor_load_total = 0.0;
                for (int load_index = 0;
                     load_index < static_cast<int>(data_.loads.size());
                     ++load_index) {
                    predictor_load_total +=
                        data_.loads[load_index].pd_nominal *
                        predictor_state.demand_factor[load_index];
                }
                output.fixed_jacobian_predictor_trace.push_back({
                    {"iteration", predictor_iteration},
                    {"branch_outage_low_rank_update",
                     branch_outage_low_rank_update},
                    {"generation_total", predictor_generation_total},
                    {"load_total", predictor_load_total},
                    {"validation", predictor_validation.to_json()},
                });
                retain_best_intermediate(
                    predictor_state, predictor_validation,
                    "fixed_base_jacobian_predictor");
                if (predictor_validation.max_residual <=
                    options_.validation_tolerance) {
                    output.converged = true;
                    output.feasible = true;
                    output.fixed_jacobian_predictor_selected = true;
                    output.solve.status = 0;
                    output.solve.objective = predictor_objective;
                    output.solve.iterations = predictor_iteration;
                    output.solve.state = std::move(predictor_state);
                    output.validation = predictor_validation;
                    output.wall_seconds = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - wall_start).count();
                    output.solve.wall_seconds = output.wall_seconds;
                    return output;
                }
                if (predictor_iteration ==
                    kMaximumFixedJacobianIterations) {
                    break;
                }
                if (predictor_iteration == 0 && outaged_branch >= 0) {
                    branch_outage_low_rank_update =
                        predictor_cache_->configure_branch_outage_update(
                            data_, base_state_, outaged_branch);
                    output.fixed_jacobian_predictor_trace.back()[
                        "branch_outage_low_rank_prepared"] =
                        branch_outage_low_rank_update;
                }

                std::vector<double> p_spec(
                    static_cast<std::size_t>(nb), 0.0);
                std::vector<double> q_spec(
                    static_cast<std::size_t>(nb), 0.0);
                std::vector<double> generated_p_by_bus(
                    static_cast<std::size_t>(nb), 0.0);
                std::vector<double> generated_q_by_bus(
                    static_cast<std::size_t>(nb), 0.0);
                std::vector<double> load_p_by_bus(
                    static_cast<std::size_t>(nb), 0.0);
                std::vector<double> load_q_by_bus(
                    static_cast<std::size_t>(nb), 0.0);
                std::vector<double> p_balance_by_bus(
                    static_cast<std::size_t>(nb), 0.0);
                std::vector<double> q_balance_by_bus(
                    static_cast<std::size_t>(nb), 0.0);
                std::vector<double> p_slack_target(
                    static_cast<std::size_t>(nb), 0.0);
                std::vector<double> q_slack_target(
                    static_cast<std::size_t>(nb), 0.0);
                for (int bus = 0; bus < nb; ++bus) {
                    for (int generator : active_at_bus[bus]) {
                        generated_p_by_bus[bus] +=
                            predictor_state.pg[generator];
                        generated_q_by_bus[bus] +=
                            predictor_state.qg[generator];
                    }
                    for (int load : data_.buses[bus].loads) {
                        load_p_by_bus[bus] +=
                            data_.loads[load].pd_nominal *
                            predictor_state.demand_factor[load];
                        load_q_by_bus[bus] +=
                            data_.loads[load].qd_nominal *
                            predictor_state.demand_factor[load];
                    }
                    p_balance_by_bus[bus] = p_network[bus] -
                        generated_p_by_bus[bus] + load_p_by_bus[bus];
                    q_balance_by_bus[bus] = q_network[bus] -
                        generated_q_by_bus[bus] + load_q_by_bus[bus];
                    p_slack_target[bus] = std::clamp(
                        p_balance_by_bus[bus], -0.49, 0.49);
                    q_slack_target[bus] = std::clamp(
                        q_balance_by_bus[bus], -0.49, 0.49);
                }
                for (int bus = 0; bus < nb; ++bus) {
                    p_slack_target[bus] = std::clamp(
                        p_balance_by_bus[bus], -0.49, 0.49);
                    q_slack_target[bus] = std::clamp(
                        q_balance_by_bus[bus], -0.49, 0.49);
                }
                const auto distribute_slack_sum = [](
                    const std::vector<int>& component,
                    const std::vector<double>& balance,
                    std::vector<double>& target) {
                    double difference = 0.0;
                    for (int bus : component) {
                        difference += balance[bus] - target[bus];
                    }
                    for (int pass = 0;
                         pass < 4 && std::abs(difference) > 1e-10;
                         ++pass) {
                        double room = 0.0;
                        for (int bus : component) {
                            room += difference > 0.0
                                ? std::max(0.0, 0.49 - target[bus])
                                : std::max(0.0, target[bus] + 0.49);
                        }
                        if (room <= 1e-12) {
                            break;
                        }
                        const double prior_difference = difference;
                        double applied = 0.0;
                        for (int bus : component) {
                            const double individual_room =
                                prior_difference > 0.0
                                ? std::max(0.0, 0.49 - target[bus])
                                : std::max(0.0, target[bus] + 0.49);
                            const double change =
                                prior_difference * individual_room / room;
                            target[bus] = std::clamp(
                                target[bus] + change, -0.49, 0.49);
                            applied += change;
                        }
                        difference -= applied;
                    }
                };
                for (const auto& component : components) {
                    distribute_slack_sum(
                        component, p_balance_by_bus, p_slack_target);
                    distribute_slack_sum(
                        component, q_balance_by_bus, q_slack_target);
                }
                for (int bus = 0; bus < nb; ++bus) {
                    p_spec[bus] = generated_p_by_bus[bus] -
                        load_p_by_bus[bus] + p_slack_target[bus];
                    q_spec[bus] = generated_q_by_bus[bus] -
                        load_q_by_bus[bus] + q_slack_target[bus];
                }
                // A full chord step is excellent for the first outage
                // correction but often overshoots once the nonlinear state is
                // near the soft balance band.  Reuse the same resident
                // factorization for a short deterministic backtracking line
                // search and accept only a validation-improving step.
                const AcState correction_reference = predictor_state;
                AcState selected_correction = predictor_state;
                ValidationReport selected_validation = predictor_validation;
                double selected_damping = 0.0;
                std::string selected_correction_mode = "none";
                constexpr std::array<double, 5> kDampingCandidates{
                    1.0, 0.5, 0.25, 0.125, 0.0625};
                const auto strongly_improving_full_damping = [&]
                    (double damping, const ValidationReport& validation) {
                    return branch_outage_low_rank_update && damping == 1.0 &&
                        validation.max_residual <=
                            0.5 * predictor_validation.max_residual;
                };
                const bool active_block_dominant =
                    predictor_validation.worst_category ==
                        "active_balance" &&
                    predictor_validation.max_active_balance_residual > 1e-4;
                int worst_active_flow_branch = -1;
                bool worst_active_flow_from_side = true;
                double worst_active_flow_excess = 0.0;
                int worst_reactive_flow_branch = -1;
                bool worst_reactive_flow_from_side = true;
                double worst_reactive_flow_excess = 0.0;
                for (int branch_index = 0;
                     branch_index < static_cast<int>(data_.branches.size());
                     ++branch_index) {
                    if (branch_index == outaged_branch ||
                        data_.branches[branch_index].status == 0) {
                        continue;
                    }
                    const double rating =
                        data_.branches[branch_index].rate_c;
                    const double from_excess =
                        std::abs(predictor_state.pf[branch_index]) - rating;
                    const double to_excess =
                        std::abs(predictor_state.pt[branch_index]) - rating;
                    const double reactive_from_excess =
                        std::abs(predictor_state.qf[branch_index]) - rating;
                    const double reactive_to_excess =
                        std::abs(predictor_state.qt[branch_index]) - rating;
                    if (from_excess > worst_active_flow_excess) {
                        worst_active_flow_excess = from_excess;
                        worst_active_flow_branch = branch_index;
                        worst_active_flow_from_side = true;
                    }
                    if (to_excess > worst_active_flow_excess) {
                        worst_active_flow_excess = to_excess;
                        worst_active_flow_branch = branch_index;
                        worst_active_flow_from_side = false;
                    }
                    if (reactive_from_excess >
                        worst_reactive_flow_excess) {
                        worst_reactive_flow_excess =
                            reactive_from_excess;
                        worst_reactive_flow_branch = branch_index;
                        worst_reactive_flow_from_side = true;
                    }
                    if (reactive_to_excess >
                        worst_reactive_flow_excess) {
                        worst_reactive_flow_excess = reactive_to_excess;
                        worst_reactive_flow_branch = branch_index;
                        worst_reactive_flow_from_side = false;
                    }
                }
                const bool active_flow_bound_dominant =
                    predictor_validation.worst_category ==
                        "variable_bound" &&
                    worst_active_flow_branch >= 0 &&
                    worst_active_flow_excess > 1e-8 &&
                    worst_active_flow_excess >=
                        worst_reactive_flow_excess - 1e-10;
                const bool reactive_flow_bound_dominant =
                    predictor_validation.worst_category ==
                        "variable_bound" &&
                    worst_reactive_flow_branch >= 0 &&
                    worst_reactive_flow_excess > 1e-8 &&
                    worst_reactive_flow_excess >
                        worst_active_flow_excess + 1e-10;
                const auto project_trial_reactive_and_validate =
                    [&](AcState& trial) {
                    rebuild_contingency_state_derived_fields(
                        data_, base_state_, commitment_, *contingency, trial);
                    std::vector<double> trial_p_network(
                        static_cast<std::size_t>(nb), 0.0);
                    std::vector<double> trial_q_network(
                        static_cast<std::size_t>(nb), 0.0);
                    for (int branch_index = 0;
                         branch_index < static_cast<int>(data_.branches.size());
                         ++branch_index) {
                        if (branch_index == outaged_branch ||
                            data_.branches[branch_index].status == 0) {
                            continue;
                        }
                        const auto& branch = data_.branches[branch_index];
                        trial_p_network[branch.from] +=
                            trial.pf[branch_index];
                        trial_p_network[branch.to] +=
                            trial.pt[branch_index];
                        trial_q_network[branch.from] +=
                            trial.qf[branch_index];
                        trial_q_network[branch.to] +=
                            trial.qt[branch_index];
                    }
                    for (int shunt_index = 0;
                         shunt_index < static_cast<int>(data_.shunts.size());
                         ++shunt_index) {
                        const auto& shunt = data_.shunts[shunt_index];
                        trial_p_network[shunt.bus] +=
                            shunt.gs * trial.vm[shunt.bus] *
                            trial.vm[shunt.bus];
                        trial_q_network[shunt.bus] -=
                            effective_shunt_susceptance(
                                data_, trial, shunt_index) *
                            trial.vm[shunt.bus] *
                            trial.vm[shunt.bus];
                    }
                    std::vector<double> trial_load_active(
                        data_.loads.size(), 0.0);
                    for (int load = 0;
                         load < static_cast<int>(data_.loads.size());
                         ++load) {
                        trial_load_active[load] =
                            data_.loads[load].pd_nominal *
                            trial.demand_factor[load];
                    }
                    for (int bus = 0; bus < nb; ++bus) {
                        double generated_p = 0.0;
                        double loaded_p = 0.0;
                        for (int generator : active_at_bus[bus]) {
                            generated_p += trial.pg[generator];
                        }
                        for (int load : data_.buses[bus].loads) {
                            loaded_p += trial_load_active[load];
                        }
                        double p_balance =
                            trial_p_network[bus] - generated_p + loaded_p;
                        double p_excess = std::copysign(
                            std::max(
                                0.0, std::abs(p_balance) - 0.49),
                            p_balance);
                        if (std::abs(p_excess) <= 1e-12) {
                            continue;
                        }
                        if (!active_at_bus[bus].empty()) {
                            double total_lower = 0.0;
                            double total_upper = 0.0;
                            for (int generator : active_at_bus[bus]) {
                                total_lower += p_lower[generator];
                                total_upper += p_upper[generator];
                            }
                            const double target_generation = std::clamp(
                                generated_p + p_excess,
                                total_lower, total_upper);
                            const auto preferred_pg = trial.pg;
                            auto proposed_pg = trial.pg;
                            if (allocate_total(
                                    active_at_bus[bus], p_lower, p_upper,
                                    preferred_pg, target_generation,
                                    proposed_pg)) {
                                trial.pg = std::move(proposed_pg);
                                generated_p = 0.0;
                                for (int generator : active_at_bus[bus]) {
                                    generated_p += trial.pg[generator];
                                }
                                p_balance = trial_p_network[bus] -
                                    generated_p + loaded_p;
                                p_excess = std::copysign(
                                    std::max(
                                        0.0,
                                        std::abs(p_balance) - 0.49),
                                    p_balance);
                            }
                        }
                        if (std::abs(p_excess) <= 1e-12 ||
                            data_.buses[bus].loads.empty()) {
                            continue;
                        }
                        double total_lower = 0.0;
                        double total_upper = 0.0;
                        for (int load : data_.buses[bus].loads) {
                            total_lower +=
                                predictor_load_power_lower[load];
                            total_upper +=
                                predictor_load_power_upper[load];
                        }
                        const double target_load = std::clamp(
                            loaded_p - p_excess,
                            total_lower, total_upper);
                        const auto preferred_load = trial_load_active;
                        if (!allocate_total(
                                data_.buses[bus].loads,
                                predictor_load_power_lower,
                                predictor_load_power_upper,
                                preferred_load, target_load,
                                trial_load_active)) {
                            continue;
                        }
                        for (int load : data_.buses[bus].loads) {
                            if (std::abs(
                                    data_.loads[load].pd_nominal) > 1e-12) {
                                trial.demand_factor[load] =
                                    trial_load_active[load] /
                                    data_.loads[load].pd_nominal;
                            }
                        }
                    }
                    for (int bus = 0; bus < nb; ++bus) {
                        if (active_at_bus[bus].empty()) {
                            continue;
                        }
                        double generated_q = 0.0;
                        double load_at_bus_q = 0.0;
                        for (int generator : active_at_bus[bus]) {
                            generated_q += trial.qg[generator];
                        }
                        for (int load : data_.buses[bus].loads) {
                            load_at_bus_q += data_.loads[load].qd_nominal *
                                trial.demand_factor[load];
                        }
                        const double q_balance =
                            trial_q_network[bus] - generated_q +
                            load_at_bus_q;
                        auto projected_qg = trial.qg;
                        if (allocate_total(
                                active_at_bus[bus], q_lower, q_upper,
                                initial_state.qg,
                                generated_q + q_balance -
                                    std::clamp(q_balance, -0.49, 0.49),
                                projected_qg)) {
                            trial.qg = std::move(projected_qg);
                        }
                    }
                    rebuild_contingency_state_derived_fields(
                        data_, base_state_, commitment_, *contingency, trial);
                    return validate_state(
                        data_, ModelMode::ContingencySoft, trial,
                        commitment_, direct_context);
                };
                const auto try_damped_corrections = [&]
                    (FixedJacobianPredictorCache& cache) {
                    if (active_flow_bound_dominant) {
                        const auto& branch =
                            data_.branches[worst_active_flow_branch];
                        const int from = branch.from;
                        const int to = branch.to;
                        const double denominator =
                            branch.r * branch.r + branch.x * branch.x;
                        const double g = denominator > 1e-20
                            ? branch.r / denominator : 0.0;
                        const double b = denominator > 1e-20
                            ? -branch.x / denominator : 0.0;
                        const double tm2 = branch.tap * branch.tap;
                        const double tr =
                            branch.tap * std::cos(branch.shift);
                        const double ti =
                            branch.tap * std::sin(branch.shift);
                        const double angle =
                            correction_reference.va[from] -
                            correction_reference.va[to];
                        const double voltage_product =
                            correction_reference.vm[from] *
                            correction_reference.vm[to];
                        const double flow = worst_active_flow_from_side
                            ? correction_reference
                                  .pf[worst_active_flow_branch]
                            : correction_reference
                                  .pt[worst_active_flow_branch];
                        double derivative = 0.0;
                        if (worst_active_flow_from_side) {
                            const double cosine_coefficient =
                                ((-g * tr + b * ti) / tm2) *
                                voltage_product;
                            const double sine_coefficient =
                                ((-b * tr - g * ti) / tm2) *
                                voltage_product;
                            derivative =
                                -cosine_coefficient * std::sin(angle) +
                                sine_coefficient * std::cos(angle);
                        } else {
                            const double cosine_coefficient =
                                ((-g * tr - b * ti) / tm2) *
                                voltage_product;
                            const double sine_coefficient =
                                ((-b * tr + g * ti) / tm2) *
                                voltage_product;
                            derivative =
                                -cosine_coefficient * std::sin(angle) -
                                sine_coefficient * std::cos(angle);
                        }
                        if (std::abs(derivative) > 1e-10) {
                            const double target_flow = std::copysign(
                                std::max(0.0, branch.rate_c - 1e-4),
                                flow);
                            const double raw_angle_change = std::clamp(
                                (target_flow - flow) / derivative,
                                -0.05, 0.05);
                            for (const double damping :
                                 kDampingCandidates) {
                                auto trial = correction_reference;
                                const double angle_change =
                                    damping * raw_angle_change;
                                if (slack[from]) {
                                    trial.va[to] -= angle_change;
                                } else if (slack[to]) {
                                    trial.va[from] += angle_change;
                                } else {
                                    trial.va[from] += 0.5 * angle_change;
                                    trial.va[to] -= 0.5 * angle_change;
                                }
                                const auto trial_validation =
                                    project_trial_reactive_and_validate(trial);
                                if (trial_validation.max_residual + 1e-10 <
                                    selected_validation.max_residual) {
                                    selected_correction = std::move(trial);
                                    selected_validation = trial_validation;
                                    selected_damping = damping;
                                    selected_correction_mode =
                                        "active_branch_flow_angle";
                                    if (strongly_improving_full_damping(
                                            damping, trial_validation)) {
                                        break;
                                    }
                                }
                            }
                        }
                        const double target_flow = std::copysign(
                            std::max(0.0, branch.rate_c - 1e-4), flow);
                        auto full_redispatch_pg =
                            correction_reference.pg;
                        auto full_redispatch_p_spec = p_spec;
                        AcState best_redispatch = correction_reference;
                        ValidationReport best_redispatch_validation =
                            predictor_validation;
                        double best_redispatch_damping = 0.0;
                        if (cache.apply_active_flow_redispatch(
                                data_, worst_active_flow_branch,
                                worst_active_flow_from_side,
                                target_flow - flow,
                                p_lower, p_upper,
                                correction_reference,
                                full_redispatch_pg,
                                full_redispatch_p_spec)) {
                            for (const double damping :
                                 kDampingCandidates) {
                                auto trial = correction_reference;
                                auto trial_p_spec = p_spec;
                                for (int generator = 0;
                                     generator < ng; ++generator) {
                                    trial.pg[generator] += damping *
                                        (full_redispatch_pg[generator] -
                                         correction_reference.pg[generator]);
                                }
                                for (int bus = 0; bus < nb; ++bus) {
                                    trial_p_spec[bus] += damping *
                                        (full_redispatch_p_spec[bus] -
                                         p_spec[bus]);
                                }
                                if (!cache.apply_active_correction(
                                        data_, trial_p_spec, p_network,
                                        trial.va, 1.0)) {
                                    continue;
                                }
                                const auto trial_validation =
                                    project_trial_reactive_and_validate(
                                        trial);
                                if (trial_validation.max_residual + 1e-10 <
                                    best_redispatch_validation.max_residual) {
                                    best_redispatch = std::move(trial);
                                    best_redispatch_validation =
                                        trial_validation;
                                    best_redispatch_damping = damping;
                                    if (strongly_improving_full_damping(
                                            damping, trial_validation)) {
                                        break;
                                    }
                                }
                            }
                        }
                        // A network-wide redispatch is generally a more
                        // stable continuation step than changing only the two
                        // endpoint angles.  Prefer any independently
                        // improving redispatch; retain the endpoint move as a
                        // fallback when redispatch cannot improve the exact
                        // residual.
                        const bool prefer_redispatch_over_endpoint_angle =
                            active_flow_bound_dominant &&
                            (worst_reactive_flow_excess <= 1e-6 ||
                             (predictor_iteration >= 6 &&
                              predictor_validation.max_residual > 0.02));
                        if (best_redispatch_damping != 0.0 &&
                            (prefer_redispatch_over_endpoint_angle ||
                             best_redispatch_validation.max_residual +
                                     1e-10 <
                                 selected_validation.max_residual)) {
                            selected_correction =
                                std::move(best_redispatch);
                            selected_validation =
                                best_redispatch_validation;
                            selected_damping = best_redispatch_damping;
                            selected_correction_mode =
                                "active_branch_flow_redispatch";
                        }
                        if (selected_damping != 0.0 &&
                            active_flow_bound_dominant) {
                            return;
                        }
                    }
                    for (const double damping : kDampingCandidates) {
                        if (!active_block_dominant) {
                            break;
                        }
                        auto trial = correction_reference;
                        if (!cache.apply_active_correction(
                                data_, p_spec, p_network,
                                trial.va, damping)) {
                            continue;
                        }
                        const auto trial_validation =
                            project_trial_reactive_and_validate(trial);
                        if (trial_validation.max_residual + 1e-10 <
                            selected_validation.max_residual) {
                            selected_correction = std::move(trial);
                            selected_validation = trial_validation;
                            selected_damping = damping;
                            selected_correction_mode = "active_angle";
                            if (strongly_improving_full_damping(
                                    damping, trial_validation)) {
                                break;
                            }
                        }
                    }
                    if (active_block_dominant) {
                        return;
                    }
                    for (const double damping : kDampingCandidates) {
                        if (active_block_dominant) {
                            break;
                        }
                        auto trial = correction_reference;
                        if (!cache.apply_local_reactive_least_squares(
                                data_, q_balance_by_bus,
                                trial.vm, damping)) {
                            continue;
                        }
                        const auto trial_validation =
                            project_trial_reactive_and_validate(trial);
                        if (trial_validation.max_residual + 1e-10 <
                            selected_validation.max_residual) {
                            selected_correction = std::move(trial);
                            selected_validation = trial_validation;
                            selected_damping = damping;
                            selected_correction_mode =
                                "local_reactive_least_squares";
                            if (strongly_improving_full_damping(
                                    damping, trial_validation)) {
                                break;
                            }
                        }
                    }
                    for (const double damping : kDampingCandidates) {
                        if (active_block_dominant) {
                            break;
                        }
                        auto trial = correction_reference;
                        if (!cache.apply_local_reactive_least_squares(
                                data_, q_balance_by_bus,
                                trial.vm, damping)) {
                            continue;
                        }
                        rebuild_contingency_state_derived_fields(
                            data_, base_state_, commitment_, *contingency,
                            trial);
                        std::vector<double> trial_p_network(
                            static_cast<std::size_t>(nb), 0.0);
                        for (int branch_index = 0;
                             branch_index <
                                 static_cast<int>(data_.branches.size());
                             ++branch_index) {
                            if (branch_index == outaged_branch ||
                                data_.branches[branch_index].status == 0) {
                                continue;
                            }
                            const auto& branch =
                                data_.branches[branch_index];
                            trial_p_network[branch.from] +=
                                trial.pf[branch_index];
                            trial_p_network[branch.to] +=
                                trial.pt[branch_index];
                        }
                        for (const auto& shunt : data_.shunts) {
                            trial_p_network[shunt.bus] +=
                                shunt.gs * trial.vm[shunt.bus] *
                                trial.vm[shunt.bus];
                        }
                        if (!cache.apply_active_correction(
                                data_, p_spec, trial_p_network,
                                trial.va, 1.0)) {
                            continue;
                        }
                        const auto trial_validation =
                            project_trial_reactive_and_validate(trial);
                        if (trial_validation.max_residual + 1e-10 <
                            selected_validation.max_residual) {
                            selected_correction = std::move(trial);
                            selected_validation = trial_validation;
                            selected_damping = damping;
                            selected_correction_mode =
                                "local_reactive_least_squares_then_"
                                "active_angle";
                            if (strongly_improving_full_damping(
                                    damping, trial_validation)) {
                                break;
                            }
                        }
                    }
                    if (selected_damping != 0.0) {
                        return;
                    }
                    for (const double damping : kDampingCandidates) {
                        if (active_block_dominant) {
                            break;
                        }
                        auto trial = correction_reference;
                        if (!cache.apply_reactive_band_correction(
                                data_, q_balance_by_bus,
                                trial.vm, damping)) {
                            continue;
                        }
                        const auto trial_validation =
                            project_trial_reactive_and_validate(trial);
                        if (trial_validation.max_residual + 1e-10 <
                            selected_validation.max_residual) {
                            selected_correction = std::move(trial);
                            selected_validation = trial_validation;
                            selected_damping = damping;
                            selected_correction_mode =
                                "reactive_feasibility_band";
                            if (strongly_improving_full_damping(
                                    damping, trial_validation)) {
                                break;
                            }
                        }
                    }
                    for (const double damping : kDampingCandidates) {
                        if (active_block_dominant) {
                            break;
                        }
                        auto trial = correction_reference;
                        if (!cache.apply_reactive_correction(
                                data_, q_spec, q_network,
                                trial.vm, damping)) {
                            continue;
                        }
                        const auto trial_validation =
                            project_trial_reactive_and_validate(trial);
                        if (trial_validation.max_residual + 1e-10 <
                            selected_validation.max_residual) {
                            selected_correction = std::move(trial);
                            selected_validation = trial_validation;
                            selected_damping = damping;
                            selected_correction_mode = "reactive_voltage";
                            if (strongly_improving_full_damping(
                                    damping, trial_validation)) {
                                break;
                            }
                        }
                    }
                    for (const double damping : kDampingCandidates) {
                        if (active_block_dominant) {
                            break;
                        }
                        auto trial = correction_reference;
                        if (!cache.apply_reactive_correction(
                                data_, q_spec, q_network,
                                trial.vm, damping)) {
                            continue;
                        }
                        rebuild_contingency_state_derived_fields(
                            data_, base_state_, commitment_, *contingency,
                            trial);
                        std::vector<double> trial_p_network(
                            static_cast<std::size_t>(nb), 0.0);
                        for (int branch_index = 0;
                             branch_index <
                                 static_cast<int>(data_.branches.size());
                             ++branch_index) {
                            if (branch_index == outaged_branch ||
                                data_.branches[branch_index].status == 0) {
                                continue;
                            }
                            const auto& branch =
                                data_.branches[branch_index];
                            trial_p_network[branch.from] +=
                                trial.pf[branch_index];
                            trial_p_network[branch.to] +=
                                trial.pt[branch_index];
                        }
                        for (const auto& shunt : data_.shunts) {
                            trial_p_network[shunt.bus] +=
                                shunt.gs * trial.vm[shunt.bus] *
                                trial.vm[shunt.bus];
                        }
                        if (!cache.apply_active_correction(
                                data_, p_spec, trial_p_network,
                                trial.va, 1.0)) {
                            continue;
                        }
                        const auto trial_validation =
                            project_trial_reactive_and_validate(trial);
                        if (trial_validation.max_residual + 1e-10 <
                            selected_validation.max_residual) {
                            selected_correction = std::move(trial);
                            selected_validation = trial_validation;
                            selected_damping = damping;
                            selected_correction_mode =
                                "reactive_voltage_then_active_angle";
                            if (strongly_improving_full_damping(
                                    damping, trial_validation)) {
                                break;
                            }
                        }
                    }
                    for (const double damping : kDampingCandidates) {
                        if (active_block_dominant) {
                            break;
                        }
                        auto trial = correction_reference;
                        if (!cache.apply_correction(
                                data_, p_spec, q_spec, p_network, q_network,
                                trial.vm, trial.va, damping)) {
                            continue;
                        }
                        const auto trial_validation =
                            project_trial_reactive_and_validate(trial);
                        if (trial_validation.max_residual + 1e-10 <
                            selected_validation.max_residual) {
                            selected_correction = std::move(trial);
                            selected_validation = trial_validation;
                            selected_damping = damping;
                            selected_correction_mode =
                                "coupled_angle_voltage";
                            if (strongly_improving_full_damping(
                                    damping, trial_validation)) {
                                break;
                            }
                        }
                    }
                };
                // The resident base-case Jacobian is intentionally cheap and
                // handles most outages.  A low-impedance branch outage can,
                // however, remove a dominant Jacobian term while the stale
                // factorization continues making small but misleading
                // progress.  Refactor only after four clearly slow steps;
                // this keeps the common path resident while giving hard
                // branch outages the correct local derivatives.
                const bool slow_branch_outage =
                    outaged_branch >= 0 && predictor_iteration >= 4 &&
                    predictor_validation.max_residual > 0.1 &&
                    !output.fixed_jacobian_predictor_trace.empty() &&
                    predictor_validation.max_residual >
                        0.25 * output.fixed_jacobian_predictor_trace.front()
                            .at("validation")
                            .at("max_residual")
                            .get<double>();
                if (!contingency_predictor_cache && slow_branch_outage) {
                    contingency_predictor_cache =
                        std::make_unique<FixedJacobianPredictorCache>(
                            data_, correction_reference, commitment_,
                            outaged_branch);
                    output.fixed_jacobian_predictor_preparation_seconds +=
                        contingency_predictor_cache->preparation_seconds;
                    output.fixed_jacobian_predictor_trace.back()[
                        "contingency_specific_refactorization"] = true;
                    output.fixed_jacobian_predictor_trace.back()[
                        "contingency_specific_refactorization_reason"] =
                        "slow_branch_outage_progress";
                }
                if (contingency_predictor_cache) {
                    try_damped_corrections(*contingency_predictor_cache);
                } else {
                    try_damped_corrections(*predictor_cache_);
                    if (selected_damping == 0.0) {
                        contingency_predictor_cache =
                            std::make_unique<FixedJacobianPredictorCache>(
                                data_, correction_reference, commitment_,
                                outaged_branch);
                        output.fixed_jacobian_predictor_preparation_seconds +=
                            contingency_predictor_cache->preparation_seconds;
                        output.fixed_jacobian_predictor_trace.back()[
                            "contingency_specific_refactorization"] = true;
                        if (contingency_predictor_cache->valid) {
                            try_damped_corrections(
                                *contingency_predictor_cache);
                        }
                    }
                }
                if (selected_damping == 0.0) {
                    int worst_q_bus = -1;
                    double worst_q_excess = 0.0;
                    for (int bus = 0; bus < nb; ++bus) {
                        const double excess =
                            std::abs(q_balance_by_bus[bus]) - 0.49;
                        if (excess > worst_q_excess) {
                            worst_q_excess = excess;
                            worst_q_bus = bus;
                        }
                    }
                    std::vector<int> shunt_candidate_buses;
                    if (worst_q_bus >= 0) {
                        shunt_candidate_buses.push_back(worst_q_bus);
                        for (int branch_index :
                             data_.buses[worst_q_bus].branches_from) {
                            if (branch_index != outaged_branch &&
                                data_.branches[branch_index].status != 0) {
                                shunt_candidate_buses.push_back(
                                    data_.branches[branch_index].to);
                            }
                        }
                        for (int branch_index :
                             data_.buses[worst_q_bus].branches_to) {
                            if (branch_index != outaged_branch &&
                                data_.branches[branch_index].status != 0) {
                                shunt_candidate_buses.push_back(
                                    data_.branches[branch_index].from);
                            }
                        }
                        std::sort(
                            shunt_candidate_buses.begin(),
                            shunt_candidate_buses.end());
                        shunt_candidate_buses.erase(
                            std::unique(
                                shunt_candidate_buses.begin(),
                                shunt_candidate_buses.end()),
                            shunt_candidate_buses.end());
                    }
                    if (reactive_flow_bound_dominant) {
                        const auto& reactive_branch =
                            data_.branches[worst_reactive_flow_branch];
                        shunt_candidate_buses.push_back(
                            reactive_branch.from);
                        shunt_candidate_buses.push_back(
                            reactive_branch.to);
                        std::sort(
                            shunt_candidate_buses.begin(),
                            shunt_candidate_buses.end());
                        shunt_candidate_buses.erase(
                            std::unique(
                                shunt_candidate_buses.begin(),
                                shunt_candidate_buses.end()),
                            shunt_candidate_buses.end());
                        output.fixed_jacobian_predictor_trace.back()[
                            "worst_reactive_flow_branch"] =
                            reactive_branch.source_key;
                        output.fixed_jacobian_predictor_trace.back()[
                            "worst_reactive_flow_side"] =
                            worst_reactive_flow_from_side
                            ? "from" : "to";
                        output.fixed_jacobian_predictor_trace.back()[
                            "worst_reactive_flow_excess"] =
                            worst_reactive_flow_excess;
                    }
                    auto* shunt_correction_cache =
                        contingency_predictor_cache &&
                            contingency_predictor_cache->valid
                        ? contingency_predictor_cache.get()
                        : predictor_cache_.get();
                    auto shunt_trial_trace = nlohmann::json::array();
                    output.fixed_jacobian_predictor_trace.back()[
                        "worst_q_bus"] = worst_q_bus >= 0
                        ? data_.buses[worst_q_bus].bus_i : -1;
                    output.fixed_jacobian_predictor_trace.back()[
                        "worst_q_balance"] = worst_q_bus >= 0
                        ? q_balance_by_bus[worst_q_bus] : 0.0;
                    auto reactive_flow_voltage_trace =
                        nlohmann::json::array();
                    if (reactive_flow_bound_dominant) {
                        const auto& reactive_branch =
                            data_.branches[worst_reactive_flow_branch];
                        constexpr std::array<double, 6>
                            kReactiveFlowVoltageChanges{
                                -0.0005, -0.00025, -0.0001,
                                 0.0001,  0.00025,  0.0005};
                        for (int candidate_bus : {
                                 reactive_branch.from,
                                 reactive_branch.to}) {
                            for (double voltage_change :
                                 kReactiveFlowVoltageChanges) {
                                const double proposed_voltage =
                                    correction_reference.vm[candidate_bus] +
                                    voltage_change;
                                if (proposed_voltage <
                                        data_.buses[candidate_bus].vmin -
                                            1e-12 ||
                                    proposed_voltage >
                                        data_.buses[candidate_bus].vmax +
                                            1e-12) {
                                    continue;
                                }
                                auto trial = correction_reference;
                                trial.vm[candidate_bus] = proposed_voltage;
                                const auto trial_validation =
                                    project_trial_reactive_and_validate(
                                        trial);
                                reactive_flow_voltage_trace.push_back({
                                    {"candidate_bus",
                                     data_.buses[candidate_bus].bus_i},
                                    {"voltage_change", voltage_change},
                                    {"validation",
                                     trial_validation.to_json()},
                                });
                                if (trial_validation.max_residual +
                                        1e-10 <
                                    selected_validation.max_residual) {
                                    selected_correction =
                                        std::move(trial);
                                    selected_validation =
                                        trial_validation;
                                    selected_damping = voltage_change;
                                    selected_correction_mode =
                                        "reactive_branch_flow_voltage";
                                }
                            }
                        }
                    }
                    output.fixed_jacobian_predictor_trace.back()[
                        "reactive_flow_voltage_trials"] =
                        std::move(reactive_flow_voltage_trace);
                    auto local_reactive_trace = nlohmann::json::array();
                    if (shunt_correction_cache != nullptr &&
                        worst_q_bus >= 0) {
                        const double desired_network_q_change =
                            std::clamp(
                                q_balance_by_bus[worst_q_bus],
                                -0.48, 0.48) -
                            q_balance_by_bus[worst_q_bus];
                        constexpr std::array<double, 7>
                            kLocalReactiveDamping{
                                1.0, 0.5, 0.25, 0.125, 0.0625,
                                0.03125, 0.015625};
                        for (double direction : {1.0, -1.0}) {
                        for (double damping : kLocalReactiveDamping) {
                            auto trial = correction_reference;
                            if (!shunt_correction_cache->
                                    apply_local_reactive_correction(
                                        data_, worst_q_bus,
                                        direction *
                                            desired_network_q_change,
                                        trial.vm, damping)) {
                                continue;
                            }
                            const auto trial_validation =
                                project_trial_reactive_and_validate(trial);
                            local_reactive_trace.push_back({
                                {"damping", damping},
                                {"direction", direction},
                                {"desired_network_q_change",
                                 direction *
                                     desired_network_q_change},
                                {"validation", trial_validation.to_json()},
                            });
                            if (trial_validation.max_residual + 1e-10 <
                                selected_validation.max_residual) {
                                selected_correction = std::move(trial);
                                selected_validation = trial_validation;
                                selected_damping = damping;
                                selected_correction_mode =
                                    "localized_reactive_voltage";
                            }
                        }
                        }
                    }
                    output.fixed_jacobian_predictor_trace.back()[
                        "local_reactive_trials"] =
                        std::move(local_reactive_trace);
                    auto local_voltage_trace = nlohmann::json::array();
                    constexpr std::array<double, 4>
                        kLocalVoltageChanges{
                            -0.005, -0.001, 0.001, 0.005};
                    for (int candidate_bus : shunt_candidate_buses) {
                        for (double voltage_change :
                             kLocalVoltageChanges) {
                            const double proposed_voltage =
                                correction_reference.vm[candidate_bus] +
                                voltage_change;
                            if (proposed_voltage <
                                    data_.buses[candidate_bus].vmin - 1e-12 ||
                                proposed_voltage >
                                    data_.buses[candidate_bus].vmax + 1e-12) {
                                continue;
                            }
                            auto trial = correction_reference;
                            trial.vm[candidate_bus] = proposed_voltage;
                            const auto trial_validation =
                                project_trial_reactive_and_validate(trial);
                            local_voltage_trace.push_back({
                                {"candidate_bus",
                                 data_.buses[candidate_bus].bus_i},
                                {"voltage_change", voltage_change},
                                {"validation", trial_validation.to_json()},
                            });
                            if (trial_validation.max_residual + 1e-10 <
                                selected_validation.max_residual) {
                                selected_correction = std::move(trial);
                                selected_validation = trial_validation;
                                selected_damping = voltage_change;
                                selected_correction_mode =
                                    "one_hop_voltage_coordinate";
                            }
                        }
                    }
                    output.fixed_jacobian_predictor_trace.back()[
                        "local_voltage_trials"] =
                        std::move(local_voltage_trace);
                    for (int candidate_bus : shunt_candidate_buses) {
                    for (int shunt_index :
                         data_.buses[candidate_bus].shunts) {
                        const auto& shunt = data_.shunts[shunt_index];
                        if (!shunt.dispatchable ||
                            correction_reference.shunt_steps.size() !=
                                data_.shunts.size() ||
                            correction_reference.shunt_bs.size() !=
                                data_.shunts.size()) {
                            continue;
                        }
                        for (int block = 0;
                             block < static_cast<int>(
                                 shunt.block_maximum_steps.size());
                             ++block) {
                            for (const int step_change : {-1, 1}) {
                                const int current_step =
                                    correction_reference
                                        .shunt_steps[shunt_index][block];
                                const int proposed_step =
                                    current_step + step_change;
                                if (proposed_step < 0 ||
                                    proposed_step >
                                        shunt.block_maximum_steps[block]) {
                                    continue;
                                }
                                const double delta_bs =
                                    static_cast<double>(step_change) *
                                    shunt.block_susceptance[block];
                                nlohmann::json shunt_trial = {
                                    {"candidate_bus",
                                     data_.buses[candidate_bus].bus_i},
                                    {"shunt", shunt.source_key},
                                    {"block", block},
                                    {"current_step", current_step},
                                    {"proposed_step", proposed_step},
                                    {"delta_bs", delta_bs},
                                };
                                // q_balance includes -bs*|V|^2, so positive
                                // excess is relieved by increasing bs and
                                // negative excess by decreasing it.
                                if (worst_q_bus < 0 ||
                                    (candidate_bus == worst_q_bus &&
                                     q_balance_by_bus[worst_q_bus] *
                                         delta_bs <= 0.0)) {
                                    shunt_trial["skipped"] =
                                        "opposite_reactive_direction";
                                    shunt_trial_trace.push_back(
                                        std::move(shunt_trial));
                                    continue;
                                }
                                auto trial = correction_reference;
                                trial.shunt_steps[shunt_index][block] =
                                    proposed_step;
                                double trial_bs = 0.0;
                                for (int candidate_block = 0;
                                     candidate_block < static_cast<int>(
                                         shunt.block_maximum_steps.size());
                                     ++candidate_block) {
                                    trial_bs += static_cast<double>(
                                        trial.shunt_steps[shunt_index]
                                            [candidate_block]) *
                                        shunt.block_susceptance
                                            [candidate_block];
                                }
                                trial.shunt_bs[shunt_index] = trial_bs;
                                auto trial_validation =
                                    project_trial_reactive_and_validate(trial);
                                shunt_trial["direct_validation"] =
                                    trial_validation.to_json();
                                if (trial_validation.max_residual + 1e-10 <
                                    selected_validation.max_residual) {
                                    // Keep trial intact for the optional
                                    // neighbor-shunt/reactive correction that
                                    // follows. Moving it here leaves its state
                                    // vectors empty and caused the fast worker
                                    // to dereference a moved-from candidate.
                                    selected_correction = trial;
                                    selected_validation = trial_validation;
                                    selected_damping = 1.0;
                                    selected_correction_mode =
                                        "single_switched_shunt_step";
                                }
                                std::vector<double> trial_q_network(
                                    static_cast<std::size_t>(nb), 0.0);
                                for (int branch_index = 0;
                                     branch_index < static_cast<int>(
                                         data_.branches.size());
                                     ++branch_index) {
                                    if (branch_index == outaged_branch ||
                                        data_.branches[branch_index].status ==
                                            0) {
                                        continue;
                                    }
                                    const auto& branch =
                                        data_.branches[branch_index];
                                    trial_q_network[branch.from] +=
                                        trial.qf[branch_index];
                                    trial_q_network[branch.to] +=
                                        trial.qt[branch_index];
                                }
                                for (int candidate_shunt = 0;
                                     candidate_shunt < static_cast<int>(
                                         data_.shunts.size());
                                     ++candidate_shunt) {
                                    const auto& network_shunt =
                                        data_.shunts[candidate_shunt];
                                    trial_q_network[network_shunt.bus] -=
                                        effective_shunt_susceptance(
                                            data_, trial,
                                            candidate_shunt) *
                                        trial.vm[network_shunt.bus] *
                                        trial.vm[network_shunt.bus];
                                }
                                // The discrete shunt step changes both the
                                // network injection and the component-wide
                                // amount of admissible reactive-balance slack.
                                // Rebuild the voltage-correction target from
                                // this trial state instead of reusing the
                                // target computed before the shunt moved.
                                std::vector<double> trial_q_spec(
                                    static_cast<std::size_t>(nb), 0.0);
                                std::vector<double> trial_q_balance(
                                    static_cast<std::size_t>(nb), 0.0);
                                std::vector<double> trial_q_slack_target(
                                    static_cast<std::size_t>(nb), 0.0);
                                std::vector<double> trial_generated_q(
                                    static_cast<std::size_t>(nb), 0.0);
                                std::vector<double> trial_load_q(
                                    static_cast<std::size_t>(nb), 0.0);
                                for (int bus = 0; bus < nb; ++bus) {
                                    for (int generator : active_at_bus[bus]) {
                                        trial_generated_q[bus] +=
                                            trial.qg[generator];
                                    }
                                    for (int load : data_.buses[bus].loads) {
                                        trial_load_q[bus] +=
                                            data_.loads[load].qd_nominal *
                                            trial.demand_factor[load];
                                    }
                                    trial_q_balance[bus] =
                                        trial_q_network[bus] -
                                        trial_generated_q[bus] +
                                        trial_load_q[bus];
                                    trial_q_slack_target[bus] = std::clamp(
                                        trial_q_balance[bus], -0.49, 0.49);
                                }
                                for (const auto& component : components) {
                                    distribute_slack_sum(
                                        component, trial_q_balance,
                                        trial_q_slack_target);
                                }
                                for (int bus = 0; bus < nb; ++bus) {
                                    trial_q_spec[bus] =
                                        trial_generated_q[bus] -
                                        trial_load_q[bus] +
                                        trial_q_slack_target[bus];
                                }
                                constexpr std::array<double, 3>
                                    kShuntCorrectionDamping{
                                        1.0, 0.5, 0.25};
                                shunt_trial["reactive_corrections"] =
                                    nlohmann::json::array();
                                for (double damping :
                                     kShuntCorrectionDamping) {
                                    auto corrected_trial = trial;
                                    if (!shunt_correction_cache->
                                            apply_reactive_correction(
                                                data_, trial_q_spec,
                                                trial_q_network,
                                                corrected_trial.vm,
                                                damping)) {
                                        continue;
                                    }
                                    const auto corrected_validation =
                                        project_trial_reactive_and_validate(
                                            corrected_trial);
                                    shunt_trial["reactive_corrections"]
                                        .push_back({
                                            {"damping", damping},
                                            {"validation",
                                             corrected_validation.to_json()},
                                        });
                                    if (corrected_validation.max_residual +
                                            1e-10 <
                                        selected_validation.max_residual) {
                                        selected_correction =
                                            std::move(corrected_trial);
                                        selected_validation =
                                            corrected_validation;
                                        selected_damping = damping;
                                        selected_correction_mode =
                                            "neighbor_switched_shunt_plus_"
                                            "reactive_voltage";
                                    }
                                }
                                shunt_trial_trace.push_back(
                                    std::move(shunt_trial));
                            }
                        }
                    }
                }
                    auto paired_voltage_trace = nlohmann::json::array();
                    if (selected_damping == 0.0 &&
                        shunt_candidate_buses.size() > 1) {
                        constexpr std::array<double, 2>
                            kPairedVoltageChanges{-0.001, 0.001};
                        for (std::size_t first = 0;
                             first < shunt_candidate_buses.size(); ++first) {
                            for (std::size_t second = first + 1;
                                 second < shunt_candidate_buses.size();
                                 ++second) {
                                const int first_bus =
                                    shunt_candidate_buses[first];
                                const int second_bus =
                                    shunt_candidate_buses[second];
                                for (double first_change :
                                     kPairedVoltageChanges) {
                                    for (double second_change :
                                         kPairedVoltageChanges) {
                                        const double first_voltage =
                                            correction_reference.vm[first_bus] +
                                            first_change;
                                        const double second_voltage =
                                            correction_reference.vm[second_bus] +
                                            second_change;
                                        if (first_voltage <
                                                data_.buses[first_bus].vmin -
                                                    1e-12 ||
                                            first_voltage >
                                                data_.buses[first_bus].vmax +
                                                    1e-12 ||
                                            second_voltage <
                                                data_.buses[second_bus].vmin -
                                                    1e-12 ||
                                            second_voltage >
                                                data_.buses[second_bus].vmax +
                                                    1e-12) {
                                            continue;
                                        }
                                        auto trial = correction_reference;
                                        trial.vm[first_bus] = first_voltage;
                                        trial.vm[second_bus] = second_voltage;
                                        const auto trial_validation =
                                            project_trial_reactive_and_validate(
                                                trial);
                                        paired_voltage_trace.push_back({
                                            {"first_bus",
                                             data_.buses[first_bus].bus_i},
                                            {"first_change", first_change},
                                            {"second_bus",
                                             data_.buses[second_bus].bus_i},
                                            {"second_change", second_change},
                                            {"validation",
                                             trial_validation.to_json()},
                                        });
                                        if (trial_validation.max_residual +
                                                1e-10 <
                                            selected_validation.max_residual) {
                                            selected_correction =
                                                std::move(trial);
                                            selected_validation =
                                                trial_validation;
                                            selected_damping = 1.0;
                                            selected_correction_mode =
                                                "paired_one_hop_voltage_"
                                                "coordinate";
                                        }
                                    }
                                }
                            }
                        }
                    }
                    output.fixed_jacobian_predictor_trace.back()[
                        "paired_voltage_trials"] =
                        std::move(paired_voltage_trace);
                    output.fixed_jacobian_predictor_trace.back()[
                        "shunt_trials"] = std::move(shunt_trial_trace);
                }
                output.fixed_jacobian_predictor_trace.back()[
                    "selected_damping"] = selected_damping;
                output.fixed_jacobian_predictor_trace.back()[
                    "selected_correction_mode"] =
                    selected_correction_mode;
                output.fixed_jacobian_predictor_trace.back()[
                    "selected_next_validation"] =
                    selected_validation.to_json();
                const bool initial_active_feasibility_repair =
                    active_feasibility_repair_attempts == 0 &&
                    selected_damping == 0.0 &&
                    (predictor_validation.max_active_balance_residual >
                         1e-4 ||
                     worst_active_flow_excess > 1e-6);
                const bool late_security_stagnation_repair =
                    active_feasibility_repair_attempts == 1 &&
                    predictor_iteration >= 28 &&
                    predictor_validation.worst_category ==
                        "variable_bound" &&
                    predictor_validation.max_variable_bound_violation >
                        1e-4 &&
                    selected_validation.max_residual >=
                        0.99 * predictor_validation.max_residual;
                if (outaged_generator >= 0 &&
                    (initial_active_feasibility_repair ||
                     late_security_stagnation_repair)) {
                    ++active_feasibility_repair_attempts;
                    const double active_balance_limit =
                        initial_active_feasibility_repair ? 0.25 : 0.20;
                    const double active_angle_trust =
                        initial_active_feasibility_repair ? 0.05 : 0.02;
                    const double active_voltage_trust =
                        initial_active_feasibility_repair ? 0.005 : 0.002;
                    auto active_repair =
                        solve_linearized_active_feasibility_repair(
                            data_, correction_reference, commitment_,
                            *direct_context, active_balance_limit,
                            active_angle_trust, 5.0,
                            active_voltage_trust, false);
                    auto active_repair_json =
                        active_repair.to_json(false);
                    if (active_repair.success) {
                        auto active_trial = std::move(active_repair.state);
                        const auto active_validation =
                            project_trial_reactive_and_validate(active_trial);
                        active_repair_json["nonlinear_validation"] =
                            active_validation.to_json();
                        if (active_validation.max_residual + 1e-10 <
                            selected_validation.max_residual) {
                            selected_correction =
                                std::move(active_trial);
                            selected_validation = active_validation;
                            selected_damping = 1.0;
                            selected_correction_mode =
                                "linearized_active_feasibility_repair";
                        }
                    }
                    output.fixed_jacobian_predictor_trace.back()[
                        "active_feasibility_repair"] =
                        std::move(active_repair_json);
                    output.fixed_jacobian_predictor_trace.back()[
                        "selected_damping"] = selected_damping;
                    output.fixed_jacobian_predictor_trace.back()[
                        "selected_correction_mode"] =
                        selected_correction_mode;
                    output.fixed_jacobian_predictor_trace.back()[
                        "selected_next_validation"] =
                        selected_validation.to_json();
                }
                if (selected_damping == 0.0) {
                    break;
                }
                predictor_state = std::move(selected_correction);
            }
        }
    }

    // The first pass of the two-stage large-case scheduler is only a
    // classifier.  Once the reusable fixed-Jacobian predictor has failed
    // independent validation, continuing into the legacy Newton loops spends
    // seconds on a candidate that the exact corrective stage must solve
    // anyway.  Return the best predictor candidate as the explicitly
    // unverified fallback seed; solve_loaded_contingency still records this as
    // a failed screen and never accepts it as a feasible contingency result.
    if (!base_mode && options_.fixed_jacobian_screen_only) {
        output.solve.status = 2;
        output.solve.iterations = output.fixed_jacobian_predictor_iterations;
        output.solve.state = best_intermediate_state;
        output.solve.objective = rebuild_contingency_state_derived_fields(
            data_, base_state_, commitment_, *contingency,
            output.solve.state);
        output.validation = validate_state(
            data_, ModelMode::ContingencySoft, output.solve.state,
            commitment_, direct_context);
        output.best_intermediate_candidate_selected = true;
        output.best_intermediate_candidate_source = best_intermediate_source;
        output.best_intermediate_candidate_validation = output.validation;
        output.feasible = false;
        output.failure_reason =
            "fixed-Jacobian screen requires exact corrective fallback; "
            "independent validation failed: "
            + output.validation.worst_category + " at "
            + output.validation.worst_identity;
        output.wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        output.solve.wall_seconds = output.wall_seconds;
        return output;
    }

    const bool use_predictor_newton_seed =
        best_intermediate_source == "fixed_base_jacobian_predictor" &&
        best_intermediate_validation.max_residual <
            output.direct_candidate_validation.max_residual;
    const AcState& newton_seed = use_predictor_newton_seed
        ? best_intermediate_state : initial_state;
    const AcState& balance_seed = use_predictor_newton_seed
        ? best_intermediate_state : direct_state;
    if (use_predictor_newton_seed) {
        pg = newton_seed.pg;
        qg = newton_seed.qg;
    }
    const YRows ybus = build_ybus(data_, outaged_branch, &newton_seed);
    std::vector<double> vm = newton_seed.vm;
    std::vector<double> va = newton_seed.va;
    for (int i = 0; i < nb; ++i) {
        vm[i] = std::clamp(vm[i], data_.buses[i].vmin, data_.buses[i].vmax);
        if (data_.buses[i].type == 3) {
            va[i] = 0.0;
        }
    }

    std::vector<double> load_p(nb, 0.0), load_q(nb, 0.0);
    for (int i = 0; i < static_cast<int>(data_.loads.size()); ++i) {
        const auto& load = data_.loads[i];
        load_p[load.bus] += load.pd_nominal * newton_seed.demand_factor[i];
        load_q[load.bus] += load.qd_nominal * newton_seed.demand_factor[i];
    }
    std::vector<double> fixed_q_bus(nb, 0.0);
    std::vector<double> active_slack_target(nb, 0.0);
    std::vector<double> reactive_slack_target(nb, 0.0);
    for (int bus = 0; bus < nb; ++bus) {
        double p_balance = 0.0;
        double q_balance = 0.0;
        for (int branch : data_.buses[bus].branches_from) {
            p_balance += balance_seed.pf[branch];
            q_balance += balance_seed.qf[branch];
        }
        for (int branch : data_.buses[bus].branches_to) {
            p_balance += balance_seed.pt[branch];
            q_balance += balance_seed.qt[branch];
        }
        for (int gen : data_.buses[bus].generators) {
            p_balance -= balance_seed.pg[gen];
            q_balance -= balance_seed.qg[gen];
        }
        for (int load : data_.buses[bus].loads) {
                p_balance += data_.loads[load].pd_nominal
                * balance_seed.demand_factor[load];
            q_balance += data_.loads[load].qd_nominal
                * balance_seed.demand_factor[load];
        }
        for (int shunt : data_.buses[bus].shunts) {
            const double vm2 = balance_seed.vm[bus] * balance_seed.vm[bus];
            p_balance += data_.shunts[shunt].gs * vm2;
            q_balance -= effective_shunt_susceptance(
                data_, balance_seed, shunt) * vm2;
        }
        // Both GO2 base and corrective models permit bounded nodal imbalance.
        // Keep the direct candidate's signed balance when it is already inside
        // the source limit and correct only the excess.  The independent full
        // validator still checks the exact source bound before acceptance.
        active_slack_target[bus] = std::clamp(p_balance, -0.49, 0.49);
        reactive_slack_target[bus] = std::clamp(q_balance, -0.49, 0.49);
    }
    for (int bus = 0; bus < nb; ++bus) {
        for (int gen : active_at_bus[bus]) {
            fixed_q_bus[bus] += qg[gen];
        }
    }

    const auto evaluate_newton_candidate = [&](bool clamp_voltage) {
        AcState candidate = newton_seed;
        candidate.vm = vm;
        candidate.va = va;
        candidate.pg = pg;
        candidate.qg = qg;
        if (clamp_voltage) {
            for (int bus = 0; bus < nb; ++bus) {
                candidate.vm[bus] = std::clamp(
                    candidate.vm[bus],
                    data_.buses[bus].vmin, data_.buses[bus].vmax);
            }
        }
        compute_branch_flows(data_, outaged_branch, !base_mode, candidate);
        for (int bus = 0; bus < nb; ++bus) {
            if (active_at_bus[bus].empty()) {
                continue;
            }
            double p_balance = 0.0;
            double q_balance = 0.0;
            for (int branch : data_.buses[bus].branches_from) {
                p_balance += candidate.pf[branch];
                q_balance += candidate.qf[branch];
            }
            for (int branch : data_.buses[bus].branches_to) {
                p_balance += candidate.pt[branch];
                q_balance += candidate.qt[branch];
            }
            double current_pg = 0.0;
            double current_qg = 0.0;
            for (int gen : data_.buses[bus].generators) {
                p_balance -= candidate.pg[gen];
                q_balance -= candidate.qg[gen];
            }
            for (int gen : active_at_bus[bus]) {
                current_pg += candidate.pg[gen];
                current_qg += candidate.qg[gen];
            }
            for (int load : data_.buses[bus].loads) {
                p_balance += data_.loads[load].pd_nominal
                    * candidate.demand_factor[load];
                q_balance += data_.loads[load].qd_nominal
                    * candidate.demand_factor[load];
            }
            for (int shunt : data_.buses[bus].shunts) {
                const double vm2 = candidate.vm[bus] * candidate.vm[bus];
                p_balance += data_.shunts[shunt].gs * vm2;
                q_balance -= effective_shunt_susceptance(
                    data_, candidate, shunt) * vm2;
            }
            const double p_target = std::clamp(p_balance, -0.49, 0.49);
            const double q_target = std::clamp(q_balance, -0.49, 0.49);
            auto proposed_pg = candidate.pg;
            if (allocate_total(
                    active_at_bus[bus], p_lower, p_upper,
                    initial_state.pg,
                    current_pg + p_balance - p_target,
                    proposed_pg)) {
                candidate.pg = std::move(proposed_pg);
            }
            auto proposed_qg = candidate.qg;
            if (allocate_total(
                    active_at_bus[bus], q_lower, q_upper,
                    initial_state.qg,
                    current_qg + q_balance - q_target,
                    proposed_qg)) {
                candidate.qg = std::move(proposed_qg);
            }
        }
        const auto balance = nodal_balance_slack_seed(
            data_, candidate, 0.5, 1e-7);
        candidate.p_delta = balance.active;
        candidate.q_delta = balance.reactive;
        candidate.gen_lambda.clear();
        for (int i = 0; i < ng; ++i) {
            if (!active[i]) {
                continue;
            }
            const auto points = active_pwl_points(
                data_.generators[i].cost, data_.generators[i].ncost,
                p_lower[i], p_upper[i]);
            append_weights(candidate.gen_lambda, points, candidate.pg[i]);
        }
        candidate.load_lambda.clear();
        for (int i = 0; i < static_cast<int>(data_.loads.size()); ++i) {
            const auto& load = data_.loads[i];
            const auto points = active_pwl_points(
                load.cost, load.ncost, load.pd_min, load.pd_max);
            append_weights(
                candidate.load_lambda, points,
                load.pd_nominal * candidate.demand_factor[i]);
        }
        auto validation = validate_state(
            data_,
            base_mode ? ModelMode::BaseSoft : ModelMode::ContingencySoft,
            candidate, commitment_, direct_context);
        return std::pair<AcState, ValidationReport>{
            std::move(candidate), std::move(validation)};
    };

    if (!base_mode) {
        const auto initial_vm = vm;
        const auto initial_va = va;
        std::vector<bool> active_only_pq(nb, false);
        std::vector<double> active_only_p_spec(nb, 0.0);
        std::vector<double> unused_q_spec(nb, 0.0);
        for (int bus = 0; bus < nb; ++bus) {
            active_only_p_spec[bus] =
                -load_p[bus] + active_slack_target[bus];
            for (int gen : active_at_bus[bus]) {
                active_only_p_spec[bus] += pg[gen];
            }
        }
        output.active_only_newton_attempted = true;
        const auto active_only = run_newton(
            data_, ybus, slack, active_only_pq,
            active_only_p_spec, unused_q_spec,
            options_.max_newton_iterations, options_.newton_tolerance,
            vm, va);
        output.active_only_newton_converged = active_only.converged;
        output.active_only_newton_iterations = active_only.iterations;
        output.newton_iterations += active_only.iterations;
        auto [active_only_state, active_only_validation] =
            evaluate_newton_candidate(false);
        output.active_only_newton_validation = active_only_validation;
        retain_best_intermediate(
            active_only_state, active_only_validation,
            "active_only_newton");
        if (active_only_validation.max_residual <=
            options_.validation_tolerance) {
            output.converged = active_only.converged;
            output.feasible = true;
            output.active_only_newton_selected = true;
            output.solve.status = active_only.converged ? 0 : 1;
            output.solve.objective = 0.0;
            output.solve.iterations = output.newton_iterations;
            output.solve.state = std::move(active_only_state);
            output.validation = active_only_validation;
            output.wall_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - wall_start).count();
            output.solve.wall_seconds = output.wall_seconds;
            return output;
        }
        if (active_only_validation.max_residual <
            output.direct_candidate_validation.max_residual) {
            vm = active_only_state.vm;
            va = active_only_state.va;
            pg = active_only_state.pg;
            qg = active_only_state.qg;
            std::fill(fixed_q_bus.begin(), fixed_q_bus.end(), 0.0);
            for (int bus = 0; bus < nb; ++bus) {
                for (int gen : active_at_bus[bus]) {
                    fixed_q_bus[bus] += qg[gen];
                }
            }
        } else {
            vm = initial_vm;
            va = initial_va;
        }

        const auto reactive_initial_vm = vm;
        const auto reactive_initial_va = va;
        std::vector<double> current_p_network;
        std::vector<double> current_q_network;
        network_injections(
            ybus, vm, va, current_p_network, current_q_network);
        std::vector<double> unused_p_spec(nb, 0.0);
        std::vector<double> reactive_q_spec(nb, 0.0);
        for (int bus = 0; bus < nb; ++bus) {
            double generated_q = 0.0;
            for (int gen : active_at_bus[bus]) {
                generated_q += qg[gen];
            }
            const double q_balance = current_q_network[bus]
                - generated_q + load_q[bus];
            const double q_target = std::clamp(q_balance, -0.49, 0.49);
            reactive_q_spec[bus] = generated_q - load_q[bus] + q_target;
        }
        std::vector<bool> fixed_angle(nb, true);
        output.reactive_only_newton_attempted = true;
        const auto reactive_only = run_newton(
            data_, ybus, fixed_angle, pq,
            unused_p_spec, reactive_q_spec,
            options_.max_newton_iterations, options_.newton_tolerance,
            vm, va);
        output.reactive_only_newton_converged = reactive_only.converged;
        output.reactive_only_newton_iterations = reactive_only.iterations;
        output.newton_iterations += reactive_only.iterations;
        auto [reactive_state, reactive_validation] =
            evaluate_newton_candidate(false);
        if (reactive_validation.worst_category == "variable_bound") {
            auto [clamped_state, clamped_validation] =
                evaluate_newton_candidate(true);
            if (clamped_validation.max_residual <
                reactive_validation.max_residual) {
                reactive_state = std::move(clamped_state);
                reactive_validation = clamped_validation;
            }
        }
        for (int scaling_pass = 0;
             scaling_pass < 6 &&
             reactive_validation.max_residual > options_.validation_tolerance &&
             (reactive_validation.worst_category == "variable_bound" ||
              reactive_validation.worst_category == "flow_limit");
             ++scaling_pass) {
            double maximum_component_ratio = 1.0;
            for (int i = 0; i < static_cast<int>(data_.branches.size()); ++i) {
                const double rating = branch_rating(i);
                if (i == outaged_branch || data_.branches[i].status == 0 ||
                    rating <= 1e-12) {
                    continue;
                }
                maximum_component_ratio = std::max(
                    maximum_component_ratio,
                    std::max({std::abs(reactive_state.pf[i]),
                              std::abs(reactive_state.qf[i]),
                              std::abs(reactive_state.pt[i]),
                              std::abs(reactive_state.qt[i])})
                        / rating);
            }
            if (maximum_component_ratio <= 1.0 + 1e-12) {
                break;
            }
            const double scale = std::min(
                0.9995, 0.9999 / std::sqrt(maximum_component_ratio));
            vm = reactive_state.vm;
            va = reactive_state.va;
            pg = reactive_state.pg;
            qg = reactive_state.qg;
            bool changed = false;
            for (int bus = 0; bus < nb; ++bus) {
                const double reduced = std::max(
                    data_.buses[bus].vmin, vm[bus] * scale);
                changed = changed || reduced < vm[bus] - 1e-12;
                vm[bus] = reduced;
            }
            if (!changed) {
                break;
            }
            auto [scaled_state, scaled_validation] =
                evaluate_newton_candidate(false);
            if (scaled_validation.max_residual + 1e-12 <
                reactive_validation.max_residual) {
                reactive_state = std::move(scaled_state);
                reactive_validation = scaled_validation;
            } else {
                break;
            }
        }
        output.reactive_only_newton_validation = reactive_validation;
        retain_best_intermediate(
            reactive_state, reactive_validation,
            "reactive_only_newton");
        if (reactive_validation.max_residual <=
            options_.validation_tolerance) {
            output.converged = reactive_only.converged;
            output.feasible = true;
            output.reactive_only_newton_selected = true;
            output.solve.status = reactive_only.converged ? 0 : 1;
            output.solve.objective = 0.0;
            output.solve.iterations = output.newton_iterations;
            output.solve.state = std::move(reactive_state);
            output.validation = reactive_validation;
            output.wall_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - wall_start).count();
            output.solve.wall_seconds = output.wall_seconds;
            return output;
        }
        if (reactive_validation.max_residual <
            output.active_only_newton_validation.max_residual) {
            vm = reactive_state.vm;
            va = reactive_state.va;
            pg = reactive_state.pg;
            qg = reactive_state.qg;
            std::fill(fixed_q_bus.begin(), fixed_q_bus.end(), 0.0);
            for (int bus = 0; bus < nb; ++bus) {
                for (int gen : active_at_bus[bus]) {
                    fixed_q_bus[bus] += qg[gen];
                }
            }
        } else {
            vm = reactive_initial_vm;
            va = reactive_initial_va;
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
                p_spec[bus] = -load_p[bus] + active_slack_target[bus];
                q_spec[bus] = fixed_q_bus[bus] - load_q[bus]
                    + reactive_slack_target[bus];
                for (int gen : active_at_bus[bus]) {
                    p_spec[bus] += pg[gen];
                }
            }
            if (active_pass == 0 && q_pass == 0) {
                std::vector<double> initial_p;
                std::vector<double> initial_q;
                network_injections(ybus, vm, va, initial_p, initial_q);
                std::vector<int> angle_index(nb, -1);
                std::vector<int> voltage_index(nb, -1);
                int angle_count = 0;
                int voltage_count = 0;
                for (int bus = 0; bus < nb; ++bus) {
                    if (!slack[bus]) {
                        angle_index[bus] = angle_count++;
                    }
                    if (pq[bus]) {
                        voltage_index[bus] = voltage_count++;
                    }
                }
                output.initial_newton_residual = mismatch_norm(
                    p_spec, q_spec, initial_p, initial_q,
                    angle_index, voltage_index);
            }
            const NewtonResult newton = run_newton(
                data_, ybus, slack, pq, p_spec, q_spec,
                options_.max_newton_iterations, options_.newton_tolerance,
                vm, va);
            output.newton_iterations += newton.iterations;
            auto [newton_state, newton_validation] =
                evaluate_newton_candidate(false);
            if (newton_validation.worst_category == "variable_bound") {
                auto [clamped_state, clamped_validation] =
                    evaluate_newton_candidate(true);
                if (clamped_validation.max_residual <
                    newton_validation.max_residual) {
                    newton_state = std::move(clamped_state);
                    newton_validation = clamped_validation;
                }
            }
            output.newton_candidate_validation = newton_validation;
            retain_best_intermediate(
                newton_state, newton_validation,
                "full_newton");
            if (newton_validation.max_residual <=
                options_.validation_tolerance) {
                output.converged = newton.converged;
                output.feasible = true;
                output.newton_candidate_selected = true;
                output.solve.status = newton.converged ? 0 : 1;
                output.solve.objective = 0.0;
                output.solve.iterations = output.newton_iterations;
                output.solve.state = std::move(newton_state);
                output.validation = newton_validation;
                output.wall_seconds = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - wall_start).count();
                output.solve.wall_seconds = output.wall_seconds;
                return output;
            }
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
                    const double target = vm[bus] < data_.buses[bus].vmin
                        ? 0.49 : -0.49;
                    vm[bus] = projected;
                    if (std::abs(reactive_slack_target[bus] - target) > 1e-12) {
                        // Base and corrective models both permit a bounded
                        // local reactive-balance slack.  Keep this bus in the
                        // Newton Q equations and solve to an interior 0.49-p.u.
                        // target before treating its voltage as fixed.  This
                        // avoids concentrating an arbitrary Q mismatch at a
                        // projected voltage-bound bus.
                        reactive_slack_target[bus] = target;
                    } else {
                        pq[bus] = false;
                    }
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
                const double required = q_network[bus] + load_q[bus]
                    - reactive_slack_target[bus];
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
            const double required_slack = p_network[slack_bus]
                + load_p[slack_bus] - active_slack_target[slack_bus];
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
                required_p = p_network[bus] + load_p[bus]
                    - active_slack_target[bus];
            }
            const double required_q = pq[bus]
                ? fixed_q_bus[bus]
                : q_network[bus] + load_q[bus]
                    - reactive_slack_target[bus];
            auto proposed_pg = pg;
            auto proposed_qg = qg;
            if (!allocate_total(
                    active_at_bus[bus], p_lower, p_upper,
                    initial_state.pg, required_p, proposed_pg) ||
                !allocate_total(
                    active_at_bus[bus], q_lower, q_upper,
                    initial_state.qg, required_q, proposed_qg)) {
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

    AcState state = initial_state;
    state.vm = std::move(vm);
    state.va = std::move(va);
    state.pg = std::move(pg);
    state.qg = std::move(qg);
    compute_branch_flows(data_, outaged_branch, !base_mode, state);
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
            q_balance -= effective_shunt_susceptance(
                data_, state, shunt) * state.vm[bus] * state.vm[bus];
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
        const double previous = base_mode
            ? load.pd_prev
            : load.pd_nominal * base_state_.demand_factor[i];
        const double factor_lower = std::abs(load.pd_nominal) <= 1e-12
            ? load.tmin
            : std::max(load.tmin,
                (previous - (base_mode ? load.prdmax * data_.delta_r
                                       : load.prdmaxctg * data_.delta_r_ctg))
                    / load.pd_nominal);
        const double factor_upper = std::abs(load.pd_nominal) <= 1e-12
            ? load.tmax
            : std::min(load.tmax,
                (previous + (base_mode ? load.prumax * data_.delta_r
                                       : load.prumaxctg * data_.delta_r_ctg))
                    / load.pd_nominal);
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
            load.qd_nominal * initial_state.demand_factor[i];
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
        compute_branch_flows(data_, outaged_branch, !base_mode, state);
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
                    q_balance[bus] -= effective_shunt_susceptance(
                        data_, state, shunt)
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
                        initial_state.qg, current + q_excess, proposed)) {
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
                    initial_state.pg, current + p_excess, proposed)) {
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

    std::optional<ContingencyContext> context;
    if (!base_mode) {
        context = ContingencyContext{};
        context->base_state = base_state_;
        context->outaged_generator = outaged_generator;
        context->outaged_branch = outaged_branch;
    }
    const auto validate_candidate = [&](const AcState& candidate) {
        return validate_state(
            data_,
            base_mode ? ModelMode::BaseSoft : ModelMode::ContingencySoft,
            candidate, commitment_, context);
    };
    output.solve.status = 0;
    output.solve.objective = 0.0;
    output.solve.iterations = output.newton_iterations;
    output.solve.state = std::move(state);
    output.converged = !approximate_candidate;
    output.validation = validate_candidate(output.solve.state);
    for (int scaling_pass = 0;
         scaling_pass < 6 && output.validation.max_residual > options_.validation_tolerance;
         ++scaling_pass) {
        if (output.validation.worst_category != "variable_bound" &&
            output.validation.worst_category != "flow_limit") {
            break;
        }
        double maximum_component_ratio = 1.0;
        for (int i = 0; i < static_cast<int>(data_.branches.size()); ++i) {
            const double rating = branch_rating(i);
            if (i == outaged_branch || data_.branches[i].status == 0 ||
                rating <= 1e-12) {
                continue;
            }
            maximum_component_ratio = std::max(maximum_component_ratio,
                std::max({std::abs(output.solve.state.pf[i]),
                          std::abs(output.solve.state.qf[i]),
                          std::abs(output.solve.state.pt[i]),
                          std::abs(output.solve.state.qt[i])})
                    / rating);
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
        output.validation = validate_candidate(output.solve.state);
    }

    // Voltage projection used to clear branch-component bounds can disturb
    // nodal balance.  A conventional single-slack power flow concentrates the
    // resulting active mismatch at one reference bus and can exceed the
    // model's per-bus soft-imbalance bound.  Re-solve the AC equations with a
    // distributed active slack across buses in each connected component. The
    // corrective model explicitly permits up to 0.5 p.u. of active imbalance
    // at each bus. Select bounded imbalance targets near the current state and
    // solve one scalar correction per component so nonlinear losses remain
    // consistent, rather than forcing the entire loss mismatch onto one bus.
    // The original state remains the incumbent and every candidate must
    // improve independent validation.
    if (options_.distributed_balance_polish &&
        output.validation.max_residual > options_.validation_tolerance &&
        (output.validation.worst_category == "active_balance" ||
         output.validation.worst_category == "reactive_balance")) {
        output.distributed_balance_polish_attempted = true;
        const AcState original_state = output.solve.state;
        const ValidationReport original_validation = output.validation;
        AcState best_state = original_state;
        ValidationReport best_validation = original_validation;

        std::vector<bool> polish_pq = pq;

        state = original_state;
        for (int repair_pass = 0; repair_pass < 8; ++repair_pass) {
            std::vector<double> p_spec(nb, 0.0);
            std::vector<double> q_spec(nb, 0.0);
            for (int bus = 0; bus < nb; ++bus) {
                for (int gen : data_.buses[bus].generators) {
                    p_spec[bus] += state.pg[gen];
                    q_spec[bus] += state.qg[gen];
                }
                for (int load : data_.buses[bus].loads) {
                    p_spec[bus] -= data_.loads[load].pd_nominal
                        * state.demand_factor[load];
                    q_spec[bus] -= data_.loads[load].qd_nominal
                        * state.demand_factor[load];
                }
            }
            std::vector<double> p_network;
            std::vector<double> q_network;
            network_injections(
                ybus, state.vm, state.va, p_network, q_network);
            std::vector<double> active_slack_weights(nb, 0.0);
            bool weights_valid = true;
            for (const auto& component : components) {
                std::vector<int> participant_buses = component;
                double total_balance = 0.0;
                for (int bus : component) {
                    total_balance += p_network[bus] - p_spec[bus];
                }
                if (participant_buses.empty()) {
                    weights_valid = false;
                    break;
                }
                if (std::abs(total_balance) <= 1e-10) {
                    const double weight = 1.0 / participant_buses.size();
                    for (int bus : participant_buses) {
                        active_slack_weights[bus] = weight;
                    }
                } else {
                    std::vector<double> target_lower(nb, 0.0);
                    std::vector<double> target_upper(nb, 0.0);
                    std::vector<double> target_preferred(nb, 0.0);
                    std::vector<double> target_values(nb, 0.0);
                    for (int bus : participant_buses) {
                        target_lower[bus] = -0.49;
                        target_upper[bus] = 0.49;
                        target_preferred[bus] = std::clamp(
                            p_network[bus] - p_spec[bus], -0.49, 0.49);
                    }
                    if (!allocate_total(
                            participant_buses,
                            target_lower, target_upper, target_preferred,
                            total_balance, target_values)) {
                        weights_valid = false;
                        break;
                    }
                    for (int bus : participant_buses) {
                        active_slack_weights[bus] = target_values[bus];
                    }
                }
            }
            if (!weights_valid) {
                output.distributed_balance_polish_failure_reason =
                    "cannot allocate bounded distributed active slack";
                break;
            }
            const NewtonResult polish = run_distributed_active_newton(
                data_, ybus, slack, polish_pq, component_of,
                active_slack_weights, p_spec, q_spec,
                options_.max_newton_iterations, options_.newton_tolerance,
                state.vm, state.va);
            output.distributed_balance_polish_iterations += polish.iterations;
            if (!polish.converged) {
                output.distributed_balance_polish_failure_reason =
                    polish.failure_reason;
                break;
            }

            bool projected_voltage = false;
            for (int bus = 0; bus < nb; ++bus) {
                if (!polish_pq[bus]) {
                    continue;
                }
                const double projected = std::clamp(
                    state.vm[bus],
                    data_.buses[bus].vmin, data_.buses[bus].vmax);
                if (std::abs(projected - state.vm[bus]) > 1e-9) {
                    state.vm[bus] = projected;
                    polish_pq[bus] = false;
                    projected_voltage = true;
                    ++output.distributed_balance_voltage_projections;
                }
            }
            if (projected_voltage) {
                continue;
            }

            bool dispatch_valid = true;
            network_injections(ybus, state.vm, state.va, p_network, q_network);
            for (int bus = 0; bus < nb; ++bus) {
                if (active_at_bus[bus].empty() || polish_pq[bus]) {
                    continue;
                }
                const double required_q = q_network[bus] + load_q[bus];
                double lower = 0.0;
                double upper = 0.0;
                for (int gen : active_at_bus[bus]) {
                    lower += q_lower[gen];
                    upper += q_upper[gen];
                }
                auto proposed = state.qg;
                if (!allocate_total(
                        active_at_bus[bus], q_lower, q_upper,
                        state.qg, std::clamp(required_q, lower, upper),
                        proposed)) {
                    dispatch_valid = false;
                    break;
                }
                state.qg = std::move(proposed);
            }
            if (!dispatch_valid) {
                output.distributed_balance_polish_failure_reason =
                    "distributed reactive dispatch exceeds generator bounds";
                break;
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
                    const double offset = state.va[angle_reference];
                    for (int bus : component) {
                        state.va[bus] -= offset;
                    }
                }
            }
            refresh_network_fields();
            const ValidationReport polished_validation =
                validate_candidate(state);
            output.distributed_balance_polish_validation = polished_validation;
            if (polished_validation.max_residual + 1e-12 <
                best_validation.max_residual) {
                best_state = state;
                best_validation = polished_validation;
            }
            if (polished_validation.max_residual <=
                options_.validation_tolerance) {
                break;
            }
            if (polished_validation.worst_category == "reactive_balance") {
                network_injections(
                    ybus, state.vm, state.va, p_network, q_network);
                bool released_reactive_bus = false;
                for (int bus = 0; bus < nb; ++bus) {
                    if (active_at_bus[bus].empty() || polish_pq[bus]) {
                        continue;
                    }
                    double generated_q = 0.0;
                    for (int gen : active_at_bus[bus]) {
                        generated_q += state.qg[gen];
                    }
                    const double q_balance =
                        q_network[bus] + load_q[bus] - generated_q;
                    if (std::abs(q_balance) <= 0.45) {
                        continue;
                    }
                    const auto& source_bus = data_.buses[bus];
                    if (state.vm[bus] <= source_bus.vmin + 1e-7 ||
                        state.vm[bus] >= source_bus.vmax - 1e-7) {
                        continue;
                    }
                    polish_pq[bus] = true;
                    released_reactive_bus = true;
                }
                if (released_reactive_bus) {
                    continue;
                }
                break;
            }
            if (polished_validation.worst_category != "variable_bound" &&
                polished_validation.worst_category != "flow_limit") {
                break;
            }

            double maximum_component_ratio = 1.0;
            for (int i = 0; i < static_cast<int>(data_.branches.size()); ++i) {
                const double rating = branch_rating(i);
                if (i == outaged_branch || data_.branches[i].status == 0 ||
                    rating <= 1e-12) {
                    continue;
                }
                maximum_component_ratio = std::max(maximum_component_ratio,
                    std::max({std::abs(state.pf[i]), std::abs(state.qf[i]),
                              std::abs(state.pt[i]), std::abs(state.qt[i])})
                        / rating);
            }
            const double scale = std::min(
                0.9995, 0.9999 / std::sqrt(maximum_component_ratio));
            bool changed = false;
            for (int bus = 0; bus < nb; ++bus) {
                const double reduced = std::max(
                    data_.buses[bus].vmin, state.vm[bus] * scale);
                changed = changed || reduced < state.vm[bus] - 1e-12;
                state.vm[bus] = reduced;
            }
            if (!changed) {
                break;
            }
            refresh_network_fields();
        }

        output.newton_iterations += output.distributed_balance_polish_iterations;
        if (best_validation.max_residual + 1e-12 <
            original_validation.max_residual) {
            output.solve.state = std::move(best_state);
            output.validation = best_validation;
            output.distributed_balance_polish_selected = true;
        } else {
            state = original_state;
        }
    }
    output.best_intermediate_candidate_validation =
        best_intermediate_validation;
    if (best_intermediate_validation.max_residual + 1e-12 <
        output.validation.max_residual) {
        output.solve.state = std::move(best_intermediate_state);
        output.validation = best_intermediate_validation;
        output.best_intermediate_candidate_selected = true;
        output.best_intermediate_candidate_source =
            best_intermediate_source;
    }
    output.feasible = output.validation.max_residual <= options_.validation_tolerance;
    if (!output.feasible) {
        const std::string validation_failure =
            "independent validation failed: "
            + output.validation.worst_category + " at "
            + output.validation.worst_identity;
        output.failure_reason = output.failure_reason.empty()
            ? validation_failure
            : output.failure_reason + "; " + validation_failure;
    } else if (approximate_candidate) {
        output.failure_reason.clear();
    }
    output.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_start).count();
    output.solve.wall_seconds = output.wall_seconds;
    return output;
}

}  // namespace gravityx
