#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include "gravityx/active_feasibility_repair.hpp"
#include "gravityx/fast_power_flow.hpp"
#include "gravityx/linearized_ac_seed.hpp"
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
    // Callers frequently allocate the generators or loads at one bus while
    // the full case contains thousands of controls.  Keep this operation
    // transactional over only the touched indices so callers do not need to
    // copy the entire case-wide vector before every local projection.
    std::vector<double> original;
    original.reserve(indices.size());
    for (int index : indices) {
        original.push_back(values[index]);
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
    if (std::abs(target - current) <= 1e-7) {
        return true;
    }
    for (std::size_t position = 0; position < indices.size(); ++position) {
        values[indices[position]] = original[position];
    }
    return false;
}

bool allocate_total_in_order(
    const std::vector<int>& ordered_indices,
    const std::vector<double>& lower,
    const std::vector<double>& upper,
    const std::vector<double>& preferred,
    double target,
    std::vector<double>& values) {
    double total_lower = 0.0;
    double total_upper = 0.0;
    for (int index : ordered_indices) {
        total_lower += lower[index];
        total_upper += upper[index];
    }
    if (target < total_lower - kAllocationTolerance ||
        target > total_upper + kAllocationTolerance) {
        return false;
    }
    std::vector<double> original;
    original.reserve(ordered_indices.size());
    double current = 0.0;
    for (int index : ordered_indices) {
        original.push_back(values[index]);
        values[index] = std::clamp(
            preferred[index], lower[index], upper[index]);
        current += values[index];
    }
    double difference = target - current;
    for (int index : ordered_indices) {
        if (std::abs(difference) <= kAllocationTolerance) {
            return true;
        }
        const double adjustment = difference > 0.0
            ? std::min(difference, upper[index] - values[index])
            : -std::min(-difference, values[index] - lower[index]);
        values[index] += adjustment;
        difference -= adjustment;
    }
    if (std::abs(difference) <= 1e-7) {
        return true;
    }
    for (std::size_t position = 0;
         position < ordered_indices.size(); ++position) {
        values[ordered_indices[position]] = original[position];
    }
    return false;
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

double pwl_value_at(
    const std::vector<PwlPoint>& points,
    double target) {
    const auto weights = interpolation_weights(points, target);
    double value = 0.0;
    for (std::size_t point = 0; point < points.size(); ++point) {
        value += weights[point] * points[point].cost;
    }
    return value;
}

struct EconomicMeritTarget {
    bool feasible{};
    bool curve_order_valid{true};
    std::vector<double> pg;
    std::vector<double> load_power;
    double initial_market_surplus{};
    double target_market_surplus{};
    double generation_movement{};
    double load_movement{};
    int component_count{};
};

// Solve the network-free portion of the corrective PWL economics exactly.
// For convex generator costs and concave load benefits, the component problem
//
//   max  benefit(D) - cost(G)       subject to G - D = fixed net injection
//
// is a one-dimensional convex resource-allocation problem.  Starting from the
// minimum feasible demand, generator segments are consumed from cheapest to
// dearest and load segments from highest to lowest benefit.  Paired increments
// continue while marginal benefit exceeds marginal cost.  The resulting
// controls are only a target; nonlinear AC repair and the complete validator
// decide whether any part of it may be accepted.
EconomicMeritTarget build_component_economic_merit_target(
    const CaseData& data,
    const std::vector<std::vector<int>>& components,
    const std::vector<bool>& active,
    const std::vector<double>& generator_lower,
    const std::vector<double>& generator_upper,
    const std::vector<double>& load_lower,
    const std::vector<double>& load_upper,
    const std::vector<double>& reference_pg,
    const std::vector<double>& reference_load_power) {
    EconomicMeritTarget output;
    const int nb = static_cast<int>(data.buses.size());
    const int ng = static_cast<int>(data.generators.size());
    const int nd = static_cast<int>(data.loads.size());
    if (active.size() != data.generators.size() ||
        generator_lower.size() != data.generators.size() ||
        generator_upper.size() != data.generators.size() ||
        load_lower.size() != data.loads.size() ||
        load_upper.size() != data.loads.size() ||
        reference_pg.size() != data.generators.size() ||
        reference_load_power.size() != data.loads.size()) {
        return output;
    }

    std::vector<int> component_of_bus(static_cast<std::size_t>(nb), -1);
    for (int component = 0;
         component < static_cast<int>(components.size()); ++component) {
        for (int bus : components[component]) {
            if (bus < 0 || bus >= nb || component_of_bus[bus] >= 0) {
                return output;
            }
            component_of_bus[bus] = component;
        }
    }
    if (std::any_of(
            component_of_bus.begin(), component_of_bus.end(),
            [](int component) { return component < 0; })) {
        return output;
    }

    output.pg.assign(data.generators.size(), 0.0);
    output.load_power = reference_load_power;

    struct Segment {
        int device{-1};
        int ordinal{};
        double remaining{};
        double slope{};
    };
    const auto append_segments = [&output](
        const std::vector<PwlPoint>& points,
        int device,
        double lower,
        double upper,
        bool convex,
        std::vector<Segment>& segments) {
        if (!std::isfinite(lower) || !std::isfinite(upper) ||
            lower > upper + 1e-10 || points.size() < 2) {
            output.curve_order_valid = false;
            return;
        }
        double covered = 0.0;
        double prior_slope = convex
            ? -std::numeric_limits<double>::infinity()
            : std::numeric_limits<double>::infinity();
        int ordinal = 0;
        for (std::size_t point = 1; point < points.size(); ++point) {
            const double x0 = points[point - 1].mw;
            const double x1 = points[point].mw;
            const double width = x1 - x0;
            if (!std::isfinite(width) || width <= 1e-14) {
                output.curve_order_valid = false;
                return;
            }
            const double slope =
                (points[point].cost - points[point - 1].cost) / width;
            if (!std::isfinite(slope) ||
                (convex && slope + 1e-8 < prior_slope) ||
                (!convex && slope > prior_slope + 1e-8)) {
                output.curve_order_valid = false;
                return;
            }
            prior_slope = slope;
            const double left = std::max(lower, x0);
            const double right = std::min(upper, x1);
            if (right > left + 1e-14) {
                segments.push_back({device, ordinal++, right - left, slope});
                covered += right - left;
            }
        }
        if (std::abs(covered - (upper - lower)) >
            1e-7 * std::max(1.0, std::abs(upper - lower))) {
            output.curve_order_valid = false;
        }
    };

    std::vector<std::vector<int>> generators_by_component(components.size());
    std::vector<std::vector<int>> loads_by_component(components.size());
    for (int generator = 0; generator < ng; ++generator) {
        if (!active[generator]) {
            continue;
        }
        const int component =
            component_of_bus[data.generators[generator].bus];
        generators_by_component[component].push_back(generator);
    }
    for (int load = 0; load < nd; ++load) {
        const int component = component_of_bus[data.loads[load].bus];
        loads_by_component[component].push_back(load);
    }

    for (int component = 0;
         component < static_cast<int>(components.size()); ++component) {
        std::vector<Segment> supply;
        std::vector<Segment> demand;
        double generation_lower_total = 0.0;
        double generation_upper_total = 0.0;
        double load_lower_total = 0.0;
        double load_upper_total = 0.0;
        double reference_net_injection = 0.0;

        for (int generator : generators_by_component[component]) {
            const double lower = generator_lower[generator];
            const double upper = generator_upper[generator];
            if (lower > upper + 1e-10) {
                return output;
            }
            output.pg[generator] = lower;
            generation_lower_total += lower;
            generation_upper_total += upper;
            reference_net_injection += reference_pg[generator];
            const auto& source = data.generators[generator];
            append_segments(
                active_pwl_points(
                    source.cost, source.ncost, lower, upper),
                generator, lower, upper, true, supply);
        }
        for (int load : loads_by_component[component]) {
            const double lower = load_lower[load];
            const double upper = load_upper[load];
            if (lower > upper + 1e-10) {
                return output;
            }
            output.load_power[load] = lower;
            load_lower_total += lower;
            load_upper_total += upper;
            reference_net_injection -= reference_load_power[load];
            const auto& source = data.loads[load];
            append_segments(
                active_pwl_points(
                    source.cost, source.ncost,
                    source.pd_min, source.pd_max),
                load, lower, upper, false, demand);
        }
        if (!output.curve_order_valid) {
            return output;
        }

        const double feasible_load_lower = std::max(
            load_lower_total,
            generation_lower_total - reference_net_injection);
        const double feasible_load_upper = std::min(
            load_upper_total,
            generation_upper_total - reference_net_injection);
        if (feasible_load_lower > feasible_load_upper + 1e-7) {
            return output;
        }
        const double selected_load_total = std::clamp(
            feasible_load_lower, load_lower_total, load_upper_total);
        const double selected_generation_total =
            selected_load_total + reference_net_injection;

        std::stable_sort(
            supply.begin(), supply.end(),
            [](const Segment& left, const Segment& right) {
                if (std::abs(left.slope - right.slope) > 1e-12) {
                    return left.slope < right.slope;
                }
                if (left.device != right.device) {
                    return left.device < right.device;
                }
                return left.ordinal < right.ordinal;
            });
        std::stable_sort(
            demand.begin(), demand.end(),
            [](const Segment& left, const Segment& right) {
                if (std::abs(left.slope - right.slope) > 1e-12) {
                    return left.slope > right.slope;
                }
                if (left.device != right.device) {
                    return left.device < right.device;
                }
                return left.ordinal < right.ordinal;
            });

        const auto consume = [](std::vector<Segment>& segments,
                                double amount,
                                std::vector<double>& values) {
            for (auto& segment : segments) {
                if (amount <= 1e-10) {
                    break;
                }
                const double used = std::min(amount, segment.remaining);
                values[segment.device] += used;
                segment.remaining -= used;
                amount -= used;
            }
            return amount <= 1e-7;
        };
        if (!consume(
                supply,
                selected_generation_total - generation_lower_total,
                output.pg) ||
            !consume(
                demand,
                selected_load_total - load_lower_total,
                output.load_power)) {
            return output;
        }

        std::size_t supply_position = 0;
        std::size_t demand_position = 0;
        while (true) {
            while (supply_position < supply.size() &&
                   supply[supply_position].remaining <= 1e-12) {
                ++supply_position;
            }
            while (demand_position < demand.size() &&
                   demand[demand_position].remaining <= 1e-12) {
                ++demand_position;
            }
            if (supply_position >= supply.size() ||
                demand_position >= demand.size() ||
                demand[demand_position].slope <=
                    supply[supply_position].slope + 1e-10) {
                break;
            }
            const double paired = std::min(
                supply[supply_position].remaining,
                demand[demand_position].remaining);
            output.pg[supply[supply_position].device] += paired;
            output.load_power[demand[demand_position].device] += paired;
            supply[supply_position].remaining -= paired;
            demand[demand_position].remaining -= paired;
        }

        double selected_net = 0.0;
        for (int generator : generators_by_component[component]) {
            selected_net += output.pg[generator];
        }
        for (int load : loads_by_component[component]) {
            selected_net -= output.load_power[load];
        }
        if (std::abs(selected_net - reference_net_injection) > 1e-7) {
            return output;
        }
    }

    for (int generator = 0; generator < ng; ++generator) {
        if (!active[generator]) {
            continue;
        }
        const auto& source = data.generators[generator];
        const auto points = active_pwl_points(
            source.cost, source.ncost,
            generator_lower[generator], generator_upper[generator]);
        output.initial_market_surplus -=
            pwl_value_at(points, reference_pg[generator]);
        output.target_market_surplus -=
            pwl_value_at(points, output.pg[generator]);
        output.generation_movement +=
            std::abs(output.pg[generator] - reference_pg[generator]);
    }
    for (int load = 0; load < nd; ++load) {
        const auto& source = data.loads[load];
        const auto points = active_pwl_points(
            source.cost, source.ncost,
            source.pd_min, source.pd_max);
        output.initial_market_surplus +=
            pwl_value_at(points, reference_load_power[load]);
        output.target_market_surplus +=
            pwl_value_at(points, output.load_power[load]);
        output.load_movement +=
            std::abs(output.load_power[load] - reference_load_power[load]);
    }
    output.feasible = std::isfinite(output.initial_market_surplus) &&
        std::isfinite(output.target_market_surplus) &&
        output.target_market_surplus + 1e-7 >=
            output.initial_market_surplus;
    return output;
}

std::vector<std::vector<int>> capacity_bounded_network_clusters(
    const CaseData& data,
    int maximum_buses,
    int outaged_branch = -1) {
    const int nb = static_cast<int>(data.buses.size());
    maximum_buses = std::clamp(maximum_buses, 1, std::max(1, nb));
    std::vector<int> parent(static_cast<std::size_t>(nb));
    std::vector<int> size(static_cast<std::size_t>(nb), 1);
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
    std::vector<int> order;
    order.reserve(data.branches.size());
    for (int branch = 0;
         branch < static_cast<int>(data.branches.size()); ++branch) {
        const auto& source = data.branches[branch];
        if (branch != outaged_branch && source.present &&
            source.status != 0 && source.from >= 0 && source.from < nb &&
            source.to >= 0 && source.to < nb) {
            order.push_back(branch);
        }
    }
    std::stable_sort(
        order.begin(), order.end(),
        [&](int left, int right) {
            const double left_capacity = std::max(
                0.0, data.branches[left].rate_c);
            const double right_capacity = std::max(
                0.0, data.branches[right].rate_c);
            if (std::abs(left_capacity - right_capacity) > 1e-12) {
                return left_capacity > right_capacity;
            }
            return left < right;
        });
    for (int branch : order) {
        int left = find_root(data.branches[branch].from);
        int right = find_root(data.branches[branch].to);
        if (left == right || size[left] + size[right] > maximum_buses) {
            continue;
        }
        if (left > right) {
            std::swap(left, right);
        }
        parent[right] = left;
        size[left] += size[right];
    }
    std::vector<int> root_to_component(static_cast<std::size_t>(nb), -1);
    std::vector<std::vector<int>> components;
    components.reserve(static_cast<std::size_t>(nb));
    for (int bus = 0; bus < nb; ++bus) {
        const int root = find_root(bus);
        if (root_to_component[root] < 0) {
            root_to_component[root] = static_cast<int>(components.size());
            components.emplace_back();
        }
        components[root_to_component[root]].push_back(bus);
    }
    return components;
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
        double from_g_self = branch.flow_from_g_self;
        double from_b_self = branch.flow_from_b_self;
        double to_g_self = branch.flow_to_g_self;
        double to_b_self = branch.flow_to_b_self;
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
            const double tap_squared = branch.tap * branch.tap;
            const double tap_real = branch.tap * std::cos(branch.shift);
            const double tap_imag = branch.tap * std::sin(branch.shift);
            from_g_self = branch.transformer
                ? g / tap_squared + branch.g_fr
                : (g + branch.g_fr) / tap_squared;
            from_b_self = branch.transformer
                ? b / tap_squared + branch.b_fr
                : (b + branch.b_fr) / tap_squared;
            to_g_self = g + branch.g_to;
            to_b_self = b + branch.b_to;
            from_cross_cos =
                (-g * tap_real + b * tap_imag) / tap_squared;
            from_cross_sin =
                (-b * tap_real - g * tap_imag) / tap_squared;
            to_cross_cos =
                (-g * tap_real - b * tap_imag) / tap_squared;
            to_cross_sin =
                (-b * tap_real + g * tap_imag) / tap_squared;
        }
        const int f = branch.from;
        const int t = branch.to;
        const double voltage_product = state.vm[f] * state.vm[t];
        const double angle_delta = state.va[f] - state.va[t];
        const double cross_cos_ft = voltage_product * std::cos(angle_delta);
        const double cross_sin_ft = voltage_product * std::sin(angle_delta);
        state.pf[i] = from_g_self * state.vm[f] * state.vm[f]
            + from_cross_cos * cross_cos_ft
            + from_cross_sin * cross_sin_ft;
        state.qf[i] = -from_b_self * state.vm[f] * state.vm[f]
            - from_cross_sin * cross_cos_ft
            + from_cross_cos * cross_sin_ft;
        // cos(-x) == cos(x) and sin(-x) == -sin(x); reuse the same
        // trigonometric pair for the opposite terminal.
        state.pt[i] = to_g_self * state.vm[t] * state.vm[t]
            + to_cross_cos * cross_cos_ft
            - to_cross_sin * cross_sin_ft;
        state.qt[i] = -to_b_self * state.vm[t] * state.vm[t]
            - to_cross_sin * cross_cos_ft
            - to_cross_cos * cross_sin_ft;

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

BalanceSlackSeed nodal_balance_slack_seed_from_network(
    const CaseData& data,
    const AcState& state,
    const std::vector<double>& p_network,
    const std::vector<double>& q_network,
    double upper_bound,
    double interior_margin) {
    const std::size_t nb = data.buses.size();
    if (p_network.size() != nb || q_network.size() != nb ||
        state.pg.size() != data.generators.size() ||
        state.qg.size() != data.generators.size() ||
        state.demand_factor.size() != data.loads.size() ||
        !std::isfinite(upper_bound) || upper_bound < 0.0 ||
        !std::isfinite(interior_margin) || interior_margin < 0.0) {
        throw std::runtime_error(
            "cannot seed nodal slacks from invalid network injections");
    }

    BalanceSlackSeed seed;
    seed.active.resize(nb);
    seed.reactive.resize(nb);
    for (std::size_t i = 0; i < nb; ++i) {
        const auto& bus = data.buses[i];
        double p = p_network[i];
        double q = q_network[i];
        for (int generator : bus.generators) {
            p -= state.pg[generator];
            q -= state.qg[generator];
        }
        for (int load : bus.loads) {
            p += data.loads[load].pd_nominal * state.demand_factor[load];
            q += data.loads[load].qd_nominal * state.demand_factor[load];
        }
        if (!std::isfinite(p) || !std::isfinite(q)) {
            throw std::runtime_error(
                "cannot seed nodal slacks from nonfinite network injections");
        }
        seed.active[i] =
            std::min(upper_bound, std::abs(p) + interior_margin);
        seed.reactive[i] =
            std::min(upper_bound, std::abs(q) + interior_margin);
    }
    return seed;
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

std::vector<unsigned char> bridge_branches(const CaseData& data) {
    struct IncidentEdge {
        int neighbor{-1};
        int branch{-1};
    };
    struct DfsFrame {
        int bus{-1};
        std::size_t next_edge{};
    };

    const int nb = static_cast<int>(data.buses.size());
    std::vector<std::vector<IncidentEdge>> adjacency(
        static_cast<std::size_t>(nb));
    for (int branch_index = 0;
         branch_index < static_cast<int>(data.branches.size());
         ++branch_index) {
        const auto& branch = data.branches[branch_index];
        if (branch.status == 0) {
            continue;
        }
        adjacency[branch.from].push_back({branch.to, branch_index});
        adjacency[branch.to].push_back({branch.from, branch_index});
    }

    std::vector<int> discovery(static_cast<std::size_t>(nb), -1);
    std::vector<int> low(static_cast<std::size_t>(nb), -1);
    std::vector<int> parent_bus(static_cast<std::size_t>(nb), -1);
    std::vector<int> parent_branch(static_cast<std::size_t>(nb), -1);
    std::vector<unsigned char> bridge(data.branches.size(), 0);
    std::vector<DfsFrame> stack;
    stack.reserve(static_cast<std::size_t>(nb));
    int next_discovery = 0;

    for (int root = 0; root < nb; ++root) {
        if (discovery[root] >= 0) {
            continue;
        }
        discovery[root] = next_discovery;
        low[root] = next_discovery;
        ++next_discovery;
        stack.push_back({root, 0});
        while (!stack.empty()) {
            auto& frame = stack.back();
            const int bus = frame.bus;
            if (frame.next_edge < adjacency[bus].size()) {
                const auto edge = adjacency[bus][frame.next_edge++];
                if (edge.branch == parent_branch[bus]) {
                    continue;
                }
                if (discovery[edge.neighbor] < 0) {
                    parent_bus[edge.neighbor] = bus;
                    parent_branch[edge.neighbor] = edge.branch;
                    discovery[edge.neighbor] = next_discovery;
                    low[edge.neighbor] = next_discovery;
                    ++next_discovery;
                    stack.push_back({edge.neighbor, 0});
                    continue;
                }
                low[bus] = std::min(low[bus], discovery[edge.neighbor]);
                continue;
            }

            stack.pop_back();
            const int parent = parent_bus[bus];
            if (parent < 0) {
                continue;
            }
            const int branch = parent_branch[bus];
            if (low[bus] > discovery[parent]) {
                bridge[branch] = 1;
            }
            low[parent] = std::min(low[parent], low[bus]);
        }
    }
    return bridge;
}

void normalize_source_reference_angles(
    const CaseData& data,
    const std::vector<std::vector<int>>& components,
    std::vector<double>& va) {
    if (va.size() != data.buses.size()) {
        throw std::runtime_error(
            "cannot normalize source reference angles with invalid dimensions");
    }
    // The Newton and fixed-Jacobian systems may use any available-generator
    // bus as their numerical slack after a generator outage.  The source
    // model's angle constraint is different: every type-3 bus must remain at
    // zero radians.  A component-wide angle translation leaves all branch
    // angle differences, AC flows, and operating controls unchanged, so
    // restore that source gauge before independently validating a candidate.
    for (const auto& component : components) {
        int source_reference = -1;
        for (int bus : component) {
            if (data.buses[bus].type == 3) {
                source_reference = bus;
                break;
            }
        }
        if (source_reference < 0) {
            continue;
        }
        const double offset = va[source_reference];
        if (offset == 0.0) {
            continue;
        }
        for (int bus : component) {
            va[bus] -= offset;
        }
    }
}

}  // namespace

namespace {

AcState interpolate_corrective_candidate(
    const AcState& from,
    const AcState& to,
    double step) {
    if (!std::isfinite(step) || step < 0.0 || step > 1.0) {
        throw std::runtime_error(
            "corrective candidate interpolation step is outside [0, 1]");
    }
    AcState candidate = from;
    const auto blend = [step](
        const std::vector<double>& left,
        const std::vector<double>& right,
        std::vector<double>& target,
        const char* identity) {
        if (left.size() != right.size()) {
            throw std::runtime_error(
                std::string("corrective candidate dimension mismatch: ") +
                identity);
        }
        target.resize(left.size());
        for (std::size_t i = 0; i < left.size(); ++i) {
            target[i] = left[i] + step * (right[i] - left[i]);
        }
    };
    blend(from.vm, to.vm, candidate.vm, "vm");
    blend(from.va, to.va, candidate.va, "va");
    blend(from.pg, to.pg, candidate.pg, "pg");
    blend(from.qg, to.qg, candidate.qg, "qg");
    blend(
        from.demand_factor, to.demand_factor,
        candidate.demand_factor, "demand_factor");
    // Switched-shunt controls are discrete.  The AC equations for `to` were
    // built with this exact shunt state, so interpolation must transfer it as
    // a unit instead of retaining an unrelated state from `from`.
    candidate.shunt_bs = to.shunt_bs;
    candidate.shunt_steps = to.shunt_steps;
    return candidate;
}

}  // namespace

void run_fast_power_flow_topology_cache_regression() {
    CaseData data;
    data.buses.resize(3);
    data.buses[0].type = 3;
    data.buses[1].type = 1;
    data.buses[2].type = 3;
    data.branches.resize(3);
    for (int branch = 0; branch < 2; ++branch) {
        data.branches[branch].from = 0;
        data.branches[branch].to = 1;
        data.branches[branch].status = 1;
    }
    data.branches[2].from = 1;
    data.branches[2].to = 2;
    data.branches[2].status = 1;

    const auto bridges = bridge_branches(data);
    if (bridges.size() != 3 || bridges[0] != 0 || bridges[1] != 0 ||
        bridges[2] == 0) {
        throw std::runtime_error(
            "fast power-flow topology cache misclassified a parallel edge "
            "or bridge");
    }
    const auto parallel_outage_components = connected_components(data, 0);
    if (parallel_outage_components.size() != 1 ||
        parallel_outage_components.front().size() != 3) {
        throw std::runtime_error(
            "fast power-flow topology cache split a parallel outage");
    }
    const auto bridge_outage_components = connected_components(data, 2);
    if (bridge_outage_components.size() != 2) {
        throw std::runtime_error(
            "fast power-flow topology cache did not split a bridge outage");
    }
    std::vector<double> angles{0.25, 0.40, -0.30};
    normalize_source_reference_angles(
        data, bridge_outage_components, angles);
    if (std::abs(angles[0]) > 1e-12 ||
        std::abs(angles[1] - 0.15) > 1e-12 ||
        std::abs(angles[2]) > 1e-12) {
        throw std::runtime_error(
            "fast power-flow topology cache normalized the wrong island");
    }

    AcState interpolation_from;
    interpolation_from.vm = {1.0, 2.0};
    interpolation_from.va = {0.0, 0.2};
    interpolation_from.pg = {1.0};
    interpolation_from.qg = {2.0};
    interpolation_from.demand_factor = {0.8};
    interpolation_from.shunt_bs = {0.1};
    interpolation_from.shunt_steps = {{1, 0}};
    AcState interpolation_to = interpolation_from;
    interpolation_to.vm = {1.4, 2.4};
    interpolation_to.va = {0.4, 0.6};
    interpolation_to.pg = {3.0};
    interpolation_to.qg = {4.0};
    interpolation_to.demand_factor = {1.0};
    interpolation_to.shunt_bs = {0.3};
    interpolation_to.shunt_steps = {{0, 1}};
    const auto interpolation = interpolate_corrective_candidate(
        interpolation_from, interpolation_to, 0.25);
    if (std::abs(interpolation.vm[0] - 1.1) > 1e-12 ||
        std::abs(interpolation.va[1] - 0.3) > 1e-12 ||
        std::abs(interpolation.pg[0] - 1.5) > 1e-12 ||
        std::abs(interpolation.qg[0] - 2.5) > 1e-12 ||
        std::abs(interpolation.demand_factor[0] - 0.85) > 1e-12 ||
        interpolation.shunt_bs != interpolation_to.shunt_bs ||
        interpolation.shunt_steps != interpolation_to.shunt_steps) {
        throw std::runtime_error(
            "corrective candidate interpolation dropped a continuous or "
            "discrete control");
    }

    data.branches[0].rate_c = 3.0;
    data.branches[1].rate_c = 2.0;
    data.branches[2].rate_c = 1.0;
    const auto bounded_clusters =
        capacity_bounded_network_clusters(data, 2);
    if (bounded_clusters.size() != 2 ||
        bounded_clusters[0] != std::vector<int>({0, 1}) ||
        bounded_clusters[1] != std::vector<int>({2})) {
        throw std::runtime_error(
            "capacity-bounded economic clusters are not deterministic or "
            "connected");
    }
    const auto singleton_clusters =
        capacity_bounded_network_clusters(data, 1);
    if (singleton_clusters.size() != 3) {
        throw std::runtime_error(
            "capacity-bounded economic clusters violated their size cap");
    }

    CaseData economic_data;
    economic_data.buses.resize(1);
    economic_data.generators.resize(2);
    economic_data.generators[0].bus = 0;
    economic_data.generators[0].ncost = 2;
    economic_data.generators[0].cost = {0.0, 0.0, 2.0, 20.0};
    economic_data.generators[1].bus = 0;
    economic_data.generators[1].ncost = 2;
    economic_data.generators[1].cost = {0.0, 0.0, 2.0, 60.0};
    economic_data.loads.resize(1);
    economic_data.loads[0].bus = 0;
    economic_data.loads[0].pd_min = 0.0;
    economic_data.loads[0].pd_max = 3.0;
    economic_data.loads[0].ncost = 4;
    economic_data.loads[0].cost = {
        0.0, 0.0,
        1.0, 100.0,
        2.0, 150.0,
        3.0, 160.0,
    };
    const auto economic_target = build_component_economic_merit_target(
        economic_data, {{0}}, {true, true},
        {0.0, 0.0}, {2.0, 2.0},
        {0.0}, {3.0},
        {1.0, 0.0}, {1.0});
    if (!economic_target.feasible ||
        std::abs(economic_target.pg[0] - 2.0) > 1e-12 ||
        std::abs(economic_target.pg[1]) > 1e-12 ||
        std::abs(economic_target.load_power[0] - 2.0) > 1e-12 ||
        std::abs(economic_target.initial_market_surplus - 90.0) > 1e-12 ||
        std::abs(economic_target.target_market_surplus - 130.0) > 1e-12) {
        throw std::runtime_error(
            "component economic merit dispatch did not match the exact "
            "convex PWL solution");
    }

    economic_data.generators[0].ncost = 3;
    economic_data.generators[0].cost = {
        0.0, 0.0,
        1.0, 20.0,
        2.0, 30.0,
    };
    const auto nonconvex_target =
        build_component_economic_merit_target(
            economic_data, {{0}}, {true, true},
            {0.0, 0.0}, {2.0, 2.0},
            {0.0}, {3.0},
            {1.0, 0.0}, {1.0});
    if (nonconvex_target.feasible ||
        nonconvex_target.curve_order_valid) {
        throw std::runtime_error(
            "component economic merit dispatch accepted a nonconvex "
            "generator curve");
    }
}

struct FastContingencyPowerFlow::EconomicMeritTargetCache {
    bool initialized{};
    bool feasible{};
    std::vector<double> pg;
    std::vector<double> load_power;
    double initial_market_surplus{};
    double target_market_surplus{};
    double generation_movement{};
    double load_movement{};
    int component_count{};
};

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
    SparseMatrix jacobian;
    Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>> factorization;
    bool active_valid{};
    SparseMatrix active_jacobian;
    Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>>
        active_factorization;
    bool reactive_valid{};
    SparseMatrix reactive_jacobian;
    Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>>
        reactive_factorization;
    LowRankUpdate full_outage_update;
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
        jacobian.resize(dimension, dimension);
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
        full_outage_update = {};
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

        std::vector<Triplet> full_difference_entries;
        const std::array<int, 2> terminal_buses{from, to};
        const std::array<int, 4> flow_rows{
            angle_index[from],
            angle_count + voltage_index[from],
            angle_index[to],
            angle_count + voltage_index[to],
        };
        for (int flow = 0; flow < 4; ++flow) {
            const int row = flow_rows[flow];
            if (row < 0) {
                continue;
            }
            for (int variable_terminal = 0; variable_terminal < 2;
                 ++variable_terminal) {
                const int bus = terminal_buses[variable_terminal];
                const int voltage_column =
                    angle_count + voltage_index[bus];
                full_difference_entries.emplace_back(
                    row, voltage_column,
                    -derivative(flow, 2 * variable_terminal));
                const int angle_column = angle_index[bus];
                if (angle_column >= 0) {
                    full_difference_entries.emplace_back(
                        row, angle_column,
                        -derivative(flow, 1 + 2 * variable_terminal));
                }
            }
        }
        SparseMatrix full_difference(
            angle_count + voltage_count,
            angle_count + voltage_count);
        full_difference.setFromTriplets(
            full_difference_entries.begin(),
            full_difference_entries.end());

        std::vector<Triplet> active_difference_entries;
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

        const bool full_ready = prepare_low_rank_update(
            full_difference, factorization, full_outage_update);
        const bool active_ready = prepare_low_rank_update(
            active_difference, active_factorization, active_outage_update);
        const bool reactive_ready = prepare_low_rank_update(
            reactive_difference, reactive_factorization,
            reactive_outage_update);
        if (!full_ready || !active_ready || !reactive_ready) {
            full_outage_update = {};
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
        const Eigen::VectorXd step = solve_with_update(
            factorization, mismatch, full_outage_update);
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
        std::vector<double> generator_sensitivity(
            data.generators.size(), 0.0);
        for (int generator = 0;
             generator < static_cast<int>(data.generators.size());
             ++generator) {
            const int bus = data.generators[generator].bus;
            if (angle_index[bus] >= 0) {
                generator_sensitivity[generator] =
                    injection_sensitivity[angle_index[bus]];
            }
        }

        struct SensitivityExtreme {
            int generator{-1};
            double value{};
        };
        SensitivityExtreme minimum_down;
        SensitivityExtreme second_minimum_down;
        SensitivityExtreme maximum_down;
        SensitivityExtreme second_maximum_down;
        for (int down = 0;
             down < static_cast<int>(data.generators.size()); ++down) {
            if (pg[down] - p_lower[down] <= 1e-8) {
                continue;
            }
            const double sensitivity = generator_sensitivity[down];
            if (minimum_down.generator < 0 ||
                sensitivity < minimum_down.value) {
                second_minimum_down = minimum_down;
                minimum_down = {down, sensitivity};
            } else if (second_minimum_down.generator < 0 ||
                       sensitivity < second_minimum_down.value) {
                second_minimum_down = {down, sensitivity};
            }
            if (maximum_down.generator < 0 ||
                sensitivity > maximum_down.value) {
                second_maximum_down = maximum_down;
                maximum_down = {down, sensitivity};
            } else if (second_maximum_down.generator < 0 ||
                       sensitivity > second_maximum_down.value) {
                second_maximum_down = {down, sensitivity};
            }
        }
        int selected_up = -1;
        int selected_down = -1;
        double selected_coefficient = 0.0;
        for (int up = 0;
             up < static_cast<int>(data.generators.size()); ++up) {
            if (p_upper[up] - pg[up] <= 1e-8) {
                continue;
            }
            const auto& primary_down = desired_flow_change > 0.0
                ? minimum_down : maximum_down;
            const auto& alternate_down = desired_flow_change > 0.0
                ? second_minimum_down : second_maximum_down;
            const auto& chosen_down = primary_down.generator != up
                ? primary_down : alternate_down;
            if (chosen_down.generator < 0) {
                continue;
            }
            const double coefficient = generator_sensitivity[up] -
                chosen_down.value;
            if (std::abs(coefficient) <= 1e-10 ||
                desired_flow_change / coefficient <= 0.0) {
                continue;
            }
            if (selected_up < 0 ||
                std::abs(coefficient) >
                    std::abs(selected_coefficient)) {
                selected_up = up;
                selected_down = chosen_down.generator;
                selected_coefficient = coefficient;
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

    bool apply_distributed_active_flow_redispatch(
        const CaseData& data,
        int branch_index,
        bool from_side,
        double desired_flow_change,
        const std::vector<double>& p_lower,
        const std::vector<double>& p_upper,
        const std::vector<double>& load_power_lower,
        const std::vector<double>& load_power_upper,
        const AcState& reference_state,
        std::vector<double>& pg,
        std::vector<double>& load_power,
        std::vector<double>& p_spec) {
        if (!active_valid || branch_index < 0 ||
            branch_index >= static_cast<int>(data.branches.size()) ||
            std::abs(desired_flow_change) <= 1e-10 ||
            pg.size() != data.generators.size() ||
            p_lower.size() != data.generators.size() ||
            p_upper.size() != data.generators.size() ||
            load_power.size() != data.loads.size() ||
            load_power_lower.size() != data.loads.size() ||
            load_power_upper.size() != data.loads.size() ||
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
            derivative = -cosine_coefficient * std::sin(angle) +
                sine_coefficient * std::cos(angle);
        } else {
            const double cosine_coefficient =
                ((-g * tr - b * ti) / tm2) * voltage_product;
            const double sine_coefficient =
                ((-b * tr + g * ti) / tm2) * voltage_product;
            derivative = -cosine_coefficient * std::sin(angle) -
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

        const int generator_count =
            static_cast<int>(data.generators.size());
        const int load_count = static_cast<int>(data.loads.size());
        const int control_count = generator_count + load_count;
        std::vector<double> sensitivity(
            static_cast<std::size_t>(control_count), 0.0);
        std::vector<double> lower_change(
            static_cast<std::size_t>(control_count), 0.0);
        std::vector<double> upper_change(
            static_cast<std::size_t>(control_count), 0.0);
        std::vector<double> change(
            static_cast<std::size_t>(control_count), 0.0);
        std::vector<unsigned char> free(
            static_cast<std::size_t>(control_count), 0);
        for (int generator = 0; generator < generator_count; ++generator) {
            const int bus = data.generators[generator].bus;
            if (angle_index[bus] >= 0) {
                sensitivity[generator] =
                    injection_sensitivity[angle_index[bus]];
            }
            lower_change[generator] = std::max(
                p_lower[generator] - pg[generator], -1.0);
            upper_change[generator] = std::min(
                p_upper[generator] - pg[generator], 1.0);
            free[generator] =
                upper_change[generator] - lower_change[generator] > 1e-10;
        }
        // A positive load-control value below is an injection increase, so it
        // corresponds to reducing the load's consumed active power.  This
        // keeps the common balance equality sum(change) == 0 for generators
        // and loads while preserving each source t/ramp bound exactly.
        for (int load = 0; load < load_count; ++load) {
            const int control = generator_count + load;
            const int bus = data.loads[load].bus;
            if (angle_index[bus] >= 0) {
                sensitivity[control] =
                    injection_sensitivity[angle_index[bus]];
            }
            lower_change[control] = std::max(
                load_power[load] - load_power_upper[load], -1.0);
            upper_change[control] = std::min(
                load_power[load] - load_power_lower[load], 1.0);
            free[control] =
                upper_change[control] - lower_change[control] > 1e-10;
        }

        bool solved = false;
        for (int active_set_round = 0;
             active_set_round < 32; ++active_set_round) {
            double fixed_sum = 0.0;
            double fixed_flow = 0.0;
            double free_count = 0.0;
            double sensitivity_sum = 0.0;
            double sensitivity_square_sum = 0.0;
            for (int control = 0; control < control_count; ++control) {
                if (free[control]) {
                    free_count += 1.0;
                    sensitivity_sum += sensitivity[control];
                    sensitivity_square_sum +=
                        sensitivity[control] * sensitivity[control];
                } else {
                    fixed_sum += change[control];
                    fixed_flow +=
                        sensitivity[control] * change[control];
                }
            }
            if (free_count < 2.0) {
                break;
            }
            const double determinant =
                free_count * sensitivity_square_sum -
                sensitivity_sum * sensitivity_sum;
            if (determinant <= 1e-16) {
                break;
            }
            const double required_sum = -fixed_sum;
            const double required_flow =
                desired_flow_change - fixed_flow;
            const double constant_multiplier =
                (required_sum * sensitivity_square_sum -
                 required_flow * sensitivity_sum) / determinant;
            const double sensitivity_multiplier =
                (free_count * required_flow -
                 sensitivity_sum * required_sum) / determinant;
            bool clipped = false;
            for (int control = 0; control < control_count; ++control) {
                if (!free[control]) {
                    continue;
                }
                const double proposed = constant_multiplier +
                    sensitivity_multiplier * sensitivity[control];
                if (proposed < lower_change[control] - 1e-10) {
                    change[control] = lower_change[control];
                    free[control] = 0;
                    clipped = true;
                } else if (proposed >
                           upper_change[control] + 1e-10) {
                    change[control] = upper_change[control];
                    free[control] = 0;
                    clipped = true;
                } else {
                    change[control] = proposed;
                }
            }
            if (!clipped) {
                solved = true;
                break;
            }
        }
        if (!solved) {
            return false;
        }
        double balance_check = 0.0;
        double flow_check = 0.0;
        for (int control = 0; control < control_count; ++control) {
            balance_check += change[control];
            flow_check += sensitivity[control] * change[control];
        }
        if (std::abs(balance_check) > 1e-7 ||
            std::abs(flow_check - desired_flow_change) > 1e-6) {
            return false;
        }
        for (int generator = 0; generator < generator_count; ++generator) {
            if (std::abs(change[generator]) <= 1e-12) {
                continue;
            }
            pg[generator] += change[generator];
            p_spec[data.generators[generator].bus] += change[generator];
        }
        for (int load = 0; load < load_count; ++load) {
            const int control = generator_count + load;
            if (std::abs(change[control]) <= 1e-12) {
                continue;
            }
            load_power[load] -= change[control];
            p_spec[data.loads[load].bus] += change[control];
        }
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
        double damping = 1.0,
        bool enable_full_active_set_row_space = false,
        double target_band = 0.49) {
        if (!reactive_valid ||
            q_balance.size() != data.buses.size() ||
            !std::isfinite(target_band) || target_band < 0.0 ||
            target_band > 0.49) {
            return false;
        }
        std::vector<std::pair<double, int>> ranked_violations;
        for (int bus = 0; bus < static_cast<int>(data.buses.size()); ++bus) {
            const double excess =
                std::abs(q_balance[bus]) - target_band;
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
        constexpr std::size_t kLegacyMaximumViolationRows = 32;
        constexpr std::size_t kFullActiveSetThreshold =
            kLegacyMaximumViolationRows;
        constexpr std::size_t kFullMaximumViolationRows = 256;
        const bool use_full_active_set_row_space =
            enable_full_active_set_row_space &&
            ranked_violations.size() > kFullActiveSetThreshold;
        const std::size_t maximum_violation_rows =
            use_full_active_set_row_space ?
            kFullMaximumViolationRows :
            kLegacyMaximumViolationRows;
        if (ranked_violations.size() > maximum_violation_rows) {
            ranked_violations.resize(maximum_violation_rows);
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
        const std::size_t maximum_control_buses =
            use_full_active_set_row_space ? 1024 : 256;
        if (control_buses.size() > maximum_control_buses) {
            control_buses.resize(maximum_control_buses);
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
                std::clamp(
                    q_balance[bus], -target_band, target_band) -
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
        Eigen::VectorXd step;
        if (use_full_active_set_row_space) {
            // When more than 32 balance rows are simultaneously active, the
            // legacy truncation can chase a moving subset indefinitely.  The
            // full violated-row count is still much smaller than the number
            // of nearby voltage controls, so solve the regularized
            // minimum-norm system in row space and map it back to controls.
            Eigen::MatrixXd normal =
                sensitivity * sensitivity.transpose();
            const double maximum_diagonal =
                std::max(1e-12, normal.diagonal().maxCoeff());
            normal.diagonal().array() += 1e-6 * maximum_diagonal;
            const Eigen::VectorXd dual_step = normal.ldlt().solve(target);
            step = sensitivity.transpose() * dual_step;
        } else {
            // Preserve the proven legacy trajectory for the common small
            // active set; changing its regularization space can add predictor
            // iterations without improving feasibility.
            Eigen::MatrixXd normal =
                sensitivity.transpose() * sensitivity;
            const double maximum_diagonal =
                std::max(1e-12, normal.diagonal().maxCoeff());
            normal.diagonal().array() += 1e-6 * maximum_diagonal;
            const Eigen::VectorXd right_hand_side =
                sensitivity.transpose() * target;
            step = normal.ldlt().solve(right_hand_side);
        }
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
        {"local_balance_candidate_attempted",
         local_balance_candidate_attempted},
        {"local_balance_candidate_selected",
         local_balance_candidate_selected},
        {"local_balance_backtracking_attempts",
         local_balance_backtracking_attempts},
        {"local_balance_selected_step", local_balance_selected_step},
        {"local_balance_candidate_validation",
         local_balance_candidate_validation.to_json()},
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
        {"economic_balance_polish_attempted",
         economic_balance_polish_attempted},
        {"economic_merit_dispatch_attempted",
         economic_merit_dispatch_attempted},
        {"economic_merit_dispatch_cache_hit",
         economic_merit_dispatch_cache_hit},
        {"economic_merit_dispatch_feasible",
         economic_merit_dispatch_feasible},
        {"economic_merit_dispatch_applied",
         economic_merit_dispatch_applied},
        {"economic_merit_dispatch_surplus_gain",
         economic_merit_dispatch_surplus_gain},
        {"economic_merit_dispatch_generation_movement",
         economic_merit_dispatch_generation_movement},
        {"economic_merit_dispatch_load_movement",
         economic_merit_dispatch_load_movement},
        {"economic_merit_dispatch_component_count",
         economic_merit_dispatch_component_count},
        {"economic_merit_fallback_attempted",
         economic_merit_fallback_attempted},
        {"economic_merit_fallback_selected",
         economic_merit_fallback_selected},
        {"economic_merit_candidate_objective",
         economic_merit_candidate_objective},
        {"economic_merit_anchored_retry_attempted",
         economic_merit_anchored_retry_attempted},
        {"economic_merit_anchored_retry_selected",
         economic_merit_anchored_retry_selected},
        {"economic_merit_anchored_retry_objective",
         economic_merit_anchored_retry_objective},
        {"economic_balance_polish_threshold_passed",
         economic_balance_polish_threshold_passed},
        {"economic_balance_polish_objective_threshold",
         economic_balance_polish_objective_threshold},
        {"economic_linearized_objective_threshold",
         economic_linearized_objective_threshold},
        {"economic_balance_polish_time_limit_reached",
         economic_balance_polish_time_limit_reached},
        {"economic_balance_polish_wall_seconds",
         economic_balance_polish_wall_seconds},
        {"economic_exact_newton_attempted",
         economic_exact_newton_attempted},
        {"economic_exact_newton_converged",
         economic_exact_newton_converged},
        {"economic_exact_newton_selected",
         economic_exact_newton_selected},
        {"economic_exact_newton_iterations",
         economic_exact_newton_iterations},
        {"economic_exact_newton_q_limit_switches",
         economic_exact_newton_q_limit_switches},
        {"economic_exact_newton_wall_seconds",
         economic_exact_newton_wall_seconds},
        {"economic_exact_newton_objective",
         economic_exact_newton_objective},
        {"economic_exact_newton_failure_reason",
         economic_exact_newton_failure_reason},
        {"economic_exact_newton_validation",
         economic_exact_newton_validation.to_json()},
        {"economic_balance_polish_selected",
         economic_balance_polish_selected},
        {"economic_balance_polish_iterations",
         economic_balance_polish_iterations},
        {"economic_balance_polish_backtracking_attempts",
         economic_balance_polish_backtracking_attempts},
        {"economic_balance_polish_objective_before",
         economic_balance_polish_objective_before},
        {"economic_balance_polish_objective_after",
         economic_balance_polish_objective_after},
        {"economic_balance_polish_active_slack_before",
         economic_balance_polish_active_slack_before},
        {"economic_balance_polish_active_slack_after",
         economic_balance_polish_active_slack_after},
        {"economic_balance_polish_reactive_slack_before",
         economic_balance_polish_reactive_slack_before},
        {"economic_balance_polish_reactive_slack_after",
         economic_balance_polish_reactive_slack_after},
        {"economic_balance_polish_validation",
         economic_balance_polish_validation.to_json()},
        {"economic_balance_polish_trace",
         economic_balance_polish_trace},
        {"newton_candidate_selected", newton_candidate_selected},
        {"newton_candidate_validation", newton_candidate_validation.to_json()},
        {"active_only_newton_attempted", active_only_newton_attempted},
        {"active_only_newton_selected", active_only_newton_selected},
        {"active_only_newton_converged", active_only_newton_converged},
        {"active_only_newton_iterations", active_only_newton_iterations},
        {"active_only_newton_validation", active_only_newton_validation.to_json()},
        {"active_only_backtracking_attempts",
         active_only_backtracking_attempts},
        {"active_only_selected_step", active_only_selected_step},
        {"active_only_backtracking_validation",
         active_only_backtracking_validation.to_json()},
        {"reactive_only_newton_attempted", reactive_only_newton_attempted},
        {"reactive_only_newton_selected", reactive_only_newton_selected},
        {"reactive_only_newton_converged", reactive_only_newton_converged},
        {"reactive_only_newton_iterations", reactive_only_newton_iterations},
        {"reactive_only_newton_validation", reactive_only_newton_validation.to_json()},
        {"reactive_only_backtracking_attempts",
         reactive_only_backtracking_attempts},
        {"reactive_only_selected_step", reactive_only_selected_step},
        {"reactive_only_trace", reactive_only_trace},
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

static double rebuild_contingency_economic_fields(
    const CaseData& data,
    const AcState& base_state,
    const std::vector<int>& commitment,
    const Contingency& contingency,
    AcState& state) {
    if (commitment.size() != data.generators.size() ||
        base_state.pg.size() != data.generators.size() ||
        state.pg.size() != data.generators.size() ||
        state.demand_factor.size() != data.loads.size()) {
        throw std::runtime_error(
            "cannot rebuild dimensionally invalid contingency economics");
    }
    const int outaged_generator =
        contingency.type == ContingencyType::Generator
        ? contingency.component : -1;

    double objective = 0.0;
    state.gen_lambda.clear();
    for (int i = 0; i < static_cast<int>(data.generators.size()); ++i) {
        const auto& gen = data.generators[i];
        const bool active = commitment[i] == 1 && i != outaged_generator;
        if (!active) {
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

static double rebuild_contingency_state_fields(
    const CaseData& data,
    const AcState& base_state,
    const std::vector<int>& commitment,
    const Contingency& contingency,
    AcState& state,
    double balance_slack_upper,
    bool rebuild_economic_fields) {
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

    for (int i = 0; i < static_cast<int>(data.generators.size()); ++i) {
        const bool active = commitment[i] == 1 && i != outaged_generator;
        if (!active) {
            state.pg[i] = 0.0;
            state.qg[i] = 0.0;
        }
    }

    compute_branch_flows(data, outaged_branch, true, state);
    const auto balance = nodal_balance_slack_seed(
        data, state, balance_slack_upper, 1e-7);
    state.p_delta = balance.active;
    state.q_delta = balance.reactive;
    if (!rebuild_economic_fields) {
        return 0.0;
    }
    return rebuild_contingency_economic_fields(
        data, base_state, commitment, contingency, state);
}

double rebuild_contingency_state_derived_fields(
    const CaseData& data,
    const AcState& base_state,
    const std::vector<int>& commitment,
    const Contingency& contingency,
    AcState& state,
    double balance_slack_upper) {
    return rebuild_contingency_state_fields(
        data, base_state, commitment, contingency, state,
        balance_slack_upper, true);
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
    const auto source_components = connected_components(data, -1);
    normalize_source_reference_angles(data, source_components, state.va);

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
      options_(options),
      base_components_(connected_components(data, -1)),
      bridge_branch_(bridge_branches(data)) {
    if (commitment_.size() != data_.generators.size()) {
        throw std::runtime_error("fast power flow commitment has wrong length");
    }
    if (base_state_.vm.size() != data_.buses.size() ||
        base_state_.pg.size() != data_.generators.size()) {
        throw std::runtime_error("fast power flow base state has wrong dimensions");
    }
    if (!std::isfinite(options_.balance_cleanup_fraction) ||
        options_.balance_cleanup_fraction <= 0.0 ||
        options_.balance_cleanup_fraction > 1.0) {
        throw std::runtime_error(
            "fast power flow balance cleanup fraction must be in (0, 1]");
    }
    if (options_.economic_merit_cluster_max_buses <= 0 ||
        !std::isfinite(
            options_.economic_merit_fallback_objective_threshold)) {
        throw std::runtime_error(
            "fast power flow economic merit options are invalid");
    }
}

FastContingencyPowerFlow::~FastContingencyPowerFlow() = default;

FastPowerFlowResult FastContingencyPowerFlow::solve(
    const Contingency& contingency) const {
    return solve_with_merit_fallback(contingency);
}

FastPowerFlowResult FastContingencyPowerFlow::solve(
    const Contingency& contingency,
    const AcState& initial_state) const {
    return solve_with_merit_fallback(contingency, &initial_state);
}

FastPowerFlowResult FastContingencyPowerFlow::solve_with_merit_fallback(
    const Contingency& contingency,
    const AcState* initial_state) const {
    const auto wall_start = std::chrono::steady_clock::now();
    auto candidate = solve_impl(&contingency, initial_state);
    candidate.economic_merit_candidate_objective =
        candidate.solve.objective;
    const bool branch_merit_candidate =
        options_.economic_merit_dispatch &&
        contingency.type == ContingencyType::Branch;
    const bool fallback_required = branch_merit_candidate &&
        (!candidate.feasible || !std::isfinite(candidate.solve.objective) ||
         candidate.solve.objective <
             options_.economic_merit_fallback_objective_threshold);
    if (!fallback_required) {
        return candidate;
    }

    candidate.economic_merit_fallback_attempted = true;
    if (!economic_merit_fallback_) {
        auto fallback_options = options_;
        fallback_options.economic_merit_dispatch = false;
        economic_merit_fallback_ =
            std::make_unique<FastContingencyPowerFlow>(
                data_, base_state_, commitment_, fallback_options);
    }
    auto fallback = economic_merit_fallback_->solve_impl(
        &contingency, initial_state);
    const bool anchored_retry_attempted =
        options_.economic_merit_retry_from_fallback &&
        fallback.feasible && std::isfinite(fallback.solve.objective);
    FastPowerFlowResult anchored;
    if (anchored_retry_attempted) {
        // The first merit attempt can spend its entire continuation path
        // removing source-authorized overload slack. Retry from the verified
        // legacy state so every accepted continuation step must improve an
        // already secure, high-objective incumbent.
        anchored = solve_impl(&contingency, &fallback.solve.state);
        anchored.economic_merit_candidate_objective =
            anchored.solve.objective;
    }

    const FastPowerFlowResult* selected = &candidate;
    const auto prefer = [&](const FastPowerFlowResult& proposal) {
        return proposal.feasible &&
            std::isfinite(proposal.solve.objective) &&
            (!selected->feasible ||
             proposal.solve.objective > selected->solve.objective + 1e-9);
    };
    if (prefer(fallback)) {
        selected = &fallback;
    }
    if (anchored_retry_attempted && prefer(anchored)) {
        selected = &anchored;
    }

    FastPowerFlowResult output = *selected;
    const FastPowerFlowResult& merit_diagnostics =
        anchored_retry_attempted ? anchored : candidate;
    output.economic_merit_dispatch_attempted =
        merit_diagnostics.economic_merit_dispatch_attempted;
    output.economic_merit_dispatch_cache_hit =
        merit_diagnostics.economic_merit_dispatch_cache_hit;
    output.economic_merit_dispatch_feasible =
        merit_diagnostics.economic_merit_dispatch_feasible;
    output.economic_merit_dispatch_applied =
        merit_diagnostics.economic_merit_dispatch_applied;
    output.economic_merit_dispatch_surplus_gain =
        merit_diagnostics.economic_merit_dispatch_surplus_gain;
    output.economic_merit_dispatch_generation_movement =
        merit_diagnostics.economic_merit_dispatch_generation_movement;
    output.economic_merit_dispatch_load_movement =
        merit_diagnostics.economic_merit_dispatch_load_movement;
    output.economic_merit_dispatch_component_count =
        merit_diagnostics.economic_merit_dispatch_component_count;
    output.economic_merit_candidate_objective =
        candidate.solve.objective;
    output.economic_merit_fallback_attempted = true;
    output.economic_merit_fallback_selected = selected == &fallback;
    output.economic_merit_anchored_retry_attempted =
        anchored_retry_attempted;
    output.economic_merit_anchored_retry_selected =
        anchored_retry_attempted && selected == &anchored;
    output.economic_merit_anchored_retry_objective =
        anchored_retry_attempted ? anchored.solve.objective : 0.0;
    output.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_start).count();
    output.solve.wall_seconds = output.wall_seconds;
    return output;
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
    const auto maximum_thermal_ratio = [&](const AcState& state) {
        double maximum = 1.0;
        for (int i = 0; i < static_cast<int>(data_.branches.size()); ++i) {
            const double rating = branch_rating(i);
            if (i == outaged_branch || data_.branches[i].status == 0 ||
                rating <= 1e-12) {
                continue;
            }
            const auto& branch = data_.branches[i];
            const double slack = std::clamp(
                state.sm_slack[i], 0.0, data_.sm_vio_limit);
            const double from_scale = branch.transformer
                ? 1.0 + slack : state.vm[branch.from] + slack;
            const double to_scale = branch.transformer
                ? 1.0 + slack : state.vm[branch.to] + slack;
            const double from_limit = rating * std::max(1e-12, from_scale);
            const double to_limit = rating * std::max(1e-12, to_scale);
            maximum = std::max({
                maximum,
                std::hypot(state.pf[i], state.qf[i]) / from_limit,
                std::hypot(state.pt[i], state.qt[i]) / to_limit,
            });
        }
        return maximum;
    };
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
        direct_context->borrow_base_state(base_state_);
        direct_context->outaged_generator = outaged_generator;
        direct_context->outaged_branch = outaged_branch;
    }
    std::vector<std::vector<int>> bridge_outage_components;
    const std::vector<std::vector<int>>* components_pointer =
        &base_components_;
    if (outaged_branch >= 0 &&
        bridge_branch_[static_cast<std::size_t>(outaged_branch)] != 0) {
        bridge_outage_components =
            connected_components(data_, outaged_branch);
        components_pointer = &bridge_outage_components;
    }
    const auto& components = *components_pointer;
    AcState direct_state = initial_state;
    direct_state.pg = pg;
    direct_state.qg = qg;
    normalize_source_reference_angles(data_, components, direct_state.va);
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
    const bool direct_balance_cleanup_pending =
        base_mode &&
        ((options_.minimize_active_balance_slack &&
          std::accumulate(
              direct_state.p_delta.begin(),
              direct_state.p_delta.end(), 0.0) > 1e-9) ||
         (options_.minimize_reactive_balance_slack &&
          std::accumulate(
              direct_state.q_delta.begin(),
              direct_state.q_delta.end(), 0.0) > 1e-9));
    const bool contingency_economic_cleanup_pending =
        !base_mode && nb >= 16000 && options_.economic_balance_polish;
    if (output.direct_candidate_validation.max_residual <=
            options_.validation_tolerance &&
        !direct_balance_cleanup_pending &&
        !contingency_economic_cleanup_pending) {
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
    if (supplied_candidate_direct_only &&
        !contingency_economic_cleanup_pending) {
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

    std::vector<int> generator_outage_distance(nb, -1);
    if (outaged_generator >= 0) {
        const int outage_bus =
            data_.generators[outaged_generator].bus;
        std::queue<int> frontier;
        generator_outage_distance[outage_bus] = 0;
        frontier.push(outage_bus);
        while (!frontier.empty()) {
            const int bus = frontier.front();
            frontier.pop();
            const auto visit_branch = [&](int branch_index) {
                const auto& branch = data_.branches[branch_index];
                if (!branch.present || branch.status == 0) {
                    return;
                }
                const int neighbor = branch.from == bus
                    ? branch.to : branch.from;
                if (generator_outage_distance[neighbor] >= 0) {
                    return;
                }
                generator_outage_distance[neighbor] =
                    generator_outage_distance[bus] + 1;
                frontier.push(neighbor);
            };
            for (int branch : data_.buses[bus].branches_from) {
                visit_branch(branch);
            }
            for (int branch : data_.buses[bus].branches_to) {
                visit_branch(branch);
            }
        }
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
            std::vector<int> allocation_participants = participants;
            if (outaged_generator >= 0) {
                std::stable_sort(
                    allocation_participants.begin(),
                    allocation_participants.end(),
                    [&](int left, int right) {
                        const int left_distance =
                            generator_outage_distance[
                                data_.generators[left].bus];
                        const int right_distance =
                            generator_outage_distance[
                                data_.generators[right].bus];
                        const int normalized_left = left_distance >= 0
                            ? left_distance
                            : std::numeric_limits<int>::max();
                        const int normalized_right = right_distance >= 0
                            ? right_distance
                            : std::numeric_limits<int>::max();
                        if (normalized_left != normalized_right) {
                            return normalized_left < normalized_right;
                        }
                        return left < right;
                    });
            }
            const double prior_generation_total = active_target;
            active_target = std::clamp(
                active_target, active_lower, active_upper);
            if (outaged_generator >= 0) {
                allocate_total_in_order(
                    allocation_participants, p_lower, p_upper,
                    initial_state.pg, active_target, pg);
            } else {
                allocate_total(
                    participants, p_lower, p_upper,
                    initial_state.pg, active_target, pg);
            }
            reactive_target = std::clamp(
                reactive_target, reactive_lower, reactive_upper);
            if (outaged_generator >= 0) {
                allocate_total_in_order(
                    allocation_participants, q_lower, q_upper,
                    initial_state.qg, reactive_target, qg);
            } else {
                allocate_total(
                    participants, q_lower, q_upper,
                    initial_state.qg, reactive_target, qg);
            }
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

    std::vector<double> economic_merit_target_pg;
    std::vector<double> economic_merit_target_demand_factor;
    if (!base_mode && outaged_branch >= 0 &&
        options_.economic_merit_dispatch && nb >= 16000) {
        output.economic_merit_dispatch_attempted = true;
        const bool reusable_branch_target = outaged_branch >= 0 &&
            bridge_branch_[static_cast<std::size_t>(outaged_branch)] == 0;
        EconomicMeritTarget merit_target;
        if (reusable_branch_target && economic_merit_target_cache_ &&
            economic_merit_target_cache_->initialized) {
            output.economic_merit_dispatch_cache_hit = true;
            merit_target.feasible =
                economic_merit_target_cache_->feasible;
            merit_target.pg = economic_merit_target_cache_->pg;
            merit_target.load_power =
                economic_merit_target_cache_->load_power;
            merit_target.initial_market_surplus =
                economic_merit_target_cache_->initial_market_surplus;
            merit_target.target_market_surplus =
                economic_merit_target_cache_->target_market_surplus;
            merit_target.generation_movement =
                economic_merit_target_cache_->generation_movement;
            merit_target.load_movement =
                economic_merit_target_cache_->load_movement;
            output.economic_merit_dispatch_component_count =
                economic_merit_target_cache_->component_count;
        } else {
            const auto injection_components =
                capacity_bounded_network_clusters(
                    data_, options_.economic_merit_cluster_max_buses,
                    reusable_branch_target ? -1 : outaged_branch);
            output.economic_merit_dispatch_component_count =
                static_cast<int>(injection_components.size());
            // Keep responsive load fixed by default: active load movement also
            // changes reactive demand and can invalidate an otherwise useful
            // active-power merit target. It is exposed only as an explicit
            // diagnostic option. Cluster size therefore controls generator
            // sharing without silently changing the source load profile.
            const std::vector<double> fixed_load_power =
                predictor_load_power;
            const auto& merit_load_lower =
                options_.economic_merit_allow_load_movement
                ? predictor_load_power_lower : fixed_load_power;
            const auto& merit_load_upper =
                options_.economic_merit_allow_load_movement
                ? predictor_load_power_upper : fixed_load_power;
            merit_target = build_component_economic_merit_target(
                data_, injection_components, active, p_lower, p_upper,
                merit_load_lower, merit_load_upper,
                pg, predictor_load_power);
            if (reusable_branch_target) {
                if (!economic_merit_target_cache_) {
                    economic_merit_target_cache_ =
                        std::make_unique<EconomicMeritTargetCache>();
                }
                economic_merit_target_cache_->initialized = true;
                economic_merit_target_cache_->feasible =
                    merit_target.feasible;
                economic_merit_target_cache_->pg = merit_target.pg;
                economic_merit_target_cache_->load_power =
                    merit_target.load_power;
                economic_merit_target_cache_->initial_market_surplus =
                    merit_target.initial_market_surplus;
                economic_merit_target_cache_->target_market_surplus =
                    merit_target.target_market_surplus;
                economic_merit_target_cache_->generation_movement =
                    merit_target.generation_movement;
                economic_merit_target_cache_->load_movement =
                    merit_target.load_movement;
                economic_merit_target_cache_->component_count =
                    output.economic_merit_dispatch_component_count;
            }
        }
        output.economic_merit_dispatch_feasible = merit_target.feasible;
        output.economic_merit_dispatch_surplus_gain = data_.delta_ctg *
            (merit_target.target_market_surplus -
             merit_target.initial_market_surplus);
        output.economic_merit_dispatch_generation_movement =
            merit_target.generation_movement;
        output.economic_merit_dispatch_load_movement =
            merit_target.load_movement;
        if (merit_target.feasible &&
            output.economic_merit_dispatch_surplus_gain > 1e-9 &&
            merit_target.pg.size() == data_.generators.size() &&
            merit_target.load_power.size() == data_.loads.size()) {
            economic_merit_target_pg = std::move(merit_target.pg);
            economic_merit_target_demand_factor =
                predictor_demand_factor;
            for (int load_index = 0;
                 load_index < static_cast<int>(data_.loads.size());
                 ++load_index) {
                if (std::abs(data_.loads[load_index].pd_nominal) > 1e-12) {
                    economic_merit_target_demand_factor[load_index] =
                        merit_target.load_power[load_index] /
                            data_.loads[load_index].pd_nominal;
                }
            }
            output.economic_merit_dispatch_applied = true;
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
            // Preserve the economically preferred contingency controls before
            // the feasibility predictor starts making local bus-by-bus
            // corrections.  Those local corrections are excellent at finding
            // a secure point, but they can move generation and load far away
            // from the verified base dispatch.  The refreshed-Jacobian Newton
            // rescue below should solve the AC equations for this preferred
            // control target, using the secure predictor only as its voltage
            // initialization.
            AcState economic_target_state = predictor_state;
            if (output.economic_merit_dispatch_applied) {
                economic_target_state.pg = economic_merit_target_pg;
                economic_target_state.demand_factor =
                    economic_merit_target_demand_factor;
            }
            ValidationReport predictor_validation =
                output.direct_candidate_validation;
            std::unique_ptr<FixedJacobianPredictorCache>
                contingency_predictor_cache;
            int active_feasibility_repair_attempts = 0;
            int ac_feasibility_repair_attempts = 0;
            int proactive_local_search_attempts = 0;
            int consecutive_full_local_reactive_active_angle = 0;
            std::string prior_selected_correction_mode;
            double prior_selected_damping = 0.0;
            std::vector<std::pair<int, double>> coordinate_history;
            double initial_predictor_validation_residual =
                std::numeric_limits<double>::infinity();
            // These work vectors are reused across predictor iterations and
            // candidate line-search trials. On large cases, allocating and
            // zero-initializing fresh case-sized vectors for every candidate
            // creates avoidable allocator contention across resident workers.
            std::vector<double> p_network(
                static_cast<std::size_t>(nb), 0.0);
            std::vector<double> q_network(
                static_cast<std::size_t>(nb), 0.0);
            std::vector<double> predictor_load_reactive(
                data_.loads.size(), 0.0);
            std::vector<double> predictor_load_active(
                data_.loads.size(), 0.0);
            std::vector<double> trial_p_network(
                static_cast<std::size_t>(nb), 0.0);
            std::vector<double> trial_q_network(
                static_cast<std::size_t>(nb), 0.0);
            std::vector<double> trial_load_active(
                data_.loads.size(), 0.0);
            // Difficult generator outages can enter a late, monotonically
            // shrinking cycle between active- and reactive-flow limits after
            // the global feasibility repair.  Keep the inexpensive local
            // corrections alive long enough to cross the exact validator's
            // tolerance instead of stopping at an arbitrary iteration edge.
            // Productive paths still return immediately and true stalls still
            // break when no improving candidate exists, so this ceiling only
            // extends the late monotone active-flow alternation.
            constexpr int kMaximumFixedJacobianIterations = 224;
            for (int predictor_iteration = 0;
                 predictor_iteration <= kMaximumFixedJacobianIterations;
                 ++predictor_iteration) {
                const auto projection_validation_start =
                    std::chrono::steady_clock::now();
                normalize_source_reference_angles(
                    data_, components, predictor_state.va);
                // The predictor changes only the physical controls below.
                // Recompute branch flows once here, then seed nodal slacks
                // after the local P/Q projections have finished.  Calling
                // rebuild_contingency_state_fields at both points performed
                // an otherwise unused full nodal-balance pass before those
                // projections on every predictor iteration.
                compute_branch_flows(
                    data_, outaged_branch, true, predictor_state);

                std::fill(p_network.begin(), p_network.end(), 0.0);
                std::fill(q_network.begin(), q_network.end(), 0.0);
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
                    static_cast<void>(allocate_total(
                        active_at_bus[bus], q_lower, q_upper,
                        initial_state.qg,
                        generated_q + q_balance -
                            std::clamp(q_balance, -0.49, 0.49),
                        predictor_state.qg));
                }

                // A PQ load is also a source-authorized corrective control.
                // If local Q generation cannot put a bus inside the soft
                // balance band, move only the loads at that bus and respect
                // the exact contingency ramp and tmin/tmax bounds computed
                // above. This is especially important for voltage-limited
                // load buses where a voltage correction has no remaining
                // feasible direction.
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
                        if (allocate_total(
                                active_at_bus[bus], p_lower, p_upper,
                                predictor_state.pg, target_generation,
                                predictor_state.pg)) {
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
                    if (!allocate_total(
                            data_.buses[bus].loads,
                            predictor_load_power_lower,
                            predictor_load_power_upper,
                            predictor_load_active, target_load,
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

                // Rank nonfinal predictor states using the same physical
                // source constraints as the full validator.  Economic PWL
                // fields and branch Ohm residuals are deterministic rebuilds
                // of Pg/load and the branch flows already computed above, so
                // recreating and rescanning them for every infeasible trial
                // adds no acceptance protection.  A candidate that crosses
                // the physical tolerance is rebuilt economically and gets
                // exactly the omitted PWL and Ohm-law checks before it can be
                // accepted.  The already-certified physical state is not
                // recomputed or rescanned.
                const auto predictor_balance =
                    nodal_balance_slack_seed_from_network(
                        data_, predictor_state, p_network, q_network,
                        0.5, 1e-7);
                predictor_state.p_delta = predictor_balance.active;
                predictor_state.q_delta = predictor_balance.reactive;
                predictor_validation =
                    validate_rebuilt_contingency_predictor(
                        data_, predictor_state, commitment_,
                        *direct_context);
                double predictor_objective = 0.0;
                if (predictor_validation.max_residual <=
                    options_.validation_tolerance) {
                    predictor_objective =
                        rebuild_contingency_economic_fields(
                            data_, base_state_, commitment_, *contingency,
                            predictor_state);
                    predictor_validation =
                        validate_rebuilt_contingency_economic_and_ohms(
                            data_, predictor_state, commitment_,
                            *direct_context, predictor_validation);
                }
                output.fixed_jacobian_predictor_iterations =
                    predictor_iteration;
                output.fixed_jacobian_predictor_validation =
                    predictor_validation;
                if (predictor_iteration == 0) {
                    initial_predictor_validation_residual =
                        predictor_validation.max_residual;
                }
                const double projection_validation_seconds =
                    std::chrono::duration<double>(
                        std::chrono::steady_clock::now() -
                        projection_validation_start).count();
                if (options_.capture_diagnostics) {
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
                        {"projection_validation_seconds",
                         projection_validation_seconds},
                    });
                }
                retain_best_intermediate(
                    predictor_state, predictor_validation,
                    "fixed_base_jacobian_predictor");
                if (predictor_validation.max_residual <=
                    options_.validation_tolerance) {
                    output.economic_balance_polish_objective_threshold =
                        options_
                            .economic_balance_polish_objective_threshold;
                    output.economic_linearized_objective_threshold =
                        options_.economic_linearized_objective_threshold;
                    output.economic_balance_polish_threshold_passed =
                        predictor_objective <=
                        options_
                            .economic_balance_polish_objective_threshold;
                    if (options_.economic_balance_polish &&
                        output.economic_balance_polish_threshold_passed &&
                        options_.max_economic_balance_polish_iterations > 0) {
                        output.economic_balance_polish_attempted = true;
                        const auto economic_polish_start =
                            std::chrono::steady_clock::now();
                        const auto economic_polish_elapsed = [&]() {
                            return std::chrono::duration<double>(
                                std::chrono::steady_clock::now() -
                                economic_polish_start).count();
                        };
                        const auto economic_polish_budget_exhausted = [&]() {
                            if (!std::isfinite(
                                    options_
                                        .economic_balance_polish_wall_seconds)) {
                                return false;
                            }
                            const bool exhausted =
                                economic_polish_elapsed() >= std::max(
                                    0.0,
                                    options_
                                        .economic_balance_polish_wall_seconds);
                            output.economic_balance_polish_time_limit_reached =
                                output.economic_balance_polish_time_limit_reached ||
                                exhausted;
                            return exhausted;
                        };
                        const auto bounded_economic_solver_seconds =
                            [&](double requested_seconds) {
                                if (!std::isfinite(
                                        options_
                                            .economic_balance_polish_wall_seconds)) {
                                    return requested_seconds;
                                }
                                return std::max(
                                    0.0,
                                    std::min(
                                        requested_seconds,
                                        options_
                                                .economic_balance_polish_wall_seconds -
                                            economic_polish_elapsed()));
                            };
                        const auto slack_sum = [](const std::vector<double>& values) {
                            return std::accumulate(
                                values.begin(), values.end(), 0.0);
                        };
                        output.economic_balance_polish_objective_before =
                            predictor_objective;
                        output.economic_balance_polish_active_slack_before =
                            slack_sum(predictor_state.p_delta);
                        output.economic_balance_polish_reactive_slack_before =
                            slack_sum(predictor_state.q_delta);

                        AcState polished_state = predictor_state;
                        double polished_objective = predictor_objective;
                        ValidationReport polished_validation =
                            predictor_validation;
                        bool exact_newton_reactive_trigger = false;
                        bool outage_update_ready = true;
                        if (outaged_branch >= 0) {
                            outage_update_ready =
                                predictor_cache_->configure_branch_outage_update(
                                    data_, polished_state, outaged_branch);
                        }
                        if (!outage_update_ready) {
                            output.economic_balance_polish_trace.push_back({
                                {"selected", false},
                                {"reason", "branch outage low-rank update unavailable"},
                            });
                        } else {
                            const YRows polish_ybus = build_ybus(
                                data_, outaged_branch, &polished_state);
                            // After the cached network correction, project
                            // controls onto the resulting injections.  This
                            // leaves every
                            // network flow unchanged, so it cannot create a new
                            // thermal or angle violation.  At generator buses,
                            // bounded Pg/Qg absorb the local mismatch.  At
                            // load-only buses, source-authorized corrective
                            // demand movement absorbs as much active mismatch
                            // as its ramp and quantity bounds allow.  Any
                            // remainder stays in the explicit paid balance
                            // variables and the full validator remains the
                            // acceptance gate.
                            const auto apply_local_injection_projection = [&]() {
                            AcState local_projection = polished_state;
                            std::vector<double> local_p_network;
                            std::vector<double> local_q_network;
                            network_injections(
                                polish_ybus, local_projection.vm,
                                local_projection.va, local_p_network,
                                local_q_network);
                            std::vector<double> local_load_power(
                                data_.loads.size(), 0.0);
                            for (int load = 0;
                                 load < static_cast<int>(data_.loads.size());
                                 ++load) {
                                local_load_power[load] =
                                    data_.loads[load].pd_nominal *
                                    local_projection.demand_factor[load];
                            }
                            for (int bus = 0; bus < nb; ++bus) {
                                if (active_at_bus[bus].empty()) {
                                    std::vector<int> adjustable_loads;
                                    double total_lower = 0.0;
                                    double total_upper = 0.0;
                                    for (int load : data_.buses[bus].loads) {
                                        if (std::abs(
                                                data_.loads[load].pd_nominal) <=
                                            1e-12) {
                                            continue;
                                        }
                                        adjustable_loads.push_back(load);
                                        total_lower +=
                                            predictor_load_power_lower[load];
                                        total_upper +=
                                            predictor_load_power_upper[load];
                                    }
                                    if (!adjustable_loads.empty()) {
                                        if (adjustable_loads.size() == 1) {
                                            const int load =
                                                adjustable_loads.front();
                                            const auto& source_load =
                                                data_.loads[load];
                                            const double factor_first =
                                                predictor_load_power_lower[load] /
                                                source_load.pd_nominal;
                                            const double factor_second =
                                                predictor_load_power_upper[load] /
                                                source_load.pd_nominal;
                                            double feasible_factor_lower =
                                                std::min(
                                                    factor_first,
                                                    factor_second);
                                            double feasible_factor_upper =
                                                std::max(
                                                    factor_first,
                                                    factor_second);
                                            const auto restrict_factor =
                                                [&](double offset,
                                                    double slope) {
                                                    constexpr double
                                                        kSlackLimit = 0.5;
                                                    if (std::abs(slope) <=
                                                        1e-14) {
                                                        if (std::abs(offset) >
                                                            kSlackLimit +
                                                                1e-9) {
                                                            feasible_factor_upper =
                                                                feasible_factor_lower -
                                                                1.0;
                                                        }
                                                        return;
                                                    }
                                                    const double first =
                                                        (-kSlackLimit - offset) /
                                                        slope;
                                                    const double second =
                                                        (kSlackLimit - offset) /
                                                        slope;
                                                    feasible_factor_lower =
                                                        std::max(
                                                            feasible_factor_lower,
                                                            std::min(
                                                                first,
                                                                second));
                                                    feasible_factor_upper =
                                                        std::min(
                                                            feasible_factor_upper,
                                                            std::max(
                                                                first,
                                                                second));
                                                };
                                            restrict_factor(
                                                local_p_network[bus],
                                                source_load.pd_nominal);
                                            restrict_factor(
                                                local_q_network[bus],
                                                source_load.qd_nominal);
                                            if (feasible_factor_lower <=
                                                feasible_factor_upper + 1e-12) {
                                                std::vector<double>
                                                    factor_candidates{
                                                        std::clamp(
                                                            local_projection
                                                                .demand_factor[load],
                                                            feasible_factor_lower,
                                                            feasible_factor_upper),
                                                        feasible_factor_lower,
                                                        feasible_factor_upper};
                                                factor_candidates.push_back(
                                                    std::clamp(
                                                        -local_p_network[bus] /
                                                            source_load.pd_nominal,
                                                        feasible_factor_lower,
                                                        feasible_factor_upper));
                                                if (std::abs(
                                                        source_load.qd_nominal) >
                                                    1e-14) {
                                                    factor_candidates.push_back(
                                                        std::clamp(
                                                            -local_q_network[bus] /
                                                                source_load
                                                                    .qd_nominal,
                                                            feasible_factor_lower,
                                                            feasible_factor_upper));
                                                }
                                                const auto load_points =
                                                    active_pwl_points(
                                                        source_load.cost,
                                                        source_load.ncost,
                                                        predictor_load_power_lower[
                                                            load],
                                                        predictor_load_power_upper[
                                                            load]);
                                                for (const auto& point :
                                                     load_points) {
                                                    factor_candidates.push_back(
                                                        std::clamp(
                                                            point.mw /
                                                                source_load
                                                                    .pd_nominal,
                                                            feasible_factor_lower,
                                                            feasible_factor_upper));
                                                }
                                                double best_factor =
                                                    factor_candidates.front();
                                                double best_local_objective =
                                                    -std::numeric_limits<double>::
                                                        infinity();
                                                for (double factor :
                                                     factor_candidates) {
                                                    const double load_power =
                                                        source_load.pd_nominal *
                                                        factor;
                                                    const auto weights =
                                                        interpolation_weights(
                                                            load_points,
                                                            load_power);
                                                    double load_benefit = 0.0;
                                                    for (std::size_t point = 0;
                                                         point <
                                                             load_points.size();
                                                         ++point) {
                                                        load_benefit +=
                                                            load_points[point].cost *
                                                            weights[point];
                                                    }
                                                    const double active_mismatch =
                                                        local_p_network[bus] +
                                                        load_power;
                                                    const double reactive_mismatch =
                                                        local_q_network[bus] +
                                                        source_load.qd_nominal *
                                                            factor;
                                                    const double local_objective =
                                                        load_benefit -
                                                        data_
                                                                .p_delta_cost_approx *
                                                            std::abs(
                                                                active_mismatch) -
                                                        data_
                                                                .q_delta_cost_approx *
                                                            std::abs(
                                                                reactive_mismatch);
                                                    if (local_objective >
                                                        best_local_objective +
                                                            1e-9) {
                                                        best_local_objective =
                                                            local_objective;
                                                        best_factor = factor;
                                                    }
                                                }
                                                local_projection
                                                    .demand_factor[load] =
                                                    best_factor;
                                                local_load_power[load] =
                                                    source_load.pd_nominal *
                                                    best_factor;
                                            }
                                        } else {
                                        std::vector<double> prior_load_power;
                                        prior_load_power.reserve(
                                            adjustable_loads.size());
                                        double prior_load_p = 0.0;
                                        double prior_load_q = 0.0;
                                        for (int load : adjustable_loads) {
                                            prior_load_power.push_back(
                                                local_load_power[load]);
                                            prior_load_p +=
                                                local_load_power[load];
                                            prior_load_q +=
                                                data_.loads[load].qd_nominal *
                                                local_projection
                                                    .demand_factor[load];
                                        }
                                        const double target_load = std::clamp(
                                            -local_p_network[bus],
                                            total_lower, total_upper);
                                        static_cast<void>(allocate_total(
                                            adjustable_loads,
                                            predictor_load_power_lower,
                                            predictor_load_power_upper,
                                            local_load_power, target_load,
                                            local_load_power));
                                        double proposed_load_p = 0.0;
                                        double proposed_load_q = 0.0;
                                        for (int load : adjustable_loads) {
                                            proposed_load_p +=
                                                local_load_power[load];
                                            proposed_load_q +=
                                                data_.loads[load].qd_nominal *
                                                local_load_power[load] /
                                                data_.loads[load].pd_nominal;
                                        }
                                        // A single global backtracking step is
                                        // needlessly limited by the worst
                                        // load-only bus.  Bound each bus's
                                        // interpolation independently so both
                                        // of its official 0.5-p.u. P/Q balance
                                        // variables remain sufficient.
                                        double feasible_step_lower = 0.0;
                                        double feasible_step_upper = 1.0;
                                        const auto restrict_balance_step =
                                            [&](double prior_mismatch,
                                                double proposed_mismatch) {
                                                constexpr double kSlackLimit =
                                                    0.5;
                                                const double direction =
                                                    proposed_mismatch -
                                                    prior_mismatch;
                                                if (std::abs(direction) <=
                                                    1e-14) {
                                                    if (std::abs(
                                                            prior_mismatch) >
                                                        kSlackLimit + 1e-9) {
                                                        feasible_step_upper =
                                                            -1.0;
                                                    }
                                                    return;
                                                }
                                                const double first =
                                                    (-kSlackLimit -
                                                     prior_mismatch) /
                                                    direction;
                                                const double second =
                                                    (kSlackLimit -
                                                     prior_mismatch) /
                                                    direction;
                                                feasible_step_lower = std::max(
                                                    feasible_step_lower,
                                                    std::min(first, second));
                                                feasible_step_upper = std::min(
                                                    feasible_step_upper,
                                                    std::max(first, second));
                                            };
                                        restrict_balance_step(
                                            local_p_network[bus] +
                                                prior_load_p,
                                            local_p_network[bus] +
                                                proposed_load_p);
                                        restrict_balance_step(
                                            local_q_network[bus] +
                                                prior_load_q,
                                            local_q_network[bus] +
                                                proposed_load_q);
                                        const double bus_step =
                                            feasible_step_upper + 1e-12 >=
                                                    feasible_step_lower
                                            ? std::clamp(
                                                  feasible_step_upper,
                                                  0.0, 1.0)
                                            : 0.0;
                                        for (std::size_t position = 0;
                                             position <
                                                 adjustable_loads.size();
                                             ++position) {
                                            const int load =
                                                adjustable_loads[position];
                                            local_load_power[load] =
                                                prior_load_power[position] +
                                                bus_step *
                                                    (local_load_power[load] -
                                                     prior_load_power[position]);
                                            local_projection
                                                .demand_factor[load] =
                                                local_load_power[load] /
                                                data_.loads[load].pd_nominal;
                                        }
                                        }
                                    }
                                }

                                if (active_at_bus[bus].empty()) {
                                    continue;
                                }
                                double load_p = 0.0;
                                double load_q = 0.0;
                                for (int load : data_.buses[bus].loads) {
                                    load_p += data_.loads[load].pd_nominal *
                                        local_projection.demand_factor[load];
                                    load_q += data_.loads[load].qd_nominal *
                                        local_projection.demand_factor[load];
                                }
                                double p_bus_lower = 0.0;
                                double p_bus_upper = 0.0;
                                double q_bus_lower = 0.0;
                                double q_bus_upper = 0.0;
                                for (int generator : active_at_bus[bus]) {
                                    p_bus_lower += p_lower[generator];
                                    p_bus_upper += p_upper[generator];
                                    q_bus_lower += q_lower[generator];
                                    q_bus_upper += q_upper[generator];
                                }
                                const double target_pg = std::clamp(
                                    local_p_network[bus] + load_p,
                                    p_bus_lower, p_bus_upper);
                                const double target_qg = std::clamp(
                                    local_q_network[bus] + load_q,
                                    q_bus_lower, q_bus_upper);
                                static_cast<void>(allocate_total(
                                    active_at_bus[bus], p_lower, p_upper,
                                    polished_state.pg, target_pg,
                                    local_projection.pg));
                                static_cast<void>(allocate_total(
                                    active_at_bus[bus], q_lower, q_upper,
                                    polished_state.qg, target_qg,
                                    local_projection.qg));
                            }
                            ensure_shunt_control_state(
                                data_, local_projection);
                            for (int bus = 0; bus < nb; ++bus) {
                                if (data_.buses[bus].shunts.empty()) {
                                    continue;
                                }
                                double load_q = 0.0;
                                double generation_q = 0.0;
                                for (int load : data_.buses[bus].loads) {
                                    load_q += data_.loads[load].qd_nominal *
                                        local_projection.demand_factor[load];
                                }
                                for (int generator : active_at_bus[bus]) {
                                    generation_q +=
                                        local_projection.qg[generator];
                                }
                                double q_balance =
                                    local_q_network[bus] - generation_q +
                                    load_q;
                                for (int coordinate_pass = 0;
                                     coordinate_pass < 512;
                                     ++coordinate_pass) {
                                    double best_absolute =
                                        std::abs(q_balance);
                                    double best_balance = q_balance;
                                    int best_shunt = -1;
                                    int best_block = -1;
                                    int best_step_change = 0;
                                    for (int shunt_index :
                                         data_.buses[bus].shunts) {
                                        const auto& shunt =
                                            data_.shunts[shunt_index];
                                        if (!shunt.dispatchable ||
                                            shunt_index >= static_cast<int>(
                                                local_projection
                                                    .shunt_steps.size())) {
                                            continue;
                                        }
                                        for (int block = 0;
                                             block < static_cast<int>(
                                                 shunt
                                                     .block_maximum_steps.size());
                                             ++block) {
                                            if (block >= static_cast<int>(
                                                    local_projection
                                                        .shunt_steps[shunt_index]
                                                        .size())) {
                                                continue;
                                            }
                                            const int current_step =
                                                local_projection
                                                    .shunt_steps[shunt_index]
                                                                [block];
                                            for (int step_change : {-1, 1}) {
                                                const int proposed_step =
                                                    current_step + step_change;
                                                if (proposed_step < 0 ||
                                                    proposed_step >
                                                        shunt
                                                            .block_maximum_steps
                                                                [block]) {
                                                    continue;
                                                }
                                                const double delta_bs =
                                                    static_cast<double>(
                                                        step_change) *
                                                    shunt
                                                        .block_susceptance[block];
                                                const double proposed_balance =
                                                    q_balance - delta_bs *
                                                        local_projection.vm[bus] *
                                                        local_projection.vm[bus];
                                                if (std::abs(proposed_balance) +
                                                        1e-12 <
                                                    best_absolute) {
                                                    best_absolute =
                                                        std::abs(
                                                            proposed_balance);
                                                    best_balance =
                                                        proposed_balance;
                                                    best_shunt = shunt_index;
                                                    best_block = block;
                                                    best_step_change =
                                                        step_change;
                                                }
                                            }
                                        }
                                    }
                                    if (best_shunt < 0) {
                                        break;
                                    }
                                    local_projection
                                        .shunt_steps[best_shunt][best_block] +=
                                        best_step_change;
                                    local_projection.shunt_bs[best_shunt] +=
                                        static_cast<double>(best_step_change) *
                                        data_.shunts[best_shunt]
                                            .block_susceptance[best_block];
                                    q_balance = best_balance;
                                }
                                if (!active_at_bus[bus].empty()) {
                                    double current_qg = 0.0;
                                    double total_lower = 0.0;
                                    double total_upper = 0.0;
                                    for (int generator :
                                         active_at_bus[bus]) {
                                        current_qg +=
                                            local_projection.qg[generator];
                                        total_lower += q_lower[generator];
                                        total_upper += q_upper[generator];
                                    }
                                    const double target_qg = std::clamp(
                                        current_qg + q_balance,
                                        total_lower, total_upper);
                                    static_cast<void>(allocate_total(
                                        active_at_bus[bus], q_lower, q_upper,
                                        local_projection.qg, target_qg,
                                        local_projection.qg));
                                }
                            }
                            const std::array<double, 11>
                                local_projection_steps{
                                    1.0, 0.75, 0.5, 0.25, 0.125,
                                    0.0625, 0.03125, 0.015625,
                                    0.0078125, 0.00390625, 0.001953125};
                            bool local_projection_selected = false;
                            double local_projection_selected_step = 0.0;
                            double local_projection_objective =
                                polished_objective;
                            AcState local_projection_selected_state =
                                polished_state;
                            ValidationReport local_projection_validation =
                                polished_validation;
                            ValidationReport local_projection_full_validation;
                            double local_projection_full_objective =
                                polished_objective;
                            for (double step : local_projection_steps) {
                                AcState candidate = polished_state;
                                if (step == 1.0) {
                                    candidate.shunt_steps =
                                        local_projection.shunt_steps;
                                    candidate.shunt_bs =
                                        local_projection.shunt_bs;
                                }
                                for (int generator = 0; generator < ng;
                                     ++generator) {
                                    candidate.pg[generator] += step *
                                        (local_projection.pg[generator] -
                                         polished_state.pg[generator]);
                                    candidate.qg[generator] += step *
                                        (local_projection.qg[generator] -
                                         polished_state.qg[generator]);
                                }
                                for (int load = 0;
                                     load <
                                         static_cast<int>(data_.loads.size());
                                     ++load) {
                                    candidate.demand_factor[load] += step *
                                        (local_projection.demand_factor[load] -
                                         polished_state.demand_factor[load]);
                                }
                                const double candidate_objective =
                                    rebuild_contingency_state_derived_fields(
                                        data_, base_state_, commitment_,
                                        *contingency, candidate);
                                const auto candidate_validation =
                                    validate_state(
                                        data_, ModelMode::ContingencySoft,
                                        candidate, commitment_,
                                        direct_context);
                                if (step == 1.0) {
                                    local_projection_full_objective =
                                        candidate_objective;
                                    local_projection_full_validation =
                                        candidate_validation;
                                }
                                if (candidate_validation.max_residual <=
                                        options_.validation_tolerance &&
                                    candidate_objective >
                                        polished_objective + 1e-9) {
                                    local_projection_selected = true;
                                    local_projection_selected_step = step;
                                    local_projection_objective =
                                        candidate_objective;
                                    local_projection_selected_state =
                                        std::move(candidate);
                                    local_projection_validation =
                                        candidate_validation;
                                    break;
                                }
                            }
                            output.economic_balance_polish_trace.push_back({
                                {"phase", "local_injection_projection"},
                                {"selected", local_projection_selected},
                                {"selected_step",
                                 local_projection_selected_step},
                                {"objective_before", polished_objective},
                                {"objective_after",
                                 local_projection_objective},
                                {"full_step_objective",
                                 local_projection_full_objective},
                                {"active_slack_after",
                                 slack_sum(
                                     local_projection_selected_state.p_delta)},
                                {"reactive_slack_after",
                                 slack_sum(
                                     local_projection_selected_state.q_delta)},
                                {"validation",
                                 local_projection_validation.to_json()},
                                {"full_step_validation",
                                 local_projection_full_validation.to_json()},
                            });
                            if (local_projection_selected) {
                                polished_state =
                                    std::move(local_projection_selected_state);
                                polished_objective =
                                    local_projection_objective;
                                polished_validation =
                                    local_projection_validation;
                            }
                            };
                            std::vector<std::vector<int>>
                                component_generators(components.size());
                            for (int component = 0;
                                 component < static_cast<int>(components.size());
                                 ++component) {
                                for (int bus : components[component]) {
                                    component_generators[component].insert(
                                        component_generators[component].end(),
                                        active_at_bus[bus].begin(),
                                        active_at_bus[bus].end());
                                }
                            }
                            const std::array<double, 10> backtracking_steps{
                                1.0, 0.5, 0.25, 0.125, 0.0625,
                                0.03125, 0.015625, 0.0078125,
                                0.00390625, 0.001953125};
                            for (int polish_iteration = 1;
                                 polish_iteration <=
                                     options_.max_economic_balance_polish_iterations &&
                                 !economic_polish_budget_exhausted();
                                 ++polish_iteration) {
                                output.economic_balance_polish_iterations =
                                    polish_iteration;
                                AcState raw_trial = polished_state;
                                std::vector<double> polish_p_network;
                                std::vector<double> polish_q_network;
                                network_injections(
                                    polish_ybus, raw_trial.vm, raw_trial.va,
                                    polish_p_network, polish_q_network);
                                std::vector<double> p_spec(nb, 0.0);
                                std::vector<double> q_spec(nb, 0.0);
                                const auto rebuild_specs = [&]() {
                                    std::fill(p_spec.begin(), p_spec.end(), 0.0);
                                    std::fill(q_spec.begin(), q_spec.end(), 0.0);
                                    for (int generator = 0; generator < ng;
                                         ++generator) {
                                        if (!active[generator]) {
                                            continue;
                                        }
                                        const int bus =
                                            data_.generators[generator].bus;
                                        p_spec[bus] += raw_trial.pg[generator];
                                        q_spec[bus] += raw_trial.qg[generator];
                                    }
                                    for (int load = 0;
                                         load < static_cast<int>(data_.loads.size());
                                         ++load) {
                                        const int bus = data_.loads[load].bus;
                                        p_spec[bus] -=
                                            data_.loads[load].pd_nominal *
                                            raw_trial.demand_factor[load];
                                        q_spec[bus] -=
                                            data_.loads[load].qd_nominal *
                                            raw_trial.demand_factor[load];
                                    }
                                };
                                rebuild_specs();

                                // The angle Jacobian omits one active equation
                                // per connected component.  Absorb that exact
                                // component loss mismatch into source-authorized
                                // corrective generation before applying the
                                // cached angle correction.
                                for (int component = 0;
                                     component < static_cast<int>(components.size());
                                     ++component) {
                                    const auto& generators =
                                        component_generators[component];
                                    if (generators.empty()) {
                                        continue;
                                    }
                                    double mismatch = 0.0;
                                    double current_generation = 0.0;
                                    double total_lower = 0.0;
                                    double total_upper = 0.0;
                                    for (int bus : components[component]) {
                                        mismatch += polish_p_network[bus] -
                                            p_spec[bus];
                                    }
                                    for (int generator : generators) {
                                        current_generation +=
                                            raw_trial.pg[generator];
                                        total_lower += p_lower[generator];
                                        total_upper += p_upper[generator];
                                    }
                                    const double target = std::clamp(
                                        current_generation + mismatch,
                                        total_lower, total_upper);
                                    static_cast<void>(allocate_total(
                                        generators, p_lower, p_upper,
                                        polished_state.pg, target,
                                        raw_trial.pg));
                                }

                                // Match local reactive injections with Qg where
                                // available.  The voltage correction handles the
                                // remaining network-wide reactive mismatch.
                                rebuild_specs();
                                for (int bus = 0; bus < nb; ++bus) {
                                    if (active_at_bus[bus].empty()) {
                                        continue;
                                    }
                                    double current_qg = 0.0;
                                    double total_lower = 0.0;
                                    double total_upper = 0.0;
                                    for (int generator : active_at_bus[bus]) {
                                        current_qg += raw_trial.qg[generator];
                                        total_lower += q_lower[generator];
                                        total_upper += q_upper[generator];
                                    }
                                    const double q_mismatch =
                                        polish_q_network[bus] - q_spec[bus];
                                    const double target = std::clamp(
                                        current_qg + q_mismatch,
                                        total_lower, total_upper);
                                    static_cast<void>(allocate_total(
                                        active_at_bus[bus], q_lower, q_upper,
                                        polished_state.qg, target,
                                        raw_trial.qg));
                                }
                                rebuild_specs();

                                bool corrected =
                                    predictor_cache_->apply_correction(
                                        data_, p_spec, q_spec,
                                        polish_p_network, polish_q_network,
                                        raw_trial.vm, raw_trial.va);
                                // Keep the decoupled factors as a conservative
                                // fallback if the coupled solve is unavailable.
                                // Every resulting candidate still goes through
                                // the exact nonlinear rebuild and validator.
                                if (!corrected) {
                                    const bool active_corrected =
                                        predictor_cache_
                                            ->apply_active_correction(
                                                data_, p_spec,
                                                polish_p_network,
                                                raw_trial.va);
                                    if (active_corrected) {
                                        network_injections(
                                            polish_ybus, raw_trial.vm,
                                            raw_trial.va, polish_p_network,
                                            polish_q_network);
                                    }
                                    const bool reactive_corrected =
                                        predictor_cache_
                                            ->apply_reactive_correction(
                                                data_, q_spec,
                                                polish_q_network,
                                                raw_trial.vm);
                                    corrected = active_corrected ||
                                        reactive_corrected;
                                }
                                if (!corrected) {
                                    output.economic_balance_polish_trace.push_back({
                                        {"iteration", polish_iteration},
                                        {"selected", false},
                                        {"reason", "cached Jacobian corrections unavailable"},
                                    });
                                    break;
                                }
                                normalize_source_reference_angles(
                                    data_, components, raw_trial.va);

                                bool selected = false;
                                double selected_step = 0.0;
                                double selected_objective = polished_objective;
                                AcState selected_state = polished_state;
                                ValidationReport selected_validation =
                                    polished_validation;
                                bool first_rejection_recorded = false;
                                double first_rejected_step = 0.0;
                                double first_rejected_objective = 0.0;
                                ValidationReport first_rejected_validation;
                                double trial_active_slack =
                                    slack_sum(polished_state.p_delta);
                                double trial_reactive_slack =
                                    slack_sum(polished_state.q_delta);
                                for (double step : backtracking_steps) {
                                    if (economic_polish_budget_exhausted()) {
                                        break;
                                    }
                                    if (step < 1.0) {
                                        ++output
                                            .economic_balance_polish_backtracking_attempts;
                                    }
                                    AcState candidate = polished_state;
                                    for (int bus = 0; bus < nb; ++bus) {
                                        candidate.vm[bus] += step *
                                            (raw_trial.vm[bus] -
                                             polished_state.vm[bus]);
                                        candidate.va[bus] += step *
                                            (raw_trial.va[bus] -
                                             polished_state.va[bus]);
                                    }
                                    for (int generator = 0; generator < ng;
                                         ++generator) {
                                        candidate.pg[generator] += step *
                                            (raw_trial.pg[generator] -
                                             polished_state.pg[generator]);
                                        candidate.qg[generator] += step *
                                            (raw_trial.qg[generator] -
                                             polished_state.qg[generator]);
                                    }
                                    normalize_source_reference_angles(
                                        data_, components, candidate.va);
                                    const double candidate_objective =
                                        rebuild_contingency_state_derived_fields(
                                            data_, base_state_, commitment_,
                                            *contingency, candidate);
                                    const auto candidate_validation =
                                        validate_state(
                                            data_, ModelMode::ContingencySoft,
                                            candidate, commitment_,
                                            direct_context);
                                    if (candidate_validation.max_residual >
                                            options_.validation_tolerance &&
                                        !first_rejection_recorded) {
                                        first_rejection_recorded = true;
                                        first_rejected_step = step;
                                        first_rejected_objective =
                                            candidate_objective;
                                        first_rejected_validation =
                                            candidate_validation;
                                        exact_newton_reactive_trigger =
                                            candidate_validation.worst_category ==
                                            "reactive_balance";
                                    }
                                    trial_active_slack =
                                        slack_sum(candidate.p_delta);
                                    trial_reactive_slack =
                                        slack_sum(candidate.q_delta);
                                    if (candidate_validation.max_residual <=
                                            options_.validation_tolerance &&
                                        candidate_objective >
                                            polished_objective + 1e-9) {
                                        selected = true;
                                        selected_step = step;
                                        selected_objective =
                                            candidate_objective;
                                        selected_state = std::move(candidate);
                                        selected_validation =
                                            candidate_validation;
                                        break;
                                    }
                                }
                                nlohmann::json polish_trace = {
                                    {"iteration", polish_iteration},
                                    {"selected", selected},
                                    {"selected_step", selected_step},
                                    {"objective_before", polished_objective},
                                    {"objective_after", selected_objective},
                                    {"active_slack_after", trial_active_slack},
                                    {"reactive_slack_after", trial_reactive_slack},
                                    {"validation", selected_validation.to_json()},
                                };
                                if (first_rejection_recorded) {
                                    polish_trace["first_rejected_step"] =
                                        first_rejected_step;
                                    polish_trace["first_rejected_objective"] =
                                        first_rejected_objective;
                                    polish_trace["first_rejected_validation"] =
                                        first_rejected_validation.to_json();
                                }
                                output.economic_balance_polish_trace.push_back(
                                    std::move(polish_trace));
                                if (!selected) {
                                    break;
                                }
                                polished_state = std::move(selected_state);
                                polished_objective = selected_objective;
                                polished_validation = selected_validation;
                                if (slack_sum(polished_state.p_delta) +
                                        slack_sum(polished_state.q_delta) <=
                                    1e-7) {
                                    break;
                                }
                            }
                            apply_local_injection_projection();

                            // Feasibility repair deliberately stops once
                            // every nodal Q mismatch fits inside the source
                            // 0.5-p.u. soft band.  Economically, thousands of
                            // individually admissible mismatches can still
                            // carry a large aggregate penalty.  Reuse the
                            // already-factorized reactive Jacobian to drive
                            // the largest in-band residuals toward zero.  An
                            // active-angle correction restores the changed P
                            // injections, and only a fully rebuilt, validated,
                            // objective-improving state can replace the
                            // incumbent.
                            for (int zero_balance_round = 1;
                                 zero_balance_round <= options_
                                     .max_economic_reactive_zero_balance_rounds &&
                                 polished_objective <= options_
                                     .economic_reactive_zero_balance_objective_threshold &&
                                 !economic_polish_budget_exhausted();
                                 ++zero_balance_round) {
                                const YRows current_ybus = build_ybus(
                                    data_, outaged_branch, &polished_state);
                                std::vector<double> current_p_network;
                                std::vector<double> current_q_network;
                                network_injections(
                                    current_ybus, polished_state.vm,
                                    polished_state.va, current_p_network,
                                    current_q_network);
                                std::vector<double> current_p_spec(
                                    static_cast<std::size_t>(nb), 0.0);
                                std::vector<double> current_q_balance =
                                    current_q_network;
                                for (int generator = 0; generator < ng;
                                     ++generator) {
                                    if (!active[generator]) {
                                        continue;
                                    }
                                    const int bus =
                                        data_.generators[generator].bus;
                                    current_p_spec[bus] +=
                                        polished_state.pg[generator];
                                    current_q_balance[bus] -=
                                        polished_state.qg[generator];
                                }
                                for (int load = 0;
                                     load < static_cast<int>(
                                         data_.loads.size());
                                     ++load) {
                                    const int bus = data_.loads[load].bus;
                                    const double factor =
                                        polished_state.demand_factor[load];
                                    current_p_spec[bus] -=
                                        data_.loads[load].pd_nominal * factor;
                                    current_q_balance[bus] +=
                                        data_.loads[load].qd_nominal * factor;
                                }
                                double q_slack_before = 0.0;
                                for (double value : current_q_balance) {
                                    q_slack_before += std::abs(value);
                                }
                                if (q_slack_before <= options_
                                        .economic_reactive_zero_balance_trigger_slack) {
                                    break;
                                }

                                // The voltage direction is independent of the
                                // damping scalar.  Build/factor it once, then
                                // interpolate candidate voltages below.  Keep
                                // the bounded 256-row active set: target_band=0
                                // otherwise promotes thousands of harmless
                                // floating-point residuals into an unbounded
                                // dense least-squares problem.
                                std::vector<double> proposed_vm =
                                    polished_state.vm;
                                if (!predictor_cache_
                                        ->apply_local_reactive_least_squares(
                                            data_, current_q_balance,
                                            proposed_vm, 1.0, true, 0.0)) {
                                    break;
                                }

                                bool selected_zero_balance = false;
                                double selected_zero_balance_step = 0.0;
                                double selected_zero_balance_objective =
                                    polished_objective;
                                AcState selected_zero_balance_state =
                                    polished_state;
                                ValidationReport selected_zero_balance_validation =
                                    polished_validation;
                                constexpr std::array<double, 8>
                                    kZeroBalanceDamping{
                                        1.0, 0.5, 0.25, 0.125,
                                        0.0625, 0.03125, 0.015625,
                                        0.0078125};
                                for (double damping : kZeroBalanceDamping) {
                                    if (economic_polish_budget_exhausted()) {
                                        break;
                                    }
                                    AcState candidate = polished_state;
                                    for (int bus = 0; bus < nb; ++bus) {
                                        candidate.vm[bus] = std::clamp(
                                            polished_state.vm[bus] +
                                                damping *
                                                    (proposed_vm[bus] -
                                                     polished_state.vm[bus]),
                                            data_.buses[bus].vmin,
                                            data_.buses[bus].vmax);
                                    }
                                    std::vector<double> candidate_p_network;
                                    std::vector<double> candidate_q_network;
                                    network_injections(
                                        current_ybus, candidate.vm,
                                        candidate.va, candidate_p_network,
                                        candidate_q_network);
                                    static_cast<void>(candidate_q_network);
                                    if (!predictor_cache_
                                            ->apply_active_correction(
                                                data_, current_p_spec,
                                                candidate_p_network,
                                                candidate.va, 1.0)) {
                                        continue;
                                    }
                                    normalize_source_reference_angles(
                                        data_, components, candidate.va);
                                    const double candidate_objective =
                                        rebuild_contingency_state_derived_fields(
                                            data_, base_state_, commitment_,
                                            *contingency, candidate);
                                    const auto candidate_validation =
                                        validate_state(
                                            data_,
                                            ModelMode::ContingencySoft,
                                            candidate, commitment_,
                                            direct_context);
                                    if (candidate_validation.max_residual <=
                                            options_.validation_tolerance &&
                                        candidate_objective >
                                            selected_zero_balance_objective +
                                                1e-9) {
                                        selected_zero_balance = true;
                                        selected_zero_balance_step = damping;
                                        selected_zero_balance_objective =
                                            candidate_objective;
                                        selected_zero_balance_state =
                                            std::move(candidate);
                                        selected_zero_balance_validation =
                                            candidate_validation;
                                    }
                                }
                                output.economic_balance_polish_trace.push_back({
                                    {"phase", "reactive_zero_balance"},
                                    {"round", zero_balance_round},
                                    {"selected", selected_zero_balance},
                                    {"selected_step",
                                     selected_zero_balance_step},
                                    {"objective_before", polished_objective},
                                    {"objective_after",
                                     selected_zero_balance_objective},
                                    {"reactive_slack_before", q_slack_before},
                                    {"validation",
                                     selected_zero_balance_validation.to_json()},
                                });
                                if (!selected_zero_balance) {
                                    break;
                                }
                                polished_state =
                                    std::move(selected_zero_balance_state);
                                polished_objective =
                                    selected_zero_balance_objective;
                                polished_validation =
                                    selected_zero_balance_validation;
                                apply_local_injection_projection();
                            }
                        }
                        // The cached Jacobian is deliberately cheap, but a
                        // small subset of severe outages can stall at an
                        // uncontrolled shunt/PQ bus.  For only that economic
                        // tail, refresh the actual sparse AC Jacobian and
                        // solve the zero-balance equations for a few bounded
                        // Newton steps.  This is still only a candidate: all
                        // generator/load/ramp limits are reconstructed below,
                        // and the complete nonlinear contingency validator is
                        // the sole acceptance gate.
                        const bool merit_exact_newton_target =
                            outaged_branch >= 0 &&
                            output.economic_merit_dispatch_applied;
                        if (options_.economic_exact_newton_rescue &&
                            (merit_exact_newton_target ||
                             exact_newton_reactive_trigger ||
                             slack_sum(polished_state.p_delta) +
                                     slack_sum(polished_state.q_delta) >
                                 options_
                                     .economic_linearized_trigger_slack) &&
                            (merit_exact_newton_target ||
                             polished_objective <= options_
                                 .economic_exact_newton_objective_threshold) &&
                            options_.economic_exact_newton_max_iterations > 0 &&
                            !economic_polish_budget_exhausted()) {
                            output.economic_exact_newton_attempted = true;
                            const auto exact_newton_start =
                                std::chrono::steady_clock::now();
                            // Generator outages need the feasibility-polished
                            // controls: their component redispatch can land on
                            // several simultaneous corrective limits, and the
                            // existing route is already robust there.  A branch
                            // outage leaves every generator available, so use
                            // the preserved economic target and retain the
                            // secure predictor only as the voltage start.
                            const AcState& exact_control_reference =
                                outaged_branch >= 0
                                ? economic_target_state : polished_state;
                            const YRows exact_ybus = build_ybus(
                                data_, outaged_branch,
                                &exact_control_reference);
                            std::vector<double> exact_p_spec(nb, 0.0);
                            std::vector<double> exact_q_spec(nb, 0.0);
                            std::vector<double> exact_load_p(nb, 0.0);
                            std::vector<double> exact_load_q(nb, 0.0);
                            for (int generator = 0; generator < ng;
                                 ++generator) {
                                if (!active[generator]) {
                                    continue;
                                }
                                const int bus =
                                    data_.generators[generator].bus;
                                exact_p_spec[bus] +=
                                    exact_control_reference.pg[generator];
                                exact_q_spec[bus] +=
                                    exact_control_reference.qg[generator];
                            }
                            for (int load = 0;
                                 load < static_cast<int>(data_.loads.size());
                                 ++load) {
                                const int bus = data_.loads[load].bus;
                                const double load_p =
                                    data_.loads[load].pd_nominal *
                                    exact_control_reference.demand_factor[load];
                                const double load_q =
                                    data_.loads[load].qd_nominal *
                                    exact_control_reference.demand_factor[load];
                                exact_load_p[bus] += load_p;
                                exact_load_q[bus] += load_q;
                                exact_p_spec[bus] -= load_p;
                                exact_q_spec[bus] -= load_q;
                            }

                            std::vector<double> exact_p_network;
                            std::vector<double> exact_q_network;
                            network_injections(
                                exact_ybus, polished_state.vm,
                                polished_state.va, exact_p_network,
                                exact_q_network);
                            std::vector<double> component_mismatch(
                                components.size(), 0.0);
                            for (int bus = 0; bus < nb; ++bus) {
                                component_mismatch[component_of[bus]] +=
                                    exact_p_network[bus] -
                                    exact_p_spec[bus];
                            }
                            std::vector<double> active_slack_weights(
                                static_cast<std::size_t>(nb), 0.0);
                            for (int bus = 0; bus < nb; ++bus) {
                                const bool upward =
                                    component_mismatch[component_of[bus]] >=
                                    0.0;
                                for (int generator : active_at_bus[bus]) {
                                    active_slack_weights[bus] += upward
                                        ? std::max(
                                              0.0, p_upper[generator] -
                                                  exact_control_reference.pg[generator])
                                        : std::max(
                                              0.0,
                                              exact_control_reference.pg[generator] -
                                                  p_lower[generator]);
                                }
                            }
                            for (int component = 0;
                                 component <
                                     static_cast<int>(components.size());
                                 ++component) {
                                double weight_sum = 0.0;
                                for (int bus : components[component]) {
                                    weight_sum += active_slack_weights[bus];
                                }
                                if (weight_sum > 1e-12) {
                                    continue;
                                }
                                for (int bus : components[component]) {
                                    if (!active_at_bus[bus].empty()) {
                                        active_slack_weights[bus] = 1.0;
                                        break;
                                    }
                                }
                            }

                            std::vector<double> exact_vm =
                                polished_state.vm;
                            std::vector<double> exact_va =
                                polished_state.va;
                            std::vector<bool> exact_pq = pq;
                            NewtonResult exact_newton;
                            bool exact_active_set_complete = false;
                            constexpr int kExactQActiveSetPasses = 12;
                            for (int active_set_pass = 0;
                                 active_set_pass < kExactQActiveSetPasses;
                                 ++active_set_pass) {
                                exact_newton = run_distributed_active_newton(
                                    data_, exact_ybus, slack, exact_pq,
                                    component_of, active_slack_weights,
                                    exact_p_spec, exact_q_spec,
                                    options_
                                        .economic_exact_newton_max_iterations,
                                    std::min(
                                        1e-8,
                                        options_.validation_tolerance * 0.1),
                                    exact_vm, exact_va);
                                output.economic_exact_newton_iterations +=
                                    exact_newton.iterations;
                                if (!exact_newton.converged) {
                                    output.economic_exact_newton_failure_reason =
                                        exact_newton.failure_reason;
                                    break;
                                }
                                network_injections(
                                    exact_ybus, exact_vm, exact_va,
                                    exact_p_network, exact_q_network);
                                bool added_q_limit = false;
                                for (int bus = 0; bus < nb; ++bus) {
                                    if (active_at_bus[bus].empty() ||
                                        exact_pq[bus]) {
                                        continue;
                                    }
                                    double total_lower = 0.0;
                                    double total_upper = 0.0;
                                    for (int generator : active_at_bus[bus]) {
                                        total_lower += q_lower[generator];
                                        total_upper += q_upper[generator];
                                    }
                                    const double required_q =
                                        exact_q_network[bus] +
                                        exact_load_q[bus];
                                    if (required_q <
                                            total_lower - kAllocationTolerance ||
                                        required_q >
                                            total_upper + kAllocationTolerance) {
                                        exact_q_spec[bus] = std::clamp(
                                            required_q, total_lower,
                                            total_upper) - exact_load_q[bus];
                                        exact_pq[bus] = true;
                                        added_q_limit = true;
                                        ++output
                                            .economic_exact_newton_q_limit_switches;
                                    }
                                }
                                if (!added_q_limit) {
                                    exact_active_set_complete = true;
                                    break;
                                }
                            }
                            output.economic_exact_newton_converged =
                                exact_active_set_complete;
                            if (exact_newton.converged &&
                                !exact_active_set_complete &&
                                output.economic_exact_newton_failure_reason.empty()) {
                                output.economic_exact_newton_failure_reason =
                                    "exact Newton reactive active-set limit";
                            }
                            if (exact_active_set_complete) {
                                AcState exact_state = exact_control_reference;
                                exact_state.vm = std::move(exact_vm);
                                exact_state.va = std::move(exact_va);
                                network_injections(
                                    exact_ybus, exact_state.vm,
                                    exact_state.va, exact_p_network,
                                    exact_q_network);
                                bool controls_reconstructed = true;
                                for (int bus = 0; bus < nb; ++bus) {
                                    if (active_at_bus[bus].empty()) {
                                        continue;
                                    }
                                    const double required_p =
                                        exact_p_network[bus] +
                                        exact_load_p[bus];
                                    if (!allocate_total(
                                            active_at_bus[bus], p_lower,
                                            p_upper, exact_control_reference.pg,
                                            required_p, exact_state.pg)) {
                                        double total_lower = 0.0;
                                        double total_upper = 0.0;
                                        for (int generator :
                                             active_at_bus[bus]) {
                                            total_lower += p_lower[generator];
                                            total_upper += p_upper[generator];
                                        }
                                        output.economic_exact_newton_failure_reason =
                                            "exact Newton active controls exceed "
                                            "source bounds at " +
                                            data_.buses[bus].source_key +
                                            "; target=" +
                                            std::to_string(required_p) +
                                            "; range=[" +
                                            std::to_string(total_lower) + "," +
                                            std::to_string(total_upper) + "]";
                                        controls_reconstructed = false;
                                        break;
                                    }
                                    const double required_q =
                                        exact_q_network[bus] +
                                        exact_load_q[bus];
                                    if (!allocate_total(
                                            active_at_bus[bus], q_lower,
                                            q_upper, exact_control_reference.qg,
                                            required_q, exact_state.qg)) {
                                        double total_lower = 0.0;
                                        double total_upper = 0.0;
                                        for (int generator :
                                             active_at_bus[bus]) {
                                            total_lower += q_lower[generator];
                                            total_upper += q_upper[generator];
                                        }
                                        output.economic_exact_newton_failure_reason =
                                            "exact Newton reactive controls exceed "
                                            "source bounds at " +
                                            data_.buses[bus].source_key +
                                            "; target=" +
                                            std::to_string(required_q) +
                                            "; range=[" +
                                            std::to_string(total_lower) + "," +
                                            std::to_string(total_upper) + "]";
                                        controls_reconstructed = false;
                                        break;
                                    }
                                }
                                if (controls_reconstructed) {
                                    normalize_source_reference_angles(
                                        data_, components, exact_state.va);
                                    const std::array<double, 11>
                                        exact_line_search_steps{
                                            1.0, 0.75, 0.5, 0.25, 0.125,
                                            0.0625, 0.03125, 0.015625,
                                            0.0078125, 0.00390625,
                                            0.001953125};
                                    for (double step :
                                         exact_line_search_steps) {
                                        AcState candidate =
                                            interpolate_corrective_candidate(
                                                polished_state, exact_state,
                                                step);
                                        normalize_source_reference_angles(
                                            data_, components, candidate.va);
                                        const double candidate_objective =
                                            rebuild_contingency_state_derived_fields(
                                                data_, base_state_, commitment_,
                                                *contingency, candidate);
                                        const auto candidate_validation =
                                            validate_state(
                                                data_,
                                                ModelMode::ContingencySoft,
                                                candidate, commitment_,
                                                direct_context);
                                        output.economic_exact_newton_objective =
                                            candidate_objective;
                                        output.economic_exact_newton_validation =
                                            candidate_validation;
                                        if (candidate_validation.max_residual <=
                                                options_.validation_tolerance &&
                                            candidate_objective >
                                                polished_objective + 1e-9) {
                                            output.economic_exact_newton_selected =
                                                true;
                                            polished_state =
                                                std::move(candidate);
                                            polished_objective =
                                                candidate_objective;
                                            polished_validation =
                                                candidate_validation;
                                            break;
                                        }
                                    }

                                    // A full zero-balance Newton point can be
                                    // blocked by one voltage or thermal bound,
                                    // forcing the direct line search to retain
                                    // nearly all of the incumbent imbalance.
                                    // Repair those explicit security defects
                                    // around the full Newton point with the
                                    // existing source-bounded linear model,
                                    // then subject every blend to the unchanged
                                    // nonlinear contingency validator.
                                    const double security_repair_seconds =
                                        bounded_economic_solver_seconds(1.0);
                                    if (security_repair_seconds > 1e-3) {
                                        AcState security_reference =
                                            exact_state;
                                        int projected_voltage_count = 0;
                                        double maximum_voltage_projection =
                                            0.0;
                                        for (int bus = 0; bus < nb; ++bus) {
                                            const double projected =
                                                std::clamp(
                                                    security_reference.vm[bus],
                                                    data_.buses[bus].vmin,
                                                    data_.buses[bus].vmax);
                                            const double movement = std::abs(
                                                projected -
                                                security_reference.vm[bus]);
                                            if (movement > 1e-12) {
                                                ++projected_voltage_count;
                                                maximum_voltage_projection =
                                                    std::max(
                                                        maximum_voltage_projection,
                                                        movement);
                                                security_reference.vm[bus] =
                                                    projected;
                                            }
                                        }
                                        rebuild_contingency_state_derived_fields(
                                            data_, base_state_, commitment_,
                                            *contingency,
                                            security_reference);
                                        const auto try_security_repair = [&] (
                                            const ActiveFeasibilityRepairResult&
                                                repair,
                                            nlohmann::json& trace) {
                                            bool first_rejection_recorded = false;
                                            if (repair.success) {
                                                const std::array<double, 11>
                                                    repair_line_search_steps{
                                                        1.0, 0.75, 0.5, 0.25,
                                                        0.125, 0.0625, 0.03125,
                                                        0.015625, 0.0078125,
                                                        0.00390625, 0.001953125};
                                                for (double step :
                                                     repair_line_search_steps) {
                                                    AcState candidate =
                                                        interpolate_corrective_candidate(
                                                            polished_state,
                                                            repair.state,
                                                            step);
                                                    normalize_source_reference_angles(
                                                        data_, components,
                                                        candidate.va);
                                                    const double
                                                        candidate_objective =
                                                            rebuild_contingency_state_derived_fields(
                                                                data_, base_state_,
                                                                commitment_,
                                                                *contingency,
                                                                candidate);
                                                    const auto
                                                        candidate_validation =
                                                            validate_state(
                                                                data_,
                                                                ModelMode::ContingencySoft,
                                                                candidate,
                                                                commitment_,
                                                                direct_context);
                                                    if (!first_rejection_recorded &&
                                                        candidate_validation
                                                                .max_residual >
                                                            options_
                                                                .validation_tolerance) {
                                                        first_rejection_recorded =
                                                            true;
                                                        trace[
                                                            "first_rejected_step"] =
                                                            step;
                                                        trace[
                                                            "first_rejected_objective"] =
                                                            candidate_objective;
                                                        trace[
                                                            "first_rejected_validation"] =
                                                            candidate_validation
                                                                .to_json();
                                                    }
                                                    if (candidate_validation
                                                                .max_residual <=
                                                            options_
                                                                .validation_tolerance &&
                                                        candidate_objective >
                                                            polished_objective +
                                                                1e-9) {
                                                        polished_state =
                                                            std::move(candidate);
                                                        polished_objective =
                                                            candidate_objective;
                                                        polished_validation =
                                                            candidate_validation;
                                                        trace["selected"] = true;
                                                        trace["selected_step"] =
                                                            step;
                                                        trace["validation"] =
                                                            candidate_validation
                                                                .to_json();
                                                        break;
                                                    }
                                                }
                                            }
                                            trace["objective_after"] =
                                                polished_objective;
                                            output.economic_balance_polish_trace
                                                .push_back(std::move(trace));
                                        };

                                        // Preserve the fast ordinary security
                                        // repair first. If its full nonlinear
                                        // point remains blocked only by a
                                        // shrinking branch/angle residual,
                                        // relinearize around that point a few
                                        // times. This is substantially smaller
                                        // than Phase I and retains the exact
                                        // Newton point's zero-balance economics.
                                        AcState security_iteration_reference =
                                            security_reference;
                                        double prior_full_residual =
                                            std::numeric_limits<double>::infinity();
                                        constexpr int
                                            kMaximumSequentialSecurityRepairs = 4;
                                        for (int security_round = 1;
                                             security_round <=
                                                 kMaximumSequentialSecurityRepairs &&
                                             !economic_polish_budget_exhausted();
                                             ++security_round) {
                                            const double round_seconds =
                                                security_round == 1
                                                ? security_repair_seconds
                                                : bounded_economic_solver_seconds(
                                                      0.75);
                                            if (round_seconds <= 1e-3) {
                                                break;
                                            }
                                            const auto security_repair =
                                                solve_linearized_active_feasibility_repair(
                                                    data_,
                                                    security_iteration_reference,
                                                    commitment_, *direct_context,
                                                    0.5, 0.05, round_seconds,
                                                    // The full Newton point can
                                                    // legitimately overshoot a
                                                    // source voltage bound by
                                                    // more than the old
                                                    // 0.02-p.u. radius. A wider
                                                    // solver trust region lets
                                                    // the LP move it back inside
                                                    // the unchanged source bound;
                                                    // the nonlinear validator
                                                    // remains the acceptance gate.
                                                    0.10, true, true, false,
                                                    false);
                                            nlohmann::json security_trace = {
                                                {"phase",
                                                 security_round == 1
                                                     ? "exact_newton_security_repair"
                                                     : "sequential_exact_newton_security_repair"},
                                                {"round", security_round},
                                                {"repair",
                                                 security_repair.to_json(false)},
                                                {"selected", false},
                                                {"selected_step", 0.0},
                                                {"objective_before",
                                                 polished_objective},
                                                {"projected_voltage_count",
                                                 projected_voltage_count},
                                                {"maximum_voltage_projection",
                                                 maximum_voltage_projection},
                                            };
                                            bool continue_sequential_repair =
                                                false;
                                            AcState next_security_reference;
                                            if (security_repair.success) {
                                                next_security_reference =
                                                    security_repair.state;
                                                normalize_source_reference_angles(
                                                    data_, components,
                                                    next_security_reference.va);
                                                const double full_objective =
                                                    rebuild_contingency_state_derived_fields(
                                                        data_, base_state_,
                                                        commitment_,
                                                        *contingency,
                                                        next_security_reference);
                                                const auto full_validation =
                                                    validate_state(
                                                        data_,
                                                        ModelMode::ContingencySoft,
                                                        next_security_reference,
                                                        commitment_,
                                                        direct_context);
                                                security_trace[
                                                    "full_candidate_objective"] =
                                                    full_objective;
                                                security_trace[
                                                    "full_candidate_validation"] =
                                                    full_validation.to_json();
                                                const bool security_only_block =
                                                    full_validation.worst_category ==
                                                        "flow_limit" ||
                                                    full_validation.worst_category ==
                                                        "angle_limit" ||
                                                    (full_validation.worst_category ==
                                                         "variable_bound" &&
                                                     full_validation.worst_identity
                                                             .rfind("branch:", 0) ==
                                                         0);
                                                continue_sequential_repair =
                                                    full_objective >
                                                        polished_objective + 1e-9 &&
                                                    full_validation.max_residual >
                                                        options_
                                                            .validation_tolerance &&
                                                    security_only_block &&
                                                    full_validation.max_residual +
                                                            1e-8 <
                                                        prior_full_residual;
                                                prior_full_residual =
                                                    full_validation.max_residual;
                                            }
                                            try_security_repair(
                                                security_repair, security_trace);
                                            if (!continue_sequential_repair) {
                                                break;
                                            }
                                            security_iteration_reference =
                                                std::move(
                                                    next_security_reference);
                                        }

                                        // A projected exact-Newton point is the
                                        // distinctive failure mode in which the
                                        // ordinary repair can retain costly P/Q
                                        // slack. Only for that case, spend a
                                        // second short Phase-I solve. Selection
                                        // remains strictly objective-monotone
                                        // and independently nonlinear-verified.
                                        if (projected_voltage_count > 0 &&
                                            !economic_polish_budget_exhausted()) {
                                            const double balance_seconds =
                                                bounded_economic_solver_seconds(
                                                    1.0);
                                            if (balance_seconds > 1e-3) {
                                                const auto balance_repair =
                                                    solve_linearized_active_feasibility_repair(
                                                        data_, security_reference,
                                                        commitment_,
                                                        *direct_context,
                                                        0.5, 0.05,
                                                        balance_seconds,
                                                        0.10, true, true,
                                                        false, true);
                                                nlohmann::json balance_trace = {
                                                    {"phase",
                                                     "projected_exact_newton_balance_phase_one"},
                                                    {"repair",
                                                     balance_repair.to_json(false)},
                                                    {"selected", false},
                                                    {"selected_step", 0.0},
                                                    {"objective_before",
                                                     polished_objective},
                                                    {"projected_voltage_count",
                                                     projected_voltage_count},
                                                    {"maximum_voltage_projection",
                                                     maximum_voltage_projection},
                                                };
                                                try_security_repair(
                                                    balance_repair,
                                                    balance_trace);
                                            }
                                        }

                                    }
                                }
                            }
                            output.economic_exact_newton_wall_seconds =
                                std::chrono::duration<double>(
                                    std::chrono::steady_clock::now() -
                                    exact_newton_start).count();
                        }
                        for (int linearized_round = 1;
                             linearized_round <=
                                 options_.max_economic_linearized_polish_rounds &&
                             polished_objective <=
                                 options_.economic_linearized_objective_threshold &&
                             slack_sum(polished_state.p_delta) +
                                     slack_sum(polished_state.q_delta) >
                                 options_.economic_linearized_trigger_slack &&
                             !economic_polish_budget_exhausted();
                             ++linearized_round) {
                            const double repair_seconds =
                                bounded_economic_solver_seconds(
                                    options_.economic_linearized_polish_seconds);
                            if (repair_seconds <= 1e-3) {
                                output.economic_balance_polish_time_limit_reached =
                                    true;
                                break;
                            }
                            const auto repair =
                                solve_linearized_active_feasibility_repair(
                                    data_, polished_state, commitment_,
                                    *direct_context,
                                    0.5, 0.05,
                                    repair_seconds,
                                    0.01, true, true, false, true);
                            nlohmann::json repair_trace = {
                                {"phase", "linearized_zero_balance"},
                                {"round", linearized_round},
                                {"repair", repair.to_json(false)},
                                {"selected", false},
                                {"selected_step", 0.0},
                                {"objective_before", polished_objective},
                            };
                            if (!repair.success) {
                                repair_trace["objective_after"] =
                                    polished_objective;
                                output.economic_balance_polish_trace.push_back(
                                    std::move(repair_trace));
                                break;
                            }
                            const std::array<double, 11> line_search_steps{
                                1.0, 0.75, 0.5, 0.25, 0.125, 0.0625,
                                0.03125, 0.015625, 0.0078125,
                                0.00390625, 0.001953125};
                            bool selected = false;
                            double selected_step = 0.0;
                            AcState selected_state = polished_state;
                            double selected_objective = polished_objective;
                            ValidationReport selected_validation =
                                polished_validation;
                            ValidationReport first_rejected_validation;
                            bool first_rejection_recorded = false;
                            double first_rejected_step = 0.0;
                            for (double step : line_search_steps) {
                                if (economic_polish_budget_exhausted()) {
                                    break;
                                }
                                if (step < 1.0) {
                                    ++output
                                        .economic_balance_polish_backtracking_attempts;
                                }
                                AcState candidate = polished_state;
                                for (int bus = 0; bus < nb; ++bus) {
                                    candidate.vm[bus] += step *
                                        (repair.state.vm[bus] -
                                         polished_state.vm[bus]);
                                    candidate.va[bus] += step *
                                        (repair.state.va[bus] -
                                         polished_state.va[bus]);
                                }
                                for (int generator = 0; generator < ng;
                                     ++generator) {
                                    candidate.pg[generator] += step *
                                        (repair.state.pg[generator] -
                                         polished_state.pg[generator]);
                                    candidate.qg[generator] += step *
                                        (repair.state.qg[generator] -
                                         polished_state.qg[generator]);
                                }
                                for (int load = 0;
                                     load < static_cast<int>(data_.loads.size());
                                     ++load) {
                                    candidate.demand_factor[load] += step *
                                        (repair.state.demand_factor[load] -
                                         polished_state.demand_factor[load]);
                                }
                                normalize_source_reference_angles(
                                    data_, components, candidate.va);
                                const double candidate_objective =
                                    rebuild_contingency_state_derived_fields(
                                        data_, base_state_, commitment_,
                                        *contingency, candidate);
                                const auto candidate_validation =
                                    validate_state(
                                        data_, ModelMode::ContingencySoft,
                                        candidate, commitment_,
                                        direct_context);
                                if (candidate_validation.max_residual >
                                        options_.validation_tolerance &&
                                    !first_rejection_recorded) {
                                    first_rejection_recorded = true;
                                    first_rejected_step = step;
                                    first_rejected_validation =
                                        candidate_validation;
                                }
                                if (candidate_validation.max_residual <=
                                        options_.validation_tolerance &&
                                    candidate_objective >
                                        polished_objective + 1e-9) {
                                    selected = true;
                                    selected_step = step;
                                    selected_state = std::move(candidate);
                                    selected_objective = candidate_objective;
                                    selected_validation =
                                        candidate_validation;
                                    break;
                                }
                            }
                            repair_trace["selected"] = selected;
                            repair_trace["selected_step"] = selected_step;
                            repair_trace["objective_after"] =
                                selected_objective;
                            repair_trace["active_slack_after"] =
                                slack_sum(selected_state.p_delta);
                            repair_trace["reactive_slack_after"] =
                                slack_sum(selected_state.q_delta);
                            repair_trace["validation"] =
                                selected_validation.to_json();
                            if (first_rejection_recorded) {
                                repair_trace["first_rejected_step"] =
                                    first_rejected_step;
                                repair_trace["first_rejected_validation"] =
                                    first_rejected_validation.to_json();
                            }
                            output.economic_balance_polish_trace.push_back(
                                std::move(repair_trace));
                            if (!selected) {
                                break;
                            }
                            polished_state = std::move(selected_state);
                            polished_objective = selected_objective;
                            polished_validation = selected_validation;
                        }
                        // Phase II: from the verified Phase-I point, solve a
                        // short compact LP whose coefficients are the actual
                        // source PWL generator/load marginal economics plus
                        // the source balance-slack penalties.  The complete
                        // primal point is supplied to HiGHS, and every trial
                        // is rebuilt with the nonlinear GO2 objective and all
                        // source contingency constraints before acceptance.
                        for (int phase_two_round = 1;
                             phase_two_round <= options_
                                 .max_economic_linearized_phase_two_rounds &&
                             polished_objective <=
                                 options_.economic_linearized_objective_threshold &&
                             !economic_polish_budget_exhausted();
                             ++phase_two_round) {
                            const double phase_two_seconds =
                                bounded_economic_solver_seconds(
                                    options_.economic_linearized_phase_two_seconds);
                            if (phase_two_seconds <= 1e-3) {
                                output.economic_balance_polish_time_limit_reached =
                                    true;
                                break;
                            }
                            const auto phase_two = solve_linearized_ac_seed(
                                data_, polished_state, commitment_, 0.5,
                                direct_context, false, true,
                                phase_two_seconds,
                                true, false, {}, false, 0.01, 0.05, true);
                            nlohmann::json phase_two_trace = {
                                {"phase", "linearized_actual_economic"},
                                {"round", phase_two_round},
                                {"linear_model", phase_two.to_json(false)},
                                {"selected", false},
                                {"selected_step", 0.0},
                                {"objective_before", polished_objective},
                                {"candidates", nlohmann::json::array()},
                            };
                            if (!phase_two.success) {
                                phase_two_trace["objective_after"] =
                                    polished_objective;
                                output.economic_balance_polish_trace.push_back(
                                    std::move(phase_two_trace));
                                break;
                            }
                            const std::array<double, 11> line_search_steps{
                                1.0, 0.75, 0.5, 0.25, 0.125, 0.0625,
                                0.03125, 0.015625, 0.0078125,
                                0.00390625, 0.001953125};
                            bool selected = false;
                            double selected_step = 0.0;
                            AcState selected_state = polished_state;
                            double selected_objective = polished_objective;
                            ValidationReport selected_validation =
                                polished_validation;
                            for (double step : line_search_steps) {
                                if (economic_polish_budget_exhausted()) {
                                    break;
                                }
                                if (step < 1.0) {
                                    ++output
                                        .economic_balance_polish_backtracking_attempts;
                                }
                                AcState candidate = polished_state;
                                for (int bus = 0; bus < nb; ++bus) {
                                    candidate.vm[bus] += step *
                                        (phase_two.state.vm[bus] -
                                         polished_state.vm[bus]);
                                    candidate.va[bus] += step *
                                        (phase_two.state.va[bus] -
                                         polished_state.va[bus]);
                                }
                                for (int generator = 0; generator < ng;
                                     ++generator) {
                                    candidate.pg[generator] += step *
                                        (phase_two.state.pg[generator] -
                                         polished_state.pg[generator]);
                                    candidate.qg[generator] += step *
                                        (phase_two.state.qg[generator] -
                                         polished_state.qg[generator]);
                                }
                                for (int load = 0;
                                     load < static_cast<int>(data_.loads.size());
                                     ++load) {
                                    candidate.demand_factor[load] += step *
                                        (phase_two.state.demand_factor[load] -
                                         polished_state.demand_factor[load]);
                                }
                                normalize_source_reference_angles(
                                    data_, components, candidate.va);
                                const double candidate_objective =
                                    rebuild_contingency_state_derived_fields(
                                        data_, base_state_, commitment_,
                                        *contingency, candidate);
                                const auto candidate_validation =
                                    validate_state(
                                        data_, ModelMode::ContingencySoft,
                                        candidate, commitment_,
                                        direct_context);
                                const bool accepted =
                                    candidate_validation.max_residual <=
                                        options_.validation_tolerance &&
                                    candidate_objective >
                                        selected_objective + 1e-9;
                                phase_two_trace["candidates"].push_back({
                                    {"step", step},
                                    {"objective", candidate_objective},
                                    {"accepted", accepted},
                                    {"validation",
                                     candidate_validation.to_json()},
                                });
                                if (accepted) {
                                    selected = true;
                                    selected_step = step;
                                    selected_state = std::move(candidate);
                                    selected_objective = candidate_objective;
                                    selected_validation =
                                        candidate_validation;
                                }
                            }
                            phase_two_trace["selected"] = selected;
                            phase_two_trace["selected_step"] =
                                selected_step;
                            phase_two_trace["objective_after"] =
                                selected_objective;
                            phase_two_trace["active_slack_after"] =
                                slack_sum(selected_state.p_delta);
                            phase_two_trace["reactive_slack_after"] =
                                slack_sum(selected_state.q_delta);
                            phase_two_trace["validation"] =
                                selected_validation.to_json();
                            output.economic_balance_polish_trace.push_back(
                                std::move(phase_two_trace));
                            if (!selected) {
                                break;
                            }
                            polished_state = std::move(selected_state);
                            polished_objective = selected_objective;
                            polished_validation = selected_validation;
                        }
                        output.economic_balance_polish_objective_after =
                            polished_objective;
                        output.economic_balance_polish_active_slack_after =
                            slack_sum(polished_state.p_delta);
                        output.economic_balance_polish_reactive_slack_after =
                            slack_sum(polished_state.q_delta);
                        output.economic_balance_polish_validation =
                            polished_validation;
                        output.economic_balance_polish_wall_seconds =
                            economic_polish_elapsed();
                        if (polished_objective > predictor_objective + 1e-9 &&
                            polished_validation.max_residual <=
                                options_.validation_tolerance) {
                            output.economic_balance_polish_selected = true;
                            predictor_state = std::move(polished_state);
                            predictor_objective = polished_objective;
                            predictor_validation = polished_validation;
                        }
                    }
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
                const auto correction_search_start =
                    std::chrono::steady_clock::now();
                if (predictor_iteration == 0 && outaged_branch >= 0 &&
                    !data_.branches[outaged_branch].transformer) {
                    branch_outage_low_rank_update =
                        predictor_cache_->configure_branch_outage_update(
                            data_, base_state_, outaged_branch);
                    if (options_.capture_diagnostics) {
                        output.fixed_jacobian_predictor_trace.back()[
                            "branch_outage_low_rank_prepared"] =
                            branch_outage_low_rank_update;
                    }
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
                if (options_.capture_diagnostics) {
                    output.fixed_jacobian_predictor_trace.back()[
                        "reactive_active_set_size"] = std::count_if(
                            q_balance_by_bus.begin(), q_balance_by_bus.end(),
                            [](double value) {
                                return std::abs(value) > 0.49 + 1e-8;
                            });
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
                int selected_coordinate_bus = -1;
                constexpr std::array<double, 5> kDampingCandidates{
                    1.0, 0.5, 0.25, 0.125, 0.0625};
                const auto strongly_improving_full_damping = [&]
                    (double damping, const ValidationReport& validation,
                     double sufficient_decrease_ratio = 0.99) {
                    // Use the damping loop as a backtracking line search, not
                    // an exhaustive best-of-five search.  Once the full step
                    // has achieved verified sufficient decrease, smaller
                    // steps in the same correction family need not be built
                    // and exhaustively validated.  This gate is independent
                    // of the branch low-rank path, so generator outages get
                    // the same shortcut.  Requiring the full step to halve
                    // the residual avoids changing basins when a smaller
                    // damping is only modestly better.  Callers can retain
                    // that conservative half-residual gate for correction
                    // families with known basin sensitivity.  A weaker full
                    // step still executes the complete deterministic search.
                    return damping == 1.0 &&
                        validation.max_residual <=
                            sufficient_decrease_ratio *
                                predictor_validation.max_residual;
                };
                const bool active_block_dominant =
                    predictor_validation.worst_category ==
                        "active_balance" &&
                    predictor_validation.max_active_balance_residual > 1e-4;
                int rejection_priority_balance_bus = -1;
                double rejection_priority_balance_residual = -1.0;
                if (predictor_validation.worst_category ==
                        "active_balance" ||
                    predictor_validation.worst_category ==
                        "reactive_balance") {
                    const bool active_priority =
                        predictor_validation.worst_category ==
                            "active_balance";
                    for (int bus = 0; bus < nb; ++bus) {
                        const double balance = active_priority
                            ? p_balance_by_bus[bus]
                            : q_balance_by_bus[bus];
                        const double slack_value = active_priority
                            ? predictor_state.p_delta[bus]
                            : predictor_state.q_delta[bus];
                        const double residual = std::max(
                            0.0, std::abs(balance) - slack_value);
                        if (residual >
                            rejection_priority_balance_residual) {
                            rejection_priority_balance_residual = residual;
                            rejection_priority_balance_bus = bus;
                        }
                    }
                }
                int worst_active_flow_branch = -1;
                bool worst_active_flow_from_side = true;
                double worst_active_flow_excess = 0.0;
                double worst_active_flow_target = 0.0;
                int worst_reactive_flow_branch = -1;
                bool worst_reactive_flow_from_side = true;
                double worst_reactive_flow_excess = 0.0;
                double worst_reactive_flow_target = 0.0;
                int rejection_priority_branch = -1;
                double rejection_priority_branch_residual = -1.0;
                const std::string sm_slack_suffix = ":sm_slack";
                const bool apparent_slack_dominant =
                    predictor_validation.worst_category ==
                        "variable_bound" &&
                    predictor_validation.worst_identity.size() >=
                        sm_slack_suffix.size() &&
                    predictor_validation.worst_identity.compare(
                        predictor_validation.worst_identity.size() -
                            sm_slack_suffix.size(),
                        sm_slack_suffix.size(), sm_slack_suffix) == 0;
                for (int branch_index = 0;
                     branch_index < static_cast<int>(data_.branches.size());
                     ++branch_index) {
                    if (branch_index == outaged_branch ||
                        data_.branches[branch_index].status == 0) {
                        continue;
                    }
                    const double rating =
                        data_.branches[branch_index].rate_c;
                    if (rating <= 1e-12) {
                        continue;
                    }
                    const auto& branch = data_.branches[branch_index];
                    const auto consider_terminal = [&]
                        (double active_flow, double reactive_flow,
                         double voltage, bool from_side) {
                        const double allowed_magnitude = rating *
                            (branch.transformer
                                 ? 1.0 + data_.sm_vio_limit
                                 : voltage + data_.sm_vio_limit);
                        const double active_limit = std::sqrt(std::max(
                            0.0,
                            allowed_magnitude * allowed_magnitude -
                                reactive_flow * reactive_flow));
                        const double reactive_limit = std::sqrt(std::max(
                            0.0,
                            allowed_magnitude * allowed_magnitude -
                                active_flow * active_flow));
                        const double active_excess =
                            std::abs(active_flow) - active_limit;
                        const double reactive_excess =
                            std::abs(reactive_flow) - reactive_limit;
                        if (active_excess > worst_active_flow_excess) {
                            worst_active_flow_excess = active_excess;
                            worst_active_flow_branch = branch_index;
                            worst_active_flow_from_side = from_side;
                            worst_active_flow_target = std::copysign(
                                std::max(0.0, active_limit - 1e-4),
                                active_flow);
                        }
                        if (reactive_excess >
                            worst_reactive_flow_excess) {
                            worst_reactive_flow_excess = reactive_excess;
                            worst_reactive_flow_branch = branch_index;
                            worst_reactive_flow_from_side = from_side;
                            worst_reactive_flow_target = std::copysign(
                                std::max(0.0, reactive_limit - 1e-4),
                                reactive_flow);
                        }
                    };
                    consider_terminal(
                        predictor_state.pf[branch_index],
                        predictor_state.qf[branch_index],
                        predictor_state.vm[branch.from], true);
                    consider_terminal(
                        predictor_state.pt[branch_index],
                        predictor_state.qt[branch_index],
                        predictor_state.vm[branch.to], false);
                    double priority_residual = -1.0;
                    if (predictor_validation.worst_category ==
                            "variable_bound") {
                        const auto bound_excess = [](
                            double value, double lower, double upper) {
                            return std::max(
                                std::max(0.0, lower - value),
                                std::max(0.0, value - upper));
                        };
                        priority_residual = bound_excess(
                            predictor_state.sm_slack[branch_index],
                            0.0, data_.sm_vio_limit);
                    } else if (predictor_validation.worst_category ==
                                   "flow_limit") {
                        const double from_squared =
                            predictor_state.pf[branch_index] *
                                predictor_state.pf[branch_index] +
                            predictor_state.qf[branch_index] *
                                predictor_state.qf[branch_index];
                        const double to_squared =
                            predictor_state.pt[branch_index] *
                                predictor_state.pt[branch_index] +
                            predictor_state.qt[branch_index] *
                                predictor_state.qt[branch_index];
                        const double from_limit = rating * rating *
                            std::pow(
                                branch.transformer
                                    ? 1.0 +
                                        predictor_state.sm_slack[branch_index]
                                    : predictor_state.vm[branch.from] +
                                        predictor_state.sm_slack[branch_index],
                                2);
                        const double to_limit = rating * rating *
                            std::pow(
                                branch.transformer
                                    ? 1.0 +
                                        predictor_state.sm_slack[branch_index]
                                    : predictor_state.vm[branch.to] +
                                        predictor_state.sm_slack[branch_index],
                                2);
                        priority_residual = std::max(
                            std::max(0.0, from_squared - from_limit),
                            std::max(0.0, to_squared - to_limit));
                    } else if (predictor_validation.worst_category ==
                                   "angle") {
                        const double start_delta =
                            base_state_.va[branch.from] -
                            base_state_.va[branch.to];
                        if (start_delta >= branch.angmin &&
                            start_delta <= branch.angmax) {
                            const double angle =
                                predictor_state.va[branch.from] -
                                predictor_state.va[branch.to];
                            priority_residual = std::max(
                                std::max(0.0, angle - branch.angmax),
                                std::max(0.0, branch.angmin - angle));
                        }
                    }
                    if (priority_residual >
                        rejection_priority_branch_residual) {
                        rejection_priority_branch_residual =
                            priority_residual;
                        rejection_priority_branch = branch_index;
                    }
                }
                const bool active_flow_blocks_balance_repair =
                    predictor_validation.worst_category ==
                        "active_balance" &&
                    predictor_validation.max_active_balance_residual >
                        1e-8 &&
                    worst_active_flow_excess >=
                        0.25 * predictor_validation
                            .max_active_balance_residual;
                const bool active_flow_bound_dominant =
                    (predictor_validation.worst_category ==
                         "variable_bound" ||
                     active_flow_blocks_balance_repair) &&
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
                const bool full_reactive_active_set_allowed =
                    outaged_branch >= 0 &&
                    !data_.branches[outaged_branch].transformer &&
                    predictor_validation.worst_category ==
                        "reactive_balance" &&
                    predictor_validation.max_variable_bound_violation <=
                        1e-5 &&
                    predictor_validation.max_flow_limit_violation <= 1e-5;
                if (options_.capture_diagnostics) {
                    output.fixed_jacobian_predictor_trace.back()[
                        "full_reactive_active_set_allowed"] =
                        full_reactive_active_set_allowed;
                    output.fixed_jacobian_predictor_trace.back()[
                        "worst_active_flow_excess"] =
                        worst_active_flow_excess;
                    output.fixed_jacobian_predictor_trace.back()[
                        "worst_reactive_flow_excess"] =
                        worst_reactive_flow_excess;
                }
                const auto project_trial_reactive_and_validate =
                    [&](AcState& trial) {
                    normalize_source_reference_angles(
                        data_, components, trial.va);
                    // Correction trials need fresh nonlinear branch flows,
                    // but their nodal slacks are only meaningful after the
                    // local active/reactive redispatch below.  Avoid the
                    // redundant pre-projection nodal-balance pass.
                    compute_branch_flows(
                        data_, outaged_branch, true, trial);
                    std::fill(
                        trial_p_network.begin(), trial_p_network.end(), 0.0);
                    std::fill(
                        trial_q_network.begin(), trial_q_network.end(), 0.0);
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
                            if (allocate_total(
                                    active_at_bus[bus], p_lower, p_upper,
                                    trial.pg, target_generation,
                                    trial.pg)) {
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
                        if (!allocate_total(
                                data_.buses[bus].loads,
                                predictor_load_power_lower,
                                predictor_load_power_upper,
                                trial_load_active, target_load,
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
                        static_cast<void>(allocate_total(
                            active_at_bus[bus], q_lower, q_upper,
                            initial_state.qg,
                            generated_q + q_balance -
                                std::clamp(q_balance, -0.49, 0.49),
                            trial.qg));
                    }
                    const auto trial_balance =
                        nodal_balance_slack_seed_from_network(
                            data_, trial, trial_p_network, trial_q_network,
                            0.5, 1e-7);
                    trial.p_delta = trial_balance.active;
                    trial.q_delta = trial_balance.reactive;
                    if (options_.capture_diagnostics) {
                        return validate_rebuilt_contingency_trial(
                            data_, trial, commitment_, *direct_context);
                    }
                    return
                        validate_rebuilt_contingency_trial_until_rejected(
                            data_, trial, commitment_, *direct_context,
                            selected_validation.max_residual,
                            rejection_priority_balance_bus,
                            rejection_priority_branch);
                };
                const auto try_damped_corrections = [&]
                    (FixedJacobianPredictorCache& cache) {
                    const auto try_local_reactive =
                        [&](double damping) {
                        auto trial = correction_reference;
                        if (!cache.apply_local_reactive_least_squares(
                                data_, q_balance_by_bus,
                                trial.vm, damping,
                                full_reactive_active_set_allowed)) {
                            return false;
                        }
                        const auto trial_validation =
                            project_trial_reactive_and_validate(trial);
                        if (trial_validation.max_residual + 1e-10 >=
                            selected_validation.max_residual) {
                            return false;
                        }
                        selected_correction = std::move(trial);
                        selected_validation = trial_validation;
                        selected_damping = damping;
                        selected_correction_mode =
                            "local_reactive_least_squares";
                        return true;
                    };
                    const auto try_local_reactive_then_active =
                        [&](double damping) {
                        auto trial = correction_reference;
                        if (!cache.apply_local_reactive_least_squares(
                                data_, q_balance_by_bus,
                                trial.vm, damping,
                                full_reactive_active_set_allowed)) {
                            return false;
                        }
                        compute_branch_flows(
                            data_, outaged_branch, true, trial);
                        std::fill(
                            trial_p_network.begin(),
                            trial_p_network.end(), 0.0);
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
                            return false;
                        }
                        const auto trial_validation =
                            project_trial_reactive_and_validate(trial);
                        if (trial_validation.max_residual + 1e-10 >=
                            selected_validation.max_residual) {
                            return false;
                        }
                        selected_correction = std::move(trial);
                        selected_validation = trial_validation;
                        selected_damping = damping;
                        selected_correction_mode =
                            "local_reactive_least_squares_then_"
                            "active_angle";
                        return true;
                    };
                    const auto try_reactive_band =
                        [&](double damping,
                            bool require_continuation_decrease) {
                        auto trial = correction_reference;
                        if (!cache.apply_reactive_band_correction(
                                data_, q_balance_by_bus,
                                trial.vm, damping)) {
                            return false;
                        }
                        const auto trial_validation =
                            project_trial_reactive_and_validate(trial);
                        if (require_continuation_decrease &&
                            trial_validation.max_residual >
                                0.999 *
                                    predictor_validation.max_residual) {
                            return false;
                        }
                        if (trial_validation.max_residual + 1e-10 >=
                            selected_validation.max_residual) {
                            return false;
                        }
                        selected_correction = std::move(trial);
                        selected_validation = trial_validation;
                        selected_damping = damping;
                        selected_correction_mode =
                            "reactive_feasibility_band";
                        return true;
                    };
                    const auto try_reactive_voltage_then_active = [&]() {
                        for (const double damping : kDampingCandidates) {
                            auto trial = correction_reference;
                            if (!cache.apply_reactive_correction(
                                    data_, q_spec, q_network,
                                    trial.vm, damping)) {
                                continue;
                            }
                            compute_branch_flows(
                                data_, outaged_branch, true, trial);
                            std::fill(
                                trial_p_network.begin(),
                                trial_p_network.end(), 0.0);
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
                    };

                    // Once this combined correction has selected a full step
                    // twice consecutively, that full step usually remains
                    // strongly improving at the next nonlinear state.  Try it
                    // as a continuation step
                    // before the ten-candidate generic search, but accept it
                    // early once it gives a material exact-residual reduction.
                    // Requiring five percent here made slowly contracting but
                    // strictly monotone reactive repairs repeat the complete
                    // ten-candidate search for more than a hundred iterations.
                    // A 0.1-percent gate still rejects numerical stagnation;
                    // final acceptance continues to require the complete
                    // independent validator below the configured tolerance.
                    if (!active_block_dominant &&
                        !active_flow_bound_dominant &&
                        prior_selected_correction_mode ==
                            "reactive_feasibility_band" &&
                        prior_selected_damping > 0.0 &&
                        try_reactive_band(
                            prior_selected_damping, true)) {
                        if (options_.capture_diagnostics) {
                            output.fixed_jacobian_predictor_trace.back()[
                                "continuation_correction_selected"] = true;
                        }
                        return;
                    }
                    if (!active_block_dominant &&
                        !active_flow_bound_dominant &&
                        prior_selected_correction_mode ==
                            "local_reactive_least_squares" &&
                        std::abs(prior_selected_damping - 1.0) <= 1e-12 &&
                        try_local_reactive(prior_selected_damping) &&
                        selected_validation.max_residual <=
                            0.999 * predictor_validation.max_residual) {
                        if (options_.capture_diagnostics) {
                            output.fixed_jacobian_predictor_trace.back()[
                                "continuation_correction_selected"] = true;
                        }
                        return;
                    }
                    if (!active_block_dominant &&
                        !active_flow_bound_dominant &&
                        consecutive_full_local_reactive_active_angle >= 1 &&
                        try_local_reactive_then_active(1.0) &&
                        selected_validation.max_residual <=
                            0.999 * predictor_validation.max_residual) {
                        if (options_.capture_diagnostics) {
                            output.fixed_jacobian_predictor_trace.back()[
                                "continuation_correction_selected"] = true;
                        }
                        return;
                    }
                    // When the exact validator identifies a reactive-flow
                    // bound as the dominant residual, the generic search
                    // commonly rejects several complete correction families
                    // before selecting this coupled voltage/angle step. Try
                    // that same fully rebuilt and validated family first.
                    // A material-decrease gate protects the nonlinear basin;
                    // a weak or failed priority attempt is discarded and the
                    // original deterministic search remains unchanged.
                    if (reactive_flow_bound_dominant) {
                        try_reactive_voltage_then_active();
                        if (selected_correction_mode ==
                                "reactive_voltage_then_active_angle" &&
                            selected_validation.max_residual <=
                                0.99 * predictor_validation.max_residual) {
                            if (options_.capture_diagnostics) {
                                output.fixed_jacobian_predictor_trace.back()[
                                    "reactive_flow_priority_correction_"
                                    "selected"] = true;
                            }
                            return;
                        }
                        selected_correction = correction_reference;
                        selected_validation = predictor_validation;
                        selected_damping = 0.0;
                        selected_correction_mode = "none";
                    }
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
                        const double target_flow =
                            worst_active_flow_target;
                        const bool prefer_redispatch_over_endpoint_angle =
                            worst_reactive_flow_excess <= 1e-6 ||
                            (predictor_iteration >= 5 &&
                             predictor_validation.max_residual > 0.02);
                        const auto try_endpoint_angle = [&]() {
                            if (std::abs(derivative) <= 1e-10) {
                                return;
                            }
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
                        };
                        if (!prefer_redispatch_over_endpoint_angle) {
                            try_endpoint_angle();
                        }
                        auto full_redispatch_pg =
                            correction_reference.pg;
                        auto full_redispatch_p_spec = p_spec;
                        AcState best_redispatch = correction_reference;
                        ValidationReport best_redispatch_validation =
                            predictor_validation;
                        double best_redispatch_damping = 0.0;
                        bool best_redispatch_distributed = false;
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
                                    if (damping == 1.0 &&
                                        prior_selected_correction_mode ==
                                            "active_branch_flow_redispatch" &&
                                        std::abs(
                                            prior_selected_damping - 1.0) <=
                                            1e-12 &&
                                        trial_validation.max_residual <=
                                            0.9999 *
                                                predictor_validation
                                                    .max_residual) {
                                        if (options_.capture_diagnostics) {
                                            output
                                                .fixed_jacobian_predictor_trace
                                                .back()[
                                                    "redispatch_continuation_"
                                                    "selected"] = true;
                                        }
                                        break;
                                    }
                                }
                            }
                        }
                        if (best_redispatch_damping == 0.0) {
                            auto distributed_pg = correction_reference.pg;
                            const int distributed_load_count =
                                static_cast<int>(data_.loads.size());
                            std::vector<double> distributed_load_power(
                                static_cast<std::size_t>(
                                    distributed_load_count), 0.0);
                            for (int load = 0;
                                 load < distributed_load_count; ++load) {
                                distributed_load_power[load] =
                                    data_.loads[load].pd_nominal *
                                    correction_reference
                                        .demand_factor[load];
                            }
                            const auto reference_load_power =
                                distributed_load_power;
                            auto distributed_p_spec = p_spec;
                            if (cache.apply_distributed_active_flow_redispatch(
                                    data_, worst_active_flow_branch,
                                    worst_active_flow_from_side,
                                    target_flow - flow,
                                    p_lower, p_upper,
                                    predictor_load_power_lower,
                                    predictor_load_power_upper,
                                    correction_reference,
                                    distributed_pg,
                                    distributed_load_power,
                                    distributed_p_spec)) {
                                for (const double damping :
                                     kDampingCandidates) {
                                    auto trial = correction_reference;
                                    auto trial_p_spec = p_spec;
                                    for (int generator = 0;
                                         generator < ng; ++generator) {
                                        trial.pg[generator] += damping *
                                            (distributed_pg[generator] -
                                             correction_reference
                                                 .pg[generator]);
                                    }
                                    for (int bus = 0; bus < nb; ++bus) {
                                        trial_p_spec[bus] += damping *
                                            (distributed_p_spec[bus] -
                                             p_spec[bus]);
                                    }
                                    for (int load = 0;
                                         load < distributed_load_count;
                                         ++load) {
                                        if (std::abs(
                                                data_.loads[load]
                                                    .pd_nominal) <= 1e-12) {
                                            continue;
                                        }
                                        const double trial_load_power =
                                            reference_load_power[load] +
                                            damping *
                                                (distributed_load_power[load] -
                                                 reference_load_power[load]);
                                        trial.demand_factor[load] =
                                            trial_load_power /
                                            data_.loads[load].pd_nominal;
                                    }
                                    if (!cache.apply_active_correction(
                                            data_, trial_p_spec, p_network,
                                            trial.va, 1.0)) {
                                        continue;
                                    }
                                    const auto trial_validation =
                                        project_trial_reactive_and_validate(
                                            trial);
                                    if (trial_validation.max_residual +
                                            1e-10 <
                                        best_redispatch_validation
                                            .max_residual) {
                                        best_redispatch = std::move(trial);
                                        best_redispatch_validation =
                                            trial_validation;
                                        best_redispatch_damping = damping;
                                        best_redispatch_distributed = true;
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
                                best_redispatch_distributed
                                ? "active_branch_flow_distributed_redispatch"
                                : "active_branch_flow_redispatch";
                        } else if (prefer_redispatch_over_endpoint_angle) {
                            // Endpoint-angle candidates are irrelevant when
                            // an improving redispatch exists.  Retain them as
                            // the exact fallback when redispatch cannot move
                            // this flow constraint.
                            try_endpoint_angle();
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
                                    damping, trial_validation, 0.5)) {
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
                        if (try_local_reactive(damping)) {
                            if (strongly_improving_full_damping(
                                    damping, selected_validation)) {
                                break;
                            }
                        }
                    }
                    for (const double damping : kDampingCandidates) {
                        if (active_block_dominant) {
                            break;
                        }
                        if (try_local_reactive_then_active(damping) &&
                            strongly_improving_full_damping(
                                damping, selected_validation)) {
                            break;
                        }
                    }
                    if (selected_damping != 0.0) {
                        return;
                    }
                    for (const double damping : kDampingCandidates) {
                        if (active_block_dominant) {
                            break;
                        }
                        if (try_reactive_band(damping, false)) {
                            if (strongly_improving_full_damping(
                                    damping, selected_validation)) {
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
                    if (!active_block_dominant) {
                        try_reactive_voltage_then_active();
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
                const bool prior_one_hop_coordinate =
                    prior_selected_correction_mode ==
                        "one_hop_voltage_coordinate" ||
                    prior_selected_correction_mode ==
                        "one_hop_voltage_coordinate_continuation";
                if (!coordinate_history.empty() &&
                    predictor_validation.worst_category ==
                        "reactive_balance") {
                    auto history_trace = nlohmann::json::array();
                    for (const auto& [candidate_bus, voltage_change] :
                         coordinate_history) {
                        if (candidate_bus < 0 || candidate_bus >= nb ||
                            std::abs(voltage_change) < 1e-12) {
                            continue;
                        }
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
                        if (options_.capture_diagnostics) {
                            history_trace.push_back({
                                {"candidate_bus",
                                 data_.buses[candidate_bus].bus_i},
                                {"voltage_change", voltage_change},
                                {"validation",
                                 trial_validation.to_json()},
                            });
                        }
                        if (trial_validation.max_residual + 1e-10 <
                            selected_validation.max_residual) {
                            selected_correction = std::move(trial);
                            selected_validation = trial_validation;
                            selected_damping = voltage_change;
                            selected_correction_mode =
                                "one_hop_voltage_coordinate_continuation";
                            selected_coordinate_bus = candidate_bus;
                        }
                    }
                    if (options_.capture_diagnostics) {
                        output.fixed_jacobian_predictor_trace.back()[
                            "coordinate_history_trials"] =
                            std::move(history_trace);
                    }
                    // The complete neighborhood search accepts any exact
                    // residual decrease. Reuse a recent coordinate only when
                    // that same complete validation shows a measurable
                    // decrease; otherwise discard it and retain the original
                    // search below.
                    if (selected_validation.max_residual >
                        0.9999 * predictor_validation.max_residual) {
                        selected_correction = correction_reference;
                        selected_validation = predictor_validation;
                        selected_damping = 0.0;
                        selected_correction_mode = "none";
                        selected_coordinate_bus = -1;
                    }
                }

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
                    predictor_validation.max_residual >
                        0.25 * initial_predictor_validation_residual;
                if (selected_damping == 0.0 &&
                    !contingency_predictor_cache && slow_branch_outage) {
                    contingency_predictor_cache =
                        std::make_unique<FixedJacobianPredictorCache>(
                            data_, correction_reference, commitment_,
                            outaged_branch);
                    output.fixed_jacobian_predictor_preparation_seconds +=
                        contingency_predictor_cache->preparation_seconds;
                    if (options_.capture_diagnostics) {
                        output.fixed_jacobian_predictor_trace.back()[
                            "contingency_specific_refactorization"] = true;
                        output.fixed_jacobian_predictor_trace.back()[
                            "contingency_specific_refactorization_reason"] =
                            "slow_branch_outage_progress";
                    }
                }
                if (selected_damping == 0.0) {
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
                            if (options_.capture_diagnostics) {
                                output.fixed_jacobian_predictor_trace.back()[
                                    "contingency_specific_refactorization"] =
                                    true;
                            }
                            if (contingency_predictor_cache->valid) {
                                try_damped_corrections(
                                    *contingency_predictor_cache);
                            }
                        }
                    }
                }
                if (selected_damping == 0.0 && prior_one_hop_coordinate &&
                    std::abs(prior_selected_damping) >= 1e-12) {
                    int continuation_worst_q_bus = -1;
                    double continuation_worst_q_excess = 0.0;
                    for (int bus = 0; bus < nb; ++bus) {
                        const double excess =
                            std::abs(q_balance_by_bus[bus]) - 0.49;
                        if (excess > continuation_worst_q_excess) {
                            continuation_worst_q_excess = excess;
                            continuation_worst_q_bus = bus;
                        }
                    }
                    std::vector<int> continuation_candidate_buses;
                    if (continuation_worst_q_bus >= 0) {
                        continuation_candidate_buses.push_back(
                            continuation_worst_q_bus);
                        for (int branch_index :
                             data_.buses[continuation_worst_q_bus]
                                 .branches_from) {
                            if (branch_index != outaged_branch &&
                                data_.branches[branch_index].status != 0) {
                                continuation_candidate_buses.push_back(
                                    data_.branches[branch_index].to);
                            }
                        }
                        for (int branch_index :
                             data_.buses[continuation_worst_q_bus]
                                 .branches_to) {
                            if (branch_index != outaged_branch &&
                                data_.branches[branch_index].status != 0) {
                                continuation_candidate_buses.push_back(
                                    data_.branches[branch_index].from);
                            }
                        }
                        std::sort(
                            continuation_candidate_buses.begin(),
                            continuation_candidate_buses.end());
                        continuation_candidate_buses.erase(
                            std::unique(
                                continuation_candidate_buses.begin(),
                                continuation_candidate_buses.end()),
                            continuation_candidate_buses.end());
                    }
                    auto continuation_trials = nlohmann::json::array();
                    AcState continuation_best = correction_reference;
                    ValidationReport continuation_best_validation =
                        predictor_validation;
                    double continuation_best_damping = 0.0;
                    int continuation_best_bus = -1;
                    const std::array<double, 2>
                        continuation_voltage_changes{
                            prior_selected_damping,
                            std::copysign(0.005, prior_selected_damping),
                    };
                    for (int candidate_bus :
                         continuation_candidate_buses) {
                        for (std::size_t change_index = 0;
                             change_index <
                                 continuation_voltage_changes.size();
                             ++change_index) {
                            const double voltage_change =
                                continuation_voltage_changes[change_index];
                            if (change_index > 0 &&
                                std::abs(
                                    voltage_change -
                                    continuation_voltage_changes[0]) <=
                                1e-12) {
                                continue;
                            }
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
                            if (options_.capture_diagnostics) {
                                continuation_trials.push_back({
                                    {"candidate_bus",
                                     data_.buses[candidate_bus].bus_i},
                                    {"voltage_change", voltage_change},
                                    {"validation", trial_validation.to_json()},
                                });
                            }
                            if (trial_validation.max_residual + 1e-10 <
                                continuation_best_validation.max_residual) {
                                continuation_best = std::move(trial);
                                continuation_best_validation =
                                    trial_validation;
                                continuation_best_damping = voltage_change;
                                continuation_best_bus = candidate_bus;
                            }
                        }
                    }
                    if (options_.capture_diagnostics) {
                        output.fixed_jacobian_predictor_trace.back()[
                            "coordinate_continuation_trials"] =
                            std::move(continuation_trials);
                    }
                    if (continuation_best_validation.max_residual <=
                        0.99 * predictor_validation.max_residual) {
                        selected_correction =
                            std::move(continuation_best);
                        selected_validation =
                            continuation_best_validation;
                        selected_damping = continuation_best_damping;
                        selected_correction_mode =
                            "one_hop_voltage_coordinate_continuation";
                        selected_coordinate_bus = continuation_best_bus;
                    }
                }
                // Once an apparent-power slack bound is only improving by
                // microscopic amounts, accepting that nominal step can
                // suppress the localized voltage search for dozens of
                // iterations. Permit one proactive local search at the first
                // verified stall. The existing candidate remains available,
                // and every localized trial is still rebuilt and validated
                // against the complete nonlinear contingency model.
                const bool proactive_apparent_slack_local_search =
                    selected_damping != 0.0 &&
                    proactive_local_search_attempts == 0 &&
                    active_feasibility_repair_attempts > 0 &&
                    predictor_iteration >= 8 &&
                    apparent_slack_dominant &&
                    predictor_validation.max_residual >
                        options_.validation_tolerance &&
                    selected_validation.max_residual >=
                        0.99 * predictor_validation.max_residual;
                if (selected_damping == 0.0 ||
                    proactive_apparent_slack_local_search) {
                    if (proactive_apparent_slack_local_search) {
                        ++proactive_local_search_attempts;
                        if (options_.capture_diagnostics) {
                            output.fixed_jacobian_predictor_trace.back()[
                                "proactive_apparent_slack_local_search"] =
                                true;
                        }
                    }
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
                        if (options_.capture_diagnostics) {
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
                    }
                    auto* shunt_correction_cache =
                        contingency_predictor_cache &&
                            contingency_predictor_cache->valid
                        ? contingency_predictor_cache.get()
                        : predictor_cache_.get();
                    auto shunt_trial_trace = nlohmann::json::array();
                    if (options_.capture_diagnostics) {
                        output.fixed_jacobian_predictor_trace.back()[
                            "worst_q_bus"] = worst_q_bus >= 0
                            ? data_.buses[worst_q_bus].bus_i : -1;
                        output.fixed_jacobian_predictor_trace.back()[
                            "worst_q_balance"] = worst_q_bus >= 0
                            ? q_balance_by_bus[worst_q_bus] : 0.0;
                    }
                    auto reactive_flow_voltage_trace =
                        nlohmann::json::array();
                    if (reactive_flow_bound_dominant) {
                        const auto& reactive_branch =
                            data_.branches[worst_reactive_flow_branch];
                        const int from = reactive_branch.from;
                        const int to = reactive_branch.to;
                        const double impedance_denominator =
                            reactive_branch.r * reactive_branch.r +
                            reactive_branch.x * reactive_branch.x;
                        const double g = impedance_denominator > 1e-20
                            ? reactive_branch.r / impedance_denominator
                            : 0.0;
                        const double b = impedance_denominator > 1e-20
                            ? -reactive_branch.x / impedance_denominator
                            : 0.0;
                        const double tm2 =
                            reactive_branch.tap * reactive_branch.tap;
                        const double tr = reactive_branch.tap *
                            std::cos(reactive_branch.shift);
                        const double ti = reactive_branch.tap *
                            std::sin(reactive_branch.shift);
                        const double vf = correction_reference.vm[from];
                        const double vt = correction_reference.vm[to];
                        const double angle =
                            correction_reference.va[from] -
                            correction_reference.va[to];
                        const double cosine = std::cos(angle);
                        const double sine = std::sin(angle);
                        const double flow = worst_reactive_flow_from_side
                            ? correction_reference
                                  .qf[worst_reactive_flow_branch]
                            : correction_reference
                                  .qt[worst_reactive_flow_branch];
                        double derivative_from = 0.0;
                        double derivative_to = 0.0;
                        if (worst_reactive_flow_from_side) {
                            const double from_b_self =
                                reactive_branch.transformer
                                ? b / tm2 + reactive_branch.b_fr
                                : (b + reactive_branch.b_fr) / tm2;
                            const double cosine_coefficient =
                                (-b * tr - g * ti) / tm2;
                            const double sine_coefficient =
                                (-g * tr + b * ti) / tm2;
                            derivative_from =
                                -2.0 * from_b_self * vf -
                                cosine_coefficient * vt * cosine +
                                sine_coefficient * vt * sine;
                            derivative_to =
                                -cosine_coefficient * vf * cosine +
                                sine_coefficient * vf * sine;
                        } else {
                            const double cosine_coefficient =
                                (-b * tr + g * ti) / tm2;
                            const double sine_coefficient =
                                (-g * tr - b * ti) / tm2;
                            derivative_from =
                                -cosine_coefficient * vt * cosine -
                                sine_coefficient * vt * sine;
                            derivative_to =
                                -2.0 * (b + reactive_branch.b_to) * vt -
                                cosine_coefficient * vf * cosine -
                                sine_coefficient * vf * sine;
                        }
                        const double gradient_square =
                            derivative_from * derivative_from +
                            derivative_to * derivative_to;
                        if (shunt_correction_cache != nullptr &&
                            shunt_correction_cache->reactive_valid &&
                            worst_q_bus >= 0) {
                            const int balance_row =
                                shunt_correction_cache
                                    ->voltage_index[worst_q_bus];
                            const int from_column =
                                shunt_correction_cache
                                    ->voltage_index[from];
                            const int to_column =
                                shunt_correction_cache
                                    ->voltage_index[to];
                            const double balance_from_derivative =
                                shunt_correction_cache->reactive_jacobian
                                    .coeff(balance_row, from_column);
                            const double balance_to_derivative =
                                shunt_correction_cache->reactive_jacobian
                                    .coeff(balance_row, to_column);
                            const double two_by_two_determinant =
                                derivative_from * balance_to_derivative -
                                derivative_to * balance_from_derivative;
                            if (std::abs(two_by_two_determinant) > 1e-12) {
                                const double target_flow =
                                    worst_reactive_flow_target;
                                const double desired_flow_change =
                                    target_flow - flow;
                                const double desired_balance_change =
                                    std::clamp(
                                        q_balance_by_bus[worst_q_bus],
                                        -0.48, 0.48) -
                                    q_balance_by_bus[worst_q_bus];
                                const double raw_from_change = std::clamp(
                                    (desired_flow_change *
                                         balance_to_derivative -
                                     derivative_to *
                                         desired_balance_change) /
                                        two_by_two_determinant,
                                    -0.02, 0.02);
                                const double raw_to_change = std::clamp(
                                    (derivative_from *
                                         desired_balance_change -
                                     desired_flow_change *
                                         balance_from_derivative) /
                                        two_by_two_determinant,
                                    -0.02, 0.02);
                                for (double damping : kDampingCandidates) {
                                    const double proposed_from = std::clamp(
                                        vf + damping * raw_from_change,
                                        data_.buses[from].vmin,
                                        data_.buses[from].vmax);
                                    const double proposed_to = std::clamp(
                                        vt + damping * raw_to_change,
                                        data_.buses[to].vmin,
                                        data_.buses[to].vmax);
                                    if (std::abs(proposed_from - vf) <=
                                            1e-12 &&
                                        std::abs(proposed_to - vt) <=
                                            1e-12) {
                                        continue;
                                    }
                                    auto trial = correction_reference;
                                    trial.vm[from] = proposed_from;
                                    trial.vm[to] = proposed_to;
                                    const auto trial_validation =
                                        project_trial_reactive_and_validate(
                                            trial);
                                    if (options_.capture_diagnostics) {
                                        reactive_flow_voltage_trace.push_back({
                                            {"targeted_flow_balance_step", true},
                                            {"balance_bus",
                                             data_.buses[worst_q_bus].bus_i},
                                            {"from_bus",
                                             data_.buses[from].bus_i},
                                            {"to_bus", data_.buses[to].bus_i},
                                            {"from_voltage_change",
                                             proposed_from - vf},
                                            {"to_voltage_change",
                                             proposed_to - vt},
                                            {"damping", damping},
                                            {"validation",
                                             trial_validation.to_json()},
                                        });
                                    }
                                    if (trial_validation.max_residual +
                                            1e-10 <
                                        selected_validation.max_residual) {
                                        selected_correction =
                                            std::move(trial);
                                        selected_validation =
                                            trial_validation;
                                        selected_damping = damping;
                                        selected_correction_mode =
                                            "reactive_branch_flow_balance";
                                    }
                                }
                            }
                        }
                        if (gradient_square > 1e-16) {
                            const double target_flow = std::copysign(
                                std::max(
                                    0.0, reactive_branch.rate_c - 1e-4),
                                flow);
                            const double desired_change = target_flow - flow;
                            const double raw_from_change = std::clamp(
                                desired_change * derivative_from /
                                    gradient_square,
                                -0.02, 0.02);
                            const double raw_to_change = std::clamp(
                                desired_change * derivative_to /
                                    gradient_square,
                                -0.02, 0.02);
                            for (double damping : kDampingCandidates) {
                                const double proposed_from = std::clamp(
                                    vf + damping * raw_from_change,
                                    data_.buses[from].vmin,
                                    data_.buses[from].vmax);
                                const double proposed_to = std::clamp(
                                    vt + damping * raw_to_change,
                                    data_.buses[to].vmin,
                                    data_.buses[to].vmax);
                                if (std::abs(proposed_from - vf) <= 1e-12 &&
                                    std::abs(proposed_to - vt) <= 1e-12) {
                                    continue;
                                }
                                auto trial = correction_reference;
                                trial.vm[from] = proposed_from;
                                trial.vm[to] = proposed_to;
                                const auto trial_validation =
                                    project_trial_reactive_and_validate(
                                        trial);
                                if (options_.capture_diagnostics) {
                                    reactive_flow_voltage_trace.push_back({
                                        {"targeted_gradient_step", true},
                                        {"from_bus", data_.buses[from].bus_i},
                                        {"to_bus", data_.buses[to].bus_i},
                                        {"from_voltage_change",
                                         proposed_from - vf},
                                        {"to_voltage_change", proposed_to - vt},
                                        {"damping", damping},
                                        {"validation",
                                         trial_validation.to_json()},
                                    });
                                }
                                const double
                                    maximum_preserved_reactive_residual =
                                        std::max(
                                            1e-3,
                                            1.25 * predictor_validation
                                                .max_reactive_balance_residual +
                                                1e-5);
                                if (trial_validation
                                        .max_reactive_balance_residual >
                                    maximum_preserved_reactive_residual) {
                                    continue;
                                }
                                if (trial_validation.max_residual + 1e-10 <
                                    selected_validation.max_residual) {
                                    selected_correction = std::move(trial);
                                    selected_validation = trial_validation;
                                    selected_damping = damping;
                                    selected_correction_mode =
                                        "reactive_branch_flow_gradient";
                                }
                            }
                        }
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
                                if (options_.capture_diagnostics) {
                                    reactive_flow_voltage_trace.push_back({
                                        {"candidate_bus",
                                         data_.buses[candidate_bus].bus_i},
                                        {"voltage_change", voltage_change},
                                        {"validation",
                                         trial_validation.to_json()},
                                    });
                                }
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
                    if (options_.capture_diagnostics) {
                        output.fixed_jacobian_predictor_trace.back()[
                            "reactive_flow_voltage_trials"] =
                            std::move(reactive_flow_voltage_trace);
                    }
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
                            if (options_.capture_diagnostics) {
                                local_reactive_trace.push_back({
                                    {"damping", damping},
                                    {"direction", direction},
                                    {"desired_network_q_change",
                                     direction *
                                         desired_network_q_change},
                                    {"validation", trial_validation.to_json()},
                                });
                            }
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
                    if (options_.capture_diagnostics) {
                        output.fixed_jacobian_predictor_trace.back()[
                            "local_reactive_trials"] =
                            std::move(local_reactive_trace);
                    }
                    auto local_voltage_trace = nlohmann::json::array();
                    constexpr std::array<double, 4>
                        kLocalVoltageChanges{
                            -0.005, -0.001, 0.001, 0.005};
                    bool strong_local_voltage_coordinate_selected = false;
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
                            if (options_.capture_diagnostics) {
                                local_voltage_trace.push_back({
                                    {"candidate_bus",
                                     data_.buses[candidate_bus].bus_i},
                                    {"voltage_change", voltage_change},
                                    {"validation", trial_validation.to_json()},
                                });
                            }
                            if (trial_validation.max_residual + 1e-10 <
                                selected_validation.max_residual) {
                                selected_correction = std::move(trial);
                                selected_validation = trial_validation;
                                selected_damping = voltage_change;
                                selected_correction_mode =
                                    "one_hop_voltage_coordinate";
                                selected_coordinate_bus = candidate_bus;
                                // Treat this coordinate loop as a guarded
                                // line search.  A fully rebuilt candidate
                                // that removes at least ten percent of the
                                // exact dominant residual is already a
                                // useful nonlinear step; evaluating every
                                // remaining neighbor and voltage magnitude
                                // only changes which improving point is
                                // chosen.  We retain the exhaustive search
                                // whenever no candidate clears this
                                // conservative gate.
                                if (trial_validation.max_residual <=
                                    0.90 * predictor_validation.max_residual) {
                                    strong_local_voltage_coordinate_selected =
                                        true;
                                    break;
                                }
                            }
                        }
                        if (strong_local_voltage_coordinate_selected) {
                            break;
                        }
                    }
                    if (options_.capture_diagnostics) {
                        output.fixed_jacobian_predictor_trace.back()[
                            "local_voltage_trials"] =
                            std::move(local_voltage_trace);
                        output.fixed_jacobian_predictor_trace.back()[
                            "strong_local_voltage_coordinate_selected"] =
                            strong_local_voltage_coordinate_selected;
                    }
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
                                nlohmann::json shunt_trial;
                                if (options_.capture_diagnostics) {
                                    shunt_trial = {
                                        {"candidate_bus",
                                         data_.buses[candidate_bus].bus_i},
                                        {"shunt", shunt.source_key},
                                        {"block", block},
                                        {"current_step", current_step},
                                        {"proposed_step", proposed_step},
                                        {"delta_bs", delta_bs},
                                    };
                                }
                                // q_balance includes -bs*|V|^2, so positive
                                // excess is relieved by increasing bs and
                                // negative excess by decreasing it.
                                if (worst_q_bus < 0 ||
                                    (candidate_bus == worst_q_bus &&
                                     q_balance_by_bus[worst_q_bus] *
                                         delta_bs <= 0.0)) {
                                    if (options_.capture_diagnostics) {
                                        shunt_trial["skipped"] =
                                            "opposite_reactive_direction";
                                        shunt_trial_trace.push_back(
                                            std::move(shunt_trial));
                                    }
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
                                if (options_.capture_diagnostics) {
                                    shunt_trial["direct_validation"] =
                                        trial_validation.to_json();
                                }
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
                                std::fill(
                                    trial_q_network.begin(),
                                    trial_q_network.end(), 0.0);
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
                                if (options_.capture_diagnostics) {
                                    shunt_trial["reactive_corrections"] =
                                        nlohmann::json::array();
                                }
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
                                    if (options_.capture_diagnostics) {
                                        shunt_trial["reactive_corrections"]
                                            .push_back({
                                                {"damping", damping},
                                                {"validation",
                                                 corrected_validation.to_json()},
                                            });
                                    }
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
                                if (options_.capture_diagnostics) {
                                    shunt_trial_trace.push_back(
                                        std::move(shunt_trial));
                                }
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
                                        if (options_.capture_diagnostics) {
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
                                        }
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
                    if (options_.capture_diagnostics) {
                        output.fixed_jacobian_predictor_trace.back()[
                            "paired_voltage_trials"] =
                            std::move(paired_voltage_trace);
                        output.fixed_jacobian_predictor_trace.back()[
                            "shunt_trials"] = std::move(shunt_trial_trace);
                    }
                }
                if (options_.capture_diagnostics) {
                    output.fixed_jacobian_predictor_trace.back()[
                        "selected_damping"] = selected_damping;
                    output.fixed_jacobian_predictor_trace.back()[
                        "selected_correction_mode"] =
                        selected_correction_mode;
                    output.fixed_jacobian_predictor_trace.back()[
                        "selected_coordinate_bus"] =
                        selected_coordinate_bus >= 0
                        ? data_.buses[selected_coordinate_bus].bus_i : -1;
                    output.fixed_jacobian_predictor_trace.back()[
                        "selected_next_validation"] =
                        selected_validation.to_json();
                }
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
                const bool branch_security_stagnation_repair =
                    outaged_branch >= 0 &&
                    active_feasibility_repair_attempts == 0 &&
                    predictor_iteration >= 12 &&
                    predictor_validation.worst_category ==
                        "variable_bound" &&
                    worst_active_flow_excess > 1e-4 &&
                    selected_validation.max_residual >=
                        0.99 * predictor_validation.max_residual;
                // Generator outages can enter the same active-flow-bound
                // stall while still finding microscopic redispatch steps.
                // Requiring selected_damping == 0 meant those tiny steps
                // suppressed Phase I indefinitely; difficult 19k cases then
                // consumed more than one hundred predictor iterations before
                // the repair was first attempted.  Trigger the identical
                // source-bounded LP when verified progress has fallen below
                // one percent and an active terminal bound is still blocking
                // the generator-outage correction.
                const bool generator_security_stagnation_repair =
                    outaged_generator >= 0 &&
                    active_feasibility_repair_attempts == 0 &&
                    predictor_iteration >= 12 &&
                    predictor_validation.worst_category ==
                        "variable_bound" &&
                    worst_active_flow_excess > 1e-4 &&
                    selected_validation.max_residual >=
                        0.99 * predictor_validation.max_residual;
                if (initial_active_feasibility_repair ||
                    branch_security_stagnation_repair ||
                    generator_security_stagnation_repair ||
                    late_security_stagnation_repair) {
                    ++active_feasibility_repair_attempts;
                    // Branch outages can also stall at an active-flow box
                    // bound while every nodal balance is already inside its
                    // source-authorized band.  The same source-bounded LP is
                    // valid for either outage type.  Give its first branch
                    // invocation the broader trust region; every returned
                    // candidate is still rebuilt through the nonlinear AC
                    // equations and independently validated before use.
                    const bool broad_active_repair =
                        initial_active_feasibility_repair ||
                        branch_security_stagnation_repair ||
                        generator_security_stagnation_repair;
                    const double active_balance_limit =
                        outaged_branch >= 0
                            ? 0.49
                            : (broad_active_repair ? 0.25 : 0.20);
                    const double active_angle_trust =
                        broad_active_repair ? 0.05 : 0.02;
                    const double active_voltage_trust =
                        broad_active_repair ? 0.005 : 0.002;
                    const double active_repair_time_limit =
                        5.0;
                    auto active_repair =
                        solve_linearized_active_feasibility_repair(
                            data_, correction_reference, commitment_,
                            *direct_context, active_balance_limit,
                            active_angle_trust, active_repair_time_limit,
                            active_voltage_trust,
                            outaged_branch >= 0,
                            outaged_branch >= 0);
                    nlohmann::json active_repair_json;
                    if (options_.capture_diagnostics) {
                        active_repair_json = active_repair.to_json(false);
                    }
                    if (active_repair.success) {
                        auto active_trial = std::move(active_repair.state);
                        const auto active_validation =
                            project_trial_reactive_and_validate(active_trial);
                        if (options_.capture_diagnostics) {
                            active_repair_json["nonlinear_validation"] =
                                active_validation.to_json();
                        }
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
                    if (options_.capture_diagnostics) {
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
                }
                // The inexpensive active-only Phase I above is sufficient for
                // most generator outages, but it cannot move voltage or
                // reactive generation.  A small number of large-case outages
                // therefore reach a genuine final stall at a reactive-balance
                // or apparent-flow bound even though a nearby corrective AC
                // point exists.  Before escalating to the much more expensive
                // nonlinear fallback, solve one current-violation-only
                // linearized AC Phase I.  This retains every source ramp,
                // generator, load, voltage, angle, and security bound.  The
                // returned point is never trusted directly: nonlinear branch
                // flows are rebuilt and the independent contingency validator
                // must show strict improvement before the predictor can use it.
                const bool final_generator_ac_feasibility_repair =
                    outaged_generator >= 0 &&
                    ac_feasibility_repair_attempts == 0 &&
                    selected_damping == 0.0 &&
                    predictor_validation.max_residual >
                        options_.validation_tolerance;
                if (final_generator_ac_feasibility_repair) {
                    ++ac_feasibility_repair_attempts;
                    auto ac_repair =
                        solve_linearized_active_feasibility_repair(
                            data_, correction_reference, commitment_,
                            *direct_context, 0.49, 0.05, 5.0, 0.01,
                            true, true);
                    nlohmann::json ac_repair_json;
                    if (options_.capture_diagnostics) {
                        ac_repair_json = ac_repair.to_json(false);
                    }
                    if (ac_repair.success) {
                        auto ac_trial = std::move(ac_repair.state);
                        const auto ac_validation =
                            project_trial_reactive_and_validate(ac_trial);
                        if (options_.capture_diagnostics) {
                            ac_repair_json["nonlinear_validation"] =
                                ac_validation.to_json();
                        }
                        if (ac_validation.max_residual + 1e-10 <
                            selected_validation.max_residual) {
                            selected_correction = std::move(ac_trial);
                            selected_validation = ac_validation;
                            selected_damping = 1.0;
                            selected_correction_mode =
                                "linearized_ac_feasibility_repair";
                        }
                    }
                    if (options_.capture_diagnostics) {
                        output.fixed_jacobian_predictor_trace.back()[
                            "ac_feasibility_repair"] =
                            std::move(ac_repair_json);
                        output.fixed_jacobian_predictor_trace.back()[
                            "selected_damping"] = selected_damping;
                        output.fixed_jacobian_predictor_trace.back()[
                            "selected_correction_mode"] =
                            selected_correction_mode;
                        output.fixed_jacobian_predictor_trace.back()[
                            "selected_next_validation"] =
                            selected_validation.to_json();
                    }
                }
                if (options_.capture_diagnostics) {
                    output.fixed_jacobian_predictor_trace.back()[
                        "correction_search_seconds"] =
                        std::chrono::duration<double>(
                            std::chrono::steady_clock::now() -
                            correction_search_start).count();
                }
                if (selected_damping == 0.0) {
                    break;
                }
                if (selected_correction_mode ==
                        "local_reactive_least_squares_then_active_angle" &&
                    selected_damping == 1.0) {
                    ++consecutive_full_local_reactive_active_angle;
                } else {
                    consecutive_full_local_reactive_active_angle = 0;
                }
                prior_selected_correction_mode = selected_correction_mode;
                prior_selected_damping = selected_damping;
                if ((selected_correction_mode ==
                         "one_hop_voltage_coordinate" ||
                     selected_correction_mode ==
                         "one_hop_voltage_coordinate_continuation") &&
                    selected_coordinate_bus >= 0) {
                    coordinate_history.erase(
                        std::remove_if(
                            coordinate_history.begin(),
                            coordinate_history.end(),
                            [&](const auto& item) {
                                return item.first == selected_coordinate_bus &&
                                    std::abs(
                                        item.second - selected_damping) <=
                                        1e-12;
                            }),
                        coordinate_history.end());
                    coordinate_history.insert(
                        coordinate_history.begin(),
                        {selected_coordinate_bus, selected_damping});
                    constexpr std::size_t kCoordinateHistoryLimit = 4;
                    if (coordinate_history.size() >
                        kCoordinateHistoryLimit) {
                        coordinate_history.resize(
                            kCoordinateHistoryLimit);
                    }
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
        normalize_source_reference_angles(
            data_, components, output.solve.state.va);
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
        active_slack_target[bus] =
            base_mode && options_.minimize_active_balance_slack
            ? std::clamp(
                (1.0 - options_.balance_cleanup_fraction) * p_balance,
                -0.49, 0.49)
            : std::clamp(p_balance, -0.49, 0.49);
        reactive_slack_target[bus] =
            base_mode && options_.minimize_reactive_balance_slack
            ? std::clamp(
                (1.0 - options_.balance_cleanup_fraction) * q_balance,
                -0.49, 0.49)
            : std::clamp(q_balance, -0.49, 0.49);
    }
    for (int bus = 0; bus < nb; ++bus) {
        for (int gen : active_at_bus[bus]) {
            fixed_q_bus[bus] += qg[gen];
        }
    }

    const auto evaluate_newton_candidate = [&] (
        bool clamp_voltage,
        bool reallocate_generation) {
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
        if (reallocate_generation) {
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
                const double p_target =
                    base_mode && options_.minimize_active_balance_slack
                    ? active_slack_target[bus]
                    : std::clamp(p_balance, -0.49, 0.49);
                const double q_target =
                    base_mode && options_.minimize_reactive_balance_slack
                    ? reactive_slack_target[bus]
                    : std::clamp(q_balance, -0.49, 0.49);
                static_cast<void>(allocate_total(
                    active_at_bus[bus], p_lower, p_upper,
                    initial_state.pg,
                    current_pg + p_balance - p_target,
                    candidate.pg));
                static_cast<void>(allocate_total(
                    active_at_bus[bus], q_lower, q_upper,
                    initial_state.qg,
                    current_qg + q_balance - q_target,
                    candidate.qg));
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

    if (base_mode && options_.minimize_active_balance_slack &&
        !options_.skip_balance_cleanup_prepasses) {
        output.local_balance_candidate_attempted = true;
        AcState local_candidate = direct_state;
        std::vector<double> load_power(data_.loads.size(), 0.0);
        std::vector<double> load_power_lower(data_.loads.size(), 0.0);
        std::vector<double> load_power_upper(data_.loads.size(), 0.0);
        for (int load_index = 0;
             load_index < static_cast<int>(data_.loads.size());
             ++load_index) {
            const auto& load = data_.loads[load_index];
            load_power[load_index] =
                load.pd_nominal * local_candidate.demand_factor[load_index];
            if (std::abs(load.pd_nominal) <= 1e-12) {
                load_power_lower[load_index] = load_power[load_index];
                load_power_upper[load_index] = load_power[load_index];
                continue;
            }
            const double factor_lower = std::max(
                load.tmin,
                (load.pd_prev - data_.delta_r * load.prdmax) /
                    load.pd_nominal);
            const double factor_upper = std::min(
                load.tmax,
                (load.pd_prev + data_.delta_r * load.prumax) /
                    load.pd_nominal);
            load_power_lower[load_index] = std::min(
                load.pd_nominal * factor_lower,
                load.pd_nominal * factor_upper);
            load_power_upper[load_index] = std::max(
                load.pd_nominal * factor_lower,
                load.pd_nominal * factor_upper);
        }
        const auto signed_active_balance = [&](int bus) {
            double balance = 0.0;
            for (int branch : data_.buses[bus].branches_from) {
                balance += direct_state.pf[branch];
            }
            for (int branch : data_.buses[bus].branches_to) {
                balance += direct_state.pt[branch];
            }
            for (int gen : data_.buses[bus].generators) {
                balance -= local_candidate.pg[gen];
            }
            for (int load_index : data_.buses[bus].loads) {
                balance += load_power[load_index];
            }
            for (int shunt : data_.buses[bus].shunts) {
                const double vm2 = direct_state.vm[bus] *
                    direct_state.vm[bus];
                balance += data_.shunts[shunt].gs * vm2;
            }
            return balance;
        };
        for (int bus = 0; bus < nb; ++bus) {
            double balance = signed_active_balance(bus);
            if (std::abs(balance) <= 1e-10) {
                continue;
            }
            if (!active_at_bus[bus].empty()) {
                double current = 0.0;
                double total_lower = 0.0;
                double total_upper = 0.0;
                for (int gen : active_at_bus[bus]) {
                    current += local_candidate.pg[gen];
                    total_lower += p_lower[gen];
                    total_upper += p_upper[gen];
                }
                const double target = std::clamp(
                    current + balance, total_lower, total_upper);
                if (allocate_total(
                        active_at_bus[bus], p_lower, p_upper,
                        local_candidate.pg, target,
                        local_candidate.pg)) {
                    balance -= target - current;
                }
            }
            if (std::abs(balance) <= 1e-10 ||
                data_.buses[bus].loads.empty()) {
                continue;
            }
            double current = 0.0;
            double total_lower = 0.0;
            double total_upper = 0.0;
            for (int load_index : data_.buses[bus].loads) {
                current += load_power[load_index];
                total_lower += load_power_lower[load_index];
                total_upper += load_power_upper[load_index];
            }
            const double target = std::clamp(
                current - balance, total_lower, total_upper);
            if (allocate_total(
                    data_.buses[bus].loads,
                    load_power_lower, load_power_upper,
                    load_power, target, load_power)) {
                for (int load_index : data_.buses[bus].loads) {
                    const double nominal =
                        data_.loads[load_index].pd_nominal;
                    if (std::abs(nominal) > 1e-12) {
                        local_candidate.demand_factor[load_index] =
                            load_power[load_index] / nominal;
                    }
                }
            }
        }
        const auto signed_reactive_balance = [&](
            const AcState& candidate,
            int bus) {
            double balance = 0.0;
            for (int branch : data_.buses[bus].branches_from) {
                balance += direct_state.qf[branch];
            }
            for (int branch : data_.buses[bus].branches_to) {
                balance += direct_state.qt[branch];
            }
            for (int gen : data_.buses[bus].generators) {
                balance -= candidate.qg[gen];
            }
            for (int load_index : data_.buses[bus].loads) {
                balance += data_.loads[load_index].qd_nominal *
                    candidate.demand_factor[load_index];
            }
            for (int shunt : data_.buses[bus].shunts) {
                const double vm2 = direct_state.vm[bus] *
                    direct_state.vm[bus];
                balance -= effective_shunt_susceptance(
                    data_, candidate, shunt) * vm2;
            }
            return balance;
        };
        constexpr double kReactiveSlackInterior = 0.499999;
        for (int bus = 0; bus < nb; ++bus) {
            if (!active_at_bus[bus].empty() ||
                data_.buses[bus].loads.empty()) {
                continue;
            }
            const double original_q =
                signed_reactive_balance(direct_state, bus);
            const double proposed_q =
                signed_reactive_balance(local_candidate, bus);
            double step = 1.0;
            if (proposed_q > kReactiveSlackInterior &&
                proposed_q > original_q + 1e-14) {
                step = (kReactiveSlackInterior - original_q) /
                    (proposed_q - original_q);
            } else if (proposed_q < -kReactiveSlackInterior &&
                       proposed_q < original_q - 1e-14) {
                step = (-kReactiveSlackInterior - original_q) /
                    (proposed_q - original_q);
            }
            step = std::clamp(step, 0.0, 1.0);
            if (step >= 1.0) {
                continue;
            }
            for (int load_index : data_.buses[bus].loads) {
                local_candidate.demand_factor[load_index] =
                    direct_state.demand_factor[load_index] + step *
                    (local_candidate.demand_factor[load_index] -
                     direct_state.demand_factor[load_index]);
            }
        }
        const auto project_local_qg = [&](int bus) {
            if (active_at_bus[bus].empty()) {
                return;
            }
            const double q_balance =
                signed_reactive_balance(local_candidate, bus);
            double current = 0.0;
            double total_lower = 0.0;
            double total_upper = 0.0;
            for (int gen : active_at_bus[bus]) {
                current += local_candidate.qg[gen];
                total_lower += q_lower[gen];
                total_upper += q_upper[gen];
            }
            const double target = std::clamp(
                current + q_balance, total_lower, total_upper);
            static_cast<void>(allocate_total(
                active_at_bus[bus], q_lower, q_upper,
                local_candidate.qg, target, local_candidate.qg));
        };
        for (int bus = 0; bus < nb; ++bus) {
            project_local_qg(bus);
        }
        ensure_shunt_control_state(data_, local_candidate);
        for (int bus = 0; bus < nb; ++bus) {
            if (data_.buses[bus].shunts.empty()) {
                continue;
            }
            double q_balance =
                signed_reactive_balance(local_candidate, bus);
            for (int coordinate_pass = 0;
                 coordinate_pass < 512;
                 ++coordinate_pass) {
                const double current_absolute = std::abs(q_balance);
                double best_absolute = current_absolute;
                double best_balance = q_balance;
                int best_shunt = -1;
                int best_block = -1;
                int best_step_change = 0;
                for (int shunt_index : data_.buses[bus].shunts) {
                    const auto& shunt = data_.shunts[shunt_index];
                    if (!shunt.dispatchable ||
                        shunt_index >= static_cast<int>(
                            local_candidate.shunt_steps.size())) {
                        continue;
                    }
                    for (int block = 0;
                         block < static_cast<int>(
                             shunt.block_maximum_steps.size());
                         ++block) {
                        if (block >= static_cast<int>(
                                local_candidate
                                    .shunt_steps[shunt_index].size())) {
                            continue;
                        }
                        const int current_step =
                            local_candidate.shunt_steps[shunt_index][block];
                        for (const int step_change : {-1, 1}) {
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
                            const double proposed_balance = q_balance -
                                delta_bs * direct_state.vm[bus] *
                                    direct_state.vm[bus];
                            const double proposed_absolute =
                                std::abs(proposed_balance);
                            if (proposed_absolute + 1e-12 < best_absolute) {
                                best_absolute = proposed_absolute;
                                best_balance = proposed_balance;
                                best_shunt = shunt_index;
                                best_block = block;
                                best_step_change = step_change;
                            }
                        }
                    }
                }
                if (best_shunt < 0) {
                    break;
                }
                local_candidate.shunt_steps[best_shunt][best_block] +=
                    best_step_change;
                local_candidate.shunt_bs[best_shunt] +=
                    static_cast<double>(best_step_change) *
                    data_.shunts[best_shunt]
                        .block_susceptance[best_block];
                q_balance = best_balance;
            }
            project_local_qg(bus);
        }
        const AcState initial_local_state = direct_state;
        AcState direct_objective_state = direct_state;
        const double direct_objective = rebuild_base_state_derived_fields(
            data_, commitment_, direct_objective_state, 0.5);
        double local_objective = rebuild_base_state_derived_fields(
            data_, commitment_, local_candidate, 0.5);
        auto local_validation = validate_state(
            data_, ModelMode::BaseSoft,
            local_candidate, commitment_);
        output.local_balance_candidate_validation = local_validation;
        double selected_step = 1.0;
        if (local_validation.max_residual >
                options_.validation_tolerance ||
            local_objective <= direct_objective + 1e-9) {
            const AcState full_local_candidate = local_candidate;
            selected_step = 0.0;
            for (const double step :
                 std::array<double, 12>{
                     0.75, 0.5, 0.25, 0.125, 0.0625, 0.03125,
                     0.015625, 0.0078125, 0.00390625,
                     0.001953125, 0.0009765625, 0.00048828125}) {
                ++output.local_balance_backtracking_attempts;
                AcState trial = initial_local_state;
                for (int gen = 0; gen < ng; ++gen) {
                    trial.pg[gen] += step *
                        (full_local_candidate.pg[gen] -
                         initial_local_state.pg[gen]);
                    trial.qg[gen] += step *
                        (full_local_candidate.qg[gen] -
                         initial_local_state.qg[gen]);
                }
                for (int load_index = 0;
                     load_index < static_cast<int>(data_.loads.size());
                     ++load_index) {
                    trial.demand_factor[load_index] += step *
                        (full_local_candidate.demand_factor[load_index] -
                         initial_local_state.demand_factor[load_index]);
                }
                const double trial_objective =
                    rebuild_base_state_derived_fields(
                        data_, commitment_, trial, 0.5);
                const auto trial_validation = validate_state(
                    data_, ModelMode::BaseSoft, trial, commitment_);
                output.local_balance_candidate_validation = trial_validation;
                if (trial_validation.max_residual <=
                        options_.validation_tolerance &&
                    trial_objective > direct_objective + 1e-9) {
                    local_candidate = std::move(trial);
                    local_validation = trial_validation;
                    local_objective = trial_objective;
                    selected_step = step;
                    break;
                }
            }
        }
        if (selected_step > 0.0) {
            output.converged = true;
            output.feasible = true;
            output.local_balance_candidate_selected = true;
            output.local_balance_selected_step = selected_step;
            output.solve.status = 0;
            output.solve.objective = local_objective;
            output.solve.iterations = 0;
            output.solve.state = std::move(local_candidate);
            output.validation = local_validation;
            output.wall_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - wall_start).count();
            output.solve.wall_seconds = output.wall_seconds;
            return output;
        }
    }

    const double direct_active_slack_sum = std::accumulate(
        direct_state.p_delta.begin(), direct_state.p_delta.end(), 0.0);
    const double direct_reactive_slack_sum = std::accumulate(
        direct_state.q_delta.begin(), direct_state.q_delta.end(), 0.0);
    if (base_mode && options_.minimize_reactive_balance_slack &&
        !options_.skip_balance_cleanup_prepasses &&
        options_.balance_cleanup_fraction >= 1.0 - 1e-12 &&
        direct_active_slack_sum <= 0.05 &&
        direct_reactive_slack_sum > 1e-8) {
        output.reactive_only_newton_attempted = true;
        const auto reactive_initial_vm = vm;
        const auto reactive_initial_va = va;
        std::vector<bool> fixed_angles(nb, true);
        std::vector<double> unused_p_spec(nb, 0.0);
        std::vector<double> zero_balance_q_spec(nb, 0.0);
        for (int bus = 0; bus < nb; ++bus) {
            zero_balance_q_spec[bus] = -load_q[bus];
            for (int gen : active_at_bus[bus]) {
                zero_balance_q_spec[bus] += qg[gen];
            }
        }
        const auto reactive_only = run_newton(
            data_, ybus, fixed_angles, pq,
            unused_p_spec, zero_balance_q_spec,
            options_.max_newton_iterations,
            options_.newton_tolerance, vm, va);
        output.reactive_only_newton_converged = reactive_only.converged;
        output.reactive_only_newton_iterations = reactive_only.iterations;
        output.newton_iterations += reactive_only.iterations;

        AcState direct_objective_state = direct_state;
        const double direct_objective = rebuild_base_state_derived_fields(
            data_, commitment_, direct_objective_state, 0.5);
        bool selected = false;
        double best_objective = direct_objective;
        double best_step = 0.0;
        AcState best_state = direct_state;
        ValidationReport best_validation =
            output.direct_candidate_validation;
        const std::array<double, 13> steps{
            1.0, 0.75, 0.5, 0.25, 0.125, 0.0625, 0.03125,
            0.015625, 0.0078125, 0.00390625,
            0.001953125, 0.0009765625, 0.00048828125};
        for (const double step : steps) {
            ++output.reactive_only_backtracking_attempts;
            AcState trial = direct_state;
            for (int bus = 0; bus < nb; ++bus) {
                trial.vm[bus] = std::clamp(
                    reactive_initial_vm[bus] + step *
                        (vm[bus] - reactive_initial_vm[bus]),
                    data_.buses[bus].vmin,
                    data_.buses[bus].vmax);
                trial.va[bus] = reactive_initial_va[bus];
            }
            compute_branch_flows(
                data_, outaged_branch, false, trial);
            for (int bus = 0; bus < nb; ++bus) {
                if (active_at_bus[bus].empty()) {
                    continue;
                }
                double q_balance = 0.0;
                for (int branch : data_.buses[bus].branches_from) {
                    q_balance += trial.qf[branch];
                }
                for (int branch : data_.buses[bus].branches_to) {
                    q_balance += trial.qt[branch];
                }
                double current_qg = 0.0;
                for (int gen : data_.buses[bus].generators) {
                    q_balance -= trial.qg[gen];
                }
                for (int gen : active_at_bus[bus]) {
                    current_qg += trial.qg[gen];
                }
                for (int load_index : data_.buses[bus].loads) {
                    q_balance += data_.loads[load_index].qd_nominal *
                        trial.demand_factor[load_index];
                }
                for (int shunt : data_.buses[bus].shunts) {
                    const double vm2 = trial.vm[bus] * trial.vm[bus];
                    q_balance -= effective_shunt_susceptance(
                        data_, trial, shunt) * vm2;
                }
                double total_lower = 0.0;
                double total_upper = 0.0;
                for (int gen : active_at_bus[bus]) {
                    total_lower += q_lower[gen];
                    total_upper += q_upper[gen];
                }
                const double target = std::clamp(
                    current_qg + q_balance,
                    total_lower, total_upper);
                static_cast<void>(allocate_total(
                    active_at_bus[bus], q_lower, q_upper,
                    trial.qg, target, trial.qg));
            }
            double trial_objective = rebuild_base_state_derived_fields(
                data_, commitment_, trial, 0.5);
            auto trial_validation = validate_state(
                data_, ModelMode::BaseSoft, trial, commitment_);
            if (step == 1.0) {
                output.reactive_only_newton_validation =
                    trial_validation;
            }
            nlohmann::json trace = {
                {"step", step},
                {"objective_before_active_repair", trial_objective},
                {"validation_before_active_repair",
                 trial_validation.to_json()},
                {"active_repair_attempted", false},
            };
            FastPowerFlowOptions active_repair_options;
            active_repair_options.minimize_active_balance_slack = true;
            active_repair_options.minimize_reactive_balance_slack = false;
            active_repair_options.max_newton_iterations = 20;
            active_repair_options.max_active_redispatch_passes = 8;
            FastContingencyPowerFlow active_repair(
                data_, trial, commitment_, active_repair_options);
            auto repaired = active_repair.solve_base();
            trace["active_repair_attempted"] = true;
            trace["active_repair_feasible"] = repaired.feasible;
            trace["active_repair_wall_seconds"] = repaired.wall_seconds;
            trace["active_repair_validation"] =
                repaired.validation.to_json();
            if (repaired.feasible) {
                repaired.solve.objective = rebuild_base_state_derived_fields(
                    data_, commitment_, repaired.solve.state, 0.5);
                repaired.validation = validate_state(
                    data_, ModelMode::BaseSoft,
                    repaired.solve.state, commitment_);
                if (repaired.validation.max_residual <=
                        options_.validation_tolerance &&
                    repaired.solve.objective > trial_objective + 1e-9) {
                    trial = std::move(repaired.solve.state);
                    trial_objective = repaired.solve.objective;
                    trial_validation = repaired.validation;
                }
            }
            trace["objective_after_active_repair"] = trial_objective;
            trace["validation_after_active_repair"] =
                trial_validation.to_json();
            output.reactive_only_trace.push_back(std::move(trace));
            if (trial_validation.max_residual <=
                    options_.validation_tolerance &&
                trial_objective > best_objective + 1e-9) {
                selected = true;
                best_objective = trial_objective;
                best_step = step;
                best_state = std::move(trial);
                best_validation = trial_validation;
            }
        }
        vm = reactive_initial_vm;
        va = reactive_initial_va;
        if (selected) {
            output.converged = reactive_only.converged;
            output.feasible = true;
            output.reactive_only_newton_selected = true;
            output.reactive_only_selected_step = best_step;
            output.solve.status = reactive_only.converged ? 0 : 1;
            output.solve.objective = best_objective;
            output.solve.iterations = output.newton_iterations;
            output.solve.state = std::move(best_state);
            output.validation = best_validation;
            output.wall_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - wall_start).count();
            output.solve.wall_seconds = output.wall_seconds;
            return output;
        }
    }

    if (base_mode && options_.minimize_active_balance_slack &&
        !options_.skip_balance_cleanup_prepasses &&
        !(options_.minimize_reactive_balance_slack &&
          direct_active_slack_sum <= 0.05)) {
        const auto cleanup_initial_vm = vm;
        const auto cleanup_initial_va = va;
        const auto cleanup_initial_pg = pg;
        std::vector<bool> active_only_pq(nb, false);
        std::vector<double> active_only_p_spec(nb, 0.0);
        std::vector<double> unused_q_spec(nb, 0.0);
        for (int bus = 0; bus < nb; ++bus) {
            active_only_p_spec[bus] = -load_p[bus];
            for (int gen : active_at_bus[bus]) {
                active_only_p_spec[bus] += pg[gen];
            }
        }
        std::vector<double> initial_p_network;
        std::vector<double> initial_q_network;
        network_injections(
            ybus, vm, va, initial_p_network, initial_q_network);
        std::vector<double> active_slack_weights(nb, 0.0);
        bool active_slack_weights_valid = true;
        for (const auto& component : components) {
            double component_mismatch = 0.0;
            for (int bus : component) {
                component_mismatch +=
                    initial_p_network[bus] - active_only_p_spec[bus];
            }
            double total_room = 0.0;
            for (int bus : component) {
                for (int gen : active_at_bus[bus]) {
                    total_room += component_mismatch >= 0.0
                        ? std::max(0.0, p_upper[gen] - pg[gen])
                        : std::max(0.0, pg[gen] - p_lower[gen]);
                }
            }
            if (total_room <= 1e-12 ||
                total_room + kAllocationTolerance <
                    std::abs(component_mismatch)) {
                active_slack_weights_valid = false;
                break;
            }
            for (int bus : component) {
                double bus_room = 0.0;
                for (int gen : active_at_bus[bus]) {
                    bus_room += component_mismatch >= 0.0
                        ? std::max(0.0, p_upper[gen] - pg[gen])
                        : std::max(0.0, pg[gen] - p_lower[gen]);
                }
                active_slack_weights[bus] = bus_room / total_room;
            }
        }
        output.active_only_newton_attempted = true;
        const auto active_only = active_slack_weights_valid
            ? run_distributed_active_newton(
                data_, ybus, slack, active_only_pq, component_of,
                active_slack_weights,
                active_only_p_spec, unused_q_spec,
                options_.max_newton_iterations, options_.newton_tolerance,
                vm, va)
            : NewtonResult{
                false, 0,
                "insufficient active range for zero-balance cleanup"};
        output.active_only_newton_converged = active_only.converged;
        output.active_only_newton_iterations = active_only.iterations;
        output.newton_iterations += active_only.iterations;
        bool generation_allocation_valid = active_only.converged;
        if (generation_allocation_valid) {
            std::vector<double> solved_p_network;
            std::vector<double> solved_q_network;
            network_injections(
                ybus, vm, va, solved_p_network, solved_q_network);
            for (int bus = 0; bus < nb; ++bus) {
                if (active_at_bus[bus].empty()) {
                    continue;
                }
                double current_pg = 0.0;
                for (int gen : active_at_bus[bus]) {
                    current_pg += pg[gen];
                }
                const double bus_adjustment =
                    solved_p_network[bus] - active_only_p_spec[bus];
                if (!allocate_total(
                        active_at_bus[bus], p_lower, p_upper,
                        cleanup_initial_pg,
                        current_pg + bus_adjustment,
                        pg)) {
                    generation_allocation_valid = false;
                    break;
                }
            }
        }
        auto [active_only_state, active_only_validation] =
            evaluate_newton_candidate(false, false);
        if (!generation_allocation_valid) {
            active_only_validation.max_residual =
                std::numeric_limits<double>::infinity();
            active_only_validation.worst_category =
                "active_generation_allocation";
            active_only_validation.worst_identity = "base";
        }
        output.active_only_newton_validation = active_only_validation;
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
        if (generation_allocation_valid) {
            const auto zero_balance_vm = vm;
            const auto zero_balance_va = va;
            const auto zero_balance_pg = pg;
            for (double step = 0.5;
                 step >= 1.0 / 4096.0;
                 step *= 0.5) {
                ++output.active_only_backtracking_attempts;
                vm = cleanup_initial_vm;
                va = cleanup_initial_va;
                pg = cleanup_initial_pg;
                for (int bus = 0; bus < nb; ++bus) {
                    vm[bus] += step *
                        (zero_balance_vm[bus] - cleanup_initial_vm[bus]);
                    va[bus] += step *
                        (zero_balance_va[bus] - cleanup_initial_va[bus]);
                }
                for (int gen = 0; gen < ng; ++gen) {
                    pg[gen] += step *
                        (zero_balance_pg[gen] - cleanup_initial_pg[gen]);
                }
                auto [trial_state, trial_validation] =
                    evaluate_newton_candidate(false, false);
                output.active_only_backtracking_validation = trial_validation;
                if (trial_validation.max_residual >
                    options_.validation_tolerance) {
                    continue;
                }
                output.converged = active_only.converged;
                output.feasible = true;
                output.active_only_newton_selected = true;
                output.active_only_selected_step = step;
                output.solve.status = active_only.converged ? 0 : 1;
                output.solve.objective = 0.0;
                output.solve.iterations = output.newton_iterations;
                output.solve.state = std::move(trial_state);
                output.validation = trial_validation;
                output.wall_seconds = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - wall_start).count();
                output.solve.wall_seconds = output.wall_seconds;
                return output;
            }
        }
        vm = cleanup_initial_vm;
        va = cleanup_initial_va;
        pg = cleanup_initial_pg;
    }

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
            evaluate_newton_candidate(false, true);
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
            evaluate_newton_candidate(false, true);
        if (reactive_validation.worst_category == "variable_bound") {
            auto [clamped_state, clamped_validation] =
                evaluate_newton_candidate(true, true);
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
            const double thermal_ratio =
                maximum_thermal_ratio(reactive_state);
            if (thermal_ratio <= 1.0 + 1e-12) {
                break;
            }
            const double scale = std::min(
                0.9995, 0.9999 / std::sqrt(thermal_ratio));
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
                evaluate_newton_candidate(false, true);
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
    const auto requested_balance_cleanup_complete =
        [&](const AcState& candidate) {
            // Rebuilt states deliberately add a 1e-7 interior margin to
            // every nodal slack variable.  Treat exactly that bookkeeping
            // floor as zero; otherwise a large case can never satisfy a
            // nominal zero-slack cleanup request.
            const double active_cleanup_floor =
                1e-7 * static_cast<double>(candidate.p_delta.size()) +
                1e-8;
            const double reactive_cleanup_floor =
                1e-7 * static_cast<double>(candidate.q_delta.size()) +
                1e-8;
            const double active_slack = std::accumulate(
                candidate.p_delta.begin(),
                candidate.p_delta.end(), 0.0);
            const double reactive_slack = std::accumulate(
                candidate.q_delta.begin(),
                candidate.q_delta.end(), 0.0);
            return (!base_mode ||
                    !options_.minimize_active_balance_slack ||
                    active_slack <= active_cleanup_floor) &&
                (!base_mode ||
                 !options_.minimize_reactive_balance_slack ||
                 reactive_slack <= reactive_cleanup_floor);
        };
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
                evaluate_newton_candidate(false, true);
            if (newton_validation.worst_category == "variable_bound") {
                auto [clamped_state, clamped_validation] =
                    evaluate_newton_candidate(true, true);
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
                    options_.validation_tolerance &&
                requested_balance_cleanup_complete(newton_state)) {
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

    normalize_source_reference_angles(data_, components, va);
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
        context->borrow_base_state(base_state_);
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
        const double thermal_ratio =
            maximum_thermal_ratio(output.solve.state);
        const double scale = std::min(0.9995,
            0.9999 / std::sqrt(thermal_ratio));
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
    const bool large_base_has_better_validated_incumbent =
        base_mode && nb >= 16000 &&
        !options_.minimize_active_balance_slack &&
        !options_.minimize_reactive_balance_slack &&
        std::isfinite(best_intermediate_validation.max_residual) &&
        best_intermediate_validation.max_residual + 1e-12 <
            output.validation.max_residual;
    if (options_.distributed_balance_polish &&
        !large_base_has_better_validated_incumbent &&
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

            normalize_source_reference_angles(
                data_, components, state.va);
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

            const double thermal_ratio = maximum_thermal_ratio(state);
            const double scale = std::min(
                0.9995, 0.9999 / std::sqrt(thermal_ratio));
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
    if (!options_.minimize_active_balance_slack &&
        !options_.minimize_reactive_balance_slack &&
        best_intermediate_validation.max_residual + 1e-12 <
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
