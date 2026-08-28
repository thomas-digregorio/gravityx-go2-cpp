#include <highs/Highs.h>

#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include "gravityx/algorithm.hpp"
#include "gravityx/fast_power_flow.hpp"
#include "gravityx/state_io.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
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

struct BranchFlowDerivatives {
    std::vector<double> pf_angle;
    std::vector<double> qf_angle;
    std::vector<double> pt_angle;
    std::vector<double> qt_angle;
    std::vector<double> pf_vf;
    std::vector<double> pf_vt;
    std::vector<double> qf_vf;
    std::vector<double> qf_vt;
    std::vector<double> pt_vf;
    std::vector<double> pt_vt;
    std::vector<double> qt_vf;
    std::vector<double> qt_vt;
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

BranchFlowDerivatives branch_flow_derivatives(
    const CaseData& data,
    const AcState& reference) {
    const std::size_t count = data.branches.size();
    BranchFlowDerivatives output{
        std::vector<double>(count, 0.0),
        std::vector<double>(count, 0.0),
        std::vector<double>(count, 0.0),
        std::vector<double>(count, 0.0),
        std::vector<double>(count, 0.0),
        std::vector<double>(count, 0.0),
        std::vector<double>(count, 0.0),
        std::vector<double>(count, 0.0),
        std::vector<double>(count, 0.0),
        std::vector<double>(count, 0.0),
        std::vector<double>(count, 0.0),
        std::vector<double>(count, 0.0),
    };
    for (int index = 0;
         index < static_cast<int>(data.branches.size()); ++index) {
        const auto& branch = data.branches[index];
        if (!branch.present || branch.status == 0) {
            continue;
        }
        const auto coefficient = flow_coefficients(branch);
        const double vf = reference.vm[branch.from];
        const double vt = reference.vm[branch.to];
        const double angle =
            reference.va[branch.from] - reference.va[branch.to];
        const double sine = std::sin(angle);
        const double cosine = std::cos(angle);
        const double voltage_product = vf * vt;
        output.pf_angle[index] = voltage_product *
            (-coefficient.from_cross_cos * sine +
             coefficient.from_cross_sin * cosine);
        output.qf_angle[index] = voltage_product *
            (coefficient.from_cross_sin * sine +
             coefficient.from_cross_cos * cosine);
        output.pt_angle[index] = voltage_product *
            (-coefficient.to_cross_cos * sine -
             coefficient.to_cross_sin * cosine);
        output.qt_angle[index] = voltage_product *
            (coefficient.to_cross_sin * sine -
             coefficient.to_cross_cos * cosine);

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
        output.pf_vf[index] =
            2.0 * coefficient.from_g_self * vf + pf_cross * vt;
        output.pf_vt[index] = pf_cross * vf;
        output.qf_vf[index] =
            -2.0 * coefficient.from_b_self * vf + qf_cross * vt;
        output.qf_vt[index] = qf_cross * vf;
        output.pt_vf[index] = pt_cross * vt;
        output.pt_vt[index] =
            2.0 * coefficient.to_g_self * vt + pt_cross * vf;
        output.qt_vf[index] = qt_cross * vt;
        output.qt_vt[index] =
            -2.0 * coefficient.to_b_self * vt + qt_cross * vf;
    }
    return output;
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

// Eliminate the linearized AC state equations with one artificial active-power
// slack per connected island. Reactive balance and voltage-magnitude response
// are retained at every bus, while bounded reactive-generation movements are
// explicit reduced-LP controls. The artificial slacks are constrained to zero
// in the reduced LP, so this is an exact algebraic reformulation of the
// selected linearized equations rather than a relaxation. Adjoint solves map
// generator/load movements directly into the small set of monitored network
// rows without materializing a dense sensitivity matrix.
class ReducedActiveNetwork {
public:
    ReducedActiveNetwork(
        const CaseData& data,
        const std::vector<int>& commitment,
        const AcState& reference,
        const std::vector<int>& component_reference,
        const BranchFlowDerivatives& derivative,
        bool pv_pq_partition)
        : nb_(static_cast<int>(data.buses.size())),
          ng_(static_cast<int>(data.generators.size())),
          nd_(static_cast<int>(data.loads.size())),
          component_count_(static_cast<int>(component_reference.size())),
          pv_pq_partition_(pv_pq_partition),
          angle_column_(static_cast<std::size_t>(nb_), -1),
          voltage_column_(static_cast<std::size_t>(nb_), -1),
          reactive_row_(static_cast<std::size_t>(nb_), -1),
          pv_bus_(static_cast<std::size_t>(nb_), 0),
          control_entries_(static_cast<std::size_t>(
              ng_ + nd_ + (pv_pq_partition ? 0 : ng_))) {
        if (nb_ <= 0 || component_count_ <= 0 ||
            commitment.size() != data.generators.size() ||
            reference.qg.size() != data.generators.size() ||
            reference.vm.size() != data.buses.size() ||
            reference.va.size() != data.buses.size() ||
            derivative.pf_angle.size() != data.branches.size() ||
            derivative.qf_angle.size() != data.branches.size() ||
            derivative.pt_angle.size() != data.branches.size() ||
            derivative.qt_angle.size() != data.branches.size()) {
            failure_reason_ = "invalid reduced active-network dimensions";
            return;
        }
        std::vector<unsigned char> is_reference(
            static_cast<std::size_t>(nb_), 0);
        for (int reference : component_reference) {
            if (reference < 0 || reference >= nb_ ||
                is_reference[reference] != 0) {
                failure_reason_ =
                    "invalid reduced active-network component reference";
                return;
            }
            is_reference[reference] = 1;
        }
        for (int bus = 0; bus < nb_; ++bus) {
            if (is_reference[bus] == 0) {
                angle_column_[bus] = angle_count_++;
            }
        }
        if (pv_pq_partition_) {
            constexpr double kReactiveRoomTolerance = 1e-6;
            for (int bus = 0; bus < nb_; ++bus) {
                bool has_committed_generator = false;
                double reactive_lower = 0.0;
                double reactive_upper = 0.0;
                double reactive_value = 0.0;
                for (int generator : data.buses[bus].generators) {
                    if (commitment[generator] == 0) {
                        continue;
                    }
                    has_committed_generator = true;
                    reactive_lower += data.generators[generator].qmin;
                    reactive_upper += data.generators[generator].qmax;
                    reactive_value += reference.qg[generator];
                }
                if (has_committed_generator &&
                    reactive_value >
                        reactive_lower + kReactiveRoomTolerance &&
                    reactive_value <
                        reactive_upper - kReactiveRoomTolerance) {
                    pv_bus_[bus] = 1;
                }
            }
        }
        for (int bus = 0; bus < nb_; ++bus) {
            if (pv_bus_[bus] != 0) {
                continue;
            }
            voltage_column_[bus] = angle_count_ + voltage_count_;
            reactive_row_[bus] = nb_ + voltage_count_;
            ++voltage_count_;
        }
        const int component_slack_offset = angle_count_ + voltage_count_;
        dimension_ = component_slack_offset + component_count_;
        if (nb_ + voltage_count_ != dimension_) {
            failure_reason_ =
                "reduced active-network basis is not square";
            return;
        }

        std::vector<Eigen::Triplet<double>> entries;
        entries.reserve(data.branches.size() * 16 +
                        static_cast<std::size_t>(component_count_));
        const auto append_angle = [&](int row, int bus, double value) {
            const int column = angle_column_[bus];
            if (column >= 0 && std::abs(value) > 1e-14) {
                entries.emplace_back(row, column, value);
            }
        };
        const auto append_voltage = [&](int row, int bus, double value) {
            const int column = voltage_column_[bus];
            if (column >= 0 && std::abs(value) > 1e-14) {
                entries.emplace_back(row, column, value);
            }
        };
        for (int index = 0;
             index < static_cast<int>(data.branches.size()); ++index) {
            const auto& branch = data.branches[index];
            if (!branch.present || branch.status == 0) {
                continue;
            }
            append_angle(
                branch.from, branch.from,
                -derivative.pf_angle[index]);
            append_angle(
                branch.from, branch.to,
                derivative.pf_angle[index]);
            append_voltage(
                branch.from, branch.from,
                -derivative.pf_vf[index]);
            append_voltage(
                branch.from, branch.to,
                -derivative.pf_vt[index]);
            append_angle(
                branch.to, branch.from,
                -derivative.pt_angle[index]);
            append_angle(
                branch.to, branch.to,
                derivative.pt_angle[index]);
            append_voltage(
                branch.to, branch.from,
                -derivative.pt_vf[index]);
            append_voltage(
                branch.to, branch.to,
                -derivative.pt_vt[index]);

            if (reactive_row_[branch.from] >= 0) {
                const int row = reactive_row_[branch.from];
                append_angle(
                    row, branch.from,
                    -derivative.qf_angle[index]);
                append_angle(
                    row, branch.to,
                    derivative.qf_angle[index]);
                append_voltage(
                    row, branch.from,
                    -derivative.qf_vf[index]);
                append_voltage(
                    row, branch.to,
                    -derivative.qf_vt[index]);
            }
            if (reactive_row_[branch.to] >= 0) {
                const int row = reactive_row_[branch.to];
                append_angle(
                    row, branch.from,
                    -derivative.qt_angle[index]);
                append_angle(
                    row, branch.to,
                    derivative.qt_angle[index]);
                append_voltage(
                    row, branch.from,
                    -derivative.qt_vf[index]);
                append_voltage(
                    row, branch.to,
                    -derivative.qt_vt[index]);
            }
        }
        for (int bus = 0; bus < nb_; ++bus) {
            double gs = 0.0;
            double bs = 0.0;
            for (int shunt : data.buses[bus].shunts) {
                gs += data.shunts[shunt].gs;
                bs += effective_shunt_susceptance(
                    data, reference, shunt);
            }
            append_voltage(
                bus, bus, -2.0 * gs * reference.vm[bus]);
            if (reactive_row_[bus] >= 0) {
                append_voltage(
                    reactive_row_[bus], bus,
                    2.0 * bs * reference.vm[bus]);
            }
        }
        for (int component = 0; component < component_count_; ++component) {
            entries.emplace_back(
                component_reference[component],
                component_slack_offset + component, 1.0);
        }
        Eigen::SparseMatrix<double> basis(dimension_, dimension_);
        basis.setFromTriplets(entries.begin(), entries.end());
        basis.makeCompressed();
        basis_factorization_.analyzePattern(basis);
        basis_factorization_.factorize(basis);
        if (basis_factorization_.info() != Eigen::Success) {
            failure_reason_ =
                "reduced active-network basis factorization failed";
            return;
        }
        Eigen::SparseMatrix<double> transpose = basis.transpose();
        transpose.makeCompressed();
        adjoint_factorization_.analyzePattern(transpose);
        adjoint_factorization_.factorize(transpose);
        if (adjoint_factorization_.info() != Eigen::Success) {
            failure_reason_ =
                "reduced active-network adjoint factorization failed";
            return;
        }

        for (int generator = 0; generator < ng_; ++generator) {
            control_entries_[generator].emplace_back(
                data.generators[generator].bus, 1.0);
            if (!pv_pq_partition_) {
                control_entries_[ng_ + nd_ + generator].emplace_back(
                    reactive_row_[data.generators[generator].bus], 1.0);
            }
        }
        for (int load = 0; load < nd_; ++load) {
            const int bus = data.loads[load].bus;
            control_entries_[ng_ + load].emplace_back(
                bus, -data.loads[load].pd_nominal);
            if (reactive_row_[bus] >= 0) {
                control_entries_[ng_ + load].emplace_back(
                    reactive_row_[bus], -data.loads[load].qd_nominal);
            }
        }
        valid_ = true;
    }

    bool valid() const { return valid_; }
    const std::string& failure_reason() const { return failure_reason_; }
    int control_count() const {
        return ng_ + nd_ + (pv_pq_partition_ ? 0 : ng_);
    }
    const std::vector<unsigned char>& pv_bus_mask() const {
        return pv_bus_;
    }

    bool append_component_balance_rows(std::vector<SparseRow>& rows) {
        for (int component = 0; component < component_count_; ++component) {
            Eigen::VectorXd functional =
                Eigen::VectorXd::Zero(dimension_);
            functional[angle_count_ + voltage_count_ + component] = 1.0;
            SparseRow row;
            row.lower = 0.0;
            row.upper = 0.0;
            if (!append_response(row, functional)) {
                return false;
            }
            normalize_row(row);
            rows.push_back(std::move(row));
        }
        return true;
    }

    bool append_angle_response(
        SparseRow& row,
        const std::vector<std::pair<int, double>>& bus_coefficients) {
        return append_state_response(row, bus_coefficients, {});
    }

    bool append_state_response(
        SparseRow& row,
        const std::vector<std::pair<int, double>>& angle_coefficients,
        const std::vector<std::pair<int, double>>& voltage_coefficients) {
        Eigen::VectorXd functional =
            Eigen::VectorXd::Zero(dimension_);
        for (const auto& [bus, coefficient] : angle_coefficients) {
            if (bus < 0 || bus >= nb_ || !std::isfinite(coefficient)) {
                return false;
            }
            const int column = angle_column_[bus];
            if (column >= 0) {
                functional[column] += coefficient;
            }
        }
        for (const auto& [bus, coefficient] : voltage_coefficients) {
            if (bus < 0 || bus >= nb_ || !std::isfinite(coefficient)) {
                return false;
            }
            const int column = voltage_column_[bus];
            if (column >= 0) {
                functional[column] += coefficient;
            }
        }
        return append_response(row, functional);
    }

    bool reconstruct(
        const std::vector<double>& control_values,
        std::vector<double>& angle_change,
        std::vector<double>& voltage_change,
        std::vector<double>& component_slack) {
        if (!valid_ ||
            control_values.size() !=
                static_cast<std::size_t>(control_count())) {
            return false;
        }
        Eigen::VectorXd right_hand_side =
            Eigen::VectorXd::Zero(dimension_);
        for (int control = 0; control < control_count(); ++control) {
            for (const auto& [row, coefficient] :
                 control_entries_[control]) {
                right_hand_side[row] -=
                    coefficient * control_values[control];
            }
        }
        const Eigen::VectorXd solution =
            basis_factorization_.solve(right_hand_side);
        if (basis_factorization_.info() != Eigen::Success ||
            !solution.allFinite()) {
            return false;
        }
        angle_change.assign(static_cast<std::size_t>(nb_), 0.0);
        for (int bus = 0; bus < nb_; ++bus) {
            if (angle_column_[bus] >= 0) {
                angle_change[bus] = solution[angle_column_[bus]];
            }
        }
        voltage_change.assign(static_cast<std::size_t>(nb_), 0.0);
        for (int bus = 0; bus < nb_; ++bus) {
            if (voltage_column_[bus] >= 0) {
                voltage_change[bus] = solution[voltage_column_[bus]];
            }
        }
        component_slack.resize(
            static_cast<std::size_t>(component_count_));
        for (int component = 0; component < component_count_; ++component) {
            component_slack[component] =
                solution[angle_count_ + voltage_count_ + component];
        }
        return true;
    }

private:
    bool append_response(SparseRow& row, const Eigen::VectorXd& functional) {
        if (!valid_ || functional.size() != dimension_) {
            return false;
        }
        const Eigen::VectorXd adjoint =
            adjoint_factorization_.solve(functional);
        if (adjoint_factorization_.info() != Eigen::Success ||
            !adjoint.allFinite()) {
            return false;
        }
        for (int control = 0; control < control_count(); ++control) {
            double response = 0.0;
            for (const auto& [row_index, coefficient] :
                 control_entries_[control]) {
                response -= adjoint[row_index] * coefficient;
            }
            append(row, control, response);
        }
        return true;
    }

    int nb_{};
    int ng_{};
    int nd_{};
    int component_count_{};
    int angle_count_{};
    int voltage_count_{};
    int dimension_{};
    bool pv_pq_partition_{};
    bool valid_{};
    std::string failure_reason_;
    std::vector<int> angle_column_;
    std::vector<int> voltage_column_;
    std::vector<int> reactive_row_;
    std::vector<unsigned char> pv_bus_;
    std::vector<std::vector<std::pair<int, double>>> control_entries_;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> basis_factorization_;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> adjoint_factorization_;
};

double pv_reactive_requirement(
    const CaseData& data,
    const AcState& state,
    int bus) {
    double required = 0.0;
    for (int branch : data.buses[bus].branches_from) {
        required += state.qf[branch];
    }
    for (int branch : data.buses[bus].branches_to) {
        required += state.qt[branch];
    }
    for (int load : data.buses[bus].loads) {
        required += data.loads[load].qd_nominal *
            state.demand_factor[load];
    }
    for (int shunt : data.buses[bus].shunts) {
        required -= effective_shunt_susceptance(
            data, state, shunt) * state.vm[bus] * state.vm[bus];
    }
    return required;
}

double linearized_pv_reactive_requirement_change(
    const CaseData& data,
    const AcState& reference,
    const BranchFlowDerivatives& derivative,
    const std::vector<double>& angle_change,
    const std::vector<double>& voltage_change,
    const std::vector<double>& control_values,
    int bus) {
    const int ng = static_cast<int>(data.generators.size());
    const int nd = static_cast<int>(data.loads.size());
    if (angle_change.size() != data.buses.size() ||
        voltage_change.size() != data.buses.size() ||
        control_values.size() < static_cast<std::size_t>(ng + nd)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double change = 0.0;
    for (int branch_index : data.buses[bus].branches_from) {
        const auto& branch = data.branches[branch_index];
        change += derivative.qf_angle[branch_index] *
            (angle_change[branch.from] - angle_change[branch.to]);
        change += derivative.qf_vf[branch_index] *
            voltage_change[branch.from];
        change += derivative.qf_vt[branch_index] *
            voltage_change[branch.to];
    }
    for (int branch_index : data.buses[bus].branches_to) {
        const auto& branch = data.branches[branch_index];
        change += derivative.qt_angle[branch_index] *
            (angle_change[branch.from] - angle_change[branch.to]);
        change += derivative.qt_vf[branch_index] *
            voltage_change[branch.from];
        change += derivative.qt_vt[branch_index] *
            voltage_change[branch.to];
    }
    for (int load : data.buses[bus].loads) {
        change += data.loads[load].qd_nominal *
            control_values[ng + load];
    }
    double shunt_susceptance = 0.0;
    for (int shunt : data.buses[bus].shunts) {
        shunt_susceptance += effective_shunt_susceptance(
            data, reference, shunt);
    }
    change -= 2.0 * shunt_susceptance * reference.vm[bus] *
        voltage_change[bus];
    return change;
}

bool append_pv_reactive_capability_row(
    const CaseData& data,
    const std::vector<int>& commitment,
    const AcState& reference,
    const BranchFlowDerivatives& derivative,
    ReducedActiveNetwork& reduced,
    int bus,
    SparseRow& row) {
    if (bus < 0 || bus >= static_cast<int>(data.buses.size()) ||
        commitment.size() != data.generators.size() ||
        reference.qf.size() != data.branches.size() ||
        reference.qt.size() != data.branches.size()) {
        return false;
    }
    std::vector<std::pair<int, double>> angle_coefficients;
    std::vector<std::pair<int, double>> voltage_coefficients;
    angle_coefficients.reserve(
        2 * (data.buses[bus].branches_from.size() +
             data.buses[bus].branches_to.size()));
    voltage_coefficients.reserve(angle_coefficients.capacity() + 1);
    for (int branch_index : data.buses[bus].branches_from) {
        const auto& branch = data.branches[branch_index];
        angle_coefficients.emplace_back(
            branch.from, derivative.qf_angle[branch_index]);
        angle_coefficients.emplace_back(
            branch.to, -derivative.qf_angle[branch_index]);
        voltage_coefficients.emplace_back(
            branch.from, derivative.qf_vf[branch_index]);
        voltage_coefficients.emplace_back(
            branch.to, derivative.qf_vt[branch_index]);
    }
    for (int branch_index : data.buses[bus].branches_to) {
        const auto& branch = data.branches[branch_index];
        angle_coefficients.emplace_back(
            branch.from, derivative.qt_angle[branch_index]);
        angle_coefficients.emplace_back(
            branch.to, -derivative.qt_angle[branch_index]);
        voltage_coefficients.emplace_back(
            branch.from, derivative.qt_vf[branch_index]);
        voltage_coefficients.emplace_back(
            branch.to, derivative.qt_vt[branch_index]);
    }
    double shunt_susceptance = 0.0;
    for (int shunt : data.buses[bus].shunts) {
        shunt_susceptance += effective_shunt_susceptance(
            data, reference, shunt);
    }
    voltage_coefficients.emplace_back(
        bus, -2.0 * shunt_susceptance * reference.vm[bus]);
    if (!reduced.append_state_response(
            row, angle_coefficients, voltage_coefficients)) {
        return false;
    }
    const int demand_offset = static_cast<int>(data.generators.size());
    for (int load : data.buses[bus].loads) {
        append(
            row, demand_offset + load,
            data.loads[load].qd_nominal);
    }

    double lower = 0.0;
    double upper = 0.0;
    for (int generator : data.buses[bus].generators) {
        if (commitment[generator] == 0) {
            continue;
        }
        lower += data.generators[generator].qmin;
        upper += data.generators[generator].qmax;
    }
    const double required = pv_reactive_requirement(data, reference, bus);
    // Keep the verified zero-movement incumbent feasible even when it uses a
    // source-authorized local Q-balance slack.  Once inside aggregate Q
    // capability, do not let a linearized economic direction leave it.
    lower = std::min(lower, required);
    upper = std::max(upper, required);
    row.lower = lower - required;
    row.upper = upper - required;
    normalize_row(row);
    return true;
}

struct ReactiveGenerationRecovery {
    int pv_bus_count{};
    int saturated_bus_count{};
    double maximum_unserved_requirement{};
};

// Once voltage and angle have been proposed, reactive generation at a PV bus
// is a dependent power-flow quantity. Recover it exactly when the aggregate
// source Q limits permit. This keeps free PV-bus Q movement from being
// mistaken for reactive-balance slack by the nonlinear candidate checker.
ReactiveGenerationRecovery recover_pv_bus_reactive_generation(
    const CaseData& data,
    const std::vector<int>& commitment,
    AcState& state,
    const std::vector<unsigned char>* pv_bus_mask = nullptr) {
    if (pv_bus_mask != nullptr &&
        pv_bus_mask->size() != data.buses.size()) {
        throw std::runtime_error(
            "PV reactive-recovery mask has invalid dimensions");
    }
    ReactiveGenerationRecovery output;
    for (int generator = 0;
         generator < static_cast<int>(data.generators.size()); ++generator) {
        if (commitment[generator] == 0) {
            state.qg[generator] = 0.0;
        }
    }
    for (int bus = 0;
         bus < static_cast<int>(data.buses.size()); ++bus) {
        if (pv_bus_mask != nullptr && (*pv_bus_mask)[bus] == 0) {
            continue;
        }
        std::vector<int> active_generators;
        double lower_total = 0.0;
        double upper_total = 0.0;
        for (int generator : data.buses[bus].generators) {
            if (commitment[generator] == 0) {
                continue;
            }
            active_generators.push_back(generator);
            lower_total += data.generators[generator].qmin;
            upper_total += data.generators[generator].qmax;
        }
        if (active_generators.empty()) {
            continue;
        }
        ++output.pv_bus_count;
        double required_total = 0.0;
        for (int branch : data.buses[bus].branches_from) {
            required_total += state.qf[branch];
        }
        for (int branch : data.buses[bus].branches_to) {
            required_total += state.qt[branch];
        }
        for (int load : data.buses[bus].loads) {
            required_total += data.loads[load].qd_nominal *
                state.demand_factor[load];
        }
        for (int shunt : data.buses[bus].shunts) {
            required_total -= effective_shunt_susceptance(
                data, state, shunt) * state.vm[bus] * state.vm[bus];
        }
        const double served_total = std::clamp(
            required_total, lower_total, upper_total);
        const double unmet = std::abs(required_total - served_total);
        output.maximum_unserved_requirement = std::max(
            output.maximum_unserved_requirement, unmet);
        if (unmet > 1e-8) {
            ++output.saturated_bus_count;
        }

        double current_total = 0.0;
        for (int generator : active_generators) {
            state.qg[generator] = std::clamp(
                state.qg[generator],
                data.generators[generator].qmin,
                data.generators[generator].qmax);
            current_total += state.qg[generator];
        }
        double remaining = served_total - current_total;
        for (int generator : active_generators) {
            const auto& source = data.generators[generator];
            const double movement = remaining >= 0.0
                ? std::min(remaining, source.qmax - state.qg[generator])
                : std::max(remaining, source.qmin - state.qg[generator]);
            state.qg[generator] += movement;
            remaining -= movement;
            if (std::abs(remaining) <= 1e-12) {
                break;
            }
        }
        output.maximum_unserved_requirement = std::max(
            output.maximum_unserved_requirement,
            std::abs(remaining));
    }
    return output;
}

std::vector<int> select_component_reference_buses(
    const CaseData& data,
    const std::vector<int>& component_of_bus,
    int component_count) {
    if (component_count <= 0 ||
        component_of_bus.size() != data.buses.size()) {
        throw std::runtime_error(
            "invalid active-network component reference dimensions");
    }
    std::vector<int> references(
        static_cast<std::size_t>(component_count), -1);
    std::vector<unsigned char> has_source_reference(
        static_cast<std::size_t>(component_count), 0);
    for (int bus = 0; bus < static_cast<int>(data.buses.size()); ++bus) {
        const int component = component_of_bus[bus];
        if (component < 0 || component >= component_count) {
            throw std::runtime_error(
                "invalid active-network component assignment");
        }
        if (references[component] < 0) {
            references[component] = bus;
        }
        // The GO2 model fixes every source type-3 bus angle to zero. Use the
        // first such bus as the eliminated gauge for its component instead of
        // an arbitrary first-index bus. Additional type-3 buses are retained
        // as explicit zero-response rows below.
        if (data.buses[bus].type == 3 &&
            has_source_reference[component] == 0) {
            references[component] = bus;
            has_source_reference[component] = 1;
        }
    }
    for (int reference : references) {
        if (reference < 0) {
            throw std::runtime_error(
                "active-network component has no reference bus");
        }
    }
    return references;
}

}  // namespace

void run_active_network_reduction_regression() {
    CaseData data;
    data.buses.resize(3);
    data.branches.resize(2);
    for (int bus = 0; bus < 3; ++bus) {
        data.buses[bus].index = bus;
        data.buses[bus].vmin = 0.9;
        data.buses[bus].vmax = 1.1;
    }
    for (int branch = 0; branch < 2; ++branch) {
        data.branches[branch].present = true;
        data.branches[branch].status = 1;
        data.branches[branch].from = branch;
        data.branches[branch].to = branch + 1;
        data.branches[branch].r = 0.0;
        data.branches[branch].x = branch == 0 ? 0.1 : 0.125;
        data.branches[branch].tap = 1.0;
    }
    data.generators.resize(1);
    data.generators[0].bus = 0;
    data.generators[0].qmin = -20.0;
    data.generators[0].qmax = 20.0;
    data.buses[0].generators = {0};
    data.loads.resize(1);
    data.loads[0].bus = 1;
    data.loads[0].pd_nominal = 2.0;
    data.loads[0].qd_nominal = 0.8;
    data.buses[1].type = 3;

    AcState reference;
    reference.vm = {1.02, 0.98, 1.01};
    reference.va = {0.07, -0.03, 0.01};
    reference.qg = {0.0};
    const auto derivative = branch_flow_derivatives(data, reference);
    const auto exact_flows = [&](int branch_index, const AcState& state) {
        const auto& branch = data.branches[branch_index];
        const auto coefficient = flow_coefficients(branch);
        const double vf = state.vm[branch.from];
        const double vt = state.vm[branch.to];
        const double angle =
            state.va[branch.from] - state.va[branch.to];
        const double cosine = std::cos(angle);
        const double sine = std::sin(angle);
        return std::array<double, 4>{
            coefficient.from_g_self * vf * vf + vf * vt *
                (coefficient.from_cross_cos * cosine +
                 coefficient.from_cross_sin * sine),
            -coefficient.from_b_self * vf * vf + vf * vt *
                (-coefficient.from_cross_sin * cosine +
                 coefficient.from_cross_cos * sine),
            coefficient.to_g_self * vt * vt + vf * vt *
                (coefficient.to_cross_cos * cosine -
                 coefficient.to_cross_sin * sine),
            -coefficient.to_b_self * vt * vt + vf * vt *
                (-coefficient.to_cross_sin * cosine -
                 coefficient.to_cross_cos * sine),
        };
    };
    constexpr double kDifferenceStep = 1e-6;
    for (int branch_index = 0; branch_index < 2; ++branch_index) {
        const auto& branch = data.branches[branch_index];
        const auto finite_difference = [&](bool angle, int bus) {
            AcState lower = reference;
            AcState upper = reference;
            auto& lower_values = angle ? lower.va : lower.vm;
            auto& upper_values = angle ? upper.va : upper.vm;
            lower_values[bus] -= kDifferenceStep;
            upper_values[bus] += kDifferenceStep;
            const auto lower_flow = exact_flows(branch_index, lower);
            const auto upper_flow = exact_flows(branch_index, upper);
            std::array<double, 4> result{};
            for (int flow = 0; flow < 4; ++flow) {
                result[flow] =
                    (upper_flow[flow] - lower_flow[flow]) /
                    (2.0 * kDifferenceStep);
            }
            return result;
        };
        const auto angle_from = finite_difference(true, branch.from);
        const auto angle_to = finite_difference(true, branch.to);
        const auto voltage_from = finite_difference(false, branch.from);
        const auto voltage_to = finite_difference(false, branch.to);
        const std::array<double, 4> analytic_angle{
            derivative.pf_angle[branch_index],
            derivative.qf_angle[branch_index],
            derivative.pt_angle[branch_index],
            derivative.qt_angle[branch_index],
        };
        const std::array<double, 4> analytic_voltage_from{
            derivative.pf_vf[branch_index],
            derivative.qf_vf[branch_index],
            derivative.pt_vf[branch_index],
            derivative.qt_vf[branch_index],
        };
        const std::array<double, 4> analytic_voltage_to{
            derivative.pf_vt[branch_index],
            derivative.qf_vt[branch_index],
            derivative.pt_vt[branch_index],
            derivative.qt_vt[branch_index],
        };
        for (int flow = 0; flow < 4; ++flow) {
            if (std::abs(angle_from[flow] - analytic_angle[flow]) > 1e-8 ||
                std::abs(angle_to[flow] + analytic_angle[flow]) > 1e-8 ||
                std::abs(
                    voltage_from[flow] -
                    analytic_voltage_from[flow]) > 1e-8 ||
                std::abs(
                    voltage_to[flow] -
                    analytic_voltage_to[flow]) > 1e-8) {
                throw std::runtime_error(
                    "active-network reduction regression: derivative failed");
            }
        }
    }
    const auto component_references = select_component_reference_buses(
        data, {0, 0, 0}, 1);
    if (component_references != std::vector<int>{1}) {
        throw std::runtime_error(
            "active-network reduction regression: source reference not selected");
    }
    ReducedActiveNetwork reduced(
        data, {1}, reference, component_references, derivative, false);
    if (!reduced.valid() || reduced.control_count() != 3) {
        throw std::runtime_error(
            "active-network reduction regression: construction failed: " +
            reduced.failure_reason());
    }
    std::vector<SparseRow> component_rows;
    if (!reduced.append_component_balance_rows(component_rows) ||
        component_rows.size() != 1) {
        throw std::runtime_error(
            "active-network reduction regression: component row failed");
    }
    // +0.5 generator movement balances a +0.25 load-factor movement at a
    // 2.0-pu load. The associated 0.2-pu reactive movement exercises the PQ
    // voltage rows while the component active slack remains exactly zero.
    const std::vector<double> balanced_controls{0.5, 0.25, 0.0};
    std::vector<double> angle_change;
    std::vector<double> voltage_change;
    std::vector<double> component_slack;
    if (!reduced.reconstruct(
            balanced_controls, angle_change, voltage_change,
            component_slack) ||
        component_slack.size() != 1 ||
        std::abs(component_slack[0]) > 1e-12) {
        throw std::runtime_error(
            "active-network reduction regression: balance reconstruction failed");
    }
    if (std::abs(angle_change[1]) > 1e-12) {
        throw std::runtime_error(
            "active-network reduction regression: source reference moved");
    }
    ReducedActiveNetwork pv_reduced(
        data, {1}, reference, component_references, derivative, true);
    std::vector<double> pv_angle_change;
    std::vector<double> pv_voltage_change;
    std::vector<double> pv_component_slack;
    if (!pv_reduced.valid() || pv_reduced.control_count() != 2 ||
        pv_reduced.pv_bus_mask()[0] == 0 ||
        !pv_reduced.reconstruct(
            {0.5, 0.25}, pv_angle_change, pv_voltage_change,
            pv_component_slack) ||
        std::abs(pv_angle_change[1]) > 1e-12 ||
        std::abs(pv_voltage_change[0]) > 1e-12 ||
        std::abs(pv_component_slack[0]) > 1e-12) {
        throw std::runtime_error(
            "active-network reduction regression: PV/PQ partition failed: " +
            pv_reduced.failure_reason());
    }
    double component_activity = 0.0;
    for (const auto& [column, coefficient] : component_rows[0].entries) {
        component_activity += coefficient * balanced_controls[column];
    }
    if (std::abs(component_activity) > 1e-12) {
        throw std::runtime_error(
            "active-network reduction regression: component sensitivity failed");
    }
    SparseRow angle_difference;
    if (!reduced.append_angle_response(
            angle_difference, {{1, 1.0}, {2, -1.0}})) {
        throw std::runtime_error(
            "active-network reduction regression: adjoint solve failed");
    }
    double adjoint_angle_difference = 0.0;
    for (const auto& [column, coefficient] : angle_difference.entries) {
        adjoint_angle_difference += coefficient * balanced_controls[column];
    }
    const double reconstructed_angle_difference =
        angle_change[1] - angle_change[2];
    if (std::abs(
            adjoint_angle_difference -
            reconstructed_angle_difference) > 1e-12) {
        throw std::runtime_error(
            "active-network reduction regression: adjoint identity failed");
    }
    SparseRow voltage_difference;
    if (!reduced.append_state_response(
            voltage_difference, {}, {{1, 1.0}, {2, -1.0}})) {
        throw std::runtime_error(
            "active-network reduction regression: voltage adjoint failed");
    }
    double adjoint_voltage_difference = 0.0;
    for (const auto& [column, coefficient] : voltage_difference.entries) {
        adjoint_voltage_difference +=
            coefficient * balanced_controls[column];
    }
    const double reconstructed_voltage_difference =
        voltage_change[1] - voltage_change[2];
    if (std::abs(reconstructed_voltage_difference) <= 1e-10 ||
        std::abs(
            adjoint_voltage_difference -
            reconstructed_voltage_difference) > 1e-12) {
        throw std::runtime_error(
            "active-network reduction regression: voltage identity failed");
    }

    data.buses[0].branches_from = {0};
    data.buses[1].branches_to = {0};
    data.buses[1].branches_from = {1};
    data.buses[2].branches_to = {1};
    data.buses[1].loads = {0};
    AcState reactive_state = reference;
    reactive_state.pg = {0.0};
    reactive_state.qg = {0.0};
    reactive_state.demand_factor = {1.0};
    reactive_state.qf.resize(2);
    reactive_state.qt.resize(2);
    for (int branch = 0; branch < 2; ++branch) {
        const auto flows = exact_flows(branch, reactive_state);
        reactive_state.qf[branch] = flows[1];
        reactive_state.qt[branch] = flows[3];
    }
    const std::vector<double> pv_controls{0.5, 0.25};
    const double linearized_reactive_change =
        linearized_pv_reactive_requirement_change(
            data, reactive_state, derivative,
            pv_angle_change, pv_voltage_change, pv_controls, 0);
    const auto reactive_requirement_at_step = [&](double step) {
        AcState perturbed = reactive_state;
        for (int bus = 0; bus < 3; ++bus) {
            perturbed.va[bus] += step * pv_angle_change[bus];
            perturbed.vm[bus] += step * pv_voltage_change[bus];
        }
        perturbed.demand_factor[0] += step * pv_controls[1];
        for (int branch = 0; branch < 2; ++branch) {
            const auto flows = exact_flows(branch, perturbed);
            perturbed.qf[branch] = flows[1];
            perturbed.qt[branch] = flows[3];
        }
        return pv_reactive_requirement(data, perturbed, 0);
    };
    const double finite_difference_reactive_change =
        (reactive_requirement_at_step(kDifferenceStep) -
         reactive_requirement_at_step(-kDifferenceStep)) /
        (2.0 * kDifferenceStep);
    SparseRow reactive_capability;
    if (!std::isfinite(linearized_reactive_change) ||
        std::abs(
            linearized_reactive_change -
            finite_difference_reactive_change) > 1e-7 ||
        !append_pv_reactive_capability_row(
            data, {1}, reactive_state, derivative,
            pv_reduced, 0, reactive_capability)) {
        throw std::runtime_error(
            "active-network reduction regression: PV Q sensitivity failed");
    }
    double reactive_row_activity = 0.0;
    for (const auto& [column, coefficient] :
         reactive_capability.entries) {
        reactive_row_activity += coefficient * pv_controls[column];
    }
    const double normalized_width =
        reactive_capability.upper - reactive_capability.lower;
    const double reactive_row_scale = 40.0 / normalized_width;
    if (!std::isfinite(reactive_row_scale) ||
        std::abs(
            reactive_row_scale * reactive_row_activity -
            linearized_reactive_change) > 1e-10) {
        throw std::runtime_error(
            "active-network reduction regression: PV Q row identity failed");
    }
    const auto recovery = recover_pv_bus_reactive_generation(
        data, {1}, reactive_state);
    if (recovery.pv_bus_count != 1 ||
        recovery.saturated_bus_count != 0 ||
        recovery.maximum_unserved_requirement > 1e-12 ||
        std::abs(reactive_state.qg[0] - reactive_state.qf[0]) > 1e-12) {
        throw std::runtime_error(
            "active-network reduction regression: PV Q recovery failed");
    }
}

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
        {"angle_response_row_count", angle_response_row_count},
        {"voltage_response_row_count", voltage_response_row_count},
        {"reactive_capability_row_count",
         reactive_capability_row_count},
        {"trust_region_constraint_generation_passes",
         trust_region_constraint_generation_passes},
        {"voltage_constraint_generation_passes",
         voltage_constraint_generation_passes},
        {"reactive_capability_constraint_generation_passes",
         reactive_capability_constraint_generation_passes},
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
    if (options.contingency_profile) {
        CaseData corrective_profile = data;
        corrective_profile.delta = data.delta_ctg;
        corrective_profile.delta_r = data.delta_r_ctg;
        for (auto& generator : corrective_profile.generators) {
            generator.prumax = generator.prumaxctg;
            generator.prdmax = generator.prdmaxctg;
        }
        for (auto& load : corrective_profile.loads) {
            load.prumax = load.prumaxctg;
            load.prdmax = load.prdmaxctg;
        }
        for (auto& branch : corrective_profile.branches) {
            branch.rate_a = branch.rate_c;
        }
        auto corrective_options = options;
        corrective_options.contingency_profile = false;
        const char* presolve_override =
            std::getenv("GRAVITYX_ACTIVE_ECONOMIC_PRESOLVE");
        corrective_options.simplex_presolve =
            presolve_override != nullptr &&
            std::string(presolve_override) != "0";
        const char* crash_basis_override =
            std::getenv("GRAVITYX_ACTIVE_ECONOMIC_CRASH_BASIS");
        corrective_options.use_simplex_crash_basis =
            crash_basis_override != nullptr &&
            std::string(crash_basis_override) != "0";
        const char* compact_override =
            std::getenv("GRAVITYX_ACTIVE_ECONOMIC_COMPACT_SIGNED");
        corrective_options.compact_signed_columns =
            compact_override == nullptr ||
            std::string(compact_override) != "0";
        const char* load_movement_override =
            std::getenv("GRAVITYX_ACTIVE_ECONOMIC_LOAD_MOVEMENT");
        corrective_options.freeze_load_movement =
            load_movement_override == nullptr ||
            std::string(load_movement_override) == "0";
        const char* eliminate_angles_override =
            std::getenv("GRAVITYX_ACTIVE_ECONOMIC_ELIMINATE_ANGLES");
        corrective_options.eliminate_angles =
            eliminate_angles_override != nullptr &&
            std::string(eliminate_angles_override) != "0";
        const char* reduced_pv_pq_override =
            std::getenv("GRAVITYX_ACTIVE_ECONOMIC_REDUCED_PV_PQ");
        corrective_options.reduced_pv_pq_partition =
            reduced_pv_pq_override != nullptr &&
            std::string(reduced_pv_pq_override) != "0";
        const char* trust_row_generation_override =
            std::getenv("GRAVITYX_ACTIVE_ECONOMIC_GENERATE_TRUST_ROWS");
        corrective_options.generate_trust_region_rows =
            trust_row_generation_override != nullptr &&
            std::string(trust_row_generation_override) != "0";
        const char* maximum_rounds_override =
            std::getenv("GRAVITYX_ACTIVE_ECONOMIC_MAX_ROUNDS");
        if (maximum_rounds_override != nullptr) {
            corrective_options.maximum_rounds =
                std::max(1, std::stoi(maximum_rounds_override));
        }
        const char* trust_rows_override =
            std::getenv("GRAVITYX_ACTIVE_ECONOMIC_TRUST_ROWS_PER_PASS");
        if (trust_rows_override != nullptr) {
            corrective_options.maximum_trust_region_rows_per_pass =
                std::max(1, std::stoi(trust_rows_override));
        }
        const char* angle_trust_override =
            std::getenv("GRAVITYX_ACTIVE_ECONOMIC_ANGLE_TRUST_RADIUS");
        if (angle_trust_override != nullptr) {
            corrective_options.angle_trust_radius =
                std::stod(angle_trust_override);
        }
        const char* voltage_trust_override =
            std::getenv("GRAVITYX_ACTIVE_ECONOMIC_VOLTAGE_TRUST_RADIUS");
        if (voltage_trust_override != nullptr) {
            corrective_options.voltage_trust_radius =
                std::stod(voltage_trust_override);
        }
        const char* sparse_full_ac_override =
            std::getenv("GRAVITYX_ACTIVE_ECONOMIC_SPARSE_FULL_AC");
        corrective_options.sparse_full_ac_linearization =
            sparse_full_ac_override != nullptr &&
            std::string(sparse_full_ac_override) != "0";
        return refine_fixed_commitment_active_network_economic(
            corrective_profile, commitment, incumbent,
            corrective_options);
    }
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
        !std::isfinite(options.voltage_trust_radius) ||
        options.voltage_trust_radius <= 0.0 ||
        !std::isfinite(
            options.maximum_candidate_repair_balance_slack) ||
        options.maximum_candidate_repair_balance_slack <= 0.0 ||
        !std::isfinite(options.thermal_row_utilization_threshold) ||
        options.thermal_row_utilization_threshold < 0.0 ||
        options.thermal_row_utilization_threshold > 1.0 ||
        options.maximum_rounds <= 0 ||
        options.maximum_candidate_trials <= 0 ||
        options.maximum_candidate_repair_newton_iterations <= 0 ||
        options.maximum_voltage_rows_per_pass <= 0 ||
        options.maximum_trust_region_rows_per_pass <= 0) {
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
    const char* trace_override =
        std::getenv("GRAVITYX_ACTIVE_ECONOMIC_TRACE");
    const bool trace_enabled = trace_override != nullptr &&
        std::string(trace_override) != "0";
    const auto trace = [&](const std::string& stage) {
        if (!trace_enabled) {
            return;
        }
        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        std::cerr << "GRAVITYX_ACTIVE_ECONOMIC_TRACE "
                  << elapsed << " " << stage << std::endl;
    };

    const int nb = static_cast<int>(data.buses.size());
    const int ng = static_cast<int>(data.generators.size());
    const int nd = static_cast<int>(data.loads.size());
    // The original base diagnostic uses nonnegative up/down movement pairs.
    // The common corrective profile instead uses one signed column per
    // control. This halves the column count and removes the severe paired-
    // column degeneracy that prevented HiGHS from retaining a feasible
    // time-limited terminal point on the 19k case.
    const bool compact_signed_columns = options.compact_signed_columns;
    const bool eliminate_angles =
        options.eliminate_angles && compact_signed_columns;
    const bool reduced_pv_pq =
        eliminate_angles && options.reduced_pv_pq_partition;
    const bool sparse_full_ac =
        options.sparse_full_ac_linearization &&
        compact_signed_columns && !eliminate_angles;
    const int pg_signed_offset = 0;
    const int demand_signed_offset = pg_signed_offset + ng;
    const int qg_signed_offset = demand_signed_offset + nd;
    const int angle_signed_offset = qg_signed_offset +
        (((eliminate_angles && !reduced_pv_pq) || sparse_full_ac)
            ? ng : 0);
    const int voltage_signed_offset = angle_signed_offset + nb;
    const int pg_up_offset = 0;
    const int pg_down_offset = pg_up_offset + ng;
    const int demand_up_offset = pg_down_offset + ng;
    const int demand_down_offset = demand_up_offset + nd;
    const int angle_up_offset = demand_down_offset + nd;
    const int angle_down_offset = angle_up_offset + nb;
    const int column_count = compact_signed_columns
        ? (eliminate_angles ? angle_signed_offset
                            : (sparse_full_ac
                                ? voltage_signed_offset + nb
                                : angle_signed_offset + nb))
        : angle_down_offset + nb;
    output.column_count = column_count;
    const auto append_delta = [&](
        SparseRow& row,
        int signed_column,
        int up_column,
        int down_column,
        double coefficient) {
        if (compact_signed_columns) {
            append(row, signed_column, coefficient);
        } else {
            append(row, up_column, coefficient);
            append(row, down_column, -coefficient);
        }
    };

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
    std::vector<int> root_to_component(static_cast<std::size_t>(nb), -1);
    std::vector<int> component_of_bus(static_cast<std::size_t>(nb), -1);
    for (int bus = 0; bus < nb; ++bus) {
        const int root = find_root(bus);
        if (root_to_component[root] < 0) {
            root_to_component[root] = output.component_count++;
        }
        component_of_bus[bus] = root_to_component[root];
    }
    const std::vector<int> component_reference =
        select_component_reference_buses(
            data, component_of_bus, output.component_count);
    std::vector<unsigned char> monitored_angle_branch(
        data.branches.size(), 0);
    std::vector<unsigned char> monitored_voltage_bus(
        static_cast<std::size_t>(nb), 0);
    std::vector<unsigned char> monitored_reactive_bus(
        static_cast<std::size_t>(nb), 0);

    for (int round = 1; round <= options.maximum_rounds; ++round) {
        trace("round=" + std::to_string(round) + ":start");
        const double elapsed_before_round = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - wall_start).count();
        const double remaining =
            options.time_limit_seconds - elapsed_before_round;
        if (remaining <= 0.05) {
            output.time_limit_reached = true;
            break;
        }
        const SolveResult reference = output.selected;
        if (eliminate_angles) {
            constexpr double kActiveVoltageBoundTolerance = 1e-8;
            for (int bus = 0; bus < nb; ++bus) {
                if (reference.state.vm[bus] - data.buses[bus].vmin <=
                        kActiveVoltageBoundTolerance ||
                    data.buses[bus].vmax - reference.state.vm[bus] <=
                        kActiveVoltageBoundTolerance) {
                    monitored_voltage_bus[bus] = 1;
                }
            }
        }

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
            const auto points = active_pwl_points(
                source.cost, source.ncost,
                physical_lower, physical_upper);
            const double slope = data.delta *
                pwl_slope(points, reference.state.pg[generator]);
            if (compact_signed_columns) {
                lower[pg_signed_offset + generator] =
                    physical_lower - reference.state.pg[generator];
                upper[pg_signed_offset + generator] =
                    physical_upper - reference.state.pg[generator];
                cost[pg_signed_offset + generator] = slope;
            } else {
                upper[pg_up_offset + generator] = std::max(
                    0.0,
                    physical_upper - reference.state.pg[generator]);
                upper[pg_down_offset + generator] = std::max(
                    0.0,
                    reference.state.pg[generator] - physical_lower);
                cost[pg_up_offset + generator] = slope;
                cost[pg_down_offset + generator] = -slope;
            }
            if ((eliminate_angles && !reduced_pv_pq) ||
                sparse_full_ac) {
                if (reference.state.qg[generator] < source.qmin - 1e-8 ||
                    reference.state.qg[generator] > source.qmax + 1e-8) {
                    output.status =
                        "reference_generator_outside_reactive_bounds";
                    break;
                }
                lower[qg_signed_offset + generator] =
                    source.qmin - reference.state.qg[generator];
                upper[qg_signed_offset + generator] =
                    source.qmax - reference.state.qg[generator];
            }
        }
        if (output.status == "reference_generator_outside_source_bounds" ||
            output.status ==
                "reference_generator_outside_reactive_bounds") {
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
            const auto points = active_pwl_points(
                source.cost, source.ncost,
                source.pd_min, source.pd_max);
            const double slope = -data.delta *
                source.pd_nominal * pwl_slope(
                    points, source.pd_nominal *
                        reference.state.demand_factor[load]);
            if (compact_signed_columns) {
                lower[demand_signed_offset + load] =
                    options.freeze_load_movement
                    ? 0.0
                    : factor_lower -
                        reference.state.demand_factor[load];
                upper[demand_signed_offset + load] =
                    options.freeze_load_movement
                    ? 0.0
                    : factor_upper -
                        reference.state.demand_factor[load];
                cost[demand_signed_offset + load] = slope;
            } else {
                if (!options.freeze_load_movement) {
                    upper[demand_up_offset + load] = std::max(
                        0.0,
                        factor_upper -
                            reference.state.demand_factor[load]);
                    upper[demand_down_offset + load] = std::max(
                        0.0,
                        reference.state.demand_factor[load] - factor_lower);
                }
                cost[demand_up_offset + load] = slope;
                cost[demand_down_offset + load] = -slope;
            }
        }
        if (output.status == "reference_load_outside_source_bounds") {
            break;
        }
        if (!eliminate_angles) {
            for (int bus = 0; bus < nb; ++bus) {
                if (compact_signed_columns) {
                    lower[angle_signed_offset + bus] =
                        -options.angle_trust_radius;
                    upper[angle_signed_offset + bus] =
                        options.angle_trust_radius;
                } else {
                    upper[angle_up_offset + bus] =
                        options.angle_trust_radius;
                    upper[angle_down_offset + bus] =
                        options.angle_trust_radius;
                }
                if (sparse_full_ac) {
                    lower[voltage_signed_offset + bus] = std::max(
                        -options.voltage_trust_radius,
                        data.buses[bus].vmin - reference.state.vm[bus]);
                    upper[voltage_signed_offset + bus] = std::min(
                        options.voltage_trust_radius,
                        data.buses[bus].vmax - reference.state.vm[bus]);
                }
            }
            for (int bus : component_reference) {
                if (compact_signed_columns) {
                    lower[angle_signed_offset + bus] = 0.0;
                    upper[angle_signed_offset + bus] = 0.0;
                } else {
                    upper[angle_up_offset + bus] = 0.0;
                    upper[angle_down_offset + bus] = 0.0;
                }
            }
            for (int bus = 0; bus < nb; ++bus) {
                if (data.buses[bus].type != 3) {
                    continue;
                }
                if (compact_signed_columns) {
                    lower[angle_signed_offset + bus] = 0.0;
                    upper[angle_signed_offset + bus] = 0.0;
                } else {
                    upper[angle_up_offset + bus] = 0.0;
                    upper[angle_down_offset + bus] = 0.0;
                }
            }
        }
        double objective_scale = 1.0;
        for (double coefficient : cost) {
            objective_scale = std::max(
                objective_scale, std::abs(coefficient));
        }
        for (double& coefficient : cost) {
            coefficient /= objective_scale;
        }

        const auto flow_derivative =
            branch_flow_derivatives(data, reference.state);
        const auto& dpf = flow_derivative.pf_angle;
        const auto& dqf = flow_derivative.qf_angle;
        const auto& dpt = flow_derivative.pt_angle;
        const auto& dqt = flow_derivative.qt_angle;

        std::vector<SparseRow> rows;
        std::unique_ptr<ReducedActiveNetwork> reduced_network;
        bool reduced_row_failure = false;
        int reduced_pv_bus_count = 0;
        if (eliminate_angles) {
            reduced_network = std::make_unique<ReducedActiveNetwork>(
                data, commitment, reference.state,
                component_reference, flow_derivative,
                reduced_pv_pq);
            if (!reduced_network->valid()) {
                output.status = reduced_network->failure_reason();
                break;
            }
            reduced_pv_bus_count = static_cast<int>(std::count(
                reduced_network->pv_bus_mask().begin(),
                reduced_network->pv_bus_mask().end(),
                static_cast<unsigned char>(1)));
            trace("round=" + std::to_string(round) +
                  ":reduced_factorization_complete");
            rows.reserve(static_cast<std::size_t>(
                output.component_count + data.branches.size()));
            if (!reduced_network->append_component_balance_rows(rows)) {
                output.status =
                    "reduced active-network component rows failed";
                break;
            }
            int reactive_capability_rows = 0;
            if (reduced_pv_pq) {
                for (int bus = 0; bus < nb; ++bus) {
                    if (monitored_reactive_bus[bus] == 0 ||
                        reduced_network->pv_bus_mask()[bus] == 0) {
                        continue;
                    }
                    SparseRow capability;
                    if (!append_pv_reactive_capability_row(
                            data, commitment, reference.state,
                            flow_derivative, *reduced_network, bus,
                            capability)) {
                        output.status =
                            "reduced active-network reactive capability row failed";
                        reduced_row_failure = true;
                        break;
                    }
                    rows.push_back(std::move(capability));
                    ++reactive_capability_rows;
                }
            }
            if (reduced_row_failure) {
                break;
            }
            output.reactive_capability_row_count =
                reactive_capability_rows;
            for (int bus = 0; bus < nb; ++bus) {
                if (data.buses[bus].type != 3 ||
                    bus == component_reference[component_of_bus[bus]]) {
                    continue;
                }
                SparseRow reference_angle;
                reference_angle.lower = 0.0;
                reference_angle.upper = 0.0;
                if (!reduced_network->append_angle_response(
                        reference_angle, {{bus, 1.0}})) {
                    output.status =
                        "reduced active-network reference-angle row failed";
                    reduced_row_failure = true;
                    break;
                }
                normalize_row(reference_angle);
                rows.push_back(std::move(reference_angle));
            }
            if (reduced_row_failure) {
                break;
            }
            int angle_response_rows = 0;
            for (int branch_index = 0;
                 branch_index < static_cast<int>(data.branches.size());
                 ++branch_index) {
                if (monitored_angle_branch[branch_index] == 0 ||
                    !data.branches[branch_index].present ||
                    data.branches[branch_index].status == 0) {
                    continue;
                }
                const auto& branch = data.branches[branch_index];
                SparseRow angle;
                angle.lower = -options.angle_trust_radius;
                angle.upper = options.angle_trust_radius;
                if (!reduced_network->append_angle_response(
                        angle,
                        {{branch.from, 1.0}, {branch.to, -1.0}})) {
                    output.status =
                        "reduced active-network angle trust row failed";
                    reduced_row_failure = true;
                    break;
                }
                normalize_row(angle);
                rows.push_back(std::move(angle));
                ++angle_response_rows;
            }
            if (reduced_row_failure) {
                break;
            }
            output.angle_response_row_count = angle_response_rows;
            int voltage_rows = 0;
            for (int bus = 0; bus < nb; ++bus) {
                if (monitored_voltage_bus[bus] == 0) {
                    continue;
                }
                SparseRow voltage;
                voltage.lower = std::max(
                    -options.voltage_trust_radius,
                    data.buses[bus].vmin - reference.state.vm[bus]);
                voltage.upper = std::min(
                    options.voltage_trust_radius,
                    data.buses[bus].vmax - reference.state.vm[bus]);
                if (!reduced_network->append_state_response(
                        voltage, {}, {{bus, 1.0}})) {
                    output.status =
                        "reduced active-network voltage row failed";
                    reduced_row_failure = true;
                    break;
                }
                normalize_row(voltage);
                rows.push_back(std::move(voltage));
                ++voltage_rows;
            }
            if (reduced_row_failure) {
                break;
            }
            output.voltage_response_row_count = voltage_rows;
        } else {
            rows.reserve(static_cast<std::size_t>(
                (sparse_full_ac ? 2 * nb : nb) +
                data.branches.size()));
            for (int bus = 0; bus < nb; ++bus) {
                SparseRow balance;
                SparseRow reactive_balance;
                for (int generator : data.buses[bus].generators) {
                    append_delta(
                        balance, pg_signed_offset + generator,
                        pg_up_offset + generator,
                        pg_down_offset + generator, 1.0);
                    if (sparse_full_ac) {
                        append(
                            reactive_balance,
                            qg_signed_offset + generator, 1.0);
                    }
                }
                for (int load : data.buses[bus].loads) {
                    append_delta(
                        balance, demand_signed_offset + load,
                        demand_up_offset + load,
                        demand_down_offset + load,
                        -data.loads[load].pd_nominal);
                    if (sparse_full_ac) {
                        append(
                            reactive_balance,
                            demand_signed_offset + load,
                            -data.loads[load].qd_nominal);
                    }
                }
                for (int branch_index : data.buses[bus].branches_from) {
                    if (!data.branches[branch_index].present ||
                        data.branches[branch_index].status == 0) {
                        continue;
                    }
                    const auto& branch = data.branches[branch_index];
                    append_delta(
                        balance, angle_signed_offset + branch.from,
                        angle_up_offset + branch.from,
                        angle_down_offset + branch.from,
                        -dpf[branch_index]);
                    append_delta(
                        balance, angle_signed_offset + branch.to,
                        angle_up_offset + branch.to,
                        angle_down_offset + branch.to,
                        dpf[branch_index]);
                    if (sparse_full_ac) {
                        append(
                            balance, voltage_signed_offset + branch.from,
                            -flow_derivative.pf_vf[branch_index]);
                        append(
                            balance, voltage_signed_offset + branch.to,
                            -flow_derivative.pf_vt[branch_index]);
                        append(
                            reactive_balance,
                            angle_signed_offset + branch.from,
                            -dqf[branch_index]);
                        append(
                            reactive_balance,
                            angle_signed_offset + branch.to,
                            dqf[branch_index]);
                        append(
                            reactive_balance,
                            voltage_signed_offset + branch.from,
                            -flow_derivative.qf_vf[branch_index]);
                        append(
                            reactive_balance,
                            voltage_signed_offset + branch.to,
                            -flow_derivative.qf_vt[branch_index]);
                    }
                }
                for (int branch_index : data.buses[bus].branches_to) {
                    if (!data.branches[branch_index].present ||
                        data.branches[branch_index].status == 0) {
                        continue;
                    }
                    const auto& branch = data.branches[branch_index];
                    append_delta(
                        balance, angle_signed_offset + branch.from,
                        angle_up_offset + branch.from,
                        angle_down_offset + branch.from,
                        -dpt[branch_index]);
                    append_delta(
                        balance, angle_signed_offset + branch.to,
                        angle_up_offset + branch.to,
                        angle_down_offset + branch.to,
                        dpt[branch_index]);
                    if (sparse_full_ac) {
                        append(
                            balance, voltage_signed_offset + branch.from,
                            -flow_derivative.pt_vf[branch_index]);
                        append(
                            balance, voltage_signed_offset + branch.to,
                            -flow_derivative.pt_vt[branch_index]);
                        append(
                            reactive_balance,
                            angle_signed_offset + branch.from,
                            -dqt[branch_index]);
                        append(
                            reactive_balance,
                            angle_signed_offset + branch.to,
                            dqt[branch_index]);
                        append(
                            reactive_balance,
                            voltage_signed_offset + branch.from,
                            -flow_derivative.qt_vf[branch_index]);
                        append(
                            reactive_balance,
                            voltage_signed_offset + branch.to,
                            -flow_derivative.qt_vt[branch_index]);
                    }
                }
                if (sparse_full_ac) {
                    double conductance = 0.0;
                    double susceptance = 0.0;
                    for (int shunt : data.buses[bus].shunts) {
                        conductance += data.shunts[shunt].gs;
                        susceptance += effective_shunt_susceptance(
                            data, reference.state, shunt);
                    }
                    append(
                        balance, voltage_signed_offset + bus,
                        -2.0 * conductance * reference.state.vm[bus]);
                    append(
                        reactive_balance, voltage_signed_offset + bus,
                        2.0 * susceptance * reference.state.vm[bus]);
                }
                balance.lower = 0.0;
                balance.upper = 0.0;
                normalize_row(balance);
                rows.push_back(std::move(balance));
                if (sparse_full_ac) {
                    reactive_balance.lower = 0.0;
                    reactive_balance.upper = 0.0;
                    normalize_row(reactive_balance);
                    rows.push_back(std::move(reactive_balance));
                }
            }
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
                if (eliminate_angles) {
                    if (!reduced_network->append_angle_response(
                            angle_limit,
                            {{branch.from, 1.0},
                             {branch.to, -1.0}})) {
                        output.status =
                            "reduced active-network angle row failed";
                        reduced_row_failure = true;
                        break;
                    }
                } else {
                    append_delta(
                        angle_limit, angle_signed_offset + branch.from,
                        angle_up_offset + branch.from,
                        angle_down_offset + branch.from, 1.0);
                    append_delta(
                        angle_limit, angle_signed_offset + branch.to,
                        angle_up_offset + branch.to,
                        angle_down_offset + branch.to, -1.0);
                }
                normalize_row(angle_limit);
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
                double voltage_from_derivative = magnitude > 1e-12
                    ? (p * dp_vf + q * dq_vf) / magnitude
                    : 0.0;
                double voltage_to_derivative = magnitude > 1e-12
                    ? (p * dp_vt + q * dq_vt) / magnitude
                    : 0.0;
                // Non-transformer ratings are RATE_A times the terminal
                // voltage magnitude. Move that first-order right-hand-side
                // response to the left with the flow-magnitude derivative.
                if (!branch.transformer) {
                    if (from_terminal) {
                        voltage_from_derivative -= branch.rate_a;
                    } else {
                        voltage_to_derivative -= branch.rate_a;
                    }
                }
                SparseRow thermal;
                thermal.upper = std::max(limit, magnitude) - magnitude;
                if (eliminate_angles) {
                    if (!reduced_network->append_state_response(
                            thermal,
                            {{branch.from, derivative},
                             {branch.to, -derivative}},
                            {{branch.from, voltage_from_derivative},
                             {branch.to, voltage_to_derivative}})) {
                        output.status =
                            "reduced active-network thermal row failed";
                        reduced_row_failure = true;
                        return;
                    }
                } else {
                    append_delta(
                        thermal, angle_signed_offset + branch.from,
                        angle_up_offset + branch.from,
                        angle_down_offset + branch.from, derivative);
                    append_delta(
                        thermal, angle_signed_offset + branch.to,
                        angle_up_offset + branch.to,
                        angle_down_offset + branch.to, -derivative);
                    if (sparse_full_ac) {
                        append(
                            thermal,
                            voltage_signed_offset + branch.from,
                            voltage_from_derivative);
                        append(
                            thermal,
                            voltage_signed_offset + branch.to,
                            voltage_to_derivative);
                    }
                }
                normalize_row(thermal);
                rows.push_back(std::move(thermal));
                ++thermal_rows;
            };
            add_thermal_row(
                from_magnitude, from_limit,
                reference.state.pf[index], reference.state.qf[index],
                dpf[index], dqf[index],
                flow_derivative.pf_vf[index],
                flow_derivative.qf_vf[index],
                flow_derivative.pf_vt[index],
                flow_derivative.qf_vt[index], true);
            add_thermal_row(
                to_magnitude, to_limit,
                reference.state.pt[index], reference.state.qt[index],
                dpt[index], dqt[index],
                flow_derivative.pt_vf[index],
                flow_derivative.qt_vf[index],
                flow_derivative.pt_vt[index],
                flow_derivative.qt_vt[index], false);
            if (reduced_row_failure) {
                break;
            }
        }
        if (reduced_row_failure) {
            break;
        }
        trace("round=" + std::to_string(round) +
              ":rows_complete");

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
        int highs_threads = 1;
        const char* threads_override =
            std::getenv("GRAVITYX_ACTIVE_ECONOMIC_THREADS");
        if (threads_override != nullptr) {
            highs_threads = std::max(1, std::stoi(threads_override));
        }
        highs.setOptionValue("threads", highs_threads);
        const char* solver_override =
            std::getenv("GRAVITYX_ACTIVE_ECONOMIC_SOLVER");
        const std::string solver = solver_override != nullptr
            ? std::string(solver_override) : "simplex";
        highs.setOptionValue(
            "presolve",
            solver == "simplex" && !options.simplex_presolve
                ? "off" : "on");
        highs.setOptionValue("solver", solver);
        if (solver == "simplex") {
            int simplex_strategy = 4;
            const char* strategy_override =
                std::getenv("GRAVITYX_ACTIVE_ECONOMIC_SIMPLEX_STRATEGY");
            if (strategy_override != nullptr) {
                simplex_strategy = std::stoi(strategy_override);
            }
            highs.setOptionValue("simplex_strategy", simplex_strategy);
            highs.setOptionValue(
                "primal_simplex_bound_perturbation_multiplier", 0.0);
        }
        highs.setOptionValue("primal_feasibility_tolerance", 1e-9);
        highs.setOptionValue("dual_feasibility_tolerance", 1e-9);
        const char* iteration_limit = std::getenv(
            "GRAVITYX_ACTIVE_ECONOMIC_SIMPLEX_ITERATION_LIMIT");
        if (solver == "simplex" && iteration_limit != nullptr) {
            const HighsInt parsed_limit = static_cast<HighsInt>(
                std::stoll(iteration_limit));
            if (parsed_limit > 0) {
                highs.setOptionValue(
                    "simplex_iteration_limit", parsed_limit);
            }
        }
        highs.setOptionValue(
            "time_limit", std::max(
                0.05, std::min(0.5 * remaining, remaining - 0.02)));
        const HighsStatus add_vars_status =
            highs.addVars(column_count, lower.data(), upper.data());
        const HighsStatus change_cost_status = add_vars_status !=
                HighsStatus::kError
            ? highs.changeColsCost(
                  0, column_count - 1, cost.data())
            : HighsStatus::kError;
        const HighsStatus add_rows_status =
            add_vars_status != HighsStatus::kError &&
                change_cost_status != HighsStatus::kError
            ? highs.addRows(
                  static_cast<HighsInt>(rows.size()), row_lower.data(),
                  row_upper.data(), static_cast<HighsInt>(indices.size()),
                  starts.data(), indices.data(), values.data())
            : HighsStatus::kError;
        const bool model_loaded =
            add_vars_status != HighsStatus::kError &&
            change_cost_status != HighsStatus::kError &&
            add_rows_status != HighsStatus::kError;
        if (!model_loaded) {
            output.status =
                "model_construction_failed:add_vars=" +
                std::to_string(static_cast<int>(add_vars_status)) +
                ",change_cost=" +
                std::to_string(static_cast<int>(change_cost_status)) +
                ",add_rows=" +
                std::to_string(static_cast<int>(add_rows_status));
            break;
        }
        std::vector<double> primal_start(
            static_cast<std::size_t>(column_count), 0.0);
        std::vector<HighsInt> start_indices(
            static_cast<std::size_t>(column_count));
        std::iota(start_indices.begin(), start_indices.end(), HighsInt{0});
        bool primal_start_attempted = false;
        HighsStatus start_status = HighsStatus::kOk;
        if (!options.simplex_presolve) {
            primal_start_attempted = true;
            start_status = highs.setSolution(
                column_count, start_indices.data(), primal_start.data());
        }
        HighsStatus basis_status = HighsStatus::kOk;
        bool basis_attempted = false;
        if (solver == "simplex" && !options.simplex_presolve &&
            !sparse_full_ac) {
            basis_attempted = true;
            HighsBasis basis;
            basis.alien = true;
            basis.useful = true;
            basis.col_status.assign(
                static_cast<std::size_t>(column_count),
                HighsBasisStatus::kLower);
            basis.row_status.assign(
                rows.size(), HighsBasisStatus::kBasic);
            if (compact_signed_columns) {
                for (int column = 0; column < column_count; ++column) {
                    if (std::abs(lower[column]) <= 1e-14) {
                        basis.col_status[column] =
                            HighsBasisStatus::kLower;
                    } else if (std::abs(upper[column]) <= 1e-14) {
                        basis.col_status[column] =
                            HighsBasisStatus::kUpper;
                    } else {
                        basis.col_status[column] =
                            HighsBasisStatus::kZero;
                    }
                }
            }
            if (options.use_simplex_crash_basis && !sparse_full_ac) {
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
                        const int angle_column = compact_signed_columns
                            ? angle_signed_offset + bus
                            : angle_up_offset + bus;
                        basis.col_status[angle_column] =
                            HighsBasisStatus::kBasic;
                    }
                }
                std::vector<unsigned char> injection_chosen(
                    static_cast<std::size_t>(output.component_count), 0);
                const auto choose_injection = [&](int component, int column) {
                    if (injection_chosen[component] == 0 &&
                        (lower[column] < -1e-12 ||
                         upper[column] > 1e-12)) {
                        basis.col_status[column] = HighsBasisStatus::kBasic;
                        injection_chosen[component] = 1;
                    }
                };
                for (int generator = 0; generator < ng; ++generator) {
                    const int component = component_of_bus[
                        data.generators[generator].bus];
                    if (compact_signed_columns) {
                        choose_injection(
                            component, pg_signed_offset + generator);
                    } else {
                        choose_injection(
                            component, pg_up_offset + generator);
                        choose_injection(
                            component, pg_down_offset + generator);
                    }
                }
                for (int load = 0; load < nd; ++load) {
                    const int component =
                        component_of_bus[data.loads[load].bus];
                    if (compact_signed_columns) {
                        choose_injection(
                            component, demand_signed_offset + load);
                    } else {
                        choose_injection(
                            component, demand_up_offset + load);
                        choose_injection(
                            component, demand_down_offset + load);
                    }
                }
                const bool complete_basis = std::all_of(
                    injection_chosen.begin(), injection_chosen.end(),
                    [](unsigned char value) { return value != 0; });
                basis_status = complete_basis
                    ? highs.setBasis(
                          basis, "active_network_tree_crash")
                    : HighsStatus::kError;
            } else {
                // Every movement variable is at its lower bound and every row
                // is basic. This is the exact zero-movement feasible point;
                // primal simplex can improve it immediately without first
                // destroying and then rediscovering feasibility.
                basis_status = highs.setBasis(
                    basis, "active_network_zero_movement");
            }
        }

        const auto solver_start = std::chrono::steady_clock::now();
        const HighsStatus run_status = highs.run();
        const double solver_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - solver_start).count();
        output.solver_wall_seconds += solver_seconds;
        trace("round=" + std::to_string(round) +
              ":solver_complete");
        const auto model_status = highs.getModelStatus();
        const auto& solution = highs.getSolution();
        const auto& info = highs.getInfo();
        const bool round_optimal =
            model_status == HighsModelStatus::kOptimal;
        output.all_solver_rounds_optimal =
            output.all_solver_rounds_optimal && round_optimal;
        output.time_limit_reached = output.time_limit_reached ||
            model_status == HighsModelStatus::kTimeLimit;
        double audited_column_violation =
            std::numeric_limits<double>::infinity();
        double audited_row_violation =
            std::numeric_limits<double>::infinity();
        int audited_worst_row = -1;
        double audited_worst_row_activity =
            std::numeric_limits<double>::quiet_NaN();
        double terminal_canonicalization_scale = 1.0;
        std::vector<double> terminal_values;
        const bool terminal_vector_available = solution.value_valid &&
            solution.col_value.size() ==
                static_cast<std::size_t>(column_count);
        if (terminal_vector_available) {
            terminal_values = solution.col_value;
            audited_column_violation = 0.0;
            for (int column = 0; column < column_count; ++column) {
                const double value = solution.col_value[column];
                if (!std::isfinite(value)) {
                    audited_column_violation =
                        std::numeric_limits<double>::infinity();
                    break;
                }
                audited_column_violation = std::max({
                    audited_column_violation,
                    lower[column] - value,
                    value - upper[column],
                });
                terminal_values[column] = std::clamp(
                    value, lower[column], upper[column]);
            }
            if (eliminate_angles &&
                (output.angle_response_row_count > 0 ||
                 output.voltage_response_row_count > 0)) {
                terminal_canonicalization_scale = 0.9999;
                for (int column = 0; column < column_count; ++column) {
                    terminal_values[column] = std::clamp(
                        terminal_canonicalization_scale *
                            terminal_values[column],
                        lower[column], upper[column]);
                }
            }
            audited_row_violation = 0.0;
            for (int row_index = 0;
                 row_index < static_cast<int>(rows.size()); ++row_index) {
                const auto& row = rows[row_index];
                long double activity = 0.0L;
                for (const auto& [column, coefficient] : row.entries) {
                    activity += static_cast<long double>(coefficient) *
                        static_cast<long double>(terminal_values[column]);
                }
                if (!std::isfinite(activity)) {
                    audited_row_violation =
                        std::numeric_limits<double>::infinity();
                    break;
                }
                if (std::isfinite(row.lower)) {
                    const double violation = static_cast<double>(
                        static_cast<long double>(row.lower) - activity);
                    if (violation > audited_row_violation) {
                        audited_row_violation = violation;
                        audited_worst_row = row_index;
                        audited_worst_row_activity =
                            static_cast<double>(activity);
                    }
                }
                if (std::isfinite(row.upper)) {
                    const double violation = static_cast<double>(
                        activity - static_cast<long double>(row.upper));
                    if (violation > audited_row_violation) {
                        audited_row_violation = violation;
                        audited_worst_row = row_index;
                        audited_worst_row_activity =
                            static_cast<double>(activity);
                    }
                }
            }
        }
        constexpr double kProposalColumnCanonicalizationTolerance = 1e-5;
        const bool independently_feasible_terminal_solution =
            terminal_vector_available &&
            audited_column_violation <=
                kProposalColumnCanonicalizationTolerance &&
            audited_row_violation <= 1e-7;
        // HiGHS can leave a primal-simplex basis exactly feasible at a time
        // limit while reporting an unknown primal_solution_status.  The
        // original sparse rows and bounds are authoritative, so preserve that
        // terminal vector only after this independent full audit.
        const bool round_feasible = info.valid &&
            independently_feasible_terminal_solution;
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
            {"angle_response_row_count",
             output.angle_response_row_count},
            {"voltage_response_row_count",
             output.voltage_response_row_count},
            {"reactive_capability_row_count",
             output.reactive_capability_row_count},
            {"sparse_full_ac_linearization", sparse_full_ac},
            {"reduced_pv_pq_partition", reduced_pv_pq},
            {"generate_trust_region_rows",
             options.generate_trust_region_rows},
            {"reduced_pv_bus_count", reduced_pv_bus_count},
            {"objective_scale", objective_scale},
            {"solver", solver},
            {"solver_threads", highs_threads},
            {"primal_start_attempted", primal_start_attempted},
            {"primal_start_status", static_cast<int>(start_status)},
            {"basis_attempted", basis_attempted},
            {"basis_status", static_cast<int>(basis_status)},
            {"run_status", static_cast<int>(run_status)},
            {"model_status", static_cast<int>(model_status)},
            {"status", output.status},
            {"linear_objective", info.objective_function_value},
            {"optimal", round_optimal},
            {"solution_value_valid", solution.value_valid},
            {"info_valid", info.valid},
            {"primal_solution_status", info.primal_solution_status},
            {"max_primal_infeasibility", info.max_primal_infeasibility},
            {"audited_column_violation", audited_column_violation},
            {"audited_row_violation", audited_row_violation},
            {"audited_worst_row", audited_worst_row},
            {"audited_worst_row_activity",
             audited_worst_row_activity},
            {"audited_worst_row_lower",
             audited_worst_row >= 0
                ? rows[audited_worst_row].lower
                : std::numeric_limits<double>::quiet_NaN()},
            {"audited_worst_row_upper",
             audited_worst_row >= 0
                ? rows[audited_worst_row].upper
                : std::numeric_limits<double>::quiet_NaN()},
            {"terminal_canonicalization_scale",
             terminal_canonicalization_scale},
            {"independently_feasible_terminal_solution",
             independently_feasible_terminal_solution},
            {"simplex_iterations", info.simplex_iteration_count},
            {"ipm_iterations", info.ipm_iteration_count},
            {"solver_wall_seconds", solver_seconds},
        });
        if (!round_feasible) {
            break;
        }
        ++output.rounds_completed;

        double reduced_step_scale = 1.0;
        std::vector<double> reduced_angle_change;
        std::vector<double> reduced_voltage_change;
        if (eliminate_angles) {
            std::vector<double> component_slack;
            if (!reduced_network->reconstruct(
                    terminal_values, reduced_angle_change,
                    reduced_voltage_change,
                    component_slack)) {
                output.status =
                    "reduced active-network reconstruction failed";
                break;
            }
            double maximum_component_slack = 0.0;
            for (double value : component_slack) {
                maximum_component_slack = std::max(
                    maximum_component_slack, std::abs(value));
            }
            if (maximum_component_slack > 1e-5) {
                output.status =
                    "reduced active-network component balance failed";
                break;
            }
            double maximum_angle_change = 0.0;
            for (double value : reduced_angle_change) {
                maximum_angle_change = std::max(
                    maximum_angle_change, std::abs(value));
            }
            double maximum_angle_difference_change = 0.0;
            for (const auto& branch : data.branches) {
                if (!branch.present || branch.status == 0) {
                    continue;
                }
                maximum_angle_difference_change = std::max(
                    maximum_angle_difference_change,
                    std::abs(
                        reduced_angle_change[branch.from] -
                        reduced_angle_change[branch.to]));
            }
            double maximum_voltage_change = 0.0;
            for (double value : reduced_voltage_change) {
                maximum_voltage_change = std::max(
                    maximum_voltage_change, std::abs(value));
            }
            output.rounds.back()["maximum_component_slack"] =
                maximum_component_slack;
            output.rounds.back()["maximum_reduced_angle_change"] =
                maximum_angle_change;
            output.rounds.back()[
                "maximum_reduced_angle_difference_change"] =
                maximum_angle_difference_change;
            output.rounds.back()["maximum_reduced_voltage_change"] =
                maximum_voltage_change;

            // Voltage magnitude is eliminated at a reduced-space PV bus, so
            // its aggregate reactive generation is a dependent quantity.
            // Generate only capability rows violated by the current LP
            // direction, then re-solve before attempting nonlinear repair.
            // This is a physical Q-limit constraint, not an optional trust-
            // region heuristic.
            if (reduced_pv_pq) {
                struct ReactiveCapabilityViolation {
                    double severity{};
                    int bus{};
                };
                std::vector<ReactiveCapabilityViolation> violations;
                violations.reserve(static_cast<std::size_t>(nb));
                bool reactive_linearization_failure = false;
                for (int bus = 0; bus < nb; ++bus) {
                    if (reduced_network->pv_bus_mask()[bus] == 0 ||
                        monitored_reactive_bus[bus] != 0) {
                        continue;
                    }
                    double lower = 0.0;
                    double upper = 0.0;
                    for (int generator : data.buses[bus].generators) {
                        if (commitment[generator] == 0) {
                            continue;
                        }
                        lower += data.generators[generator].qmin;
                        upper += data.generators[generator].qmax;
                    }
                    const double required = pv_reactive_requirement(
                        data, reference.state, bus);
                    // Preserve an independently verified incumbent that uses
                    // source-authorized local Q slack outside aggregate
                    // capability, exactly as the generated row does.
                    lower = std::min(lower, required);
                    upper = std::max(upper, required);
                    const double change =
                        linearized_pv_reactive_requirement_change(
                            data, reference.state, flow_derivative,
                            reduced_angle_change, reduced_voltage_change,
                            terminal_values, bus);
                    if (!std::isfinite(required) ||
                        !std::isfinite(change)) {
                        output.status =
                            "reduced active-network reactive capability screening failed";
                        reactive_linearization_failure = true;
                        break;
                    }
                    const double predicted = required + change;
                    const double violation = std::max(
                        lower - predicted, predicted - upper);
                    if (violation > 1e-8) {
                        const double scale = std::max(
                            {1.0, std::abs(lower), std::abs(upper)});
                        violations.push_back({violation / scale, bus});
                    }
                }
                if (reactive_linearization_failure) {
                    break;
                }
                std::sort(
                    violations.begin(), violations.end(),
                    [](const auto& left, const auto& right) {
                        if (left.severity != right.severity) {
                            return left.severity > right.severity;
                        }
                        return left.bus < right.bus;
                    });
                int new_reactive_rows = 0;
                for (const auto& violation : violations) {
                    if (new_reactive_rows >=
                        options.maximum_trust_region_rows_per_pass) {
                        break;
                    }
                    monitored_reactive_bus[violation.bus] = 1;
                    ++new_reactive_rows;
                }
                output.rounds.back()[
                    "violated_reactive_capability_rows"] =
                    violations.size();
                output.rounds.back()["new_reactive_capability_rows"] =
                    new_reactive_rows;
                if (new_reactive_rows > 0) {
                    ++output
                        .reactive_capability_constraint_generation_passes;
                    trace("round=" + std::to_string(round) +
                          ":reactive_capability_rows_added=" +
                          std::to_string(new_reactive_rows) +
                          ":violated=" +
                          std::to_string(violations.size()));
                    continue;
                }
            }

            if (options.generate_trust_region_rows) {
            struct StateResponseViolation {
                double severity{};
                int type{};  // 0 = branch-angle difference, 1 = voltage
                int bus{};   // branch index for type 0, bus index for type 1
            };
            std::vector<StateResponseViolation> violations;
            violations.reserve(static_cast<std::size_t>(nb));
            for (int branch_index = 0;
                 branch_index < static_cast<int>(data.branches.size());
                 ++branch_index) {
                const auto& branch = data.branches[branch_index];
                if (!branch.present || branch.status == 0 ||
                    monitored_angle_branch[branch_index] != 0) {
                    continue;
                }
                const double difference_change =
                    reduced_angle_change[branch.from] -
                    reduced_angle_change[branch.to];
                if (std::abs(difference_change) >
                    options.angle_trust_radius + 1e-8) {
                    violations.push_back({
                        std::abs(difference_change) /
                                options.angle_trust_radius -
                            1.0,
                        0,
                        branch_index,
                    });
                }
            }
            for (int bus = 0; bus < nb; ++bus) {
                const double lower_change =
                    std::max(
                        -options.voltage_trust_radius,
                        data.buses[bus].vmin - reference.state.vm[bus]);
                const double upper_change =
                    std::min(
                        options.voltage_trust_radius,
                        data.buses[bus].vmax - reference.state.vm[bus]);
                const double violation = std::max(
                    lower_change - reduced_voltage_change[bus],
                    reduced_voltage_change[bus] - upper_change);
                if (violation > 1e-8 &&
                    monitored_voltage_bus[bus] == 0) {
                    violations.push_back({
                        violation / options.voltage_trust_radius,
                        1,
                        bus,
                    });
                }
            }
            std::sort(
                violations.begin(), violations.end(),
                [](const auto& left, const auto& right) {
                    if (left.severity != right.severity) {
                        return left.severity > right.severity;
                    }
                    if (left.type != right.type) {
                        return left.type < right.type;
                    }
                    return left.bus < right.bus;
                });
            int new_angle_rows = 0;
            int new_voltage_rows = 0;
            int new_state_rows = 0;
            for (const auto& violation : violations) {
                if (new_state_rows >=
                    options.maximum_trust_region_rows_per_pass) {
                    break;
                }
                if (violation.type == 0) {
                    monitored_angle_branch[violation.bus] = 1;
                    ++new_angle_rows;
                    ++new_state_rows;
                } else if (new_voltage_rows <
                    options.maximum_voltage_rows_per_pass) {
                    monitored_voltage_bus[violation.bus] = 1;
                    ++new_voltage_rows;
                    ++new_state_rows;
                }
            }
            if (new_state_rows > 0) {
                ++output.trust_region_constraint_generation_passes;
                if (new_voltage_rows > 0) {
                    ++output.voltage_constraint_generation_passes;
                }
                output.rounds.back()["new_angle_response_rows"] =
                    new_angle_rows;
                output.rounds.back()["new_voltage_response_rows"] =
                    new_voltage_rows;
                output.rounds.back()["violated_state_response_rows"] =
                    violations.size();
                trace("round=" + std::to_string(round) +
                      ":state_rows_added=" +
                      std::to_string(new_state_rows) +
                      ":violated=" +
                      std::to_string(violations.size()));
                continue;
            }
            }
            if (maximum_angle_difference_change >
                options.angle_trust_radius) {
                reduced_step_scale =
                    options.angle_trust_radius /
                    maximum_angle_difference_change;
            }
            if (maximum_voltage_change > options.voltage_trust_radius) {
                reduced_step_scale = std::min(
                    reduced_step_scale,
                    options.voltage_trust_radius /
                        maximum_voltage_change);
            }
            output.rounds.back()["trust_region_step_scale"] =
                reduced_step_scale;
            for (int bus = 0; bus < nb; ++bus) {
                const double change = reduced_voltage_change[bus];
                if (change > 1e-14 &&
                    reduced_step_scale * change >
                        data.buses[bus].vmax -
                            reference.state.vm[bus]) {
                    reduced_step_scale = std::min(
                        reduced_step_scale,
                        std::max(
                            0.0, data.buses[bus].vmax -
                                reference.state.vm[bus]) / change);
                } else if (change < -1e-14 &&
                    reduced_step_scale * change <
                        data.buses[bus].vmin -
                            reference.state.vm[bus]) {
                    reduced_step_scale = std::min(
                        reduced_step_scale,
                        std::max(
                            0.0, reference.state.vm[bus] -
                                data.buses[bus].vmin) / -change);
                }
            }
            output.rounds.back()["selected_step_scale"] =
                reduced_step_scale;
        }

        AcState target = reference.state;
        for (int generator = 0; generator < ng; ++generator) {
            target.pg[generator] += reduced_step_scale *
                (compact_signed_columns
                ? terminal_values[pg_signed_offset + generator]
                : terminal_values[pg_up_offset + generator] -
                    terminal_values[pg_down_offset + generator]);
            if ((eliminate_angles && !reduced_pv_pq) ||
                sparse_full_ac) {
                target.qg[generator] += reduced_step_scale *
                    terminal_values[qg_signed_offset + generator];
            }
        }
        for (int load = 0; load < nd; ++load) {
            target.demand_factor[load] += reduced_step_scale *
                (compact_signed_columns
                ? terminal_values[demand_signed_offset + load]
                : terminal_values[demand_up_offset + load] -
                    terminal_values[demand_down_offset + load]);
        }
        if (eliminate_angles) {
            for (int bus = 0; bus < nb; ++bus) {
                target.va[bus] +=
                    reduced_step_scale * reduced_angle_change[bus];
                target.vm[bus] +=
                    reduced_step_scale * reduced_voltage_change[bus];
            }
        } else {
            for (int bus = 0; bus < nb; ++bus) {
                target.va[bus] += compact_signed_columns
                    ? terminal_values[angle_signed_offset + bus]
                    : terminal_values[angle_up_offset + bus] -
                        terminal_values[angle_down_offset + bus];
                if (sparse_full_ac) {
                    target.vm[bus] +=
                        terminal_values[voltage_signed_offset + bus];
                }
            }
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
                proposal.state.qg[generator] += fraction *
                    (target.qg[generator] - reference.state.qg[generator]);
            }
            for (int load = 0; load < nd; ++load) {
                proposal.state.demand_factor[load] += fraction *
                    (target.demand_factor[load] -
                     reference.state.demand_factor[load]);
            }
            for (int bus = 0; bus < nb; ++bus) {
                proposal.state.va[bus] += fraction *
                    (target.va[bus] - reference.state.va[bus]);
                proposal.state.vm[bus] += fraction *
                    (target.vm[bus] - reference.state.vm[bus]);
            }
            // The first rebuild supplies exact nonlinear branch flows. Then
            // recover the dependent PV-bus reactive outputs and rebuild the
            // objective/slacks against that physically consistent Q state.
            rebuild_base_state_derived_fields(
                data, commitment, proposal.state);
            ReactiveGenerationRecovery reactive_recovery;
            if (reduced_pv_pq) {
                reactive_recovery = recover_pv_bus_reactive_generation(
                    data, commitment, proposal.state,
                    &reduced_network->pv_bus_mask());
            } else if (!eliminate_angles && !sparse_full_ac) {
                reactive_recovery = recover_pv_bus_reactive_generation(
                    data, commitment, proposal.state);
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
            const double raw_thermal_slack = std::accumulate(
                proposal.state.sm_slack.begin(),
                proposal.state.sm_slack.end(), 0.0);
            const double raw_economic_objective = raw_objective +
                data.delta * (
                    data.p_delta_cost_approx * active_slack +
                    data.q_delta_cost_approx * reactive_slack +
                    data.sm_cost_approx * raw_thermal_slack);
            const auto raw_pg = proposal.state.pg;
            const auto raw_demand_factor =
                proposal.state.demand_factor;
            bool repair_attempted = false;
            bool repair_skipped_balance_slack_guard = false;
            double repair_seconds = 0.0;
            int repair_newton_iterations = 0;
            int repair_active_redispatch_passes = 0;
            int repair_reactive_limit_passes = 0;
            double repair_initial_newton_residual = 0.0;
            bool repair_newton_candidate_selected = false;
            std::string repair_failure_reason;
            if (active_slack + reactive_slack > 1e-8 ||
                !validated_candidate_is_feasible(
                    proposal, proposal_validation,
                    options.validation_tolerance)) {
                repair_skipped_balance_slack_guard =
                    active_slack + reactive_slack >
                    options.maximum_candidate_repair_balance_slack;
                if (!repair_skipped_balance_slack_guard) {
                    repair_attempted = true;
                    trace("round=" + std::to_string(round) +
                          ":trial=" + std::to_string(trial + 1) +
                          ":repair_start");
                    FastPowerFlowOptions fast_options;
                    fast_options.distributed_balance_polish = false;
                    fast_options.minimize_active_balance_slack = true;
                    fast_options.minimize_reactive_balance_slack = true;
                    fast_options.skip_balance_cleanup_prepasses =
                        eliminate_angles || sparse_full_ac;
                    fast_options.max_newton_iterations =
                        options.maximum_candidate_repair_newton_iterations;
                    fast_options.max_active_redispatch_passes = 6;
                    fast_options.max_reactive_limit_passes = 4;
                    FastContingencyPowerFlow repair(
                        data, proposal.state, commitment, fast_options);
                    auto repaired = repair.solve_base();
                    repair_seconds = repaired.wall_seconds;
                    repair_newton_iterations = repaired.newton_iterations;
                    repair_active_redispatch_passes =
                        repaired.active_redispatch_passes;
                    repair_reactive_limit_passes =
                        repaired.reactive_limit_passes;
                    repair_initial_newton_residual =
                        repaired.initial_newton_residual;
                    repair_newton_candidate_selected =
                        repaired.newton_candidate_selected;
                    repair_failure_reason = repaired.failure_reason;
                    repaired.solve.objective =
                        rebuild_base_state_derived_fields(
                            data, commitment, repaired.solve.state);
                    repaired.validation = validate_state(
                        data, ModelMode::BaseSoft,
                        repaired.solve.state, commitment);
                    proposal = std::move(repaired.solve);
                    proposal_validation = repaired.validation;
                    trace("round=" + std::to_string(round) +
                          ":trial=" + std::to_string(trial + 1) +
                          ":repair_complete");
                }
            }
            const double candidate_active_slack = std::accumulate(
                proposal.state.p_delta.begin(),
                proposal.state.p_delta.end(), 0.0);
            const double candidate_reactive_slack = std::accumulate(
                proposal.state.q_delta.begin(),
                proposal.state.q_delta.end(), 0.0);
            const double candidate_thermal_slack = std::accumulate(
                proposal.state.sm_slack.begin(),
                proposal.state.sm_slack.end(), 0.0);
            const double candidate_economic_objective =
                proposal.objective + data.delta * (
                    data.p_delta_cost_approx * candidate_active_slack +
                    data.q_delta_cost_approx * candidate_reactive_slack +
                    data.sm_cost_approx * candidate_thermal_slack);
            double repair_generation_movement = 0.0;
            double repair_maximum_generation_movement = 0.0;
            for (int generator = 0; generator < ng; ++generator) {
                const double movement = std::abs(
                    proposal.state.pg[generator] - raw_pg[generator]);
                repair_generation_movement += movement;
                repair_maximum_generation_movement = std::max(
                    repair_maximum_generation_movement, movement);
            }
            double repair_load_factor_movement = 0.0;
            for (int load = 0; load < nd; ++load) {
                repair_load_factor_movement += std::abs(
                    proposal.state.demand_factor[load] -
                    raw_demand_factor[load]);
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
                {"raw_thermal_slack", raw_thermal_slack},
                {"raw_economic_objective", raw_economic_objective},
                {"pv_reactive_bus_count",
                 reactive_recovery.pv_bus_count},
                {"pv_reactive_saturated_bus_count",
                 reactive_recovery.saturated_bus_count},
                {"pv_reactive_maximum_unserved_requirement",
                 reactive_recovery.maximum_unserved_requirement},
                {"raw_validation", raw_validation.to_json()},
                {"repair_attempted", repair_attempted},
                {"repair_skipped_balance_slack_guard",
                 repair_skipped_balance_slack_guard},
                {"repair_seconds", repair_seconds},
                {"repair_newton_iterations", repair_newton_iterations},
                {"repair_active_redispatch_passes",
                 repair_active_redispatch_passes},
                {"repair_reactive_limit_passes",
                 repair_reactive_limit_passes},
                {"repair_initial_newton_residual",
                 repair_initial_newton_residual},
                {"repair_newton_candidate_selected",
                 repair_newton_candidate_selected},
                {"repair_failure_reason", repair_failure_reason},
                {"repair_generation_movement",
                 repair_generation_movement},
                {"repair_maximum_generation_movement",
                 repair_maximum_generation_movement},
                {"repair_load_factor_movement",
                 repair_load_factor_movement},
                {"candidate_objective", proposal.objective},
                {"candidate_active_balance_slack",
                 candidate_active_slack},
                {"candidate_reactive_balance_slack",
                 candidate_reactive_slack},
                {"candidate_thermal_slack", candidate_thermal_slack},
                {"candidate_economic_objective",
                 candidate_economic_objective},
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
