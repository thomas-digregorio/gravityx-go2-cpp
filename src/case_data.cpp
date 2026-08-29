#include "gravityx/case_data.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <utility>

namespace gravityx {
namespace {

using json = nlohmann::json;

template <typename T>
T value_or(const json& object, const char* key, T fallback) {
    const auto it = object.find(key);
    return it == object.end() || it->is_null() ? fallback : it->get<T>();
}

template <typename Entry, typename Builder>
std::vector<Entry> sorted_entries(const json& object, Builder&& builder) {
    std::vector<Entry> result;
    result.reserve(object.size());
    for (auto it = object.begin(); it != object.end(); ++it) {
        result.push_back(builder(it.key(), it.value()));
    }
    std::sort(result.begin(), result.end(), [](const Entry& left, const Entry& right) {
        return left.index < right.index;
    });
    return result;
}

std::vector<double> doubles(const json& array) {
    std::vector<double> values;
    values.reserve(array.size());
    for (const auto& value : array) {
        values.push_back(value.get<double>());
    }
    return values;
}

std::string source_identifier(const json& item, std::size_t position) {
    const auto& value = item.at("source_id").at(position);
    if (value.is_string()) {
        return value.get<std::string>();
    }
    if (value.is_number_integer()) {
        return std::to_string(value.get<long long>());
    }
    if (value.is_number_unsigned()) {
        return std::to_string(value.get<unsigned long long>());
    }
    throw std::runtime_error("source identifier is not a string or integer");
}

int source_integer(const json& item, std::size_t position) {
    const auto& value = item.at("source_id").at(position);
    if (value.is_number_integer() || value.is_number_unsigned()) {
        return value.get<int>();
    }
    if (value.is_string()) {
        return std::stoi(value.get<std::string>());
    }
    throw std::runtime_error("source identifier is not an integer");
}

}  // namespace

void refresh_branch_flow_coefficients(Branch& branch) {
    const double denominator = branch.r * branch.r + branch.x * branch.x;
    const double g = denominator > 1e-20
        ? branch.r / denominator : 0.0;
    const double b = denominator > 1e-20
        ? -branch.x / denominator : 0.0;
    if (!std::isfinite(branch.tap) || std::abs(branch.tap) <= 1e-12) {
        throw std::runtime_error(
            "branch has an invalid tap ratio: " + branch.source_key);
    }
    const double tap_squared = branch.tap * branch.tap;
    const double tap_real = branch.tap * std::cos(branch.shift);
    const double tap_imag = branch.tap * std::sin(branch.shift);
    branch.flow_from_g_self = branch.transformer
        ? g / tap_squared + branch.g_fr
        : (g + branch.g_fr) / tap_squared;
    branch.flow_from_b_self = branch.transformer
        ? b / tap_squared + branch.b_fr
        : (b + branch.b_fr) / tap_squared;
    branch.flow_to_g_self = g + branch.g_to;
    branch.flow_to_b_self = b + branch.b_to;
    branch.flow_from_cross_cos =
        (-g * tap_real + b * tap_imag) / tap_squared;
    branch.flow_from_cross_sin =
        (-b * tap_real - g * tap_imag) / tap_squared;
    branch.flow_to_cross_cos =
        (-g * tap_real - b * tap_imag) / tap_squared;
    branch.flow_to_cross_sin =
        (-b * tap_real + g * tap_imag) / tap_squared;
    branch.flow_coefficients_valid = true;
}

std::pair<int, int> transformer_step_bounds(const Branch& branch) {
    if (!branch.transformer) {
        throw std::runtime_error(
            "transformer step requested for a non-transformer branch: " +
            branch.source_key);
    }
    const int mode = std::abs(branch.control_mode);
    const int steps = mode == 1 ? branch.tm_steps
        : mode == 3 ? branch.ta_steps : 0;
    if (steps < 3 || steps % 2 == 0) {
        throw std::runtime_error(
            "transformer control positions must be an odd source count: " +
            branch.source_key);
    }
    const int maximum_step = (steps - 1) / 2;
    return {-maximum_step, maximum_step};
}

double transformer_step_value(const Branch& branch, int step) {
    const auto [lower_step, upper_step] = transformer_step_bounds(branch);
    if (step < lower_step || step > upper_step) {
        throw std::runtime_error(
            "transformer step is outside the source range: " +
            branch.source_key);
    }
    const int mode = std::abs(branch.control_mode);
    const int steps = mode == 1 ? branch.tm_steps : branch.ta_steps;
    const double minimum = mode == 1 ? branch.tm_min : branch.ta_min;
    const double maximum = mode == 1 ? branch.tm_max : branch.ta_max;
    if (!std::isfinite(minimum) || !std::isfinite(maximum) ||
        maximum < minimum) {
        throw std::runtime_error(
            "transformer has an invalid source control range: " +
            branch.source_key);
    }
    const double increment = (maximum - minimum) /
        static_cast<double>(steps - 1);
    const double midpoint = 0.5 * (minimum + maximum);
    return midpoint + increment * static_cast<double>(step);
}

void set_transformer_step(Branch& branch, int step) {
    const double value = transformer_step_value(branch, step);
    const int mode = std::abs(branch.control_mode);
    if (mode == 1) {
        branch.tap = value;
        branch.tm_step = step;
    } else {
        branch.shift = value;
        branch.ta_step = step;
    }
    refresh_branch_flow_coefficients(branch);
}

CaseData CaseData::load(const std::string& path) {
    std::ifstream stream(path);
    if (!stream) {
        throw std::runtime_error("cannot open normalized case JSON: " + path);
    }
    json root;
    stream >> root;

    CaseData data;
    data.name = root.at("name").get<std::string>();
    data.base_mva = root.at("baseMVA").get<double>();
    data.delta = root.at("delta").get<double>();
    data.delta_r = root.at("deltar").get<double>();
    data.delta_ctg = root.at("deltactg").get<double>();
    data.delta_r_ctg = root.at("deltarctg").get<double>();
    data.sm_vio_limit = root.at("sm_vio_limit").get<double>();
    data.sm_cost_approx = root.at("sm_cost_approx").get<double>();
    data.p_delta_cost_approx = root.at("p_delta_cost_approx").get<double>();
    data.q_delta_cost_approx = root.at("q_delta_cost_approx").get<double>();

    data.buses = sorted_entries<Bus>(root.at("bus"), [](const std::string& key, const json& item) {
        Bus bus;
        bus.source_key = key;
        bus.index = item.at("index").get<int>();
        bus.bus_i = item.at("bus_i").get<int>();
        bus.type = item.at("bus_type").get<int>();
        bus.vmin = item.at("vmin").get<double>();
        bus.vmax = item.at("vmax").get<double>();
        bus.vm_start = value_or(item, "vm_start", value_or(item, "vm", 1.0));
        bus.va_start = value_or(item, "va_start", value_or(item, "va", 0.0));
        bus.present = item.value("present", true);
        return bus;
    });
    for (int position = 0; position < static_cast<int>(data.buses.size()); ++position) {
        data.bus_position.emplace(data.buses[position].bus_i, position);
    }

    const auto bus_pos = [&data](int bus_i) {
        const auto it = data.bus_position.find(bus_i);
        if (it == data.bus_position.end()) {
            throw std::runtime_error("component references missing bus " + std::to_string(bus_i));
        }
        return it->second;
    };

    data.generators = sorted_entries<Generator>(root.at("gen"), [&bus_pos](const std::string& key, const json& item) {
        Generator gen;
        gen.source_key = key;
        gen.index = item.at("index").get<int>();
        gen.bus = bus_pos(item.at("gen_bus").get<int>());
        gen.status_prev = item.at("status_prev").get<int>();
        gen.suqual = item.at("suqual").get<int>();
        gen.sdqual = item.at("sdqual").get<int>();
        gen.pg_start = value_or(item, "pg_start", value_or(item, "pg", 0.0));
        gen.qg_start = value_or(item, "qg_start", value_or(item, "qg", 0.0));
        gen.pg_prev = item.at("pg_prev").get<double>();
        gen.pmin = item.at("pmin").get<double>();
        gen.pmax = item.at("pmax").get<double>();
        gen.qmin = item.at("qmin").get<double>();
        gen.qmax = item.at("qmax").get<double>();
        gen.prumax = item.at("prumax").get<double>();
        gen.prdmax = item.at("prdmax").get<double>();
        gen.prumaxctg = item.at("prumaxctg").get<double>();
        gen.prdmaxctg = item.at("prdmaxctg").get<double>();
        gen.oncost = item.at("oncost").get<double>();
        gen.sucost = item.at("sucost").get<double>();
        gen.sdcost = item.at("sdcost").get<double>();
        gen.ncost = item.at("ncost").get<int>();
        gen.cost = doubles(item.at("cost"));
        gen.present = item.value("present", true);
        gen.source_bus = source_integer(item, 1);
        gen.source_id = source_identifier(item, 2);
        return gen;
    });

    data.loads = sorted_entries<Load>(root.at("load"), [&bus_pos](const std::string& key, const json& item) {
        Load load;
        load.source_key = key;
        load.index = item.at("index").get<int>();
        load.bus = bus_pos(item.at("load_bus").get<int>());
        load.pd_nominal = item.at("pd_nominal").get<double>();
        load.qd_nominal = item.at("qd_nominal").get<double>();
        load.pd_prev = item.at("pd_prev").get<double>();
        load.qd_prev = item.at("qd_prev").get<double>();
        load.pd_min = item.at("pd_min").get<double>();
        load.pd_max = item.at("pd_max").get<double>();
        load.tmin = item.at("tmin").get<double>();
        load.tmax = item.at("tmax").get<double>();
        load.prumax = item.at("prumax").get<double>();
        load.prdmax = item.at("prdmax").get<double>();
        load.prumaxctg = item.at("prumaxctg").get<double>();
        load.prdmaxctg = item.at("prdmaxctg").get<double>();
        const double pd = value_or(item, "pd", load.pd_prev);
        load.z_start = std::abs(load.pd_nominal) > 1e-12 ? pd / load.pd_nominal : 1.0;
        load.ncost = item.at("ncost").get<int>();
        load.cost = doubles(item.at("cost"));
        load.present = item.value("present", true);
        load.source_bus = source_integer(item, 1);
        load.source_id = source_identifier(item, 2);
        return load;
    });

    data.shunts = sorted_entries<Shunt>(root.at("shunt"), [&bus_pos](const std::string& key, const json& item) {
        Shunt shunt;
        shunt.source_key = key;
        shunt.index = item.at("index").get<int>();
        shunt.bus = bus_pos(item.at("shunt_bus").get<int>());
        shunt.gs = item.at("gs").get<double>();
        shunt.bs = item.at("bs").get<double>();
        shunt.dispatchable = item.value("dispatchable", false);
        if (item.contains("xst")) {
            shunt.steps = item.at("xst").get<std::vector<int>>();
        }
        if (item.contains("blocks")) {
            for (const auto& block : item.at("blocks")) {
                if (!block.is_array() || block.size() != 2) {
                    throw std::runtime_error(
                        "invalid switched-shunt block: " +
                        shunt.source_key);
                }
                shunt.block_maximum_steps.push_back(
                    block.at(0).get<int>());
                shunt.block_susceptance.push_back(
                    block.at(1).get<double>());
            }
        }
        if (shunt.steps.size() < shunt.block_maximum_steps.size()) {
            shunt.steps.resize(shunt.block_maximum_steps.size(), 0);
        }
        shunt.present = item.value("present", true);
        shunt.source_bus = source_integer(item, 1);
        return shunt;
    });

    data.branches = sorted_entries<Branch>(root.at("branch"), [&bus_pos](const std::string& key, const json& item) {
        Branch branch;
        branch.source_key = key;
        branch.index = item.at("index").get<int>();
        branch.status = item.value(
            "br_status", item.value("status_prev", 1));
        if (branch.status != 0 && branch.status != 1) {
            throw std::runtime_error(
                "branch status is not binary: " + branch.source_key);
        }
        branch.from = bus_pos(item.at("f_bus").get<int>());
        branch.to = bus_pos(item.at("t_bus").get<int>());
        branch.transformer = item.at("transformer").get<bool>();
        branch.r = item.at("br_r").get<double>();
        branch.x = item.at("br_x").get<double>();
        branch.g_fr = item.at("g_fr").get<double>();
        branch.b_fr = item.at("b_fr").get<double>();
        branch.g_to = item.at("g_to").get<double>();
        branch.b_to = item.at("b_to").get<double>();
        branch.tap = item.at("tap").get<double>();
        branch.shift = item.at("shift").get<double>();
        branch.angmin = item.at("angmin").get<double>();
        branch.angmax = item.at("angmax").get<double>();
        branch.rate_a = item.at("rate_a").get<double>();
        branch.rate_b = item.value("rate_b", branch.rate_a);
        branch.rate_c = item.value("rate_c", branch.rate_a);
        branch.present = item.value("present", true);
        branch.source_from = source_integer(item, 1);
        branch.source_to = source_integer(item, 2);
        branch.source_id = source_identifier(
            item, branch.transformer ? 4 : 3);
        branch.control_mode = item.value("control_mode", 0);
        branch.tm_min = item.value("tm_min", branch.tap);
        branch.tm_max = item.value("tm_max", branch.tap);
        branch.tm_steps = item.value("tm_steps", 1);
        branch.tm_step = item.value("tm_step", 0);
        branch.ta_min = item.value("ta_min", branch.shift);
        branch.ta_max = item.value("ta_max", branch.shift);
        branch.ta_steps = item.value("ta_steps", 1);
        branch.ta_step = item.value("ta_step", 0);
        refresh_branch_flow_coefficients(branch);
        return branch;
    });

    std::unordered_map<int, int> generator_position;
    std::unordered_map<int, int> branch_position;
    for (int i = 0; i < static_cast<int>(data.generators.size()); ++i) {
        generator_position.emplace(data.generators[i].index, i);
    }
    for (int i = 0; i < static_cast<int>(data.branches.size()); ++i) {
        branch_position.emplace(data.branches[i].index, i);
    }
    const auto append_contingencies = [&data](
        const json& array,
        ContingencyType type,
        const std::unordered_map<int, int>& positions) {
        for (const auto& item : array) {
            const int source_index = item.at("idx").get<int>();
            const auto position = positions.find(source_index);
            if (position == positions.end()) {
                throw std::runtime_error("contingency references a missing component index " +
                                         std::to_string(source_index));
            }
            data.contingencies.push_back({
                item.at("label").get<std::string>(), type, source_index, position->second});
        }
    };
    append_contingencies(root.at("gen_contingencies"), ContingencyType::Generator, generator_position);
    append_contingencies(root.at("branch_contingencies"), ContingencyType::Branch, branch_position);
    std::sort(data.contingencies.begin(), data.contingencies.end(),
              [](const Contingency& left, const Contingency& right) {
                  return left.label < right.label;
              });

    for (int i = 0; i < static_cast<int>(data.generators.size()); ++i) {
        data.buses[data.generators[i].bus].generators.push_back(i);
    }
    for (int i = 0; i < static_cast<int>(data.loads.size()); ++i) {
        data.buses[data.loads[i].bus].loads.push_back(i);
    }
    for (int i = 0; i < static_cast<int>(data.shunts.size()); ++i) {
        data.buses[data.shunts[i].bus].shunts.push_back(i);
    }
    for (int i = 0; i < static_cast<int>(data.branches.size()); ++i) {
        data.buses[data.branches[i].from].branches_from.push_back(i);
        data.buses[data.branches[i].to].branches_to.push_back(i);
    }

    return data;
}

std::vector<PwlPoint> active_pwl_points(
    const std::vector<double>& flat_cost,
    int ncost,
    double pmin,
    double pmax,
    double tolerance) {
    if (ncost < 2 || static_cast<int>(flat_cost.size()) != 2 * ncost) {
        throw std::runtime_error("piecewise-linear curve must contain at least two complete points");
    }
    if (!(pmin <= pmax) || !std::isfinite(pmin) || !std::isfinite(pmax)) {
        throw std::runtime_error("piecewise-linear curve requires finite ordered bounds");
    }

    std::vector<PwlPoint> points;
    points.reserve(ncost);
    for (int i = 0; i < ncost; ++i) {
        points.push_back({flat_cost[2 * i], flat_cost[2 * i + 1]});
    }

    std::size_t first = 0;
    for (std::size_t i = 0; i + 1 < points.size(); ++i) {
        first = i;
        if (pmin <= points[i + 1].mw) {
            break;
        }
    }

    std::size_t last = points.size() - 1;
    for (std::size_t distance = 1; distance < points.size(); ++distance) {
        const std::size_t preceding = points.size() - 1 - distance;
        last = points.size() - distance;
        if (pmax >= points[preceding].mw) {
            break;
        }
    }
    if (first >= last) {
        throw std::runtime_error("piecewise-linear active range has fewer than two points");
    }

    std::vector<PwlPoint> active(points.begin() + static_cast<std::ptrdiff_t>(first),
                                 points.begin() + static_cast<std::ptrdiff_t>(last + 1));

    if (active.front().mw > pmin) {
        const double x0 = pmin - tolerance;
        const double dx = active[1].mw - active[0].mw;
        const double slope = std::abs(dx) > 0.0 ? (active[1].cost - active[0].cost) / dx : NAN;
        const double y0 = std::isnan(slope)
            ? active.front().cost
            : active[1].cost - slope * (active[1].mw - x0);
        active.front() = {x0, y0};
    }

    if (active.back().mw < pmax) {
        const std::size_t n = active.size();
        const double x3 = pmax + tolerance;
        const double dx = active[n - 1].mw - active[n - 2].mw;
        const double slope = std::abs(dx) > 0.0 ? (active[n - 1].cost - active[n - 2].cost) / dx : NAN;
        const double y3 = std::isnan(slope)
            ? active.back().cost
            : slope * (x3 - active[n - 2].mw) + active[n - 2].cost;
        active.back() = {x3, y3};
    }

    return active;
}

double branch_terminal_component_bound(
    const CaseData& data,
    const Branch& branch,
    double rating,
    bool from_terminal) {
    if (!std::isfinite(rating) || rating < 0.0) {
        throw std::runtime_error(
            "branch terminal component bound requires a finite nonnegative rating");
    }
    const double voltage_upper = branch.transformer
        ? 1.0
        : data.buses[from_terminal ? branch.from : branch.to].vmax;
    return rating * (voltage_upper + data.sm_vio_limit);
}

}  // namespace gravityx
