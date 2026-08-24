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

}  // namespace

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
        return load;
    });

    data.shunts = sorted_entries<Shunt>(root.at("shunt"), [&bus_pos](const std::string& key, const json& item) {
        Shunt shunt;
        shunt.source_key = key;
        shunt.index = item.at("index").get<int>();
        shunt.bus = bus_pos(item.at("shunt_bus").get<int>());
        shunt.gs = item.at("gs").get<double>();
        shunt.bs = item.at("bs").get<double>();
        return shunt;
    });

    data.branches = sorted_entries<Branch>(root.at("branch"), [&bus_pos](const std::string& key, const json& item) {
        Branch branch;
        branch.source_key = key;
        branch.index = item.at("index").get<int>();
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

}  // namespace gravityx
