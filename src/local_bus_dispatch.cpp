#include "gravityx/local_bus_dispatch.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace gravityx {
namespace {

constexpr double balance_limit = 0.5 - 1e-7;

double curve_value(const std::vector<PwlPoint>& curve, double x) {
    if (x <= curve.front().mw) return curve.front().cost;
    if (x >= curve.back().mw) return curve.back().cost;
    for (std::size_t i = 1; i < curve.size(); ++i) {
        if (x <= curve[i].mw) {
            const double width = curve[i].mw - curve[i - 1].mw;
            if (width <= 1e-14) return curve[i].cost;
            const double t = (x - curve[i - 1].mw) / width;
            return (1.0 - t) * curve[i - 1].cost + t * curve[i].cost;
        }
    }
    throw std::runtime_error("local dispatch PWL lookup failed");
}

bool clip_balance_interval(double current, double mismatch, double slope,
                           double& lower, double& upper) {
    if (std::abs(slope) <= 1e-15) return std::abs(mismatch) <= balance_limit;
    const double a = current + (-balance_limit - mismatch) / slope;
    const double b = current + ( balance_limit - mismatch) / slope;
    lower = std::max(lower, std::min(a, b));
    upper = std::min(upper, std::max(a, b));
    return lower <= upper;
}

struct CoordinateResult { double value{}, gain{}; };

// Exhaustive kink enumeration for this ONE-dimensional piecewise-linear
// subproblem, without a convexity assumption or convexifying source curves.
CoordinateResult coordinate_minimum(
    const std::vector<PwlPoint>& curve, double power_scale, double cost_sign,
    double current, double lower, double upper, double rp, double rq,
    double p_slope, double q_slope, double kp, double kq) {
    CoordinateResult result{current, 0.0};
    if (!clip_balance_interval(current, rp, p_slope, lower, upper) ||
        !clip_balance_interval(current, rq, q_slope, lower, upper)) return result;
    const auto cost = [&](double value) {
        const double change = value - current;
        return cost_sign * curve_value(curve, power_scale * value) +
            kp * std::abs(rp + p_slope * change) +
            kq * std::abs(rq + q_slope * change);
    };
    const double original = cost(current);
    double best = original;
    const auto consider = [&](double x) {
        if (!std::isfinite(x) || x < lower || x > upper) return;
        const double value = cost(x);
        if (std::isfinite(value) && value < best - 1e-9) {
            best = value;
            result.value = x;
        }
    };
    consider(lower); consider(upper);
    if (std::abs(p_slope) > 1e-15) consider(current - rp / p_slope);
    if (std::abs(q_slope) > 1e-15) consider(current - rq / q_slope);
    if (std::abs(power_scale) > 1e-15) {
        for (const auto& point : curve) consider(point.mw / power_scale);
    }
    result.gain = original - best;
    return result;
}

void finite_vector(const std::vector<double>& values) {
    for (double value : values) if (!std::isfinite(value))
        throw std::runtime_error("nonfinite local dispatch input");
}

}  // namespace

LocalBusDispatchResult improve_fixed_network_bus_dispatch(
    const CaseData& data, const AcState& original_base,
    const std::vector<int>& commitment, const Contingency& contingency,
    const std::vector<double>& network_p, const std::vector<double>& network_q,
    AcState& candidate, int passes) {
    const auto ng = data.generators.size(), nd = data.loads.size(), nb = data.buses.size();
    if (passes < 1 || passes > 2 || commitment.size() != ng ||
        original_base.pg.size() != ng || original_base.demand_factor.size() != nd ||
        candidate.pg.size() != ng || candidate.qg.size() != ng ||
        candidate.demand_factor.size() != nd || network_p.size() != nb || network_q.size() != nb ||
        !std::isfinite(data.delta_ctg) || data.delta_ctg <= 0.0 ||
        !std::isfinite(data.p_delta_cost_approx) || data.p_delta_cost_approx < 0.0 ||
        !std::isfinite(data.q_delta_cost_approx) || data.q_delta_cost_approx < 0.0) {
        throw std::runtime_error("invalid local dispatch dimensions or policy");
    }
    finite_vector(original_base.pg); finite_vector(original_base.demand_factor);
    finite_vector(candidate.pg); finite_vector(candidate.qg); finite_vector(candidate.demand_factor);
    finite_vector(network_p); finite_vector(network_q);
    const int outaged_generator = contingency.type == ContingencyType::Generator
        ? contingency.component : -1;
    if (outaged_generator >= static_cast<int>(ng) ||
        (contingency.type == ContingencyType::Generator && outaged_generator < 0)) {
        throw std::runtime_error("invalid local dispatch generator outage");
    }
    std::vector<double> lower_g(ng), upper_g(ng), lower_t(nd), upper_t(nd);
    std::vector<std::vector<PwlPoint>> gen_curves(ng), load_curves(nd);
    for (std::size_t g = 0; g < ng; ++g) {
        if (commitment[g] != 0 && commitment[g] != 1)
            throw std::runtime_error("local dispatch requires fixed binary commitment");
        if (commitment[g] == 0 || static_cast<int>(g) == outaged_generator) continue;
        const auto& gen = data.generators[g];
        lower_g[g] = std::max(gen.pmin, original_base.pg[g] - data.delta_r_ctg * gen.prdmaxctg);
        upper_g[g] = std::min(gen.pmax, original_base.pg[g] + data.delta_r_ctg * gen.prumaxctg);
        if (lower_g[g] > upper_g[g] || gen.qmin > gen.qmax)
            throw std::runtime_error("inconsistent source local generation bounds");
        gen_curves[g] = active_pwl_points(gen.cost, gen.ncost, lower_g[g], upper_g[g]);
    }
    for (std::size_t l = 0; l < nd; ++l) {
        const auto& load = data.loads[l];
        lower_t[l] = load.tmin; upper_t[l] = load.tmax;
        if (std::abs(load.pd_nominal) > 1e-12) {
            const double base_power = load.pd_nominal * original_base.demand_factor[l];
            // Match the source model's factor/ramp intersection exactly.
            lower_t[l] = std::max(lower_t[l], (base_power - data.delta_r_ctg * load.prdmaxctg) / load.pd_nominal);
            upper_t[l] = std::min(upper_t[l], (base_power + data.delta_r_ctg * load.prumaxctg) / load.pd_nominal);
        }
        if (lower_t[l] > upper_t[l]) continue; // No candidate from an inconsistent interval.
        load_curves[l] = active_pwl_points(load.cost, load.ncost, load.pd_min, load.pd_max);
        // The unchanged PWL power equation also requires the power to lie in
        // the source curve's represented interval (including its existing
        // source-format boundary extension). Never flatten an extrapolation.
        if (std::abs(load.pd_nominal) > 1e-12) {
            const double a = load_curves[l].front().mw / load.pd_nominal;
            const double b = load_curves[l].back().mw / load.pd_nominal;
            lower_t[l] = std::max(lower_t[l], std::min(a, b));
            upper_t[l] = std::min(upper_t[l], std::max(a, b));
            if (lower_t[l] > upper_t[l]) load_curves[l].clear();
        }
    }
    auto rp = network_p, rq = network_q;
    for (std::size_t g = 0; g < ng; ++g) {
        rp[data.generators[g].bus] -= candidate.pg[g];
        rq[data.generators[g].bus] -= candidate.qg[g];
    }
    for (std::size_t l = 0; l < nd; ++l) {
        rp[data.loads[l].bus] += data.loads[l].pd_nominal * candidate.demand_factor[l];
        rq[data.loads[l].bus] += data.loads[l].qd_nominal * candidate.demand_factor[l];
    }
    LocalBusDispatchResult result;
    for (int pass = 0; pass < passes; ++pass) {
        const int changes_before = result.active_generation_changes + result.reactive_generation_changes + result.load_changes;
        for (std::size_t g = 0; g < ng; ++g) {
            if (gen_curves[g].empty()) continue;
            const auto& gen = data.generators[g];
            const int bus = gen.bus;
            const double q = std::clamp(candidate.qg[g] + rq[bus], gen.qmin, gen.qmax);
            const double q_change = q - candidate.qg[g];
            const double q_gain = data.q_delta_cost_approx * (std::abs(rq[bus]) - std::abs(rq[bus] - q_change));
            if (q_gain > 1e-9) {
                candidate.qg[g] = q; rq[bus] -= q_change;
                ++result.reactive_generation_changes;
                result.predicted_native_gain += data.delta_ctg * q_gain;
            }
            const auto best = coordinate_minimum(gen_curves[g], 1.0, 1.0,
                candidate.pg[g], lower_g[g], upper_g[g], rp[bus], rq[bus],
                -1.0, 0.0, data.p_delta_cost_approx, data.q_delta_cost_approx);
            if (best.value != candidate.pg[g]) {
                rp[bus] -= best.value - candidate.pg[g]; candidate.pg[g] = best.value;
                ++result.active_generation_changes;
                result.predicted_native_gain += data.delta_ctg * best.gain;
            }
        }
        for (std::size_t l = 0; l < nd; ++l) {
            if (load_curves[l].empty()) continue;
            const auto& load = data.loads[l];
            const int bus = load.bus;
            const auto best = coordinate_minimum(load_curves[l], load.pd_nominal, -1.0,
                candidate.demand_factor[l], lower_t[l], upper_t[l], rp[bus], rq[bus],
                load.pd_nominal, load.qd_nominal, data.p_delta_cost_approx, data.q_delta_cost_approx);
            if (best.value != candidate.demand_factor[l]) {
                const double change = best.value - candidate.demand_factor[l];
                rp[bus] += load.pd_nominal * change; rq[bus] += load.qd_nominal * change;
                candidate.demand_factor[l] = best.value;
                ++result.load_changes;
                result.predicted_native_gain += data.delta_ctg * best.gain;
            }
        }
        if (changes_before == result.active_generation_changes + result.reactive_generation_changes + result.load_changes) break;
    }
    return result;
}

void run_local_bus_dispatch_regression() {
    const auto require = [](bool condition, const char* message) {
        if (!condition) throw std::runtime_error(message);
    };
    // A dense one-dimensional oracle also covers nonconvex source slopes,
    // signed reactive demand, simultaneous P/Q kinks and active slack caps.
    for (const std::vector<PwlPoint>& curve : {
            std::vector<PwlPoint>{{0, 0}, {0.7, 1}, {1.3, 9}, {2, 13}},
            std::vector<PwlPoint>{{0, 0}, {0.7, 8}, {1.3, 3}, {2, 7}}}) {
        for (double sign : {-1.0, 1.0}) for (double p_slope : {-1.0, 1.0})
        for (double q_slope : {-0.4, 0.0, 0.3}) {
            const double current = 1.0, rp = 0.2, rq = -0.1;
            const auto best = coordinate_minimum(curve, 1.0, sign, current,
                0.4, 1.6, rp, rq, p_slope, q_slope, 12.0, 20.0);
            const auto cost = [&](double x) { return sign * curve_value(curve, x) +
                12.0 * std::abs(rp + p_slope * (x - current)) + 20.0 * std::abs(rq + q_slope * (x - current)); };
            for (int i = 0; i <= 12000; ++i) {
                const double x = 0.4 + i * 0.0001;
                if (std::abs(rp + p_slope * (x - current)) <= balance_limit &&
                    std::abs(rq + q_slope * (x - current)) <= balance_limit) {
                    require(cost(best.value) <= cost(x) + 1e-9, "local coordinate missed PWL minimum");
                }
            }
            require(best.gain >= 0.0 && best.value >= 0.4 && best.value <= 1.6,
                "local coordinate left bounds or lost objective");
        }
    }
    CaseData data;
    data.buses.resize(2); data.delta_ctg = 1.0; data.delta_r_ctg = 1.0;
    data.p_delta_cost_approx = 1000; data.q_delta_cost_approx = 1000;
    Generator gen;
    gen.bus = 0; gen.pmin = 0.5; gen.pmax = 2.0; gen.qmin = -0.2; gen.qmax = 0.2;
    gen.prumaxctg = 0.05; gen.prdmaxctg = 0.03; gen.ncost = 3;
    gen.cost = {0.0, 0.0, 1.0, 10.0, 2.0, 30.0};
    data.generators = {gen, gen, gen};
    Load load;
    load.bus = 1; load.pd_nominal = 1.0; load.qd_nominal = -0.25;
    load.tmin = 0.7; load.tmax = 1.3; load.pd_min = 0.7; load.pd_max = 1.3;
    load.prumaxctg = 0.08; load.prdmaxctg = 0.06; load.ncost = 2;
    load.cost = {0.0, 0.0, 2.0, 20.0}; data.loads = {load};
    AcState base; base.pg = {1.0, 0.0, 1.0}; base.qg = {0.0, 0.0, 0.0}; base.demand_factor = {1.0};
    base.vm = {1.0, 0.99}; base.va = {0.0, -0.01}; base.pf = {0.4}; base.shunt_steps = {{1}};
    auto candidate = base; candidate.pg[0] = 0.98; candidate.pg[2] = 0.0;
    Contingency outage; outage.type = ContingencyType::Generator; outage.component = 2;
    const auto result = improve_fixed_network_bus_dispatch(data, base, {1,0,1}, outage,
        {1.1,-0.8}, {0.1,0.25}, candidate);
    require(result.changed() && result.predicted_native_gain > 0, "local dispatch did not improve fixture");
    require(candidate.pg[0] >= 0.97 && candidate.pg[0] <= 1.05,
        "local generation changed original-base ramp anchor");
    require(candidate.pg[0] >= gen.pmin && candidate.qg[0] >= gen.qmin && candidate.qg[0] <= gen.qmax,
        "local dispatch broke PMIN or Q bounds");
    require(candidate.pg[1] == 0 && candidate.pg[2] == 0 && candidate.qg[1] == 0 && candidate.qg[2] == 0,
        "local dispatch activated an unavailable or outaged generator");
    require(candidate.demand_factor[0] >= 0.94 && candidate.demand_factor[0] <= 1.08,
        "local load changed source corrective ramp bounds");
    require(candidate.vm == base.vm && candidate.va == base.va && candidate.pf == base.pf &&
        candidate.shunt_steps == base.shunt_steps && base.pg[0] == 1.0 && base.demand_factor[0] == 1.0,
        "local dispatch mutated network controls or original base");
    auto pmin_base = base; pmin_base.pg[0] = 0.51;
    auto pmin_trial = candidate; pmin_trial.pg[0] = 0.52;
    improve_fixed_network_bus_dispatch(data, pmin_base, {1,0,1}, outage,
        {0.4,-0.8}, {0.1,0.25}, pmin_trial);
    require(pmin_trial.pg[0] == gen.pmin && gen.pmin == 0.5,
        "positive source PMIN was not the exact limiting bound");
    auto negative_q = candidate;
    improve_fixed_network_bus_dispatch(data, base, {1,0,1}, outage,
        {1.1,-0.8}, {-0.4,0.25}, negative_q);
    require(negative_q.qg[0] == gen.qmin, "negative Q bound not respected");
    const std::vector<PwlPoint> zero_power_curve{{0.0,0.0},{1.0,1.0}};
    const auto reactive_only = coordinate_minimum(zero_power_curve, 0.0, -1.0,
        1.0, 0.5, 1.5, 0.1, 0.2, 0.0, 0.5, 12.0, 20.0);
    require(std::abs(reactive_only.value - 0.6) < 1e-12,
        "zero-active-power reactive coordinate was mishandled");
    auto invalid = candidate; invalid.pg[0] = std::numeric_limits<double>::quiet_NaN();
    bool rejected = false;
    try { improve_fixed_network_bus_dispatch(data, base, {1,0,1}, outage, {1.1,-0.8}, {0.1,0.25}, invalid); }
    catch (const std::runtime_error&) { rejected = true; }
    require(rejected, "nonfinite local dispatch state accepted");
}

}  // namespace gravityx
