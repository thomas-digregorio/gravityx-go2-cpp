#pragma once

#include <nlohmann/json.hpp>

#include <string>
#include <utility>
#include <unordered_map>
#include <vector>

namespace gravityx {

struct PwlPoint {
    double mw{};
    double cost{};
};

struct Bus {
    std::string source_key;
    int index{};
    int bus_i{};
    int type{};
    double vmin{};
    double vmax{};
    double vm_start{};
    double va_start{};
    std::vector<int> generators;
    std::vector<int> loads;
    std::vector<int> shunts;
    std::vector<int> branches_from;
    std::vector<int> branches_to;
    bool present{true};
};

struct Generator {
    std::string source_key;
    int index{};
    int bus{};
    int status_prev{};
    int suqual{};
    int sdqual{};
    double pg_start{};
    double qg_start{};
    double pg_prev{};
    double pmin{};
    double pmax{};
    double qmin{};
    double qmax{};
    double prumax{};
    double prdmax{};
    double prumaxctg{};
    double prdmaxctg{};
    double oncost{};
    double sucost{};
    double sdcost{};
    int ncost{};
    std::vector<double> cost;
    bool present{true};
    int source_bus{};
    std::string source_id;
};

struct Load {
    std::string source_key;
    int index{};
    int bus{};
    double pd_nominal{};
    double qd_nominal{};
    double pd_prev{};
    double qd_prev{};
    double pd_min{};
    double pd_max{};
    double tmin{};
    double tmax{};
    double prumax{};
    double prdmax{};
    double prumaxctg{};
    double prdmaxctg{};
    double z_start{};
    int ncost{};
    std::vector<double> cost;
    bool present{true};
    int source_bus{};
    std::string source_id;
};

struct Shunt {
    std::string source_key;
    int index{};
    int bus{};
    double gs{};
    double bs{};
    bool dispatchable{};
    std::vector<int> steps;
    std::vector<int> block_maximum_steps;
    std::vector<double> block_susceptance;
    bool present{true};
    int source_bus{};
};

struct Branch {
    std::string source_key;
    int index{};
    int status{1};
    int from{};
    int to{};
    bool transformer{};
    double r{};
    double x{};
    double g_fr{};
    double b_fr{};
    double g_to{};
    double b_to{};
    double tap{};
    double shift{};
    // Immutable AC-flow coefficients derived from r/x, tap, phase shift,
    // and terminal admittances at ingest. Fast contingency screening rebuilds
    // branch flows many times for the same network, so recomputing these
    // quantities (including sin/cos of the fixed phase shift) in every trial
    // is pure overhead. Hand-built test fixtures may leave the cache invalid;
    // flow routines retain the canonical formula as a fallback.
    bool flow_coefficients_valid{false};
    double flow_from_g_self{};
    double flow_from_b_self{};
    double flow_to_g_self{};
    double flow_to_b_self{};
    double flow_from_cross_cos{};
    double flow_from_cross_sin{};
    double flow_to_cross_cos{};
    double flow_to_cross_sin{};
    double angmin{};
    double angmax{};
    double rate_a{};
    double rate_b{};
    double rate_c{};
    bool present{true};
    int source_from{};
    int source_to{};
    std::string source_id;
    int control_mode{};
    double tm_min{};
    double tm_max{};
    int tm_steps{1};
    int tm_step{};
    double ta_min{};
    double ta_max{};
    int ta_steps{1};
    int ta_step{};
};

enum class ContingencyType {
    Generator,
    Branch,
};

struct Contingency {
    std::string label;
    ContingencyType type{};
    int source_index{};
    int component{};
};

struct CaseData {
    std::string name;
    double base_mva{};
    double delta{};
    double delta_r{};
    double delta_ctg{};
    double delta_r_ctg{};
    double sm_vio_limit{};
    double sm_cost_approx{};
    double p_delta_cost_approx{};
    double q_delta_cost_approx{};

    std::vector<Bus> buses;
    std::vector<Generator> generators;
    std::vector<Load> loads;
    std::vector<Shunt> shunts;
    std::vector<Branch> branches;
    std::vector<Contingency> contingencies;

    std::unordered_map<int, int> bus_position;

    static CaseData load(const std::string& path);
};

std::vector<PwlPoint> active_pwl_points(
    const std::vector<double>& flat_cost,
    int ncost,
    double pmin,
    double pmax,
    double tolerance = 1e-2);

// A redundant scalar bound implied by the official GO2 apparent-current
// constraint and the source voltage/soft-margin bounds.  It is intentionally
// wider than RATE_A/RATE_C: an individual P or Q component may exceed the
// rating while the permitted apparent-current slack remains feasible.
double branch_terminal_component_bound(
    const CaseData& data,
    const Branch& branch,
    double rating,
    bool from_terminal);

// GO Challenge 2 controllable transformers use an integer position centered
// on zero.  The physical tap ratio or phase shift is the midpoint of the
// source range plus the integer position times the uniform source step size.
// These helpers retain that exact source discretization and refresh every
// cached AC-flow coefficient after a position change.
std::pair<int, int> transformer_step_bounds(const Branch& branch);
double transformer_step_value(const Branch& branch, int step);
void refresh_branch_flow_coefficients(Branch& branch);
void set_transformer_step(Branch& branch, int step);

}  // namespace gravityx
