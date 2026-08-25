#pragma once

#include <nlohmann/json.hpp>

#include <string>
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
    double angmin{};
    double angmax{};
    double rate_a{};
    double rate_b{};
    double rate_c{};
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

}  // namespace gravityx
