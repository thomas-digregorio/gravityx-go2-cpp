#include "gravityx/ac_model.hpp"
#include "gravityx/active_feasibility_repair.hpp"
#include "gravityx/algorithm.hpp"
#include "gravityx/case_data.hpp"
#include "gravityx/fast_power_flow.hpp"
#include "gravityx/linearized_ac_seed.hpp"
#include "gravityx/solution_writer.hpp"
#include "gravityx/state_io.hpp"
#include "gravityx/validation.hpp"

#include <gravity/solver.h>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>

using namespace gravity;

namespace {

void reject_onedrive(const std::string& path) {
    std::string normalized = path;
    std::transform(normalized.begin(), normalized.end(), normalized.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (normalized.find("onedrive") != std::string::npos) {
        throw std::runtime_error("refusing OneDrive path: " + path);
    }
}

nlohmann::json read_json_file(const std::string& path) {
    reject_onedrive(path);
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open JSON input: " + path);
    }
    nlohmann::json value;
    input >> value;
    return value;
}

void write_json_file(const std::string& path, const nlohmann::json& value) {
    reject_onedrive(path);
    const std::filesystem::path output(path);
    if (output.has_parent_path()) {
        std::filesystem::create_directories(output.parent_path());
    }
    const auto temporary = output.string() + ".tmp";
    {
        std::ofstream stream(temporary, std::ios::trunc);
        if (!stream) {
            throw std::runtime_error("cannot open JSON output: " + temporary);
        }
        const char* pretty_json = std::getenv("GRAVITYX_PRETTY_JSON");
        const int indentation = pretty_json != nullptr &&
                std::string(pretty_json) != "0"
            ? 2 : -1;
        stream << value.dump(indentation) << '\n';
        if (!stream) {
            throw std::runtime_error("failed while writing JSON output: " + temporary);
        }
    }
    std::filesystem::rename(temporary, output);
}

void require_near(double actual, double expected, double tolerance, const std::string& name) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        throw std::runtime_error("component test failed: " + name);
    }
}

bool bounded_fast_newton_rescue_candidate(
    const gravityx::ValidationReport& validation) {
    return validation.max_residual <= 0.35 &&
        (validation.worst_category == "active_balance" ||
         validation.worst_category == "reactive_balance");
}

bool bounded_fast_candidate_repair_candidate(
    const gravityx::ValidationReport& validation) {
    return bounded_fast_newton_rescue_candidate(validation) ||
        (validation.max_residual <= 0.1 &&
         validation.worst_category == "variable_bound");
}

struct PassivePocket {
    int boundary_branch{-1};
    int outside_bus{-1};
    std::vector<int> buses;
};

std::vector<PassivePocket> find_small_passive_outage_pockets(
    const gravityx::CaseData& data,
    int outaged_branch,
    std::size_t maximum_pocket_buses = 64,
    nlohmann::json* diagnostics = nullptr) {
    std::vector<PassivePocket> pockets;
    if (outaged_branch < 0 ||
        outaged_branch >= static_cast<int>(data.branches.size())) {
        return pockets;
    }
    const auto& outage = data.branches[outaged_branch];
    const int nb = static_cast<int>(data.buses.size());
    const int nl = static_cast<int>(data.branches.size());
    const auto incident_branches = [&](int bus) {
        std::vector<int> result = data.buses[bus].branches_from;
        result.insert(
            result.end(), data.buses[bus].branches_to.begin(),
            data.buses[bus].branches_to.end());
        return result;
    };
    const auto other_bus = [&](int branch, int bus) {
        return data.branches[branch].from == bus
            ? data.branches[branch].to : data.branches[branch].from;
    };

    // Only bridges close to an endpoint of the newly opened branch can bound
    // the passive pocket created by that outage.  Limit the discovery radius
    // and every bridge-side search so this remains negligible on the largest
    // cases.
    constexpr int kDiscoveryDepth = 4;
    std::vector<int> depth(static_cast<std::size_t>(nb), -1);
    std::vector<unsigned char> candidate_seen(
        static_cast<std::size_t>(nl), 0);
    std::vector<int> candidate_branches;
    std::queue<int> discovery;
    for (int endpoint : {outage.from, outage.to}) {
        if (depth[endpoint] < 0) {
            depth[endpoint] = 0;
            discovery.push(endpoint);
        }
    }
    while (!discovery.empty()) {
        const int bus = discovery.front();
        discovery.pop();
        for (int branch : incident_branches(bus)) {
            if (branch == outaged_branch ||
                data.branches[branch].status == 0) {
                continue;
            }
            if (!candidate_seen[branch]) {
                candidate_seen[branch] = 1;
                candidate_branches.push_back(branch);
            }
            const int neighbor = other_bus(branch, bus);
            if (depth[bus] < kDiscoveryDepth && depth[neighbor] < 0) {
                depth[neighbor] = depth[bus] + 1;
                discovery.push(neighbor);
            }
        }
    }

    std::vector<int> visited(static_cast<std::size_t>(nb), 0);
    int visit_token = 0;
    int rejected_oversized_or_nonbridge = 0;
    int rejected_without_outage_endpoint = 0;
    int rejected_nonpassive = 0;
    int rejected_nonunit_ratio = 0;
    nlohmann::json rejected_nonpassive_details =
        nlohmann::json::array();
    std::vector<unsigned char> accepted_side(
        static_cast<std::size_t>(2 * nl), 0);
    for (int boundary_branch : candidate_branches) {
        const auto& boundary = data.branches[boundary_branch];
        for (int side = 0; side < 2; ++side) {
            const int start = side == 0 ? boundary.from : boundary.to;
            const int outside = side == 0 ? boundary.to : boundary.from;
            const int side_key = 2 * boundary_branch + side;
            if (accepted_side[side_key]) {
                continue;
            }
            ++visit_token;
            std::queue<int> queue;
            std::vector<int> component;
            visited[start] = visit_token;
            queue.push(start);
            bool oversized = false;
            while (!queue.empty() && !oversized) {
                const int bus = queue.front();
                queue.pop();
                component.push_back(bus);
                if (component.size() > maximum_pocket_buses) {
                    oversized = true;
                    break;
                }
                for (int branch : incident_branches(bus)) {
                    if (branch == outaged_branch ||
                        branch == boundary_branch ||
                        data.branches[branch].status == 0) {
                        continue;
                    }
                    const int neighbor = other_bus(branch, bus);
                    if (visited[neighbor] != visit_token) {
                        visited[neighbor] = visit_token;
                        queue.push(neighbor);
                    }
                }
            }
            if (oversized || visited[outside] == visit_token) {
                ++rejected_oversized_or_nonbridge;
                continue;
            }
            const bool contains_outage_endpoint =
                visited[outage.from] == visit_token ||
                visited[outage.to] == visit_token;
            if (!contains_outage_endpoint) {
                ++rejected_without_outage_endpoint;
                continue;
            }
            bool passive = true;
            bool unit_ratio = true;
            for (int bus : component) {
                const auto& source_bus = data.buses[bus];
                // Loads and source-dispatchable shunts remain subject to their
                // exact corrective bounds below.  Do not apply this candidate
                // transformation to generation or a reference bus.
                if (!source_bus.generators.empty() || source_bus.type == 3) {
                    rejected_nonpassive_details.push_back({
                        {"boundary_branch",
                         data.branches[boundary_branch].source_key},
                        {"bus", source_bus.source_key},
                        {"generator_count", source_bus.generators.size()},
                        {"load_count", source_bus.loads.size()},
                        {"shunt_count", source_bus.shunts.size()},
                        {"bus_type", source_bus.type},
                        {"component_size", component.size()},
                    });
                    passive = false;
                    break;
                }
                for (int branch : incident_branches(bus)) {
                    if (branch == outaged_branch ||
                        data.branches[branch].status == 0) {
                        continue;
                    }
                    const auto& source_branch = data.branches[branch];
                    if (std::abs(source_branch.tap - 1.0) > 1e-9 ||
                        std::abs(source_branch.shift) > 1e-9) {
                        unit_ratio = false;
                        break;
                    }
                }
                if (!unit_ratio) {
                    break;
                }
            }
            if (!passive) {
                ++rejected_nonpassive;
                continue;
            }
            if (!unit_ratio) {
                ++rejected_nonunit_ratio;
                continue;
            }
            std::sort(component.begin(), component.end());
            accepted_side[side_key] = 1;
            pockets.push_back({
                boundary_branch,
                outside,
                std::move(component),
            });
        }
    }
    if (diagnostics != nullptr) {
        *diagnostics = {
            {"candidate_branch_count", candidate_branches.size()},
            {"accepted_pocket_count", pockets.size()},
            {"rejected_oversized_or_nonbridge",
             rejected_oversized_or_nonbridge},
            {"rejected_without_outage_endpoint",
             rejected_without_outage_endpoint},
            {"rejected_nonpassive", rejected_nonpassive},
            {"rejected_nonpassive_details",
             std::move(rejected_nonpassive_details)},
            {"rejected_nonunit_ratio", rejected_nonunit_ratio},
        };
    }
    return pockets;
}

struct PassivePocketRepair {
    std::optional<gravityx::FastPowerFlowResult> result;
    nlohmann::json diagnostics;
};

std::optional<PassivePocketRepair> try_passive_outage_pocket_repair(
    const gravityx::CaseData& data,
    const gravityx::AcState& base_state,
    const std::vector<int>& commitment,
    const gravityx::Contingency& contingency,
    const gravityx::ContingencyContext& context,
    const gravityx::AcState& reference,
    const gravityx::ValidationReport& reference_validation) {
    if (contingency.type != gravityx::ContingencyType::Branch) {
        return std::nullopt;
    }
    const auto wall_start = std::chrono::steady_clock::now();
    nlohmann::json search_diagnostics;
    const auto pockets = find_small_passive_outage_pockets(
        data, contingency.component, 64, &search_diagnostics);
    nlohmann::json candidates = nlohmann::json::array();
    int rejected_voltage_bounds = 0;
    std::optional<gravityx::FastPowerFlowResult> best;
    for (const auto& pocket : pockets) {
        const double boundary_vm = reference.vm[pocket.outside_bus];
        bool voltage_bounds_allow_equalization = true;
        for (int bus : pocket.buses) {
            if (boundary_vm < data.buses[bus].vmin - 1e-12 ||
                boundary_vm > data.buses[bus].vmax + 1e-12) {
                voltage_bounds_allow_equalization = false;
                break;
            }
        }
        if (!voltage_bounds_allow_equalization) {
            ++rejected_voltage_bounds;
            continue;
        }
        auto candidate = reference;
        nlohmann::json shunt_changes = nlohmann::json::array();
        if (candidate.shunt_steps.size() == data.shunts.size() &&
            candidate.shunt_bs.size() == data.shunts.size()) {
            for (int bus : pocket.buses) {
                for (int shunt_index : data.buses[bus].shunts) {
                    const auto& shunt = data.shunts[shunt_index];
                    if (!shunt.dispatchable ||
                        shunt.block_maximum_steps.size() !=
                            shunt.block_susceptance.size()) {
                        continue;
                    }
                    std::size_t configuration_count = 1;
                    for (int maximum : shunt.block_maximum_steps) {
                        if (maximum < 0 ||
                            configuration_count >
                                4096 / static_cast<std::size_t>(maximum + 1)) {
                            configuration_count = 0;
                            break;
                        }
                        configuration_count *=
                            static_cast<std::size_t>(maximum + 1);
                    }
                    if (configuration_count == 0 ||
                        configuration_count > 4096) {
                        continue;
                    }
                    const auto original_steps =
                        candidate.shunt_steps[shunt_index];
                    std::vector<int> best_steps = original_steps;
                    double best_bs = candidate.shunt_bs[shunt_index];
                    int best_movement = std::numeric_limits<int>::max();
                    for (std::size_t code = 0;
                         code < configuration_count; ++code) {
                        std::size_t remaining = code;
                        std::vector<int> steps(
                            shunt.block_maximum_steps.size(), 0);
                        double bs = 0.0;
                        int movement = 0;
                        for (std::size_t block = 0;
                             block < steps.size(); ++block) {
                            const int radix =
                                shunt.block_maximum_steps[block] + 1;
                            steps[block] = static_cast<int>(
                                remaining % static_cast<std::size_t>(radix));
                            remaining /= static_cast<std::size_t>(radix);
                            bs += static_cast<double>(steps[block]) *
                                shunt.block_susceptance[block];
                            if (block < original_steps.size()) {
                                movement += std::abs(
                                    steps[block] - original_steps[block]);
                            }
                        }
                        if (std::abs(bs) + 1e-12 < std::abs(best_bs) ||
                            (std::abs(std::abs(bs) - std::abs(best_bs)) <=
                                 1e-12 &&
                             movement < best_movement)) {
                            best_steps = std::move(steps);
                            best_bs = bs;
                            best_movement = movement;
                        }
                    }
                    if (best_steps != original_steps) {
                        shunt_changes.push_back({
                            {"shunt", shunt.source_key},
                            {"before_bs", candidate.shunt_bs[shunt_index]},
                            {"after_bs", best_bs},
                            {"before_steps", original_steps},
                            {"after_steps", best_steps},
                        });
                        candidate.shunt_steps[shunt_index] =
                            std::move(best_steps);
                        candidate.shunt_bs[shunt_index] = best_bs;
                    }
                }
            }
        }
        for (int bus : pocket.buses) {
            candidate.vm[bus] = boundary_vm;
            candidate.va[bus] = reference.va[pocket.outside_bus];
        }
        const double objective =
            gravityx::rebuild_contingency_state_derived_fields(
                data, base_state, commitment, contingency, candidate);
        const auto validation = gravityx::validate_state(
            data, gravityx::ModelMode::ContingencySoft,
            candidate, commitment, context);
        nlohmann::json bus_keys = nlohmann::json::array();
        for (int bus : pocket.buses) {
            bus_keys.push_back(data.buses[bus].source_key);
        }
        candidates.push_back({
            {"boundary_branch",
             data.branches[pocket.boundary_branch].source_key},
            {"outside_bus", data.buses[pocket.outside_bus].source_key},
            {"pocket_buses", std::move(bus_keys)},
            {"shunt_changes", std::move(shunt_changes)},
            {"validation", validation.to_json()},
        });
        const double best_residual = best
            ? best->validation.max_residual
            : reference_validation.max_residual;
        if (validation.max_residual + 1e-12 >= best_residual) {
            continue;
        }
        gravityx::FastPowerFlowResult result;
        result.converged = false;
        result.feasible = validation.max_residual <= 1e-5;
        result.solve.status = result.feasible ? 0 : 1;
        result.solve.objective = objective;
        result.solve.state = std::move(candidate);
        result.validation = validation;
        result.failure_reason = result.feasible
            ? std::string{}
            : "passive outage-pocket equalization requires further repair";
        best = std::move(result);
    }
    const double wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_start).count();
    if (best) {
        best->wall_seconds = wall_seconds;
        best->solve.wall_seconds = wall_seconds;
    }
    return PassivePocketRepair{
        std::move(best),
        {
            {"search", std::move(search_diagnostics)},
            {"rejected_voltage_bounds", rejected_voltage_bounds},
            {"candidate_count", candidates.size()},
            {"candidates", std::move(candidates)},
            {"wall_seconds", wall_seconds},
        },
    };
}

int run_component_tests() {
    gravityx::run_fast_power_flow_topology_cache_regression();
    const auto points = gravityx::active_pwl_points(
        {0.0, 0.0, 10.0, 100.0, 20.0, 300.0}, 3, 5.0, 15.0);
    if (points.size() != 3) {
        throw std::runtime_error("component test failed: clipped PWL point count");
    }
    require_near(points[0].mw, 0.0, 1e-12, "lower PWL bracket MW");
    require_near(points[0].cost, 0.0, 1e-12, "lower PWL bracket cost");
    require_near(points[1].mw, 10.0, 1e-12, "interior PWL MW");
    require_near(points[2].mw, 20.0, 1e-12, "upper PWL bracket MW");
    require_near(points[2].cost, 300.0, 1e-12, "upper PWL bracket cost");

    gravityx::AcState source;
    source.vm = {1.0, 1.01};
    source.va = {0.0, -0.02};
    source.pg = {0.5};
    source.commitment = {1.0};
    source.gen_lambda = {0.25, 0.75};
    const auto restored = gravityx::ac_state_from_json(gravityx::ac_state_to_json(source));
    if (restored.vm != source.vm || restored.va != source.va || restored.pg != source.pg ||
        restored.commitment != source.commitment || restored.gen_lambda != source.gen_lambda) {
        throw std::runtime_error("component test failed: AC state JSON round trip");
    }
    source.qg = {0.1};
    source.demand_factor = {0.95};
    source.pf = {0.4};
    const auto submission_state = gravityx::ac_submission_state_to_json(source);
    if (submission_state.size() != 5 || !submission_state.contains("vm") ||
        !submission_state.contains("va") || !submission_state.contains("pg") ||
        !submission_state.contains("qg") ||
        !submission_state.contains("demand_factor") ||
        submission_state.contains("pf") || submission_state.contains("gen_lambda")) {
        throw std::runtime_error(
            "component test failed: compact submission state schema");
    }

    gravityx::SolveResult nonconverged_feasible;
    nonconverged_feasible.status = -2;
    nonconverged_feasible.objective = 42.0;
    gravityx::ValidationReport feasible_validation;
    feasible_validation.max_residual = 1e-8;
    if (!gravityx::validated_candidate_is_feasible(
            nonconverged_feasible, feasible_validation, 1e-5)) {
        throw std::runtime_error(
            "component test failed: validated nonconverged candidate was rejected");
    }
    feasible_validation.max_residual = 1e-4;
    if (gravityx::validated_candidate_is_feasible(
            nonconverged_feasible, feasible_validation, 1e-5)) {
        throw std::runtime_error(
            "component test failed: invalid nonconverged candidate was accepted");
    }
    gravityx::ValidationReport rescue_validation;
    rescue_validation.max_residual = 0.15;
    rescue_validation.worst_category = "reactive_balance";
    if (!bounded_fast_newton_rescue_candidate(rescue_validation)) {
        throw std::runtime_error(
            "component test failed: reactive fast-Newton rescue was not eligible");
    }
    rescue_validation.worst_category = "flow_limit";
    if (bounded_fast_newton_rescue_candidate(rescue_validation)) {
        throw std::runtime_error(
            "component test failed: non-balance fast-Newton rescue was eligible");
    }
    rescue_validation.worst_category = "active_balance";
    rescue_validation.max_residual = 0.36;
    if (bounded_fast_newton_rescue_candidate(rescue_validation)) {
        throw std::runtime_error(
            "component test failed: oversized fast-Newton rescue was eligible");
    }
    rescue_validation.worst_category = "variable_bound";
    rescue_validation.max_residual = 0.05;
    if (bounded_fast_newton_rescue_candidate(rescue_validation) ||
        !bounded_fast_candidate_repair_candidate(rescue_validation)) {
        throw std::runtime_error(
            "component test failed: bounded security repair routing failed");
    }
    feasible_validation.max_residual = 1e-8;
    nonconverged_feasible.objective = std::numeric_limits<double>::quiet_NaN();
    if (gravityx::validated_candidate_is_feasible(
            nonconverged_feasible, feasible_validation, 1e-5)) {
        throw std::runtime_error(
            "component test failed: nonfinite candidate was accepted");
    }

    gravityx::CaseData passive_pocket_case;
    passive_pocket_case.buses.resize(4);
    passive_pocket_case.generators.resize(1);
    passive_pocket_case.loads.resize(1);
    passive_pocket_case.buses[0].generators = {0};
    passive_pocket_case.buses[3].loads = {0};
    passive_pocket_case.branches.resize(3);
    for (auto& branch : passive_pocket_case.branches) {
        branch.status = 1;
        branch.tap = 1.0;
        branch.shift = 0.0;
    }
    passive_pocket_case.branches[0].from = 0;
    passive_pocket_case.branches[0].to = 1;
    passive_pocket_case.branches[1].from = 1;
    passive_pocket_case.branches[1].to = 2;
    passive_pocket_case.branches[2].from = 2;
    passive_pocket_case.branches[2].to = 3;
    passive_pocket_case.buses[0].branches_from = {0};
    passive_pocket_case.buses[1].branches_to = {0};
    passive_pocket_case.buses[1].branches_from = {1};
    passive_pocket_case.buses[2].branches_to = {1};
    passive_pocket_case.buses[2].branches_from = {2};
    passive_pocket_case.buses[3].branches_to = {2};
    const auto passive_pockets = find_small_passive_outage_pockets(
        passive_pocket_case, 2, 8);
    const auto expected_pocket = std::find_if(
        passive_pockets.begin(), passive_pockets.end(),
        [](const PassivePocket& pocket) {
            return pocket.boundary_branch == 0 &&
                pocket.outside_bus == 0 &&
                pocket.buses == std::vector<int>({1, 2});
        });
    if (expected_pocket == passive_pockets.end()) {
        throw std::runtime_error(
            "component test failed: passive outage pocket was not found");
    }

    gravityx::CaseData seed_case;
    seed_case.buses.resize(2);
    seed_case.buses[0].generators = {0};
    seed_case.buses[0].branches_from = {0};
    seed_case.buses[1].loads = {0};
    seed_case.buses[1].branches_to = {0};
    seed_case.generators.resize(1);
    seed_case.generators[0].bus = 0;
    seed_case.loads.resize(1);
    seed_case.loads[0].bus = 1;
    seed_case.loads[0].pd_nominal = 0.39;
    seed_case.loads[0].qd_nominal = 0.095;
    seed_case.branches.resize(1);
    seed_case.branches[0].from = 0;
    seed_case.branches[0].to = 1;

    gravityx::AcState balanced;
    balanced.vm = {1.0, 1.0};
    balanced.pg = {0.4};
    balanced.qg = {0.1};
    balanced.demand_factor = {1.0};
    balanced.pf = {0.4};
    balanced.qf = {0.1};
    balanced.pt = {-0.39};
    balanced.qt = {-0.095};
    const auto balanced_seed = gravityx::nodal_balance_slack_seed(
        seed_case, balanced, 0.5, 0.0);
    require_near(balanced_seed.active[0], 0.0, 1e-12, "balanced active seed at from bus");
    require_near(balanced_seed.reactive[1], 0.0, 1e-12, "balanced reactive seed at to bus");

    auto branch_outage = balanced;
    branch_outage.pf[0] = 0.0;
    branch_outage.qf[0] = 0.0;
    branch_outage.pt[0] = 0.0;
    branch_outage.qt[0] = 0.0;
    const auto outage_seed = gravityx::nodal_balance_slack_seed(
        seed_case, branch_outage, 0.5, 0.0);
    require_near(outage_seed.active[0], 0.4, 1e-12, "branch-outage active seed at from bus");
    require_near(outage_seed.active[1], 0.39, 1e-12, "branch-outage active seed at to bus");
    require_near(outage_seed.reactive[0], 0.1, 1e-12, "branch-outage reactive seed at from bus");
    require_near(outage_seed.reactive[1], 0.095, 1e-12, "branch-outage reactive seed at to bus");

    auto generator_outage = balanced;
    generator_outage.pg[0] = 0.0;
    generator_outage.qg[0] = 0.0;
    const auto generator_seed = gravityx::nodal_balance_slack_seed(
        seed_case, generator_outage, 0.5, 0.0);
    require_near(generator_seed.active[0], 0.4, 1e-12, "generator-outage active seed");
    require_near(generator_seed.reactive[0], 0.1, 1e-12, "generator-outage reactive seed");

    gravityx::CaseData writer_case;
    writer_case.buses.resize(2);
    writer_case.buses[0].index = 101;
    writer_case.buses[0].bus_i = 101;
    writer_case.buses[1].index = 202;
    writer_case.buses[1].bus_i = 202;
    gravityx::Load writer_load;
    writer_load.index = 3;
    writer_load.bus = 1;
    writer_load.source_bus = 202;
    writer_load.source_id = "L1";
    writer_case.loads.push_back(writer_load);
    gravityx::Generator writer_generator;
    writer_generator.index = 7;
    writer_generator.bus = 0;
    writer_generator.source_bus = 101;
    writer_generator.source_id = "G1";
    writer_case.generators.push_back(writer_generator);
    gravityx::Branch writer_line;
    writer_line.index = 11;
    writer_line.from = 0;
    writer_line.to = 1;
    writer_line.source_from = 101;
    writer_line.source_to = 202;
    writer_line.source_id = "1";
    writer_line.status = 1;
    writer_case.branches.push_back(writer_line);
    gravityx::Branch writer_transformer = writer_line;
    writer_transformer.index = 12;
    writer_transformer.transformer = true;
    writer_transformer.source_id = "T1";
    writer_transformer.control_mode = 3;
    writer_transformer.ta_step = 2;
    writer_case.branches.push_back(writer_transformer);
    gravityx::Shunt writer_shunt;
    writer_shunt.index = 5;
    writer_shunt.bus = 0;
    writer_shunt.source_bus = 101;
    writer_shunt.dispatchable = true;
    writer_shunt.steps = {0, 0};
    writer_case.shunts.push_back(writer_shunt);
    gravityx::AcState writer_state;
    writer_state.vm = {1.0, 0.5};
    writer_state.va = {0.0, -0.125};
    writer_state.demand_factor = {1.0};
    writer_state.pg = {0.75};
    writer_state.qg = {0.25};
    writer_state.shunt_steps = {{2, 0}};
    gravityx::Contingency writer_outage;
    writer_outage.label = "line-outage";
    writer_outage.type = gravityx::ContingencyType::Branch;
    writer_outage.source_index = 11;
    writer_outage.component = 0;
    const std::string expected_solution =
        "--bus section\n"
        "i, v, theta\n"
        "101, 1, 0\n"
        "202, 0.5, -0.125\n"
        "--load section\n"
        "i, id, t\n"
        "202, L1, 1\n"
        "--generator section\n"
        "i, id, p, q, x\n"
        "101, G1, 0.75, 0.25, 1\n"
        "--line section\n"
        "iorig, idest, id, x\n"
        "--transformer section\n"
        "iorig, idest, id, x, xst\n"
        "101, 202, T1, 1, 2\n"
        "--switched shunt section\n"
        "i, xst1, xst2, xst3, xst4, xst5, xst6, xst7, xst8\n"
        "101, 2, 0\n";
    if (gravityx::go_solution_text(
            writer_case, writer_state, {1}, &writer_outage) !=
        expected_solution) {
        throw std::runtime_error(
            "component test failed: official solution text schema");
    }
    writer_state.vm[0] = 1.123456789012345;
    const auto precision_text = gravityx::go_solution_text(
        writer_case, writer_state, {1}, &writer_outage);
    if (precision_text.find("101, 1.12345678901, 0\n") ==
        std::string::npos) {
        throw std::runtime_error(
            "component test failed: bounded solution precision");
    }
    std::cout << "component tests passed\n";
    return 0;
}

int run_parallel_circuit_regression() {
    gravityx::CaseData data;
    data.name = "parallel-circuit-regression";
    data.base_mva = 100.0;
    data.delta = 1.0;
    data.delta_r = 1.0;
    data.delta_ctg = 1.0;
    data.delta_r_ctg = 1.0;
    data.sm_vio_limit = 0.2;
    data.sm_cost_approx = 1e5;
    data.p_delta_cost_approx = 1e5;
    data.q_delta_cost_approx = 1e5;
    data.buses = {
        {"1", 1, 1, 3, 0.9, 1.1, 1.0, 0.0},
        {"2", 2, 2, 1, 0.9, 1.1, 1.0, -0.05},
    };
    gravityx::Generator generator;
    generator.source_key = "1";
    generator.index = 1;
    generator.bus = 0;
    generator.status_prev = 1;
    generator.suqual = 1;
    generator.sdqual = 1;
    generator.pg_start = 0.5;
    generator.pg_prev = 0.5;
    generator.pmin = 0.0;
    generator.pmax = 2.0;
    generator.qmin = -2.0;
    generator.qmax = 2.0;
    generator.prumax = generator.prdmax = 2.0;
    generator.prumaxctg = generator.prdmaxctg = 2.0;
    generator.ncost = 2;
    generator.cost = {0.0, 0.0, 2.0, 20.0};
    data.generators.push_back(generator);

    gravityx::Load load;
    load.source_key = "1";
    load.index = 1;
    load.bus = 1;
    load.pd_nominal = load.pd_prev = 0.5;
    load.qd_nominal = load.qd_prev = 0.1;
    load.pd_min = 0.0;
    load.pd_max = 1.0;
    load.tmin = load.tmax = load.z_start = 1.0;
    load.prumax = load.prdmax = 1.0;
    load.prumaxctg = load.prdmaxctg = 1.0;
    load.ncost = 2;
    load.cost = {0.0, 0.0, 1.0, 1000.0};
    data.loads.push_back(load);

    gravityx::Branch branch;
    branch.from = 0;
    branch.to = 1;
    branch.r = 0.01;
    branch.x = 0.1;
    branch.tap = 1.0;
    branch.angmin = -1.0;
    branch.angmax = 1.0;
    branch.rate_a = 2.0;
    branch.rate_b = 2.0;
    branch.rate_c = 2.0;
    branch.source_key = "1";
    branch.index = 1;
    data.branches.push_back(branch);
    branch.source_key = "2";
    branch.index = 2;
    data.branches.push_back(branch);

    data.buses[0].generators = {0};
    data.buses[1].loads = {0};
    data.buses[0].branches_from = {0, 1};
    data.buses[1].branches_to = {0, 1};

    const auto source_base = gravityx::build_validated_source_base(data, {1});
    if (!source_base.feasible || source_base.validation.max_residual > 1e-5) {
        throw std::runtime_error(
            "validated source-base regression failed with residual " +
            std::to_string(source_base.validation.max_residual));
    }
    auto inactive_branch_data = data;
    inactive_branch_data.branches[0].status = 0;
    // Deliberately make the inactive circuit electrically extreme.  Every
    // model layer must still keep it at exactly zero flow and exclude it from
    // the network equations.
    inactive_branch_data.branches[0].r = 1e-10;
    inactive_branch_data.branches[0].x = 1e-10;
    const auto inactive_source = gravityx::build_validated_source_base(
        inactive_branch_data, {1});
    if (!inactive_source.feasible ||
        inactive_source.validation.max_residual > 1e-5 ||
        std::abs(inactive_source.solve.state.pf[0]) > 1e-12 ||
        std::abs(inactive_source.solve.state.qf[0]) > 1e-12 ||
        std::abs(inactive_source.solve.state.pt[0]) > 1e-12 ||
        std::abs(inactive_source.solve.state.qt[0]) > 1e-12 ||
        std::abs(inactive_source.solve.state.sm_slack[0]) > 1e-12) {
        throw std::runtime_error(
            "inactive-branch source-topology regression failed");
    }
    gravityx::AcModel inactive_model(
        inactive_branch_data, gravityx::ModelMode::BaseSoft, {1});
    const auto inactive_solve = inactive_model.solve(0, 1e-7);
    const auto inactive_validation = gravityx::validate_state(
        inactive_branch_data, gravityx::ModelMode::BaseSoft,
        inactive_solve.state, {1});
    if ((inactive_solve.status != 0 && inactive_solve.status != 1) ||
        inactive_validation.max_residual > 1e-5 ||
        std::abs(inactive_solve.state.pf[0]) > 1e-12 ||
        std::abs(inactive_solve.state.qf[0]) > 1e-12 ||
        std::abs(inactive_solve.state.pt[0]) > 1e-12 ||
        std::abs(inactive_solve.state.qt[0]) > 1e-12 ||
        std::abs(inactive_solve.state.sm_slack[0]) > 1e-12) {
        throw std::runtime_error(
            "inactive-branch exact-topology regression failed");
    }
    const auto linear_seed = gravityx::solve_linearized_ac_seed(
        data, source_base.solve.state, {1});
    if (!linear_seed.success) {
        throw std::runtime_error(
            "linearized AC seed regression failed: " + linear_seed.status);
    }
    const auto extended_linear_seed = gravityx::solve_linearized_ac_seed(
        data, source_base.solve.state, {1}, 0.49, std::nullopt, false, false,
        90.0);
    if (!extended_linear_seed.success ||
        std::abs(extended_linear_seed.time_limit_seconds - 90.0) > 1e-12) {
        throw std::runtime_error(
            "extended linearized AC time-limit regression failed");
    }
    constexpr double kTestVoltageTrustRadius = 0.04;
    constexpr double kTestAngleTrustRadius = 0.10;
    const auto lightweight_seed = gravityx::solve_linearized_ac_seed(
        data, source_base.solve.state, {1}, 0.49, std::nullopt, false, true,
        60.0, false, false, {}, false,
        kTestVoltageTrustRadius, kTestAngleTrustRadius);
    if (!lightweight_seed.success ||
        std::abs(lightweight_seed.voltage_trust_radius -
                 kTestVoltageTrustRadius) > 1e-12 ||
        std::abs(lightweight_seed.angle_trust_radius -
                 kTestAngleTrustRadius) > 1e-12) {
        throw std::runtime_error(
            "lightweight linearized AC seed regression failed: " +
            lightweight_seed.status);
    }
    for (int i = 0; i < static_cast<int>(data.buses.size()); ++i) {
        if (std::abs(lightweight_seed.state.vm[i] - source_base.solve.state.vm[i]) >
                kTestVoltageTrustRadius + 1e-9 ||
            std::abs(lightweight_seed.state.va[i] - source_base.solve.state.va[i]) >
                kTestAngleTrustRadius + 1e-9) {
            throw std::runtime_error(
                "lightweight linearized AC seed left its trust region");
        }
    }
    auto outside_voltage_reference = source_base.solve.state;
    outside_voltage_reference.vm[0] = data.buses[0].vmax + 0.2;
    const auto projected_voltage_seed = gravityx::solve_linearized_ac_seed(
        data, outside_voltage_reference, {1}, 0.49, std::nullopt, false, true,
        60.0, false, false, {}, false,
        kTestVoltageTrustRadius, kTestAngleTrustRadius);
    if (projected_voltage_seed.projected_reference_voltage_count != 1 ||
        std::abs(projected_voltage_seed.maximum_reference_voltage_projection -
                 0.2) > 1e-10 ||
        projected_voltage_seed.row_count <= 0 ||
        projected_voltage_seed.column_count <= 0 ||
        projected_voltage_seed.status.empty() ||
        (projected_voltage_seed.success &&
         (projected_voltage_seed.state.vm[0] < data.buses[0].vmin - 1e-9 ||
          projected_voltage_seed.state.vm[0] > data.buses[0].vmax + 1e-9))) {
        throw std::runtime_error(
            "out-of-bound reference-voltage projection regression failed: " +
            projected_voltage_seed.status);
    }
    const auto balance_only_base_seed = gravityx::solve_linearized_ac_seed(
        data, source_base.solve.state, {1}, 0.25, std::nullopt,
        true, true, 60.0, true, true);
    if (!balance_only_base_seed.success ||
        !balance_only_base_seed.projected_balance_slack ||
        !balance_only_base_seed.branch_security_rows_omitted ||
        !balance_only_base_seed.feasibility_only ||
        balance_only_base_seed.branch_security_subset_count != 0 ||
        balance_only_base_seed.row_count >= linear_seed.row_count ||
        balance_only_base_seed.column_count >= linear_seed.column_count) {
        throw std::runtime_error(
            "balance-only base Phase-I regression failed");
    }
    gravityx::AcModel model(data, gravityx::ModelMode::BaseSoft, {1});
    const auto solve = model.solve(0, 1e-7);
    const auto validation = gravityx::validate_state(
        data, gravityx::ModelMode::BaseSoft, solve.state, {1});
    if ((solve.status != 0 && solve.status != 1) || !std::isfinite(solve.objective) ||
        validation.max_residual > 1e-5) {
        throw std::runtime_error("parallel-circuit symbolic-DAG regression failed");
    }
    gravityx::Contingency branch_contingency;
    branch_contingency.label = "parallel-outage";
    branch_contingency.type = gravityx::ContingencyType::Branch;
    branch_contingency.source_index = 1;
    branch_contingency.component = 0;
    gravityx::FastContingencyPowerFlow fast_screen(
        data, solve.state, {1});
    const auto fast_result = fast_screen.solve(branch_contingency);
    if (!fast_result.feasible || fast_result.validation.max_residual > 1e-5) {
        throw std::runtime_error(
            "validated fast power-flow contingency regression failed: "
            + fast_result.failure_reason);
    }
    auto deliberately_bad_rolling_candidate = solve.state;
    deliberately_bad_rolling_candidate.va[1] += 0.5;
    const auto direct_rolling_screen = fast_screen.screen_candidate(
        branch_contingency, deliberately_bad_rolling_candidate);
    if (direct_rolling_screen.feasible ||
        direct_rolling_screen.fixed_jacobian_predictor_attempted ||
        direct_rolling_screen.direct_candidate_validation.max_residual <=
            1e-5 ||
        direct_rolling_screen.direct_candidate_validation
            .worst_identity.empty()) {
        throw std::runtime_error(
            "rolling corrective candidate was not screened directly");
    }
    auto globally_shifted_candidate = fast_result.solve.state;
    for (double& angle : globally_shifted_candidate.va) {
        angle += 0.123;
    }
    const auto normalized_gauge_screen = fast_screen.screen_candidate(
        branch_contingency, globally_shifted_candidate);
    if (!normalized_gauge_screen.feasible ||
        normalized_gauge_screen.validation.max_residual > 1e-5 ||
        normalized_gauge_screen.validation.max_reference_angle_residual >
            1e-12) {
        throw std::runtime_error(
            "source reference-angle gauge normalization regression failed: " +
            normalized_gauge_screen.failure_reason);
    }
    auto contingency_rating_data = data;
    contingency_rating_data.branches[1].rate_a = 0.1;
    contingency_rating_data.branches[1].rate_c = 2.0;
    if (std::hypot(fast_result.solve.state.pf[1],
                   fast_result.solve.state.qf[1]) <= 0.1) {
        throw std::runtime_error(
            "contingency-rating regression fixture does not exceed RATE_A");
    }
    auto base_rating_state = fast_result.solve.state;
    gravityx::rebuild_base_state_derived_fields(
        contingency_rating_data, {1}, base_rating_state);
    if (base_rating_state.sm_slack[1] <= 0.0) {
        throw std::runtime_error(
            "base-case RATE_A regression failed to create thermal slack");
    }
    gravityx::ContingencyContext branch_context;
    branch_context.outaged_branch = branch_contingency.component;
    branch_context.base_state = solve.state;
    const auto require_same_validation_numbers = [&] (
        const gravityx::ValidationReport& actual,
        const gravityx::ValidationReport& expected,
        const std::string& label) {
        require_near(
            actual.max_variable_bound_violation,
            expected.max_variable_bound_violation, 1e-12,
            label + " variable bounds");
        require_near(
            actual.max_pwl_sum_residual,
            expected.max_pwl_sum_residual, 1e-12,
            label + " PWL sum");
        require_near(
            actual.max_pwl_power_residual,
            expected.max_pwl_power_residual, 1e-12,
            label + " PWL power");
        require_near(
            actual.max_reference_angle_residual,
            expected.max_reference_angle_residual, 1e-12,
            label + " reference angle");
        require_near(
            actual.max_generator_residual,
            expected.max_generator_residual, 1e-12,
            label + " generator");
        require_near(
            actual.max_load_ramp_violation,
            expected.max_load_ramp_violation, 1e-12,
            label + " load ramp");
        require_near(
            actual.max_active_balance_residual,
            expected.max_active_balance_residual, 1e-12,
            label + " active balance");
        require_near(
            actual.max_reactive_balance_residual,
            expected.max_reactive_balance_residual, 1e-12,
            label + " reactive balance");
        require_near(
            actual.max_ohms_residual,
            expected.max_ohms_residual, 1e-12,
            label + " Ohm law");
        require_near(
            actual.max_angle_violation,
            expected.max_angle_violation, 1e-12,
            label + " angle");
        require_near(
            actual.max_flow_limit_violation,
            expected.max_flow_limit_violation, 1e-12,
            label + " flow limit");
        require_near(
            actual.max_residual, expected.max_residual, 1e-12,
            label + " overall residual");
    };
    const auto full_rebuilt_validation = gravityx::validate_state(
        data, gravityx::ModelMode::ContingencySoft,
        fast_result.solve.state, {1}, branch_context);
    const auto physical_rebuilt_validation =
        gravityx::validate_rebuilt_contingency_predictor(
            data, fast_result.solve.state, {1}, branch_context);
    const auto completed_rebuilt_validation =
        gravityx::validate_rebuilt_contingency_economic_and_ohms(
            data, fast_result.solve.state, {1}, branch_context,
            physical_rebuilt_validation);
    require_same_validation_numbers(
        completed_rebuilt_validation, full_rebuilt_validation,
        "split rebuilt validation");

    auto bad_pwl_state = fast_result.solve.state;
    if (bad_pwl_state.gen_lambda.empty()) {
        throw std::runtime_error(
            "split rebuilt validation regression has no generator PWL");
    }
    bad_pwl_state.gen_lambda[0] += 0.125;
    const auto bad_pwl_physical =
        gravityx::validate_rebuilt_contingency_predictor(
            data, bad_pwl_state, {1}, branch_context);
    const auto bad_pwl_completed =
        gravityx::validate_rebuilt_contingency_economic_and_ohms(
            data, bad_pwl_state, {1}, branch_context,
            bad_pwl_physical);
    const auto bad_pwl_full = gravityx::validate_state(
        data, gravityx::ModelMode::ContingencySoft,
        bad_pwl_state, {1}, branch_context);
    require_same_validation_numbers(
        bad_pwl_completed, bad_pwl_full,
        "split rebuilt PWL rejection");
    if (bad_pwl_completed.max_pwl_sum_residual <= 0.1) {
        throw std::runtime_error(
            "split rebuilt validation did not reject a PWL violation");
    }

    auto bad_ohms_state = fast_result.solve.state;
    bad_ohms_state.pf[1] += 0.25;
    const auto bad_ohms_physical =
        gravityx::validate_rebuilt_contingency_predictor(
            data, bad_ohms_state, {1}, branch_context);
    const auto bad_ohms_completed =
        gravityx::validate_rebuilt_contingency_economic_and_ohms(
            data, bad_ohms_state, {1}, branch_context,
            bad_ohms_physical);
    const auto bad_ohms_full = gravityx::validate_state(
        data, gravityx::ModelMode::ContingencySoft,
        bad_ohms_state, {1}, branch_context);
    require_same_validation_numbers(
        bad_ohms_completed, bad_ohms_full,
        "split rebuilt Ohm rejection");
    if (bad_ohms_completed.max_ohms_residual <= 0.2) {
        throw std::runtime_error(
            "split rebuilt validation did not reject an Ohm-law violation");
    }
    const auto lightweight_trial_validation =
        gravityx::validate_rebuilt_contingency_trial(
            data, fast_result.solve.state, {1}, branch_context);
    require_near(
        lightweight_trial_validation.max_variable_bound_violation,
        fast_result.validation.max_variable_bound_violation, 1e-12,
        "lightweight trial variable bounds");
    require_near(
        lightweight_trial_validation.max_active_balance_residual,
        fast_result.validation.max_active_balance_residual, 1e-12,
        "lightweight trial active balance");
    require_near(
        lightweight_trial_validation.max_reactive_balance_residual,
        fast_result.validation.max_reactive_balance_residual, 1e-12,
        "lightweight trial reactive balance");
    require_near(
        lightweight_trial_validation.max_angle_violation,
        fast_result.validation.max_angle_violation, 1e-12,
        "lightweight trial angle limits");
    require_near(
        lightweight_trial_validation.max_flow_limit_violation,
        fast_result.validation.max_flow_limit_violation, 1e-12,
        "lightweight trial flow limits");
    if (lightweight_trial_validation.max_pwl_sum_residual != 0.0 ||
        lightweight_trial_validation.max_pwl_power_residual != 0.0 ||
        lightweight_trial_validation.max_ohms_residual != 0.0) {
        throw std::runtime_error(
            "lightweight trial validation repeated rebuilt invariants");
    }
    const auto improving_cutoff_validation =
        gravityx::validate_rebuilt_contingency_trial_until_rejected(
            data, fast_result.solve.state, {1}, branch_context,
            lightweight_trial_validation.max_residual + 1e-6);
    require_same_validation_numbers(
        improving_cutoff_validation, lightweight_trial_validation,
        "lightweight cutoff complete improving report");
    auto rejected_cutoff_state = fast_result.solve.state;
    rejected_cutoff_state.vm[0] = data.buses[0].vmax + 0.2;
    const auto rejected_cutoff_validation =
        gravityx::validate_rebuilt_contingency_trial_until_rejected(
            data, rejected_cutoff_state, {1}, branch_context, 0.1);
    if (rejected_cutoff_validation.max_residual + 1e-10 < 0.1) {
        throw std::runtime_error(
            "lightweight cutoff validation did not prove rejection");
    }
    auto predictor_identity_state = fast_result.solve.state;
    predictor_identity_state.vm[0] = data.buses[0].vmax + 0.2;
    gravityx::rebuild_contingency_state_derived_fields(
        data, solve.state, {1}, branch_contingency,
        predictor_identity_state);
    const auto identity_free_trial_validation =
        gravityx::validate_rebuilt_contingency_trial(
            data, predictor_identity_state, {1}, branch_context);
    const auto predictor_identity_validation =
        gravityx::validate_rebuilt_contingency_predictor(
            data, predictor_identity_state, {1}, branch_context);
    require_near(
        predictor_identity_validation.max_residual,
        identity_free_trial_validation.max_residual, 1e-12,
        "predictor validation residual identity capture");
    if (!identity_free_trial_validation.worst_identity.empty() ||
        predictor_identity_validation.worst_identity.empty() ||
        predictor_identity_validation.worst_category.empty()) {
        throw std::runtime_error(
            "lightweight predictor validation did not preserve residual "
            "routing identity");
    }
    auto active_repair_reference = fast_result.solve.state;
    active_repair_reference.va[1] += 0.35;
    gravityx::rebuild_contingency_state_derived_fields(
        data, solve.state, {1}, branch_contingency,
        active_repair_reference);
    const auto active_repair =
        gravityx::solve_linearized_active_feasibility_repair(
            data, active_repair_reference, {1}, branch_context,
            0.49, 0.5, 5.0, 0.1);
    if (!active_repair.success || active_repair.row_count < 2 ||
        active_repair.branch_security_row_count <= 0 ||
        !active_repair.finite_solution_values ||
        active_repair.maximum_column_violation > 1e-7 ||
        active_repair.maximum_linearized_violation > 1e-7 ||
        active_repair.maximum_angle_change > 0.5 + 1e-9) {
        throw std::runtime_error(
            "active feasibility repair regression failed: " +
            active_repair.status);
    }
    const auto screened_active_repair =
        gravityx::solve_linearized_active_feasibility_repair(
            data, active_repair_reference, {1}, branch_context,
            0.49, 0.5, 5.0, 0.1, true, true);
    if (!screened_active_repair.success ||
        !screened_active_repair.include_reactive ||
        !screened_active_repair.current_security_rows_only ||
        screened_active_repair.simplex_strategy != 4 ||
        !screened_active_repair.finite_solution_values ||
        screened_active_repair.maximum_column_violation > 1e-7 ||
        screened_active_repair.maximum_linearized_violation > 1e-7 ||
        screened_active_repair.maximum_voltage_change > 0.1 + 1e-9 ||
        screened_active_repair.branch_security_row_count >
            active_repair.branch_security_row_count) {
        throw std::runtime_error(
            "screened active feasibility repair regression failed: " +
            screened_active_repair.status);
    }
    auto voltage_bound_reference = fast_result.solve.state;
    voltage_bound_reference.vm[1] = data.buses[1].vmin - 0.05;
    gravityx::rebuild_contingency_state_derived_fields(
        data, solve.state, {1}, branch_contingency,
        voltage_bound_reference);
    const auto voltage_bound_repair =
        gravityx::solve_linearized_active_feasibility_repair(
            data, voltage_bound_reference, {1}, branch_context,
            0.49, 0.5, 5.0, 0.1, true, true);
    if (!voltage_bound_repair.success ||
        voltage_bound_repair.state.vm[1] < data.buses[1].vmin - 1e-9 ||
        voltage_bound_repair.maximum_voltage_change < 0.05 - 1e-9) {
        throw std::runtime_error(
            "active feasibility repair did not enforce an out-of-bound "
            "reference voltage: " + voltage_bound_repair.status);
    }
    auto active_repair_state = active_repair.state;
    gravityx::rebuild_contingency_state_derived_fields(
        data, solve.state, {1}, branch_contingency,
        active_repair_state);
    const auto active_repair_validation = gravityx::validate_state(
        data, gravityx::ModelMode::ContingencySoft,
        active_repair_state, {1}, branch_context);
    if (active_repair_validation.max_generator_residual > 1e-5 ||
        active_repair_validation.max_load_ramp_violation > 1e-5 ||
        active_repair_validation.max_reference_angle_residual > 1e-5) {
        throw std::runtime_error(
            "active feasibility repair source-bound validation failed with residual " +
            std::to_string(active_repair_validation.max_residual));
    }
    const auto full_contingency_seed = gravityx::solve_linearized_ac_seed(
        data, solve.state, {1}, 0.49, branch_context);
    const auto balance_only_contingency_seed =
        gravityx::solve_linearized_ac_seed(
            data, solve.state, {1}, 0.49, branch_context,
            true, false, 60.0, true, true);
    const auto subset_security_contingency_seed =
        gravityx::solve_linearized_ac_seed(
            data, solve.state, {1}, 0.49, branch_context,
            true, false, 60.0, true, true, {1});
    if (!full_contingency_seed.success ||
        !balance_only_contingency_seed.success ||
        !subset_security_contingency_seed.success ||
        !full_contingency_seed.model_preflight_passed ||
        !full_contingency_seed.model_construction_success ||
        full_contingency_seed.add_vars_status != 0 ||
        full_contingency_seed.change_cols_cost_status != 0 ||
        full_contingency_seed.add_rows_status != 0 ||
        !balance_only_contingency_seed.projected_balance_slack ||
        !balance_only_contingency_seed.branch_security_rows_omitted ||
        !balance_only_contingency_seed.feasibility_only ||
        subset_security_contingency_seed.branch_security_subset_count != 1 ||
        std::abs(balance_only_contingency_seed.ipm_optimality_tolerance -
                 1e-4) > 1e-12 ||
        std::abs(full_contingency_seed.ipm_optimality_tolerance -
                 1e-8) > 1e-14 ||
        balance_only_contingency_seed.row_count >=
            subset_security_contingency_seed.row_count ||
        subset_security_contingency_seed.row_count >
            full_contingency_seed.row_count) {
        throw std::runtime_error(
            "balance-only contingency Phase-I regression failed");
    }
    auto nonfinite_linear_reference = solve.state;
    nonfinite_linear_reference.vm[0] =
        std::numeric_limits<double>::quiet_NaN();
    const auto nonfinite_linear_seed = gravityx::solve_linearized_ac_seed(
        data, nonfinite_linear_reference, {1}, 0.49, branch_context);
    if (nonfinite_linear_seed.success ||
        nonfinite_linear_seed.model_preflight_passed ||
        nonfinite_linear_seed.model_load_failure_call != "preflight" ||
        nonfinite_linear_seed.model_preflight_failure.empty()) {
        throw std::runtime_error(
            "linearized model preflight did not reject non-finite input");
    }
    const auto contingency_rating_validation = gravityx::validate_state(
        contingency_rating_data, gravityx::ModelMode::ContingencySoft,
        fast_result.solve.state, {1}, branch_context);
    if (contingency_rating_validation.max_residual > 1e-5) {
        throw std::runtime_error(
            "contingency RATE_C regression failed with residual " +
            std::to_string(contingency_rating_validation.max_residual));
    }
    auto linear_reference = solve.state;
    gravityx::ValidationReport linear_validation;
    bool linear_contingency_feasible = false;
    for (int round = 0; round < 3; ++round) {
        const auto contingency_seed = gravityx::solve_linearized_ac_seed(
            data, linear_reference, {1}, 0.49, branch_context);
        if (!contingency_seed.success) {
            throw std::runtime_error(
                "linearized contingency regression failed: "
                + contingency_seed.status);
        }
        linear_reference = contingency_seed.state;
        gravityx::rebuild_contingency_state_derived_fields(
            data, solve.state, {1}, branch_contingency, linear_reference);
        linear_validation = gravityx::validate_state(
            data, gravityx::ModelMode::ContingencySoft,
            linear_reference, {1}, branch_context);
        if (linear_validation.max_residual <= 1e-5) {
            linear_contingency_feasible = true;
            break;
        }
    }
    if (!linear_contingency_feasible) {
        throw std::runtime_error(
            "linearized contingency validation regression failed with residual "
            + std::to_string(linear_validation.max_residual));
    }
    auto nonfinite_state = solve.state;
    nonfinite_state.vm[0] = std::numeric_limits<double>::quiet_NaN();
    bool nonfinite_rejected = false;
    try {
        static_cast<void>(gravityx::validate_state(
            data, gravityx::ModelMode::BaseSoft, nonfinite_state, {1}));
    } catch (const std::runtime_error&) {
        nonfinite_rejected = true;
    }
    if (!nonfinite_rejected) {
        throw std::runtime_error("nonfinite-state validation regression failed");
    }
    std::cout << "parallel-circuit regression passed with max residual "
              << validation.max_residual << '\n';
    return 0;
}

nlohmann::json state_summary(const gravityx::SolveResult& result) {
    const auto& state = result.state;
    nlohmann::json summary = {
        {"status", result.status},
        {"objective", result.objective},
        {"wall_seconds", result.wall_seconds},
    };
    if (!state.commitment.empty()) {
        int committed = 0;
        int fractional = 0;
        for (double value : state.commitment) {
            committed += value >= 0.5 ? 1 : 0;
            fractional += value > 1e-6 && value < 1.0 - 1e-6 ? 1 : 0;
        }
        summary["committed_at_half"] = committed;
        summary["fractional_commitments"] = fractional;
    }
    if (!state.p_delta.empty()) {
        summary["max_p_delta"] = *std::max_element(state.p_delta.begin(), state.p_delta.end());
        summary["max_q_delta"] = *std::max_element(state.q_delta.begin(), state.q_delta.end());
    }
    summary["max_sm_slack"] = state.sm_slack.empty()
        ? 0.0
        : *std::max_element(state.sm_slack.begin(), state.sm_slack.end());
    return summary;
}

int run_smoke() {
    Model model("Gravity framework smoke test");
    var<double> x("x", -10.0, 10.0);
    x.initialize_all(0.0);
    model.add_var(x);
    model.set_objective(min((x - 3.0) * (x - 3.0)));

    solver nlp(model, ipopt);
    const auto status = nlp.run(0, false, 1e-9, "mumps");
    const double value = x.getvalue();
    std::cout << "status=" << status << " x=" << value
              << " objective=" << model._obj_val << '\n';
    if (!std::isfinite(value) || std::abs(value - 3.0) > 1e-5) {
        throw std::runtime_error("Gravity/Ipopt smoke solution is incorrect");
    }
    return 0;
}

int run_inspect(const std::string& path) {
    reject_onedrive(path);
    const auto data = gravityx::CaseData::load(path);
    int transformers = 0;
    int source_active_branches = 0;
    int initially_on = 0;
    int startup_eligible = 0;
    int shutdown_eligible = 0;
    for (const auto& branch : data.branches) {
        transformers += branch.transformer ? 1 : 0;
        source_active_branches += branch.status;
    }
    for (const auto& gen : data.generators) {
        initially_on += gen.status_prev;
        startup_eligible += gen.status_prev == 0 && gen.suqual == 1 ? 1 : 0;
        shutdown_eligible += gen.status_prev == 1 && gen.sdqual == 1 ? 1 : 0;
    }
    const nlohmann::json result = {
        {"name", data.name},
        {"buses", data.buses.size()},
        {"generators", data.generators.size()},
        {"loads", data.loads.size()},
        {"shunts", data.shunts.size()},
        {"branches", data.branches.size()},
        {"source_active_branches", source_active_branches},
        {"source_inactive_branches",
         static_cast<int>(data.branches.size()) - source_active_branches},
        {"contingencies", data.contingencies.size()},
        {"transformers", transformers},
        {"lines", static_cast<int>(data.branches.size()) - transformers},
        {"initially_on", initially_on},
        {"startup_eligible", startup_eligible},
        {"shutdown_eligible", shutdown_eligible},
    };
    std::cout << result.dump(2) << '\n';
    return 0;
}

int run_ac_model(const std::string& command, const std::string& path, int print_level) {
    reject_onedrive(path);
    const auto data = gravityx::CaseData::load(path);
    const auto mode = command == "solve-base"
        ? gravityx::ModelMode::BaseSoft
        : gravityx::ModelMode::UnitCommitmentRelaxation;
    gravityx::AcModel model(data, mode);
    const auto result = model.solve(print_level, 1e-6);
    auto summary = state_summary(result);
    summary["validation"] = gravityx::validate_state(data, mode, result.state).to_json();
    if (data.buses.size() <= 50) {
        summary["state"] = {
            {"vm", result.state.vm},
            {"va", result.state.va},
            {"pg", result.state.pg},
            {"qg", result.state.qg},
            {"demand_factor", result.state.demand_factor},
            {"commitment", result.state.commitment},
            {"gen_lambda", result.state.gen_lambda},
            {"load_lambda", result.state.load_lambda},
        };
    }
    std::cout << summary.dump(2) << '\n';
    return result.status == 0 ? 0 : 1;
}

int run_ibr(const std::string& path, int print_level) {
    reject_onedrive(path);
    const auto data = gravityx::CaseData::load(path);
    gravityx::IbrOptions options;
    options.print_level = print_level;
    const auto result = gravityx::run_iterative_batch_rounding(data, options);
    std::cout << result.to_json(data.buses.size() <= 50).dump(2) << '\n';
    return result.success ? 0 : 1;
}

int run_ibr_json(
    const std::string& path,
    const std::string& output_path,
    int print_level,
    bool source_status_only = false) {
    reject_onedrive(path);
    reject_onedrive(output_path);
    const auto data = gravityx::CaseData::load(path);
    gravityx::IbrOptions options;
    options.print_level = print_level;
    options.source_status_only = source_status_only;
    const auto result = gravityx::run_iterative_batch_rounding(data, options);
    const auto json = result.to_json(true);
    write_json_file(output_path, json);
    std::cout << nlohmann::json({
        {"output", output_path},
        {"success", result.success},
        {"candidate_accepted", result.candidate_accepted},
        {"wall_seconds", result.wall_seconds},
        {"selected_commitment", result.commitment},
    }).dump(2) << '\n';
    return result.success ? 0 : 1;
}

int run_validated_source_base_json(
    const std::string& path,
    const std::string& output_path,
    bool allow_exact_fallback = true,
    bool allow_large_base_newton_restart = true) {
    reject_onedrive(path);
    reject_onedrive(output_path);
    const auto command_start = std::chrono::steady_clock::now();
    const char* base_log_value = std::getenv("GRAVITYX_BASE_LOG");
    const bool base_log_enabled = base_log_value != nullptr &&
        std::string(base_log_value) != "0";
    const auto log_base_phase = [&](const std::string& phase,
                                    const nlohmann::json& detail) {
        if (!base_log_enabled) {
            return;
        }
        auto message = detail;
        message["phase"] = phase;
        message["elapsed_seconds"] = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - command_start).count();
        std::cerr << "GRAVITYX_BASE_PHASE " << message.dump() << '\n';
        std::cerr.flush();
    };
    const auto data = gravityx::CaseData::load(path);
    log_base_phase("case_loaded", {
        {"buses", data.buses.size()},
        {"generators", data.generators.size()},
        {"loads", data.loads.size()},
        {"branches", data.branches.size()},
        {"contingencies", data.contingencies.size()},
    });
    std::vector<int> commitment;
    commitment.reserve(data.generators.size());
    for (const auto& generator : data.generators) {
        commitment.push_back(generator.status_prev);
    }
    const auto source = gravityx::build_validated_source_base(data, commitment);
    log_base_phase("source_candidate", {
        {"feasible", source.feasible},
        {"max_residual", source.validation.max_residual},
        {"worst_category", source.validation.worst_category},
        {"wall_seconds", source.wall_seconds},
    });
    gravityx::SolveResult selected_solve = source.solve;
    gravityx::ValidationReport selected_validation = source.validation;
    bool success = source.feasible;
    double wall_seconds = source.wall_seconds;
    nlohmann::json repair_json = nlohmann::json::array();
    nlohmann::json linearized_repair_json = nlohmann::json::array();
    nlohmann::json exact_repair_json = nullptr;
    bool base_optimization_performed = false;
    std::string base_method = "independently_validated_source_operating_point";
    if (!success) {
        gravityx::FastContingencyPowerFlow repair(
            data, source.solve.state, commitment);
        auto repaired = repair.solve_base();
        log_base_phase("source_voltage_repair", {
            {"feasible", repaired.feasible},
            {"converged", repaired.converged},
            {"max_residual", repaired.validation.max_residual},
            {"worst_category", repaired.validation.worst_category},
            {"wall_seconds", repaired.wall_seconds},
        });
        auto source_repair_json = repaired.to_json();
        source_repair_json["start"] = "source_voltage";
        repair_json.push_back(std::move(source_repair_json));
        wall_seconds += repaired.wall_seconds;
        if (repaired.validation.max_residual < selected_validation.max_residual) {
            selected_solve = repaired.solve;
            selected_validation = repaired.validation;
        }
        success = repaired.feasible;
        if (!success && !(data.buses.size() >= 16000 &&
                          selected_validation.max_residual <= 0.25)) {
            auto flat_state = source.solve.state;
            for (int i = 0; i < static_cast<int>(data.buses.size()); ++i) {
                flat_state.vm[i] = std::clamp(
                    1.0, data.buses[i].vmin, data.buses[i].vmax);
                flat_state.va[i] = 0.0;
            }
            gravityx::FastContingencyPowerFlow flat_repair(
                data, flat_state, commitment);
            const auto flat = flat_repair.solve_base();
            log_base_phase("flat_voltage_repair", {
                {"feasible", flat.feasible},
                {"converged", flat.converged},
                {"max_residual", flat.validation.max_residual},
                {"worst_category", flat.validation.worst_category},
                {"wall_seconds", flat.wall_seconds},
            });
            auto flat_json = flat.to_json();
            flat_json["start"] = "flat_voltage";
            repair_json.push_back(std::move(flat_json));
            wall_seconds += flat.wall_seconds;
            if (flat.validation.max_residual < selected_validation.max_residual) {
                selected_solve = flat.solve;
                selected_validation = flat.validation;
            }
            success = flat.feasible;
        }
        if (!success && !(data.buses.size() >= 16000 &&
                          selected_validation.max_residual <= 0.25)) {
            const std::vector<double> load_fractions{
                0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.5, 0.4, 0.25, 0.0};
            auto continuation_state = selected_solve.state;
            double continuation_residual = selected_validation.max_residual;
            for (double fraction : load_fractions) {
                auto seed = continuation_state;
                double total_load = 0.0;
                for (int i = 0; i < static_cast<int>(data.loads.size()); ++i) {
                    const auto& load = data.loads[i];
                    double lower = load.tmin;
                    if (std::abs(load.pd_nominal) > 1e-12) {
                        lower = std::max(
                            lower,
                            (load.pd_prev - load.prdmax * data.delta_r)
                                / load.pd_nominal);
                    }
                    seed.demand_factor[i] = lower + fraction
                        * (source.solve.state.demand_factor[i] - lower);
                    total_load += load.pd_nominal * seed.demand_factor[i];
                }

                std::vector<double> lower(data.generators.size(), 0.0);
                std::vector<double> upper(data.generators.size(), 0.0);
                double lower_sum = 0.0;
                double upper_sum = 0.0;
                double current_sum = 0.0;
                for (int i = 0; i < static_cast<int>(data.generators.size()); ++i) {
                    const auto& gen = data.generators[i];
                    if (commitment[i] == 0) {
                        seed.pg[i] = 0.0;
                        seed.qg[i] = 0.0;
                        continue;
                    }
                    const double previous = gen.status_prev == 0
                        ? gen.pmin : gen.pg_prev;
                    lower[i] = std::max(
                        gen.pmin, previous - data.delta_r * gen.prdmax);
                    upper[i] = std::min(
                        gen.pmax, previous + data.delta_r * gen.prumax);
                    seed.pg[i] = std::clamp(seed.pg[i], lower[i], upper[i]);
                    lower_sum += lower[i];
                    upper_sum += upper[i];
                    current_sum += seed.pg[i];
                }
                const double target = std::clamp(
                    total_load * 1.02, lower_sum, upper_sum);
                for (int pass = 0; pass < 8; ++pass) {
                    const double difference = target - current_sum;
                    if (std::abs(difference) <= 1e-9) {
                        break;
                    }
                    double room = 0.0;
                    for (int i = 0; i < static_cast<int>(data.generators.size()); ++i) {
                        if (commitment[i] == 0) {
                            continue;
                        }
                        room += difference > 0.0
                            ? std::max(0.0, upper[i] - seed.pg[i])
                            : std::max(0.0, seed.pg[i] - lower[i]);
                    }
                    if (room <= 1e-12) {
                        break;
                    }
                    current_sum = 0.0;
                    for (int i = 0; i < static_cast<int>(data.generators.size()); ++i) {
                        if (commitment[i] == 0) {
                            continue;
                        }
                        const double individual_room = difference > 0.0
                            ? std::max(0.0, upper[i] - seed.pg[i])
                            : std::max(0.0, seed.pg[i] - lower[i]);
                        seed.pg[i] = std::clamp(
                            seed.pg[i] + difference * individual_room / room,
                            lower[i], upper[i]);
                        current_sum += seed.pg[i];
                    }
                }
                gravityx::FastContingencyPowerFlow reduced_repair(
                    data, seed, commitment);
                const auto reduced = reduced_repair.solve_base();
                const double prior_continuation_residual =
                    continuation_residual;
                log_base_phase("reduced_load_repair", {
                    {"load_fraction_from_lower_to_source", fraction},
                    {"feasible", reduced.feasible},
                    {"converged", reduced.converged},
                    {"max_residual", reduced.validation.max_residual},
                    {"worst_category", reduced.validation.worst_category},
                    {"wall_seconds", reduced.wall_seconds},
                });
                auto reduced_json = reduced.to_json();
                reduced_json["start"] = "reduced_load_continuation";
                reduced_json["load_fraction_from_lower_to_source"] = fraction;
                reduced_json["target_total_generation"] = target;
                repair_json.push_back(std::move(reduced_json));
                wall_seconds += reduced.wall_seconds;
                if (reduced.converged ||
                    reduced.validation.max_residual < continuation_residual) {
                    continuation_state = reduced.solve.state;
                    continuation_residual = reduced.validation.max_residual;
                }
                if (reduced.validation.max_residual <
                    selected_validation.max_residual) {
                    selected_solve = reduced.solve;
                    selected_validation = reduced.validation;
                }
                if (reduced.feasible) {
                    selected_solve = reduced.solve;
                    selected_validation = reduced.validation;
                    success = true;
                    break;
                }
                // The 19k-bus source cases reach their best continuation
                // center near the 70% load point, then spend roughly another
                // minute walking toward lower loads while the independently
                // measured residual worsens.  Stop at the first non-improving
                // point.  This changes only seed selection; no candidate is
                // accepted without the unchanged exact validator.
                if (data.buses.size() >= 16000 &&
                    reduced.validation.max_residual + 1e-8 >=
                        prior_continuation_residual) {
                    log_base_phase("reduced_load_stagnation_stop", {
                        {"load_fraction_from_lower_to_source", fraction},
                        {"prior_residual", prior_continuation_residual},
                        {"current_residual",
                         reduced.validation.max_residual},
                    });
                    break;
                }
                // On the 8k-bus cases this continuation is only a way to find
                // a useful center for the linearized feasibility LP.  Once a
                // sub-0.1 p.u. candidate exists, the lower-load points were
                // measured to move away from feasibility while consuming most
                // of twenty seconds.  The candidate is not accepted here: it
                // still goes through the linearized solve, sparse Newton, and
                // the unchanged independent validator below.
                if (data.buses.size() >= 8000 &&
                    selected_validation.max_residual <= 0.1) {
                    break;
                }
            }
        }
        base_method = "validated_sparse_newton_source_point_repair";
    }
    if (!success && data.buses.size() >= 16000 && !allow_exact_fallback) {
        std::vector<int> dynamic_security_branches;
        std::vector<unsigned char> dynamic_security_selected(
            data.branches.size(), 0);
        const auto collect_violated_base_security_branches =
            [&](const gravityx::AcState& state) {
                nlohmann::json added = nlohmann::json::array();
                if (state.pf.size() != data.branches.size() ||
                    state.qf.size() != data.branches.size() ||
                    state.pt.size() != data.branches.size() ||
                    state.qt.size() != data.branches.size() ||
                    state.sm_slack.size() != data.branches.size()) {
                    return added;
                }
                constexpr double kSecurityCollectionTolerance = 1e-5;
                for (std::size_t i = 0; i < data.branches.size(); ++i) {
                    const auto& branch = data.branches[i];
                    if (branch.status == 0) {
                        continue;
                    }
                    const double rating = branch.rate_a;
                    const double box_violation = std::max({
                        std::abs(state.pf[i]) - rating,
                        std::abs(state.qf[i]) - rating,
                        std::abs(state.pt[i]) - rating,
                        std::abs(state.qt[i]) - rating,
                    });
                    const double source_delta =
                        data.buses[branch.from].va_start -
                        data.buses[branch.to].va_start;
                    double angle_violation = 0.0;
                    if (source_delta >= branch.angmin &&
                        source_delta <= branch.angmax) {
                        const double angle =
                            state.va[branch.from] - state.va[branch.to];
                        angle_violation = std::max(
                            angle - branch.angmax,
                            branch.angmin - angle);
                    }
                    const double slack = state.sm_slack[i];
                    const double from_scale = branch.transformer
                        ? 1.0 + slack : state.vm[branch.from] + slack;
                    const double to_scale = branch.transformer
                        ? 1.0 + slack : state.vm[branch.to] + slack;
                    const double apparent_violation = std::max(
                        state.pf[i] * state.pf[i] +
                            state.qf[i] * state.qf[i] -
                            rating * rating * from_scale * from_scale,
                        state.pt[i] * state.pt[i] +
                            state.qt[i] * state.qt[i] -
                            rating * rating * to_scale * to_scale);
                    if (std::max({box_violation, angle_violation,
                                  apparent_violation}) <=
                            kSecurityCollectionTolerance ||
                        dynamic_security_selected[i]) {
                        continue;
                    }
                    dynamic_security_selected[i] = 1;
                    dynamic_security_branches.push_back(
                        static_cast<int>(i));
                    added.push_back({
                        {"component_position", static_cast<int>(i)},
                        {"source_key", branch.source_key},
                        {"box_violation", std::max(0.0, box_violation)},
                        {"angle_violation", std::max(0.0, angle_violation)},
                        {"apparent_flow_violation",
                         std::max(0.0, apparent_violation)},
                    });
                }
                std::sort(dynamic_security_branches.begin(),
                          dynamic_security_branches.end());
                return added;
            };

        auto linear_reference = selected_solve.state;
        double linear_reference_residual = selected_validation.max_residual;
        // Flow variables retain explicit component bounds independently of
        // the apparent-flow soft slack.  Seed the restricted Phase I with the
        // exact branches that violate those bounds (or angle/apparent-flow
        // constraints), then add newly exposed rows after nonlinear repair.
        const nlohmann::json initial_dynamic_security_branches =
            collect_violated_base_security_branches(linear_reference);
        log_base_phase("dynamic_base_phase_one_initialized", {
            {"initial_security_branch_count",
             dynamic_security_branches.size()},
            {"reference_residual", linear_reference_residual},
        });
        constexpr int kMaximumDynamicBaseRounds = 16;
        for (int round = 1; round <= kMaximumDynamicBaseRounds; ++round) {
            const double phase_one_time_limit_seconds =
                dynamic_security_branches.empty() ? 60.0 : 90.0;
            // A zero-objective targeted Phase-I model may choose any feasible
            // point in its trust box.  When only a small security violation
            // remains, the original 0.01/0.006 box is much larger than the
            // required correction and can create a large nonlinear reactive
            // mismatch.  Contract the box with the independently measured
            // residual while retaining the original radii for coarse repair.
            const bool use_adaptive_security_trust =
                !dynamic_security_branches.empty() &&
                linear_reference_residual <= 0.02;
            const double voltage_trust_radius =
                use_adaptive_security_trust
                ? std::clamp(
                    2.0 * linear_reference_residual, 0.001, 0.01)
                : -1.0;
            const double angle_trust_radius =
                use_adaptive_security_trust
                ? std::clamp(
                    0.25 * linear_reference_residual, 0.0001, 0.006)
                : -1.0;
            log_base_phase("dynamic_base_phase_one_start", {
                {"round", round},
                {"reference_residual", linear_reference_residual},
                {"security_branch_count", dynamic_security_branches.size()},
                {"time_limit_seconds", phase_one_time_limit_seconds},
                {"adaptive_security_trust", use_adaptive_security_trust},
                {"voltage_trust_radius", voltage_trust_radius},
                {"angle_trust_radius", angle_trust_radius},
            });
            bool lightweight = true;
            auto linear = gravityx::solve_linearized_ac_seed(
                data, linear_reference, commitment, 0.49, std::nullopt,
                true, lightweight, phase_one_time_limit_seconds, true, true,
                dynamic_security_branches, false,
                voltage_trust_radius, angle_trust_radius);
            if (!linear.success && linear.status == "Infeasible") {
                // A local voltage/angle trust box may not contain a balance
                // repair when the source point is several p.u. infeasible.
                // Retry the same reduced Phase-I rows without that box.
                lightweight = false;
                linear = gravityx::solve_linearized_ac_seed(
                    data, linear_reference, commitment, 0.49, std::nullopt,
                    true, lightweight, phase_one_time_limit_seconds,
                    true, true,
                    dynamic_security_branches, false,
                    voltage_trust_radius, angle_trust_radius);
            }
            wall_seconds += linear.wall_seconds;
            auto round_json = linear.to_json(false);
            round_json["round"] = round;
            round_json["dynamic_base_phase_one"] = true;
            round_json["lightweight_large_base_seed"] = lightweight;
            round_json["initial_dynamic_security_branches"] =
                initial_dynamic_security_branches;
            round_json["active_dynamic_security_branch_positions"] =
                dynamic_security_branches;
            log_base_phase("dynamic_base_phase_one_complete", {
                {"round", round},
                {"success", linear.success},
                {"status", linear.status},
                {"rows", linear.row_count},
                {"columns", linear.column_count},
                {"nonzeros", linear.nonzero_count},
                {"iterations", linear.iterations},
                {"wall_seconds", linear.wall_seconds},
                {"lightweight_large_base_seed", lightweight},
                {"security_branch_count", dynamic_security_branches.size()},
            });
            if (!linear.success) {
                round_json["nonlinear_repair"] = nullptr;
                linearized_repair_json.push_back(std::move(round_json));
                break;
            }

            auto linear_state = linear.state;
            const double linear_objective =
                gravityx::rebuild_base_state_derived_fields(
                    data, commitment, linear_state);
            const auto linear_validation = gravityx::validate_state(
                data, gravityx::ModelMode::BaseSoft,
                linear_state, commitment);
            round_json["linear_objective"] = linear_objective;
            round_json["linear_validation"] = linear_validation.to_json();
            const auto security_branches_from_full_linear_point =
                collect_violated_base_security_branches(linear_state);
            round_json["new_dynamic_security_branches_from_full_linear_point"] =
                security_branches_from_full_linear_point;
            round_json["active_dynamic_security_branch_positions_after_full_linear_point"] =
                dynamic_security_branches;
            if (linear_validation.max_residual <
                selected_validation.max_residual) {
                selected_solve.status = 0;
                selected_solve.objective = linear_objective;
                selected_solve.wall_seconds = linear.wall_seconds;
                selected_solve.iterations = linear.iterations;
                selected_solve.state = linear_state;
                selected_validation = linear_validation;
            }
            if (linear_validation.max_residual <= 1e-5) {
                selected_solve.status = 0;
                selected_solve.objective = linear_objective;
                selected_solve.wall_seconds = linear.wall_seconds;
                selected_solve.iterations = linear.iterations;
                selected_solve.state = std::move(linear_state);
                selected_validation = linear_validation;
                success = true;
                base_method = "highs_dynamic_security_base_phase_one";
                round_json["nonlinear_repair"] = nullptr;
                linearized_repair_json.push_back(std::move(round_json));
                break;
            }

            auto repair_seed = linear_state;
            auto repair_seed_validation = linear_validation;
            double repair_seed_objective = linear_objective;
            double repair_seed_fraction = 1.0;
            const bool full_step_is_security_feasible =
                linear_validation.max_variable_bound_violation <= 1e-5 &&
                linear_validation.max_angle_violation <= 1e-5 &&
                linear_validation.max_flow_limit_violation <= 1e-5;
            nlohmann::json blend_candidates = nlohmann::json::array({{
                {"fraction", 1.0},
                {"objective", linear_objective},
                {"validation", linear_validation.to_json()},
            }});
            for (const double fraction : {0.75, 0.5, 0.25}) {
                auto blended = linear_reference;
                const auto interpolate = [fraction](
                    const std::vector<double>& reference_values,
                    const std::vector<double>& linear_values,
                    std::vector<double>& output_values) {
                    output_values.resize(reference_values.size());
                    for (std::size_t i = 0; i < reference_values.size(); ++i) {
                        output_values[i] = reference_values[i] + fraction *
                            (linear_values[i] - reference_values[i]);
                    }
                };
                interpolate(
                    linear_reference.vm, linear_state.vm, blended.vm);
                interpolate(
                    linear_reference.va, linear_state.va, blended.va);
                interpolate(
                    linear_reference.pg, linear_state.pg, blended.pg);
                interpolate(
                    linear_reference.qg, linear_state.qg, blended.qg);
                interpolate(
                    linear_reference.demand_factor,
                    linear_state.demand_factor,
                    blended.demand_factor);
                const double blended_objective =
                    gravityx::rebuild_base_state_derived_fields(
                        data, commitment, blended);
                const auto blended_validation = gravityx::validate_state(
                    data, gravityx::ModelMode::BaseSoft,
                    blended, commitment);
                blend_candidates.push_back({
                    {"fraction", fraction},
                    {"objective", blended_objective},
                    {"validation", blended_validation.to_json()},
                });
                if (!full_step_is_security_feasible &&
                    blended_validation.max_residual + 1e-12 <
                     repair_seed_validation.max_residual) {
                    repair_seed = std::move(blended);
                    repair_seed_validation = blended_validation;
                    repair_seed_objective = blended_objective;
                    repair_seed_fraction = fraction;
                }
            }
            round_json["blend_candidates"] = std::move(blend_candidates);
            round_json["selected_blend_fraction"] = repair_seed_fraction;
            round_json["full_step_is_security_feasible"] =
                full_step_is_security_feasible;
            round_json["selected_blend_validation"] =
                repair_seed_validation.to_json();
            log_base_phase("dynamic_base_blend_selected", {
                {"round", round},
                {"fraction", repair_seed_fraction},
                {"max_residual", repair_seed_validation.max_residual},
                {"worst_category", repair_seed_validation.worst_category},
                {"worst_identity", repair_seed_validation.worst_identity},
            });
            if (repair_seed_validation.max_residual <
                selected_validation.max_residual) {
                selected_solve.status = 0;
                selected_solve.objective = repair_seed_objective;
                selected_solve.wall_seconds = linear.wall_seconds;
                selected_solve.iterations = linear.iterations;
                selected_solve.state = repair_seed;
                selected_validation = repair_seed_validation;
            }
            if (repair_seed_validation.max_residual <= 1e-5) {
                success = true;
                base_method =
                    "highs_dynamic_security_base_phase_one_blended";
                round_json["nonlinear_repair"] = nullptr;
                linearized_repair_json.push_back(std::move(round_json));
                break;
            }

            gravityx::FastContingencyPowerFlow nonlinear_repair(
                data, repair_seed, commitment);
            auto nonlinear = nonlinear_repair.solve_base();
            wall_seconds += nonlinear.wall_seconds;
            round_json["nonlinear_repair"] = nonlinear.to_json();
            const auto new_security_branches =
                collect_violated_base_security_branches(
                    nonlinear.solve.state);
            round_json["new_dynamic_security_branches"] =
                new_security_branches;
            round_json["active_dynamic_security_branch_positions_after_repair"] =
                dynamic_security_branches;
            linearized_repair_json.push_back(std::move(round_json));
            if (nonlinear.validation.max_residual <
                selected_validation.max_residual) {
                selected_solve = nonlinear.solve;
                selected_validation = nonlinear.validation;
            }
            if (nonlinear.feasible) {
                selected_solve = nonlinear.solve;
                selected_validation = nonlinear.validation;
                success = true;
                base_method =
                    "highs_dynamic_security_base_phase_one_plus_validated_sparse_newton";
                break;
            }
            if (nonlinear.converged ||
                nonlinear.validation.max_residual <
                    linear_reference_residual) {
                linear_reference = nonlinear.solve.state;
                linear_reference_residual =
                    nonlinear.validation.max_residual;
            } else {
                if (use_adaptive_security_trust &&
                    linear_reference_residual <= 0.01) {
                    log_base_phase("dynamic_base_phase_one_stagnation_stop", {
                        {"round", round},
                        {"reference_residual", linear_reference_residual},
                        {"nonlinear_residual",
                         nonlinear.validation.max_residual},
                    });
                    break;
                }
                linear_reference = repair_seed;
                linear_reference_residual =
                    repair_seed_validation.max_residual;
            }
        }
    }
    if (!success && data.buses.size() < 16000) {
        auto linear_reference = selected_solve.state;
        double linear_reference_residual = selected_validation.max_residual;
        for (int round = 1; round <= 4; ++round) {
            const bool lightweight_large_base_seed =
                data.buses.size() >= 16000 ||
                (data.buses.size() >= 8000 && round >= 2 &&
                 !allow_large_base_newton_restart);
            log_base_phase("linearized_seed_start", {
                {"round", round},
                {"reference_residual", linear_reference_residual},
                {"lightweight_large_base_seed", lightweight_large_base_seed},
            });
            const auto linear = gravityx::solve_linearized_ac_seed(
                data, linear_reference, commitment, 0.49, std::nullopt,
                false, lightweight_large_base_seed);
            log_base_phase("linearized_seed_complete", {
                {"round", round},
                {"success", linear.success},
                {"status", linear.status},
                {"rows", linear.row_count},
                {"columns", linear.column_count},
                {"nonzeros", linear.nonzero_count},
                {"iterations", linear.iterations},
                {"wall_seconds", linear.wall_seconds},
                {"lightweight_large_base_seed", lightweight_large_base_seed},
            });
            auto round_json = linear.to_json(false);
            round_json["round"] = round;
            round_json["lightweight_large_base_seed"] =
                lightweight_large_base_seed;
            wall_seconds += linear.wall_seconds;
            if (!linear.success) {
                round_json["nonlinear_repair"] = nullptr;
                linearized_repair_json.push_back(std::move(round_json));
                break;
            }
            auto linear_state = linear.state;
            const double linear_objective =
                gravityx::rebuild_base_state_derived_fields(
                    data, commitment, linear_state);
            const auto linear_validation = gravityx::validate_state(
                data, gravityx::ModelMode::BaseSoft,
                linear_state, commitment);
            round_json["linear_objective"] = linear_objective;
            round_json["linear_validation"] = linear_validation.to_json();
            if (linear_validation.max_residual <
                selected_validation.max_residual) {
                selected_solve.status = 0;
                selected_solve.objective = linear_objective;
                selected_solve.wall_seconds = linear.wall_seconds;
                selected_solve.iterations = linear.iterations;
                selected_solve.state = linear_state;
                selected_validation = linear_validation;
            }
            if (linear_validation.max_residual <= 1e-5) {
                selected_solve.status = 0;
                selected_solve.objective = linear_objective;
                selected_solve.wall_seconds = linear.wall_seconds;
                selected_solve.iterations = linear.iterations;
                selected_solve.state = std::move(linear_state);
                selected_validation = linear_validation;
                success = true;
                base_method = "highs_sequential_linearized_ac";
                round_json["nonlinear_repair"] = nullptr;
                linearized_repair_json.push_back(std::move(round_json));
                break;
            }
            gravityx::FastContingencyPowerFlow nonlinear_repair(
                data, linear_state, commitment);
            auto nonlinear = nonlinear_repair.solve_base();
            wall_seconds += nonlinear.wall_seconds;
            round_json["nonlinear_repair"] = nonlinear.to_json();
            nlohmann::json nonlinear_restarts = nlohmann::json::array();
            // A failed Newton pass can leave a much better voltage/dispatch
            // center than the LP candidate that started it.  Restarting from
            // that center costs about a second at 8k buses and can avoid a
            // second 45-second LP.  Keep at most two strictly improving
            // restarts; none is accepted without the same full validator.
            for (int restart = 1;
                 allow_large_base_newton_restart &&
                 data.buses.size() >= 8000 && restart <= 2 &&
                 !nonlinear.feasible &&
                 nonlinear.validation.max_residual <= 0.25;
                 ++restart) {
                gravityx::FastContingencyPowerFlow restart_repair(
                    data, nonlinear.solve.state, commitment);
                const auto restarted = restart_repair.solve_base();
                wall_seconds += restarted.wall_seconds;
                auto restart_json = restarted.to_json();
                restart_json["restart"] = restart;
                nonlinear_restarts.push_back(std::move(restart_json));
                if (restarted.feasible ||
                    restarted.validation.max_residual + 1e-12 <
                        nonlinear.validation.max_residual) {
                    nonlinear = restarted;
                } else {
                    break;
                }
            }
            round_json["nonlinear_restarts"] =
                std::move(nonlinear_restarts);
            linearized_repair_json.push_back(std::move(round_json));
            if (nonlinear.validation.max_residual <
                selected_validation.max_residual) {
                selected_solve = nonlinear.solve;
                selected_validation = nonlinear.validation;
            }
            if (nonlinear.feasible) {
                selected_solve = nonlinear.solve;
                selected_validation = nonlinear.validation;
                success = true;
                base_method = "highs_linearized_ac_plus_validated_sparse_newton";
                break;
            }
            if (nonlinear.converged ||
                nonlinear.validation.max_residual < linear_reference_residual) {
                linear_reference = nonlinear.solve.state;
                linear_reference_residual = nonlinear.validation.max_residual;
            } else {
                linear_reference = linear.state;
            }
        }
    }
    if (!success && allow_exact_fallback) {
        const double exact_initial_residual = selected_validation.max_residual;
        log_base_phase("exact_fallback_start", {
            {"initial_max_residual", exact_initial_residual},
            {"initial_worst_category", selected_validation.worst_category},
            {"initial_worst_identity", selected_validation.worst_identity},
        });
        const auto exact_start = std::chrono::steady_clock::now();
        gravityx::AcModel exact_model(
            data, gravityx::ModelMode::BaseSoft, commitment);
        exact_model.initialize_from(selected_solve.state);
        const auto exact = exact_model.solve(0, 1e-6);
        const auto exact_validation = gravityx::validate_state(
            data, gravityx::ModelMode::BaseSoft, exact.state, commitment);
        const double exact_total_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - exact_start).count();
        wall_seconds += exact_total_seconds;
        success = (exact.status == 0 || exact.status == 1) &&
            gravityx::validated_candidate_is_feasible(
                exact, exact_validation, 1e-5);
        selected_solve = exact;
        selected_validation = exact_validation;
        exact_repair_json = {
            {"success", success},
            {"total_wall_seconds", exact_total_seconds},
            {"solve", gravityx::solve_result_to_json(exact, false)},
            {"validation", exact_validation.to_json()},
            {"initial_max_residual", exact_initial_residual},
        };
        base_optimization_performed = true;
        base_method = "gravity_ipopt_from_validated_fast_ac_seed";
        log_base_phase("exact_fallback_complete", {
            {"success", success},
            {"solve_status", exact.status},
            {"iterations", exact.iterations},
            {"max_residual", exact_validation.max_residual},
            {"worst_category", exact_validation.worst_category},
            {"worst_identity", exact_validation.worst_identity},
            {"wall_seconds", exact_total_seconds},
        });
    }
    gravityx::IbrResult result;
    result.success = success;
    result.wall_seconds = wall_seconds;
    result.base = selected_solve;
    result.base_validation = selected_validation;
    result.commitment = commitment;
    result.selected_state = selected_solve.state;
    auto output = result.to_json(true);
    output["base_method"] = base_method;
    output["base_optimization_performed"] = base_optimization_performed;
    output["exact_fallback_allowed"] = allow_exact_fallback;
    output["source_candidate"] = source.to_json();
    output["source_repair"] = repair_json;
    output["linearized_repair"] = linearized_repair_json;
    output["exact_repair"] = exact_repair_json;
    write_json_file(output_path, output);
    std::cout << nlohmann::json({
        {"output", output_path},
        {"success", success},
        {"base_method", base_method},
        {"wall_seconds", wall_seconds},
        {"max_residual", selected_validation.max_residual},
        {"worst_category", selected_validation.worst_category},
        {"worst_identity", selected_validation.worst_identity},
    }).dump(2) << '\n';
    return success ? 0 : 1;
}

struct BasePoint {
    std::vector<int> commitment;
    gravityx::AcState state;
};

struct CorrectiveSeed {
    std::string label;
    gravityx::AcState state;
};

struct ContingencyComputation {
    nlohmann::json result;
    gravityx::AcState state;
};

BasePoint load_base_point(const std::string& base_result_path) {
    const auto base_json = read_json_file(base_result_path);
    if (!base_json.value("success", false)) {
        throw std::runtime_error("base IBR result was not successful");
    }
    if (!base_json.contains("commitment") || !base_json.contains("selected_state")) {
        throw std::runtime_error("base IBR result lacks commitment or selected_state");
    }
    return {
        base_json.at("commitment").get<std::vector<int>>(),
        gravityx::ac_state_from_json(base_json.at("selected_state")),
    };
}

bool solve_loaded_contingency(
    const gravityx::CaseData& data,
    const BasePoint& base,
    const std::string& label,
    const std::string& output_path,
    int print_level,
    std::unique_ptr<gravityx::AcModel>* reusable_model = nullptr,
    bool acceptable_termination = false,
    gravityx::FastContingencyPowerFlow* fast_power_flow = nullptr,
    bool fast_only = false,
    bool linearized_fallback = false,
    bool linearized_only = false,
    const gravityx::AcState* precomputed_fast_state = nullptr,
    const gravityx::AcState* rolling_corrective_seed = nullptr,
    const std::string* rolling_corrective_seed_label = nullptr,
    const std::vector<CorrectiveSeed>* corrective_seed_bank = nullptr,
    std::optional<ContingencyComputation>* completed_computation = nullptr,
    bool persist_result = true) {
    reject_onedrive(output_path);
    const auto match = std::find_if(
        data.contingencies.begin(), data.contingencies.end(),
        [&label](const gravityx::Contingency& item) { return item.label == label; });
    if (match == data.contingencies.end()) {
        throw std::runtime_error("unknown contingency label: " + label);
    }
    gravityx::ContingencyContext context;
    context.base_state = base.state;
    if (match->type == gravityx::ContingencyType::Generator) {
        context.outaged_generator = match->component;
    } else {
        context.outaged_branch = match->component;
    }

    auto complete = [&](nlohmann::json output,
                        const gravityx::AcState& state,
                        bool persist_full_state = false) {
        if (persist_result) {
            output["solve"]["state"] = persist_full_state
                ? gravityx::ac_state_to_json(state)
                : gravityx::ac_submission_state_to_json(state);
            write_json_file(output_path, output);
        }
        if (completed_computation != nullptr) {
            completed_computation->emplace(ContingencyComputation{
                std::move(output),
                state,
            });
        }
    };

    std::optional<gravityx::FastPowerFlowResult> fast_result;
    const bool reused_fast_screen_reference = precomputed_fast_state != nullptr;
    bool rolling_seed_fast_screen_selected = false;
    bool bounded_fast_newton_rescue_selected = false;
    bool passive_pocket_repair_selected = false;
    nlohmann::json passive_pocket_repair_diagnostics = nullptr;
    bool bounded_fast_linearized_repair_selected = false;
    bool bounded_fast_postlinear_newton_selected = false;
    std::optional<gravityx::ActiveFeasibilityRepairResult>
        bounded_fast_linearized_repair;
    nlohmann::json bounded_fast_linearized_repair_attempts =
        nlohmann::json::array();
    std::optional<gravityx::FastPowerFlowResult>
        bounded_fast_postlinear_newton;
    std::optional<std::string> selected_direct_seed_label;
    if (precomputed_fast_state) {
        fast_result.emplace();
        fast_result->solve.state = *precomputed_fast_state;
        fast_result->failure_reason = "precomputed_fast_screen_reference";
    } else if (fast_power_flow) {
        const bool rolling_seed_dimensions_match =
            rolling_corrective_seed != nullptr &&
            rolling_corrective_seed->vm.size() == data.buses.size() &&
            rolling_corrective_seed->va.size() == data.buses.size() &&
            rolling_corrective_seed->pg.size() == data.generators.size() &&
            rolling_corrective_seed->qg.size() == data.generators.size() &&
            rolling_corrective_seed->demand_factor.size() == data.loads.size();
        if (rolling_seed_dimensions_match) {
            auto translated_rolling_seed = *rolling_corrective_seed;
            if (rolling_corrective_seed_label != nullptr &&
                match->type == gravityx::ContingencyType::Generator) {
                const auto prior_match = std::find_if(
                    data.contingencies.begin(), data.contingencies.end(),
                    [&](const gravityx::Contingency& item) {
                        return item.label == *rolling_corrective_seed_label;
                    });
                if (prior_match != data.contingencies.end() &&
                    prior_match->type ==
                        gravityx::ContingencyType::Generator &&
                    data.generators[prior_match->component].bus ==
                        data.generators[match->component].bus) {
                    const int prior_outage = prior_match->component;
                    const int current_outage = match->component;
                    translated_rolling_seed.pg[prior_outage] =
                        translated_rolling_seed.pg[current_outage];
                    translated_rolling_seed.qg[prior_outage] =
                        translated_rolling_seed.qg[current_outage];
                    translated_rolling_seed.pg[current_outage] = 0.0;
                    translated_rolling_seed.qg[current_outage] = 0.0;
                }
            }
            fast_result = fast_power_flow->screen_candidate(
                *match, translated_rolling_seed);
            rolling_seed_fast_screen_selected = fast_result->feasible;
            if (rolling_seed_fast_screen_selected &&
                rolling_corrective_seed_label != nullptr) {
                selected_direct_seed_label =
                    *rolling_corrective_seed_label;
            }
        }
        if (!rolling_seed_fast_screen_selected) {
            fast_result = fast_power_flow->solve(*match);
        }
        // A failed predictor is rare and the exact corrective fallback is
        // expensive on the largest cases.  Before escalating, direct-screen
        // the small resident bank of previously verified corrective states.
        // This does not transfer feasibility: screen_candidate rebuilds the
        // requested outage and applies the complete independent validator.
        if (!fast_result->feasible && corrective_seed_bank != nullptr) {
            for (const auto& seed : *corrective_seed_bank) {
                if (rolling_corrective_seed_label != nullptr &&
                    seed.label == *rolling_corrective_seed_label) {
                    continue;
                }
                auto bank_screen = fast_power_flow->screen_candidate(
                    *match, seed.state);
                if (!bank_screen.feasible) {
                    continue;
                }
                fast_result = std::move(bank_screen);
                selected_direct_seed_label = seed.label;
                break;
            }
        }
        if (fast_only && !fast_result->feasible &&
            bounded_fast_candidate_repair_candidate(
                fast_result->validation)) {
            const double prior_screen_seconds = fast_result->wall_seconds;
            const auto rescue_start = std::chrono::steady_clock::now();
            gravityx::FastPowerFlowOptions rescue_options;
            rescue_options.max_newton_iterations = 12;
            rescue_options.max_active_redispatch_passes = 4;
            rescue_options.max_reactive_limit_passes = 4;
            gravityx::FastContingencyPowerFlow rescue_solver(
                data, base.state, base.commitment, rescue_options);
            auto rescue_result = rescue_solver.solve(
                *match, fast_result->solve.state);
            const double combined_wall_seconds = prior_screen_seconds +
                std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - rescue_start).count();
            rescue_result.wall_seconds = combined_wall_seconds;
            rescue_result.solve.wall_seconds = combined_wall_seconds;
            // Retain a strictly better validated rescue even when it has not
            // crossed the feasibility tolerance.  This preserves the actual
            // rescue diagnostics and gives any exact fallback the best known
            // state, without weakening the independent feasibility gate.
            if (rescue_result.feasible ||
                rescue_result.validation.max_residual + 1e-12 <
                    fast_result->validation.max_residual) {
                fast_result = std::move(rescue_result);
                bounded_fast_newton_rescue_selected = fast_result->feasible;
            }
            if (!fast_result->feasible &&
                fast_result->validation.max_residual <= 0.2 &&
                match->type == gravityx::ContingencyType::Branch) {
                auto pocket_repair = try_passive_outage_pocket_repair(
                    data, base.state, base.commitment, *match, context,
                    fast_result->solve.state, fast_result->validation);
                if (pocket_repair) {
                    passive_pocket_repair_diagnostics =
                        std::move(pocket_repair->diagnostics);
                    if (pocket_repair->result) {
                        const double combined_pocket_wall_seconds =
                            fast_result->wall_seconds +
                            pocket_repair->result->wall_seconds;
                        pocket_repair->result->wall_seconds =
                            combined_pocket_wall_seconds;
                        pocket_repair->result->solve.wall_seconds =
                            combined_pocket_wall_seconds;
                        if (pocket_repair->result->feasible ||
                            pocket_repair->result->validation.max_residual +
                                    1e-12 <
                                fast_result->validation.max_residual) {
                            fast_result =
                                std::move(*pocket_repair->result);
                            passive_pocket_repair_selected =
                                fast_result->feasible;
                        }
                    }
                }
            }
            if (!fast_result->feasible &&
                fast_result->validation.worst_category == "variable_bound" &&
                fast_result->validation.max_residual <= 0.1) {
                double wall_after_linearized_repair =
                    fast_result->wall_seconds;
                auto linearized_reference = fast_result->solve.state;
                bool linearized_repair_produced_candidate = false;
                constexpr int kMaximumCompactRepairRounds = 3;
                for (int compact_round = 1;
                     compact_round <= kMaximumCompactRepairRounds;
                     ++compact_round) {
                    bounded_fast_linearized_repair =
                        gravityx::solve_linearized_active_feasibility_repair(
                            data, linearized_reference, base.commitment,
                            context, 0.45, 0.15, 5.0, 0.1, true, true);
                    wall_after_linearized_repair +=
                        bounded_fast_linearized_repair->wall_seconds;
                    auto compact_attempt =
                        bounded_fast_linearized_repair->to_json(false);
                    compact_attempt["round"] = compact_round;
                    if (!bounded_fast_linearized_repair->success) {
                        compact_attempt["exact_validation"] = nullptr;
                        bounded_fast_linearized_repair_attempts.push_back(
                            std::move(compact_attempt));
                        bounded_fast_linearized_repair =
                            gravityx::solve_linearized_active_feasibility_repair(
                                data, linearized_reference, base.commitment,
                                context, 0.49, 0.5, 5.0, 0.2, true, true);
                        wall_after_linearized_repair +=
                            bounded_fast_linearized_repair->wall_seconds;
                        compact_attempt =
                            bounded_fast_linearized_repair->to_json(false);
                        compact_attempt["round"] = compact_round;
                        compact_attempt["wide_trust_retry"] = true;
                        if (!bounded_fast_linearized_repair->success) {
                            compact_attempt["exact_validation"] = nullptr;
                            bounded_fast_linearized_repair_attempts.push_back(
                                std::move(compact_attempt));
                            break;
                        }
                    }
                    linearized_repair_produced_candidate = true;
                    auto linearized_state =
                        bounded_fast_linearized_repair->state;
                    const double linearized_objective =
                        gravityx::rebuild_contingency_state_derived_fields(
                            data, base.state, base.commitment, *match,
                            linearized_state);
                    const auto linearized_validation =
                        gravityx::validate_state(
                            data, gravityx::ModelMode::ContingencySoft,
                            linearized_state, base.commitment, context);
                    compact_attempt["exact_objective"] =
                        linearized_objective;
                    compact_attempt["exact_validation"] =
                        linearized_validation.to_json();
                    const bool improved =
                        linearized_validation.max_residual + 1e-12 <
                        fast_result->validation.max_residual;
                    compact_attempt["selected_as_best"] = improved;
                    bounded_fast_linearized_repair_attempts.push_back(
                        std::move(compact_attempt));
                    linearized_reference = linearized_state;
                    if (improved) {
                        gravityx::FastPowerFlowResult linearized_result;
                        linearized_result.converged = false;
                        linearized_result.feasible =
                            linearized_validation.max_residual <= 1e-5;
                        linearized_result.wall_seconds =
                            wall_after_linearized_repair;
                        linearized_result.failure_reason =
                            linearized_result.feasible
                            ? std::string{}
                            : "compact linearized feasibility repair requires "
                              "nonlinear polish";
                        linearized_result.solve.status =
                            linearized_result.feasible ? 0 : 1;
                        linearized_result.solve.objective =
                            linearized_objective;
                        linearized_result.solve.wall_seconds =
                            wall_after_linearized_repair;
                        linearized_result.solve.state = linearized_state;
                        linearized_result.validation = linearized_validation;
                        fast_result = std::move(linearized_result);
                        bounded_fast_linearized_repair_selected =
                            fast_result->feasible;
                    }
                    if (fast_result->feasible) {
                        break;
                    }
                }
                if (!fast_result->feasible &&
                    linearized_repair_produced_candidate) {
                    gravityx::FastPowerFlowOptions postlinear_options;
                    postlinear_options.max_newton_iterations = 20;
                    postlinear_options.max_active_redispatch_passes = 6;
                    postlinear_options.max_reactive_limit_passes = 6;
                    gravityx::FastContingencyPowerFlow postlinear_solver(
                        data, base.state, base.commitment,
                        postlinear_options);
                    bounded_fast_postlinear_newton = postlinear_solver.solve(
                        *match, linearized_reference);
                    const double combined_postlinear_wall_seconds =
                        wall_after_linearized_repair +
                        bounded_fast_postlinear_newton->wall_seconds;
                    bounded_fast_postlinear_newton->wall_seconds =
                        combined_postlinear_wall_seconds;
                    bounded_fast_postlinear_newton->solve.wall_seconds =
                        combined_postlinear_wall_seconds;
                    if (bounded_fast_postlinear_newton->feasible ||
                        bounded_fast_postlinear_newton->validation
                                .max_residual + 1e-12 <
                            fast_result->validation.max_residual) {
                        fast_result = *bounded_fast_postlinear_newton;
                        bounded_fast_postlinear_newton_selected =
                            fast_result->feasible;
                    }
                    constexpr int kMaximumPostNewtonRecenters = 2;
                    for (int recenter_round = 1;
                         !fast_result->feasible &&
                         fast_result->validation.max_residual <= 0.1 &&
                         recenter_round <= kMaximumPostNewtonRecenters;
                         ++recenter_round) {
                        const double residual_before_recenter =
                            fast_result->validation.max_residual;
                        bounded_fast_linearized_repair =
                            gravityx::solve_linearized_active_feasibility_repair(
                                data, fast_result->solve.state,
                                base.commitment, context,
                                0.49, 0.5, 5.0, 0.2, true, true);
                        const double wall_after_recenter =
                            fast_result->wall_seconds +
                            bounded_fast_linearized_repair->wall_seconds;
                        auto recenter_attempt =
                            bounded_fast_linearized_repair->to_json(false);
                        recenter_attempt["phase"] =
                            "post_newton_recenter";
                        recenter_attempt["round"] = recenter_round;
                        if (!bounded_fast_linearized_repair->success) {
                            recenter_attempt["exact_validation"] = nullptr;
                            bounded_fast_linearized_repair_attempts.push_back(
                                std::move(recenter_attempt));
                            break;
                        }
                        auto recentered_state =
                            bounded_fast_linearized_repair->state;
                        const double recentered_objective =
                            gravityx::rebuild_contingency_state_derived_fields(
                                data, base.state, base.commitment, *match,
                                recentered_state);
                        const auto recentered_validation =
                            gravityx::validate_state(
                                data, gravityx::ModelMode::ContingencySoft,
                                recentered_state, base.commitment, context);
                        recenter_attempt["exact_objective"] =
                            recentered_objective;
                        recenter_attempt["exact_validation"] =
                            recentered_validation.to_json();
                        const bool recentered_improved =
                            recentered_validation.max_residual + 1e-12 <
                            fast_result->validation.max_residual;
                        recenter_attempt["selected_as_best"] =
                            recentered_improved;
                        if (recentered_improved) {
                            gravityx::FastPowerFlowResult recentered_result;
                            recentered_result.converged = false;
                            recentered_result.feasible =
                                recentered_validation.max_residual <= 1e-5;
                            recentered_result.wall_seconds =
                                wall_after_recenter;
                            recentered_result.solve.status =
                                recentered_result.feasible ? 0 : 1;
                            recentered_result.solve.objective =
                                recentered_objective;
                            recentered_result.solve.wall_seconds =
                                wall_after_recenter;
                            recentered_result.solve.state = recentered_state;
                            recentered_result.validation =
                                recentered_validation;
                            fast_result = std::move(recentered_result);
                            bounded_fast_linearized_repair_selected =
                                fast_result->feasible;
                        }
                        if (!fast_result->feasible) {
                            gravityx::FastPowerFlowOptions recentered_options;
                            recentered_options.max_newton_iterations = 20;
                            recentered_options.max_active_redispatch_passes = 6;
                            recentered_options.max_reactive_limit_passes = 6;
                            gravityx::FastContingencyPowerFlow
                                recentered_solver(
                                    data, base.state, base.commitment,
                                    recentered_options);
                            auto recentered_newton =
                                recentered_solver.solve(
                                    *match, recentered_state);
                            const double combined_recenter_wall_seconds =
                                wall_after_recenter +
                                recentered_newton.wall_seconds;
                            recentered_newton.wall_seconds =
                                combined_recenter_wall_seconds;
                            recentered_newton.solve.wall_seconds =
                                combined_recenter_wall_seconds;
                            recenter_attempt["nonlinear_polish"] =
                                recentered_newton.to_json();
                            if (recentered_newton.feasible ||
                                recentered_newton.validation.max_residual +
                                        1e-12 <
                                    fast_result->validation.max_residual) {
                                bounded_fast_postlinear_newton =
                                    recentered_newton;
                                fast_result = std::move(recentered_newton);
                                bounded_fast_postlinear_newton_selected =
                                    fast_result->feasible;
                            }
                        }
                        bounded_fast_linearized_repair_attempts.push_back(
                            std::move(recenter_attempt));
                        if (fast_result->feasible ||
                            fast_result->validation.max_residual + 1e-12 >=
                                residual_before_recenter) {
                            break;
                        }
                    }
                }
            }
        }
        if (fast_result->feasible) {
            const std::string solution_method =
                bounded_fast_postlinear_newton_selected
                ? "bounded_fast_linearized_repair_plus_newton"
                : bounded_fast_linearized_repair_selected
                ? "bounded_fast_linearized_feasibility_repair"
                : passive_pocket_repair_selected
                ? "small_outage_pocket_equalization"
                : bounded_fast_newton_rescue_selected
                ? "bounded_fast_newton_rescue"
                : rolling_seed_fast_screen_selected
                ? "rolling_corrective_seed_direct_screen"
                : selected_direct_seed_label
                ? "corrective_seed_bank_direct_screen"
                : fast_result->direct_candidate_selected
                ? "direct_base_state_outage_candidate"
                : fast_result->fixed_jacobian_predictor_selected
                ? "resident_fixed_jacobian_predictor"
                : "fast_newton_power_flow";
            nlohmann::json output = {
                {"success", true},
                {"solver_status_success", true},
                {"accepted_feasible_nonconverged", false},
                {"label", match->label},
                {"type", match->type == gravityx::ContingencyType::Generator ? "gen" : "branch"},
                {"source_index", match->source_index},
                {"component_position", match->component},
                {"solution_method", solution_method},
                {"fast_power_flow_screen", true},
                {"fast_screen", persist_result
                    ? fast_result->to_json()
                    : nlohmann::json({
                          {"failure_reason", fast_result->failure_reason},
                          {"wall_seconds", fast_result->wall_seconds},
                      })},
                {"bounded_fast_linearized_repair",
                 bounded_fast_linearized_repair
                     ? bounded_fast_linearized_repair->to_json(false)
                     : nlohmann::json(nullptr)},
                {"bounded_fast_linearized_repair_attempts",
                 bounded_fast_linearized_repair_attempts},
                {"bounded_fast_postlinear_newton",
                 bounded_fast_postlinear_newton
                     ? bounded_fast_postlinear_newton->to_json()
                     : nlohmann::json(nullptr)},
                {"passive_pocket_repair",
                 passive_pocket_repair_diagnostics},
                {"rolling_corrective_seed_label",
                 selected_direct_seed_label
                     ? nlohmann::json(*selected_direct_seed_label)
                     : nlohmann::json(nullptr)},
                {"resident_parametric_model", false},
                {"acceptable_termination_enabled", false},
                {"resident_model_created", false},
                {"model_preparation_wall_seconds", 0.0},
                {"solve", gravityx::solve_result_to_json(
                    fast_result->solve, false)},
                {"validation", fast_result->validation.to_json()},
            };
            complete(std::move(output), fast_result->solve.state);
            std::cout << nlohmann::json({
                {"output", output_path},
                {"success", true},
                {"solver_status_success", true},
                {"accepted_feasible_nonconverged", false},
                {"label", match->label},
                {"status", 0},
                {"objective", fast_result->solve.objective},
                {"wall_seconds", fast_result->wall_seconds},
                {"iterations", fast_result->newton_iterations},
                {"resident_reoptimization", false},
                {"model_preparation_wall_seconds", 0.0},
                {"max_residual", fast_result->validation.max_residual},
                {"solution_method", solution_method},
            }).dump(2) << '\n';
            return true;
        }
        if (fast_only) {
            nlohmann::json output = {
                {"success", false},
                {"screen_completed", true},
                {"requires_exact_fallback", true},
                {"solver_status_success", false},
                {"accepted_feasible_nonconverged", false},
                {"label", match->label},
                {"type", match->type == gravityx::ContingencyType::Generator ? "gen" : "branch"},
                {"source_index", match->source_index},
                {"component_position", match->component},
                {"solution_method", "fast_newton_screen_failed"},
                {"fast_power_flow_screen", true},
                {"fast_screen", fast_result->to_json()},
                {"bounded_fast_linearized_repair",
                 bounded_fast_linearized_repair
                     ? bounded_fast_linearized_repair->to_json(false)
                     : nlohmann::json(nullptr)},
                {"bounded_fast_linearized_repair_attempts",
                 bounded_fast_linearized_repair_attempts},
                {"bounded_fast_postlinear_newton",
                 bounded_fast_postlinear_newton
                     ? bounded_fast_postlinear_newton->to_json()
                     : nlohmann::json(nullptr)},
                {"passive_pocket_repair",
                 passive_pocket_repair_diagnostics},
                {"resident_parametric_model", false},
                {"acceptable_termination_enabled", false},
                {"resident_model_created", false},
                {"model_preparation_wall_seconds", 0.0},
                {"solve", gravityx::solve_result_to_json(
                    fast_result->solve, false)},
                {"validation", fast_result->validation.to_json()},
            };
            complete(std::move(output), fast_result->solve.state);
            return true;
        }
    }

    std::optional<gravityx::AcState> linearized_seed;
    std::optional<gravityx::FastPowerFlowResult> prelinear_fast_repair;
    std::optional<gravityx::FastPowerFlowResult> rolling_fast_repair;
    nlohmann::json corrective_seed_fast_repair_attempts =
        nlohmann::json::array();
    nlohmann::json linearized_attempts = nlohmann::json::array();
    if (linearized_fallback) {
        auto best_state = fast_result
            ? fast_result->solve.state : base.state;
        double best_objective = gravityx::rebuild_contingency_state_derived_fields(
            data, base.state, base.commitment, *match, best_state);
        auto best_validation = gravityx::validate_state(
            data, gravityx::ModelMode::ContingencySoft,
            best_state, base.commitment, context);
        if (reused_fast_screen_reference) {
            fast_result->solve.objective = best_objective;
            fast_result->validation = best_validation;
        }
        // The large-case fast solver can discover a much better intermediate
        // AC state than its terminal iterate.  After preserving that state,
        // give Newton inexpensive restarts before constructing a 121k-column
        // LP.  A resident worker may also reuse its last independently verified
        // corrective state as an initialization for the new outage.  Neither
        // path transfers feasibility: every candidate is rebuilt for this
        // outage and must pass the complete independent validator.
        if (data.buses.size() >= 16000 && fast_power_flow != nullptr &&
            fast_result.has_value() &&
            best_validation.max_residual > 0.1) {
            double fast_repair_wall_seconds = fast_result->wall_seconds;
            auto evaluate_fast_repair = [&](
                gravityx::FastPowerFlowResult& repair,
                const std::string& solution_method,
                const std::string* corrective_seed_label) {
                fast_repair_wall_seconds += repair.wall_seconds;
                if (repair.solve.state.vm.empty()) {
                    if (corrective_seed_label != nullptr) {
                        corrective_seed_fast_repair_attempts.push_back({
                            {"seed_label", *corrective_seed_label},
                            {"result", repair.to_json()},
                        });
                    }
                    return false;
                }
                auto repaired_state = repair.solve.state;
                const double repaired_objective =
                    gravityx::rebuild_contingency_state_derived_fields(
                        data, base.state, base.commitment, *match,
                        repaired_state);
                const auto repaired_validation = gravityx::validate_state(
                    data, gravityx::ModelMode::ContingencySoft,
                    repaired_state, base.commitment, context);
                repair.solve.state = repaired_state;
                repair.solve.objective = repaired_objective;
                repair.validation = repaired_validation;
                repair.feasible = repaired_validation.max_residual <= 1e-5;
                if (corrective_seed_label != nullptr) {
                    corrective_seed_fast_repair_attempts.push_back({
                        {"seed_label", *corrective_seed_label},
                        {"result", repair.to_json()},
                    });
                }
                if (repaired_validation.max_residual <
                    best_validation.max_residual) {
                    best_state = repaired_state;
                    best_objective = repaired_objective;
                    best_validation = repaired_validation;
                }
                if (!repair.feasible) {
                    return false;
                }
                auto solve = repair.solve;
                solve.wall_seconds = fast_repair_wall_seconds;
                const nlohmann::json selected_corrective_seed =
                    corrective_seed_label != nullptr
                    ? nlohmann::json(*corrective_seed_label)
                    : nlohmann::json(nullptr);
                nlohmann::json output = {
                    {"success", true},
                    {"solver_status_success", true},
                    {"accepted_feasible_nonconverged", false},
                    {"label", match->label},
                    {"type", match->type == gravityx::ContingencyType::Generator ? "gen" : "branch"},
                    {"source_index", match->source_index},
                    {"component_position", match->component},
                    {"solution_method", solution_method},
                    {"precomputed_fast_screen_reference", reused_fast_screen_reference},
                    {"fast_power_flow_screen", true},
                    {"fast_screen", fast_result->to_json()},
                    {"prelinear_fast_repair", prelinear_fast_repair ? prelinear_fast_repair->to_json() : nlohmann::json(nullptr)},
                    {"rolling_fast_repair", rolling_fast_repair ? rolling_fast_repair->to_json() : nlohmann::json(nullptr)},
                    {"rolling_corrective_seed_label", selected_corrective_seed},
                    {"corrective_seed_fast_repair_attempts", corrective_seed_fast_repair_attempts},
                    {"linearized_attempts", linearized_attempts},
                    {"resident_parametric_model", false},
                    {"acceptable_termination_enabled", false},
                    {"resident_model_created", false},
                    {"model_preparation_wall_seconds", 0.0},
                    {"solve", gravityx::solve_result_to_json(solve, false)},
                    {"validation", repaired_validation.to_json()},
                };
                complete(std::move(output), solve.state);
                std::cout << nlohmann::json({
                    {"output", output_path},
                    {"success", true},
                    {"label", match->label},
                    {"status", solve.status},
                    {"objective", repaired_objective},
                    {"wall_seconds", solve.wall_seconds},
                    {"iterations", solve.iterations},
                    {"max_residual", repaired_validation.max_residual},
                    {"solution_method", solution_method},
                    {"corrective_seed_label", selected_corrective_seed},
                }).dump(2) << '\n';
                return true;
            };

            prelinear_fast_repair =
                fast_power_flow->solve(*match, best_state);
            if (evaluate_fast_repair(
                    *prelinear_fast_repair,
                    "iterated_fast_newton_power_flow", nullptr)) {
                return true;
            }

            const bool rolling_seed_dimensions_match =
                rolling_corrective_seed != nullptr &&
                rolling_corrective_seed->vm.size() == data.buses.size() &&
                rolling_corrective_seed->va.size() == data.buses.size() &&
                rolling_corrective_seed->pg.size() == data.generators.size() &&
                rolling_corrective_seed->qg.size() == data.generators.size() &&
                rolling_corrective_seed->demand_factor.size() == data.loads.size();
            if (rolling_seed_dimensions_match) {
                rolling_fast_repair = fast_power_flow->solve(
                    *match, *rolling_corrective_seed);
                if (evaluate_fast_repair(
                        *rolling_fast_repair,
                        "rolling_corrective_seed_fast_newton",
                        rolling_corrective_seed_label)) {
                    return true;
                }
            }
            if (corrective_seed_bank != nullptr) {
                for (const auto& seed : *corrective_seed_bank) {
                    if (rolling_corrective_seed_label != nullptr &&
                        seed.label == *rolling_corrective_seed_label) {
                        continue;
                    }
                    const bool seed_dimensions_match =
                        seed.state.vm.size() == data.buses.size() &&
                        seed.state.va.size() == data.buses.size() &&
                        seed.state.pg.size() == data.generators.size() &&
                        seed.state.qg.size() == data.generators.size() &&
                        seed.state.demand_factor.size() == data.loads.size();
                    if (!seed_dimensions_match) {
                        continue;
                    }
                    auto bank_repair = fast_power_flow->solve(
                        *match, seed.state);
                    if (evaluate_fast_repair(
                            bank_repair,
                            "corrective_seed_bank_fast_newton",
                            &seed.label)) {
                        return true;
                    }
                }
            }
        }
        const bool use_near_feasible_exact_polish =
            data.buses.size() >= 16000 &&
            best_validation.max_residual <= 0.1;
        double linearized_wall_seconds = 0.0;
        int linearized_iterations = 0;
        std::string last_status = "not_run";
        int last_model_status = -1;
        if (use_near_feasible_exact_polish) {
            linearized_seed = best_state;
        } else {
        std::vector<int> dynamic_security_branches;
        std::vector<unsigned char> dynamic_security_selected(
            data.branches.size(), 0);
        const auto collect_violated_security_branches =
            [&](const gravityx::AcState& state) {
                nlohmann::json added = nlohmann::json::array();
                if (state.pf.size() != data.branches.size() ||
                    state.qf.size() != data.branches.size() ||
                    state.pt.size() != data.branches.size() ||
                    state.qt.size() != data.branches.size() ||
                    state.sm_slack.size() != data.branches.size()) {
                    return added;
                }
                constexpr double kSecurityCollectionTolerance = 1e-5;
                for (std::size_t i = 0; i < data.branches.size(); ++i) {
                    if (static_cast<int>(i) == context.outaged_branch ||
                        data.branches[i].status == 0) {
                        continue;
                    }
                    const auto& branch = data.branches[i];
                    const double rating = branch.rate_c;
                    const double box_violation = std::max({
                        std::abs(state.pf[i]) - rating,
                        std::abs(state.qf[i]) - rating,
                        std::abs(state.pt[i]) - rating,
                        std::abs(state.qt[i]) - rating,
                    });
                    const double source_delta =
                        base.state.va[branch.from] - base.state.va[branch.to];
                    double angle_violation = 0.0;
                    if (source_delta >= branch.angmin &&
                        source_delta <= branch.angmax) {
                        const double angle =
                            state.va[branch.from] - state.va[branch.to];
                        angle_violation = std::max(
                            angle - branch.angmax,
                            branch.angmin - angle);
                    }
                    const double slack = state.sm_slack[i];
                    const double from_scale = branch.transformer
                        ? 1.0 + slack : state.vm[branch.from] + slack;
                    const double to_scale = branch.transformer
                        ? 1.0 + slack : state.vm[branch.to] + slack;
                    const double apparent_violation = std::max(
                        state.pf[i] * state.pf[i] +
                            state.qf[i] * state.qf[i] -
                            rating * rating * from_scale * from_scale,
                        state.pt[i] * state.pt[i] +
                            state.qt[i] * state.qt[i] -
                            rating * rating * to_scale * to_scale);
                    if (std::max({box_violation, angle_violation,
                                  apparent_violation}) <=
                        kSecurityCollectionTolerance ||
                        dynamic_security_selected[i]) {
                        continue;
                    }
                    dynamic_security_selected[i] = 1;
                    dynamic_security_branches.push_back(
                        static_cast<int>(i));
                    added.push_back({
                        {"component_position", static_cast<int>(i)},
                        {"source_key", branch.source_key},
                        {"box_violation", std::max(0.0, box_violation)},
                        {"angle_violation", std::max(0.0, angle_violation)},
                        {"apparent_flow_violation",
                         std::max(0.0, apparent_violation)},
                    });
                }
                std::sort(dynamic_security_branches.begin(),
                          dynamic_security_branches.end());
                return added;
            };
        const auto initial_dynamic_security_branches =
            collect_violated_security_branches(best_state);
        auto reference = best_state;
        std::string reference_source = fast_result ? "fast_screen" : "base";
        constexpr int kMaximumLinearizedRounds = 3;
        for (int round = 1; round <= kMaximumLinearizedRounds; ++round) {
            // At 16k-bus scale the full security-row seed exceeds 120k rows
            // and repeatedly exhausts its 60-second LP limit before producing
            // a candidate.  Start with network balance, operating bounds, and
            // trust regions only.  This is candidate generation, not an
            // acceptance relaxation: sparse nonlinear repair and the complete
            // validator still check every source branch rating and angle row.
            const bool balance_only_phase_one =
                data.buses.size() >= 16000;
            // Keep a 0.25-p.u. interior margin below the source model's
            // 0.5-p.u. corrective balance limit.  The previous 0.49 target
            // left only 0.01 p.u. for AC linearization error; CTG_000005's
            // returned linear point was LP-feasible but 0.115 p.u. outside
            // the exact nonlinear balance envelope.
            const double linearized_balance_slack =
                balance_only_phase_one ? 0.25 : 0.49;
            const double linearized_time_limit_seconds =
                balance_only_phase_one ? 90.0 : 60.0;
            auto linear = gravityx::solve_linearized_ac_seed(
                data, reference, base.commitment,
                linearized_balance_slack, context,
                balance_only_phase_one, false,
                linearized_time_limit_seconds, balance_only_phase_one,
                balance_only_phase_one, dynamic_security_branches);
            linearized_wall_seconds += linear.wall_seconds;
            linearized_iterations += std::max(0, linear.iterations);
            auto attempt = linear.to_json(false);
            attempt["round"] = round;
            attempt["reference_source"] = reference_source;
            attempt["balance_only_phase_one"] = balance_only_phase_one;
            attempt["linearized_balance_slack"] =
                linearized_balance_slack;
            attempt["linearized_time_limit_seconds"] =
                linearized_time_limit_seconds;
            attempt["initial_dynamic_security_branches"] =
                initial_dynamic_security_branches;
            attempt["active_dynamic_security_branch_positions"] =
                dynamic_security_branches;
            if (!linear.success && data.buses.size() >= 16000 &&
                !dynamic_security_branches.empty() &&
                (linear.status == "Infeasible" ||
                 linear.status == "Unknown")) {
                // A violated security row can be unreachable inside the
                // deliberately tight targeted trust region even though the
                // complete corrective problem is feasible.  Retrying that
                // identical restricted LP cannot change the answer.  First
                // solve the elastic balance Phase I with the wider ordinary
                // trust region; its exact nonlinear image is screened below,
                // and any still-violated security rows remain in the dynamic
                // set for the next round.
                attempt["balance_phase_one_retry_scheduled"] = true;
                attempt["delayed_security_branch_positions"] =
                    dynamic_security_branches;
                linearized_attempts.push_back(std::move(attempt));
                constexpr double kLargeBalancePhaseOneTimeLimitSeconds = 90.0;
                linear = gravityx::solve_linearized_ac_seed(
                    data, reference, base.commitment,
                    0.25,
                    context, true, false,
                    kLargeBalancePhaseOneTimeLimitSeconds,
                    true, true, {}, true);
                linearized_wall_seconds += linear.wall_seconds;
                linearized_iterations += std::max(0, linear.iterations);
                attempt = linear.to_json(false);
                attempt["round"] = round;
                attempt["balance_phase_one_retry"] = true;
                attempt["delayed_security_branch_positions"] =
                    dynamic_security_branches;
                attempt["reference_source"] = reference_source;
            } else if (!linear.success && data.buses.size() >= 8000 &&
                       data.buses.size() < 16000 &&
                       (linear.status == "Infeasible" ||
                        linear.status == "Unknown")) {
                attempt["projected_balance_retry_scheduled"] = true;
                linearized_attempts.push_back(std::move(attempt));
                // Under eight-way contention the large projected-balance IPM
                // can be only a few iterations short at the historical
                // 60-second cutoff.  Let that same Phase-I formulation finish
                // instead of discarding it and rebuilding two costlier
                // reference retries.  The runner's per-worker and global
                // deadlines remain the hard outer bounds.
                constexpr double kLargeProjectedBalanceTimeLimitSeconds = 90.0;
                linear = gravityx::solve_linearized_ac_seed(
                    data, reference, base.commitment,
                    balance_only_phase_one ? 0.25 : 0.49,
                    context, true, false,
                    kLargeProjectedBalanceTimeLimitSeconds,
                    balance_only_phase_one, balance_only_phase_one,
                    dynamic_security_branches);
                linearized_wall_seconds += linear.wall_seconds;
                linearized_iterations += std::max(0, linear.iterations);
                attempt = linear.to_json(false);
                attempt["round"] = round;
                attempt["projected_balance_retry"] = true;
                attempt["reference_source"] = reference_source;
            }
            last_status = linear.status;
            last_model_status = linear.model_status;
            if (!linear.success) {
                if (fast_result && reference_source != "base") {
                    attempt["base_reference_retry_scheduled"] = true;
                    linearized_attempts.push_back(std::move(attempt));
                    reference = base.state;
                    reference_source = "base";
                    continue;
                }
                attempt["validation"] = nullptr;
                linearized_attempts.push_back(std::move(attempt));
                break;
            }

            auto candidate = linear.state;
            const auto maximum_step = [](const std::vector<double>& left,
                                         const std::vector<double>& right) {
                double maximum = 0.0;
                const auto count = std::min(left.size(), right.size());
                for (std::size_t i = 0; i < count; ++i) {
                    maximum = std::max(maximum, std::abs(left[i] - right[i]));
                }
                return maximum;
            };
            attempt["maximum_voltage_magnitude_step"] =
                maximum_step(candidate.vm, reference.vm);
            attempt["maximum_voltage_angle_step"] =
                maximum_step(candidate.va, reference.va);
            attempt["maximum_active_generation_step"] =
                maximum_step(candidate.pg, reference.pg);
            attempt["maximum_reactive_generation_step"] =
                maximum_step(candidate.qg, reference.qg);
            attempt["maximum_demand_factor_step"] =
                maximum_step(candidate.demand_factor, reference.demand_factor);
            const double objective =
                gravityx::rebuild_contingency_state_derived_fields(
                    data, base.state, base.commitment, *match, candidate);
            const auto validation = gravityx::validate_state(
                data, gravityx::ModelMode::ContingencySoft,
                candidate, base.commitment, context);
            attempt["exact_objective"] = objective;
            attempt["validation"] = validation.to_json();
            if (validation.max_residual < best_validation.max_residual) {
                best_state = candidate;
                best_objective = objective;
                best_validation = validation;
            }
            if (validation.max_residual <= 1e-5) {
                attempt["nonlinear_repair"] = nullptr;
                linearized_attempts.push_back(std::move(attempt));
                gravityx::SolveResult solve;
                solve.status = 0;
                solve.objective = objective;
                solve.wall_seconds = linearized_wall_seconds;
                solve.iterations = linearized_iterations;
                solve.state = std::move(candidate);
                nlohmann::json output = {
                    {"success", true},
                    {"solver_status_success", true},
                    {"accepted_feasible_nonconverged", false},
                    {"label", match->label},
                    {"type", match->type == gravityx::ContingencyType::Generator ? "gen" : "branch"},
                    {"source_index", match->source_index},
                    {"component_position", match->component},
                    {"solution_method", "highs_sequential_linearized_contingency"},
                    {"precomputed_fast_screen_reference", reused_fast_screen_reference},
                    {"fast_power_flow_screen", fast_power_flow != nullptr},
                    {"fast_screen", fast_result ? fast_result->to_json() : nlohmann::json(nullptr)},
                    {"prelinear_fast_repair", prelinear_fast_repair ? prelinear_fast_repair->to_json() : nlohmann::json(nullptr)},
                    {"rolling_fast_repair", rolling_fast_repair ? rolling_fast_repair->to_json() : nlohmann::json(nullptr)},
                    {"rolling_corrective_seed_label", rolling_corrective_seed_label ? nlohmann::json(*rolling_corrective_seed_label) : nlohmann::json(nullptr)},
                    {"corrective_seed_fast_repair_attempts", corrective_seed_fast_repair_attempts},
                    {"linearized_attempts", linearized_attempts},
                    {"resident_parametric_model", false},
                    {"acceptable_termination_enabled", false},
                    {"resident_model_created", false},
                    {"model_preparation_wall_seconds", 0.0},
                    {"solve", gravityx::solve_result_to_json(solve, false)},
                    {"validation", validation.to_json()},
                };
                complete(std::move(output), solve.state);
                std::cout << nlohmann::json({
                    {"output", output_path},
                    {"success", true},
                    {"label", match->label},
                    {"status", 0},
                    {"objective", objective},
                    {"wall_seconds", linearized_wall_seconds},
                    {"iterations", linearized_iterations},
                    {"linearized_rounds", round},
                    {"max_residual", validation.max_residual},
                    {"solution_method", "highs_sequential_linearized_contingency"},
                }).dump(2) << '\n';
                return true;
            }

            std::optional<gravityx::FastPowerFlowResult> nonlinear;
            if (fast_power_flow) {
                nonlinear = fast_power_flow->solve(*match, candidate);
            } else {
                gravityx::FastContingencyPowerFlow nonlinear_solver(
                    data, base.state, base.commitment);
                nonlinear = nonlinear_solver.solve(*match, candidate);
            }
            linearized_wall_seconds += nonlinear->wall_seconds;
            linearized_iterations += std::max(0, nonlinear->solve.iterations);
            attempt["nonlinear_repair"] = nonlinear->to_json();

            auto nonlinear_candidate = nonlinear->solve.state;
            double nonlinear_objective = best_objective;
            auto nonlinear_validation = nonlinear->validation;
            if (!nonlinear_candidate.vm.empty()) {
                nonlinear_objective =
                    gravityx::rebuild_contingency_state_derived_fields(
                        data, base.state, base.commitment, *match,
                        nonlinear_candidate);
                nonlinear_validation = gravityx::validate_state(
                    data, gravityx::ModelMode::ContingencySoft,
                    nonlinear_candidate, base.commitment, context);
                attempt["nonlinear_exact_objective"] = nonlinear_objective;
                attempt["nonlinear_validation"] =
                    nonlinear_validation.to_json();
                if (nonlinear_validation.max_residual <
                    best_validation.max_residual) {
                    best_state = nonlinear_candidate;
                    best_objective = nonlinear_objective;
                    best_validation = nonlinear_validation;
                }
                if (data.buses.size() >= 16000) {
                    attempt["new_dynamic_security_branches"] =
                        collect_violated_security_branches(
                            nonlinear_candidate);
                    attempt["active_dynamic_security_branch_positions_after_repair"] =
                        dynamic_security_branches;
                }
            }
            const char* linearized_log = std::getenv("GRAVITYX_HIGHS_LOG");
            if (linearized_log != nullptr &&
                std::string(linearized_log) != "0") {
                std::cerr << "GRAVITYX_LINEARIZED_ROUND_RESULT "
                          << nlohmann::json({
                                 {"label", match->label},
                                 {"round", round},
                                 {"reference_source", reference_source},
                                 {"linear_status", linear.status},
                                 {"linearized_balance_slack",
                                  linearized_balance_slack},
                                 {"linear_validation", validation.to_json()},
                                 {"linear_max_residual", validation.max_residual},
                                 {"linear_worst_category", validation.worst_category},
                                 {"linear_worst_identity", validation.worst_identity},
                                 {"nonlinear_converged", nonlinear->converged},
                                 {"nonlinear_feasible", nonlinear->feasible},
                                 {"nonlinear_failure_reason", nonlinear->failure_reason},
                                 {"nonlinear_validation",
                                  nonlinear_validation.to_json()},
                                 {"nonlinear_max_residual",
                                  nonlinear_validation.max_residual},
                                 {"nonlinear_worst_category",
                                  nonlinear_validation.worst_category},
                                 {"nonlinear_worst_identity",
                                  nonlinear_validation.worst_identity},
                                 {"nonlinear_wall_seconds", nonlinear->wall_seconds},
                                 {"nonlinear_iterations",
                                  nonlinear->solve.iterations},
                             }).dump()
                          << '\n';
            }
            linearized_attempts.push_back(std::move(attempt));
            if (!nonlinear_candidate.vm.empty() &&
                nonlinear_validation.max_residual <= 1e-5) {
                gravityx::SolveResult solve;
                solve.status = 0;
                solve.objective = nonlinear_objective;
                solve.wall_seconds = linearized_wall_seconds;
                solve.iterations = linearized_iterations;
                solve.state = std::move(nonlinear_candidate);
                nlohmann::json output = {
                    {"success", true},
                    {"solver_status_success", true},
                    {"accepted_feasible_nonconverged", false},
                    {"label", match->label},
                    {"type", match->type == gravityx::ContingencyType::Generator ? "gen" : "branch"},
                    {"source_index", match->source_index},
                    {"component_position", match->component},
                    {"solution_method", "highs_linearized_contingency_plus_fast_newton"},
                    {"precomputed_fast_screen_reference", reused_fast_screen_reference},
                    {"fast_power_flow_screen", fast_power_flow != nullptr},
                    {"fast_screen", fast_result ? fast_result->to_json() : nlohmann::json(nullptr)},
                    {"prelinear_fast_repair", prelinear_fast_repair ? prelinear_fast_repair->to_json() : nlohmann::json(nullptr)},
                    {"rolling_fast_repair", rolling_fast_repair ? rolling_fast_repair->to_json() : nlohmann::json(nullptr)},
                    {"rolling_corrective_seed_label", rolling_corrective_seed_label ? nlohmann::json(*rolling_corrective_seed_label) : nlohmann::json(nullptr)},
                    {"corrective_seed_fast_repair_attempts", corrective_seed_fast_repair_attempts},
                    {"linearized_attempts", linearized_attempts},
                    {"resident_parametric_model", false},
                    {"acceptable_termination_enabled", false},
                    {"resident_model_created", false},
                    {"model_preparation_wall_seconds", 0.0},
                    {"solve", gravityx::solve_result_to_json(solve, false)},
                    {"validation", nonlinear_validation.to_json()},
                };
                complete(std::move(output), solve.state);
                std::cout << nlohmann::json({
                    {"output", output_path},
                    {"success", true},
                    {"label", match->label},
                    {"status", 0},
                    {"objective", nonlinear_objective},
                    {"wall_seconds", linearized_wall_seconds},
                    {"iterations", linearized_iterations},
                    {"linearized_rounds", round},
                    {"max_residual", nonlinear_validation.max_residual},
                    {"solution_method", "highs_linearized_contingency_plus_fast_newton"},
                }).dump(2) << '\n';
                return true;
            }
            if (!nonlinear_candidate.vm.empty() &&
                nonlinear_validation.max_residual < validation.max_residual) {
                reference = std::move(nonlinear_candidate);
                reference_source = "iterated_nonlinear_repair";
            } else {
                reference = std::move(candidate);
                reference_source = "iterated_linear_seed";
            }
        }
        linearized_seed = best_state;
        }
        if (linearized_only) {
            gravityx::SolveResult solve;
            solve.status = last_model_status;
            solve.objective = best_objective;
            solve.wall_seconds = linearized_wall_seconds;
            solve.iterations = linearized_iterations;
            solve.state = std::move(best_state);
            nlohmann::json output = {
                {"success", false},
                {"solver_status_success", false},
                {"accepted_feasible_nonconverged", false},
                {"label", match->label},
                {"type", match->type == gravityx::ContingencyType::Generator ? "gen" : "branch"},
                {"source_index", match->source_index},
                {"component_position", match->component},
                {"solution_method", "highs_sequential_linearized_contingency_failed"},
                {"precomputed_fast_screen_reference", reused_fast_screen_reference},
                {"linearized_status", last_status},
                {"fast_power_flow_screen", fast_power_flow != nullptr},
                {"fast_screen", fast_result ? fast_result->to_json() : nlohmann::json(nullptr)},
                {"prelinear_fast_repair", prelinear_fast_repair ? prelinear_fast_repair->to_json() : nlohmann::json(nullptr)},
                {"rolling_fast_repair", rolling_fast_repair ? rolling_fast_repair->to_json() : nlohmann::json(nullptr)},
                {"rolling_corrective_seed_label", rolling_corrective_seed_label ? nlohmann::json(*rolling_corrective_seed_label) : nlohmann::json(nullptr)},
                {"corrective_seed_fast_repair_attempts", corrective_seed_fast_repair_attempts},
                {"linearized_attempts", linearized_attempts},
                {"resident_parametric_model", false},
                {"acceptable_termination_enabled", false},
                {"resident_model_created", false},
                {"model_preparation_wall_seconds", 0.0},
                {"solve", gravityx::solve_result_to_json(solve, false)},
                {"validation", best_validation.to_json()},
            };
            complete(std::move(output), solve.state, true);
            return false;
        }
    }

    const auto preparation_start = std::chrono::steady_clock::now();
    std::unique_ptr<gravityx::AcModel> fresh_model;
    gravityx::AcModel* model = nullptr;
    bool resident_model_created = false;
    if (reusable_model) {
        if (!*reusable_model) {
            *reusable_model = std::make_unique<gravityx::AcModel>(
                data, gravityx::ModelMode::ContingencySoft,
                base.commitment, context, true, acceptable_termination);
            resident_model_created = true;
        } else {
            (*reusable_model)->set_contingency(context);
        }
        model = reusable_model->get();
    } else {
        fresh_model = std::make_unique<gravityx::AcModel>(
            data, gravityx::ModelMode::ContingencySoft, base.commitment, context);
        model = fresh_model.get();
    }
    if (linearized_seed) {
        model->initialize_from(*linearized_seed);
    }
    const auto preparation_finish = std::chrono::steady_clock::now();
    const double preparation_seconds = std::chrono::duration<double>(
        preparation_finish - preparation_start).count();
    const auto solve = model->solve(print_level, 1e-6);
    const auto validation = gravityx::validate_state(
        data, gravityx::ModelMode::ContingencySoft,
        solve.state, base.commitment, context);
    const bool solver_status_success = solve.status == 0 || solve.status == 1;
    const bool success = gravityx::validated_candidate_is_feasible(
        solve, validation, 1e-5);
    const bool accepted_feasible_nonconverged = success && !solver_status_success;
    nlohmann::json output = {
        {"success", success},
        {"solver_status_success", solver_status_success},
        {"accepted_feasible_nonconverged", accepted_feasible_nonconverged},
        {"label", match->label},
        {"type", match->type == gravityx::ContingencyType::Generator ? "gen" : "branch"},
        {"source_index", match->source_index},
        {"component_position", match->component},
        {"solution_method", "ipopt_corrective_fallback"},
        {"precomputed_fast_screen_reference", reused_fast_screen_reference},
        {"fast_power_flow_screen", fast_power_flow != nullptr},
        {"fast_screen", fast_result ? fast_result->to_json() : nlohmann::json(nullptr)},
        {"prelinear_fast_repair", prelinear_fast_repair ? prelinear_fast_repair->to_json() : nlohmann::json(nullptr)},
        {"rolling_fast_repair", rolling_fast_repair ? rolling_fast_repair->to_json() : nlohmann::json(nullptr)},
        {"rolling_corrective_seed_label", rolling_corrective_seed_label ? nlohmann::json(*rolling_corrective_seed_label) : nlohmann::json(nullptr)},
        {"corrective_seed_fast_repair_attempts", corrective_seed_fast_repair_attempts},
        {"linearized_attempts", linearized_attempts},
        {"resident_parametric_model", reusable_model != nullptr},
        {"acceptable_termination_enabled", acceptable_termination},
        {"resident_model_created", resident_model_created},
        {"model_preparation_wall_seconds", preparation_seconds},
        {"solve", gravityx::solve_result_to_json(solve, false)},
        {"validation", validation.to_json()},
    };
    complete(std::move(output), solve.state);
    std::cout << nlohmann::json({
        {"output", output_path},
        {"success", success},
        {"solver_status_success", solver_status_success},
        {"accepted_feasible_nonconverged", accepted_feasible_nonconverged},
        {"label", match->label},
        {"status", solve.status},
        {"objective", solve.objective},
        {"wall_seconds", solve.wall_seconds},
        {"iterations", solve.iterations},
        {"resident_reoptimization", solve.resident_reoptimization},
        {"model_preparation_wall_seconds", preparation_seconds},
        {"max_residual", validation.max_residual},
    }).dump(2) << '\n';
    return success;
}

int solve_contingency(
    const std::string& case_path,
    const std::string& base_result_path,
    const std::string& label,
    const std::string& output_path,
    int print_level) {
    reject_onedrive(case_path);
    reject_onedrive(base_result_path);
    const auto data = gravityx::CaseData::load(case_path);
    const auto base = load_base_point(base_result_path);
    return solve_loaded_contingency(
        data, base, label, output_path, print_level) ? 0 : 1;
}

int screen_all_contingencies(
    const std::string& case_path,
    const std::string& base_result_path,
    const std::string& output_path) {
    reject_onedrive(case_path);
    reject_onedrive(base_result_path);
    reject_onedrive(output_path);
    const auto data = gravityx::CaseData::load(case_path);
    const auto base = load_base_point(base_result_path);
    gravityx::FastContingencyPowerFlow screen(
        data, base.state, base.commitment);
    nlohmann::json records = nlohmann::json::array();
    int feasible_count = 0;
    const auto started = std::chrono::steady_clock::now();
    for (const auto& contingency : data.contingencies) {
        const auto result = screen.solve(contingency);
        auto record = result.to_json();
        record["label"] = contingency.label;
        record["type"] = contingency.type == gravityx::ContingencyType::Generator
            ? "gen" : "branch";
        record["source_index"] = contingency.source_index;
        records.push_back(std::move(record));
        feasible_count += result.feasible ? 1 : 0;
    }
    const double wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    write_json_file(output_path, {
        {"contingency_count", data.contingencies.size()},
        {"feasible_count", feasible_count},
        {"fallback_count", static_cast<int>(data.contingencies.size()) - feasible_count},
        {"wall_seconds", wall_seconds},
        {"records", records},
    });
    std::cout << nlohmann::json({
        {"output", output_path},
        {"contingency_count", data.contingencies.size()},
        {"feasible_count", feasible_count},
        {"fallback_count", static_cast<int>(data.contingencies.size()) - feasible_count},
        {"wall_seconds", wall_seconds},
    }).dump(2) << '\n';
    return 0;
}

int solve_contingency_batch(
    const std::string& case_path,
    const std::string& base_result_path,
    const std::string& manifest_path,
    int print_level) {
    reject_onedrive(case_path);
    reject_onedrive(base_result_path);
    reject_onedrive(manifest_path);
    const auto data = gravityx::CaseData::load(case_path);
    const auto base = load_base_point(base_result_path);
    const auto manifest = read_json_file(manifest_path);
    if (!manifest.contains("tasks") || !manifest.at("tasks").is_array() ||
        manifest.at("tasks").empty()) {
        throw std::runtime_error("contingency batch manifest has no tasks");
    }
    int completed = 0;
    for (const auto& task : manifest.at("tasks")) {
        const auto label = task.at("label").get<std::string>();
        const auto output_path = task.at("output_path").get<std::string>();
        if (!solve_loaded_contingency(
                data, base, label, output_path, print_level)) {
            return 1;
        }
        ++completed;
        std::cout << nlohmann::json({
            {"batch_completed", completed},
            {"batch_size", manifest.at("tasks").size()},
            {"label", label},
        }).dump() << '\n';
    }
    return 0;
}

int run_contingency_worker(
    const std::string& case_path,
    const std::string& base_result_path,
    int print_level,
    bool reusable_model,
    bool acceptable_termination,
    bool fast_power_flow_screen,
    bool fast_only,
    bool linearized_fallback,
    bool linearized_only) {
    reject_onedrive(case_path);
    reject_onedrive(base_result_path);
    const auto data = gravityx::CaseData::load(case_path);
    const auto base = load_base_point(base_result_path);
    const gravityx::GoSolutionWriter solution_writer(data);
    if (fast_only && !fast_power_flow_screen) {
        throw std::runtime_error(
            "fast-only contingency worker requires fast-pf mode");
    }
    if (linearized_only && !linearized_fallback) {
        throw std::runtime_error(
            "linearized-only contingency worker requires linearized mode");
    }
    std::unique_ptr<gravityx::AcModel> resident_model;
    std::unique_ptr<gravityx::FastContingencyPowerFlow> fast_power_flow;
    if (fast_power_flow_screen) {
        gravityx::FastPowerFlowOptions fast_options;
        fast_options.fixed_jacobian_screen_only =
            fast_only && data.buses.size() >= 16000;
        const char* fast_diagnostics =
            std::getenv("GRAVITYX_FAST_PF_DIAGNOSTICS");
        fast_options.capture_diagnostics = fast_diagnostics != nullptr &&
            std::string(fast_diagnostics) != "0";
        fast_power_flow = std::make_unique<gravityx::FastContingencyPowerFlow>(
            data, base.state, base.commitment, fast_options);
    }
    std::optional<gravityx::AcState> rolling_corrective_seed;
    std::string rolling_corrective_seed_label;
    std::vector<CorrectiveSeed> corrective_seed_bank;
    std::cout << "GRAVITYX_WORKER_READY" << std::endl;
    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.empty()) {
            continue;
        }
        const auto task = nlohmann::json::parse(line);
        if (task.value("stop", false)) {
            return 0;
        }
        const auto label = task.at("label").get<std::string>();
        const auto output_path = task.at("output_path").get<std::string>();
        const std::optional<std::string> fallback_output_path =
            task.contains("fallback_output_path")
            ? std::make_optional(
                  task.at("fallback_output_path").get<std::string>())
            : std::nullopt;
        const bool remove_output_after_result =
            task.value("remove_output_after_result", false);
        if (fallback_output_path) {
            reject_onedrive(*fallback_output_path);
        }
        const std::optional<std::string> solution_path =
            task.contains("solution_path")
            ? std::make_optional(
                  task.at("solution_path").get<std::string>())
            : std::nullopt;
        if (solution_path) {
            reject_onedrive(*solution_path);
        }
        std::optional<gravityx::AcState> precomputed_fast_state;
        if (task.contains("fast_screen_path")) {
            const auto fast_screen_path =
                task.at("fast_screen_path").get<std::string>();
            const auto fast_screen_json = read_json_file(fast_screen_path);
            if (fast_screen_json.value("label", std::string()) != label ||
                !fast_screen_json.value("screen_completed", false) ||
                !fast_screen_json.value("requires_exact_fallback", false) ||
                !fast_screen_json.contains("solve") ||
                !fast_screen_json.at("solve").contains("state")) {
                throw std::runtime_error(
                    "invalid precomputed fast-screen result for " + label);
            }
            precomputed_fast_state = gravityx::ac_state_from_json(
                fast_screen_json.at("solve").at("state"));
        }
        std::optional<ContingencyComputation> completed_computation;
        const bool success = solve_loaded_contingency(
            data, base, label, output_path, print_level,
            reusable_model ? &resident_model : nullptr,
            acceptable_termination, fast_power_flow.get(), fast_only,
            linearized_fallback, linearized_only,
            precomputed_fast_state ? &*precomputed_fast_state : nullptr,
            rolling_corrective_seed ? &*rolling_corrective_seed : nullptr,
            rolling_corrective_seed ? &rolling_corrective_seed_label : nullptr,
            (linearized_fallback || fast_only)
                ? &corrective_seed_bank : nullptr,
            &completed_computation, !remove_output_after_result);
        double result_read_seconds = 0.0;
        if (success && !completed_computation) {
            throw std::runtime_error(
                "contingency worker produced no in-memory result for " + label);
        }
        bool solution_written = false;
        double solution_write_seconds = 0.0;
        if (solution_path && completed_computation &&
            completed_computation->result.value("success", false)) {
            const auto contingency = std::find_if(
                data.contingencies.begin(), data.contingencies.end(),
                [&](const gravityx::Contingency& item) {
                    return item.label == label;
                });
            if (contingency == data.contingencies.end()) {
                throw std::runtime_error(
                    "cannot write solution for unknown contingency " + label);
            }
            const auto solution_write_start =
                std::chrono::steady_clock::now();
            solution_writer.write(
                *solution_path, completed_computation->state,
                base.commitment,
                &*contingency);
            solution_write_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - solution_write_start)
                .count();
            solution_written = true;
        }
        bool fallback_result_persisted = false;
        double fallback_result_persist_seconds = 0.0;
        if (fallback_output_path && completed_computation &&
            completed_computation->result.value(
                "requires_exact_fallback", false)) {
            const auto persist_start = std::chrono::steady_clock::now();
            auto fallback_result = completed_computation->result;
            fallback_result["solve"]["state"] =
                gravityx::ac_submission_state_to_json(
                    completed_computation->state);
            write_json_file(*fallback_output_path, fallback_result);
            fallback_result_persist_seconds =
                std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - persist_start)
                    .count();
            fallback_result_persisted = true;
        }
        bool rolling_corrective_seed_updated = false;
        if (success && (linearized_fallback || fast_only) &&
            data.buses.size() >= 16000 && completed_computation) {
            const auto& result_json = completed_computation->result;
            if (result_json.value("success", false)) {
                auto verified_corrective_state = completed_computation->state;
                rolling_corrective_seed = verified_corrective_state;
                rolling_corrective_seed_label = label;
                corrective_seed_bank.erase(
                    std::remove_if(
                        corrective_seed_bank.begin(),
                        corrective_seed_bank.end(),
                        [&](const CorrectiveSeed& seed) {
                            return seed.label == label;
                        }),
                    corrective_seed_bank.end());
                corrective_seed_bank.insert(
                    corrective_seed_bank.begin(),
                    CorrectiveSeed{label, std::move(verified_corrective_state)});
                constexpr std::size_t kMaximumCorrectiveSeedBankSize = 16;
                if (corrective_seed_bank.size() >
                    kMaximumCorrectiveSeedBankSize) {
                    corrective_seed_bank.resize(
                        kMaximumCorrectiveSeedBankSize);
                }
                rolling_corrective_seed_updated = true;
            }
        }
        nlohmann::json result_summary = nullptr;
        if (completed_computation) {
            const auto& result_json = completed_computation->result;
            const auto& solve_json = result_json.at("solve");
            const auto& validation_json = result_json.at("validation");
            std::string fast_screen_failure_reason;
            if (result_json.contains("fast_screen") &&
                result_json.at("fast_screen").is_object()) {
                fast_screen_failure_reason = result_json.at("fast_screen")
                    .value("failure_reason", std::string());
            }
            result_summary = {
                {"success", result_json.value("success", false)},
                {"screen_completed",
                 result_json.value("screen_completed", false)},
                {"requires_exact_fallback",
                 result_json.value("requires_exact_fallback", false)},
                {"solver_status_success",
                 result_json.value("solver_status_success", false)},
                {"accepted_feasible_nonconverged",
                 result_json.value(
                     "accepted_feasible_nonconverged", false)},
                {"solution_method",
                 result_json.value("solution_method", std::string())},
                {"fast_power_flow_screen",
                 result_json.value("fast_power_flow_screen", false)},
                {"rolling_corrective_seed_label",
                 result_json.contains("rolling_corrective_seed_label")
                     ? result_json.at("rolling_corrective_seed_label")
                     : nlohmann::json(nullptr)},
                {"precomputed_fast_screen_reference",
                 result_json.value(
                     "precomputed_fast_screen_reference", false)},
                {"model_preparation_wall_seconds",
                 result_json.value(
                     "model_preparation_wall_seconds", 0.0)},
                {"solve", {
                    {"status", solve_json.value("status", -1)},
                    {"objective", solve_json.value("objective", 0.0)},
                    {"wall_seconds",
                     solve_json.value("wall_seconds", 0.0)},
                    {"iterations", solve_json.value("iterations", -1)},
                    {"resident_reoptimization",
                     solve_json.value("resident_reoptimization", false)},
                    {"acceptable_termination_enabled",
                     solve_json.value(
                         "acceptable_termination_enabled", false)},
                }},
                {"validation", {
                    {"max_residual",
                     validation_json.value("max_residual", 0.0)},
                }},
                {"fast_screen", {
                    {"failure_reason", fast_screen_failure_reason},
                }},
            };
        }
        bool transient_output_removed = false;
        if (remove_output_after_result) {
            std::error_code remove_error;
            if (std::filesystem::exists(output_path, remove_error)) {
                std::filesystem::remove(output_path, remove_error);
            }
            if (remove_error) {
                throw std::runtime_error(
                    "cannot remove transient contingency result " +
                    output_path + ": " + remove_error.message());
            }
            transient_output_removed =
                !std::filesystem::exists(output_path, remove_error);
            if (remove_error) {
                throw std::runtime_error(
                    "cannot verify transient contingency cleanup " +
                    output_path + ": " + remove_error.message());
            }
        }
        std::cout << "GRAVITYX_TASK_RESULT " << nlohmann::json({
            {"label", label},
            {"success", success},
            {"precomputed_fast_screen_reference",
             precomputed_fast_state.has_value()},
            {"rolling_corrective_seed_updated",
             rolling_corrective_seed_updated},
            {"rolling_corrective_seed_label",
             rolling_corrective_seed
                 ? nlohmann::json(rolling_corrective_seed_label)
                 : nlohmann::json(nullptr)},
            {"corrective_seed_bank_size", corrective_seed_bank.size()},
            {"solution_written", solution_written},
            {"result_read_seconds", result_read_seconds},
            {"solution_write_seconds", solution_write_seconds},
            {"fallback_result_persisted", fallback_result_persisted},
            {"fallback_result_persist_seconds",
             fallback_result_persist_seconds},
            {"transient_output_removed", transient_output_removed},
            {"result_summary", result_summary},
        }).dump() << std::endl;
        if (!success) {
            return 1;
        }
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "smoke") {
            return run_smoke();
        }
        if (argc == 2 && std::string(argv[1]) == "component-tests") {
            return run_component_tests();
        }
        if (argc == 2 && std::string(argv[1]) == "parallel-circuit-test") {
            return run_parallel_circuit_regression();
        }
        if (argc == 3 && std::string(argv[1]) == "inspect") {
            return run_inspect(argv[2]);
        }
        if ((argc == 3 || argc == 4) &&
            (std::string(argv[1]) == "solve-base" || std::string(argv[1]) == "solve-relax")) {
            const int print_level = argc == 4 ? std::stoi(argv[3]) : 0;
            return run_ac_model(argv[1], argv[2], print_level);
        }
        if ((argc == 3 || argc == 4) && std::string(argv[1]) == "run-ibr") {
            const int print_level = argc == 4 ? std::stoi(argv[3]) : 0;
            return run_ibr(argv[2], print_level);
        }
        if ((argc >= 4 && argc <= 6) && std::string(argv[1]) == "run-ibr-json") {
            const int print_level = argc == 5 ? std::stoi(argv[4]) : 0;
            bool source_status_only = false;
            for (int i = 4; i < argc; ++i) {
                const std::string option = argv[i];
                if (option == "source-only") {
                    source_status_only = true;
                } else if (i == 4) {
                    continue;
                } else {
                    throw std::runtime_error("unknown run-ibr-json option: " + option);
                }
            }
            return run_ibr_json(argv[2], argv[3], print_level, source_status_only);
        }
        if ((argc >= 4 && argc <= 6) &&
            std::string(argv[1]) == "validated-source-base-json") {
            bool allow_exact_fallback = true;
            bool allow_large_base_newton_restart = true;
            for (int i = 4; i < argc; ++i) {
                const std::string option = argv[i];
                if (option == "fast-only") {
                    allow_exact_fallback = false;
                } else if (option == "robust-contingency-seed") {
                    allow_large_base_newton_restart = false;
                } else {
                    throw std::runtime_error(
                        "unknown validated-source-base-json option: " + option);
                }
            }
            return run_validated_source_base_json(
                argv[2], argv[3], allow_exact_fallback,
                allow_large_base_newton_restart);
        }
        if ((argc == 6 || argc == 7) && std::string(argv[1]) == "solve-contingency") {
            const int print_level = argc == 7 ? std::stoi(argv[6]) : 0;
            return solve_contingency(argv[2], argv[3], argv[4], argv[5], print_level);
        }
        if (argc == 5 && std::string(argv[1]) == "screen-all-contingencies") {
            return screen_all_contingencies(argv[2], argv[3], argv[4]);
        }
        if ((argc == 5 || argc == 6) &&
            std::string(argv[1]) == "solve-contingency-batch") {
            const int print_level = argc == 6 ? std::stoi(argv[5]) : 0;
            return solve_contingency_batch(argv[2], argv[3], argv[4], print_level);
        }
        if ((argc >= 4 && argc <= 11) &&
            std::string(argv[1]) == "contingency-worker") {
            const int print_level = argc >= 5 ? std::stoi(argv[4]) : 0;
            bool reusable_model = false;
            bool acceptable_termination = false;
            bool fast_power_flow_screen = false;
            bool fast_only = false;
            bool linearized_fallback = false;
            bool linearized_only = false;
            for (int i = 5; i < argc; ++i) {
                const std::string option = argv[i];
                if (option == "resident") {
                    reusable_model = true;
                } else if (option == "acceptable") {
                    acceptable_termination = true;
                } else if (option == "fast-pf") {
                    fast_power_flow_screen = true;
                } else if (option == "fast-only") {
                    fast_only = true;
                } else if (option == "linearized") {
                    linearized_fallback = true;
                } else if (option == "linearized-only") {
                    linearized_fallback = true;
                    linearized_only = true;
                } else {
                    throw std::runtime_error(
                        "unknown contingency-worker option: " + option);
                }
            }
            if (acceptable_termination && !reusable_model) {
                throw std::runtime_error(
                    "acceptable contingency termination requires resident mode");
            }
            return run_contingency_worker(
                argv[2], argv[3], print_level, reusable_model,
                acceptable_termination, fast_power_flow_screen, fast_only,
                linearized_fallback, linearized_only);
        }
        std::cerr << "usage:\n"
                  << "  gravityx_go2 smoke\n"
                  << "  gravityx_go2 component-tests\n"
                  << "  gravityx_go2 parallel-circuit-test\n"
                  << "  gravityx_go2 inspect CASE.json\n"
                  << "  gravityx_go2 solve-base CASE.json [print-level]\n"
                  << "  gravityx_go2 solve-relax CASE.json [print-level]\n"
                  << "  gravityx_go2 run-ibr CASE.json [print-level]\n"
                  << "  gravityx_go2 run-ibr-json CASE.json OUTPUT.json [print-level]\n"
                  << "  gravityx_go2 validated-source-base-json CASE.json OUTPUT.json [fast-only] [robust-contingency-seed]\n"
                  << "  gravityx_go2 solve-contingency CASE.json BASE.json LABEL OUTPUT.json [print-level]\n"
                  << "  gravityx_go2 screen-all-contingencies CASE.json BASE.json OUTPUT.json\n"
                  << "  gravityx_go2 solve-contingency-batch CASE.json BASE.json MANIFEST.json [print-level]\n"
                  << "  gravityx_go2 contingency-worker CASE.json BASE.json [print-level] [resident] [acceptable] [fast-pf] [fast-only] [linearized] [linearized-only]\n";
        return 2;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
