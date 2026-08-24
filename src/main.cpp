#include "gravityx/ac_model.hpp"
#include "gravityx/algorithm.hpp"
#include "gravityx/case_data.hpp"
#include "gravityx/fast_power_flow.hpp"
#include "gravityx/linearized_ac_seed.hpp"
#include "gravityx/state_io.hpp"
#include "gravityx/validation.hpp"

#include <gravity/solver.h>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>

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
        stream << value.dump(2) << '\n';
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

int run_component_tests() {
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
    feasible_validation.max_residual = 1e-8;
    nonconverged_feasible.objective = std::numeric_limits<double>::quiet_NaN();
    if (gravityx::validated_candidate_is_feasible(
            nonconverged_feasible, feasible_validation, 1e-5)) {
        throw std::runtime_error(
            "component test failed: nonfinite candidate was accepted");
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
    const auto linear_seed = gravityx::solve_linearized_ac_seed(
        data, source_base.solve.state, {1});
    if (!linear_seed.success) {
        throw std::runtime_error(
            "linearized AC seed regression failed: " + linear_seed.status);
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
    int initially_on = 0;
    int startup_eligible = 0;
    int shutdown_eligible = 0;
    for (const auto& branch : data.branches) {
        transformers += branch.transformer ? 1 : 0;
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
    bool allow_exact_fallback = true) {
    reject_onedrive(path);
    reject_onedrive(output_path);
    const auto data = gravityx::CaseData::load(path);
    std::vector<int> commitment;
    commitment.reserve(data.generators.size());
    for (const auto& generator : data.generators) {
        commitment.push_back(generator.status_prev);
    }
    const auto source = gravityx::build_validated_source_base(data, commitment);
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
        auto source_repair_json = repaired.to_json();
        source_repair_json["start"] = "source_voltage";
        repair_json.push_back(std::move(source_repair_json));
        wall_seconds += repaired.wall_seconds;
        if (repaired.validation.max_residual < selected_validation.max_residual) {
            selected_solve = repaired.solve;
            selected_validation = repaired.validation;
        }
        success = repaired.feasible;
        if (!success) {
            auto flat_state = source.solve.state;
            for (int i = 0; i < static_cast<int>(data.buses.size()); ++i) {
                flat_state.vm[i] = std::clamp(
                    1.0, data.buses[i].vmin, data.buses[i].vmax);
                flat_state.va[i] = 0.0;
            }
            gravityx::FastContingencyPowerFlow flat_repair(
                data, flat_state, commitment);
            const auto flat = flat_repair.solve_base();
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
        if (!success) {
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
    if (!success) {
        auto linear_reference = selected_solve.state;
        double linear_reference_residual = selected_validation.max_residual;
        for (int round = 1; round <= 4; ++round) {
            const auto linear = gravityx::solve_linearized_ac_seed(
                data, linear_reference, commitment);
            auto round_json = linear.to_json(false);
            round_json["round"] = round;
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
    const gravityx::AcState* precomputed_fast_state = nullptr) {
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

    std::optional<gravityx::FastPowerFlowResult> fast_result;
    const bool reused_fast_screen_reference = precomputed_fast_state != nullptr;
    if (precomputed_fast_state) {
        fast_result.emplace();
        fast_result->solve.state = *precomputed_fast_state;
        fast_result->failure_reason = "precomputed_fast_screen_reference";
    } else if (fast_power_flow) {
        fast_result = fast_power_flow->solve(*match);
        if (fast_result->feasible) {
            const std::string solution_method =
                fast_result->direct_candidate_selected
                ? "direct_base_state_outage_candidate"
                : "fast_newton_power_flow";
            const nlohmann::json output = {
                {"success", true},
                {"solver_status_success", true},
                {"accepted_feasible_nonconverged", false},
                {"label", match->label},
                {"type", match->type == gravityx::ContingencyType::Generator ? "gen" : "branch"},
                {"source_index", match->source_index},
                {"component_position", match->component},
                {"solution_method", solution_method},
                {"fast_power_flow_screen", true},
                {"fast_screen", fast_result->to_json()},
                {"resident_parametric_model", false},
                {"acceptable_termination_enabled", false},
                {"resident_model_created", false},
                {"model_preparation_wall_seconds", 0.0},
                {"solve", gravityx::solve_result_to_submission_json(fast_result->solve)},
                {"validation", fast_result->validation.to_json()},
            };
            write_json_file(output_path, output);
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
            const nlohmann::json output = {
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
                {"resident_parametric_model", false},
                {"acceptable_termination_enabled", false},
                {"resident_model_created", false},
                {"model_preparation_wall_seconds", 0.0},
                {"solve", gravityx::solve_result_to_submission_json(fast_result->solve)},
                {"validation", fast_result->validation.to_json()},
            };
            write_json_file(output_path, output);
            return true;
        }
    }

    std::optional<gravityx::AcState> linearized_seed;
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
        auto reference = best_state;
        double linearized_wall_seconds = 0.0;
        int linearized_iterations = 0;
        std::string last_status = "not_run";
        int last_model_status = -1;
        std::string reference_source = fast_result ? "fast_screen" : "base";
        constexpr int kMaximumLinearizedRounds = 3;
        for (int round = 1; round <= kMaximumLinearizedRounds; ++round) {
            auto linear = gravityx::solve_linearized_ac_seed(
                data, reference, base.commitment, 0.49, context);
            linearized_wall_seconds += linear.wall_seconds;
            linearized_iterations += std::max(0, linear.iterations);
            auto attempt = linear.to_json(false);
            attempt["round"] = round;
            attempt["reference_source"] = reference_source;
            if (!linear.success && data.buses.size() >= 8000 &&
                (linear.status == "Infeasible" ||
                 linear.status == "Unknown")) {
                attempt["projected_balance_retry_scheduled"] = true;
                linearized_attempts.push_back(std::move(attempt));
                linear = gravityx::solve_linearized_ac_seed(
                    data, reference, base.commitment, 0.49, context, true);
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
                const nlohmann::json output = {
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
                    {"linearized_attempts", linearized_attempts},
                    {"resident_parametric_model", false},
                    {"acceptable_termination_enabled", false},
                    {"resident_model_created", false},
                    {"model_preparation_wall_seconds", 0.0},
                    {"solve", gravityx::solve_result_to_submission_json(solve)},
                    {"validation", validation.to_json()},
                };
                write_json_file(output_path, output);
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
                const nlohmann::json output = {
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
                    {"linearized_attempts", linearized_attempts},
                    {"resident_parametric_model", false},
                    {"acceptable_termination_enabled", false},
                    {"resident_model_created", false},
                    {"model_preparation_wall_seconds", 0.0},
                    {"solve", gravityx::solve_result_to_submission_json(solve)},
                    {"validation", nonlinear_validation.to_json()},
                };
                write_json_file(output_path, output);
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
        if (linearized_only) {
            gravityx::SolveResult solve;
            solve.status = last_model_status;
            solve.objective = best_objective;
            solve.wall_seconds = linearized_wall_seconds;
            solve.iterations = linearized_iterations;
            solve.state = std::move(best_state);
            const nlohmann::json output = {
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
                {"linearized_attempts", linearized_attempts},
                {"resident_parametric_model", false},
                {"acceptable_termination_enabled", false},
                {"resident_model_created", false},
                {"model_preparation_wall_seconds", 0.0},
                {"solve", gravityx::solve_result_to_json(solve, true)},
                {"validation", best_validation.to_json()},
            };
            write_json_file(output_path, output);
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
    const nlohmann::json output = {
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
        {"linearized_attempts", linearized_attempts},
        {"resident_parametric_model", reusable_model != nullptr},
        {"acceptable_termination_enabled", acceptable_termination},
        {"resident_model_created", resident_model_created},
        {"model_preparation_wall_seconds", preparation_seconds},
        {"solve", gravityx::solve_result_to_submission_json(solve)},
        {"validation", validation.to_json()},
    };
    write_json_file(output_path, output);
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
        fast_power_flow = std::make_unique<gravityx::FastContingencyPowerFlow>(
            data, base.state, base.commitment);
    }
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
        const bool success = solve_loaded_contingency(
            data, base, label, output_path, print_level,
            reusable_model ? &resident_model : nullptr,
            acceptable_termination, fast_power_flow.get(), fast_only,
            linearized_fallback, linearized_only,
            precomputed_fast_state ? &*precomputed_fast_state : nullptr);
        std::cout << "GRAVITYX_TASK_RESULT " << nlohmann::json({
            {"label", label},
            {"success", success},
            {"precomputed_fast_screen_reference",
             precomputed_fast_state.has_value()},
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
        if ((argc == 4 || argc == 5) &&
            std::string(argv[1]) == "validated-source-base-json") {
            const bool allow_exact_fallback = argc == 4 ||
                std::string(argv[4]) != "fast-only";
            return run_validated_source_base_json(
                argv[2], argv[3], allow_exact_fallback);
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
                  << "  gravityx_go2 validated-source-base-json CASE.json OUTPUT.json [fast-only]\n"
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
