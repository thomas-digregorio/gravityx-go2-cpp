#include "gravityx/ac_model.hpp"
#include "gravityx/algorithm.hpp"
#include "gravityx/case_data.hpp"
#include "gravityx/state_io.hpp"
#include "gravityx/validation.hpp"

#include <gravity/solver.h>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
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

    gravityx::AcModel model(data, gravityx::ModelMode::BaseSoft, {1});
    const auto solve = model.solve(0, 1e-7);
    const auto validation = gravityx::validate_state(
        data, gravityx::ModelMode::BaseSoft, solve.state, {1});
    if ((solve.status != 0 && solve.status != 1) || !std::isfinite(solve.objective) ||
        validation.max_residual > 1e-5) {
        throw std::runtime_error("parallel-circuit symbolic-DAG regression failed");
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

int run_ibr_json(const std::string& path, const std::string& output_path, int print_level) {
    reject_onedrive(path);
    reject_onedrive(output_path);
    const auto data = gravityx::CaseData::load(path);
    gravityx::IbrOptions options;
    options.print_level = print_level;
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
    int print_level) {
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

    gravityx::AcModel model(
        data, gravityx::ModelMode::ContingencySoft, base.commitment, context);
    const auto solve = model.solve(print_level, 1e-6);
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
        {"solve", gravityx::solve_result_to_json(solve, true)},
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
    int print_level) {
    reject_onedrive(case_path);
    reject_onedrive(base_result_path);
    const auto data = gravityx::CaseData::load(case_path);
    const auto base = load_base_point(base_result_path);
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
        const bool success = solve_loaded_contingency(
            data, base, label, output_path, print_level);
        std::cout << "GRAVITYX_TASK_RESULT " << nlohmann::json({
            {"label", label},
            {"success", success},
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
        if ((argc == 4 || argc == 5) && std::string(argv[1]) == "run-ibr-json") {
            const int print_level = argc == 5 ? std::stoi(argv[4]) : 0;
            return run_ibr_json(argv[2], argv[3], print_level);
        }
        if ((argc == 6 || argc == 7) && std::string(argv[1]) == "solve-contingency") {
            const int print_level = argc == 7 ? std::stoi(argv[6]) : 0;
            return solve_contingency(argv[2], argv[3], argv[4], argv[5], print_level);
        }
        if ((argc == 5 || argc == 6) &&
            std::string(argv[1]) == "solve-contingency-batch") {
            const int print_level = argc == 6 ? std::stoi(argv[5]) : 0;
            return solve_contingency_batch(argv[2], argv[3], argv[4], print_level);
        }
        if ((argc == 4 || argc == 5) &&
            std::string(argv[1]) == "contingency-worker") {
            const int print_level = argc == 5 ? std::stoi(argv[4]) : 0;
            return run_contingency_worker(argv[2], argv[3], print_level);
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
                  << "  gravityx_go2 solve-contingency CASE.json BASE.json LABEL OUTPUT.json [print-level]\n"
                  << "  gravityx_go2 solve-contingency-batch CASE.json BASE.json MANIFEST.json [print-level]\n"
                  << "  gravityx_go2 contingency-worker CASE.json BASE.json [print-level]\n";
        return 2;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
