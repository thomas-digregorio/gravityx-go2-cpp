#pragma once

#include "gravityx/case_data.hpp"

#include <gravity/solver.h>

#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace gravityx {

enum class ModelMode {
    BaseSoft,
    UnitCommitmentRelaxation,
    ContingencySoft,
};

struct AcState {
    std::vector<double> vm;
    std::vector<double> va;
    std::vector<double> pg;
    std::vector<double> qg;
    std::vector<double> demand_factor;
    std::vector<double> pf;
    std::vector<double> qf;
    std::vector<double> pt;
    std::vector<double> qt;
    std::vector<double> sm_slack;
    std::vector<double> p_delta;
    std::vector<double> q_delta;
    std::vector<double> commitment;
    std::vector<double> startup;
    std::vector<double> shutdown;
    std::vector<double> gen_lambda;
    std::vector<double> load_lambda;
};

struct SolveResult {
    int status{};
    double objective{};
    double wall_seconds{};
    AcState state;
};

struct ContingencyContext {
    int outaged_generator{-1};
    int outaged_branch{-1};
    AcState base_state;
};

class AcModel {
public:
    AcModel(
        const CaseData& data,
        ModelMode mode,
        std::vector<int> fixed_status = {},
        std::optional<ContingencyContext> contingency = std::nullopt);

    SolveResult solve(int print_level = 0, double tolerance = 1e-6);
    void initialize_from(const AcState& state);
    void set_commitment_bound(int generator, int status);
    void set_commitment_start(int generator, double value);

private:
    const CaseData& data_;
    ModelMode mode_;
    std::vector<int> fixed_status_;
    std::optional<ContingencyContext> contingency_;
    gravity::Model model_;

    std::unique_ptr<gravity::var<double>> vm_;
    std::unique_ptr<gravity::var<double>> va_;
    std::unique_ptr<gravity::var<double>> pg_;
    std::unique_ptr<gravity::var<double>> qg_;
    std::unique_ptr<gravity::var<double>> demand_;
    std::unique_ptr<gravity::var<double>> pf_;
    std::unique_ptr<gravity::var<double>> qf_;
    std::unique_ptr<gravity::var<double>> pt_;
    std::unique_ptr<gravity::var<double>> qt_;
    std::unique_ptr<gravity::var<double>> sm_slack_;
    std::unique_ptr<gravity::var<double>> p_delta_;
    std::unique_ptr<gravity::var<double>> q_delta_;
    std::unique_ptr<gravity::var<double>> commitment_;
    std::unique_ptr<gravity::var<double>> startup_;
    std::unique_ptr<gravity::var<double>> shutdown_;
    std::unique_ptr<gravity::var<double>> gen_lambda_;
    std::unique_ptr<gravity::var<double>> load_lambda_;

    std::vector<std::vector<PwlPoint>> gen_points_;
    std::vector<std::vector<PwlPoint>> load_points_;
    std::vector<int> gen_lambda_offset_;
    std::vector<int> load_lambda_offset_;

    void build_variables();
    void build_constraints_and_objective();
    void initialize_source_point();
    AcState capture_state() const;
};

}  // namespace gravityx
