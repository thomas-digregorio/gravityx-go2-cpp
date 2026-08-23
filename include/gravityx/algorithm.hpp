#pragma once

#include "gravityx/ac_model.hpp"
#include "gravityx/case_data.hpp"
#include "gravityx/validation.hpp"

#include <nlohmann/json.hpp>

#include <vector>

namespace gravityx {

struct IbrOptions {
    int batch_count{4};
    double threshold{0.5};
    double tolerance{1e-6};
    double validation_tolerance{1e-5};
    int print_level{0};
};

struct IbrRound {
    int round{};
    std::vector<int> batch;
    std::vector<int> proposed_status;
    bool fallback_to_prior{};
    SolveResult solve;
    ValidationReport validation;

    nlohmann::json to_json() const;
};

struct IbrResult {
    bool success{};
    bool candidate_accepted{};
    double wall_seconds{};
    double switching_cost{};
    double candidate_proxy{};
    SolveResult base;
    ValidationReport base_validation;
    std::vector<IbrRound> rounds;
    SolveResult fixed_repair;
    ValidationReport fixed_validation;
    std::vector<int> commitment;
    AcState selected_state;

    nlohmann::json to_json(bool include_state = false) const;
};

bool validated_candidate_is_feasible(
    const SolveResult& result,
    const ValidationReport& validation,
    double tolerance);

IbrResult run_iterative_batch_rounding(const CaseData& data, const IbrOptions& options = {});

}  // namespace gravityx
