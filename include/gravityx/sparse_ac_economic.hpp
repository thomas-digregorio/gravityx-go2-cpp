#pragma once

#include "gravityx/ac_model.hpp"
#include "gravityx/case_data.hpp"
#include "gravityx/validation.hpp"

#include <nlohmann/json.hpp>

#include <string>
#include <vector>

namespace gravityx {

struct SparseAcEconomicOptions {
    double time_limit_seconds{60.0};
    double tolerance{1e-7};
    double acceptable_tolerance{1e-5};
    int print_level{};
};

struct SparseAcEconomicResult {
    bool attempted{};
    bool solver_initialized{};
    bool candidate_returned{};
    bool candidate_verified{};
    bool best_intermediate_found{};
    bool improved{};
    int application_status{-99};
    int solver_return_status{-99};
    int iterations{-1};
    int variable_count{};
    int constraint_count{};
    int jacobian_nonzero_count{};
    int intermediate_callbacks{};
    int intermediate_iterates_retrieved{};
    int intermediate_verified_candidates{};
    int intermediate_capture_failures{};
    int best_intermediate_iteration{-1};
    double wall_seconds{};
    double intermediate_callback_seconds{};
    double scaled_solver_objective{};
    double incumbent_objective{};
    double candidate_objective{};
    double initial_constraint_violation{};
    double candidate_constraint_violation{};
    double best_intermediate_objective{};
    double best_intermediate_max_residual{};
    std::string status;
    std::string selected_source{"incumbent"};
    std::string intermediate_capture_error;
    SolveResult selected;
    ValidationReport selected_validation;

    nlohmann::json to_json(bool include_state = false) const;
};

SparseAcEconomicResult solve_sparse_fixed_commitment_ac_economic(
    const CaseData& data,
    const std::vector<int>& commitment,
    const SolveResult& incumbent,
    const SparseAcEconomicOptions& options = {});

}  // namespace gravityx
