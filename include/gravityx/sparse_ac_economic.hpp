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
    bool pwl_epigraph{false};
    bool exact_hessian{false};
};

struct SparseAcEconomicResult {
    bool attempted{};
    bool solver_initialized{};
    bool candidate_returned{};
    bool candidate_verified{};
    bool incumbent_verified{};
    bool best_intermediate_found{};
    bool improved{};
    bool common_corrective_reference{};
    int application_status{-99};
    int solver_return_status{-99};
    int iterations{-1};
    int variable_count{};
    int constraint_count{};
    int jacobian_nonzero_count{};
    bool exact_hessian_enabled{};
    int hessian_nonzero_count{};
    int hessian_evaluations{};
    bool pwl_epigraph_enabled{};
    int pwl_epigraph_curve_count{};
    int pwl_epigraph_row_count{};
    int pwl_original_curve_count{};
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

// Starts from this run's verified base, which remains the immutable anchor
// for every corrective bound. No actual outage is solved or certified here.
SparseAcEconomicResult solve_sparse_common_corrective_reference(
    const CaseData& data, const std::vector<int>& commitment,
    const SolveResult& original_base,
    const SparseAcEconomicOptions& options = {});

void run_sparse_ac_corrective_reference_regression(
    const CaseData& data, const std::vector<int>& commitment,
    const AcState& original_base);

// Tiny-fixture checks of equivalent costs, unchanged physical rows/bounds,
// and the new analytic cost Jacobian. Does not solve a production case.
void run_sparse_ac_pwl_epigraph_regression(
    const CaseData& data, const std::vector<int>& commitment,
    const AcState& start);

// Finite differences of the complete Lagrangian gradient on tiny fixtures.
void run_sparse_ac_exact_hessian_regression(
    const CaseData& data, const std::vector<int>& commitment,
    const AcState& start);

}  // namespace gravityx
