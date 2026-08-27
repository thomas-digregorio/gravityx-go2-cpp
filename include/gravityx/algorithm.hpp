#pragma once

#include "gravityx/ac_model.hpp"
#include "gravityx/case_data.hpp"
#include "gravityx/validation.hpp"

#include <nlohmann/json.hpp>

#include <vector>

namespace gravityx {

struct IbrOptions {
    int batch_count{4};
    bool source_status_only{};
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

struct SparseEconomicRefinementOptions {
    double time_limit_seconds{60.0};
    double validation_tolerance{1e-5};
    double objective_tolerance{1e-9};
    int maximum_rounds{24};
    int maximum_linear_economic_rounds{1};
    double linear_economic_time_limit_seconds{60.0};
    double linear_economic_voltage_trust_radius{0.02};
    double linear_economic_angle_trust_radius{0.05};
    int voltage_coordinate_bus_count{64};
    int maximum_voltage_coordinate_passes{8};
};

struct SparseEconomicRefinementResult {
    bool incumbent_verified{};
    bool attempted{};
    bool improved{};
    bool time_limit_reached{};
    double wall_seconds{};
    double incumbent_objective{};
    double selected_objective{};
    SolveResult selected;
    ValidationReport selected_validation;
    nlohmann::json rounds = nlohmann::json::array();

    nlohmann::json to_json(bool include_state = false) const;
};

bool validated_candidate_is_feasible(
    const SolveResult& result,
    const ValidationReport& validation,
    double tolerance);

bool verified_economic_candidate_improves_incumbent(
    const SolveResult& incumbent,
    const ValidationReport& incumbent_validation,
    const SolveResult& candidate,
    const ValidationReport& candidate_validation,
    double validation_tolerance,
    double objective_tolerance = 1e-9);

SparseEconomicRefinementResult refine_fixed_commitment_sparse(
    const CaseData& data,
    const std::vector<int>& commitment,
    const SolveResult& incumbent,
    const SparseEconomicRefinementOptions& options = {});

IbrResult run_iterative_batch_rounding(const CaseData& data, const IbrOptions& options = {});

}  // namespace gravityx
