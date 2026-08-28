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

struct BusInjectionCommitmentOptions {
    double time_limit_seconds{10.0};
    double mip_relative_gap{1e-3};
    double validation_tolerance{1e-5};
    double objective_tolerance{1e-9};
    bool enforce_generator_contingency_headroom{true};
};

struct BusInjectionCommitmentResult {
    bool attempted{};
    bool solver_feasible{};
    bool candidate_verified{};
    bool improved{};
    double wall_seconds{};
    double incumbent_objective{};
    double candidate_objective{};
    double incumbent_transition_cost{};
    double candidate_transition_cost{};
    double incumbent_official_proxy{};
    double candidate_official_proxy{};
    double maximum_bus_active_injection_change{};
    double maximum_bus_reactive_injection_change{};
    int online_before{};
    int online_after{};
    int shutdown_count{};
    int row_count{};
    int column_count{};
    int nonzero_count{};
    int run_status{};
    int model_status{};
    int primal_solution_status{};
    int mip_node_count{};
    double mip_gap{};
    double mip_dual_bound{};
    std::string status;
    std::vector<int> commitment;
    SolveResult selected;
    ValidationReport selected_validation;

    nlohmann::json to_json(bool include_state = false) const;
};

// A deliberately small fixed-commitment relaxation used only to propose an
// economic direction.  It preserves the source PWL curves, exact base ramp
// and PMIN/PMAX/load bounds, and active-power balance within every connected
// component.  Network equations and limits are not relaxed in the accepted
// answer: every interpolated proposal is passed through the nonlinear fast
// power flow and the complete validator, and the verified incumbent is kept
// unless a strictly better feasible point is found.
struct ComponentEconomicDispatchOptions {
    double time_limit_seconds{5.0};
    double validation_tolerance{1e-5};
    double objective_tolerance{1e-9};
    int maximum_rounds{32};
    int maximum_candidate_trials{9};
};

struct ComponentEconomicDispatchResult {
    bool incumbent_verified{};
    bool attempted{};
    bool solver_feasible{};
    bool solver_optimal{};
    bool improved{};
    bool time_limit_reached{};
    bool primal_start_attempted{};
    bool primal_start_accepted{};
    double wall_seconds{};
    double solver_wall_seconds{};
    double incumbent_objective{};
    double relaxed_market_surplus{};
    double selected_objective{};
    double selected_fraction{};
    int component_count{};
    int rounds_completed{};
    int row_count{};
    int column_count{};
    int nonzero_count{};
    int run_status{};
    int model_status{};
    int primal_solution_status{};
    int primal_start_status{};
    int simplex_iterations{};
    int ipm_iterations{};
    std::string status;
    SolveResult selected;
    ValidationReport selected_validation;
    nlohmann::json trials = nlohmann::json::array();

    nlohmann::json to_json(bool include_state = false) const;
};

struct ActiveNetworkEconomicDispatchOptions {
    double time_limit_seconds{10.0};
    double validation_tolerance{1e-5};
    double objective_tolerance{1e-9};
    double angle_trust_radius{0.02};
    double thermal_row_utilization_threshold{0.7};
    int maximum_rounds{8};
    int maximum_candidate_trials{9};
};

struct ActiveNetworkEconomicDispatchResult {
    bool incumbent_verified{};
    bool attempted{};
    bool solver_feasible{};
    bool all_solver_rounds_optimal{true};
    bool improved{};
    bool time_limit_reached{};
    double wall_seconds{};
    double solver_wall_seconds{};
    double incumbent_objective{};
    double selected_objective{};
    double selected_fraction{};
    double maximum_selected_angle_change{};
    int component_count{};
    int rounds_completed{};
    int thermal_row_count{};
    int row_count{};
    int column_count{};
    int nonzero_count{};
    int simplex_iterations{};
    int ipm_iterations{};
    std::string status;
    SolveResult selected;
    ValidationReport selected_validation;
    nlohmann::json rounds = nlohmann::json::array();
    nlohmann::json trials = nlohmann::json::array();

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

double base_commitment_transition_cost(
    const CaseData& data,
    const std::vector<int>& commitment);

// Optimize only the split of P/Q and on/off status among generators connected
// to the same bus. Exact per-bus P and Q totals are fixed to a verified
// incumbent, so the AC network point is unchanged. The returned candidate is
// nevertheless rebuilt and independently validated before it may improve the
// incumbent.
BusInjectionCommitmentResult refine_commitment_preserving_bus_injections(
    const CaseData& data,
    const std::vector<int>& incumbent_commitment,
    const SolveResult& incumbent,
    const BusInjectionCommitmentOptions& options = {});

ComponentEconomicDispatchResult refine_fixed_commitment_component_economic(
    const CaseData& data,
    const std::vector<int>& commitment,
    const SolveResult& incumbent,
    const ComponentEconomicDispatchOptions& options = {});

ActiveNetworkEconomicDispatchResult
refine_fixed_commitment_active_network_economic(
    const CaseData& data,
    const std::vector<int>& commitment,
    const SolveResult& incumbent,
    const ActiveNetworkEconomicDispatchOptions& options = {});

IbrResult run_iterative_batch_rounding(const CaseData& data, const IbrOptions& options = {});

}  // namespace gravityx
