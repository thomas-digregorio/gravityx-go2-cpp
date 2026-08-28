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

// Network-free one-period unit-commitment proposal over each connected
// component.  This is deliberately not an accepted AC or security solution:
// exact source eligibility, conditional PMIN/PMAX, ramps, transition costs,
// PWL economics, component balance, and generator-outage headroom are enforced
// in the MILP, after which the candidate must still pass nonlinear base repair
// and independent validation.  The verified incumbent is retained otherwise.
struct ComponentCommitmentOptions {
    double time_limit_seconds{10.0};
    double mip_relative_gap{1e-3};
    double validation_tolerance{1e-5};
    double objective_tolerance{1e-9};
    bool enforce_generator_contingency_headroom{true};
    int maximum_commitment_changes{};
    // Negative disables the row family.  A nonnegative value bounds the
    // absolute change in generation minus responsive load at every bus.
    double bus_active_injection_trust_radius{-1.0};
    bool initialize_near_incumbent_dispatch{true};
    // Co-optimize commitment and dispatch against a first-order AC network
    // model around the independently verified incumbent.  Voltage magnitudes
    // remain fixed, while voltage angles and generator reactive power move.
    // This is a proposal model only; nonlinear repair and validation remain
    // authoritative before a candidate can replace the incumbent.
    bool linearized_active_network{};
    bool linearized_reactive_network{};
    double angle_trust_radius{0.02};
    double voltage_trust_radius{0.02};
    double thermal_row_utilization_threshold{0.70};
    // Zero omits the L1 movement auxiliaries from the large MILP.  A positive
    // value is useful for small diagnostics or a fixed-commitment polishing
    // LP, but roughly doubles the network-state row count.
    double network_movement_penalty{0.0};
    // Zero leaves every source-eligible commitment binary free.  A positive
    // value ranks economically plausible changes and fixes all remaining
    // commitments to the verified incumbent before model construction.
    int maximum_candidate_generators{};
    // Replace per-point convex-combination lambdas by mathematically
    // equivalent power plus PWL epigraph/hypograph variables when every
    // active source curve has the required convex/concave slope ordering.
    bool compact_pwl_formulation{};
};

struct ComponentCommitmentResult {
    bool incumbent_verified{};
    bool attempted{};
    bool solver_feasible{};
    bool solver_optimal{};
    bool mip_start_attempted{};
    bool mip_start_accepted{};
    bool compact_pwl_formulation{};
    bool candidate_repair_attempted{};
    bool candidate_repair_feasible{};
    bool candidate_repair_converged{};
    bool candidate_verified{};
    bool improved{};
    bool time_limit_reached{};
    double wall_seconds{};
    double solver_wall_seconds{};
    double candidate_repair_wall_seconds{};
    double incumbent_objective{};
    double raw_candidate_objective{};
    double candidate_objective{};
    double incumbent_penalty_slack{};
    double raw_candidate_penalty_slack{};
    double candidate_penalty_slack{};
    double incumbent_transition_cost{};
    double candidate_transition_cost{};
    double incumbent_official_proxy{};
    double candidate_official_proxy{};
    double solver_objective{};
    double maximum_milp_residual{};
    double mip_start_maximum_column_violation{};
    double mip_start_maximum_row_violation{};
    double candidate_headroom_residual{};
    std::string maximum_milp_residual_identity;
    bool used_near_incumbent_dispatch{};
    int component_count{};
    int bus_active_injection_trust_rows{};
    int linearized_angle_columns{};
    int linearized_voltage_columns{};
    int linearized_reactive_generation_columns{};
    int linearized_active_balance_rows{};
    int linearized_reactive_balance_rows{};
    int linearized_reactive_capability_rows{};
    int linearized_angle_limit_rows{};
    int linearized_thermal_rows{};
    int generator_contingency_headroom_rows{};
    int candidate_generator_count{};
    int fixed_incumbent_generator_count{};
    int online_before{};
    int online_after{};
    int startup_count{};
    int shutdown_count{};
    int commitment_change_count{};
    int row_count{};
    int column_count{};
    int nonzero_count{};
    int run_status{};
    int model_status{};
    int primal_solution_status{};
    int mip_start_status{};
    int mip_start_worst_column{-1};
    int mip_start_worst_row{-1};
    int mip_node_count{};
    double mip_gap{};
    double mip_dual_bound{};
    std::string status;
    std::string candidate_repair_failure_reason;
    std::vector<int> candidate_commitment;
    std::vector<int> selected_commitment;
    SolveResult candidate;
    ValidationReport raw_candidate_validation;
    ValidationReport candidate_validation;
    SolveResult selected;
    ValidationReport selected_validation;

    nlohmann::json to_json(bool include_state = false) const;
};

struct GreedyCommitmentSearchOptions {
    double time_limit_seconds{30.0};
    double proposal_time_limit_seconds{5.0};
    double validation_tolerance{1e-5};
    double objective_tolerance{1e-9};
    int maximum_rounds{4};
    int maximum_candidates_per_round{64};
    bool enforce_generator_contingency_headroom{true};
    // The component MILP supplies a globally coherent commitment direction,
    // but each nonlinear toggle remains authoritative.  Stop a round at the
    // first verified improving toggle instead of repairing every remaining
    // alternative merely to choose the largest one-step gain.
    bool accept_first_improving_toggle{true};
};

struct GreedyCommitmentSearchResult {
    bool incumbent_verified{};
    bool proposal_attempted{};
    bool proposal_feasible{};
    bool improved{};
    bool time_limit_reached{};
    double wall_seconds{};
    double proposal_wall_seconds{};
    double candidate_repair_wall_seconds{};
    double incumbent_objective{};
    double selected_objective{};
    double incumbent_official_proxy{};
    double selected_official_proxy{};
    double incumbent_penalty_slack{};
    double selected_penalty_slack{};
    double selected_headroom_residual{};
    int proposal_change_count{};
    int candidate_pool_size{};
    int rounds_completed{};
    int candidates_attempted{};
    int candidates_precheck_rejected{};
    int candidates_repaired{};
    int candidates_verified{};
    int accepted_moves{};
    int candidate_order_refreshes{};
    int first_improvement_selections{};
    std::string status;
    std::vector<int> selected_commitment;
    SolveResult selected;
    ValidationReport selected_validation;
    nlohmann::json trials = nlohmann::json::array();

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
    double voltage_trust_radius{0.02};
    double maximum_candidate_repair_balance_slack{0.5};
    double thermal_row_utilization_threshold{0.7};
    int maximum_candidate_repair_newton_iterations{12};
    int maximum_voltage_rows_per_pass{8};
    int maximum_trust_region_rows_per_pass{32};
    int maximum_rounds{8};
    int maximum_candidate_trials{9};
    // Diagnostic/proposal profile for a common corrective target: use the
    // source contingency ramps, interval weight, and RATE_C while retaining
    // the no-outage topology. The returned state is never an accepted
    // contingency answer by itself; each outage still receives full source
    // reconstruction and validation.
    bool contingency_profile{};
    bool simplex_presolve{};
    bool use_simplex_crash_basis{true};
    bool compact_signed_columns{};
    bool freeze_load_movement{};
    bool eliminate_angles{};
    bool reduced_pv_pq_partition{};
    bool generate_trust_region_rows{};
    // Retain the complete sparse linearized P/Q balance system with explicit
    // angle, voltage, and reactive-generation movements. This puts the state
    // trust region inside the LP instead of shrinking an unconstrained
    // reduced-space direction after the solve.
    bool sparse_full_ac_linearization{};
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
    int angle_response_row_count{};
    int voltage_response_row_count{};
    int reactive_capability_row_count{};
    int trust_region_constraint_generation_passes{};
    int voltage_constraint_generation_passes{};
    int reactive_capability_constraint_generation_passes{};
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

ComponentCommitmentResult refine_component_economic_commitment(
    const CaseData& data,
    const std::vector<int>& incumbent_commitment,
    const SolveResult& incumbent,
    const ComponentCommitmentOptions& options = {});

GreedyCommitmentSearchResult refine_greedy_economic_commitment(
    const CaseData& data,
    const std::vector<int>& incumbent_commitment,
    const SolveResult& incumbent,
    const GreedyCommitmentSearchOptions& options = {});

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

void run_active_network_reduction_regression();

IbrResult run_iterative_batch_rounding(const CaseData& data, const IbrOptions& options = {});

}  // namespace gravityx
