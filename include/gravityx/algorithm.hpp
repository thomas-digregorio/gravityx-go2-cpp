#pragma once

#include "gravityx/ac_model.hpp"
#include "gravityx/case_data.hpp"
#include "gravityx/validation.hpp"

#include <nlohmann/json.hpp>

#include <limits>
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
    // Proposal consumers may need only the exact MILP commitment/dispatch.
    // Skip the nonlinear candidate repair in that case; no unverified point
    // is selected or serialized as the accepted answer.
    bool repair_candidate{true};
    // The component MILP supplies a component-balanced economic Pg/load
    // target.  Preserve that target during nonlinear repair instead of
    // replacing it with the generic feasibility redispatch.  The complete
    // nonlinear validator remains authoritative and the verified incumbent
    // is retained whenever this repair does not pass.
    bool preserve_candidate_dispatch_during_repair{};
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
    bool candidate_repair_preserved_dispatch{};
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

// A small mixed-integer commitment layer over the reduced, first-order AC
// response used by ActiveNetworkEconomicDispatch.  Only deterministic
// source-eligible toggles proposed by the exact component MILP receive binary
// columns.  All continuous controls retain source ramps and PMIN/PMAX, and
// monitored angle, voltage, thermal, aggregate-Q, and generator-outage
// headroom rows remain in the proposal model.  The mixed-integer point is
// never authoritative: nonlinear AC repair plus the complete validator must
// certify a strict objective improvement before it replaces the incumbent.
struct ReducedNetworkCommitmentOptions {
    double time_limit_seconds{30.0};
    double proposal_time_limit_seconds{5.0};
    double validation_tolerance{1e-5};
    double objective_tolerance{1e-9};
    double angle_trust_radius{0.02};
    double voltage_trust_radius{0.02};
    double generation_trust_radius{0.10};
    double load_factor_trust_radius{0.02};
    double thermal_row_utilization_threshold{0.70};
    double mip_relative_gap{1e-3};
    int maximum_rounds{8};
    int minimum_commitment_changes_per_round{1};
    int maximum_commitment_changes_per_round{1};
    int maximum_candidate_generators{96};
    int maximum_rejected_patterns_per_round{8};
    int maximum_constraint_generation_passes{6};
    int maximum_rows_per_constraint_pass{32};
    int maximum_new_shutdowns_per_bus{1};
    int maximum_startup_correction_candidates{16};
    int startup_correction_candidate_offset{};
    // When load movement is enabled, zero retains every source-eligible load
    // control.  A positive value keeps only the largest component-MILP
    // proposed MW movements, reducing the commitment model without changing
    // any selected load's source bounds, ramps, or PWL value curve.
    int maximum_movable_loads{};
    // A feasibility repair can leave an economically strong commitment with
    // an arbitrary dispatch.  Optionally spend a small, bounded amount of the
    // same global search budget on the existing fixed-commitment sequential
    // AC refinement before ranking the commitment.  Zero preserves the
    // original behavior.
    double post_repair_economic_seconds{};
    int post_repair_economic_maximum_rounds{2};
    // Before the first reduced commitment MILP, optionally move the voltage
    // targets at a small number of economically attractive multi-generator
    // PV buses so that retiring one unit does not begin from an avoidable
    // aggregate-Q-capability violation.  The proposed state is accepted only
    // after the ordinary nonlinear base power flow and complete validator
    // pass; it is a search reference, not a relaxed incumbent.
    bool prepare_reactive_headroom{};
    int reactive_headroom_maximum_buses{4};
    double reactive_headroom_maximum_voltage_step{0.005};
    double reactive_headroom_margin{0.02};
    // Large aggregate-Q deficits are structural search directions rather than
    // local voltage-conditioning opportunities.  Restrict preparation to
    // near-bound buses so a few hopelessly large adjustments cannot displace
    // the useful, inexpensive ones from the bounded preparation set.
    double reactive_headroom_maximum_residual{1.0};
    bool freeze_load_movement{true};
    bool enforce_generator_contingency_headroom{true};
    bool allow_nonimproving_continuation{};
    // Permit the exact component MILP to correct an earlier shutdown when
    // another generator at the same bus remains online.  Such a restart
    // preserves the reduced model's PV/PQ partition; source startup
    // eligibility, PMIN/PMAX, ramps, transition costs, Q capability,
    // nonlinear repair, and independent validation all remain enforced.
    bool allow_startup_corrections_at_pv_bus{};
    // Fill unused reduced-MILP candidate slots with deterministic,
    // source-eligible shutdowns and reversals of earlier internal shutdowns.
    // The exact component-MILP direction remains first in the pool.
    bool expand_candidate_pool_beyond_component_proposal{};
    bool require_one_startup_correction{};
    // In a correction beam, require one shutdown at a bus that already has a
    // newly retired source-online unit.  This repairs a sequential-search
    // artifact: several legal same-bus shutdowns can be selected together in
    // one batch, but a conservative per-bus guard can otherwise hide the last
    // member from all later rounds.  Exact aggregate-Q rows, nonlinear repair,
    // and validation remain authoritative.
    bool require_additional_same_bus_shutdown{};
    int same_bus_minimum_prior_shutdowns{1};
    bool require_mixed_reactive_shutdowns{};
    // Search only patterns containing at least one shutdown whose removal
    // changes the reference bus's aggregate reactive-capability envelope.
    // This is a beam-diversity constraint; exact AC repair and validation
    // still decide whether the resulting pattern is acceptable.
    bool require_reactive_coupled_shutdown{};
    bool enumerate_nonimproving_patterns{};
    // When exact no-good enumeration is enabled, continue past the first
    // verified improving pattern and retain the best repaired candidate seen
    // within the same bounded beam.  The monotone incumbent is updated only
    // after enumeration, never by a lower-quality alternative.
    bool enumerate_improving_patterns{};
    // After a repaired pattern is evaluated, exclude its selected
    // reactive-coupled shutdown set from the remainder of the local beam.
    // This is a heuristic search-diversity row, not a global feasibility cut;
    // the independently verified incumbent remains outside the beam.
    bool diversify_reactive_shutdowns_between_patterns{};
    // Reuse the same HiGHS model both while adding generated security rows
    // and, during exact enumeration, while adding globally valid no-good
    // rows.  Candidate bounds and the reference linearization must remain
    // unchanged for the latter reuse.
    bool reuse_solver_between_constraint_passes{};
    bool reuse_component_proposal_candidates{};
};

struct ReducedNetworkCommitmentResult {
    bool incumbent_verified{};
    bool proposal_feasible{};
    bool improved{};
    bool time_limit_reached{};
    double wall_seconds{};
    double proposal_wall_seconds{};
    double solver_wall_seconds{};
    double repair_wall_seconds{};
    double economic_polish_wall_seconds{};
    double incumbent_official_proxy{};
    double selected_official_proxy{};
    double selected_headroom_residual{};
    double frontier_official_proxy{};
    int proposal_change_count{};
    int proposal_refreshes{};
    int candidate_generator_count{};
    int rounds_completed{};
    int constraint_generation_passes{};
    int rejected_patterns{};
    int bus_shutdown_guard_rejections{};
    int accepted_moves{};
    int best_incumbent_updates{};
    int economic_polish_attempts{};
    int economic_polish_improvements{};
    int simplex_iterations{};
    int mip_nodes{};
    std::string status;
    std::vector<int> selected_commitment;
    std::vector<int> frontier_commitment;
    SolveResult selected;
    ValidationReport selected_validation;
    // The search frontier can be a verified, deliberately non-improving
    // bridge commitment whose dispatch has not yet received a full economic
    // polish.  Keep that state separate from `selected`: callers must never
    // mistake it for the monotone incumbent, but may use it as the seed for a
    // subsequent fixed-commitment optimization.
    SolveResult frontier;
    ValidationReport frontier_validation;
    nlohmann::json rounds = nlohmann::json::array();
    nlohmann::json proposals = nlohmann::json::array();
    nlohmann::json trials = nlohmann::json::array();
    nlohmann::json reactive_headroom_preparation = nullptr;

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
    // Optional per-round real-generation movement trust radius in per unit.
    // Infinity preserves the full source PMIN/PMAX and ramp interval. A
    // finite value intersects that interval only for the linear proposal;
    // nonlinear acceptance still uses the exact source bounds.
    double generation_trust_radius{
        std::numeric_limits<double>::infinity()};
    // Optional per-round dispatchable-load factor movement trust radius.
    // Infinity preserves the full source TMIN/TMAX and ramp interval. A
    // finite value intersects that interval only for the sequential linear
    // proposal; nonlinear acceptance still uses the exact source bounds.
    double load_factor_trust_radius{
        std::numeric_limits<double>::infinity()};
    // Add a continuous switched-shunt susceptance direction to the reduced
    // LP. Every nonlinear trial is projected back to exact source block
    // counts before it can be repaired, verified, or accepted.
    bool optimize_switched_shunts{};
    double switched_shunt_susceptance_trust_radius{
        std::numeric_limits<double>::infinity()};
    double switched_shunt_total_movement_budget{
        std::numeric_limits<double>::infinity()};
    double switched_shunt_movement_penalty{1e-3};
    // Optional diagonal proximal term on real-power control movement
    // (generation and dispatchable-load factor) in the sequential linearized
    // proposal objective. Zero preserves the LP. Positive values form a
    // convex QP that discourages a simultaneous bang-bang move across
    // thousands of controls; the exact GO2 objective remains the sole
    // nonlinear acceptance criterion.
    double control_proximal_weight{};
    double maximum_candidate_repair_balance_slack{0.5};
    double thermal_row_utilization_threshold{0.7};
    // Optional valid polygonal outer-approximation cuts on each side of the
    // current apparent-power tangent for branches whose source soft-margin
    // variable is near its upper bound. Zero preserves the single tangent.
    int thermal_polygon_side_cuts{};
    double thermal_polygon_angle_radians{0.02};
    double thermal_polygon_activation_margin{1e-3};
    int maximum_candidate_repair_newton_iterations{12};
    int maximum_voltage_rows_per_pass{8};
    int maximum_trust_region_rows_per_pass{32};
    // When a reduced PV/PQ direction first violates reactive capability,
    // also guard the nearest nonviolated buses up to this total monitored-row
    // target. Zero retains purely violated-row generation.
    int reactive_capability_guard_row_target{};
    int maximum_rounds{8};
    int maximum_candidate_trials{9};
    // Center the geometric nonlinear candidate schedule on the step accepted
    // in the preceding sequential-linearization round. This spends repair
    // work near the measured local trust radius instead of restarting nine
    // times from a full step on every round.
    bool adaptive_candidate_fractions{};
    // Re-linearize immediately after the first independently verified
    // improving fraction. This avoids repairing smaller points on a stale
    // direction when a deliberately local trust-region step already passes.
    bool accept_first_improving_fraction{};
    // Diagnostic/proposal profile for a common corrective target: use the
    // source contingency ramps, interval weight, and RATE_C while retaining
    // the no-outage topology. The returned state is never an accepted
    // contingency answer by itself; each outage still receives full source
    // reconstruction and validation.
    bool contingency_profile{};
    bool simplex_presolve{};
    bool use_simplex_crash_basis{true};
    // Reuse the preceding optimal simplex basis when the next sequential
    // linearization has the same row and column dimensions. The matrix and
    // bounds are still rebuilt exactly; HiGHS reoptimizes the inherited
    // basis against the new coefficients.
    bool reuse_simplex_basis_between_rounds{};
    bool compact_signed_columns{};
    bool freeze_load_movement{};
    bool eliminate_angles{};
    bool reduced_pv_pq_partition{};
    // Hold the initial PV/PQ partition fixed during one sequential
    // refinement. Aggregate reactive-capability rows remain enforced and the
    // independently verified nonlinear state remains the acceptance gate.
    bool stable_reduced_pv_pq_partition{};
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
    int reused_basis_attempts{};
    int reused_basis_accepted{};
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

ReducedNetworkCommitmentResult
refine_reduced_network_economic_commitment(
    const CaseData& data,
    const std::vector<int>& incumbent_commitment,
    const SolveResult& incumbent,
    const ReducedNetworkCommitmentOptions& options = {});

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
