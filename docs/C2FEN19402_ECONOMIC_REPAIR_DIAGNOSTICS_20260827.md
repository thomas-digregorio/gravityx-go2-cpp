# C2FEN19402 economic-repair diagnostics (2026-08-27)

This checkpoint is a targeted diagnostic milestone on the separate
`codex/c2fen19402-economic-improvement` branch. It does not replace the frozen
`c477a2d4673a87d4d569110471956e420bb919e6` six-family evidence, and it is not a
full scenario result.

## Scope and unchanged gates

- Source case: C2FEN19402 scenario 095.
- Commitment and base point: previously verified scenario-095 sparse-AC point.
- Source contingencies, source bounds, Division 1 controls, objective semantics,
  and the independent nonlinear validator are unchanged.
- A candidate can replace the incumbent only when its rebuilt GO2 objective is
  strictly higher and its complete nonlinear residual is no greater than
  `1e-5`.
- All diagnostics below are non-official, fixed-subset calibration runs.

## Base economic direction experiments

The component economic relaxation had a market-surplus ceiling of approximately
`217169.5`, but routing that direction through the complete AC equations improved
the verified base objective only from `211834.469236` to `211922.451801`
(`+87.982565`). The selected point passed at approximately machine precision and
required `13.25 s`. This route is retained as a diagnostic tool but is not yet
part of the production full-scenario path because its measured gain is too small.

## Corrective-tail root cause

The worst retained sampled contingency, `CTG_001661`, was not suffering from a
stale serialized objective. Its exact Newton state exceeded source voltage bounds
at 22 buses, with a maximum projection of `0.280299 p.u.`. The earlier linearized
security repair rejected that reference before solving.

The verified correction is:

1. Project the repair reference into the unchanged source `VMIN/VMAX` bounds.
2. Rebuild all contingency derived fields.
3. Run the original small security-repair LP first.
4. Only when voltage projection was required, run a short Phase-I LP that
   minimizes paid active/reactive balance slack.
5. Line-search both proposals against the current verified incumbent and accept
   only a strict, independently validated objective improvement.

The order matters. Replacing the original repair with Phase I caused large
generator-outage regressions because the larger LP reached its time limit before
recovering the high-value security point. The ordered cascade restores that path.

The active-feasibility repair now also audits a time-limited terminal vector
against every original row and column bound before preserving it. This prevents
a feasible HiGHS terminal point from being falsely rejected by a subsequent
canonicalization step; an unaudited terminal vector is never exposed to the
nonlinear validator.

## Five-label monotonicity check

The final bounded check used one resident worker and evaluated each selected
label once. The comparison column uses the matching label from the earlier
100-case exact-security-repair sample.

| Label | Earlier objective | Ordered-cascade objective | Change | Residual | Solve time (s) |
|---|---:|---:|---:|---:|---:|
| CTG_001661 | -132434.090 | -1933.987 | +130500.103 | 4.69e-12 | 4.743 |
| CTG_005116 | -29400.085 | 29742.942 | +59143.027 | 3.67e-6 | 3.265 |
| CTG_000000 | 91383.083 | 91374.351 | -8.732 | 8.88e-16 | 4.799 |
| CTG_000133 | 136833.489 | 136833.489 | 0.000 | 1.78e-15 | 1.768 |
| CTG_006312 | 161802.239 | 183837.195 | +22034.956 | 6.50e-6 | 4.492 |

All five contingencies passed. The aggregate fixed-subset wall time was
`21.597 s`, and the maximum residual was `6.5044e-6`. The `-8.732` change for
`CTG_000000` is a 0.0096% time-limited LP-path variation, not a security loss or
a candidate-selection regression within a run.

## Rejected paths

The following targeted experiments were not retained in the production path:

- Extending the same Phase-I solve from one to two seconds did not improve
  `CTG_001661`.
- A recentered Phase-I attempt and a nested forced-zero-balance Newton solve
  added runtime without material objective improvement.
- The diagnostic HiGHS IPM path stopped after two iterations with a large primal
  infeasibility and was worse than simplex for this repair LP.
- The component-only base economic direction produced too little verified gain
  for its runtime.

## Sequential nonlinear security repair

A later 100-label longest-first calibration reached its 180-second diagnostic
limit after 37 labels. Thirty-four completed labels were generator outages, and
their initial corrective states carried multi-million-dollar paid imbalance.
This configuration was explicitly rejected for a full run: it averaged
`4.692 s` per completed label on one resident worker and could not satisfy the
300-second system boundary without tighter routing.

Attribution on severe generator outage `CTG_000353` showed that the first full
exact-Newton security-repair candidate improved the objective from
`-2293035.173` to `29194.751` with exact P/Q balance, but one branch required
overload slack `0.144758` above its source bound. The one-shot linearization
therefore accepted only a `1/64` blend and retained most of the economic loss.

The retained correction performs at most four compact sequential security
repairs, relinearizing only while:

- the full candidate has a strictly higher objective than the verified
  incumbent;
- the only blocker is a branch flow, angle, or branch-slack bound;
- the complete nonlinear residual strictly decreases; and
- the global economic-polish deadline remains.

On the same `CTG_000353` point, the four sequential rounds reduced the nonlinear
branch residual from `0.144758` to `2.50e-10`, raised the verified objective from
`-2256039.412` to `-15107.286`, and remained within the unchanged source model.
This is a gain of approximately `2.241 million` on that contingency. The solve
took `6.072 s`; further work is required to remove the remaining paid P/Q slack
and to gate the path to the generator-outage tail before a full scenario run.

## Status

The component test suite passes 4/4. No official scenario-095 rerun has been
performed at this checkpoint. A broader representative calibration is still
required before spending the one permitted cold full-scenario run.

## Economic-target branch Newton and complete state transfer

The fixed-Jacobian feasibility predictor can alter generation, load, and
switched-shunt controls substantially before it reaches a secure state.  The
refreshed-Jacobian Newton rescue previously solved the AC equations for those
already-distorted controls.  The retained branch-outage correction now:

1. preserves the economically preferred post-outage controls before the local
   feasibility projections;
2. uses the secure predictor voltages only as the Newton initialization;
3. transfers load factors and the exact discrete switched-shunt state into the
   Newton and repair candidates; and
4. permits the existing monotone sequential security repair for branch as well
   as generator outages.

Generator outages retain their previously verified feasibility-polished Newton
reference.  An initial diagnostic showed that replacing it could miss an exact
corrective generator bound by approximately `0.001 p.u.`; that path was rejected
before retention.

The final six-label diagnostic evaluated each label once.  Every point remained
secure under the unchanged complete nonlinear validator.

| Label | Type | Prior objective | Retained objective | Change | Residual | Time (s) |
|---|---|---:|---:|---:|---:|---:|
| CTG_000353 | generator | -15107.286 | -15107.286 | 0.000 | 2.50e-10 | 5.861 |
| CTG_000111 | generator | 212203.094 | 212203.094 | 0.000 | 9.25e-16 | 2.133 |
| CTG_000545 | branch | -1263951.215 | -14624.067 | +1249327.148 | 1.31e-10 | 6.232 |
| CTG_001538 | branch | -94528.054 | 199370.469 | +293898.523 | 3.02e-8 | 3.503 |
| CTG_002190 | branch | 211729.810 | 211680.760 | -49.050 | 1.19e-15 | 2.581 |
| CTG_002330 | branch | -5552.737 | 193947.903 | +199500.640 | 4.71e-7 | 4.799 |

The tiny `CTG_002190` variation is a `0.023%` bounded-path difference.  Within
every run, candidate replacement remains strictly objective-monotone.  A new
component regression verifies that all continuous controls are interpolated
and that the target switched-shunt state is transferred atomically; 4/4 tests
pass after this correction.  These results remain targeted diagnostics, not an
official scenario score.
