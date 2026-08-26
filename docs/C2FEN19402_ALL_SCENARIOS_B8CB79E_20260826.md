# C2FEN19402 all-scenario recovery — `b8cb79e` — 2026-08-26

## Verdict

All five GO Challenge 2 Final Event Division 1 scenarios for the parsed
19,402-bus family passed from cold starts on the frozen implementation. Every
run completed its complete source-supplied contingency set, returned official
infeasibility `0.0`, passed the independent nonlinear residual checks, and
finished inside the 300-second end-to-end limit.

No solver source, acceptance threshold, contingency set, or executable was
changed between these five runs.

## Frozen boundary

- Git revision: `b8cb79ec82e3e5d2ff7eb23745cd0f424fd5e6de`
- Git tag: `c2fen19402-recovery-freeze-20260826-v2`
- Standard exact executable SHA-256:
  `e870443dd975d1373246d1bf9beacc5d1094711bcf19f458a0157f4837f4be05`
- Native executable SHA-256:
  `5689c9fcbfe2461424b07fe73d1a045ecbcbc67e2cb151c5b7cf5588f9158495`
- Scenario 010 timing-only profile SHA-256:
  `d739ae2d353e762e61c2e9de6f2660548f056821abd5e26799f658d984d260aa`
- Pre-freeze tests: standard CTest `4/4`, native CTest `4/4`, Python `46/46`.
- Run artifacts:
  `runs/c2fen19402_recovery_b8cb79e_20260826`

Scenario 006 used the standard executable for its numerically sensitive exact
base solve and the native executable for parallel contingency screening. The
other scenarios used the native executable for both stages. Scenario 010 used
the case-hash-bound `full_v6` schedule containing timing measurements for all
6,693 labels; it contains no prior primal, dual, commitment, network, or solver
state.

## Results

Times below are the serialized end-to-end boundaries in `run_summary.json`.
Streaming official evaluation overlaps contingency execution, so base and
contingency times should not be added to reconstruct the total.

| Scenario | Completed/source contingencies | Base (s) | Contingency (s) | End-to-end (s) | Official objective | Official infeasibility | Maximum independent contingency residual | Exact fallbacks | Result |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 006 | 6,693/6,693 | 27.792 | 220.796 | 256.166 | -679,351,489.807097 | 0.0 | 9.996e-6 | 0 | PASS |
| 010 | 6,693/6,693 | 24.832 | 255.285 | 286.062 | -582,029,659.645846 | 0.0 | 9.782e-6 | 0 | PASS |
| 069 | 6,620/6,620 | 26.481 | 213.204 | 246.126 | -594,813,595.851483 | 0.0 | 9.921e-6 | 0 | PASS |
| 077 | 6,584/6,584 | 23.053 | 136.754 | 171.659 | -433,016,087.214137 | 0.0 | 8.113e-6 | 1 | PASS |
| 095 | 6,579/6,579 | 6.506 | 121.374 | 139.037 | -36,294,035.586258 | 0.0 | 9.913e-6 | 0 | PASS |

Across the suite, 33,169 source contingencies were solved and independently
checked. The largest reported contingency residual was
`9.99584980010404e-6`, below the unchanged `1e-5` acceptance threshold. Each
official evaluation certificate observed exactly one base-case detail plus
every expected contingency detail, listed no infeasible labels, and reproduced
the reported objective exactly from the per-case details.

## Recovered regressions

The prior frozen regression failed scenario 006 in the native base path and
timed out scenario 010 near the tail. Scenario 006 was recovered by retaining
the standard exact build for its base solve. Scenario 010's runtime regression
was traced to repeatedly rebuilding connected components during source
reference-angle normalization inside contingency predictor/Newton iterations.

Revision `b8cb79e` computes the base connected components and bridge flags once
per resident solver. Non-bridge outages reuse that immutable topology; bridge
outages still construct their exact post-outage components. Angle
normalization, source constraints, nonlinear validation, and official scoring
remain unchanged. A parallel-edge/bridge/island normalization regression test
was added with the cache.
