# Frozen six-family Final Event regression — 2026-08-26

## Verdict

The frozen implementation passed 34 of 37 cold GO Challenge 2 Final Event 4
Division 1 scenario runs. Every reported PASS completed the complete source
contingency set, passed the independent nonlinear checks, returned official
infeasibility `0.0`, and stayed below the 300-second end-to-end boundary.

Three scenarios did not pass:

- `C2FEN02020` scenario `260` failed after 298 of 300 contingencies when the
  exact corrective worker could not construct the linearized AC HiGHS model
  for the next fallback task.
- `C2FEN19402` scenario `006` failed its base-point acceptance check with a
  `0.001561013488` p.u. variable-bound residual at `branch:6104:pt`.
- `C2FEN19402` scenario `010` solved all 6,693 contingencies but exhausted the
  end-to-end work budget before final official evaluation. It is a timeout,
  not a verified solution.

No solver source, tolerance, acceptance rule, or algorithm setting was changed
in response to those outcomes, and none of the three scientific failures was
retried.

## Frozen boundary

- Git tag: `regression-freeze-20260826`
- Git revision: `c7fa432a8922839e4efe7bbea1aa4dc1278bf216`
- Branch: `codex/gravity-framework-v2`
- Native executable SHA-256:
  `a3fe124aba539ec8de04654c67145e9b25da11ab80618cca09a7670a35cf91a2`
- Runner SHA-256:
  `09d0e8e567ad8484a372c6fb34face709c9d51d13477d170831fc62cb997751a`
- Referenced-equation evaluator SHA-256:
  `3a3e465b6cb023866b1288d9088fb74c25cfda9befa05c6ec864020fc40288c8`
- Referenced vendor evaluator SHA-256:
  `3e98ca0e5dda571fc90c2b16bc4000d652f50c02acc0c9a203d99d132aeb603e`
- Pre-run tests: standard CTest `4/4`, native CTest `4/4`, Python `47/47`.
- All output and source paths were outside OneDrive.

The common scalable path used the native executable, one exact corrective
worker, 24 fast-screen workers, complete nonlinear acceptance checks, and
streaming referenced-equation evaluation. The end-to-end limit was 300
seconds, with two seconds reserved for evaluation and two seconds for final
serialization. Validation-tail sizes were fixed by contingency-count class so
that every shard partition was valid. The 8,300-bus family used the previously
registered robust base initialization. `C2FEN19402` scenario `010` used the
versioned timing-only heavy-lane profile
`config/C2FEN19402_s010_fast_screen_heavy_v4.json`; it contains no primal,
dual, commitment, network, or solver state. No operating point was carried
between scenarios.

## Family summary

| Parsed buses | Case family | Scenarios | PASS | Failed/timed out | PASS time range (s) |
|---:|---|---:|---:|---:|---:|
| 617 | `C2FEN00617` | 5 | 5 | 0 | 11.123–14.331 |
| 2,020 | `C2FEN02020` | 5 | 4 | 1 | 11.785–167.433 |
| 4,224 | `C2FEN04200` | 7 | 7 | 0 | 123.760–125.644 |
| 8,300 | `C2FEN08300` | 7 | 7 | 0 | 37.134–238.075 |
| 16,789 | `C2FEN17700` | 8 | 8 | 0 | 43.955–51.303 |
| 19,402 | `C2FEN19402` | 5 | 3 | 2 | 147.967–271.598 |
| **Total** |  | **37** | **34** | **3** |  |

## Scenario results

Times are the serialized end-to-end boundaries in `run_status.json`. The
contingency phase overlaps streaming validation, so base and contingency times
should not be added to infer the total. A dash means that no certified final
value exists.

| Buses | Scenario | Completed/source contingencies | Base (s) | Contingency (s) | End-to-end (s) | Official objective | Max independent residual | Result |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 617 | 005 | 105/105 | 2.418 | 9.889 | 14.331 | 837,528.141490 | 2.220e-16 | PASS |
| 617 | 017 | 103/103 | 0.090 | 9.261 | 11.216 | 677,959.663942 | 1.110e-16 | PASS |
| 617 | 024 | 104/104 | 0.114 | 9.276 | 11.123 | 1,024,779.091343 | 4.441e-16 | PASS |
| 617 | 062 | 103/103 | 0.108 | 9.266 | 11.242 | 973,222.118662 | 2.220e-16 | PASS |
| 617 | 073 | 107/107 | 0.100 | 9.338 | 11.216 | 758,559.216716 | 4.441e-16 | PASS |
| 2,020 | 025 | 292/292 | 0.171 | 9.615 | 11.785 | 3,984,584.198397 | 8.882e-16 | PASS |
| 2,020 | 121 | 300/300 | 0.276 | 165.041 | 167.433 | -1,745,787.467477 | 1.872e-06 | PASS |
| 2,020 | 134 | 300/300 | 0.180 | 9.803 | 12.406 | 4,735,923.219387 | 1.776e-15 | PASS |
| 2,020 | 260 | 298/300 | 0.166 | 12.965 | 13.420 | — | — | FAIL |
| 2,020 | 262 | 300/300 | 0.260 | 44.336 | 46.362 | -3,013,880.273921 | 1.421e-14 | PASS |
| 4,224 | 009 | 455/455 | 15.268 | 106.322 | 124.087 | 1,079,986.571967 | 2.842e-14 | PASS |
| 4,224 | 010 | 455/455 | 15.156 | 106.247 | 123.760 | 1,335,684.862638 | 2.842e-14 | PASS |
| 4,224 | 014 | 455/455 | 15.373 | 107.987 | 125.644 | 1,058,163.791479 | 7.105e-15 | PASS |
| 4,224 | 055 | 455/455 | 15.339 | 107.379 | 124.853 | 2,093,982.002567 | 1.421e-14 | PASS |
| 4,224 | 056 | 455/455 | 15.336 | 107.404 | 125.024 | 2,380,103.181476 | 1.421e-14 | PASS |
| 4,224 | 057 | 455/455 | 15.586 | 107.364 | 125.275 | 2,448,049.880033 | 2.842e-14 | PASS |
| 4,224 | 060 | 455/455 | 15.445 | 107.664 | 125.415 | 2,017,696.701851 | 1.421e-14 | PASS |
| 8,300 | 003 | 607/607 | 2.936 | 33.603 | 40.108 | -13,609,513.240397 | 6.256e-06 | PASS |
| 8,300 | 012 | 607/607 | 0.720 | 32.385 | 39.300 | -24,501,716.970450 | 6.256e-06 | PASS |
| 8,300 | 013 | 607/607 | 0.620 | 32.899 | 37.134 | -24,332,327.798295 | 6.256e-06 | PASS |
| 8,300 | 022 | 607/607 | 0.638 | 33.434 | 37.431 | -22,032,420.321919 | 6.256e-06 | PASS |
| 8,300 | 043 | 607/607 | 0.673 | 31.601 | 39.628 | -26,809,846.003040 | 5.684e-14 | PASS |
| 8,300 | 052 | 607/607 | 0.744 | 35.859 | 41.830 | -24,902,976.143012 | 5.684e-14 | PASS |
| 8,300 | 166 | 603/603 | 0.832 | 233.843 | 238.075 | -22,389,627.637556 | 1.705e-13 | PASS |
| 16,789 | 019 | 238/238 | 1.236 | 41.624 | 46.927 | -77,574,836.444221 | 2.720e-06 | PASS |
| 16,789 | 020 | 238/238 | 1.214 | 41.937 | 47.502 | -77,474,813.228184 | 2.720e-06 | PASS |
| 16,789 | 021 | 238/238 | 1.275 | 41.366 | 47.486 | -77,591,853.767029 | 2.720e-06 | PASS |
| 16,789 | 089 | 238/238 | 1.273 | 44.982 | 51.303 | -77,682,806.676277 | 3.301e-06 | PASS |
| 16,789 | 094 | 238/238 | 1.399 | 33.098 | 43.955 | -77,628,200.543698 | 8.149e-06 | PASS |
| 16,789 | 106 | 238/238 | 1.451 | 39.931 | 47.973 | -77,398,785.958770 | 1.137e-13 | PASS |
| 16,789 | 107 | 238/238 | 1.274 | 37.336 | 45.264 | -77,604,085.796992 | 1.137e-13 | PASS |
| 16,789 | 115 | 238/238 | 1.538 | 41.741 | 48.028 | -80,045,119.805866 | 2.844e-06 | PASS |
| 19,402 | 006 | 0/6,693 | 28.394 | — | 28.784 | — | 1.561e-03 base | FAIL |
| 19,402 | 010 | 6,693/6,693 | 24.116 | 276.416 | >300 | — | — | TIMEOUT |
| 19,402 | 069 | 6,620/6,620 | 24.630 | 240.786 | 271.598 | -594,814,273.479574 | 9.921e-06 | PASS |
| 19,402 | 077 | 6,584/6,584 | 23.307 | 143.030 | 175.804 | -433,016,142.352673 | 8.113e-06 | PASS |
| 19,402 | 095 | 6,579/6,579 | 6.504 | 129.450 | 147.967 | -36,293,539.541442 | 9.913e-06 | PASS |

## Failure evidence

### `C2FEN02020` scenario `260`

The base succeeded. All 300 contingencies were screened; 296 passed the fast
screen and four entered the exact fallback queue. Two fallback tasks completed,
bringing the total to 298. The worker then emitted
`failed to construct the linearized AC HiGHS model`. The pending fast-screen
record for `CTG_0020` had a `0.3696233599` p.u. variable-bound residual at
`branch:1488:pf`. No complete official result exists.

### `C2FEN19402` scenario `006`

The base method `validated_sparse_newton_source_point_repair` returned
`success=false`. Its worst independent residual was `0.0015610134882357851`
p.u., category `variable_bound`, identity `branch:6104:pt`. No contingency was
started.

### `C2FEN19402` scenario `010`

The base took `24.1158495` seconds and the contingency phase took
`276.4156363` seconds. All 6,693 contingency solutions were produced, but the
controller reported `Code2 exhausted the time reserved for official
evaluation`. Because final exhaustive official evaluation did not complete,
there is no certified objective or infeasibility result.

## Launch-preflight record

The first launch command incorrectly applied the 19,402-case validation tail
sizes to the five 617-bus cases and to 2,020-bus scenario `025`. The runner
rejected those six launches because the requested completion-order shards could
not be formed. They produced no contingency result and are not counted among
the 37 scientific runs above. Their directories were retained. Replacement
launches changed only the validation shard partition to a size-valid fixed
taper; solver source, binary, mathematical checks, tolerances, and time limit
were unchanged.

## Artifact location and storage

The complete campaign is under
`runs/regression_freeze_c7fa432_20260826`. At completion it contained 168,377
files totaling approximately 107.823 GB. Windows `C:` retained approximately
151.05 GB free. Generated run artifacts remain ignored by Git; this report is
the tracked summary.
