# Final Event six-family frozen regression — `c477a2d` — 2026-08-26

## Verdict

The frozen algorithm passed **37 of 37** intended cold scenario runs across all
six requested network families: 617, 2,020, 4,224, 8,300, 16,789, and 19,402
buses. It evaluated and independently certified all **44,517 of 44,517**
source contingencies. Every run had a complete official label set, official
infeasibility `0.0`, exact agreement between the reported and detail-derived
objective, and an independent residual no larger than
`9.99584980010404e-6`, within the frozen `1e-5` acceptance tolerance.

The summed end-to-end scenario boundaries were **3,128.085 seconds
(52:08.085)**. Every individual run completed within the frozen 300-second
limit. The slowest was 19,402-bus scenario 010 at **288.611 seconds**.

## Surgical improvement and regression controls

The 19,402-bus scenario-010 problem was deadline overhead, not a security or
grid-feasibility failure. A prior run completed all 6,693 contingency solves,
but repeated copies of large immutable contingency state left insufficient
time for the last official-evaluation shard before the 300-second boundary.

The final change uses an explicit borrowed view of the immutable base state,
moves completed states rather than copying them, and defers construction of
fallback state until fallback is actually needed. Ownership and lifetime are
explicit; a context can still own a base state when required. This changes
state transport, not the power-system formulation, source contingencies,
solver tolerances, feasibility checks, or official evaluator.

The safeguards against regression were:

- component tests for owned, borrowed, and copied-borrowed contexts;
- standard and native CTest suites, both `4/4` PASS, plus Python `48/48` PASS;
- a fixed 48-contingency equivalence test with the same zero-fallback result
  and the same worst residual (`7.869767e-6`), while elapsed time fell from
  56.40 to 51.81 seconds;
- one exact frozen revision and executable hashes for the complete campaign;
- complete official evaluation plus independent residual verification for
  every source contingency in every scenario.

The formerly failing 2,020-bus scenario 260 also passed. Its old failure was a
Code 2 implementation/status-handling problem: a recoverable HiGHS warning was
treated as a fatal linear-seed construction failure. The correction preserves
valid warning-return results while genuine solver errors remain failures.
Scenario 260 now completed all 300 contingencies in 23.957 seconds with zero
official infeasibility and a `7.105e-15` maximum independent residual.

## Frozen boundary

- Algorithm revision: `c477a2d4673a87d4d569110471956e420bb919e6`
- Candidate tag: `six-family-regression-candidate-20260826-v2`
- Native executable SHA-256:
  `670e624a8eff93f1282777f7b8cbacbf6ff8373c10b4ad1f3eb6d2b08d49c9b7`
- Standard executable SHA-256:
  `03d42ae80f3ba685614db01a3c2f015bf13f38e732c5e6167b1fa5cd6a8ba598`
- One cold successful invocation per intended scenario; no warm start from a
  different scenario and no source/configuration change during the campaign.
- Hard end-to-end limit: 300 seconds per scenario.
- `C2FEN04200` parses to 4,224 buses; `C2FEN17700` parses to 16,789 buses.
- The 19,402-bus scenario-006 exact base solve used the standard executable
  and its contingency screen used the native executable. Other successful
  runs used the native executable for both boundaries.

## Family summary

| Parsed buses | Source family | Scenarios passed | Certified contingencies | End-to-end range (s) | Sum (s) | Max residual |
|---:|---|---:|---:|---:|---:|---:|
| 617 | C2FEN00617 | 5/5 | 522/522 | 11.960–14.903 | 63.333 | 4.441e-16 |
| 2,020 | C2FEN02020 | 5/5 | 1,492/1,492 | 13.443–172.187 | 276.564 | 1.872e-6 |
| 4,224 | C2FEN04200 | 7/7 | 3,185/3,185 | 129.046–137.481 | 927.028 | 2.842e-14 |
| 8,300 | C2FEN08300 | 7/7 | 4,245/4,245 | 38.262–249.185 | 494.872 | 6.256e-6 |
| 16,789 | C2FEN17700 | 8/8 | 1,904/1,904 | 34.333–40.484 | 301.502 | 8.149e-6 |
| 19,402 | C2FEN19402 | 5/5 | 33,169/33,169 | 145.064–288.611 | 1,064.786 | 9.996e-6 |
| **Total** | — | **37/37** | **44,517/44,517** | **11.960–288.611** | **3,128.085** | **9.996e-6** |

## Detailed results

Times are retained `run_summary.json` boundaries. Streaming official
evaluation overlaps contingency work, so the base and contingency columns
must not be added to reconstruct end-to-end time.

| Buses | Scenario | Contingencies | Base (s) | Contingency (s) | End-to-end (s) | Official objective | Max independent residual | Exact fallbacks | Result |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 617 | 005 | 105 | 2.654 | 10.054 | 14.903 | 837,528.141490 | 2.220e-16 | 0 | PASS |
| 617 | 017 | 103 | 0.102 | 9.895 | 11.960 | 677,959.663942 | 1.110e-16 | 0 | PASS |
| 617 | 024 | 104 | 0.096 | 9.876 | 12.131 | 1,024,779.091343 | 4.441e-16 | 0 | PASS |
| 617 | 062 | 103 | 0.103 | 9.881 | 12.250 | 973,222.118662 | 2.220e-16 | 0 | PASS |
| 617 | 073 | 107 | 0.101 | 10.014 | 12.089 | 758,559.216716 | 4.441e-16 | 0 | PASS |
| 2,020 | 025 | 292 | 2.479 | 10.533 | 15.735 | 3,984,584.198397 | 8.882e-16 | 0 | PASS |
| 2,020 | 121 | 300 | 0.261 | 169.811 | 172.187 | -1,745,787.467477 | 1.872e-6 | 43 | PASS |
| 2,020 | 134 | 300 | 0.183 | 10.752 | 13.443 | 4,735,923.219387 | 1.776e-15 | 0 | PASS |
| 2,020 | 260 | 300 | 0.169 | 21.424 | 23.957 | 5,736,074.714256 | 7.105e-15 | 4 | PASS |
| 2,020 | 262 | 300 | 0.286 | 48.656 | 51.242 | -3,013,880.273921 | 1.421e-14 | 9 | PASS |
| 4,224 | 009 | 455 | 18.825 | 115.722 | 137.481 | 1,079,986.571967 | 2.842e-14 | 15 | PASS |
| 4,224 | 010 | 455 | 16.212 | 114.613 | 134.846 | 1,335,684.862638 | 2.842e-14 | 15 | PASS |
| 4,224 | 014 | 455 | 15.882 | 111.928 | 131.549 | 1,058,163.791479 | 7.105e-15 | 15 | PASS |
| 4,224 | 055 | 455 | 15.524 | 113.561 | 132.408 | 2,093,982.002567 | 1.421e-14 | 15 | PASS |
| 4,224 | 056 | 455 | 16.951 | 111.585 | 132.120 | 2,380,103.181476 | 1.421e-14 | 15 | PASS |
| 4,224 | 057 | 455 | 15.792 | 110.314 | 129.578 | 2,448,049.880033 | 2.842e-14 | 15 | PASS |
| 4,224 | 060 | 455 | 16.038 | 109.583 | 129.046 | 2,017,696.701851 | 1.421e-14 | 15 | PASS |
| 8,300 | 003 | 607 | 0.679 | 34.373 | 39.804 | -13,609,513.240398 | 6.256e-6 | 1 | PASS |
| 8,300 | 012 | 607 | 0.643 | 34.895 | 40.278 | -24,501,716.970451 | 6.256e-6 | 1 | PASS |
| 8,300 | 013 | 607 | 0.657 | 34.100 | 42.859 | -24,332,327.798296 | 6.256e-6 | 1 | PASS |
| 8,300 | 022 | 607 | 0.653 | 34.846 | 43.119 | -22,032,420.321919 | 6.256e-6 | 1 | PASS |
| 8,300 | 043 | 607 | 0.762 | 32.250 | 41.366 | -26,809,846.003070 | 5.684e-14 | 1 | PASS |
| 8,300 | 052 | 607 | 0.700 | 32.739 | 38.262 | -24,902,976.143042 | 5.684e-14 | 1 | PASS |
| 8,300 | 166 | 603 | 0.775 | 243.315 | 249.185 | -22,389,627.637575 | 1.705e-13 | 10 | PASS |
| 16,789 | 019 | 238 | 3.587 | 29.222 | 38.210 | -77,574,836.488953 | 2.720e-6 | 0 | PASS |
| 16,789 | 020 | 238 | 1.250 | 30.102 | 37.172 | -77,474,917.225687 | 2.720e-6 | 0 | PASS |
| 16,789 | 021 | 238 | 1.326 | 30.044 | 39.672 | -77,591,749.858990 | 2.720e-6 | 0 | PASS |
| 16,789 | 089 | 238 | 1.325 | 23.991 | 36.906 | -77,682,401.183982 | 3.301e-6 | 0 | PASS |
| 16,789 | 094 | 238 | 1.275 | 24.903 | 34.333 | -77,628,200.541403 | 8.149e-6 | 0 | PASS |
| 16,789 | 106 | 238 | 1.440 | 27.736 | 37.818 | -77,398,719.037299 | 1.137e-13 | 0 | PASS |
| 16,789 | 107 | 238 | 1.307 | 28.526 | 36.907 | -77,604,143.844662 | 1.137e-13 | 0 | PASS |
| 16,789 | 115 | 238 | 1.619 | 30.788 | 40.484 | -80,046,323.928941 | 2.844e-6 | 0 | PASS |
| 19,402 | 006 | 6,693 | 28.289 | 194.970 | 233.536 | -679,351,812.814108 | 9.996e-6 | 0 | PASS |
| 19,402 | 010 | 6,693 | 31.482 | 251.305 | 288.611 | -582,029,795.284929 | 9.704e-6 | 0 | PASS |
| 19,402 | 069 | 6,620 | 31.726 | 185.629 | 227.178 | -594,813,806.702179 | 9.921e-6 | 0 | PASS |
| 19,402 | 077 | 6,584 | 36.594 | 120.497 | 170.397 | -433,016,325.710104 | 8.113e-6 | 0 | PASS |
| 19,402 | 095 | 6,579 | 18.105 | 113.958 | 145.064 | -36,295,267.953619 | 9.913e-6 | 0 | PASS |

## Setup preflight and evidence handling

The first 617-bus scenario-005 directory is a setup preflight, not one of the
37 experiments. It rejected before the base solve because the initial command
used incorrect evaluator-tail sizes. No optimization ran. The corrected cold
scenario-005 invocation is retained as `C2FEN00617_s005_cold_replacement` and
is the sole scenario-005 result in the tables.

- Main local run root:
  `runs/six_family_regression_c477a2d_20260826`
- Scenario-010 local run root:
  `runs/six_family_retest_c477a2d_20260826/C2FEN19402_s010_cold`
- Main compact provenance archive (36 intended successes plus the rejected
  preflight):
  `docs/evidence/SIX_FAMILY_REGRESSION_C477A2D_PROVENANCE_20260826.json`
- Main archive SHA-256:
  `5121a5f4d17a9ad51a79dd9505a38d2bb9930b0c28566a116a1fb141493da9d3`
- Scenario-010 provenance archive:
  `docs/evidence/SIX_FAMILY_REGRESSION_C477A2D_EXTERNAL_S010_PROVENANCE_20260826.json`
- Scenario-010 archive SHA-256:
  `a43e308f9617b7f3dc572ad423b1b137fb0df7debf04a3fc1938e5ed4a06517a`

Large generated solution payloads were deleted only after their summaries,
statuses, official certificates, logs, and hashes were retained and audited.
No source data, configuration, executable, report, or certification evidence
was deleted.
