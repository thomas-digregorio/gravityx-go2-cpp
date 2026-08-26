# Final Event six-family frozen regression — `9010777` — 2026-08-26

## Verdict

The latest frozen algorithm passed **36 of 37** cold scenario runs across all
six requested network families. It passed every 617-, 4,224-, 8,300-,
16,789-, and 19,402-bus scenario, plus four of five 2,020-bus scenarios.

The sole failure was 2,020-bus scenario 260. It completed 298 of 300 source
contingencies, then the serial corrective worker exited while constructing the
linearized AC HiGHS fallback model for the next queued fallback
(`CTG_0020`). This was a Code 2 implementation failure after 10.699 seconds,
not a timeout and not a certified grid infeasibility. It was not retried.

Across the campaign, 44,515 of 44,517 requested contingency executions
completed. The 36 successful scenarios independently certified 44,217 source
contingencies. Every successful run had a complete label set, official
infeasibility `0.0`, exact agreement between its reported objective and
detail-derived objective, and a maximum independent residual no larger than
`9.99584980010404e-6`.

## Frozen boundary

- Git tag: `six-family-regression-freeze-20260826-v1`
- Recorded revision: `901077733fbf465a5d85c61fbab32976a43f99a9`
- Algorithm source revision: `b8cb79ec82e3e5d2ff7eb23745cd0f424fd5e6de`
- Standard exact executable SHA-256:
  `e870443dd975d1373246d1bf9beacc5d1094711bcf19f458a0157f4837f4be05`
- Native executable SHA-256:
  `5689c9fcbfe2461424b07fe73d1a045ecbcbc67e2cb151c5b7cf5588f9158495`
- Scenario-010 timing-only profile SHA-256:
  `d739ae2d353e762e61c2e9de6f2660548f056821abd5e26799f658d984d260aa`
- Pre-run compiled tests: standard CTest `4/4`; native CTest `4/4`.
- One cold invocation per scenario, no retry, no executable rebuild, and no
  source or configuration change during the campaign.
- Hard end-to-end limit: 300 seconds per scenario.
- Scenario 006 in the 19,402-bus family used the frozen standard executable
  for the exact base solve and the native executable for contingency
  screening. Every other scenario used the frozen native executable for both.
- Scenario 010 in the 19,402-bus family used the frozen `full_v6` timing-only
  work-order profile. It contains no prior primal, dual, network, commitment,
  or solver state.

The repository was clean after the campaign, the freeze tag still resolved to
the recorded revision, and all three frozen hashes still matched.

## Family summary

`C2FEN04200` parses to 4,224 buses and `C2FEN17700` parses to 16,789 buses.
The timing range for the 2,020-bus family includes the failed scenario 260.

| Parsed buses | Source family | Scenarios passed | Completed/requested contingencies | End-to-end range (s) | Sum of scenario boundaries (s) |
|---:|---|---:|---:|---:|---:|
| 617 | C2FEN00617 | 5/5 | 522/522 | 11.290–14.985 | 62.114 |
| 2,020 | C2FEN02020 | 4/5 | 1,490/1,492 | 10.699–165.782 | 254.974 |
| 4,224 | C2FEN04200 | 7/7 | 3,185/3,185 | 127.949–131.532 | 905.741 |
| 8,300 | C2FEN08300 | 7/7 | 4,245/4,245 | 37.940–242.570 | 495.014 |
| 16,789 | C2FEN17700 | 8/8 | 1,904/1,904 | 37.328–44.595 | 327.449 |
| 19,402 | C2FEN19402 | 5/5 | 33,169/33,169 | 140.048–285.442 | 1,098.711 |
| **Total** | — | **36/37** | **44,515/44,517** | **10.699–285.442** | **3,144.003** |

The summed scenario boundaries are 52.400 minutes. They exclude setup,
between-run evidence handling, and artifact pruning.

## Detailed results

Times are serialized boundaries from each `run_status.json` and
`run_summary.json`. Streaming official evaluation overlaps contingency work,
so base and contingency columns should not be added to reconstruct end-to-end
time.

| Buses | Scenario | Completed/source | Base (s) | Contingency (s) | End-to-end (s) | Official objective | Max independent residual | Exact fallbacks | Result |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 617 | 005 | 105/105 | 3.564 | 9.470 | 14.985 | 837,528.141490 | 2.220e-16 | 0 | PASS |
| 617 | 017 | 103/103 | 0.109 | 9.985 | 12.213 | 677,959.663942 | 1.110e-16 | 0 | PASS |
| 617 | 024 | 104/104 | 0.094 | 9.110 | 11.290 | 1,024,779.091343 | 4.441e-16 | 0 | PASS |
| 617 | 062 | 103/103 | 0.098 | 9.233 | 11.973 | 973,222.118662 | 2.220e-16 | 0 | PASS |
| 617 | 073 | 107/107 | 0.095 | 9.205 | 11.653 | 758,559.216716 | 4.441e-16 | 0 | PASS |
| 2,020 | 025 | 292/292 | 3.427 | 9.883 | 15.738 | 3,984,584.198397 | 8.882e-16 | 0 | PASS |
| 2,020 | 121 | 300/300 | 0.233 | 162.529 | 165.782 | -1,745,787.467477 | 1.872e-6 | 43 | PASS |
| 2,020 | 134 | 300/300 | 0.179 | 9.965 | 12.391 | 4,735,923.219387 | 1.776e-15 | 0 | PASS |
| 2,020 | 260 | 298/300 | 0.182 | — | 10.699 | — | — | — | **FAIL: Code 2 worker** |
| 2,020 | 262 | 300/300 | 0.268 | 47.422 | 50.365 | -3,013,880.273921 | 1.421e-14 | 9 | PASS |
| 4,224 | 009 | 455/455 | 17.877 | 109.555 | 131.532 | 1,079,986.571967 | 2.842e-14 | 15 | PASS |
| 4,224 | 010 | 455/455 | 15.760 | 108.445 | 128.046 | 1,335,684.862638 | 2.842e-14 | 15 | PASS |
| 4,224 | 014 | 455/455 | 15.768 | 110.534 | 130.104 | 1,058,163.791479 | 7.105e-15 | 15 | PASS |
| 4,224 | 055 | 455/455 | 15.513 | 110.387 | 129.677 | 2,093,982.002567 | 1.421e-14 | 15 | PASS |
| 4,224 | 056 | 455/455 | 15.566 | 109.024 | 128.452 | 2,380,103.181476 | 1.421e-14 | 15 | PASS |
| 4,224 | 057 | 455/455 | 15.611 | 108.524 | 127.949 | 2,448,049.880033 | 2.842e-14 | 15 | PASS |
| 4,224 | 060 | 455/455 | 15.475 | 110.631 | 129.982 | 2,017,696.701851 | 1.421e-14 | 15 | PASS |
| 8,300 | 003 | 607/607 | 3.023 | 37.900 | 50.479 | -13,609,513.240398 | 6.256e-6 | 1 | PASS |
| 8,300 | 012 | 607/607 | 0.652 | 34.987 | 43.427 | -24,501,716.970451 | 6.256e-6 | 1 | PASS |
| 8,300 | 013 | 607/607 | 0.727 | 33.138 | 39.079 | -24,332,327.798296 | 6.256e-6 | 1 | PASS |
| 8,300 | 022 | 607/607 | 0.702 | 32.824 | 41.136 | -22,032,420.321919 | 6.256e-6 | 1 | PASS |
| 8,300 | 043 | 607/607 | 0.710 | 31.501 | 40.382 | -26,809,846.003070 | 5.684e-14 | 1 | PASS |
| 8,300 | 052 | 607/607 | 0.737 | 32.083 | 37.940 | -24,902,976.143042 | 5.684e-14 | 1 | PASS |
| 8,300 | 166 | 603/603 | 0.851 | 236.486 | 242.570 | -22,389,627.637575 | 1.705e-13 | 10 | PASS |
| 16,789 | 019 | 238/238 | 3.650 | 35.159 | 44.595 | -77,574,732.536182 | 2.720e-6 | 0 | PASS |
| 16,789 | 020 | 238/238 | 1.359 | 34.234 | 41.350 | -77,474,917.225687 | 2.720e-6 | 0 | PASS |
| 16,789 | 021 | 238/238 | 1.447 | 34.023 | 41.253 | -77,591,722.692500 | 2.720e-6 | 0 | PASS |
| 16,789 | 089 | 238/238 | 1.262 | 29.999 | 39.305 | -77,682,401.183982 | 3.301e-6 | 0 | PASS |
| 16,789 | 094 | 238/238 | 1.289 | 25.054 | 37.328 | -77,628,114.043451 | 8.149e-6 | 0 | PASS |
| 16,789 | 106 | 238/238 | 1.372 | 33.507 | 40.720 | -77,398,719.037299 | 1.137e-13 | 0 | PASS |
| 16,789 | 107 | 238/238 | 1.247 | 32.515 | 39.577 | -77,604,018.875521 | 1.137e-13 | 0 | PASS |
| 16,789 | 115 | 238/238 | 1.537 | 33.823 | 43.320 | -80,045,119.774033 | 2.844e-6 | 0 | PASS |
| 19,402 | 006 | 6,693/6,693 | 27.693 | 218.763 | 257.352 | -679,351,881.006607 | 9.996e-6 | 0 | PASS |
| 19,402 | 010 | 6,693/6,693 | 26.581 | 248.893 | 285.442 | -582,029,745.991480 | 9.704e-6 | 0 | PASS |
| 19,402 | 069 | 6,620/6,620 | 26.160 | 214.310 | 247.694 | -594,814,154.478509 | 9.921e-6 | 0 | PASS |
| 19,402 | 077 | 6,584/6,584 | 25.642 | 130.977 | 168.175 | -433,016,057.964859 | 8.113e-6 | 0 | PASS |
| 19,402 | 095 | 6,579/6,579 | 8.812 | 118.815 | 140.048 | -36,294,106.898381 | 9.913e-6 | 0 | PASS |

## Scenario 260 failure evidence

The 2,020-bus scenario-260 base case succeeded. The fast stage screened all
300 contingencies and accepted 298. The serial corrective worker then solved
fallbacks `CTG_0035` and `CTG_0001`, after which its log ended with:

```text
error: failed to construct the linearized AC HiGHS model
```

The next queued fallback was `CTG_0020`. Its retained fast-screen record shows
that it required exact fallback after independent validation found a branch
active-flow variable-bound residual at `branch:1488:pf`. Because the worker
exited before producing the final two cases, no official objective,
infeasibility score, or security certificate exists for scenario 260. It must
therefore remain a failure, not be labeled feasible or infeasible.

## Evidence and storage handling

- Local run root:
  `runs/six_family_regression_9010777_20260826`
- Compact provenance archive:
  `docs/evidence/SIX_FAMILY_REGRESSION_9010777_PROVENANCE_20260826.json`
- Provenance archive SHA-256:
  `7c6576c99219f66014cb3f2a7e01463cf1929628da442bca7dd5b8bf6c90a97a`
- The archive contains all 37 run-status records and hash-backed provenance.
- Per-run summaries, failure logs, contingency records, validated evaluator
  shards, and official certificates remain local.
- To avoid exhausting the laptop disk, 90,791 generated `solution_*.txt`
  files totaling approximately 116.593 GB were deleted only after each run's
  retained result/certificate was checked. No source, configuration, binary,
  report, or failure evidence was deleted. Recreating those exact solution
  text files would require rerunning the corresponding scenarios.
