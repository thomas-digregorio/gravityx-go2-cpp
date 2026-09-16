# V34: all 37 scenarios PASS

All 37 scenarios across the six requested network families passed one cold
run each on the same frozen V34 algorithm. This replaces the incomplete
V33 result of 36/37. No benchmark repetitions, source changes, omitted
contingencies, final-tolerance relaxation or deadline extension were used.

## Frozen implementation and acceptance

- Algorithm/run revision: `3b335872647296c2d0e2c8ae5acbe4359ae72f0a`.
- Executable SHA256: `353439cb4b29a27b0a75cf4a006a894c65adf7bf5d0f06f8bd193bd21ab93efa`.
- Configuration: `config/v34_all_scenarios_20260916.json`.
- Configuration SHA256: `9a50c0f540c4633dce6c6f9cd2e15fabbc6b99e2228f34644a4e167ed8a2daef`.
- Component checks before full runs: 108 Python tests and four native CTest groups passed.
- Every scenario started from its own source inputs, with no previous run's solution imported.
- All source scenario directories in the six families are covered: 37 unique cases.
- All 44,517 source contingencies, plus the 37 base cases, completed verification.
- Every official infeasibility value is 0.0.
- Largest native residual: 9.990036206664055e-6, below the unchanged 1e-5 limit.
- Largest measured run time: 268.0952593999973 seconds, below the unchanged 300-second cap.

Times cover normalized-case loading through referenced official-evaluator
equations and result serialization. Build/setup and prior normalization
are outside the registered boundary. Work cutoff remains 295 seconds,
reserving three seconds for evaluation and two for finalization.
Hash archival and this post-run evidence audit are not additional solves.

PASS certifies complete feasibility checks under the existing GO2-derived
model, which permits penalized imbalances. It does not mean zero penalized
imbalance, a global optimum, a certified optimality gap, or top-five
competition placement. Higher official objective is better. These single
runs do not establish a statistically measured speedup, and not every
scenario's objective improved.

## What fixed scenario 089

The old repair could clear a signed balance pair before discovering that
its reconstructed value exceeded a bound by around 2e-11. Rejecting that
roundoff left an artificial 0.5 row violation in the diagnostic. Reconstruction
is now atomic, keeps the exact column bounds, and retains the existing
full 1e-7 LP row audit. Nonfinite or material errors are still rejected.
Tiny tests cover both signs, fixed columns, and preservation on failure.

The bounded expanded repair also previously stopped after a non-improving
full step. It now retains the previous point and retries with quarter-sized
angle/voltage trust radii. Attempt count, neighborhood size, per-LP time
and total repair time remain explicitly bounded. The final independent
nonlinear acceptance check is unchanged. There was no separate full-case
ablation to attribute the improvement between these two fixes.

| 16,789-family scenario 089 | V33 | V34 |
|---|---:|---:|
| Result | Work-deadline timeout | PASS |
| Completed contingencies | 237/238 | 238/238 |
| End-to-end seconds | 295.1255333000008 | 84.17479409999942 |
| Complete official objective | None | 11,455,738.985055316 |
| Official infeasibility | No complete certificate | 0.0 |
| Maximum native residual for complete case | Not certified | 8.414052062200028e-6 |

V34 required no slow corrective fallback for 089. This was the first V34
cold run and counts toward the 37; it was not repeated after passing.
It was followed by all five 19,402-family cases and the other 31 smaller
cases, without changing the algorithm or configuration.

## Results by family

| Network family | Passed / tested | End-to-end range (s) | Source contingencies |
|---|---:|---:|---:|
| 617 | 5/5 | 2.98–6.03 | 522 |
| 2,020 | 5/5 | 10.70–16.47 | 1,492 |
| 4,224 | 7/7 | 49.38–57.16 | 3,185 |
| 8,300 | 7/7 | 115.73–173.58 | 4,245 |
| 16,789 | 8/8 | 71.13–84.17 | 1,904 |
| 19,402 | 5/5 | 214.32–268.10 | 33,169 |

## Every scenario

These are V34 results, not a mixture of best historical runs.
Unrounded values are in `evidence/V34_ALL_SCENARIOS_SUMMARY_20260916.json`.

| Network family | Scenario | Contingencies | End-to-end (s) | Official-evaluator objective | Maximum native residual | Result |
|---|---|---:|---:|---:|---:|---|
| 617 | 005 | 105/105 | 6.03 | 1,232,415.20 | 2.71374030e-8 | PASS |
| 617 | 017 | 103/103 | 2.98 | 1,044,179.56 | 3.58609111e-8 | PASS |
| 617 | 024 | 104/104 | 3.34 | 1,434,134.94 | 3.44472100e-8 | PASS |
| 617 | 062 | 103/103 | 3.10 | 1,346,818.18 | 3.56636893e-8 | PASS |
| 617 | 073 | 107/107 | 3.25 | 1,131,781.93 | 3.42950903e-8 | PASS |
| 2,020 | 025 | 292/292 | 10.70 | 4,039,032.52 | 1.96553808e-6 | PASS |
| 2,020 | 121 | 300/300 | 10.73 | 6,219,881.67 | 3.95226541e-6 | PASS |
| 2,020 | 134 | 300/300 | 10.90 | 4,792,821.26 | 4.44089210e-16 | PASS |
| 2,020 | 260 | 300/300 | 10.91 | 5,892,570.70 | 1.77635684e-15 | PASS |
| 2,020 | 262 | 300/300 | 16.47 | 5,652,949.73 | 9.16317877e-6 | PASS |
| 4,224 | 009 | 455/455 | 56.98 | 4,959,833.45 | 7.90348977e-8 | PASS |
| 4,224 | 010 | 455/455 | 57.16 | 5,229,164.71 | 5.84174842e-9 | PASS |
| 4,224 | 014 | 455/455 | 54.26 | 4,943,537.42 | 6.24819521e-8 | PASS |
| 4,224 | 055 | 455/455 | 49.38 | 6,008,746.35 | 4.24446012e-9 | PASS |
| 4,224 | 056 | 455/455 | 56.90 | 6,249,036.64 | 7.76761100e-9 | PASS |
| 4,224 | 057 | 455/455 | 56.18 | 6,326,034.06 | 2.01926227e-8 | PASS |
| 4,224 | 060 | 455/455 | 54.62 | 5,879,776.79 | 5.35703695e-8 | PASS |
| 8,300 | 003 | 607/607 | 152.35 | 19,524,401.69 | 9.34397091e-6 | PASS |
| 8,300 | 012 | 607/607 | 119.58 | 8,546,719.57 | 9.92694571e-6 | PASS |
| 8,300 | 013 | 607/607 | 115.73 | 8,586,981.10 | 9.24326474e-6 | PASS |
| 8,300 | 022 | 607/607 | 121.05 | 10,952,868.05 | 9.65750739e-6 | PASS |
| 8,300 | 043 | 607/607 | 138.00 | 8,930,999.63 | 9.94428342e-6 | PASS |
| 8,300 | 052 | 607/607 | 140.53 | 10,670,225.91 | 9.95185428e-6 | PASS |
| 8,300 | 166 | 603/603 | 173.58 | 9,943,280.04 | 9.82314900e-6 | PASS |
| 16,789 | 019 | 238/238 | 75.25 | 10,912,788.55 | 9.67645424e-6 | PASS |
| 16,789 | 020 | 238/238 | 73.42 | 11,154,375.47 | 6.52903053e-6 | PASS |
| 16,789 | 021 | 238/238 | 77.12 | 10,928,554.56 | 6.41552307e-6 | PASS |
| 16,789 | 089 | 238/238 | 84.17 | 11,455,738.99 | 8.41405206e-6 | PASS |
| 16,789 | 094 | 238/238 | 71.13 | 11,635,258.31 | 5.50573451e-6 | PASS |
| 16,789 | 106 | 238/238 | 82.66 | 10,572,549.61 | 9.18043304e-6 | PASS |
| 16,789 | 107 | 238/238 | 74.37 | 10,429,320.12 | 8.42078478e-6 | PASS |
| 16,789 | 115 | 238/238 | 71.15 | 12,261,185.92 | 7.49720861e-6 | PASS |
| 19,402 | 006 | 6,693/6,693 | 264.52 | 335,287.72 | 9.95223635e-6 | PASS |
| 19,402 | 010 | 6,693/6,693 | 268.10 | 471,449.39 | 9.99003621e-6 | PASS |
| 19,402 | 069 | 6,620/6,620 | 234.55 | 888,780.09 | 9.92728923e-6 | PASS |
| 19,402 | 077 | 6,584/6,584 | 214.32 | 376,772.83 | 4.16777614e-6 | PASS |
| 19,402 | 095 | 6,579/6,579 | 234.87 | 402,081.43 | 9.98385224e-6 | PASS |

Network labels are the familiar family labels. The normalized input has
8,316 buses for each 8,300-family scenario. In the 19,402 family, 006/010
have 19,402 buses, 069 has 18,916, 077 has 18,889, and 095 has 18,877.
Other displayed labels equal the audited normalized bus counts. Every
input hash matches its frozen manifest; this distinction changes only
the reporting label.

## Independent post-run audit

The read-only audit rechecked every acceptance gate against saved status,
summary and official certificate, all source RAW/JSON/CON and normalized
model hashes, exact cold commands, executable/configuration hashes, source
scenario coverage, uniqueness, per-run archive/status/adjacent-JSON hashes,
and full output-vector counts or the explicitly recorded pruning exception.

The audit returned 37 passes, zero failures, 44,517 complete source
contingencies, all hashes matching, and all cold starts. No optimization
processes remain active.

Campaign archive:
`evidence/V34_ALL_SCENARIOS_PROVENANCE_20260916.json`,
4,179,793 bytes, 37 runs, no archive diagnostics.
SHA256: `30a154253f7cb9ab93608de279e8bb5c17871777b27f38b1985307723bc62a4e`.

## Storage and exact-output retention

Only solution text payloads from five redundant V33/V34 run pairs were
removed, after complete same-input verification, archive/hash checks,
resolved local non-reparse path checks, and confirmation that workers had
stopped. Source data, code, environments, logs, status, summaries,
certificates and compact evidence were retained.

| 19,402-family scenario | Pruned solution-text run | Retained complete run | Removed text entries | Physical GiB recovered |
|---|---|---|---:|---:|
| 006 | V33 | V34 | 13,434 | 11.48 |
| 010 | V33 | V34 | 13,434 | 11.49 |
| 069 | V33 | V34 | 13,288 | 11.01 |
| 077 | V34 | V33 | 13,216 | 10.95 |
| 095 | V33 | V34 | 13,206 | 10.94 |

Total physical space recovered by those operations was
55.879825592041016 GiB; newly generated campaign outputs reused some space.
Windows C: free space after the campaign/archive was 32.24 GiB
(34,621,403,136 bytes). Individual pre/post byte counts and archive hashes
are in the five `evidence/V34_PRUNE_19402_S*_20260916.json` receipts.

Full V34 solutions remain for 36 scenarios. For 19,402-family 077, the
higher-scoring complete V33 solution (376,775.64657868724) remains, while
the V34 solution texts (376,772.8285620564) were pruned after verification.
The table above deliberately reports the actual V34 score and time.
Its compact verification record remains, but cannot reconstruct or
independently re-evaluate those exact deleted V34 text outputs.

The four pruned V33 solutions were 006, 010, 069 and 095; their complete
V34 replacements remain. Thus each scenario still has a complete solution
retained, and all V34 verification records are preserved.
