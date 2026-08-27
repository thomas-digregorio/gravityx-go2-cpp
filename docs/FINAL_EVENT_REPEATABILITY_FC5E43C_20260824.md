# Final Event repeatability audit — `fc5e43c` — 2026-08-24

## Verdict

The unchanged current implementation passed 23 of the 24 requested cold
scenario runs. Every passing run completed below the 300-second end-to-end
limit, included every source contingency, passed the independent nonlinear
checker, and returned official infeasibility `0.0`.

`C2FEN08300` scenario `166` did not pass. Its base case was feasible, all 603
contingencies were screened, and 595 contingency solutions were completed,
but eight linearized corrective cases remained when the end-to-end work
deadline was reached. The official evaluator was therefore not started. This
is an incomplete timing result, not an infeasibility finding.

| Network | Requested | Full PASS | Incomplete/failed |
|---:|---:|---:|---:|
| 617 buses | 5 | 5 | 0 |
| 2,020 buses | 5 | 5 | 0 |
| 4,224 buses (`C2FEN04200`) | 7 | 7 | 0 |
| 8,300 buses | 7 | 6 | 1 |
| **Total** | **24** | **23** | **1** |

## Frozen execution boundary

- Branch: `codex/gravity-framework-v2`
- Commit: `fc5e43c9e5245c3b17932dbe95d3130417d3172d`
- Binary SHA-256:
  `95e4ce6310ff9047fb9e1e66fffc509e52ad3a529eb4024b9dda9baa50e55f0e`
- Official source archive SHA-256:
  `8a35cf683129fdcd318bac17e3c5cf98911f4887b8d8ce52de8999276b47f6f7`
- One cold run per scenario; no cross-scenario starts and no replacement runs.
- Hard boundary: 300 seconds from normalized-case loading through independent
  checking, official evaluation, and serialization.
- Pre-run checks: CTest `4/4` and pytest `12/12` passed.
- All raw inputs and outputs were kept under local `C:\Users\thoma\Documents`
  paths, outside OneDrive.

The source archive was reacquired from the official OpenEI Challenge 2 dataset
after the earlier raw extraction had been removed. Its archive hash and every
extracted `case.raw`, `case.con`, and `case.json` hash matched the registered
evidence before the 8,300-bus runs began.

## Results

Times below are final end-to-end process-return times. The baseline columns are
the previously retained published results. Objective change is
`(repeat - baseline) / abs(baseline)`; Challenge 2 maximizes this objective.
The retained baselines were produced by earlier implementation revisions, so
these comparisons are regression checks against prior evidence, not
same-binary timing repetitions.

### 617 buses

| Scenario | Contingencies | Baseline (s) | Repeat (s) | Repeat objective | Objective change | Max residual | Result |
|---:|---:|---:|---:|---:|---:|---:|---|
| 005 | 105 | 210.426 | 177.735 | 1,325,555.700484 | +0.0000% | 6.07e-8 | PASS |
| 017 | 103 | 213.297 | 171.192 | 1,193,434.865383 | -0.0000% | 4.27e-7 | PASS |
| 024 | 104 | 231.010 | 179.666 | 1,481,664.381598 | -0.0000% | 6.13e-8 | PASS |
| 062 | 103 | 193.352 | 174.469 | 1,416,741.921925 | +0.0601% | 7.75e-8 | PASS |
| 073 | 107 | 180.638 | 170.261 | 1,264,991.446413 | -0.0000% | 1.99e-8 | PASS |

Scenario 062 used one independently feasible state returned without a full
solver convergence status. It passed the explicit residual checks and the
official evaluator. The other four scenarios had zero such acceptances.

### 2,020 buses

| Scenario | Contingencies | Baseline (s) | Repeat (s) | Repeat objective | Objective change | Max residual | Result |
|---:|---:|---:|---:|---:|---:|---:|---|
| 025 | 292 | 216.767 | 94.560 | 3,499,217.938688 | -1.0815% | 3.40e-6 | PASS |
| 121 | 300 | 232.523 | 197.490 | 5,311,446.061792 | +0.2615% | 8.28e-6 | PASS |
| 134 | 300 | 226.131 | 202.224 | 995,579.276743 | +0.5512% | 1.43e-6 | PASS |
| 260 | 300 | 270.092 | 220.430 | 3,838,550.954070 | +1.4490% | 9.97e-6 | PASS |
| 262 | 300 | 230.922 | 200.300 | 4,377,078.568384 | +0.8746% | 8.92e-6 | PASS |

### 4,224 buses

| Scenario | Contingencies | Baseline (s) | Repeat (s) | Repeat objective | Objective change | Max residual | Result |
|---:|---:|---:|---:|---:|---:|---:|---|
| 009 | 455 | 244.622 | 123.685 | 1,469,108.401892 | +14.6063% | 9.61e-6 | PASS |
| 010 | 455 | 251.301 | 123.110 | 1,738,062.803125 | +12.0723% | 9.61e-6 | PASS |
| 014 | 455 | 252.510 | 120.454 | 1,466,165.132143 | +3.0614% | 9.49e-6 | PASS |
| 055 | 455 | 255.471 | 119.004 | 2,554,080.919971 | +7.3038% | 9.60e-6 | PASS |
| 056 | 455 | 253.314 | 115.331 | 2,819,891.095862 | +7.1978% | 9.61e-6 | PASS |
| 057 | 455 | 256.517 | 127.516 | 2,807,829.658402 | +5.3680% | 9.59e-6 | PASS |
| 060 | 455 | 258.220 | 120.392 | 2,584,708.088361 | +7.8179% | 9.61e-6 | PASS |

The current path accepted 428 or 429 contingencies through the fast screen and
sent only 26 or 27 to linearized correction. The retained baselines sent 62 to
64 cases to correction. This method-partition change accounts for much of the
large timing and objective difference.

### 8,300 buses

| Scenario | Contingencies | Baseline (s) | Repeat (s) | Repeat objective | Objective change | Max residual | Result |
|---:|---:|---:|---:|---:|---:|---:|---|
| 003 | 607 | 253.021 | 81.867 | 7,654,994.397662 | +0.9389% | 5.68e-14 | PASS |
| 012 | 607 | 260.215 | 83.586 | -3,223,328.612327 | +2.0793% | 5.68e-14 | PASS |
| 013 | 607 | 270.410 | 83.171 | -3,054,141.603864 | +2.1968% | 5.68e-14 | PASS |
| 022 | 607 | 255.892 | 84.121 | -756,951.998440 | +8.3799% | 5.68e-14 | PASS |
| 043 | 607 | 280.457 | 94.835 | 3,212,002.960929 | +7.7776% | 1.14e-13 | PASS |
| 052 | 607 | 261.047 | 94.176 | 5,124,197.229243 | +4.7488% | 1.14e-13 | PASS |
| 166 | 603 | 279.228 | 276.262 | — | — | — | **INCOMPLETE** |

For scenarios 003 through 022, 605 cases passed the fast path and two used
linearized correction. Scenarios 043 and 052 used 603 fast and four
linearized cases. The retained baselines used 33 or 34 linearized cases, again
explaining both the speedup and the different feasible corrective objectives.

## Scenario 166 failure record

- Robust base: feasible; 103.129 seconds.
- Source contingencies screened: 603/603.
- Fast-path feasible: 572.
- Registered linearized fallbacks: 31.
- Completed fallbacks: 23/31.
- Completed contingency solutions: 595/603.
- Written solution files: 596, including the base case.
- Elapsed time at clean stop: 276.262 seconds.
- Error: `contingency worker 4 reached the end-to-end work deadline`.
- Official evaluator details: 0/604; evaluation was correctly not started.

The runner reserves 24 seconds for the eight-process official evaluator and
one second for final serialization. Its work deadline is therefore about 275
seconds into the 300-second boundary. Stopping with eight cases unfinished was
the correct behavior: proceeding would have made an end-to-end PASS
impossible. The result does not establish that any contingency is infeasible.

## Acceptance audit

For each of the 23 passing runs:

- `run_status.json` reports `success=true`, stage `complete`, and the
  end-to-end gate satisfied;
- every source contingency and the base case have a solution and official
  evaluator detail;
- the official detail label set is exact and complete;
- official infeasibility is exactly `0.0`;
- the independent maximum residual is no greater than `1e-5`;
- raw, contingency, scenario, and normalized-model hashes match the retained
  scenario evidence.

The largest independent residual among the passing runs was
`9.967884830430762e-6` (2,020-bus scenario 260), below the acceptance limit.

## Interpretation and outstanding issues

1. **Functional repeatability:** 23 scenarios still solve and pass under the
   current implementation. Scenario 166 has a timing regression or tail
   variability that prevents a complete repeat inside 300 seconds.
2. **Objective repeatability:** the current run often returns a different
   feasible corrective state and therefore a different official objective.
   The largest change is +14.6063% for 4,224-bus scenario 009. The evidence
   points to the newer fast/direct sparse-Newton acceptance path replacing many
   former linearized or local optimization solves.
3. **What is not shown:** because the retained baselines use earlier commits
   and each current scenario was run only once, this campaign does not measure
   same-commit stochastic reproducibility. It establishes current functional
   pass/fail status against the frozen source cases.
4. **No automatic repair:** no scenario was rerun, no tolerance was changed,
   and no source or algorithm file was modified after the failure.

## Storage maintenance

Before the campaign, the Ubuntu `ext4.vhdx` was trimmed, WSL was shut down,
and the detached VHDX was compacted. Its physical size fell from 461.072 GiB
to about 170.933 GiB, returning about 290.14 GiB to Windows. C: had about
291.22 GiB free immediately afterward and 278.04 GiB free after the retained
repeat artifacts were written. WSL was stopped again after the campaign.

All repeat artifacts are under:

`C:\Users\thoma\Documents\gravityx-go2-cpp\runs\repeatability_fc5e43c_20260824`

