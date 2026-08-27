# C2FEN08300 all-scenario regression — `eab4f29` — 2026-08-24

## Outcome

All seven public Final Event 4 `C2FEN08300` scenarios passed one cold run on
revision `eab4f29eac292b1093f7188040f4911373e6c842`. Every source contingency
and the base case produced a solution, passed the independent nonlinear
validator, and appeared in the official evaluator detail set. Official
aggregate infeasibility was `0.0` in every scenario.

The frozen executable SHA-256 was
`b73b23af4592a963269801341fb8d7af7d2c52c2af5730ed49de949ad4ee52fa`.
Before the retained runs, CTest passed 4/4 and the Python component suite
passed 12/12.

## Corrections

The mathematical model, source contingencies, feasibility tolerance, and
official evaluator were not changed. Three orchestration/solver-control
corrections addressed the scenario-166 tail:

1. Measured fallback durations now feed one global longest-first priority
   queue. Historical worker assignments no longer trap work in per-worker FIFO
   queues.
2. The projected-balance Phase-I LP on 8k-bus retries receives a logged
   90-second internal limit. Two observed cases had reached the prior
   60-second cutoff after 19 IPM iterations, while a valid solution required
   22. The per-worker and 300-second end-to-end deadlines remain hard limits.
3. The corrective pool starts with eight workers while sixteen fast screeners
   are active. After exhaustive screening finishes and those screeners exit,
   the pool expands to twelve workers, reusing four freed CPU cores.

No failed solver state or nonconverged heuristic state was accepted.

## Results

The end-to-end column uses the conservative final process-return timing. Stage
times come from the persisted run summary and can differ slightly because of
loading, orchestration, and final serialization.

| Scenario | Source contingencies | Base (s) | Code 2 (s) | Evaluation (s) | End-to-end (s) | Official objective | Official infeasibility | Max independent residual | Fast / linearized | Result |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 003 | 607 | 3.679 | 38.445 | 40.077 | **83.066** | 7,654,994.397662 | 0.0 | 5.68e-14 | 605 / 2 | PASS |
| 012 | 607 | 3.652 | 38.350 | 39.414 | **82.280** | -3,223,328.612327 | 0.0 | 5.68e-14 | 605 / 2 | PASS |
| 013 | 607 | 3.651 | 38.150 | 39.584 | **82.229** | -3,054,141.603864 | 0.0 | 5.68e-14 | 605 / 2 | PASS |
| 022 | 607 | 3.651 | 38.419 | 39.727 | **82.650** | -756,951.998440 | 0.0 | 5.68e-14 | 605 / 2 | PASS |
| 043 | 607 | 3.552 | 49.498 | 39.877 | **93.789** | 3,212,002.960929 | 0.0 | 1.14e-13 | 603 / 4 | PASS |
| 052 | 607 | 3.504 | 50.065 | 39.366 | **93.792** | 5,124,197.229243 | 0.0 | 1.14e-13 | 603 / 4 | PASS |
| 166 | 603 | 101.966 | 168.136 | 15.823 | **286.821** | 2,284,068.929973 | 0.0 | 1.14e-13 | 572 / 31 | PASS |

Scenarios 003 through 052 each produced 608 solution files and 608 evaluator
details. Scenario 166 produced 604 of each. Every evaluation certificate
reported the expected and observed counts as equal, a complete label set, and
an empty infeasible-label list.

## Scenario 166 evidence

The two Phase-I cases that had triggered expensive rebuilds completed their
projected-balance attempt optimally in 22 IPM iterations:

| Contingency | Corrective wall (s) | Internal Phase-I limit (s) | Status |
|---|---:|---:|---|
| `CTG_000017` | 73.729 | 90 | Optimal and independently feasible |
| `CTG_000323` | 74.665 | 90 | Optimal and independently feasible |

The final run screened all 603 contingencies, accepted 572 fast solutions, and
sent exactly 31 to the linearized corrective path. Once screening completed,
the corrective pool expanded from 8 to 12 workers. The two cases left
unfinished by the preceding run, `CTG_000452` and `CTG_000486`, both completed
in the retained run.

The official MPI certificate recorded 604 expected and observed details,
objective `2,284,068.9299734263`, infeasibility `0.0`, and no infeasible label.
The maximum independent contingency residual was `1.1368683772161603e-13`,
well below `1e-5`.

## Local evidence

Large run artifacts remain ignored and local under:

- `runs/validation_eab4f29_20260824/` for scenarios 003–052;
- `runs/postscreen_eab4f29_20260824/` for scenario 166.

They were not added to Git. Code, tests, and this report remain outside
OneDrive in `C:\Users\thoma\Documents\gravityx-go2-cpp`.
