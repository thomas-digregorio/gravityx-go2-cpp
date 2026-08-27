# C2FEN08300 final-scenario results — 2026-08-24

## Scope

This report records all seven public GO Challenge 2 Final Event 4
`C2FEN08300` scenarios: `003`, `012`, `013`, `022`, `043`, `052`, and `166`.
The source files define 607 contingencies in scenarios 003 through 052: 41
generator outages and 566 branch outages. Scenario 166 defines 603
contingencies: 37 generator outages and the same 566 branch outages. No source
contingency was omitted.

Each retained run was cold. It loaded its own normalized scenario and source
prior state without carrying an operating point from another scenario. The
hard 300-second boundary starts at normalized-case loading and includes base
construction, every source contingency, independent validation, complete
official evaluation, result serialization, and the final timing check.

Scenarios 003 through 052 used frozen revision
`bec6b3d7805290697fb490902da88d09321c92bb`. Scenario 166 exposed a slower
base and corrective tail. Its retained replacement used frozen revision
`c5832293b1c7fee6ec9349c77c06575f566a523d`, which adds an independently
validated lightweight second base linearization, a measured 31-case fallback
schedule, and strict parallel-evaluator bookkeeping validation. The later
paths are opt-in; the earlier retained run artifacts were not modified.

## Method

The scalable C++ path used the source commitment and the following feasibility
pipeline:

1. construct a base operating point with HiGHS sequential-linearized AC
   subproblems and accept it only after complete nonlinear validation;
2. screen every source contingency with parallel sparse-Newton AC power flow;
3. send every screen failure to an isolated HiGHS sequential-linearized
   corrective solve;
4. accept no base or corrective state unless the independent nonlinear checker
   passes all source bounds, AC equations, ramps, and applicable limits;
5. write one solution for the base case and every contingency and evaluate the
   complete set with the pinned official GO Challenge 2 evaluator.

The first six scenarios used eight fast-screen and eight corrective workers.
Scenario 166 used sixteen fast-screen workers, eight corrective workers, and a
balanced scheduling-only profile for its 31 measured fallbacks. The profile
contains labels, worker assignments, and prior wall times only—no primal,
dual, commitment, or network state. All HiGHS instances used one thread to
avoid sparse-factorization oversubscription.

The first six official evaluations were serial. Scenario 166 used eight
Microsoft MPI evaluator processes. The unchanged vendor MPI evaluator writes
correct per-case details and recomputes the top-level objective and
infeasibility, but leaves four serial accumulator fields stale. The runner
therefore required the exact 604-label set, independently checked every detail
and aggregate, archived the raw vendor summary, and reconstructed only those
four bookkeeping fields. Any missing, extra, malformed, non-finite, or
infeasible detail would have failed the run.

The internal contingency checker applied the stricter normal `RATE_A` branch
ratings, while the official evaluator applied the source corrective `RATE_C`
semantics.
The reported PASS therefore requires both the conservative internal screen and
the unchanged official evaluation; no source limit was weakened.

## Results

| Scenario | Contingencies | Base (s) | Code 2 (s) | Evaluation (s) | End-to-end (s) | Official objective | Official infeasibility | Max independent residual | Fast / linearized | Revision | Result |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| 003 | 607 | 1.271 | 213.389 | 37.716 | **253.021** | 7,583,789.940537 | 0.0 | 1.71e-13 | 574 / 33 | `bec6b3d` | PASS |
| 012 | 607 | 4.658 | 214.988 | 39.896 | **260.215** | -3,291,776.055542 | 0.0 | 1.71e-13 | 574 / 33 | `bec6b3d` | PASS |
| 013 | 607 | 4.607 | 226.861 | 38.272 | **270.410** | -3,122,743.496285 | 0.0 | 1.71e-13 | 574 / 33 | `bec6b3d` | PASS |
| 022 | 607 | 4.331 | 212.303 | 38.608 | **255.892** | -826,185.200476 | 0.0 | 1.71e-13 | 574 / 33 | `bec6b3d` | PASS |
| 043 | 607 | 16.631 | 224.179 | 38.951 | **280.457** | 2,980,214.020276 | 0.0 | 1.71e-13 | 573 / 34 | `bec6b3d` | PASS |
| 052 | 607 | 1.238 | 221.854 | 37.298 | **261.047** | 4,891,889.007695 | 0.0 | 1.71e-13 | 573 / 34 | `bec6b3d` | PASS |
| 166 | 603 | 98.064 | 168.504 | 11.795 | **279.228** | 2,284,068.929973 | 0.0 | 1.14e-13 | 572 / 31 | `c583229` | PASS |

`Fast / linearized` partitions every source contingency by the method whose
state passed independent validation. Stage times do not sum exactly to the
end-to-end value because loading, orchestration, file generation, and final
serialization are outside those stage counters. The table uses the final
process-return timing check, which is slightly more conservative than the
timing checkpoint persisted immediately before the final status rewrite.

All seven scenarios completed below 300 seconds. The slowest was scenario 043
at 280.457 seconds, leaving 19.543 seconds. Scenario 166 completed at 279.228
seconds, leaving 20.772 seconds.

## Scenario 166 correction

The best earlier cold scenario-166 attempt completed 600 of 603 contingencies
in 294.230 seconds and never reached official evaluation. Its robust base took
116.885 process seconds, and its measured fallback queues left three cases
unfinished.

The retained correction preserves all AC nodal-balance equations in its second
base linearization, fixes explicit balance slacks to zero, and places voltage
and angle variables in a local trust box. A branch-flow or angle row is omitted
only when interval evaluation proves that row redundant throughout the box.
The resulting candidate still must pass the complete nonlinear checker. On the
retained run, the base took 98.064 process seconds and its validation maximum
was `5.68e-14`.

The corrected base reduced the contingency partition to 572 fast acceptances
and 31 linearized repairs. The balanced fallback queues completed all 603
contingencies in 168.504 seconds. Eight-process official evaluation then took
11.795 seconds. No Ipopt fallback and no nonconverged acceptance was used.

## Acceptance evidence

The retained evidence proves, for every scenario:

- `run_status.json` reports `success=true`, stage `complete`, and the
  end-to-end gate satisfied;
- the recorded contingency count equals the source count;
- exactly one solution and one official evaluator detail exist for the base
  case and every source contingency;
- the objective and infeasibility maps contain that same complete label set;
- `eval_summary.json` reports aggregate `infeas=0.0`;
- every independent maximum residual is below `1e-5`;
- fast and linearized accepted counts exhaustively partition the source
  contingencies, with zero nonconverged acceptances.

For scenario 166, the parallel-evaluation certificate additionally reports 604
expected and observed details, a complete label set, zero infeasible labels,
and archived raw vendor JSON and CSV summaries.

The source and normalized-case hashes are:

| Scenario | `case.raw` SHA-256 | `case.con` SHA-256 | `case.json` SHA-256 | Normalized-case SHA-256 |
|---:|---|---|---|---|
| 003 | `d9111af99d32a626977abbb26e03d65555af4f0912e698ebc59897d16e3f9cb3` | `43bc83b5bd22a75e5cfbe4e4a946db03d95daf0a848b86bd32f7353caa67fead` | `08e2756d5555fd70bb880309fecd0b96bafdd0960c6bb15bfa3ec0095ce8adaa` | `17e819cac41c8d0f6238724adc91083cb203de8ec77bf92f4b840b309e1dbde3` |
| 012 | `d9111af99d32a626977abbb26e03d65555af4f0912e698ebc59897d16e3f9cb3` | `43bc83b5bd22a75e5cfbe4e4a946db03d95daf0a848b86bd32f7353caa67fead` | `06240bd997951b5d38dafc19cccb6c812ec7f6e8dbe626a0fbce90f27b936338` | `e2e5ad78165018aee20c271e68a184e5755d7c4544f9ccf511ce23b747809fe8` |
| 013 | `d9111af99d32a626977abbb26e03d65555af4f0912e698ebc59897d16e3f9cb3` | `43bc83b5bd22a75e5cfbe4e4a946db03d95daf0a848b86bd32f7353caa67fead` | `850499edd5b441aad28c2d8acb8538764b57d6e7d9afeb417ffe8009460b7d1d` | `bb2fcd0a66144030330d5bcf4235424eab6efcac7f475a2b961491039215dade` |
| 022 | `d9111af99d32a626977abbb26e03d65555af4f0912e698ebc59897d16e3f9cb3` | `43bc83b5bd22a75e5cfbe4e4a946db03d95daf0a848b86bd32f7353caa67fead` | `50a9553fc319d58f819a03cd6b5600a3f902f70b741dc655c4de965fcfdc7ebb` | `c0cb193252d5b485b785c0f5792f8eabdb19f0eb494caf93d0bfed44d2a1aa7b` |
| 043 | `128cbda24e42572765de3103377df181df8436eb8b4aed31d15179e4c81fae24` | `43bc83b5bd22a75e5cfbe4e4a946db03d95daf0a848b86bd32f7353caa67fead` | `c20d7e2d10832c2aeb5c9afda565f07b98fd1516aa46fb5df34b75c74fb10789` | `ab45fa62a7febab039a89a933cd02dc843cecc0d316ff1b2206c0f86bdc53f6e` |
| 052 | `128cbda24e42572765de3103377df181df8436eb8b4aed31d15179e4c81fae24` | `43bc83b5bd22a75e5cfbe4e4a946db03d95daf0a848b86bd32f7353caa67fead` | `2cb54c69c876ad6d04520d45c98692b838712b8573e9f58444022ed12546529b` | `36019208299491f836a5a4bdca03f3ac4567955652b587ec70658ecc59a2db10` |
| 166 | `cb40922ea24c74ab0524c42676e51817fecc4e6e5135e6e09e4343a172f67b13` | `149901274d9b61077d7a94041f241539f3bc2a8f4890d7ce70ea3b6a206dfd5a` | `99f06e343a14113a3284771c8b5b21fec691ccf048b6b504aea7e9328ab9ac6b` | `6ead1a631345205dfdc7d38a262ab41c0de639cbf466cb7fb8f7ba3fafefaef0` |

Raw competition inputs and run outputs remain ignored. The local retained
directories contain the complete statuses, timing summaries, worker logs,
independent validation records, official evaluator outputs, and all base and
contingency solutions. Scenario 003 is retained in the local WSL archive; the
other six are retained under the repository's ignored `runs` directory.
