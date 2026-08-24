# C2FEN02020 final-scenario results — 2026-08-24

## Scope

This report records the five public GO Challenge 2 Final Event 4
`C2FEN02020` scenarios. Each reported run was cold, used only open-source
Gravity C++ and Ipopt/MUMPS, enforced the complete source contingency list,
and had a 300-second end-to-end deadline from normalized-case loading through
official evaluation and result serialization.

The base commitment is the source prior commitment. Corrective states retain
the source generator, load, branch, transformer, ramp, voltage, and balance
semantics. Every returned state is checked independently before it is written,
and the complete solution set is then checked by the pinned official GO
Challenge 2 evaluator. No contingency is omitted and no feasibility tolerance
is relaxed.

## Results

| Scenario | Source contingencies | Workers | Base (s) | Code 2 (s) | Evaluation (s) | End-to-end (s) | Official objective | Official infeasibility | Max independent residual | Fast / exact | Revision | Result |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| 025 | 292 | 2 | 70.021 | 135.549 | 11.018 | **216.767** | 3,537,475.452987 | 0.0 | 8.64e-6 | 287 / 5 | `657958d` | PASS |
| 121 | 300 | 3 | 73.955 | 146.953 | 11.425 | **232.523** | 5,297,595.453676 | 0.0 | 8.48e-6 | 294 / 6 | `7d2eecd` | PASS |
| 134 | 300 | 3 | 73.248 | 139.497 | 13.195 | **226.131** | 990,121.476138 | 0.0 | 7.90e-6 | 296 / 4 | `7d2eecd` | PASS |
| 260 | 300 | 3 | 95.668 | 162.742 | 11.494 | **270.092** | 3,783,725.438590 | 0.0 | 9.90e-6 | 293 / 7 | `7d2eecd` | PASS |
| 262 | 300 | 3 | 76.264 | 142.503 | 11.969 | **230.922** | 4,339,127.157325 | 0.0 | 9.84e-6 | 296 / 4 | `7d2eecd` | PASS |

`Fast / exact` is the number of contingencies accepted by the independent
sparse-Newton checker versus the number sent to the exact resident Ipopt
corrective model. None of the exact fallbacks was accepted with a failed Ipopt
termination status in these five retained runs.

Scenario 025 is the already-retained cold PASS from the immediately preceding
validated fast-power-flow revision. It was not repeated after scenario 121
motivated the distributed-balance improvement. Scenarios 121, 134, 260, and
262 all use the same frozen and pushed revision, `7d2eecd`.

## Scenario-121 correction

The first scenario-121 run exhausted the end-to-end work budget after 73 of
300 contingencies. Its original fast screen accepted 282 cases and sent 18 to
Ipopt. Most rejected fast states were otherwise physical AC solutions whose
loss mismatch had been concentrated at one reference bus beyond the source
model's 0.5-p.u. per-bus balance-slack bound.

Revision `7d2eecd` adds a distributed-slack Newton correction. For each
connected component it:

1. constructs active-balance targets no larger than 0.49 p.u. per bus;
2. solves the AC equations and one scalar nonlinear-loss correction together;
3. enforces voltage and reactive-generator bounds;
4. retains the original state unless full independent validation improves;
5. sends any still-invalid state to the unchanged exact Ipopt fallback.

On the saved scenario-121 base point, exhaustive screening improved from
282/300 to 294/300 fast-valid states. The cold replacement then passed in
232.523 seconds. This is an algorithmic use of balance variables already in
the official corrective model; it does not widen their 0.5-p.u. bounds.

## Acceptance evidence

For every retained run:

- every source contingency produced a solution file;
- independent validation was at or below `1e-5`;
- official `infeas` was exactly `0.0`;
- official per-case infeasibility flags were false;
- the complete end-to-end run finished below 300 seconds.

The normalized-case SHA-256 hashes are, by scenario:

| Scenario | Normalized-case SHA-256 |
|---:|---|
| 025 | `b8e60e82747990fe5b8348619a790bafb78ec9f191143c0a83f8aed5d386fbf7` |
| 121 | `b842e35c61822017bf958b18c4a178144e3d96bae272c9c9af4218f873f18d23` |
| 134 | `d7ab0c338ee8846e61f9620b14d14fa5b5ebc36c7b153777a3bce8df9e665f25` |
| 260 | `c858ad365e95cd425f2e1a260687d648afc9316d59b89e31feff4f122302aed6` |
| 262 | `5aad89ebff838d9ef30da8cbf7b620bbf7f43a2f3bf9018825449ab49ae850bc` |

The ignored run directories contain the full manifests, worker assignments,
independent validation records, official evaluator detail, and all base and
contingency solution files.
