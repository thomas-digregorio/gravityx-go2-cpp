# C2FEN19402 economic-improvement pilot — `6f0afaa`

Date: 2026-08-27

## Outcome

Scenario 095 passed every requested gate in its sole cold run. Its official
objective improved from `-36,295,267.953619` to `210,762.417006907`, exceeded
the MSpp prior-point benchmark by `7,436.937870907`, covered all 6,579 source
contingencies, had official infeasibility `0.0`, passed the independent
`1e-5` residual gate, and completed in `239.470985400` seconds.

The next authorized cold run, scenario 006, exceeded the 300-second limit. The
controller stopped after `301.800822300` seconds with only 3,196 of 6,693
contingencies completed. Therefore it produced no complete official objective
or security certificate. Per the registered stop policy, scenarios 010, 069,
and 077 were not started. No retries or broader regression were performed.

## Frozen identities and tests

- Preserved baseline revision:
  `c477a2d4673a87d4d569110471956e420bb919e6`.
- Tested economic-improvement revision:
  `6f0afaaf7b238d67373fa0d8949c9648f9ac1faf`.
- Branch: `codex/c2fen19402-economic-improvement`.
- The retained baseline runs, report, tag, and evidence were not changed or
  rerun.
- Both standard and native C++ component suites passed 4/4, and the Python
  suite passed 49/49 before the full runs.
- The implementation continued to use the open-source CPU stack: Gravity C++,
  HiGHS, and Ipopt/MUMPS. It did not add GPU or commercial-solver work, omit
  source contingencies, change a source case, relax a tolerance, or enable
  Division 1-prohibited line/transformer switching.

The official workbook comparison comes from `Division 1 scores!D106:AQ110` in
`C:\Users\thoma\Documents\goc2-ac-score-check\final.xlsx`. Published best
values are comparison points, not proven global optima. Competitor
component-level solution artifacts were unavailable, so no competitor
component decomposition is inferred.

## Requested five-scenario comparison

Percentage improvement is `(new - frozen) / abs(frozen)`. The final column is
the requested `(new - MSpp) / (published best - MSpp)`.

| Scenario | Outcome | Frozen objective | New objective | MSpp | Best published objective (team) | Absolute improvement | Improvement | Fraction of best gain captured |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 095 | **PASS** | -36,295,267.953619 | 210,762.417007 | 203,325.479136 | 444,334.416768 (Pearl Street Technologies) | 36,506,030.370626 | 100.5807% | 3.0858% |
| 006 | **STOP — deadline** | -679,351,812.814108 | Not produced | 152,989.330341 | 459,401.669038 (Pearl Street Technologies) | N/A | N/A | N/A |
| 010 | Not run after 006 stop | -582,029,795.284929 | N/A | 275,970.553456 | 549,955.199224 (Pearl Street Technologies) | N/A | N/A | N/A |
| 069 | Not run after 006 stop | -594,813,806.702179 | N/A | 712,738.622667 | 953,207.754053 (GravityX) | N/A | N/A | N/A |
| 077 | Not run after 006 stop | -433,016,325.710104 | N/A | 155,651.050929 | 395,220.476943 (GravityX) | N/A | N/A | N/A |

Scenario 095 remains `233,571.999761` below the published best. Exceeding MSpp
is a useful first economic target, but capturing 3.0858% of the
MSpp-to-published-best interval also shows that substantial economic headroom
remains.

## Scenario 095 score decomposition

The independent detail audit required exactly one unique record for every
source contingency, required all 47 duplicated base records to be
byte-identical, and reproduced the score exactly as

`199,628.070728807 + 11,134.346278100 = 210,762.417006907`.

| Objective contribution | Base | Contingency average | Base + contingency average |
|---|---:|---:|---:|
| Load benefit | 348,005.707419 | 348,002.499602 | 696,008.207021 |
| Generator energy cost | 144,074.902239 | 144,075.119795 | 288,150.022034 |
| Generator on/no-load cost | 767.500000 | 767.441607 | 1,534.941607 |
| Startup cost | 0.000000 | 0.000000 | 0.000000 |
| Shutdown cost | 0.000000 | 0.000000 | 0.000000 |
| Active-balance penalty | 1,826.488757 | 87,814.782614 | 89,641.271371 |
| Reactive-balance penalty | 380.987621 | 103,277.782305 | 103,658.769926 |
| Total bus penalty | 2,207.476378 | 191,092.564919 | 193,300.041297 |
| Line-limit cost | 1,327.758074 | 933.027003 | 2,260.785077 |
| Transformer cost | 0.000000 | 0.000000 | 0.000000 |
| Line/transformer switching cost | 0.000000 | 0.000000 | 0.000000 |
| **Net objective** | **199,628.070729** | **11,134.346278** | **210,762.417007** |

The objective uses load benefit positively and subtracts the listed costs and
penalties.

## Why the 095 score changed

The frozen score collapse was not caused by commitment costs. All 921 source
generators were online, and the new solution retained all 921 with zero
commitment changes, startups, or shutdowns.

The dominant change was the permitted but heavily penalized bus-imbalance
term. Its aggregate contribution fell from approximately `$36.699 million` in
the frozen result to `$0.193 million`. That approximately `$36.506 million`
penalty reduction dominates the score gain. It more than offset approximately
`$22.400 thousand` less load benefit and `$2.130 thousand` more line-limit
cost; generator cost also fell by approximately `$24.956 thousand`.

The implementation made the smallest evidence-supported correction:

1. It removed the implementation-only separate branch-flow-component boxes
   while retaining the official apparent-current constraints and soft-margin
   semantics.
2. It always preserved the last independently verified feasible state as the
   fallback incumbent.
3. It applied bounded fixed-commitment market-surplus refinements to continuous
   dispatch and permitted AC controls.
4. A candidate could replace the incumbent only after the unchanged exact
   nonlinear security check passed.
5. Corrective economic work was explicitly budgeted so exhaustive screening,
   official evaluation, independent verification, and serialization could
   finish inside the global limit.

The selected base state served `123,202.547140286 MW` versus
`133,590.401993550 MW` in the source point and dispatched
`125,100.224660100 MW` versus `136,074.379946700 MW` in the source point.
All 921 generators changed real dispatch, 914 changed reactive dispatch, and
12,267 loads changed quantity. These are continuous economic/control changes,
not commitment changes.

## Security and timing evidence

| Scenario | Completed / source contingencies | Screened | Official infeasibility | Max independent residual | Fallback count | Base process | Contingency stage | Evaluation work | End-to-end |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 095 | 6,579 / 6,579 | 6,579 | 0.0 | 9.469268751949e-6 | 0 | 55.756 s | 175.354 s | 144.330 s overlapped; 5.889 s tail | **239.471 s** |
| 006 | 3,196 / 6,693 | 3,192 | Not reached | Not complete | 1 before stop | 56.306 s | 243.987 s before stop | Not complete | **301.801 s — FAIL** |
| 010 | 0 / 6,693 | 0 | Not run | Not run | N/A | N/A | N/A | N/A | N/A |
| 069 | 0 / 6,620 | 0 | Not run | Not run | N/A | N/A | N/A | N/A | N/A |
| 077 | 0 / 6,584 | 0 | Not run | Not run | N/A | N/A | N/A | N/A | N/A |

For 095, the official certificate records 6,580 unique labels (base plus all
6,579 contingencies), a complete label set, zero infeasible labels, and exact
agreement between the reported and detail-derived objective. The independent
maximum residual is below, but close to, the requested `1e-5` gate.

For 006, the base solve itself succeeded in `51.058406434` algorithm seconds
and `56.305873200` process seconds. The exact stop reason was:

`fast screen CTG_004927 finished after the end-to-end work deadline`

The last completed label recorded was `CTG_005572`. Because the full label set
was not completed, no partial objective is reported as a solution and the run
is neither certified feasible nor an infeasibility result.

## Stop decision

The registered policy required an immediate stop after the first deadline or
security failure. Consequently, 010, 069, and 077 remain deliberately unrun,
and the complete six-family regression remains unstarted pending review and
explicit approval.
