# C2FEN19402 frozen-score economic diagnosis — 2026-08-27

## Scope and evidence boundary

This is a read-only audit of the five retained runs from frozen algorithm
revision `c477a2d4673a87d4d569110471956e420bb919e6`.  No scenario was rerun and
no frozen artifact was changed.  `scripts/analyze_retained_objective.py`
independently read every retained official-evaluator detail, required one
unique detail for every source contingency, required all duplicated BASECASE
details to be byte-identical, and reproduced the official score as

```
BASECASE objective + arithmetic mean(contingency objectives).
```

The published MSpp and best scores below are aggregate workbook results only.
Matching per-component published solution artifacts are not available in the
repository, so this report does not invent a decomposition for those scores.

## Finding

The low scores are not primarily a unit-commitment-cost problem.  Every source
generator was online in all five source cases, and the frozen path preserved
that status exactly: there were zero commitment changes, zero startups, and
zero shutdowns.  Generator on-cost was only `$1,535–$1,613` in the aggregate
score.

The dominant loss is the official evaluator's permitted active/reactive bus
imbalance penalty.  It accounts for essentially the entire gap to MSpp.  The
frozen independent residual checker validated nodal balance *after subtracting
the modeled `p_delta`/`q_delta` variables*.  Therefore its near-zero residual
did not mean raw physical imbalance was near zero; it meant the remaining
imbalance was represented by an allowed, heavily penalized variable.

Four scenarios suffered a second, causal implementation problem.  The fast
path treated each branch-flow component as if it had the hard box bound
`|pf|, |qf|, |pt|, |qt| <= rating`.  The official AC formulation instead uses
the apparent-current constraint with a permitted soft margin,
`p^2 + q^2 <= (rating * (vm + sm_slack))^2` for lines (and the corresponding
transformer expression).  The source-voltage repairs for 006, 010, 069, and
077 already had zero modeled balance residual and negligible official thermal
residual, but were rejected solely by the extra component box with violations
of `0.0440–0.1184` p.u.  The subsequent Phase-I feasibility path then drove
most flexible load toward its lower bound and retained the first officially
feasible, not economically optimized, state.

Scenario 095 did not encounter that box-bound rejection.  It retained 99.90%
of prior load, but the fast Newton repair deliberately preserved permitted bus
imbalance and never ran the actual market-surplus optimization.  Its bus
penalty was still `$36.70 million`.

## Exact objective reconstruction

Amounts below are the aggregate contribution `base + contingency mean`.
Benefits enter the objective positively; all listed costs enter negatively.

| Scenario | Base objective | Contingency mean | Reproduced objective | Load benefit | Bus penalty | Generator cost | Line cost | Transformer cost |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 006 | -339,671,825.332 | -339,679,987.482 | -679,351,812.814 | 187,025.678 | 679,245,437.044 | 292,094.357 | 1,307.090 | 0.000 |
| 010 | -291,010,062.208 | -291,019,733.077 | -582,029,795.285 | 207,884.393 | 582,008,054.121 | 228,955.086 | 670.471 | 0.000 |
| 069 | -297,400,730.869 | -297,413,075.833 | -594,813,806.702 | 313,020.463 | 594,871,644.997 | 254,322.053 | 860.115 | 0.000 |
| 077 | -216,482,502.196 | -216,533,823.514 | -433,016,325.710 | 233,787.121 | 433,021,524.498 | 228,530.903 | 57.430 | 0.000 |
| 095 | -18,042,485.138 | -18,252,782.816 | -36,295,267.954 | 718,408.654 | 36,698,904.811 | 314,641.286 | 130.512 | 0.000 |

The bus penalty splits as follows:

| Scenario | Active-balance penalty | Reactive-balance penalty | Total bus penalty |
|---:|---:|---:|---:|
| 006 | 371,526,689.802 | 307,718,747.243 | 679,245,437.044 |
| 010 | 321,342,839.841 | 260,665,214.280 | 582,008,054.121 |
| 069 | 302,655,548.943 | 292,216,096.054 | 594,871,644.997 |
| 077 | 240,067,154.433 | 192,954,370.065 | 433,021,524.498 |
| 095 | 33,261,709.923 | 3,437,194.887 | 36,698,904.811 |

Generator cost further decomposes into energy and on-cost.  Startup and
shutdown costs are exactly zero in every base and contingency detail; line
cost is entirely limit-violation cost; line/transformer switching and all
transformer costs are exactly zero.  The exhaustive `total_*` scan found no
other nonzero official objective component.

## Ranked gap to the prior-point benchmark

Because MSpp component artifacts are unavailable, only the aggregate gap can
be measured.  The bus-penalty column is a confirmed frozen component; its
share of the gap demonstrates its dominance but is not represented as an MSpp
component decomposition.

| Scenario | MSpp | MSpp minus frozen | Frozen bus penalty | Bus penalty / gap | Frozen net non-bus contribution |
|---:|---:|---:|---:|---:|---:|
| 006 | 152,989.330 | 679,504,802.144 | 679,245,437.044 | 99.962% | -106,375.770 |
| 010 | 275,970.553 | 582,305,765.838 | 582,008,054.121 | 99.949% | -21,741.164 |
| 069 | 712,738.623 | 595,526,545.325 | 594,871,644.997 | 99.890% | 57,838.295 |
| 077 | 155,651.051 | 433,171,976.761 | 433,021,524.498 | 99.965% | 5,198.788 |
| 095 | 203,325.479 | 36,498,593.433 | 36,698,904.811 | 100.549% | 403,636.857 |

## Dispatch, load, commitment, and controls

| Scenario | Committed / source units | Commitment changes | Load served (MW) | Prior load (MW) | Load retained | Generation (MW) | Prior generation (MW) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 006 | 968 / 968 | 0 | 33,015.736 | 155,304.998 | 21.259% | 117,926.166 | 158,531.577 |
| 010 | 968 / 968 | 0 | 34,405.789 | 150,704.913 | 22.830% | 115,273.920 | 154,024.722 |
| 069 | 927 / 927 | 0 | 28,126.067 | 131,329.032 | 21.416% | 99,013.089 | 133,919.196 |
| 077 | 921 / 921 | 0 | 37,784.507 | 131,706.774 | 28.688% | 102,979.567 | 134,212.502 |
| 095 | 921 / 921 | 0 | 133,456.970 | 133,590.402 | 99.900% | 131,058.342 | 136,074.380 |

Confirmed control observations:

- Division 1 line and transformer switching remained zero.
- Base transformer phase-shift movement was zero.  The evaluator reported a
  maximum base tap displacement of `0.095` in every scenario; contingency tap
  movement from the selected base was zero.
- Base switched-shunt displacement was numerical zero.  A small fraction of
  contingencies changed shunts; the maximum evaluator-reported displacement
  was `2.0`.
- Base real-load movement reached the source ramp allowance (`0.4166666665`
  p.u.) in every scenario.  Generator real and reactive movements were also
  material, but their direct cost was tiny relative to the bus penalty.
- `base_optimization_performed` is `false` for all five retained runs.

## Corrective design implication

The smallest evidence-supported change is not commitment optimization.  It is:

1. remove only the implementation-only branch-component box restriction while
   retaining the official apparent-current limit, soft-margin bound, penalty,
   and complete validator;
2. preserve the fast path's verified state as an incumbent;
3. run an actual fixed-commitment market-surplus refinement from that state;
4. accept the refinement only after the unchanged complete nonlinear and
   official gates pass; otherwise serialize the incumbent;
5. let corrective solves inherit the economically repaired base, then add a
   short Phase-II economic pass only where its measured cost fits the global
   deadline.

Fast-start commitment optimization is not justified at this stage because all
source units are online and commitment-related costs are less than 0.001% of
the dominant imbalance penalty.
