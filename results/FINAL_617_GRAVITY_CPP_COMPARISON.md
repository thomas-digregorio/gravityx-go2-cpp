# Final-617 Gravity C++ cold experiment

## Outcome

The framework-faithful Gravity C++ pipeline completed the GO Challenge 2 Final
case `C2FEN00617`, scenario `005`, and passed the official evaluator for the
base case and all 105 source contingencies.  The official objective was
`1325555.7041623043` with official infeasibility `0.0`.

This is a clean-room reproduction of the public GravityX method, not the
unpublished winning binary.  The nonlinear models are genuinely formulated in
the pinned Gravity C++ framework and solved by Ipopt/MUMPS.  The winning team's
private batch ordering, rounding rules, and tuning remain unavailable and are
replaced by the deterministic choices documented in `docs/METHOD_SCOPE.md`.

## Frozen experiment identity

| Item | Frozen value |
|---|---|
| Repository revision | `c6d22a87c888221750fbd21dcb1deb706944ede1` |
| Gravity revision | `f5af33ed9572829d96c6f54b1a3ad30f53677fe7` |
| Solver | Ipopt 3.14.19 with MUMPS |
| Case | Final `C2FEN00617`, scenario `005` |
| Network | 617 buses, 94 generators, 723 lines, 130 transformers |
| Source contingencies | 105 total: 51 generator, 50 line, 4 transformer |
| Corrective workers | 22 processes, one BLAS thread per process |
| Raw input SHA-256 | `bdadd320f4ed425e46ccc1c19242c08b079af447e478fcb095379770524b5adf` |
| JSON input SHA-256 | `77c12c6218ca6aa228c6e04a862bd406740ca1a8d7d1146f9c930a950380261b` |
| CON input SHA-256 | `1e9e4c7ff07561eeb0a77b28f938a3abdf1148896bf3f453ab5c99eb1c57cd83` |
| Normalized model SHA-256 | `5f0fa87430e8246ffaad8b8296a44e7108fb1b780b0d27cd81d0a927e9f4fa67` |

The run began in a new process with no retained optimizer state.  Each
corrective contingency also used a separate process from the selected base
state.  Generated run files remain local and ignored by Git.

## Verified result

| Metric | Result |
|---|---:|
| Official total objective | 1325555.7041623043 |
| Official base objective | 645902.4212236030 |
| Official aggregate contingency objective | 679653.2829387012 |
| Official infeasibility | 0.0 |
| Evaluated cases present and feasible | 106 / 106 |
| Committed generators | 57 / 94 |
| Base generator startups | 8 |
| Independent base maximum residual | 3.56449385563451e-8 |
| Independent exhaustive maximum contingency residual | 3.57185312438446e-8 |
| Evaluator maximum line violation | 0.0 |
| Evaluator maximum transformer violation | 0.0 |
| Evaluator maximum active-balance slack | 2.04695859915005e-9 p.u. |
| Evaluator maximum reactive-balance slack | 0.0331735922172332 p.u. |

The reactive-balance value is an allowed, priced soft slack in the source GO2
formulation.  It is therefore reflected in the objective rather than reported
as hard infeasibility.  The independent residuals above measure satisfaction
of the formulated equations, including those slack variables.

The eight startups were generators `487:1`, `488:1`, `495:1`, `507:1`,
`521:1`, `538:3`, `538:5`, and `562:2`.

## Timing

| Boundary | Wall time (s) |
|---|---:|
| Base process, including C++ startup/model construction | 90.484 |
| Base algorithm inside the C++ executable | 78.069 |
| Parallel corrective-contingency boundary | 198.784 |
| Sum of 105 individual contingency solver times | 2621.747 |
| Official evaluator | 4.891 |
| Complete end to end | 294.180 |

The sum of individual contingency times is accumulated work, not elapsed
time; 22 independent workers reduced it to a 198.784-second wall boundary.
The slowest individual solve was `CTG_000046` at 27.861 solver seconds.

The base algorithm used 13.128 seconds for the source ACOPF, 16.712 seconds
for the initial continuous AC-UC relaxation, 4.763, 5.522, 5.481, and 7.222
seconds for the four batch-fixing rounds, and 15.472 seconds for the final
fixed-commitment AC repair.  No rounding round fell back to its prior state.

## Comparison with the prior Julia/PowerModels reproduction

Both experiments used the same Final-617 scenario, source contingencies,
Ipopt/MUMPS solver family, selected commitment semantics, official solution
format, and evaluator.  The earlier implementation used Julia and
PowerModelsSecurityConstrained rather than Gravity C++.

| Metric | Earlier Julia/PowerModels | Gravity C++ | Difference |
|---|---:|---:|---:|
| Official objective | 1325547.8344663763 | 1325555.7041623043 | +7.8696959279 |
| Official base objective | 645898.5266603442 | 645902.4212236030 | +3.8945632587 |
| Aggregate contingency objective | 679649.3078060322 | 679653.2829387012 | +3.9751326690 |
| Official infeasibility | 0.0 | 0.0 | 0.0 |
| Source contingencies evaluated | 105 | 105 | 0 |
| Solver pipeline wall time | 74.244 s | 289.268 s | 3.90x slower |
| End to end including evaluator | 76.352 s | 294.180 s | 3.85x slower |

The Gravity C++ score is higher by about `0.000594%`.  It scored higher in the
base case and every contingency.  Almost all of the total improvement is the
`7.869543` reduction in soft bus-imbalance cost; total generation cost was
identical and load benefit changed by only `0.000153`.

This is not a different commitment or active dispatch:

| Base decision comparison | Result |
|---|---:|
| Commitment mismatches | 0 / 94 |
| Committed units in each result | 57 |
| Maximum active-generator dispatch difference | 1.39e-17 p.u. |
| Maximum reactive-generator dispatch difference | 8.64e-6 p.u. |
| Maximum voltage-magnitude difference | 2.05e-8 p.u. |
| Maximum voltage-angle difference | 6.87e-10 radians |
| Maximum served-load-factor difference | 2.21e-8 |

The framework-faithful implementation was slower on this laptop.  The main
structural costs are Gravity's symbolic nonlinear model construction and
derivative graph handling plus construction of 105 separate corrective models
in isolated processes.  This result should not be read as a benchmark of the
original GravityX competition executable, whose private model-generation,
ordering, reuse, and tuning were not published.

## Numerical implementation gate

The first attempted full run stopped before producing a solution when the
pinned Gravity revision merged identically rendered nonlinear expressions for
parallel circuits.  That left a derivative-graph child without a value and
raised `Func get_val out of range`.  The correction attaches a fixed,
unit-valued branch/equation-specific symbolic parameter to each branch-flow
and thermal expression.  It changes expression identity without changing the
equation value or derivative.  A two-bus case with identical parallel circuits
now exercises this condition in the component suite.

After that correction, the replacement experiment above was launched cold at
the frozen revision.  All component tests and the official evaluator passed;
the aborted pre-fix directory is not a benchmark result.
