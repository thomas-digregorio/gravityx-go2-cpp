# Validation gates

## Independent formulation checks

The C++ checker recomputes, without consulting Gravity's expression graph:

- voltage, generator, load, flow, slack, and commitment bounds;
- exact piecewise-linear convex-combination sums and power recovery;
- source PMIN/PMAX and base/contingency ramp windows;
- AC active and reactive nodal balances, including shunts;
- line and transformer polar AC branch equations;
- applicable angle limits and soft thermal/current limits; and
- the exact generator or branch removed in each contingency.

Base and commitment-rounding results require a successful Gravity/Ipopt status
and an independent maximum residual no greater than `1e-5`.  A corrective
contingency point with a finite objective is retained whenever the independent
checker verifies the same `1e-5` residual limit, even if Ipopt terminates after
finding that feasible point without reporting convergence.  The raw Ipopt
status and this acceptance path are recorded, and the official evaluator must
still certify every submitted case.

The component suite also solves a two-bus AC case with two electrically
identical parallel circuits.  This guards the branch-specific symbolic tags
needed by the pinned Gravity revision to prevent collisions between otherwise
identical nonlinear DAG nodes.  Every tag is a fixed parameter with value one,
so it changes symbolic identity but not the equation value or derivative.

## Tiny semantic oracle

The official `PowerModelsSecurityConstrained.jl` formulation was used only as
a development oracle.  With the same selected base state, the following
objectives were observed:

| Model | Gravity C++ | Official Julia formulation | Absolute difference |
|---|---:|---:|---:|
| Base soft ACOPF | 593064.468772815 | 593064.4687728151 | < 1e-9 |
| Continuous AC-UC relaxation | 593064.1922037087 | 593064.192203709 | < 1e-9 |
| Generator outage `G_1_1` | -118718.210379945 | -118718.210379947 | 2.1e-9 |
| Branch outage `L_1_2_1` | 585913.891034013 | 585913.891034009 | 3.4e-9 |

The displayed voltage, angle, active-generation, and reactive-generation
vectors also matched.  The complete 14-bus fixture, including all nine source
contingencies, then passed the official evaluator with:

- official objective: `1109587.974478508`;
- official infeasibility: `0.0`;
- maximum independent C++ contingency residual: `5.06e-8`; and
- zero evaluated line and transformer limit violation.

These checks establish formulation equivalence for the exercised components;
they do not make the reproduction identical to the unpublished winning
GravityX implementation.

## Cold-run discipline

`scripts/run_experiment.py` refuses a nonempty output directory.  The base
solve starts a new process.  Corrective contingencies are assigned
deterministically, round-robin, to isolated persistent worker processes; each
worker loads the immutable case and selected base state once, then constructs
and solves a fresh Gravity/Ipopt model for each assigned contingency.  Workers
run concurrently and are limited to one BLAS thread each to avoid
oversubscription.  No contingency solution or optimizer state is reused by a
different contingency.

For Final-event Division 1 experiments, the runner enforces the historical
timing stages separately: Code1 has a 300-second wall limit and Code2 receives
two seconds per source contingency.  A failed or timed-out run leaves a
machine-readable `run_status.json`; a run is successful only after the official
evaluator confirms every submitted case is feasible.  The five-scenario
reproduction additionally enforces a stricter 300-second end-to-end limit from
normalized-case loading through official evaluation and result serialization,
with seven seconds reserved for evaluation and one second for finalization.
