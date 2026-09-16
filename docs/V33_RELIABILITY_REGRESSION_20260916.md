# V33 bounded fallback repair and frozen regression registration

V32 passed 19,402-bus scenarios 006 and 010 but missed one contingency in
069. See `V32_RELIABILITY_REGRESSION_20260916.md` for exact evidence. No V32
077/095 or smaller-network runs were performed.

## Algorithm change

Only after every existing fast path fails, try a bounded expanded local
feasibility repair before the large corrective fallback LP. Existing passing
paths do not enter this new route. Search neighborhoods are four and six
post-outage hops, capped at 2,048 control buses, with the existing one-hop
fixed-control boundary. Each neighborhood starts from the same within-run
base state, with the original immutable base retained as the ramp anchor.

At most four linearization rounds per neighborhood; each LP receives at
most two seconds of the shared 16-second construction-inclusive work budget.
The candidate trust radii are 0.15 radian and 0.05 p.u.; the existing
0.5-p.u. candidate imbalance bound is unchanged. These are search bounds,
not changes to source constraints or acceptance tolerances. All original
LP rows are audited after fixed-column substitution, and every proposed
state is rebuilt and checked against the complete nonlinear source model.
Stop on stagnation, an unsuccessful LP, a work cap, or certified feasibility.
No feasible candidate means the existing full-network fallback remains.
The external 300-second end-to-end deadline remains the hard boundary.

Diagnostics record entry, selection, elapsed time, neighborhood depth,
linearized audit residuals, and full nonlinear residuals. No cross-run
initialization, new dependency, GPU, commercial solver or source change.

## Component checks

All 108 Python tests and all four native CTest groups passed. New tiny tests
cover an initially infeasible voltage point repaired by the new LP route,
full independent nonlinear validation, unchanged initial/base states,
zero/nonfinite/negative budget handling, bounded attempt/control counts,
expanded LP row/column audits, and absence of extra work on an already
successful path. No full network was used as a component test.

Release executable SHA256:
`108eced2c6e295e56de5eb52d0b63670042ab905f8b4c1223ddaf2a5baa3af92`.
Configuration: `config/v33_reliability_regression_20260916.json`.

## Run order and acceptance

After commit/push, run 069 cold once to test the failed case, then cold
006/010/077/095 once each if 069 passes. All five must pass on this unchanged
revision before starting the smaller-network regression: all five 617-bus,
five 2,020-bus, seven 4,224-bus, seven 8,300-bus, and eight 16,789-bus
scenarios. The configured scenario lists match every local source scenario
directory in those six families.

Every run uses the source state, all source contingencies, independent
native residual <=1e-5, complete official evaluation with infeasibility zero,
and <=300-second end-to-end time. Report official objective quality
separately; PASS does not establish global optimality, zero penalized
imbalance, or a top-five competition placement. Do not substitute historical
results, rerun unchanged scenarios, or alter the algorithm during regression.

Results pending. Preserve compact hash-backed evidence for all attempts and
the protected complete solution references during any authorized cleanup.
