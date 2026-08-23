# Method and fidelity scope

## Public GravityX method implemented

1. Construct the nonlinear AC unit-commitment model in Gravity C++.
2. Solve its continuous relaxation with Ipopt.
3. Order eligible commitment variables by distance from the rounding threshold.
4. Round and fix a deterministic batch, then re-solve the reduced NLP.
5. Continue until all eligible commitments are integral.
6. Solve a fixed-commitment AC repair and retain only a valid, better candidate.
7. Produce all source contingency solutions and score with the official
   evaluator.

The corrective stage uses the selected base commitment without additional
switching.  Each source contingency receives its own AC state and corrective
generator/load movement within `deltarctg` and the source `prumaxctg` /
`prdmaxctg` limits.  Soft active/reactive balance and branch-limit penalties
match the official GO2 formulation.

## Evidence boundary

The public solver report describes iterative batch rounding, custom batch
ordering, custom rounding, and AC-loss-aware unit commitment.  It does not
publish the winning submission's source, exact batch rule, or dataset-specific
parameter values.  This project therefore distinguishes:

- **framework fidelity**: actual Gravity C++ is used to formulate and solve the
  nonlinear model;
- **method fidelity**: the publicly described relaxation/fix/re-solve loop is
  implemented;
- **unavailable details**: proprietary ordering, rounding, and tuning are
  replaced by deterministic documented choices.

No run from this repository should be described as the original GravityX
binary or an exact reproduction of unpublished code.

The reproduction also does not claim that its deterministic ordering (distance
from 0.5, then source index), four-batch partition, or 0.5 tie rule matches the
winning team's private choices.  Those are frozen, auditable substitutes for
details that were never released.
