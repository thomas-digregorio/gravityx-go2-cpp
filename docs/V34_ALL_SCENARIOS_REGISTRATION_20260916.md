# V34 all-scenario regression registration

The user clarified the completion gate: all 37 scenarios across all six
families must PASS on one frozen algorithm. V33's 36/37 is not completion.
Start with one cold 16,789-family scenario 089 run, then test every other
scenario once at this same revision. That 089 run counts toward the 37;
do not repeat it just to populate a second table. A failed scenario requires
a new evidence-backed fix, component tests, freeze, and replacement campaign,
not a changed acceptance threshold or a success claim.

## Fixes

1. Signed balance-pair reconstruction previously cleared both columns
before checking the replacement against a 1e-12 bound comparison. An
otherwise valid LP point at a 0.5 bound with around 2e-11 reconstruction
roundoff could be rejected and reported as a 0.5 row violation. Reconstruction
is now atomic: both proposed values are checked first, differences within
the existing 1e-8 solver primal tolerance are clamped to the exact existing
column bounds, and the complete existing 1e-7 expanded-row audit remains.
Larger or nonfinite errors fail without corrupting the original pair.
Fixed out-of-neighborhood balance columns remain fixed.
2. Expanded local nonlinear repair previously stopped on the first
non-improving full step. It now retains the preceding point, quarters its
angle/voltage trust radii and retries. At most 12 LP attempts per neighborhood,
four-/six-hop neighborhoods, 2,048 control-bus cap, two-second construction-
inclusive LP cap and the same shared 16-second repair budget. Retry radius
has a finite minimum. Every accepted iterate reduces independently checked
nonlinear residual; success still requires the complete original validator.
Previously successful earlier fast paths are unchanged.

No source constraint, PMIN, ramp anchor, contingency, model slack cap,
acceptance tolerance, GPU work or dependency was changed.

## Tests and freeze

108 Python component tests and all four native CTest groups passed before
any V34 full run. Tiny tests cover signed roundoff at both active bounds,
genuine and nonfinite reconstruction failures, atomic failure preservation,
fixed out-of-neighborhood pairs, shrinking retries, improving/stagnant
steps, retry bounds, unchanged base state, full LP audits and full nonlinear
validation. No full network was used as a component test.

Executable SHA256:
`353439cb4b29a27b0a75cf4a006a894c65adf7bf5d0f06f8bd193bd21ab93efa`.

Configuration: `config/v34_all_scenarios_20260916.json`.
The campaign manifest records the frozen commit, configuration, executable,
source RAW/JSON/CON and normalized-model hashes, exact command and cold start.

## Acceptance and storage

Each run has the unchanged 300-second end-to-end cap: normalized-case loading
through official-equation evaluation and result serialization. Work cutoff
is 295 seconds. Require every source contingency, complete independent
native residual <=1e-5, official infeasibility zero, and complete official
detail/label coverage. Penalized imbalances remain visible in the objective;
PASS is feasibility, not global optimality or top-five competition rank.

Preserve logs, status, summaries, certificates and hash-backed compact
evidence. For space, after a complete verified same-input replacement,
retain the better-scoring complete solution and prune only redundant solution
text payloads after archive/hash/path/liveness checks. Preserve active runs,
source inputs, code and environments. Report pruning and actual physical
space change; a compact certificate cannot reconstruct deleted outputs.

All results are pending at registration. Do not substitute V33 results.
