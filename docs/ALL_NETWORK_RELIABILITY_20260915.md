# All-network reliability work — 2026-09-15

Status: **in progress; full-suite success is not yet established.**

## Scope and acceptance

The requested network families are 617, 2,020, **4,224** (source family
`C2FEN04200`), 8,300, 16,789, and 19,402 buses. The registered suite contains
37 Final Event Division 1 scenarios. Each scenario starts cold from its own
source case. Within-scenario corrective-state reuse is unchanged; solutions
are not carried between scenarios.

PASS requires all source contingencies, the independent native residual
check at `1e-5`, zero official infeasibility, exact official label coverage,
and at most 300 seconds from normalized-case loading through verification
and serialization. It is a feasibility claim, not a global-optimality claim
or a promise of zero source-permitted penalized violations.

## Preserved implementation and model

The reliability branch starts at `6f0afaaf7b238d67373fa0d8949c9648f9ac1faf`.
The historical checkout and the separate, newer economic-improvement branch
are preserved. Work and outputs are local, outside OneDrive.

Commit `8b5c6be` adds supporting half-planes for jointly overloaded active /
reactive branch flow in the feasibility-seed LP. The previous separate
component boxes did not enforce the combined terminal-current limit. These
are candidate-generation constraints; the unchanged independent AC checker
still decides acceptance. Generator limits, PMIN, source soft-violation
caps, contingency lists, physics, and acceptance tolerances are unchanged.

The reliability configuration also uses existing staged WSL worker inputs,
transient in-memory screening output, four initial corrective workers and
eight after screening. These scheduling/I/O choices are part of the tested
configuration; timing improvements cannot be attributed solely to the
new branch-loading cuts.

Commit `7fe014d` adds a timing-only schedule for 19,402-bus scenario 010.
Previously slow groups and its unfinished contingency receive earlier
dispatch. Unfinished tasks are explicitly censored observations, **not
invented measured runtimes**. No prior primal, dual, commitment, or network
state is included in that scheduling profile. The profile and input hashes
are frozen and recorded by the suite runner.

## Targeted development runs (not the final 37-scenario campaign)

| Buses | Scenario | Revision | Seconds | Completed/source contingencies | Official objective | Result |
|---:|---:|---|---:|---:|---:|---|
| 4,224 | 009 | 8b5c6be | 97.881 | 455/455 | 4,949,906.58 | PASS |
| 8,300 | 012 | 7326430 | 147.026 | 607/607 | 9,084,823.80 | PASS |
| 19,402 | 010 | 7326430 | 295.419 | 6,692/6,693 | Not certified | Deadline during CTG_002365 |
| 19,402 | 010 | 7fe014d | 295.101 | 6,693/6,693 | Not certified | Unused corrective-worker startup hit deadline |
| 16,789 | 094 | 7fe014d | 295.289 | 236/238 | Not certified | Two generator-outage repairs unfinished |
| 19,402 | 010 | 69596b2 | 296.761 | 6,693/6,693 | 164,882.98 | PASS |
| 16,789 | 094 | cc53fca | 295.129 | 236/238 | Not certified | Same two generator-outage repairs unfinished |
| 16,789 | 094 | f471795 | 295.131 | 236/238 | Not certified | IPM policy alone did not finish the cold run |
| 16,789 | 094 | 3790e89 | 185.940 | 238/238 | 10,359,780.01 | PASS |
| 19,402 | 006 | 3790e89 | 295.351 | 6,690/6,693 | Not certified | Two unfinished screens and one corrective LP |
| 19,402 | 010 | 3790e89 | 295.600 | 6,693/6,693 | 163,091.10 | PASS |
| 19,402 | 069 | 3790e89 | 213.551 | 6,620/6,620 | 627,542.21 | PASS |
| 19,402 | 077 | 3790e89 | 184.786 | 6,584/6,584 | 84,914.92 | PASS |
| 19,402 | 095 | 3790e89 | 214.897 | 6,579/6,579 | 183,630.65 | PASS |

The failed 19,402-bus test did not establish complete security. Its logs
show a 149.23-second solver task for CTG_001697 and a late unfinished
CTG_002365. This motivated the scheduling change, not a relaxed acceptance
rule. Completed-task counts never substitute for the final official check.

The scheduling-only replacement did complete CTG_002365 (225.49 seconds of
solver work), but then started four additional corrective workers despite
there being **zero fallback tasks**. One unused worker timed out during
startup and prevented finalization. The controller now claims a task before
launching a native corrective process and caps pool expansion by actual
queued work. This removes wasted work; it does not suppress a failure from
an assigned contingency or waive final evaluation.

The controller-fix run passed with all 6,694 official labels (base plus
6,693 contingencies), zero official infeasibility, and maximum independent
contingency residual `9.703835222196755e-6`. It has only 3.24 seconds of margin
to the 300-second cap; the full-suite robustness goal is still outstanding.

The separate 16,789-bus 094 issue is not that lifecycle bug: CTG_000007 and
CTG_000008 remained in corrective repair. Their saved fast-screen states
have roughly 2.87 p.u. active-balance residual at bus 14110. Further repair
work is required before registering the final regression campaign.

The next exact computational improvement reuses sparse symbolic analysis
within each ordinary or distributed-slack Newton call. Its topology and
PQ/reference mask are fixed within that call; numerical factorization and
pivoting still occur at every step. The cache is not shared across calls or
outages. Tiny lossless and lossy AC fixtures compare cached and uncached
solutions, iteration counts and symbolic-analysis counts; voltages/angles
must agree to `1e-12`. All four CTest groups and 59 Python tests pass.
The full 094 run still timed out after this exact computational improvement;
it is not sufficient by itself. An opt-in, bounded single-contingency trace
uses its saved base/fast-screen state to identify the stalled repair phase.
That diagnostic is explicitly separate from cold scenario acceptance.

The trace identified the actual 094 bottleneck: prelinear Newton repair
finished in 2.8 seconds; the subsequent large elastic Phase-I LP stalled
under simplex. A bounded diagnostic on CTG_000007 using the existing IPM
option completed in 94.495 seconds (87.851 seconds in that LP), followed by
exact AC repair with maximum residual `1.7763568394002505e-15`. This used a
saved base and fast-screen state and is **not** a cold scenario PASS.
The numerical-method policy now selects IPM for elastic balance Phase I;
economic and non-elastic large-network LPs retain their previous policy.
Matrix construction, LP tolerances, and physical acceptance are unchanged.
The result log records the actual selected LP method, including overrides.

The cold IPM-policy replacement still timed out. Its saved base voltages,
angles, active/reactive dispatch and demand factors exactly match the prior
diagnostic input. Replaying CTG_000007 alone again passed in 95.238 seconds
(88.603 seconds reported for LP plus postlinear Newton). This points to a
runtime/context issue rather than a changed base point, but does not prove
the cause of the cold timeout. Corrective worker logging is now enabled in
the suite. The elastic LP's inner budget is raised from 90 to 180 seconds;
the outer 300-second scenario deadline and five-second finalization reserve
are unchanged. Small and explicit-security repair LP budgets are unchanged.

With that headroom, cold 16,789-bus scenario 094 passed. Its two concurrent
elastic LPs took 109.424 and 110.972 seconds, both exceeding the former
90-second inner cap. The final native contingency residual was
`7.163773372287352e-6`, with all 239 official labels and zero infeasibility.

The ensuing five-case 19,402-bus batch passed four cases. Scenario 006
stopped with 6,690/6,693 completed. Fast-screen logs contain 6,691 task
acknowledgements, but one was a failed screen requiring corrective repair;
an acknowledgement is not a feasibility certificate. CTG_000262/000263 had
no completed fast-screen result. CTG_000179's corrective elastic IPM hit
180.145 seconds after seven iterations without an accepted candidate, after
22.32 seconds of prelinear Newton work. A timing-only profile now promotes
these unfinished tasks and the nine measured screens exceeding 30 seconds.
An analytic primal/basis simplex diagnostic is separate from cold runs.

That primal/basis diagnostic did not pass: its 150-second native limit
expired without an accepted corrective candidate (155.088 seconds including
diagnostic setup/termination). It is not enabled as a production default.
Inspection then found a concrete row-generation bug: base and contingency
security collectors used the state's reconstructed `sm_slack` even when it
exceeded the source cap. This masks the very branch violation that needs a
repair constraint. Both collectors now evaluate flow excess at the allowed
slack; they do not clip or alter the candidate, source cap, or final checker.
Tiny tests cover illegal versus permitted slack, both line terminals, voltage
dependence, transformer behavior, and nonmutation of candidate state.

Failure reporting now distinguishes absent final verification evidence from
an actual observed revision or executable-hash mismatch. All missing
evidence remains a failed acceptance gate.

The security-collector correction fixed the saved CTG_000179 diagnostic in
12.085 seconds, and that outage also passed in the next cold 006 run
(`target_19402_s006_v2`, c396f66) in 50.665 seconds under concurrency. The full
run still stopped at 295.356 seconds with 6,691/6,693 contingencies completed.
Its affinity schedule places CTG_000263 first and CTG_000262 second in the
same generator-bus group. The group's worker log contains only its ready
message: 263 had not returned, so 262 was not a second concurrently stalled
solve. Both still lack certification.

A saved-base direct-LP diagnostic on 263 timed out, as did a second diagnostic
that first tried a 20-second projected-balance LP (152.531 seconds including
termination). The short probe found no feasible LP point and is removed from
the production path; it is not a claimed improvement. A 90-second fast-only
diagnostic also timed out without returning its intermediate trace. The next
diagnostic mechanism therefore adds an opt-in cooperative predictor budget:
between iterations it can return the best rebuilt, still-unverified state for
corrective repair. Default execution remains unbounded internally, subject to
the existing hard scenario deadline. Fixed-Jacobian progress is logged only
when repair logging is explicitly enabled. These diagnostics do not count as
cold scenario passes or replace any of the final 37 required checks.

The first budgeted diagnostic returned a secure CTG_000263 candidate in
65.946 seconds, with independent maximum residual `1.7763568394002505e-15`.
It also exposed a control-flow bug: the nominally bounded Newton rescue ran
the full fixed-Jacobian predictor before its capped Newton steps. The first
predictor stopped around 20 seconds, then the rescue repeated 114 predictor
iterations. A dedicated option now disables that predictor only for the
bounded Newton rescue; ordinary fast screens retain their existing search.
The next diagnostic passed in 25.451 seconds via compact linearized repair,
again with residual `1.7763568394002505e-15`. Its solve time was 23.732 seconds;
the rest was diagnostic setup and serialization. These are saved-base
diagnostics, not official cold-run times or score improvements.

The diagnosed 263 task now has an explicit 20-second predictor-stage budget
in the hash-bound schema-4 scheduling profile. This is a cooperative handoff
between iterations, not a relaxation of the global 300-second cap. The
Newton/compact repair still needs to pass complete AC validation. No other
task receives this budget; 262 remains an unmeasured queued sibling. A
task-local predictor instance prevents the override leaking into other work,
and the worker must acknowledge the requested budget. Tiny tests cover
invalid budgets, failure-status preservation, the rescue option, and profile
identity/range checks. Diagnostic stdout and stderr are now merged through a
pipe with one file writer, avoiding overlapping WSL file-output offsets.

## Storage pruning during development

After the five-case batch ended and solver processes were absent, checked
exact target containment, absence of reparse points, and archived status /
adjacent-JSON hashes for old 19,402-bus 010 targets v1, v2, and v3. Removed
only their 40,284 `solution_*.txt` file entries (74,185,604,185 logical bytes,
including duplicate/hard-linked payload entries). C: free space increased
from 49,013,534,720 to 86,038,650,880 bytes: 37,025,116,160 bytes recovered.
All archived JSON hashes still matched after pruning. Source, code, logs,
internal JSON states, summaries, certificates, and the latest five-case
batch were retained. Deletion was permanent; certificates retain the record
of verification, not the removed full solution text for direct re-evaluation.

Run evidence is under
`C:\Users\thoma\Documents\gravityx-go2-cpp\runs\reliability_20260915`.
The two `target_*_v1` controller runs have hash-backed evidence archives,
while the 4,224-bus development run retains its full summary and certificate.

## Tests and publication

The branch-loading repair passed the native tiny-fixture regression. The
scheduling and provenance update passed 57 Python tests and all four CTest
groups. These are component checks, not replacements for full-network runs.

Verified code milestones are committed locally. Publication is currently
blocked: the configured GitHub identity received HTTP 403 when pushing to
`thomas-digregorio/gravityx-go2-cpp`. No alternate credentials were attempted.

The final result must use one common frozen implementation/configuration
and report every scenario, including any failure. Earlier successful runs
will not be substituted for missing runs in that final campaign.
