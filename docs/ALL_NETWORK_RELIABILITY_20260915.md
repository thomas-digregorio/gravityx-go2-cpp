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

The separate 16,789-bus 094 issue is not that lifecycle bug: CTG_000007 and
CTG_000008 remained in corrective repair. Their saved fast-screen states
have roughly 2.87 p.u. active-balance residual at bus 14110. Further repair
work is required before registering the final regression campaign.

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
