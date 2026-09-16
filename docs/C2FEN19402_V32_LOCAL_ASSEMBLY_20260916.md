# V31 result and V32 bounded-work improvements

## V31 single cold 006 result

Frozen revision `2674c8599ed6f46f08cd80ef59860babdbe14d93` stopped at the
Code2 work deadline after 295.26296379999985 s. It completed 5,122/6,693
contingencies on time, versus V30's 870. Base algorithm time was
51.580595364 s. This is improved coverage, NOT a complete SCOPF solution,
official objective, certified speedup or satisfied top-five target.

The 5,127 terminal native records include five late records not credited by
the controller. All terminal records passed their physical checks. Methods:

| Method | Records |
|---|---:|
| Local feasible incumbent plus cached economic polish | 4,286 |
| Quality-gated local LP | 224 |
| Exact parallel direct reuse | 459 |
| Quality-gated seed-bank reuse | 57 |
| Fresh predictor/economic polish | 86 |
| Bounded linear feasibility repair | 13 |
| Linear repair plus Newton | 1 |
| Bounded Newton rescue | 1 |

All 4,286 local-incumbent polishes entered at predictor iteration zero; none
lost incumbent objective. There were 9,003 LP attempts: 8,971 Optimal and 32
Infeasible, with no iteration-limit exits. Mean reduced column/control-bus
counts were 774.379984449628 / 56.196378984782854.

Summed worker-seconds (nested, not additive and not elapsed time):

| Counter | Seconds |
|---|---:|
| All native tasks | 5,442.366621916 |
| Local route including checks and polish | 3,655.538566491 |
| LP attempts | 1,805.434452865 |
| Local-incumbent economic polish | 1,165.656905812 |
| Priority quality-seed probes | 675.137381847 |

Archive `docs/evidence/C2FEN19402_TOP5_V31_S006_WORK_DEADLINE_20260916.json`
has zero diagnostics, SHA256
`8d89e320b63ccea015c46edbe25b6450805f5780ae6d751b9ebd46047064ef32`.
Status SHA256:
`751f4e84d9b183e97c41742569b34d4dd06c2fe8c84475d07c7650c3a868aafd`.
Source inputs, executable and configuration match the frozen registration.
No other V31 scenario ran.

## Storage retention

After archive/hash/zero-diagnostic checks, process-liveness checks, resolved
local-path/no-reparse checks, and verification of the protected full 006
reference/certificate, 10,225 disposable V31 solution-text entries were
removed. Apparent size including links: 17.533 GiB. Physical free space went
from 43,135,950,848 to 52,593,684,480 bytes, recovering 8.80820083618164 GiB.
The retained status hash is unchanged and zero V31 vector entries remain.
Logs, status, partial evaluator records, base JSON and compact provenance
remain. Deleted vectors cannot be reconstructed from these records. Source
data, Git/code/dependencies, unrelated user files and protected complete
solutions remain untouched.

## V32 design

The local LP changes about 56 buses on the observed subset, yet evaluates
finite-difference derivatives across the complete 34,899-branch case before
substituting most variable increments to zero. V32 omits derivative evaluation
only when BOTH endpoint voltage/angle increments are exactly fixed at zero.
It retains all actual terminal flows, balances, source bounds, constant-row
infeasibility checks and expanded LP audits. Derivatives touching a movable
bus are unchanged. This is exact substitution in the restricted candidate
LP, not deleting a physical branch or contingency. Every nonlinear acceptance
and final official check still covers the complete original network.

The full-derivative assembly remains an explicit component-test oracle.
Counters record evaluated/skipped branch derivatives and assembly time. A
tiny connected zero-injection tail supplies fixed/fixed branches; tests must
compare reduced dimensions, objective and full nonlinear residuals against
the full-derivative oracle.

V32 also disables the BEFORE-local-repair quality seed-bank probes in the
production cached mode: V31 spent 675.137 worker-seconds to select just 57
such records. Exact equivalent parallel reuse stays first. The same-run bank
and its later failed-predictor fallback remain; no prior-run input is added.
The optional priority route and its tests remain available. This ordering
change is heuristic and may affect cost; only new complete official results
can establish an improvement. The incumbent-quality floor, local LP bounds,
two proposal rounds, time/pivot budgets and monotone bounded cleanup stay.

Source cold start, immutable base ramp anchor, source equations/PMIN/costs,
full contingency set, 1e-5 native residual gate and 300-second end-to-end
limit remain unchanged. No GPU, commercial solver or new dependency is used.

## Registration

All 108 Python tests and four native CTest groups pass. The connected-tail
fixture exercises two fixed/fixed branches and confirms identical reduced
dimensions, LP objective and full nonlinear residual against the complete
derivative oracle. The worker test confirms that disabled priority probes
perform zero probes and retain independently checked local routing. Existing
exact-parallel, immutable-anchor, constant-row rejection and source-bound
tests continue to pass. No full network was used as a development test.
Release executable SHA256:
`bbf9f3f8ce1c8259a79b8fbfba8868847e455de312bd4925387589f468ee138f`.

After component tests and a pushed freeze, run 006 once using
`config/c2fen19402_top5_local_assembly_v32.json` into
`runs/top5_cold_20260915/v32_s006`. Only a full verified pass permits the
unchanged version's one-each 010, 069, 077 and 095 runs. Preserve the current
protected full solutions and archive/hash-check terminal evidence before
pruning disposable outputs. All five top-five targets remain open.
