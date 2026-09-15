# Frozen a1bebd5: all-network reliability results

**37/37 scenarios PASS.** All six requested network families completed on one
frozen CPU implementation and configuration, with all **44,517 source-listed
contingency evaluations** and 37 base cases verified. No scenario exceeded
300 seconds. The largest independent residual was
`0.00000997666398094843`, below the unchanged `1e-5` limit.

## Frozen identity and execution boundary

- Algorithm commit: `a1bebd5e3470d2ad760951bde010a0e37c8435df`.
- Native executable SHA256: `6094158e183aded78a8e1b4f47bd09c5c332f26ed90e7205772c96656c7e7e83`.
- Configuration SHA256: `a127c128de718b35608d799a44ef74bbdebe572427cdcdcfb47c380dd3abf33b`.
- Configuration: `config/reliability_suite.json`; CPU only, no GPU work.
- Scope: GO Challenge 2 Final Event, Division 1; 37 registered scenarios.
- Each scenario starts cold from its own source case. Within-scenario state
  reuse is allowed; no solution is transferred from another scenario.
- Time runs from normalized-case loading through independent verification
  and serialization. Setup and prior normalization are excluded.
- Each scenario has the same hard 300-second end-to-end limit. Solver work
  stops early enough to reserve evaluation/finalization time.
- No algorithm or configuration changes occurred during these 37 recorded
  runs. The 006 record is `target_19402_s006_v6`, its first run on this exact
  revision. Other records are under `frozen_a1bebd5_*`. Earlier development
  attempts and historical best scores are not substituted.
- 62 Python component tests and all four native CTest groups passed before
  the frozen campaign. Component tests are not full-network repetitions.

The sum of the 37 measured scenario times is **3,759.925 seconds**
(62 minutes 39.925 seconds). This is not total project elapsed time: setup,
diagnostics, development attempts, storage cleanup and gaps between runs are
excluded.

## Network summary

| Buses | Scenarios passed | Total source contingency evaluations | Fastest case (s) | Slowest case (s) |
|---:|---:|---:|---:|---:|
| 617 | 5/5 | 522 | 16.86 | 20.52 |
| 2,020 | 5/5 | 1,492 | 52.73 | 59.38 |
| 4,224 | 7/7 | 3,185 | 79.22 | 92.66 |
| 8,300 | 7/7 | 4,245 | 90.85 | 200.99 |
| 16,789 | 8/8 | 1,904 | 74.62 | 84.72 |
| 19,402 | 5/5 | 33,169 | 184.92 | 294.77 |

## Every scenario

Objectives are the official evaluator's reported objectives, not lower-bound
proofs. The residual column is the maximum of the independent base and
contingency residuals for that scenario.

| Buses | Scenario | Source contingencies | End-to-end seconds | Official objective | Max independent residual | Result |
|---:|---:|---:|---:|---:|---:|---|
| 617 | 005 | 105 | 20.517 | 1,214,187.43 | 3.301077e-8 | PASS |
| 617 | 017 | 103 | 16.856 | 1,034,613.22 | 8.257222e-17 | PASS |
| 617 | 024 | 104 | 19.203 | 1,425,466.14 | 2.220446e-16 | PASS |
| 617 | 062 | 103 | 18.328 | 1,337,215.52 | 2.490394e-16 | PASS |
| 617 | 073 | 107 | 17.506 | 1,122,605.37 | 5.551115e-17 | PASS |
| 2,020 | 025 | 292 | 52.730 | 4,060,248.92 | 8.881784e-16 | PASS |
| 2,020 | 121 | 300 | 54.693 | 6,322,492.50 | 1.421085e-14 | PASS |
| 2,020 | 134 | 300 | 52.732 | 4,817,183.06 | 1.776357e-15 | PASS |
| 2,020 | 260 | 300 | 52.942 | 5,898,950.44 | 1.776357e-15 | PASS |
| 2,020 | 262 | 300 | 59.383 | 5,755,551.87 | 1.421085e-14 | PASS |
| 4,224 | 009 | 455 | 92.664 | 4,948,246.70 | 1.421085e-14 | PASS |
| 4,224 | 010 | 455 | 79.224 | 5,178,256.41 | 2.273737e-13 | PASS |
| 4,224 | 014 | 455 | 92.306 | 4,922,287.09 | 2.842171e-14 | PASS |
| 4,224 | 055 | 455 | 79.469 | 5,946,375.64 | 1.421085e-14 | PASS |
| 4,224 | 056 | 455 | 84.104 | 6,262,077.40 | 2.842171e-14 | PASS |
| 4,224 | 057 | 455 | 85.616 | 6,334,170.19 | 2.842171e-14 | PASS |
| 4,224 | 060 | 455 | 84.462 | 5,937,600.97 | 2.842171e-14 | PASS |
| 8,300 | 003 | 607 | 90.851 | 20,972,782.31 | 2.273737e-13 | PASS |
| 8,300 | 012 | 607 | 160.542 | 8,844,425.81 | 8.227195e-6 | PASS |
| 8,300 | 013 | 607 | 142.887 | 9,339,715.63 | 8.588511e-6 | PASS |
| 8,300 | 022 | 607 | 141.648 | 11,505,642.58 | 8.900913e-6 | PASS |
| 8,300 | 043 | 607 | 156.261 | 9,311,479.92 | 5.342846e-6 | PASS |
| 8,300 | 052 | 607 | 113.971 | 11,253,967.25 | 9.675693e-6 | PASS |
| 8,300 | 166 | 603 | 200.988 | 11,318,485.79 | 7.533006e-6 | PASS |
| 16,789 | 019 | 238 | 84.464 | 9,771,784.03 | 1.136868e-13 | PASS |
| 16,789 | 020 | 238 | 76.395 | 10,204,842.35 | 5.511475e-6 | PASS |
| 16,789 | 021 | 238 | 78.915 | 9,656,728.64 | 1.705303e-13 | PASS |
| 16,789 | 089 | 238 | 75.180 | 10,152,888.60 | 2.842171e-14 | PASS |
| 16,789 | 094 | 238 | 84.716 | 8,489,982.72 | 7.163773e-6 | PASS |
| 16,789 | 106 | 238 | 80.371 | 9,992,600.48 | 6.517924e-6 | PASS |
| 16,789 | 107 | 238 | 78.299 | 10,164,368.13 | 2.842171e-14 | PASS |
| 16,789 | 115 | 238 | 74.619 | 11,139,565.24 | 1.705303e-13 | PASS |
| 19,402 | 006 | 6,693 | 294.770 | 300.64 | 9.956144e-6 | PASS |
| 19,402 | 010 | 6,693 | 263.388 | 165,985.63 | 9.703835e-6 | PASS |
| 19,402 | 069 | 6,620 | 209.998 | 627,636.09 | 9.976664e-6 | PASS |
| 19,402 | 077 | 6,584 | 184.921 | 86,153.33 | 9.663354e-6 | PASS |
| 19,402 | 095 | 6,579 | 204.006 | 173,320.66 | 8.859175e-6 | PASS |

## What PASS does and does not mean

Every recorded run has:

1. Terminal success, with all source contingencies completed.
2. A complete official label set including the base case, without infeasible
   labels, and official infeasibility exactly zero.
3. Independent native base and contingency residuals no greater than `1e-5`.
4. A finite official objective and end-to-end time no greater than 300 seconds.
5. Matching frozen commit, executable, configuration, source RAW/JSON/CON,
   normalized-model and applicable scheduling-profile identities.

This establishes feasibility under the source GO2 formulation, including its
permitted and penalized soft violations. It is **not** a global optimality
certificate, a guarantee of zero soft violations, or a claim that economic
quality improved for every scenario. In particular, 19,402-bus scenario 006
has a weak official objective of **300.64** despite passing feasibility.
Scenario 095 is **173,320.66** in this campaign, not the historical 210,762.42
from a different run. Those historical scores must not be substituted here.

The final evidence audit reread all 37 statuses, summaries and certificates,
verified their archive hashes, rehashed the frozen inputs, independently
recounted the source contingency declarations, and checked exact registered
case coverage. It found 37 distinct expected cases, all passing. All solver
processes exited; no further solve is queued.

## The original 4,224-bus base-validation failure

In the older `6f0afaa` campaign, all seven cases stopped before contingencies.
The candidate needed about `0.255702646` of branch-loading slack at
`branch:366:sm_slack`, but the source model allowed only `0.2`. The remaining
bound violation was about `0.055702646`, far above the acceptance tolerance.

The repair LP could report an optimal approximate solution while that state
still failed the full AC check. Its separate active/reactive component bounds
were insufficient to control combined terminal loading. The reliability
implementation adds joint-flow supporting constraints during candidate
generation and still requires the independent full AC check. It does not
increase the source slack cap or relax PMIN, physical constraints, contingency
coverage, or the validation tolerance.

The old runs stopped in roughly 44-48 seconds, not at the five-minute deadline.
They demonstrated a failed repair path, not mathematical proof that the case
was infeasible. All seven now pass on the common revision above.

## Local evidence, storage and publication

The full current-run records remain under:

`C:\Users\thoma\Documents\gravityx-go2-cpp\runs\reliability_20260915`

The compact, machine-readable audited record is
[FROZEN_A1BEBD5_ALL_NETWORK_PASS_20260915.json](evidence/FROZEN_A1BEBD5_ALL_NETWORK_PASS_20260915.json).
It retains exact objectives, times, residuals, source identities, local run
paths and hashes for each status, summary, certificate and archive.

During the final campaign, older completed/failed artifacts were pruned only
between runs after process, containment, reparse-point and archived-hash checks:

- Failed 006 v4/v5: 26,376 solution-text entries permanently removed;
  24,221,732,864 physical bytes recovered. Their logs, internal states and
  archived failure evidence remain.
- Older successful 3790e89 scenario 010: 13,434 solution-text entries
  permanently removed; 12,338,941,952 physical bytes recovered. Its certificate,
  logs, internal states and archived records remain; the newer successful
  scenario 010 remains in full.
- Earlier development pruning is documented in
  [the reliability work log](ALL_NETWORK_RELIABILITY_20260915.md).

No source data, code, environments or current frozen-campaign solution text
was deleted by those operations. Archived verification records do not replace
deleted full outputs for direct re-evaluation. Free C: space was about
49.18 GB after the last run.

Verified changes are committed locally on `codex/all-network-pass-20260915`.
The configured GitHub identity previously received HTTP 403 when pushing;
publication remains blocked. No alternate credentials were used and no
successful push is claimed.
