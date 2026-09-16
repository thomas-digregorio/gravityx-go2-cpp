# V32 cold reliability and smaller-network regression

## Acceptance and freeze

The user's current request is to obtain complete feasible passes on all five
19,402-bus scenarios, then test all scenarios on the five smaller families.
Objective quality is reported separately; a feasibility PASS is not a top-five
competition placement or a global optimality certificate.

Native algorithm freeze: `234b8adf39f8464b6796ccf4af0eed8b5e874221`.
Release executable SHA256:
`bbf9f3f8ce1c8259a79b8fbfba8868847e455de312bd4925387589f468ee138f`.
All 108 Python component tests and four native CTest groups passed before
the first V32 full run. No algorithm changes or repeated same-version runs
are authorized by this registration.

Each case starts cold from its own immutable source inputs. PASS requires
all source contingencies, independently checked native residuals <=1e-5,
official infeasibility zero, complete official detail coverage, and end-to-end
time <=300 seconds. The existing work/evaluation/finalization deadlines stay
in force. Penalized imbalances permitted by GO2 remain visible in the score;
they must not be described as zero physical mismatch.

## 19,402-bus results

Configuration: `config/c2fen19402_top5_local_assembly_v32.json`.

| Scenario | Contingencies | Seconds | Official objective | Max native residual | Result |
|---|---:|---:|---:|---:|---|
| 006 | 6,693 / 6,693 | 264.9913063000058 | 332,087.2988758405 | 9.952236352528399e-6 | PASS |
| 010 | 6,693 / 6,693 | 256.152768699998 | 470,268.793058862 | 9.990036206664055e-6 | PASS |
| 069 | 6,619 / 6,620 | 295.1634219000043 | Not verified | Not complete | Work-deadline failure |
| 077 | Pending | Pending | Pending | Pending | Not started |
| 095 | Pending | Pending | Pending | Pending | Not started |

006 retains all 6,694 full solution vectors, including the base case.
The verified provenance archive SHA256 is
`bf894c74f793e4c5a36f0dac16d770a74ace164903b105fcd05dfcf9750dcb4f`;
the independent objective-component audit SHA256 is
`c96ff1cda048fdcaeee8af2c838c7d729f7c876503eedf3d02261cc502b6a832`.
The reconstructed official objective agrees exactly with the reported value.
The fifth-best published objective for 006 is 436,687.482340259, so this
verified result does not satisfy the separate top-five objective target.

## Storage protection

After verifying the replacement 006 full solution, complete official
certificate, matching raw/normalized input hashes, compact archives, absence
of active workers and resolved non-reparse paths, the superseded verified
006 run's solution-text payload was pruned using the guarded repository
helper. No source data, code, environment, logs or certificates were deleted.

Superseded run:
`runs/reliability_20260915/target_19402_s006_v6/C2FEN19402_s006_cold`.
Replacement:
`runs/top5_cold_20260915/v32_s006/C2FEN19402_s006_cold`.

13,434 text-file entries were removed (22.974 apparent GiB including linked
entries); physical C: free space increased from 40,044,257,280 to
52,347,658,240 bytes, recovering 11.45843505859375 GiB. Zero old vector
entries remain and all 6,694 replacement root vectors remain. The old
status hash is unchanged. Deleted vectors are not recoverable from compact
certificates; those preserve the verification record, not the full outputs.

Old-run pre-prune evidence SHA256:
`a17725e1c17d0c54f593634e5a39be04f0805b857550c836bf7d0f9d31166a0f`.

The same guards were applied before replacing the superseded verified 010
reference at `runs/reliability_20260915/frozen_a1bebd5_19402_s010/C2FEN19402_s010_cold`
with the full V32/010 solution. The replacement objective audit agrees
exactly, input hashes match, and its objective exceeds the older reference.
13,434 old vector entries were removed (23.041 apparent GiB); physical free
space increased from 39,814,393,856 to 52,145,266,688 bytes, recovering
11.484020233154297 GiB. All 6,694 replacement root vectors remain.
Old status/evidence are unchanged; old vectors are permanently removed.

010 compact evidence SHA256:
`3fbd23fc9d2ac372429f88856ef6610ae4114e02134eb97fbe8d2009902cad91`.
010 objective audit SHA256:
`97937bd862bc3366a8f5e3b3364408e368c83b70679e5ed7c04648f8cc4c52cc`.
Superseded 010 pre-prune archive SHA256:
`b33216a7a3f5225c9c3336a33edcf1e3241ddd5be593047e4dd338a1b056823e`.

## Smaller-network regression

Not started. V32 did not pass all five scenarios, so its 077/095 runs and
smaller-network regression were not started. After all five scenarios pass on one unchanged
algorithm, test all 32 smaller scenarios in order: 617 (5), 2,020 (5),
4,224 (7), 8,300 (7), and 16,789 (8). Source family directory identifiers
are respectively `00617`, `02020`, `04200`, `08300`, and `17700`.
No historical results will substitute for these new frozen-version tests.

## 069 diagnosis

The only unfinished case was branch outage `CTG_001525` (source branch 4840,
between buses 39948 and 43629). The initial 85-control-bus, three-hop local
LP was infeasible. The subsequent fast nonlinear search stopped with
reactive-balance residual 0.5666488630117348 at bus 39948 and variable-bound
residual 0.24627610027355845. This exceeded the bounded Newton/linear-repair
candidate threshold, so it required the full-network corrective fallback.
That fallback spent 90 seconds on a simplex feasibility LP without a valid
primal, then was interrupted during another large LP at the global work
deadline. This is algorithm failure, NOT a proof that the source outage is
infeasible. The complete original validator never accepted the candidate.

The four- and six-hop neighborhoods contain 220 and 1,151 control buses.
V33 will add a bounded wider-neighborhood repair only after all existing
fast paths fail. Previously successful paths, the source model, every
contingency, and all acceptance tolerances remain unchanged.

069 timeout archive SHA256:
`ae3f8bd75c2d8785b4b9cae04ac6c475f1252f09c1d909dd84ebdcb4cb765786`.

After archive/status verification, no-reparse/local-containment and worker
liveness checks, pruned only the failed V32/069 solution texts: 13,270 file
entries, 22.077 apparent GiB. Physical free space increased from
40,117,534,720 to 51,941,023,744 bytes (11.011482238769531 GiB recovered).
Zero failed-run vector entries remain. The prior full verified 069 reference
retains all 6,621 root vectors and its complete certificate. Failed-run logs,
base/fallback JSON, status and partial evaluator evidence remain; there is
no complete official certificate for this failed attempt. Deleted texts are
not recoverable from the compact archive. Status SHA256 remains
`b6f24ac5a77348ea151a038c1251729230915fafb3df295b9ed55c8b6d69992a`.
