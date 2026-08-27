# C2FEN19402 fixed-commitment economic pilot — stopped at time gate

Date: 2026-08-27

## Scope and frozen identities

- Frozen retained baseline: `c477a2d4673a87d4d569110471956e420bb919e6`.
- Economic-refinement implementation: `4651517464f4c4cf31e6902a186f87645548c947`.
- Native executable SHA-256:
  `f2edb09de2b71b0b71b6e6e6a39fa58880f24587a77f65582c610dee02df5926`.
- Normalized scenario 095 SHA-256:
  `11f85a7dbef5fbf9e8f34e14f7cd4df97fca9f62ea779837b3c5b6d313b3225b`.
- Raw source SHA-256:
  `bc9e819e4f34e2271e2bf9362642c01b3e7298a0c589066ddbf2c3514532daa3`.
- Contingency source SHA-256:
  `203f684a3ef1168040da9e1f683419924e65e5e75c8dd9b991b4ccc6e4125091`.
- The retained baseline run and evidence were not modified or rerun.

## Pre-run verification

Before the sole full pilot, both standard and native C++ test suites passed
4/4, the Python suite passed 49/49, and a tiny fixed-commitment fixture
accepted an independently verified market-surplus improvement. The worktree
was clean and the implementation milestone was pushed before the pilot.

The scenario 095 pilot used the retained execution pipeline: one corrective
worker, 24 fast-screen workers, two-stage screening, all 6,579 source
contingencies, C++ solution writing, overlapped referenced evaluator equations,
and independent `1e-5` verification. The only algorithmic addition was a
90-second fixed-commitment base economic-refinement request. The registered
end-to-end limit was 300 seconds, with 2 seconds each reserved for evaluation
and finalization.

## Result

| Gate or measure | Result |
|---|---:|
| Pilot scenario | 095 |
| Frozen objective | -36,295,267.953619 |
| Requested economic-refinement budget | 90.000 s |
| C++ inner deadline | 295.508 s |
| Base subprocess wall time | 305.528 s |
| Recorded end-to-end wall time | 306.020 s |
| Base solution serialized | No |
| Completed contingencies | 0 / 6,579 |
| Official evaluation reached | No |
| Independent exhaustive verification reached | No |
| New objective | Not produced |
| Outcome | **FAIL — time gate** |

The controller terminated the C++ base subprocess after its 295.508-second
inner deadline. The only base console artifact is `OUTER_TIMEOUT`; no
`base.json` was completed. The run-status evidence is retained at:

`runs/c2fen19402_economic_4651517_20260827/C2FEN19402_s095_cold/run_status.json`

Its SHA-256 is
`91611bc8c2544ec05c328fd74373993a7ff825d14a524a98c4a129836c3167ff`.

This is neither an infeasibility result nor an economic score. The run did not
reach contingency screening, official evaluation, or independent security
verification.

## Bottleneck diagnosis

The economic budget is not currently a complete wall-clock budget for the
economic phase. The implementation starts its timer, constructs a full
`AcModel`, installs the verified incumbent, and only then computes the
remaining time passed to `AcModel::solve`. Inside `solve`, Ipopt's
`max_cpu_time` constrains `OptimizeTNLP`, but model construction and the Gravity
to-Ipopt program update occur outside that solver limit.

Confirmed facts are therefore:

1. The full base subprocess consumed the global deadline.
2. No serialized base result was available for the fallback pipeline.
3. The 90-second solver limit did not bound all pre-solve construction and
   translation work.

Because subprocess output was buffered and the process was killed before the
phase record was serialized, this evidence cannot distinguish how much of the
unbounded time was in the `AcModel` constructor versus the pre-`OptimizeTNLP`
program update. Calling the whole interval "Ipopt solve time" would be
incorrect.

## Stop decision

The required pilot gate was not met. In accordance with the one-run and
stop-on-failure policy, scenarios 006, 010, 069, and 077 were not started; no
retry, tuning run, or broader six-family regression was performed.

