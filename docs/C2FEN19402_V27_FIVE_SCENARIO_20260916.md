# V27 cold five-scenario result: two pass, three work-deadline failures

Each scenario ran once, cold, with the exact source inputs, all source
contingencies, a 300-second normalized-input-to-verification limit and the
unchanged 1e-5 native residual gate. No previous solution was initialization.
The native algorithm change is `7095c57`; binary SHA256 is
`ac786b023bcc1ad7ed90789e10e40ace414be5ba1e246e1b48da27c0fd936f48`.
Config: `config/c2fen19402_top5_branch_aware_backtrack_v27.json`.
095 ran at `7b59ec0616e4a1b354867ff4e6ab5ba6a90922fa`; the other four ran at
`c6bcb2a11d8d77551e22de7d142d3c8dffc8d0fd`. The intervening commits contain
only documentation/evidence: native, controller, tests and configuration
sources are identical. All recorded binary and source-input hashes match.

| Scenario | Result | On-time contingencies / expected | End-to-end s | Accepted official objective | Published fifth-best objective | Shortfall, passing results only |
|---|---|---:|---:|---:|---:|---:|
| 006 | Work deadline | 987 / 6693 | 295.3469361999996 | None | 436687.482340259 | Not comparable |
| 010 | Work deadline | 1608 / 6693 | 295.3634312000031 | None | 534621.593420241 | Not comparable |
| 069 | Work deadline; one unresolved first screen | 2281 / 6620 | 295.3334049999976 | None | 949809.594696075 | Not comparable |
| 077 | PASS | 6584 / 6584 | 177.62342740000167 | 373417.7496024141 | 392060.669448444 | 18642.919846030 |
| 095 | PASS | 6579 / 6579 | 295.39416020000135 | 382016.64361338085 | 434319.522756719 | 52302.879143338 |

The failed runs stopped at the reserved Code2 work deadline, not because a
complete infeasibility proof was obtained. They have no accepted official
score. Their late native records are not credited. The 069 screen for
`CTG_001525` returned `fast_newton_screen_failed` with residual
0.5666488524298401 and was routed for fallback; it is NOT a feasible
contingency result or proof that the source problem is infeasible. Its
worker-level success flag indicates successful task processing, not physical
feasibility. The controller correctly withheld scenario success.

077 has 6,585/6,585 unique official detail labels including the base, official
infeasibility zero, and max native residual 7.434533043193525e-6. Its audited
score is base 191214.2314162511 plus contingency mean 182203.518186163.
All 47 base-detail copies agree. All 921 units retain source commitment;
there is no startup/shutdown. Its full solution vectors are retained.

095 has 6,580/6,580 unique official detail labels, official infeasibility zero,
and max native residual 9.90449499832513e-6. Its independent audit reproduces
382016.64361338085 exactly. This improves the previous accepted V25 campaign
result by 768.334186211. Full solution vectors are retained.

GO2 source semantics permit bounded, penalized imbalance. These are verified
source-feasible results, not zero-imbalance or globally optimal certificates.
Neither passing result reaches its fifth-place target, and the three
incomplete cases cannot be scored. The all-five goal remains open.

## Timing diagnosis and next work

| Scenario | Native records, including late | Mean task solve wall s | Mean returned solver iterations |
|---|---:|---:|---:|
| 006 | 988 | 5.554465084 | 6.228744939 |
| 010 | 1609 | 3.433470747 | 4.812305780 |
| 069 | 2288 | 2.407218761 | 3.606643357 |
| 077 | 6584 | 0.328297972 | 0.995291616 |
| 095 | 6579 | 0.787229002 | 1.339413285 |

These are summed/concurrent worker diagnostics, not extra benchmark runs or
additive end-to-end timing. Incomplete-case samples are biased toward the
registered early task order and cannot establish full-case average runtime.
For 006, economic polish consumed only 228.425895068 of 5487.811502992
summed task seconds. The first predictor path, including its checks/polish,
accounted for 5180.406326578 seconds. Thus micro-optimizing the final
branch-step proposal alone cannot resolve this deadline failure.

The economic route deliberately starts fresh from this run's optimized base;
it probes the resident verified-seed bank only after that predictor fails,
except for exactly equivalent parallel outages. This avoids the poor-quality
rolling-seed propagation observed in V1, but can spend many iterations on
each difficult outage before consulting a feasible within-run alternative.
A next experiment should measure and bound that feasibility effort, preserve
the best fully verified candidate, and consider quality-ranked within-run
seeds without reinstating unchecked or cross-run initialization. No such
policy change or replacement run is included in V27's reported results.

## Evidence

All files below are in `docs/evidence/`; every provenance archive reports
zero diagnostics. Source/input, code, config, log and status hashes are
preserved in those archives. Hashes listed below were checked before pruning.

| Scenario | Archive | SHA256 |
|---|---|---|
| 006 | C2FEN19402_TOP5_V27_S006_WORK_DEADLINE_20260916.json | 3c05abab7b3f109aedb3327e4145dde84d2ce1845b638175a626efdd5e70ef32 |
| 010 | C2FEN19402_TOP5_V27_S010_WORK_DEADLINE_20260916.json | 89d24ac85592ad57739984bf120ce051612840914fab913456af4f64c853b343 |
| 069 | C2FEN19402_TOP5_V27_S069_WORK_DEADLINE_20260916.json | bad3dc57a93ffc91a96f8d4ddbf0e1cf4ddc0ceb50603b4dd6a423effcf7e3ad |
| 077 | C2FEN19402_TOP5_V27_S077_VERIFIED_20260916.json | a6214160255f0faece02126ccaf0f6d56d72f12f2a93e253b71e1525a2498ba7 |
| 095 | C2FEN19402_TOP5_V27_S095_VERIFIED_20260916.json | 2b83b711299b2e054b5c89051ded4e787cf12f9b8a9a57458c534bb0e9d30228 |

077 independent audit SHA256:
`85dfedca3a17b67c5e42084ac5787956566adb8146572c7e974814eea8aadcb1`.
095 independent audit SHA256:
`dfd9c93da796f7c27c5ffa26b5750bbbf9f02275a7046b5aa492277594b4e5b2`.

Run-status hashes:

- 006: `73c2f47afe3302ca78684defd406bf2ea23120998f366311862ced9a7cf3c778`
- 010: `3708a72f351f1fda9397f8e3aefb613b03d123c6a2e81542c114a43d3918fec5`
- 069: `ef66991498cebfcc748164d2cae58e27461908b28bb4592d908aff464e558596`
- 077: `0fe457a807d5c48ed70a927307c90b1ee1c6766f444f912cb46c4c335a14c658`
- 095: `3b82c89cbbeeb83abee5a9d564246c6a63ed785c718180f80ecddf412cbc5e60`

Retain the prior verified full outputs for 006, 010 and 069, and the V27 full
outputs for 077 and 095. After committing/pushing this evidence, only the
three V27 Code2-timeout runs' disposable solution texts may be pruned using
the guarded policy. Their logs, failed-state JSON, partial official records
and timeout provenance must remain. That pruning loses exact re-evaluation
of deleted partial vectors, not any accepted scenario result or algorithm.

That guarded cleanup is complete. The three timeout runs' 9,663 solution-text
entries were removed, recovering 8.2671051025 GiB; their status hashes and
diagnostics remain unchanged. The superseded original 077 payload was also
archived, verified and pruned, recovering 10.9777221680 GiB. The protected
complete sets are now old-reference 006/010/069 plus V27 077/095. Their full
vectors remain available; old exact 077 and failed partial vectors cannot be
recreated from certificates alone. Detailed paths, guards, archive hashes
and physical-space readings are recorded in
`docs/CLEANUP_AND_V19_RESUMPTION_20260916.md`.
