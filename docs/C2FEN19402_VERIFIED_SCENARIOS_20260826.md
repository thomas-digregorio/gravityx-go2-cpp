# C2FEN19402 verified final-event scenarios

This note records only cold runs that completed the complete source contingency
set, passed the referenced official evaluator equations, and remained inside the
300-second end-to-end boundary. Failed and diagnostic runs are not promoted to
this table.

| Scenario | Revision | Source contingencies | Base (s) | Contingencies (s) | Serialized boundary (s) | Final controller report (s) | Official objective | Official infeasibility | Maximum independent residual | Result |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 006 | `b5f65c9aa0c59bd1c6af47fe3207ae8f11939ff6` | 6,693 | 55.483 | 219.948 | 298.978 | 298.978 | -685,035,133.344363 | 0.0 | 9.8504204e-6 | PASS |
| 010 | `52103f80d256298864e2aee9e52fe017b66ce9b5` | 6,693 | 27.293 | 265.336 | 297.671 | 297.743 | -582,029,928.559992 | 0.0 | 9.7037211e-6 | PASS |

## Scenario 010 certificate

- Local retained run: `runs/c2fen19402_52103f8_20260826/s010_cold_r24_h11`
- Base commitment count: 968
- Generated contingency solutions: 6,693 of 6,693
- Unique exact-evaluation details: 6,694 (base plus every contingency)
- Infeasible labels: none
- Fast-screen fallbacks: zero
- Timing-only heavy lane: 11 of 24 resident screen workers
- Exact validation: 47 completion-order shards; the final 96, 64, 48, 32,
  and 16-contingency shards taper the deadline tail
- Validator implementation: canonical binary numeric parser with line/string
  fallback, while retaining the referenced vendor equations, limits,
  tolerances, objective, and infeasibility decision
- Referenced vendor evaluator SHA-256:
  `3e98ca0e5dda571fc90c2b16bc4000d652f50c02acc0c9a203d99d132aeb603e`
- Normalized case SHA-256:
  `1f0b056260851e105be0c79809e87de9cba4c9981a5feb8c23598deedb89ac7b`
- RAW SHA-256:
  `e6010ea3192069203630bd4b9fedf5e510bc432ac5331dea958c855701d88eec`
- CON SHA-256:
  `9059f78af875ee1b186f64f458504d0490614a7eeb090cba6d5f10c07b4b9179`

Artifact hashes for the retained scenario 010 run:

| Artifact | SHA-256 |
|---|---|
| `run_summary.json` | `76fb40abfefa0ef6e52d8d943780456e8532f24d85c77cec00f30018f601d805` |
| `eval_summary.json` | `44dd433d2272fbf9d47482d21ccd65d783a57f6f179a61f1f6300bb142077a62` |
| `run_status.json` | `3136d3e4cf6e887a8851221dcf875555cc9deb4a29cab7f5e53c5302d5c41ac6` |
| `internal/base.json` | `f21ecfe7308ff9cb777e803a0be4251aefb608f751db025eafcc4386e9103aae` |

The serialized-boundary time is the value preserved in the result artifacts.
The final controller report is the later console measurement after those final
timing records were written; both are below 300 seconds.
