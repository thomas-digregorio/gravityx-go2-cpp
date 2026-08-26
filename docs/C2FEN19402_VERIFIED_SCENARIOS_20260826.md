# C2FEN19402 verified final-event scenarios

This note records only cold runs that completed the complete source contingency
set, passed the referenced official evaluator equations, and remained inside the
300-second end-to-end boundary. Failed and diagnostic runs are not promoted to
this table.

| Scenario | Revision | Source contingencies | Base (s) | Contingencies (s) | Serialized boundary (s) | Final controller report (s) | Official objective | Official infeasibility | Maximum independent residual | Result |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 006 | `b5f65c9aa0c59bd1c6af47fe3207ae8f11939ff6` | 6,693 | 55.483 | 219.948 | 298.978 | 298.978 | -685,035,133.344363 | 0.0 | 9.8504204e-6 | PASS |
| 010 | `52103f80d256298864e2aee9e52fe017b66ce9b5` | 6,693 | 27.293 | 265.336 | 297.671 | 297.743 | -582,029,928.559992 | 0.0 | 9.7037211e-6 | PASS |
| 069 | `1dd89725b25d7b5d8d6dc860b775440e5fb9e82c` | 6,620 | 31.087 | 237.030 | 276.073 | 276.812 | -594,811,182.974502 | 0.0 | 9.9198612e-6 | PASS |
| 077 | `0e4808a045d0455244bc554b60c62536af81e69f` | 6,584 | 26.548 | 131.702 | 169.878 | 170.605 | -433,016,702.926414 | 0.0 | 8.1132918e-6 | PASS |

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

## Scenario 069 certificate

- Local retained run:
  `runs/c2fen19402_1dd8972_20260826/s069_cold_r3_corrected`
- Base commitment count: 927 of 927 source generators
- Generated contingency solutions: 6,620 of 6,620
- Unique exact-evaluation details: 6,621 (base plus every contingency)
- Infeasible labels: none
- Fast-screen fallbacks: zero
- Accepted feasible nonconverged states: zero
- Base independent residual: `8.8817842e-16`
- Maximum independent contingency residual: `9.9198612e-6`
- Exact validation: 47 completion-order shards, with 16 persistent evaluator
  processes after the screen released its worker capacity
- Evaluator implementation SHA-256:
  `3a3e465b6cb023866b1288d9088fb74c25cfda9befa05c6ec864020fc40288c8`
- Referenced vendor evaluator SHA-256:
  `3e98ca0e5dda571fc90c2b16bc4000d652f50c02acc0c9a203d99d132aeb603e`
- Normalized case SHA-256:
  `d39ce71d2d7212116f63b38a37b60833e75d862c0849bbba91a300b0c2b0bcf0`
- RAW SHA-256:
  `65a37a85378162457e5f890964e7701ae26e0aace679c0943d511b9fbb318150`
- CON SHA-256:
  `97f9c22e564fadb4037b3b3c499044075d051555b69214deda32ebed0694835a`
- Supplemental JSON SHA-256:
  `1351c33d5efa0a5387b0a14561bb4a1e54cd22611fe34971e078d48dc3012236`

Artifact hashes for the retained scenario 069 run:

| Artifact | SHA-256 |
|---|---|
| `run_summary.json` | `2515142336fa2aaa211d874d40208cc4eee5de1a06a299afa30a4933771ca59b` |
| `eval_summary.json` | `c5fb8ae419db17ecb94d88be632d03e4916d1da134dd7886a6cc9888bb933ba1` |
| `run_status.json` | `5ed9c1a574c5083ab5302ebac3afb127a9720b0cc8141f50aaf3d43f7c18ce52` |
| `internal/base.json` | `a21ccecf2300d9b02bc71208a01c49ebc33347fb5a7610b35180007dad70ad6b` |
| `internal/official_evaluation_certificate.json` | `a1f3f65ab8f93f1d06d876fa1ba742456a9ec2108d57861c2b2decf452054fad` |

The protocol-misconfigured attempt immediately preceding this run produced no
scientific result. Revision `1dd8972` adds a preflight protocol guard, a ready
handshake, evaluator path/hash provenance, and deterministic pipe cleanup.

## Scenario 077 certificate

- Local retained run:
  `runs/c2fen19402_0e4808a_20260826/s077_cold_r2_replacement`
- Base commitment count: 921 of 921 source generators
- Generated contingency solutions: 6,584 of 6,584
- Unique exact-evaluation details: 6,585 (base plus every contingency)
- Infeasible labels: none
- Fast-screen fallbacks: zero
- Accepted feasible nonconverged states: zero
- Base independent residual: `1.7763568e-15`
- Maximum independent contingency residual: `8.1132918e-6`
- Evaluator implementation SHA-256:
  `3a3e465b6cb023866b1288d9088fb74c25cfda9befa05c6ec864020fc40288c8`
- Referenced vendor evaluator SHA-256:
  `3e98ca0e5dda571fc90c2b16bc4000d652f50c02acc0c9a203d99d132aeb603e`
- Normalized case SHA-256:
  `1764770ff954ed69d558147d847448e0914895386e644394815f6a098993dc1f`
- RAW SHA-256:
  `53ced48bec2e28f31202fa47d8e67887c7d789f79d95560731929a18b7d9bf2e`
- CON SHA-256:
  `5759d66c310d59a2b9690a3d67bcabf42eac7c3425d3b0a3becf4d1bb8b3d1ed`
- Supplemental JSON SHA-256:
  `5092b564d4207d49dbff4aae7db7a27b452b512b510ec04d9f79b8cb1b064f79`

Artifact hashes for the retained scenario 077 run:

| Artifact | SHA-256 |
|---|---|
| `run_summary.json` | `54ceca9d8553c5e1d6e7e898d39ad9db766072ba03b01ed7d14122e384de9c32` |
| `eval_summary.json` | `bd15f11225b134678c3ad8a8035b43bcdb771e851bc5efd5f3492e2ee2043b10` |
| `run_status.json` | `1354c1a7ebdfbd82b4c5ce1b99c1626579b528cbd46247ead4fb0738922015b0` |
| `internal/base.json` | `49abf3eba2a7e7aaaec97b21568c5da4be8ddd7cbea27beedbf9c7f138a9b967` |
| `internal/official_evaluation_certificate.json` | `7836c6072032d41ea053c29fa6db2fe9f88d5158e7d3544518955edeaacc273e` |

The initial scenario 077 attempt reached no contingency solve because two
lightweight LP columns had empty voltage intervals. Revision `0e4808a`
projects only the LP linearization voltages into source `VMIN/VMAX`, records
the projection, and leaves the complete nonlinear acceptance checks unchanged.
