# C2FEN04200 final-scenario results — 2026-08-24

## Scope

This report records all seven public GO Challenge 2 Final Event 4
`C2FEN04200` scenarios: `009`, `010`, `014`, `055`, `056`, `057`, and `060`.
The family name is rounded: the evaluated network has 4,224 buses, 399
generators, 1,673 loads, and 4,930 branches (2,605 lines and 2,325
transformers).

Each retained run was cold: it loaded its own normalized scenario and source
prior state without carrying a commitment or operating point from another
scenario. The source contingency file contains 455 contingencies: 12 generator
outages, 391 line outages, and 52 transformer outages. No contingency was
omitted.

The 300-second end-to-end boundary begins at normalized-case loading and ends
after independent validation, complete official evaluation, and result
serialization. All retained runs used eight fast-screen workers and eight
linearized-repair workers on frozen implementation revision
`a42a03ff72e48a2520b2c96f5a08718cb040dad0`.

## Method

The scalable C++ path used the source commitment and the following feasibility
pipeline:

1. construct a base operating point with HiGHS sequential-linearized AC
   subproblems and accept it only after nonlinear independent validation;
2. screen all 455 contingencies with parallel sparse-Newton AC power flow;
3. send every screen failure to an isolated HiGHS sequential-linearized
   corrective solve, with up to three solve-and-validate rounds;
4. accept no state unless the independent nonlinear checker reports a maximum
   residual no greater than `1e-5`;
5. write all 456 solution files and evaluate them with the pinned official GO
   Challenge 2 evaluator.

HiGHS 1.15.1 used presolve and IPM for the linearized subproblems. Each HiGHS
instance was limited to one thread so eight independent repair processes could
run concurrently without oversubscribing the laptop. The retained runs needed
no Ipopt fallback and accepted no solver-failed state under a
"feasible despite nonconvergence" exception.

## Results

| Scenario | Base (s) | Fast screen (s) | Code 2 (s) | Evaluation (s) | End-to-end (s) | Official objective | Official infeasibility | Max independent residual | Fast / linearized | Revision | Result |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| 009 | 31.632 | 14.189 | 187.702 | 25.019 | **244.622** | 1,281,874.038728 | 0.0 | 9.70e-6 | 391 / 64 | `a42a03f` | PASS |
| 010 | 34.848 | 14.424 | 190.380 | 25.791 | **251.301** | 1,550,840.040832 | 0.0 | 9.70e-6 | 391 / 64 | `a42a03f` | PASS |
| 014 | 36.068 | 14.625 | 189.592 | 26.571 | **252.510** | 1,422,612.953720 | 0.0 | 8.22e-6 | 393 / 62 | `a42a03f` | PASS |
| 055 | 35.759 | 15.154 | 192.989 | 26.425 | **255.471** | 2,380,234.479304 | 0.0 | 9.70e-6 | 391 / 64 | `a42a03f` | PASS |
| 056 | 36.528 | 15.009 | 189.213 | 27.277 | **253.314** | 2,630,550.255438 | 0.0 | 9.70e-6 | 391 / 64 | `a42a03f` | PASS |
| 057 | 36.859 | 14.661 | 192.651 | 26.725 | **256.517** | 2,664,784.345189 | 0.0 | 9.71e-6 | 392 / 63 | `a42a03f` | PASS |
| 060 | 36.696 | 15.025 | 193.950 | 27.286 | **258.220** | 2,397,290.764803 | 0.0 | 9.68e-6 | 391 / 64 | `a42a03f` | PASS |

`Fast screen` is a measured subset of `Code 2`; it must not be added to the
other timing columns. `Fast / linearized` partitions all 455 source
contingencies by the method whose resulting state passed independent nonlinear
validation. The recorded end-to-end value is the final persisted timing after
result serialization. The slowest run was scenario 060 at 258.220 seconds,
leaving 41.780 seconds under the requested limit.

## Development failure and correction

The first full scenario-009 attempt on revision `aa70b36` was not retained. It
failed on `CTG_000371` when the first HiGHS repair LP reached its 10-second
local safety limit while eight unconstrained HiGHS processes competed for CPU
threads. This was a resource-oversubscription timeout, not a mathematical
infeasibility result.

Revision `a42a03f` limited each HiGHS instance to one thread and increased the
per-subproblem safety limit to 30 seconds while retaining the same mathematical
model, validation tolerance, source data, and global 300-second boundary. A
targeted component check passed, and the cold scenario-009 replacement then
completed in 244.622 seconds. All seven reported runs use that corrected frozen
revision.

## Acceptance evidence

For each of the seven retained directories:

- `run_status.json` reports `success=true`, stage `complete`, all 455
  contingencies completed, and the 300-second end-to-end gate satisfied;
- exactly 456 solution files and 456 official evaluator-detail files exist;
- `eval_summary.json` reports `solutions_exist=true` and aggregate
  `infeas=0.0`;
- all 456 entries in `infeas_all_cases` are false;
- the independently measured maximum residual is at or below `1e-5`;
- all 455 contingencies are accounted for by the fast and linearized accepted
  counts, with zero Ipopt fallbacks and zero nonconverged acceptances.

The common source-file hashes are:

- `case.raw`: `835acc6d3f559e81d1819b216292cd2362073a321eaf38a28508b235b25ce7e1`
- `case.con`: `755f20d30b22e053d2ce58a2322719189c05834e144e9b9fc87095ddfe3b10fe`

Scenario-specific hashes are:

| Scenario | `case.json` SHA-256 | Normalized-case SHA-256 |
|---:|---|---|
| 009 | `41b3ec7e8e2bbfba940580574bc659013bc54e6bd5d9166d91508473f36fec70` | `18cf4e358e570c263cef5ca97659e4976d718685332fb5a5473769e4d0fea865` |
| 010 | `863b369d576bf607de1c9f5ef48a19a1c5a64889956fe567a9520dec4eeff369` | `b62a68eee02e6a414e47e74b2710a2cf3b6c7c43419a85ef21bcef935d9f4633` |
| 014 | `04076602930ab238f965d9229032c3efdcb73745e0fc081ff3902bdafce176c3` | `b50fc532b07debea633c937a60a2ecf7ae9a4a99136efc5333283b3c61f1c76a` |
| 055 | `e41abb0d6878a601e581d5bd0bba1865584a910123e39a6aba09e6126690d1bb` | `f5660d703e61452a341d995de4bfde4b4b551da4081c99116d224f72d242cf33` |
| 056 | `a5c7f3a807571cd1a95fc3aa3b0909105c1e9ac4cc472bfe9594a16b14a17f9a` | `eafb7e46d4fb782a0383c6f778037d7502c7f9180a92cb72f9da363fd5881bfc` |
| 057 | `d417e2d6c5cfc9f84aaca3adbc61a48e297077a2bcc4a208c4534f5f5a07168c` | `4ad2703246f419eb61f70516b6e4051e23e8e948137caac3f340c7bdfb9c9c52` |
| 060 | `3c262f835927a94bf982872b095053182a28b6a377102488453bb827bf0d6a2f` | `6ff5b0d8109c30dcebd27348de386aaca99b315b7f5ada64fcaf82e31344a738` |

The ignored local run directories retain the full run status, timing summary,
worker logs, independent validation records, official evaluator output, and
all base and contingency solution files. Raw competition inputs and run
artifacts are not committed.
