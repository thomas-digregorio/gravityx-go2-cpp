# C2FEN17700 all-scenario results — 2026-08-25

## Scope

This report records one retained cold run of every public GO Challenge 2 Final
Event 4 Division 1 `C2FEN17700` scenario: `019`, `020`, `021`, `089`, `094`,
`106`, `107`, and `115`. Although the family name is `17700`, the parsed
network contains 16,789 buses.

Every scenario defines the same 238 source contingencies: 15 generator
outages and 223 branch outages. No source contingency was omitted. Across the
eight runs, the retained evidence therefore covers 1,904 contingency cases
and eight base cases.

All retained runs used frozen revision
`44b07a595ed84edc21c0b2dc41fd6931483e0e0a`. Each scenario started cold from
its own source prior state; no operating point, corrective state, or solver
state was transferred between scenarios. Within one scenario, previously
verified corrective states could be reused only as initial guesses for later
outages. Every reused candidate was rebuilt for the current outage and had to
pass the complete independent nonlinear validator.

## Final large-case method

The final C++ path used the source commitment and the following feasibility
pipeline:

1. construct and independently validate a base AC operating point;
2. screen all 238 source contingencies with 16 parallel sparse-Newton AC power
   flow workers;
3. defer failed screens until screening is complete, then process the
   corrective queue on one resident worker;
4. initialize corrective attempts from a scenario-local bank of up to 16
   previously verified states;
5. use a Phase-I HiGHS model that starts with AC balance feasibility and adds
   only branch-security rows actually violated by the current nonlinear
   candidate;
6. rebuild each candidate with the exact nonlinear equations and accept it
   only if the independent checker passes all source bounds, AC balance,
   outage status, ramps, and applicable limits;
7. write the base and all contingency solutions and run the pinned official GO
   Challenge 2 evaluator over the complete label set.

Dynamic branch-row generation avoids repeatedly factoring the roughly
126,000-row full corrective LP while preserving the exact acceptance gate.
The scenario-local state bank is an initialization device, not a feasibility
shortcut. No Ipopt fallback or feasible-but-nonconverged acceptance was used
in any retained run.

The identical retained configuration used one corrective worker, 16 fast
screen workers, one post-screen corrective worker, a 300-second base limit, a
470-second contingency-solve limit, a 300-second Code 1 limit, a 476-second
Code 2 limit (`2 seconds x 238 contingencies`), a 3,600-second end-to-end
limit, and a 300-second evaluation reserve. Longest-first scheduling,
validated source bases, two-stage fast screening, deferred fallback, and the
linearized corrective fallback were enabled.

## Results

| Scenario | Contingencies | Base (s) | Code 2 (s) | Evaluation (s) | End-to-end (s) | Official objective | Official infeasibility | Max independent residual | Screen pass / corrective queue | Result |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 019 | 238 | 4.003 | 189.164 | 21.860 | **215.950** | 6,725,801.100325 | 0.0 | 1.14e-13 | 220 / 18 | PASS |
| 020 | 238 | 3.710 | 190.616 | 21.682 | **216.958** | 6,824,396.741104 | 0.0 | 1.14e-13 | 220 / 18 | PASS |
| 021 | 238 | 3.753 | 190.883 | 21.617 | **217.170** | 6,709,423.548195 | 0.0 | 1.14e-13 | 220 / 18 | PASS |
| 089 | 238 | 3.846 | 128.971 | 21.873 | **155.620** | 6,647,679.269698 | 0.0 | 5.68e-14 | 219 / 19 | PASS |
| 094 | 238 | 14.024 | 119.379 | 24.602 | **158.912** | -4,263,601.238114 | 0.0 | 1.14e-13 | 225 / 13 | PASS |
| 106 | 238 | 4.362 | 214.761 | 21.807 | **241.843** | 94,215.919863 | 0.0 | 1.14e-13 | 207 / 31 | PASS |
| 107 | 238 | 4.381 | 216.415 | 21.689 | **243.389** | -653,282.857771 | 0.0 | 1.14e-13 | 207 / 31 | PASS |
| 115 | 238 | 14.683 | 90.281 | 22.002 | **127.868** | 387,128.584813 | 0.0 | 5.68e-14 | 229 / 9 | PASS |

`Screen pass / corrective queue` partitions all source contingencies into
states accepted directly by the independent checker after fast screening and
states sent to the corrective pipeline. It does not imply that the screen
itself was the acceptance authority.

Stage times do not sum exactly to the end-to-end value because normalized-case
loading, orchestration, solution writing, and final serialization occur
outside the individual stage counters. The table uses the persisted final
`run_status.json` timing checkpoint; the process-return messages were about
0.2–0.3 seconds more conservative after the final status rewrite.

All eight Code 2 stages completed below 476 seconds. The slowest was scenario
107 at 216.415 seconds. The slowest persisted end-to-end run was also scenario
107 at 243.389 seconds, or approximately 4.06 minutes.

## Acceptance evidence

A final automated audit returned `AUDIT_PASS` and established, for every
scenario:

- `run_status.json` reports `success=true`, stage `complete`, and satisfied
  Code 1, Code 2, and end-to-end gates;
- exactly 238 contingencies were screened and completed;
- all 238 internal contingency records report success;
- exactly 239 solution files and 239 official evaluator details exist: one
  base case plus all 238 source contingencies;
- the official certificate reports a complete label set and no infeasible
  label;
- `eval_summary.json` reports `solutions_exist=true`, 238 contingencies, and
  aggregate `infeas=0.0`;
- the largest independent exact residual is
  `1.1368683772161603e-13`, well below the `1e-5` gate;
- every artifact records the same frozen Git revision.

The current component suite also passed 4/4 CTest tests and 13/13 Python unit
tests before the retained runs.

The official objective is the evaluator's aggregate score for the submitted
feasible base and corrective states. Negative values are valid under the GO
cost/benefit convention. These runs establish complete feasibility and
evaluation, not global economic optimality or a leaderboard comparison.

## Input identity

| Scenario | `case.raw` SHA-256 | `case.con` SHA-256 | `case.json` SHA-256 | Normalized-case SHA-256 |
|---:|---|---|---|---|
| 019 | `5c714d85a8f494c39ba736b75614464bebdfbfafd7693aa3051ed4f37bf2b4c6` | `fd79ef00b41afda792a48f2415f5bf49ce4ea1f2b6e58182665293a5814cf804` | `a2b47857c80414b5a2800d7ffd9f6dedfd417a00bab9df7de4e9ca9ec330f6e1` | `4b857f37b975790d84f41e4b9522bec9d56d4dffbcce939574d2e67c287ddd7e` |
| 020 | `5c714d85a8f494c39ba736b75614464bebdfbfafd7693aa3051ed4f37bf2b4c6` | `fd79ef00b41afda792a48f2415f5bf49ce4ea1f2b6e58182665293a5814cf804` | `107a759b925e63baedf8559cfd455206adc39ffb4ac6ae2b198d54ce8ab3926d` | `57f27c400a006f216c12941c8b709a8dced81e33285f0d2e1b950d9e309edb70` |
| 021 | `5c714d85a8f494c39ba736b75614464bebdfbfafd7693aa3051ed4f37bf2b4c6` | `fd79ef00b41afda792a48f2415f5bf49ce4ea1f2b6e58182665293a5814cf804` | `9e30ca86759036a6da9d607b42d1427c7e35d43bf75f9d03ecce81d0711bc360` | `86e9720bba17cc8410d6b5245d46640b7de6fe176dde90d6a27b2768494425c5` |
| 089 | `b986cf2fca0a17ae31693c6eaac78ac87d3b2f8716f2a6a712de1e813707dfdd` | `fd79ef00b41afda792a48f2415f5bf49ce4ea1f2b6e58182665293a5814cf804` | `88b226b77bc20f951a4d405f6355f4d98fc21aa979a931764d426b2e9fe793ea` | `aebfeb43bdb8b31a63924c6907db84742827ccf4c989cebc8c79ac1ea182e008` |
| 094 | `26929e36627eb535cedd34f8de80641eb30d01644c1a78dd9664d7d9de6402cf` | `fd79ef00b41afda792a48f2415f5bf49ce4ea1f2b6e58182665293a5814cf804` | `cd447da4cfcc0bf95ea000d903c250e157b8ac9374f3a1b87d30872e64a61167` | `e2363c546f2c62cd49883c70288556d14dc0be617ca6754ee965be146051c064` |
| 106 | `32960ba4cce17bac4450273be91252b882567d012307a83ed417fc4854cdb23a` | `fd79ef00b41afda792a48f2415f5bf49ce4ea1f2b6e58182665293a5814cf804` | `1de5f131f72b92f8ca605910fcf7a08daaf4a212d81baf7d16ecd830d28cdbeb` | `d3e65f88b0f3c1f01b97abe0fc8d67809891c8e0ec5079ae47083f1763d1a81a` |
| 107 | `32960ba4cce17bac4450273be91252b882567d012307a83ed417fc4854cdb23a` | `fd79ef00b41afda792a48f2415f5bf49ce4ea1f2b6e58182665293a5814cf804` | `4e3cc02cacc70861d32d2406c0f70eaa0ba4794f74012abb589453e8d2285ef9` | `56ccf2e395fdb3fc95ca83ed94898424e77493d6db809a842c9044b2733b3a44` |
| 115 | `a4fc6e5533e83415be55a72eea87b4d0fa5b50d0c0d587cdfb26020e112891c4` | `fd79ef00b41afda792a48f2415f5bf49ce4ea1f2b6e58182665293a5814cf804` | `5a647245aabca25f30b693ba7403ec3b32cae92558c222e78991e4338c82b2b9` | `544b88a5e689ebcddd23a7bbbf0a273374b42ae1f1ea54ad99573090ab6b9534` |

Raw competition inputs and the complete run directories remain ignored by
Git. They are retained locally under
`C:\Users\thoma\Documents\gravityx-go2-cpp\runs\c2fen17700_44b07a5_w1_all_20260824`,
off OneDrive. The tracked report records the immutable hashes, configuration,
results, and acceptance gates without publishing the competition archive or
large generated artifacts.
