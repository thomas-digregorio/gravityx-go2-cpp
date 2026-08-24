# GravityX GO Challenge 2 C++ reproduction

This repository is a clean-room, framework-faithful reproduction of the
publicly described GravityX approach for GO Competition Challenge 2. Its core
nonlinear optimization models are constructed with the actual open-source
[Gravity C++](https://github.com/coin-or/Gravity) modeling framework and solved
with open-source Ipopt. The scalable feasibility path additionally uses
open-source HiGHS for sequential-linearized AC subproblems whose returned
states must pass the same independent nonlinear and official checks.

The original competition solver source and its dataset-specific batch ordering
and rounding rules are not public.  Consequently, this project does **not**
claim bit-for-bit identity with the winning submission.  It implements the
published continuous AC-MINLP relaxation and iterative batch-rounding method,
with every reproduction-specific choice recorded in the run manifest.

## Reproducibility boundary

- Initial reproduction target: GO Challenge 2 Final `C2FEN00617`, scenario
  `005`; validated scale suites now include the public Final 617-, 2,020-, and
  4,200-, and 8,300-bus scenario families.
- Input normalization: the pinned official PowerModelsSecurityConstrained GO2
  parser exports a neutral JSON representation of RAW/JSON/CON data.
- Optimization: Gravity C++ variables, objectives, and constraints; Ipopt with
  MUMPS for continuous nonlinear solves.
- Scoring: the pinned official GO Challenge 2 evaluator.
- Working tree and all generated artifacts must remain outside OneDrive.

See `docs/METHOD_SCOPE.md` for the exact fidelity boundary and
`provenance/dependencies.json` for pinned revisions.

The completed cold Final-617 experiment and its apples-to-apples comparison
with the earlier Julia/PowerModels reproduction are recorded in
`results/FINAL_617_GRAVITY_CPP_COMPARISON.md`.

The multi-scenario scale results are recorded in:

- `docs/C2FEN00617_ABLATION_20260823.md`
- `docs/C2FEN02020_FINAL_SCENARIOS_20260824.md`
- `docs/C2FEN04200_FINAL_SCENARIOS_20260824.md`
- `docs/C2FEN08300_FINAL_SCENARIOS_20260824.md`
- `docs/FINAL_EVENT_REPEATABILITY_FC5E43C_20260824.md`

## What is implemented

- The base soft ACOPF, continuous AC unit-commitment relaxation, and every
  corrective contingency NLP are formulated with the pinned Gravity C++
  framework.
- Ipopt/MUMPS solves the continuous nonlinear models.
- HiGHS solves the sequential-linearized AC base and corrective feasibility
  subproblems used by the 4,224-bus scale path; every candidate is accepted
  only after independent nonlinear validation.
- Iterative batch rounding fixes commitment variables in four deterministic
  batches, warm-starting each structurally related relaxation.
- All source generator and branch contingencies are solved independently with
  the selected base commitment fixed and source corrective ramp limits.
- An independent C++ checker validates every physical result, after which the
  official GO Challenge 2 evaluator supplies the reported score.
- The official evaluator can run through Microsoft MPI.  Because its MPI
  branch leaves serial accumulator fields stale, the runner requires the exact
  base-plus-contingency detail-label set, checks every per-case infeasibility,
  recomputes the aggregate, archives the raw vendor summary, and normalizes
  only those bookkeeping fields.

Python is used only for process orchestration, official text-file generation,
and evaluator invocation.  The Julia tool under `tools/` is a development-only
semantic oracle and is not used by the experiment runner.

## Build

From WSL Ubuntu, in a repository path outside OneDrive:

```bash
./scripts/bootstrap_wsl.sh
```

The script checks out the exact Gravity revision recorded in
`provenance/dependencies.json`, creates a local conda-forge Ipopt/MUMPS
environment, builds Gravity and this executable, and runs the component tests.

## Input export

Raw GO2 data are intentionally ignored.  The normalized model is produced by
the pinned official parser:

```powershell
$env:JULIA_DEPOT_PATH = 'C:\Users\thoma\Documents\goc2-ac-score-check\julia-depot'
& 'C:\Users\thoma\Documents\goc2-ac-score-check\julia\julia-1.10.11\bin\julia.exe' `
  --project='C:\Users\thoma\Documents\gravityx-inspired-go2' `
  tools\export_case.jl CASE_DIR .data\final_617_005\model.json C2FEN00617 005
```

## One cold experiment

The registered Final-617 settings are in `config/final_617_005.json`.
`OUTPUT_DIR` must be absent or empty:

```powershell
& 'C:\Users\thoma\Documents\goc2-ac-score-check\evaluator-venv\Scripts\python.exe' `
  scripts\run_experiment.py `
  --case-json .data\final_617_005\model.json `
  --case-dir CASE_DIR `
  --output-dir OUTPUT_DIR `
  --workers 8
```

See `docs/VALIDATION.md` for the independent checks and tiny-case equivalence
evidence.
