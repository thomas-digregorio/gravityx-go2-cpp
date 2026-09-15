"""Bounded single-contingency trace; never an official cold scenario run."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import subprocess
import time

from run_experiment import cpp_command, reject_onedrive, sha256, to_wsl, write_json


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ("case-json", "base-json", "executable", "output"):
        parser.add_argument("--" + option, type=Path, required=True)
    parser.add_argument("--fast-screen-json", type=Path)
    parser.add_argument("--phase", choices=("corrective", "fast-screen", "direct-linearized"),
                        default="corrective")
    parser.add_argument("--label", required=True)
    parser.add_argument("--seconds", type=float, default=45)
    parser.add_argument("--linear-seed-solver", choices=("simplex", "ipm"))
    parser.add_argument("--elastic-phase-one-start", action="store_true")
    parser.add_argument("--simplex-strategy", type=int, choices=(1, 4))
    parser.add_argument("--fast-screen-seconds", type=float)
    args = parser.parse_args()
    if not re.fullmatch(r"CTG_\d{6}", args.label) or not 0 < args.seconds <= 300:
        raise ValueError("Invalid diagnostic label or time budget")
    if (args.phase == "corrective") != (args.fast_screen_json is not None):
        raise ValueError("Only corrective diagnostics require a saved fast-screen state")
    if args.fast_screen_seconds is not None and (
            args.phase != "fast-screen" or
            not 0 < args.fast_screen_seconds < args.seconds):
        raise ValueError("The fast-screen handoff budget must be inside its process limit")
    input_paths = [args.case_json, args.base_json, args.executable]
    if args.fast_screen_json is not None:
        input_paths.append(args.fast_screen_json)
    for path in [*input_paths, args.output]:
        reject_onedrive(path)
    if args.output.exists():
        raise ValueError("Diagnostic output exists; refusing duplicate trace")
    args.output.mkdir(parents=True)
    phase_arguments = {"corrective": ["fast-pf", "linearized"],
                       "fast-screen": ["fast-pf", "fast-only"],
                       "direct-linearized": ["linearized", "linearized-only"]}[args.phase]
    command = cpp_command(args.executable, "Ubuntu-24.04", [
        "contingency-worker", to_wsl(args.case_json), to_wsl(args.base_json),
        "0", *phase_arguments], args.seconds)
    command.insert(command.index("timeout"), "GRAVITYX_REPAIR_LOG=1")
    command.insert(command.index("timeout"), "GRAVITYX_HIGHS_LOG=1")
    command.insert(command.index("timeout"), "GRAVITYX_FAST_PF_DIAGNOSTICS=1")
    if args.fast_screen_seconds is not None:
        command.insert(command.index("timeout"),
                       "GRAVITYX_FAST_PF_SCREEN_SECONDS=" + str(args.fast_screen_seconds))
    if args.linear_seed_solver:
        command.insert(command.index("timeout"), "GRAVITYX_LINEAR_SEED_SOLVER=" + args.linear_seed_solver)
    if args.elastic_phase_one_start:
        command.insert(command.index("timeout"), "GRAVITYX_ELASTIC_PHASE_ONE_START=1")
    if args.simplex_strategy is not None:
        command.insert(command.index("timeout"), "GRAVITYX_LINEAR_SEED_SIMPLEX_STRATEGY=" + str(args.simplex_strategy))
    task = {"label": args.label, "output_path": to_wsl(args.output / "contingency.json")}
    if args.fast_screen_json is not None:
        task["fast_screen_path"] = to_wsl(args.fast_screen_json)
    manifest = {"purpose": "single-contingency diagnostic, not an official scenario run",
                "phase": args.phase,
                "label": args.label, "seconds_limit": args.seconds,
                "input_sha256": {str(path): sha256(path) for path in input_paths},
                "command": command, "task": task}
    write_json(args.output / "manifest.json", manifest)
    start = time.perf_counter()
    with (args.output / "trace.log").open("w") as log:
        # WSL's stdout/stderr relays can use independent offsets when both
        # inherit a disk handle. Merge into a pipe, then use one file writer.
        with subprocess.Popen(command, stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              text=True, bufsize=1) as run:
            assert run.stdin is not None and run.stdout is not None
            run.stdin.write(json.dumps(task) + '\n{"stop":true}\n')
            run.stdin.close()
            for line in run.stdout:
                log.write(line)
                log.flush()
            returncode = run.wait()
    result = {"returncode": returncode, "wall_seconds": time.perf_counter() - start,
              "timed_out": returncode == 124, "official_scenario_pass": False}
    write_json(args.output / "diagnostic_status.json", result)
    print(json.dumps(result))
    return 0 if returncode in (0, 124) else returncode


if __name__ == "__main__":
    raise SystemExit(main())
