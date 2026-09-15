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
    for option in ("case-json", "base-json", "fast-screen-json", "executable", "output"):
        parser.add_argument("--" + option, type=Path, required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--seconds", type=float, default=45)
    parser.add_argument("--linear-seed-solver", choices=("simplex", "ipm"))
    args = parser.parse_args()
    if not re.fullmatch(r"CTG_\d{6}", args.label) or not 0 < args.seconds <= 300:
        raise ValueError("Invalid diagnostic label or time budget")
    for path in (args.case_json, args.base_json, args.fast_screen_json, args.executable, args.output):
        reject_onedrive(path)
    if args.output.exists():
        raise ValueError("Diagnostic output exists; refusing duplicate trace")
    args.output.mkdir(parents=True)
    command = cpp_command(args.executable, "Ubuntu-24.04", [
        "contingency-worker", to_wsl(args.case_json), to_wsl(args.base_json),
        "0", "fast-pf", "linearized"], args.seconds)
    command.insert(command.index("timeout"), "GRAVITYX_REPAIR_LOG=1")
    command.insert(command.index("timeout"), "GRAVITYX_HIGHS_LOG=1")
    if args.linear_seed_solver:
        command.insert(command.index("timeout"), "GRAVITYX_LINEAR_SEED_SOLVER=" + args.linear_seed_solver)
    task = {"label": args.label, "output_path": to_wsl(args.output / "contingency.json"),
            "fast_screen_path": to_wsl(args.fast_screen_json)}
    manifest = {"purpose": "single-contingency diagnostic, not an official scenario run",
                "label": args.label, "seconds_limit": args.seconds,
                "input_sha256": {str(path): sha256(path) for path in
                    (args.case_json, args.base_json, args.fast_screen_json, args.executable)},
                "command": command, "task": task}
    write_json(args.output / "manifest.json", manifest)
    start = time.perf_counter()
    with (args.output / "trace.log").open("w") as log:
        run = subprocess.run(command, input=json.dumps(task) + '\n{"stop":true}\n',
                             text=True, stdout=log, stderr=subprocess.STDOUT)
    result = {"returncode": run.returncode, "wall_seconds": time.perf_counter() - start,
              "timed_out": run.returncode == 124, "official_scenario_pass": False}
    write_json(args.output / "diagnostic_status.json", result)
    print(json.dumps(result))
    return 0 if run.returncode in (0, 124) else run.returncode


if __name__ == "__main__":
    raise SystemExit(main())
