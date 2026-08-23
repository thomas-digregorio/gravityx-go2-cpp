#!/usr/bin/env python3
"""Run one cold Gravity C++ GO2 experiment and the official evaluator.

All optimization subprocesses use the pinned Gravity C++ executable.  Python
only orchestrates isolated contingency processes, writes the official text
format, records provenance, and invokes the official evaluator.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
from typing import Any


REPO = Path(__file__).resolve().parents[1]
DEFAULT_EXE = REPO / ".build" / "gravityx_go2"
DEFAULT_PYTHON = Path(
    r"C:\Users\thoma\Documents\goc2-ac-score-check\evaluator-venv\Scripts\python.exe"
)
DEFAULT_EVALUATOR = Path(
    r"C:\Users\thoma\Documents\goc2-dc-scopf-cpu\vendor\C2DataUtilities\data_utilities\evaluation.py"
)
WSL_LIBRARY_PATH = "/home/thomasdigregorio/.local/share/gravityx-go2-cpp/env/lib"


def reject_onedrive(path: Path) -> None:
    if "onedrive" in str(path.resolve()).lower():
        raise ValueError(f"refusing OneDrive path: {path}")


def to_wsl(path: Path) -> str:
    path = path.resolve()
    reject_onedrive(path)
    drive, tail = os.path.splitdrive(str(path))
    if not drive:
        return str(path).replace("\\", "/")
    return f"/mnt/{drive[0].lower()}/{tail.lstrip('\\/').replace('\\', '/')}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as stream:
        return json.load(stream)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
    temporary.replace(path)


def ordered(case: dict[str, Any], group: str) -> list[dict[str, Any]]:
    return sorted(case[group].values(), key=lambda item: int(item["index"]))


def number(value: Any) -> str:
    if isinstance(value, int):
        return str(value)
    return format(float(value), ".17g")


def write_solution(
    destination: Path,
    case: dict[str, Any],
    state: dict[str, Any],
    commitment: list[int],
    contingency: dict[str, Any] | None = None,
) -> None:
    outage_type = contingency["type"] if contingency else None
    outage_index = int(contingency["idx"]) if contingency else None
    buses = ordered(case, "bus")
    loads = ordered(case, "load")
    generators = ordered(case, "gen")
    branches = ordered(case, "branch")
    shunts = ordered(case, "shunt")

    lines: list[str] = ["--bus section", "i, v, theta"]
    for position, bus in enumerate(buses):
        if bus.get("present", True):
            lines.append(
                f"{bus['index']}, {number(state['vm'][position])}, "
                f"{number(state['va'][position])}"
            )

    lines.extend(["--load section", "i, id, t"])
    for position, load in enumerate(loads):
        if load.get("present", True):
            source = load["source_id"]
            lines.append(
                f"{source[1]}, {source[2]}, {number(state['demand_factor'][position])}"
            )

    lines.extend(["--generator section", "i, id, p, q, x"])
    for position, generator in enumerate(generators):
        if not generator.get("present", True):
            continue
        if outage_type == "gen" and int(generator["index"]) == outage_index:
            continue
        source = generator["source_id"]
        lines.append(
            f"{source[1]}, {source[2]}, {number(state['pg'][position])}, "
            f"{number(state['qg'][position])}, {int(commitment[position])}"
        )

    lines.extend(["--line section", "iorig, idest, id, x"])
    for branch in branches:
        if not branch.get("present", True) or branch["transformer"]:
            continue
        if outage_type == "branch" and int(branch["index"]) == outage_index:
            continue
        source = branch["source_id"]
        lines.append(f"{source[1]}, {source[2]}, {source[3]}, {int(branch['br_status'] != 0)}")

    lines.extend(["--transformer section", "iorig, idest, id, x, xst"])
    for branch in branches:
        if not branch.get("present", True) or not branch["transformer"]:
            continue
        if outage_type == "branch" and int(branch["index"]) == outage_index:
            continue
        source = branch["source_id"]
        control_mode = int(branch.get("control_mode", 0))
        step = 0
        if control_mode in (1, -1):
            step = int(branch["tm_step"])
        elif control_mode in (3, -3):
            step = int(branch["ta_step"])
        lines.append(
            f"{source[1]}, {source[2]}, {source[4]}, "
            f"{int(branch['br_status'] != 0)}, {step}"
        )

    lines.extend(
        [
            "--switched shunt section",
            "i, xst1, xst2, xst3, xst4, xst5, xst6, xst7, xst8",
        ]
    )
    for shunt in shunts:
        if shunt.get("present", True) and shunt.get("dispatchable", False):
            source = shunt["source_id"]
            values = [str(source[1]), *(str(int(value)) for value in shunt.get("xst", []))]
            lines.append(", ".join(values))

    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def safe_label(label: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]", "_", label)
    if not value or value in (".", ".."):
        raise ValueError(f"unsafe contingency label: {label!r}")
    return value


def run_cpp(
    executable: Path,
    distro: str,
    arguments: list[str],
    log_path: Path,
    timeout: float,
) -> subprocess.CompletedProcess[str]:
    command = [
        "wsl",
        "-d",
        distro,
        "--",
        "env",
        f"LD_LIBRARY_PATH={WSL_LIBRARY_PATH}",
        "OMP_NUM_THREADS=1",
        "OPENBLAS_NUM_THREADS=1",
        to_wsl(executable),
        *arguments,
    ]
    started = time.perf_counter()
    try:
        completed = subprocess.run(
            command,
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired as error:
        output = error.stdout or ""
        if isinstance(output, bytes):
            output = output.decode(errors="replace")
        log_path.write_text(output + "\nTIMEOUT\n", encoding="utf-8")
        raise RuntimeError(f"C++ subprocess timed out after {timeout}s: {arguments}") from error
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(completed.stdout, encoding="utf-8")
    completed.wall_seconds = time.perf_counter() - started  # type: ignore[attr-defined]
    return completed


def contingency_records(case: dict[str, Any]) -> list[dict[str, Any]]:
    records = [dict(item) for item in case["gen_contingencies"]]
    records.extend(dict(item) for item in case["branch_contingencies"])
    return sorted(records, key=lambda item: item["label"])


def git_revision() -> str | None:
    completed = subprocess.run(
        ["git", "-C", str(REPO), "rev-parse", "HEAD"],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    return completed.stdout.strip() if completed.returncode == 0 else None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-json", type=Path, required=True)
    parser.add_argument("--case-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--executable", type=Path, default=DEFAULT_EXE)
    parser.add_argument("--workers", type=int, default=22)
    parser.add_argument("--distro", default="Ubuntu-24.04")
    parser.add_argument("--base-timeout", type=float, default=900.0)
    parser.add_argument("--contingency-timeout", type=float, default=300.0)
    parser.add_argument("--python", type=Path, default=DEFAULT_PYTHON)
    parser.add_argument("--evaluator", type=Path, default=DEFAULT_EVALUATOR)
    parser.add_argument("--skip-evaluation", action="store_true")
    args = parser.parse_args()

    for path in (args.case_json, args.case_dir, args.output_dir, args.executable):
        reject_onedrive(path)
    if args.workers < 1:
        raise ValueError("workers must be positive")
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise ValueError(f"cold-run output directory is not empty: {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    case = read_json(args.case_json)
    started_utc = dt.datetime.now(dt.timezone.utc).isoformat()
    wall_start = time.perf_counter()
    internal = args.output_dir / "internal"
    base_json = internal / "base.json"

    base_start = time.perf_counter()
    base_process = run_cpp(
        args.executable,
        args.distro,
        ["run-ibr-json", to_wsl(args.case_json), to_wsl(base_json), "0"],
        internal / "base.console.log",
        args.base_timeout,
    )
    base_wall = time.perf_counter() - base_start
    if base_process.returncode != 0:
        raise RuntimeError("cold C++ base/IBR solve failed; see internal/base.console.log")
    base = read_json(base_json)
    if not base.get("success", False):
        raise RuntimeError("cold C++ base/IBR result was not successful")
    commitment = [int(value) for value in base["commitment"]]
    base_state = base["selected_state"]
    write_solution(args.output_dir / "solution_BASECASE.txt", case, base_state, commitment)

    contingencies = contingency_records(case)
    contingency_start = time.perf_counter()

    def solve_one(item: dict[str, Any]) -> dict[str, Any]:
        label = str(item["label"])
        stem = safe_label(label)
        result_path = internal / "contingencies" / f"{stem}.json"
        log_path = internal / "logs" / f"{stem}.log"
        process = run_cpp(
            args.executable,
            args.distro,
            [
                "solve-contingency",
                to_wsl(args.case_json),
                to_wsl(base_json),
                label,
                to_wsl(result_path),
                "0",
            ],
            log_path,
            args.contingency_timeout,
        )
        if process.returncode != 0:
            raise RuntimeError(f"contingency {label} failed; see {log_path}")
        result = read_json(result_path)
        if not result.get("success", False):
            raise RuntimeError(f"contingency {label} did not pass independent validation")
        write_solution(
            args.output_dir / f"solution_{label}.txt",
            case,
            result["solve"]["state"],
            commitment,
            item,
        )
        return {
            "label": label,
            "type": item["type"],
            "source_index": int(item["idx"]),
            "process_wall_seconds": process.wall_seconds,  # type: ignore[attr-defined]
            "solver_wall_seconds": result["solve"]["wall_seconds"],
            "objective": result["solve"]["objective"],
            "max_residual": result["validation"]["max_residual"],
        }

    records: list[dict[str, Any]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(solve_one, item): item for item in contingencies}
        for future in concurrent.futures.as_completed(futures):
            records.append(future.result())
            print(f"completed {len(records)}/{len(contingencies)}: {records[-1]['label']}", flush=True)
    records.sort(key=lambda item: item["label"])
    contingency_wall = time.perf_counter() - contingency_start

    evaluation_wall = 0.0
    evaluation_summary: dict[str, Any] | None = None
    if not args.skip_evaluation:
        for path in (args.python, args.evaluator):
            reject_onedrive(path)
        evaluation_start = time.perf_counter()
        completed = subprocess.run(
            [str(args.python), str(args.evaluator), "1", str(args.case_dir), str(args.output_dir)],
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        evaluation_wall = time.perf_counter() - evaluation_start
        (internal / "evaluation.console.log").write_text(completed.stdout, encoding="utf-8")
        if completed.returncode != 0:
            raise RuntimeError("official evaluator failed; see internal/evaluation.console.log")
        evaluation_summary = read_json(args.output_dir / "eval_summary.json")

    total_wall = time.perf_counter() - wall_start
    max_residual = max((item["max_residual"] for item in records), default=0.0)
    summary = {
        "method": "Gravity C++ continuous AC-UC relaxation plus deterministic iterative batch rounding",
        "exact_unpublished_gravityx_binary": False,
        "framework_faithful": True,
        "started_at_utc": started_utc,
        "finished_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "git_revision": git_revision(),
        "input_hashes": {
            "normalized_case_sha256": sha256(args.case_json),
            "raw_sha256": sha256(args.case_dir / "case.raw"),
            "con_sha256": sha256(args.case_dir / "case.con"),
            "json_sha256": sha256(args.case_dir / "case.json"),
        },
        "workers": args.workers,
        "base_process_wall_seconds": base_wall,
        "base_algorithm_wall_seconds": base["wall_seconds"],
        "contingency_parallel_wall_seconds": contingency_wall,
        "contingency_solver_seconds_sum": sum(item["solver_wall_seconds"] for item in records),
        "evaluation_wall_seconds": evaluation_wall,
        "total_wall_seconds": total_wall,
        "contingency_count": len(records),
        "max_independent_contingency_residual": max_residual,
        "base": base,
        "contingencies": records,
        "official_evaluation": evaluation_summary,
    }
    write_json(args.output_dir / "run_summary.json", summary)
    print(
        json.dumps(
            {
                "output_dir": str(args.output_dir),
                "total_wall_seconds": total_wall,
                "base_wall_seconds": base_wall,
                "contingency_wall_seconds": contingency_wall,
                "objective": evaluation_summary.get("obj") if evaluation_summary else None,
                "infeasibility": evaluation_summary.get("infeas") if evaluation_summary else None,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"error: {error}", file=sys.stderr)
        raise
