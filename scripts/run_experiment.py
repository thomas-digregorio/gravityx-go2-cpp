#!/usr/bin/env python3
"""Run one cold Gravity C++ GO2 experiment and the official evaluator.

All optimization subprocesses use the pinned Gravity C++ executable.  Python
only orchestrates isolated contingency processes, writes the official text
format, records provenance, and invokes the official evaluator.
"""

from __future__ import annotations

import argparse
import csv
import concurrent.futures
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import queue
import re
import shutil
import subprocess
import sys
import threading
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
FINALIZATION_RESERVE_SECONDS = 1.0


class CompetitionTimeout(RuntimeError):
    """Raised when a GO Competition stage exhausts its wall-clock allowance."""


def code2_time_limit(contingency_count: int, seconds_per_contingency: float) -> float:
    if contingency_count < 0:
        raise ValueError("contingency count cannot be negative")
    if seconds_per_contingency <= 0:
        raise ValueError("seconds per contingency must be positive")
    return contingency_count * seconds_per_contingency


def effective_process_timeout(
    per_process_timeout: float,
    stage_deadline: float,
    now: float | None = None,
) -> float:
    if per_process_timeout <= 0:
        raise ValueError("per-process timeout must be positive")
    remaining = stage_deadline - (time.perf_counter() if now is None else now)
    if remaining <= 0:
        raise CompetitionTimeout("competition stage deadline has expired")
    return min(per_process_timeout, remaining)


def streamed_queue_get(
    task_queue: queue.Queue[dict[str, Any]],
    screening_finished: threading.Event,
    abort: threading.Event,
    poll_seconds: float = 0.1,
    profiled_queue: queue.PriorityQueue[
        tuple[float, int, str, dict[str, Any]]
    ] | None = None,
) -> dict[str, Any] | None:
    """Wait for streamed fallback work until screening is complete or aborted.

    Measured fallback durations define one global longest-first priority queue.
    Any available worker can claim its next item, so neither stale worker
    affinity nor fast-screen completion order can strand a long corrective LP.
    """
    while not abort.is_set():
        if profiled_queue is not None:
            try:
                return profiled_queue.get_nowait()[3]
            except queue.Empty:
                pass
        try:
            return task_queue.get(timeout=poll_seconds)
        except queue.Empty:
            if screening_finished.is_set():
                return None
    return None


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
    temporary = path.with_name(
        f"{path.name}.{os.getpid()}.{threading.get_ident()}.tmp"
    )
    with temporary.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
    delays = (0.01, 0.02, 0.04, 0.08, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10)
    for attempt, delay in enumerate(delays):
        try:
            temporary.replace(path)
            return
        except PermissionError:
            if attempt == len(delays) - 1:
                raise
            time.sleep(delay)


def validate_and_normalize_evaluation_details(
    output_dir: Path,
    internal_dir: Path,
    expected_contingency_labels: set[str],
    summary: dict[str, Any],
    parallel_processes: int,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Require exhaustive official details and repair stale MPI bookkeeping."""
    expected_labels = {"BASECASE", *expected_contingency_labels}
    detail_paths = {
        path.stem.removeprefix("eval_detail_"): path
        for path in output_dir.glob("eval_detail_*.json")
    }
    actual_labels = set(detail_paths)
    missing = sorted(expected_labels - actual_labels)
    extra = sorted(actual_labels - expected_labels)
    if missing or extra:
        raise RuntimeError(
            "official evaluator detail-set mismatch: "
            f"missing={missing}, extra={extra}"
        )

    objectives: dict[str, float] = {}
    infeasibilities: dict[str, bool] = {}
    for label in sorted(expected_labels):
        detail = read_json(detail_paths[label])
        try:
            objective = float(detail["obj"]["val"])
            infeasible = bool(detail["infeas"]["val"])
        except (KeyError, TypeError, ValueError) as error:
            raise RuntimeError(
                f"invalid official evaluator detail for {label}"
            ) from error
        if not math.isfinite(objective):
            raise RuntimeError(
                f"non-finite official objective for {label}: {objective}"
            )
        objectives[label] = objective
        infeasibilities[label] = infeasible

    contingency_count = len(expected_contingency_labels)
    if int(summary.get("num_ctg", -1)) != contingency_count:
        raise RuntimeError(
            "official evaluator contingency count mismatch: "
            f"summary={summary.get('num_ctg')}, expected={contingency_count}"
        )
    expected_objective = objectives["BASECASE"]
    if contingency_count:
        expected_objective += math.fsum(
            objectives[label] for label in sorted(expected_contingency_labels)
        ) / contingency_count
    serial_cumulative_objective = objectives["BASECASE"]
    for label in sorted(expected_contingency_labels):
        serial_cumulative_objective += objectives[label] / contingency_count
    reported_objective = float(summary.get("obj", math.nan))
    if not math.isclose(
        reported_objective,
        expected_objective,
        rel_tol=1e-12,
        abs_tol=1e-6,
    ):
        raise RuntimeError(
            "official evaluator objective disagrees with its detail files: "
            f"summary={reported_objective}, details={expected_objective}"
        )

    infeasible_labels = sorted(
        label for label, infeasible in infeasibilities.items() if infeasible
    )
    expected_infeasibility = float(len(infeasible_labels))
    reported_infeasibility = float(summary.get("infeas", math.nan))
    if reported_infeasibility != expected_infeasibility:
        raise RuntimeError(
            "official evaluator infeasibility disagrees with its detail files: "
            f"summary={reported_infeasibility}, details={expected_infeasibility}"
        )

    normalized = dict(summary)
    repaired_fields: list[str] = []
    if parallel_processes > 1:
        internal_dir.mkdir(parents=True, exist_ok=True)
        write_json(internal_dir / "eval_summary.vendor_mpi.json", summary)
        vendor_csv = output_dir / "eval_summary.csv"
        if vendor_csv.exists():
            shutil.copy2(vendor_csv, internal_dir / "eval_summary.vendor_mpi.csv")

        normalized["obj_cumulative"] = serial_cumulative_objective
        normalized["obj_all_cases"] = {
            label: objectives[label] for label in sorted(expected_labels)
        }
        normalized["infeas_cumulative"] = expected_infeasibility
        normalized["infeas_all_cases"] = {
            label: infeasibilities[label] for label in sorted(expected_labels)
        }
        repaired_fields = [
            "obj_cumulative",
            "obj_all_cases",
            "infeas_cumulative",
            "infeas_all_cases",
        ]
        write_json(output_dir / "eval_summary.json", normalized)

        if vendor_csv.exists():
            with vendor_csv.open("r", encoding="utf-8", newline="") as stream:
                rows = list(csv.reader(stream))
            for row in rows:
                if row and row[0] in repaired_fields:
                    row[1:] = [str(normalized[row[0]])]
            with vendor_csv.open("w", encoding="utf-8", newline="") as stream:
                csv.writer(stream).writerows(rows)

    certificate = {
        "schema_version": 1,
        "parallel_processes": parallel_processes,
        "expected_contingency_count": contingency_count,
        "expected_detail_count": contingency_count + 1,
        "observed_detail_count": len(detail_paths),
        "complete_label_set": True,
        "objective_from_details": expected_objective,
        "serial_cumulative_objective_from_details": serial_cumulative_objective,
        "reported_objective": reported_objective,
        "infeasible_labels": infeasible_labels,
        "reported_infeasibility": reported_infeasibility,
        "repaired_vendor_mpi_bookkeeping_fields": repaired_fields,
    }
    write_json(internal_dir / "official_evaluation_certificate.json", certificate)
    return normalized, certificate


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


def cpp_command(
    executable: Path,
    distro: str,
    arguments: list[str],
    timeout: float,
) -> list[str]:
    if timeout <= 0:
        raise ValueError("C++ subprocess timeout must be positive")
    inner_timeout = max(timeout, 0.001)
    return [
        "wsl",
        "-d",
        distro,
        "--",
        "env",
        f"LD_LIBRARY_PATH={WSL_LIBRARY_PATH}",
        "OMP_NUM_THREADS=1",
        "OPENBLAS_NUM_THREADS=1",
        "timeout",
        "--signal=TERM",
        "--kill-after=5s",
        f"{inner_timeout:.3f}s",
        to_wsl(executable),
        *arguments,
    ]


def run_cpp(
    executable: Path,
    distro: str,
    arguments: list[str],
    log_path: Path,
    timeout: float,
) -> subprocess.CompletedProcess[str]:
    if timeout <= 0:
        raise ValueError("C++ subprocess timeout must be positive")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    inner_timeout = max(timeout, 0.001)
    command = cpp_command(executable, distro, arguments, inner_timeout)
    started = time.perf_counter()
    try:
        completed = subprocess.run(
            command,
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=inner_timeout + 10.0,
        )
    except subprocess.TimeoutExpired as error:
        output = error.stdout or ""
        if isinstance(output, bytes):
            output = output.decode(errors="replace")
        log_path.write_text(output + "\nOUTER_TIMEOUT\n", encoding="utf-8")
        raise CompetitionTimeout(
            f"C++ subprocess exceeded its {inner_timeout:.3f}s inner deadline: {arguments}"
        ) from error
    output = completed.stdout or ""
    timed_out = completed.returncode == 124
    if timed_out:
        output += "\nTIMEOUT\n"
    log_path.write_text(output, encoding="utf-8")
    completed.wall_seconds = time.perf_counter() - started  # type: ignore[attr-defined]
    completed.timed_out = timed_out  # type: ignore[attr-defined]
    return completed


def contingency_records(case: dict[str, Any]) -> list[dict[str, Any]]:
    records = [dict(item) for item in case["gen_contingencies"]]
    records.extend(dict(item) for item in case["branch_contingencies"])
    return sorted(records, key=lambda item: item["label"])


def longest_first_contingencies(
    case: dict[str, Any],
    base_state: dict[str, Any],
    records: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    generators = ordered(case, "gen")
    branches = ordered(case, "branch")
    generator_position = {
        int(item["index"]): position for position, item in enumerate(generators)
    }
    branch_position = {
        int(item["index"]): position for position, item in enumerate(branches)
    }
    scheduled: list[dict[str, Any]] = []
    for source in records:
        item = dict(source)
        component = int(item["idx"])
        if item["type"] == "gen":
            position = generator_position[component]
            score = math.hypot(
                float(base_state["pg"][position]),
                float(base_state["qg"][position]),
            )
        else:
            position = branch_position[component]
            score = max(
                math.hypot(
                    float(base_state["pf"][position]),
                    float(base_state["qf"][position]),
                ),
                math.hypot(
                    float(base_state["pt"][position]),
                    float(base_state["qt"][position]),
                ),
            )
        item["schedule_score_base_apparent_power"] = score
        scheduled.append(item)
    scheduled.sort(
        key=lambda item: (
            -float(item["schedule_score_base_apparent_power"]),
            str(item["label"]),
        )
    )
    for rank, item in enumerate(scheduled, start=1):
        item["schedule_rank"] = rank
    return scheduled


def load_fallback_schedule_profile(
    path: Path,
    case_sha256: str,
    contingency_labels: set[str],
    worker_count: int,
) -> tuple[dict[str, dict[str, float | int]], dict[str, Any]]:
    reject_onedrive(path)
    raw = read_json(path)
    if not isinstance(raw, dict) or raw.get("schema_version") != 1:
        raise ValueError("fallback schedule profile must use schema_version 1")
    if raw.get("case_sha256") != case_sha256:
        raise ValueError("fallback schedule profile case hash does not match")
    if int(raw.get("worker_count", -1)) != worker_count:
        raise ValueError("fallback schedule profile worker count does not match")
    assignments = raw.get("assignments")
    if not isinstance(assignments, list) or not assignments:
        raise ValueError("fallback schedule profile assignments must be nonempty")

    by_label: dict[str, dict[str, float | int]] = {}
    predicted_loads = [0.0] * worker_count
    for assignment in assignments:
        if not isinstance(assignment, dict):
            raise ValueError("fallback schedule assignment must be an object")
        label = str(assignment.get("label", ""))
        worker = int(assignment.get("worker", -1))
        predicted = float(assignment.get("predicted_wall_seconds", math.nan))
        if label not in contingency_labels:
            raise ValueError(f"unknown profiled contingency label: {label}")
        if label in by_label:
            raise ValueError(f"duplicate profiled contingency label: {label}")
        if worker < 0 or worker >= worker_count:
            raise ValueError(f"invalid profiled worker for {label}: {worker}")
        if not math.isfinite(predicted) or predicted <= 0.0:
            raise ValueError(f"invalid predicted wall time for {label}: {predicted}")
        by_label[label] = {
            "worker": worker,
            "predicted_wall_seconds": predicted,
        }
        predicted_loads[worker] += predicted

    metadata = {
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "schema_version": 1,
        "worker_count": worker_count,
        "assignment_count": len(by_label),
        "predicted_worker_load_seconds": predicted_loads,
        "worker_assignments_used_for_dispatch": False,
        "dispatch_priority": "descending predicted wall seconds",
        "uses_prior_solution_state": False,
    }
    return by_label, metadata


def apply_fallback_schedule_profile(
    contingencies: list[dict[str, Any]],
    profile: dict[str, dict[str, float | int]],
) -> list[dict[str, Any]]:
    scheduled = [dict(item) for item in contingencies]
    prior_rank = {
        str(item["label"]): int(item.get("schedule_rank", position + 1))
        for position, item in enumerate(scheduled)
    }
    scheduled.sort(
        key=lambda item: (
            0 if str(item["label"]) in profile else 1,
            -float(profile[str(item["label"])]["predicted_wall_seconds"])
            if str(item["label"]) in profile else prior_rank[str(item["label"])],
            str(item["label"]),
        )
    )
    for rank, item in enumerate(scheduled, start=1):
        item["schedule_rank"] = rank
        assignment = profile.get(str(item["label"]))
        if assignment is not None:
            item["profiled_fallback_worker"] = int(assignment["worker"])
            item["profiled_fallback_wall_seconds"] = float(
                assignment["predicted_wall_seconds"]
            )
    return scheduled


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
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--fast-workers", type=int, default=8)
    parser.add_argument("--post-screen-workers", type=int)
    parser.add_argument("--distro", default="Ubuntu-24.04")
    parser.add_argument("--base-timeout", type=float, default=900.0)
    parser.add_argument("--contingency-timeout", type=float, default=300.0)
    parser.add_argument("--code1-time-limit", type=float, default=300.0)
    parser.add_argument("--code2-seconds-per-contingency", type=float, default=2.0)
    parser.add_argument("--total-time-limit", type=float, default=300.0)
    parser.add_argument("--evaluation-reserve", type=float, default=7.0)
    parser.add_argument("--python", type=Path, default=DEFAULT_PYTHON)
    parser.add_argument("--evaluator", type=Path, default=DEFAULT_EVALUATOR)
    parser.add_argument("--evaluation-processes", type=int, default=1)
    parser.add_argument("--mpiexec", type=Path)
    parser.add_argument("--skip-evaluation", action="store_true")
    parser.add_argument("--resident-contingency-model", action="store_true")
    parser.add_argument("--ipopt-acceptable-termination", action="store_true")
    parser.add_argument("--fast-power-flow-screen", action="store_true")
    parser.add_argument("--source-status-base", action="store_true")
    parser.add_argument("--validated-source-base", action="store_true")
    parser.add_argument("--robust-contingency-base", action="store_true")
    parser.add_argument("--two-stage-contingency-screen", action="store_true")
    parser.add_argument("--linearized-contingency-fallback", action="store_true")
    parser.add_argument("--linearized-contingency-only", action="store_true")
    parser.add_argument("--longest-first-schedule", action="store_true")
    parser.add_argument("--fallback-schedule-profile", type=Path)
    args = parser.parse_args()

    if args.ipopt_acceptable_termination and not args.resident_contingency_model:
        parser.error(
            "--ipopt-acceptable-termination requires --resident-contingency-model"
        )
    if args.two_stage_contingency_screen and not args.fast_power_flow_screen:
        parser.error(
            "--two-stage-contingency-screen requires --fast-power-flow-screen"
        )
    if args.robust_contingency_base and not args.validated_source_base:
        parser.error(
            "--robust-contingency-base requires --validated-source-base"
        )
    if (args.linearized_contingency_only and
            not args.linearized_contingency_fallback):
        parser.error(
            "--linearized-contingency-only requires "
            "--linearized-contingency-fallback"
        )
    if (args.fallback_schedule_profile is not None and
            not args.two_stage_contingency_screen):
        parser.error(
            "--fallback-schedule-profile requires "
            "--two-stage-contingency-screen"
        )

    for path in (args.case_json, args.case_dir, args.output_dir, args.executable):
        reject_onedrive(path)
    if args.fallback_schedule_profile is not None:
        reject_onedrive(args.fallback_schedule_profile)
    if args.workers < 1:
        raise ValueError("workers must be positive")
    if args.fast_workers < 1:
        raise ValueError("fast workers must be positive")
    if args.post_screen_workers is None:
        args.post_screen_workers = args.workers
    if args.post_screen_workers < args.workers:
        parser.error("--post-screen-workers cannot be less than --workers")
    if (args.post_screen_workers > args.workers and
            not args.two_stage_contingency_screen):
        parser.error(
            "--post-screen-workers expansion requires "
            "--two-stage-contingency-screen"
        )
    if args.evaluation_processes < 1:
        raise ValueError("evaluation processes must be positive")
    if args.evaluation_processes > 1 and args.mpiexec is None:
        parser.error("--evaluation-processes greater than one requires --mpiexec")
    if args.mpiexec is not None:
        reject_onedrive(args.mpiexec)
    if args.code1_time_limit <= 0:
        raise ValueError("Code1 time limit must be positive")
    if args.total_time_limit <= 0:
        raise ValueError("end-to-end time limit must be positive")
    if args.evaluation_reserve <= 0:
        raise ValueError("evaluation reserve must be positive")
    if args.evaluation_reserve + FINALIZATION_RESERVE_SECONDS >= args.total_time_limit:
        raise ValueError("evaluation and finalization reserves exhaust the end-to-end limit")
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise ValueError(f"cold-run output directory is not empty: {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    started_utc = dt.datetime.now(dt.timezone.utc).isoformat()
    wall_start = time.perf_counter()
    total_deadline = wall_start + args.total_time_limit
    work_deadline = total_deadline - args.evaluation_reserve
    evaluation_deadline = total_deadline - FINALIZATION_RESERVE_SECONDS

    case = read_json(args.case_json)
    contingencies = contingency_records(case)
    case_json_sha256 = sha256(args.case_json)
    fallback_schedule: dict[str, dict[str, float | int]] = {}
    fallback_schedule_metadata: dict[str, Any] | None = None
    if args.fallback_schedule_profile is not None:
        fallback_schedule, fallback_schedule_metadata = (
            load_fallback_schedule_profile(
                args.fallback_schedule_profile,
                case_json_sha256,
                {str(item["label"]) for item in contingencies},
                min(args.workers, len(contingencies)),
            )
        )
    code2_limit = code2_time_limit(
        len(contingencies), args.code2_seconds_per_contingency
    )
    internal = args.output_dir / "internal"
    base_json = internal / "base.json"
    records: list[dict[str, Any]] = []
    run_status: dict[str, Any] = {
        "success": False,
        "stage": "initializing",
        "started_at_utc": started_utc,
        "git_revision": git_revision(),
        "competition_timing": {
            "division": 1,
            "code1_time_limit_seconds": args.code1_time_limit,
            "code2_seconds_per_contingency": args.code2_seconds_per_contingency,
            "code2_time_limit_seconds": code2_limit,
            "contingency_count": len(contingencies),
            "timing_semantics": "Code1 and Code2 are separate wall-clock stages",
        },
        "end_to_end_timing": {
            "time_limit_seconds": args.total_time_limit,
            "evaluation_reserve_seconds": args.evaluation_reserve,
            "finalization_reserve_seconds": FINALIZATION_RESERVE_SECONDS,
            "measurement_boundary": "normalized-case loading through official evaluation and result serialization",
        },
        "completed_contingency_count": 0,
        "resident_contingency_model": args.resident_contingency_model,
        "ipopt_acceptable_termination": args.ipopt_acceptable_termination,
        "fast_power_flow_screen": args.fast_power_flow_screen,
        "source_status_base": args.source_status_base,
        "validated_source_base": args.validated_source_base,
        "robust_contingency_base": args.robust_contingency_base,
        "two_stage_contingency_screen": args.two_stage_contingency_screen,
        "initial_corrective_worker_count": args.workers,
        "post_screen_corrective_worker_count": args.post_screen_workers,
        "streaming_fallback_overlap": False,
        "linearized_contingency_fallback": args.linearized_contingency_fallback,
        "linearized_contingency_only": args.linearized_contingency_only,
        "longest_first_schedule": args.longest_first_schedule,
        "fallback_schedule_profile": fallback_schedule_metadata,
        "profiled_fallback_global_priority": bool(fallback_schedule),
        "evaluation_processes": args.evaluation_processes,
    }

    def checkpoint() -> None:
        write_json(args.output_dir / "run_status.json", run_status)

    checkpoint()

    base_start = time.perf_counter()
    try:
        base_deadline = min(
            base_start + args.code1_time_limit,
            work_deadline,
        )
        base_timeout = effective_process_timeout(
            min(args.base_timeout, args.code1_time_limit),
            base_deadline,
        )
        if args.validated_source_base:
            base_arguments = [
                "validated-source-base-json",
                to_wsl(args.case_json),
                to_wsl(base_json),
                "fast-only",
            ]
            if args.robust_contingency_base:
                base_arguments.append("robust-contingency-seed")
        else:
            base_arguments = [
                "run-ibr-json", to_wsl(args.case_json), to_wsl(base_json), "0"
            ]
            if args.source_status_base:
                base_arguments.append("source-only")
        base_process = run_cpp(
            args.executable,
            args.distro,
            base_arguments,
            internal / "base.console.log",
            base_timeout,
        )
    except Exception as error:
        run_status.update(
            {
                "stage": "code1",
                "base_process_wall_seconds": time.perf_counter() - base_start,
                "code1_within_limit": False,
                "end_to_end_within_limit": time.perf_counter() <= total_deadline,
                "error": str(error),
                "finished_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                "total_wall_seconds": time.perf_counter() - wall_start,
            }
        )
        checkpoint()
        raise
    base_wall = time.perf_counter() - base_start
    run_status.update(
        {
            "stage": "code1",
            "base_process_wall_seconds": base_wall,
            "code1_within_limit": base_wall <= args.code1_time_limit,
        }
    )
    if base_process.returncode != 0:
        run_status["error"] = (
            "Code1 reached its competition or end-to-end work deadline"
            if getattr(base_process, "timed_out", False)
            else "cold C++ base/IBR solve failed"
        )
        checkpoint()
        if getattr(base_process, "timed_out", False):
            raise CompetitionTimeout(run_status["error"])
        raise RuntimeError(f"{run_status['error']}; see internal/base.console.log")
    if base_wall > args.code1_time_limit:
        run_status["error"] = "Code1 completed after its competition deadline"
        checkpoint()
        raise CompetitionTimeout(run_status["error"])
    if time.perf_counter() > work_deadline:
        run_status["error"] = "Code1 exhausted the time reserved for Code2 and evaluation"
        checkpoint()
        raise CompetitionTimeout(run_status["error"])
    base = read_json(base_json)
    if not base.get("success", False):
        run_status["error"] = "cold C++ base/IBR result was not successful"
        checkpoint()
        raise RuntimeError("cold C++ base/IBR result was not successful")
    commitment = [int(value) for value in base["commitment"]]
    base_state = base["selected_state"]
    write_solution(args.output_dir / "solution_BASECASE.txt", case, base_state, commitment)
    if args.longest_first_schedule:
        contingencies = longest_first_contingencies(case, base_state, contingencies)
        schedule_mode = "descending base apparent power, deterministic label tie-break"
    else:
        contingencies = [dict(item) for item in contingencies]
        for rank, item in enumerate(contingencies, start=1):
            item["schedule_rank"] = rank
        schedule_mode = "deterministic label order"
    if fallback_schedule:
        contingencies = apply_fallback_schedule_profile(
            contingencies, fallback_schedule
        )
        schedule_mode = (
            "global descending profiled fallback wall time, then "
            + schedule_mode
        )
    run_status.update(
        {
            "stage": "code2",
            "base_success": True,
            "base_algorithm_wall_seconds": base["wall_seconds"],
            "contingency_schedule_mode": schedule_mode,
            "contingency_schedule": [
                {
                    "label": item["label"],
                    "rank": item["schedule_rank"],
                    "score": item.get("schedule_score_base_apparent_power"),
                    "profiled_fallback_worker": item.get(
                        "profiled_fallback_worker"
                    ),
                    "profiled_fallback_wall_seconds": item.get(
                        "profiled_fallback_wall_seconds"
                    ),
                }
                for item in contingencies
            ],
        }
    )
    checkpoint()

    contingency_start = time.perf_counter()
    code2_deadline = contingency_start + code2_limit
    contingency_deadline = min(code2_deadline, work_deadline)
    contingency_deadline_name = (
        "end-to-end work deadline"
        if work_deadline <= code2_deadline
        else "Code2 competition deadline"
    )
    abort_contingencies = threading.Event()
    progress_lock = threading.Lock()
    worker_records: list[dict[str, Any]] = []
    screen_records: list[dict[str, Any]] = []
    fast_worker_count = 0
    fast_screen_wall = 0.0

    def contingency_record(
        item: dict[str, Any],
        result: dict[str, Any],
        worker_id: int,
        execution_phase: str,
    ) -> dict[str, Any]:
        return {
            "label": str(item["label"]),
            "type": item["type"],
            "source_index": int(item["idx"]),
            "schedule_rank": int(item["schedule_rank"]),
            "schedule_score_base_apparent_power": item.get(
                "schedule_score_base_apparent_power"
            ),
            "worker_id": worker_id,
            "profiled_fallback_worker": item.get("profiled_fallback_worker"),
            "profiled_fallback_wall_seconds": item.get(
                "profiled_fallback_wall_seconds"
            ),
            "profiled_fallback_stolen": (
                item.get("profiled_fallback_worker") is not None
                and int(item["profiled_fallback_worker"]) != worker_id
            ),
            "execution_phase": execution_phase,
            "solver_wall_seconds": result["solve"]["wall_seconds"],
            "objective": result["solve"]["objective"],
            "max_residual": result["validation"]["max_residual"],
            "solver_status": result["solve"]["status"],
            "solver_iterations": result["solve"].get("iterations", -1),
            "resident_reoptimization": result["solve"].get(
                "resident_reoptimization", False
            ),
            "acceptable_termination_enabled": result["solve"].get(
                "acceptable_termination_enabled", False
            ),
            "model_preparation_wall_seconds": result.get(
                "model_preparation_wall_seconds", 0.0
            ),
            "solver_status_success": result.get("solver_status_success", False),
            "accepted_feasible_nonconverged": result.get(
                "accepted_feasible_nonconverged", False
            ),
            "solution_method": result.get("solution_method", "ipopt_corrective"),
            "fast_power_flow_screen": result.get(
                "fast_power_flow_screen", False
            ),
            "fast_screen": result.get("fast_screen"),
            "precomputed_fast_screen_reference": result.get(
                "precomputed_fast_screen_reference", False
            ),
        }

    def save_secure_result(
        item: dict[str, Any],
        result: dict[str, Any],
        worker_id: int,
        execution_phase: str,
    ) -> None:
        label = str(item["label"])
        write_solution(
            args.output_dir / f"solution_{label}.txt",
            case,
            result["solve"]["state"],
            commitment,
            item,
        )
        record = contingency_record(item, result, worker_id, execution_phase)
        with progress_lock:
            records.append(record)
            run_status["completed_contingency_count"] = len(records)
            run_status["last_completed_contingency"] = label
            checkpoint()
            print(
                f"completed {len(records)}/{len(contingencies)}: "
                f"{label} on {execution_phase} worker {worker_id}",
                flush=True,
            )

    exact_contingencies = list(contingencies)
    task_queue: queue.Queue[dict[str, Any]] = queue.Queue()
    worker_count = min(args.workers, len(contingencies))
    post_screen_worker_count = min(
        args.post_screen_workers, len(contingencies)
    )
    profiled_task_queue: queue.PriorityQueue[
        tuple[float, int, str, dict[str, Any]]
    ] = queue.PriorityQueue()
    screening_finished = threading.Event()
    fast_pool: concurrent.futures.ThreadPoolExecutor | None = None
    fast_futures: dict[
        concurrent.futures.Future[dict[str, Any]], int
    ] = {}
    if args.two_stage_contingency_screen:
        fast_screen_start = time.perf_counter()
        screen_queue: queue.Queue[dict[str, Any]] = queue.Queue()
        for item in contingencies:
            screen_queue.put(item)
        fallback_items: list[dict[str, Any]] = []
        fast_worker_count = min(args.fast_workers, len(contingencies))

        def screen_worker(worker_id: int) -> dict[str, Any]:
            log_path = (
                internal / "fast_screen_worker_logs" / f"worker_{worker_id:03d}.log"
            )
            log_path.parent.mkdir(parents=True, exist_ok=True)
            task_timeout = effective_process_timeout(
                args.contingency_timeout, contingency_deadline
            )
            command = cpp_command(
                args.executable,
                args.distro,
                [
                    "contingency-worker",
                    to_wsl(args.case_json),
                    to_wsl(base_json),
                    "0",
                    "fast-pf",
                    "fast-only",
                ],
                task_timeout,
            )
            started = time.perf_counter()
            output_lines: list[str] = []
            assigned_labels: list[str] = []
            process = subprocess.Popen(
                command,
                text=True,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                bufsize=1,
            )

            def read_until(prefix: str) -> str:
                assert process.stdout is not None
                while True:
                    line = process.stdout.readline()
                    if line == "":
                        return_code = process.wait()
                        if return_code == 124:
                            raise CompetitionTimeout(
                                f"fast-screen worker {worker_id} reached the "
                                f"{contingency_deadline_name}"
                            )
                        raise RuntimeError(
                            f"fast-screen worker {worker_id} exited with status "
                            f"{return_code}; see {log_path}"
                        )
                    output_lines.append(line)
                    stripped = line.rstrip("\r\n")
                    if stripped.startswith(prefix):
                        return stripped[len(prefix):].strip()

            try:
                read_until("GRAVITYX_WORKER_READY")
                while not abort_contingencies.is_set():
                    try:
                        item = screen_queue.get_nowait()
                    except queue.Empty:
                        break
                    label = str(item["label"])
                    assigned_labels.append(label)
                    result_path = (
                        internal / "contingencies" / f"{safe_label(label)}.json"
                    )
                    assert process.stdin is not None
                    process.stdin.write(
                        json.dumps(
                            {
                                "label": label,
                                "output_path": to_wsl(result_path),
                            },
                            separators=(",", ":"),
                        )
                        + "\n"
                    )
                    process.stdin.flush()
                    acknowledgement = json.loads(
                        read_until("GRAVITYX_TASK_RESULT ")
                    )
                    if acknowledgement.get("label") != label or not acknowledgement.get(
                        "success", False
                    ):
                        raise RuntimeError(
                            f"fast screen {label} failed to execute on worker "
                            f"{worker_id}; see {log_path}"
                        )
                    if time.perf_counter() > contingency_deadline:
                        raise CompetitionTimeout(
                            f"fast screen {label} finished after the "
                            f"{contingency_deadline_name}"
                        )
                    result = read_json(result_path)
                    if result.get("success", False):
                        save_secure_result(item, result, worker_id, "fast_screen")
                    elif result.get("screen_completed", False):
                        with progress_lock:
                            fallback_items.append(item)
                            assignment = fallback_schedule.get(label)
                            if assignment is None:
                                task_queue.put(item)
                            else:
                                profiled_task_queue.put(
                                    (
                                        -float(
                                            assignment["predicted_wall_seconds"]
                                        ),
                                        int(item["schedule_rank"]),
                                        label,
                                        item,
                                    )
                                )
                    else:
                        raise RuntimeError(
                            f"fast screen {label} returned an incomplete result"
                        )
                    with progress_lock:
                        screen_records.append(
                            {
                                "label": label,
                                "worker_id": worker_id,
                                "feasible": bool(result.get("success", False)),
                                "wall_seconds": result["solve"]["wall_seconds"],
                                "max_residual": result["validation"]["max_residual"],
                                "failure_reason": (
                                    result.get("fast_screen") or {}
                                ).get("failure_reason", ""),
                            }
                        )
                        run_status["screened_contingency_count"] = len(screen_records)
                        run_status["fast_screen_fallback_count"] = len(fallback_items)
                        checkpoint()
                    screen_queue.task_done()
                if process.poll() is None:
                    assert process.stdin is not None
                    process.stdin.write('{"stop":true}\n')
                    process.stdin.flush()
                    process.stdin.close()
                    assert process.stdout is not None
                    output_lines.extend(process.stdout.readlines())
                return_code = process.wait(timeout=10.0)
                if return_code != 0:
                    raise RuntimeError(
                        f"fast-screen worker {worker_id} exited with status "
                        f"{return_code}; see {log_path}"
                    )
                return {
                    "phase": "fast_screen",
                    "worker_id": worker_id,
                    "task_count": len(assigned_labels),
                    "process_wall_seconds": time.perf_counter() - started,
                    "labels": assigned_labels,
                }
            except Exception:
                abort_contingencies.set()
                raise
            finally:
                if process.poll() is None:
                    process.terminate()
                    try:
                        process.wait(timeout=2.0)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        process.wait()
                log_path.write_text("".join(output_lines), encoding="utf-8")

        fast_pool = concurrent.futures.ThreadPoolExecutor(
            max_workers=fast_worker_count
        )
        fast_futures = {
            fast_pool.submit(screen_worker, worker_id): worker_id
            for worker_id in range(fast_worker_count)
        }
    else:
        for item in exact_contingencies:
            task_queue.put(item)
        screening_finished.set()

    if not args.two_stage_contingency_screen:
        worker_count = min(args.workers, len(exact_contingencies))

    def solve_worker(worker_id: int) -> dict[str, Any]:
        if abort_contingencies.is_set():
            raise CompetitionTimeout("Code2 was cancelled after another worker failed")
        log_path = internal / "worker_logs" / f"worker_{worker_id:03d}.log"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        task_timeout = effective_process_timeout(
            args.contingency_timeout, contingency_deadline
        )
        worker_arguments = [
            "contingency-worker",
            to_wsl(args.case_json),
            to_wsl(base_json),
            "0",
        ]
        if args.resident_contingency_model:
            worker_arguments.append("resident")
        if args.ipopt_acceptable_termination:
            worker_arguments.append("acceptable")
        # Recompute the fast-screen candidate in second-stage workers as well.
        # The linearized repair uses that candidate as its reference; dropping
        # it between stages can make an otherwise tractable trust-region LP
        # appear infeasible.
        if args.fast_power_flow_screen:
            worker_arguments.append("fast-pf")
        if args.linearized_contingency_fallback:
            worker_arguments.append("linearized")
        if args.linearized_contingency_only:
            worker_arguments.append("linearized-only")
        command = cpp_command(
            args.executable,
            args.distro,
            worker_arguments,
            task_timeout,
        )
        started = time.perf_counter()
        output_lines: list[str] = []
        assigned_labels: list[str] = []
        stolen_labels: list[str] = []
        process = subprocess.Popen(
            command,
            text=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            bufsize=1,
        )

        def read_until(prefix: str) -> str:
            assert process.stdout is not None
            while True:
                line = process.stdout.readline()
                if line == "":
                    return_code = process.wait()
                    if return_code == 124:
                        raise CompetitionTimeout(
                            f"contingency worker {worker_id} reached the "
                            f"{contingency_deadline_name}"
                        )
                    raise RuntimeError(
                        f"contingency worker {worker_id} exited with status "
                        f"{return_code}; see {log_path}"
                    )
                output_lines.append(line)
                stripped = line.rstrip("\r\n")
                if stripped.startswith(prefix):
                    return stripped[len(prefix):].strip()

        try:
            read_until("GRAVITYX_WORKER_READY")
            while not abort_contingencies.is_set():
                item = streamed_queue_get(
                    task_queue,
                    screening_finished,
                    abort_contingencies,
                    profiled_queue=(
                        profiled_task_queue if fallback_schedule else None
                    ),
                )
                if item is None:
                    break
                label = str(item["label"])
                assigned_labels.append(label)
                profiled_worker = item.get("profiled_fallback_worker")
                if (
                    profiled_worker is not None
                    and int(profiled_worker) != worker_id
                ):
                    stolen_labels.append(label)
                result_path = internal / "contingencies" / f"{safe_label(label)}.json"
                task = {
                    "label": label,
                    "output_path": to_wsl(result_path),
                }
                if args.two_stage_contingency_screen:
                    task["fast_screen_path"] = to_wsl(result_path)
                assert process.stdin is not None
                process.stdin.write(
                    json.dumps(
                        task,
                        separators=(",", ":"),
                    )
                    + "\n"
                )
                process.stdin.flush()
                acknowledgement = json.loads(
                    read_until("GRAVITYX_TASK_RESULT ")
                )
                if acknowledgement.get("label") != label:
                    raise RuntimeError(
                        f"contingency worker {worker_id} acknowledged the wrong task"
                    )
                if (args.two_stage_contingency_screen and
                        not acknowledgement.get(
                            "precomputed_fast_screen_reference", False)):
                    raise RuntimeError(
                        f"contingency worker {worker_id} did not accept the "
                        f"precomputed fast-screen reference for {label}"
                    )
                if not acknowledgement.get("success", False):
                    raise RuntimeError(
                        f"contingency {label} failed in worker {worker_id}; see {log_path}"
                    )
                if time.perf_counter() > contingency_deadline:
                    raise CompetitionTimeout(
                        f"contingency {label} finished after the "
                        f"{contingency_deadline_name}"
                    )
                result = read_json(result_path)
                if not result.get("success", False):
                    raise RuntimeError(
                        f"contingency {label} did not pass independent validation"
                    )
                save_secure_result(
                    item,
                    result,
                    worker_id,
                    (
                        "linearized_fallback"
                        if args.linearized_contingency_fallback
                        else "exact_fallback"
                    )
                    if args.two_stage_contingency_screen
                    else "combined_screen_and_fallback",
                )
            if process.poll() is None:
                assert process.stdin is not None
                process.stdin.write('{"stop":true}\n')
                process.stdin.flush()
                process.stdin.close()
                assert process.stdout is not None
                output_lines.extend(process.stdout.readlines())
            return_code = process.wait(timeout=10.0)
            if return_code != 0:
                if return_code == 124:
                    raise CompetitionTimeout(
                        f"contingency worker {worker_id} reached the "
                        f"{contingency_deadline_name}"
                    )
                raise RuntimeError(
                    f"contingency worker {worker_id} exited with status "
                    f"{return_code}; see {log_path}"
                )
            return {
                "phase": (
                    "linearized_fallback"
                    if args.linearized_contingency_fallback
                    else "exact_fallback"
                ),
                "worker_id": worker_id,
                "task_count": len(assigned_labels),
                "process_wall_seconds": time.perf_counter() - started,
                "labels": assigned_labels,
                "profiled_stolen_labels": stolen_labels,
            }
        except Exception:
            abort_contingencies.set()
            raise
        finally:
            if process.poll() is None:
                try:
                    assert process.stdin is not None
                    process.stdin.write('{"stop":true}\n')
                    process.stdin.flush()
                    process.stdin.close()
                    process.wait(timeout=2.0)
                except Exception:
                    process.terminate()
                    try:
                        process.wait(timeout=2.0)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        process.wait()
            log_path.write_text("".join(output_lines), encoding="utf-8")

    pool = concurrent.futures.ThreadPoolExecutor(
        max_workers=max(1, post_screen_worker_count)
    )
    futures: dict[concurrent.futures.Future[dict[str, Any]], int] = {}
    try:
        futures = {
            pool.submit(solve_worker, worker_id): worker_id
            for worker_id in range(worker_count)
        }
        if args.two_stage_contingency_screen:
            assert fast_pool is not None
            for future in concurrent.futures.as_completed(fast_futures):
                worker_records.append(future.result())
            fast_pool.shutdown(wait=True)
            fast_pool = None
            fast_screen_wall = time.perf_counter() - fast_screen_start
            if len(screen_records) != len(contingencies):
                raise RuntimeError(
                    f"fast-screen queue completed {len(screen_records)} of "
                    f"{len(contingencies)} tasks"
                )
            exact_contingencies = sorted(
                fallback_items, key=lambda item: int(item["schedule_rank"])
            )
            run_status["fast_screen_feasible_count"] = sum(
                bool(item["feasible"]) for item in screen_records
            )
            run_status["fast_screen_fallback_count"] = len(exact_contingencies)
            run_status["streaming_fallback_overlap"] = True
            if post_screen_worker_count > worker_count:
                for worker_id in range(worker_count, post_screen_worker_count):
                    future = pool.submit(solve_worker, worker_id)
                    futures[future] = worker_id
                print(
                    "expanded corrective worker pool after screening: "
                    f"{worker_count} -> {post_screen_worker_count}",
                    flush=True,
                )
            run_status["active_post_screen_corrective_worker_count"] = (
                post_screen_worker_count
            )
            checkpoint()
            screening_finished.set()
        for future in concurrent.futures.as_completed(futures):
            worker = future.result()
            worker_records.append(
                {
                    "phase": worker.get("phase", "exact_fallback"),
                    "worker_id": worker["worker_id"],
                    "task_count": worker["task_count"],
                    "process_wall_seconds": worker["process_wall_seconds"],
                    "labels": worker["labels"],
                    "profiled_stolen_labels": worker.get(
                        "profiled_stolen_labels", []
                    ),
                }
            )
            with progress_lock:
                run_status["completed_worker_count"] = len(worker_records)
                run_status["last_completed_worker"] = worker["worker_id"]
                checkpoint()
    except Exception as error:
        abort_contingencies.set()
        screening_finished.set()
        for future in fast_futures:
            future.cancel()
        for future in futures:
            future.cancel()
        if fast_pool is not None:
            fast_pool.shutdown(wait=True, cancel_futures=True)
        pool.shutdown(wait=True, cancel_futures=True)
        contingency_wall = time.perf_counter() - contingency_start
        run_status.update(
            {
                "stage": "code2",
                "code2_within_limit": contingency_wall <= code2_limit,
                "end_to_end_within_limit": time.perf_counter() <= total_deadline,
                "contingency_parallel_wall_seconds": contingency_wall,
                "error": str(error),
                "finished_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                "total_wall_seconds": time.perf_counter() - wall_start,
            }
        )
        checkpoint()
        raise
    else:
        pool.shutdown(wait=True)
    if len(records) != len(contingencies):
        raise RuntimeError(
            f"dynamic contingency queue completed {len(records)} of "
            f"{len(contingencies)} tasks"
        )
    records.sort(key=lambda item: item["label"])
    worker_records.sort(key=lambda item: item["worker_id"])
    contingency_wall = time.perf_counter() - contingency_start
    run_status.update(
        {
            "code2_within_limit": contingency_wall <= code2_limit,
            "contingency_parallel_wall_seconds": contingency_wall,
        }
    )
    if contingency_wall > code2_limit:
        run_status["error"] = "Code2 completed after its competition deadline"
        checkpoint()
        raise CompetitionTimeout(run_status["error"])
    if time.perf_counter() > work_deadline:
        run_status["error"] = "Code2 exhausted the time reserved for official evaluation"
        checkpoint()
        raise CompetitionTimeout(run_status["error"])
    checkpoint()

    evaluation_wall = 0.0
    evaluation_summary: dict[str, Any] | None = None
    evaluation_certificate: dict[str, Any] | None = None
    if not args.skip_evaluation:
        for path in (args.python, args.evaluator):
            reject_onedrive(path)
        if args.mpiexec is not None:
            reject_onedrive(args.mpiexec)
        stale_evaluation_outputs = [
            path
            for pattern in (
                "eval_detail_*.json",
                "eval_summary.json",
                "eval_summary.csv",
                "eval_detail.json",
                "eval_detail.csv",
                "sol_change_*.csv",
                "sol_change.csv",
            )
            for path in args.output_dir.glob(pattern)
        ]
        if stale_evaluation_outputs:
            raise RuntimeError(
                "refusing stale official-evaluation artifacts: "
                + ", ".join(sorted(path.name for path in stale_evaluation_outputs))
            )
        evaluation_start = time.perf_counter()
        evaluator_timeout = effective_process_timeout(
            args.total_time_limit,
            evaluation_deadline,
        )
        evaluator_command = [
            str(args.python),
            str(args.evaluator),
            "1",
            str(args.case_dir),
            str(args.output_dir),
        ]
        if args.evaluation_processes > 1:
            evaluator_command = [
                str(args.mpiexec),
                "-np",
                str(args.evaluation_processes),
                *evaluator_command,
            ]
        try:
            completed = subprocess.run(
                evaluator_command,
                check=False,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=evaluator_timeout,
            )
        except subprocess.TimeoutExpired as error:
            output = error.stdout or ""
            if isinstance(output, bytes):
                output = output.decode(errors="replace")
            (internal / "evaluation.console.log").write_text(
                output + "\nEND_TO_END_TIMEOUT\n", encoding="utf-8"
            )
            run_status.update(
                {
                    "stage": "evaluation",
                    "error": "official evaluator reached the end-to-end deadline",
                    "end_to_end_within_limit": False,
                    "finished_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                    "total_wall_seconds": time.perf_counter() - wall_start,
                }
            )
            checkpoint()
            raise CompetitionTimeout(run_status["error"]) from error
        evaluation_wall = time.perf_counter() - evaluation_start
        (internal / "evaluation.console.log").write_text(completed.stdout, encoding="utf-8")
        if completed.returncode != 0:
            run_status.update(
                {
                    "stage": "evaluation",
                    "error": "official evaluator failed",
                    "end_to_end_within_limit": time.perf_counter() <= total_deadline,
                    "finished_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                    "total_wall_seconds": time.perf_counter() - wall_start,
                }
            )
            checkpoint()
            raise RuntimeError("official evaluator failed; see internal/evaluation.console.log")
        evaluation_summary = read_json(args.output_dir / "eval_summary.json")
        evaluation_summary, evaluation_certificate = (
            validate_and_normalize_evaluation_details(
                args.output_dir,
                internal,
                {str(item["label"]) for item in contingencies},
                evaluation_summary,
                args.evaluation_processes,
            )
        )
        evaluation_certificate.update(
            {
                "evaluator_path": str(args.evaluator.resolve()),
                "evaluator_sha256": sha256(args.evaluator),
                "python_path": str(args.python.resolve()),
                "mpiexec_path": (
                    str(args.mpiexec.resolve()) if args.mpiexec is not None else None
                ),
                "mpiexec_sha256": (
                    sha256(args.mpiexec) if args.mpiexec is not None else None
                ),
            }
        )
        write_json(
            internal / "official_evaluation_certificate.json",
            evaluation_certificate,
        )
        infeasible_cases = [
            label
            for label, infeasible in evaluation_summary.get("infeas_all_cases", {}).items()
            if bool(infeasible)
        ]
        if (
            not bool(evaluation_summary.get("solutions_exist", False))
            or float(evaluation_summary.get("infeas", 1.0)) != 0.0
            or infeasible_cases
        ):
            run_status.update(
                {
                    "stage": "evaluation",
                    "error": "official evaluator did not certify every case feasible",
                    "end_to_end_within_limit": time.perf_counter() <= total_deadline,
                    "official_infeasibility": evaluation_summary.get("infeas"),
                    "official_infeasible_cases": infeasible_cases,
                    "finished_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                    "total_wall_seconds": time.perf_counter() - wall_start,
                }
            )
            checkpoint()
            raise RuntimeError(run_status["error"])

    if time.perf_counter() > evaluation_deadline:
        run_status.update(
            {
                "stage": "evaluation",
                "error": "official evaluation left insufficient time for result serialization",
                "end_to_end_within_limit": time.perf_counter() <= total_deadline,
                "finished_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                "total_wall_seconds": time.perf_counter() - wall_start,
            }
        )
        checkpoint()
        raise CompetitionTimeout(run_status["error"])

    total_wall = time.perf_counter() - wall_start
    max_residual = max((item["max_residual"] for item in records), default=0.0)
    summary = {
        "method": (
            "C++ source commitment with HiGHS sequential-linearized AC base, "
            "parallel sparse-Newton screening, and isolated resident Ipopt fallback"
            if args.validated_source_base and args.two_stage_contingency_screen
            else (
            "Gravity C++ source-status AC base plus validated sparse-Newton contingency screen"
            if args.source_status_base and args.fast_power_flow_screen
            else "Gravity C++ continuous AC-UC relaxation plus deterministic iterative batch rounding"
            )
        ),
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
        "workers": worker_count,
        "requested_workers": args.workers,
        "fast_workers": fast_worker_count,
        "requested_fast_workers": args.fast_workers,
        "resident_contingency_model": args.resident_contingency_model,
        "ipopt_acceptable_termination": args.ipopt_acceptable_termination,
        "fast_power_flow_screen": args.fast_power_flow_screen,
        "source_status_base": args.source_status_base,
        "validated_source_base": args.validated_source_base,
        "robust_contingency_base": args.robust_contingency_base,
        "two_stage_contingency_screen": args.two_stage_contingency_screen,
        "linearized_contingency_fallback": args.linearized_contingency_fallback,
        "linearized_contingency_only": args.linearized_contingency_only,
        "ipopt_acceptable_options": (
            {
                "acceptable_tol": 1e-3,
                "acceptable_iter": 3,
                "acceptable_constr_viol_tol": 5e-6,
                "acceptable_dual_inf_tol": 1e3,
                "acceptable_compl_inf_tol": 1e-3,
                "acceptable_obj_change_tol": 1e-7,
            }
            if args.ipopt_acceptable_termination
            else None
        ),
        "longest_first_schedule": args.longest_first_schedule,
        "contingency_schedule_mode": run_status["contingency_schedule_mode"],
        "contingency_schedule": run_status["contingency_schedule"],
        "contingency_execution_mode": (
            (
                "parallel fast-only sparse-Newton screen followed by "
                "HiGHS sequential-linearized contingency fallback queue"
                if args.linearized_contingency_fallback
                else "parallel fast-only sparse-Newton screen followed by "
                     "resident Ipopt fallback queue"
            )
            if args.two_stage_contingency_screen
            else (
            "validated sparse-Newton screen with resident Ipopt fallback per isolated worker"
            if args.fast_power_flow_screen and args.resident_contingency_model
            else (
                "resident parametric model per isolated process worker with dynamic queue"
                if args.resident_contingency_model
                else "fresh model per task in persistent isolated process workers with dynamic queue"
            )
            )
        ),
        "competition_timing": {
            **run_status["competition_timing"],
            "code1_within_limit": True,
            "code2_within_limit": True,
        },
        "end_to_end_timing": {
            **run_status["end_to_end_timing"],
            "within_limit": True,
        },
        "base_process_wall_seconds": base_wall,
        "base_algorithm_wall_seconds": base["wall_seconds"],
        "contingency_parallel_wall_seconds": contingency_wall,
        "fast_screen_parallel_wall_seconds": fast_screen_wall,
        "contingency_worker_process_seconds_sum": sum(
            item["process_wall_seconds"] for item in worker_records
        ),
        "contingency_solver_seconds_sum": sum(item["solver_wall_seconds"] for item in records),
        "contingency_model_preparation_seconds_sum": sum(
            item["model_preparation_wall_seconds"] for item in records
        ),
        "contingency_ipopt_iterations_sum": sum(
            max(0, int(item["solver_iterations"])) for item in records
        ),
        "evaluation_wall_seconds": evaluation_wall,
        "evaluation_processes": args.evaluation_processes,
        "official_evaluation_certificate": evaluation_certificate,
        "total_wall_seconds": total_wall,
        "contingency_count": len(records),
        "accepted_feasible_nonconverged_count": sum(
            bool(item["accepted_feasible_nonconverged"]) for item in records
        ),
        "fast_power_flow_accepted_count": sum(
            item["solution_method"] in {
                "fast_newton_power_flow",
                "direct_base_state_outage_candidate",
            }
            for item in records
        ),
        "direct_outage_candidate_accepted_count": sum(
            item["solution_method"] == "direct_base_state_outage_candidate"
            for item in records
        ),
        "ipopt_fallback_count": sum(
            item["solution_method"] == "ipopt_corrective_fallback" for item in records
        ),
        "linearized_contingency_accepted_count": sum(
            item["solution_method"] in {
                "highs_sequential_linearized_contingency",
                "highs_linearized_contingency_plus_fast_newton",
            }
            for item in records
        ),
        "fast_screen_fallback_count": sum(
            not bool(item["feasible"]) for item in screen_records
        ),
        "max_independent_contingency_residual": max_residual,
        "base": base,
        "contingency_workers": worker_records,
        "fast_screen_records": sorted(
            screen_records, key=lambda item: item["label"]
        ),
        "contingencies": records,
        "official_evaluation": evaluation_summary,
    }
    write_json(args.output_dir / "run_summary.json", summary)
    run_status.update(
        {
            "success": True,
            "stage": "complete",
            "finished_at_utc": summary["finished_at_utc"],
            "total_wall_seconds": total_wall,
            "official_objective": evaluation_summary.get("obj") if evaluation_summary else None,
            "official_infeasibility": evaluation_summary.get("infeas") if evaluation_summary else None,
            "completed_contingency_count": len(records),
            "accepted_feasible_nonconverged_count": summary[
                "accepted_feasible_nonconverged_count"
            ],
            "end_to_end_within_limit": True,
        }
    )
    checkpoint()
    finalized_total = time.perf_counter() - wall_start
    if finalized_total > args.total_time_limit:
        run_status.update(
            {
                "success": False,
                "stage": "finalization",
                "error": "result serialization exceeded the end-to-end deadline",
                "end_to_end_within_limit": False,
                "finished_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                "total_wall_seconds": finalized_total,
            }
        )
        checkpoint()
        raise CompetitionTimeout(run_status["error"])
    summary["total_wall_seconds"] = finalized_total
    run_status["total_wall_seconds"] = finalized_total
    write_json(args.output_dir / "run_summary.json", summary)
    checkpoint()
    reported_total = time.perf_counter() - wall_start
    if reported_total > args.total_time_limit:
        run_status.update(
            {
                "success": False,
                "stage": "finalization",
                "error": "final timing records exceeded the end-to-end deadline",
                "end_to_end_within_limit": False,
                "finished_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                "total_wall_seconds": reported_total,
            }
        )
        checkpoint()
        raise CompetitionTimeout(run_status["error"])
    print(
        json.dumps(
            {
                "output_dir": str(args.output_dir),
                "total_wall_seconds": reported_total,
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
