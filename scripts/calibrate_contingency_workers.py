#!/usr/bin/env python3
"""Calibrate resident contingency-worker concurrency on a fixed hard subset."""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
from pathlib import Path
import subprocess
import threading
import time
from typing import Any

from run_experiment import (
    AffinityScreenWorkQueue,
    DEFAULT_EXE,
    CompetitionTimeout,
    contingency_records,
    cpp_command,
    fast_screen_affinity_groups,
    load_fast_screen_heavy_profile,
    longest_first_contingencies,
    read_json,
    reject_onedrive,
    safe_label,
    sha256,
    to_wsl,
    write_json,
)


def run_trial(
    case_json: Path,
    base_json: Path,
    executable: Path,
    distro: str,
    output_dir: Path,
    records: list[dict[str, Any]],
    workers: int,
    timeout_seconds: float,
    fast_power_flow_screen: bool,
    fast_only: bool = False,
    task_groups: list[list[dict[str, Any]]] | None = None,
    linearized_fallback: bool = False,
    precomputed_fast_screen_dir: Path | None = None,
    wsl_fast_screen_scratch: bool = False,
    heavy_labels: set[str] | None = None,
    heavy_worker_count: int = 0,
    heavy_label_seconds: dict[str, float] | None = None,
    economic_contingency_polish: bool = False,
    worker_environment: dict[str, str] | None = None,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=False)
    if task_groups is None:
        task_groups = [[item] for item in records]
    if [item["label"] for group in task_groups for item in group] != [
        item["label"] for item in records
    ]:
        raise ValueError("task groups must contain every record exactly once in order")
    task_queue = AffinityScreenWorkQueue(
        task_groups,
        workers,
        heavy_labels=heavy_labels,
        heavy_worker_count=heavy_worker_count,
        heavy_label_seconds=heavy_label_seconds,
    )
    abort = threading.Event()
    deadline = time.perf_counter() + timeout_seconds
    completed: list[dict[str, Any]] = []
    completed_lock = threading.Lock()

    def worker(worker_id: int) -> dict[str, Any]:
        log_path = output_dir / "worker_logs" / f"worker_{worker_id:03d}.log"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        remaining = deadline - time.perf_counter()
        if remaining <= 0:
            raise CompetitionTimeout("calibration deadline expired before worker launch")
        worker_arguments = [
            "contingency-worker",
            to_wsl(case_json),
            to_wsl(base_json),
            "0",
            "resident",
            "acceptable",
        ]
        if fast_power_flow_screen:
            worker_arguments.append("fast-pf")
        if fast_only:
            worker_arguments.append("fast-only")
        if economic_contingency_polish:
            worker_arguments.append("economic-polish")
        if linearized_fallback:
            worker_arguments.append("linearized")
        command = cpp_command(
            executable,
            distro,
            worker_arguments,
            remaining,
            worker_environment,
        )
        output_lines: list[str] = []
        labels: list[str] = []
        affinity_split_group_count = 0
        affinity_split_contingency_count = 0
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
                            f"worker {worker_id} reached the calibration deadline"
                        )
                    raise RuntimeError(
                        f"worker {worker_id} exited with status {return_code}"
                    )
                output_lines.append(line)
                stripped = line.rstrip("\r\n")
                if stripped.startswith(prefix):
                    return stripped[len(prefix) :].strip()

        started = time.perf_counter()
        try:
            read_until("GRAVITYX_WORKER_READY")
            while not abort.is_set():
                work = task_queue.get(abort, worker_id)
                if work is None:
                    break
                work_source, group = work
                for group_position, item in enumerate(group):
                    label = str(item["label"])
                    labels.append(label)
                    result_path = (
                        output_dir / "results" / f"{safe_label(label)}.json"
                    )
                    assert process.stdin is not None
                    task = {
                        "label": label,
                        "output_path": (
                            f"/dev/shm/gravityx_calibrate_{os.getpid()}_"
                            f"{worker_id:03d}.json"
                            if wsl_fast_screen_scratch
                            else to_wsl(result_path)
                        ),
                    }
                    if wsl_fast_screen_scratch:
                        task.update(
                            {
                                "fallback_output_path": to_wsl(result_path),
                                "remove_output_after_result": True,
                                "return_result_summary": True,
                            }
                        )
                    if precomputed_fast_screen_dir is not None:
                        task["fast_screen_path"] = to_wsl(
                            precomputed_fast_screen_dir
                            / f"{safe_label(label)}.json"
                        )
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
                            f"worker {worker_id} acknowledged wrong task"
                        )
                    if not acknowledgement.get("success", False):
                        raise RuntimeError(
                            f"calibration contingency {label} failed"
                        )
                    if wsl_fast_screen_scratch:
                        result = acknowledgement.get("result_summary")
                        if not isinstance(result, dict):
                            raise RuntimeError(
                                f"calibration contingency {label} returned no "
                                "compact result summary"
                            )
                        if not acknowledgement.get(
                            "transient_output_removed", False
                        ):
                            raise RuntimeError(
                                f"calibration contingency {label} left its "
                                "transient output behind"
                            )
                        requires_fallback = bool(
                            result.get("requires_exact_fallback", False)
                        )
                        persisted = bool(
                            acknowledgement.get(
                                "fallback_result_persisted", False
                            )
                        )
                        if requires_fallback != persisted:
                            raise RuntimeError(
                                f"calibration contingency {label} persistence "
                                "does not match its fallback status"
                            )
                        if result_path.exists() != requires_fallback:
                            raise RuntimeError(
                                f"calibration contingency {label} persistent "
                                "result presence does not match fallback status"
                            )
                    else:
                        result = read_json(result_path)
                    with completed_lock:
                        completed.append(
                            {
                                "label": label,
                                "worker_id": worker_id,
                                "solver_wall_seconds": result["solve"][
                                    "wall_seconds"
                                ],
                                "model_preparation_wall_seconds": result[
                                    "model_preparation_wall_seconds"
                                ],
                                "iterations": result["solve"].get(
                                    "iterations", -1
                                ),
                                "max_residual": result["validation"][
                                    "max_residual"
                                ],
                                "solution_method": result.get(
                                    "solution_method"
                                ),
                                "secure": bool(result.get("success", False)),
                                "screen_completed": bool(
                                    result.get("screen_completed", False)
                                ),
                                "requires_exact_fallback": bool(
                                    result.get("requires_exact_fallback", False)
                                ),
                                "rolling_corrective_seed_label": result.get(
                                    "rolling_corrective_seed_label"
                                ),
                                "corrective_seed_bank_size": int(
                                    acknowledgement.get(
                                        "corrective_seed_bank_size", 0
                                    )
                                ),
                            }
                        )
                    split_count = 0
                    if result.get("requires_exact_fallback", False):
                        split_count = task_queue.requeue_remaining_as_singletons(
                            group,
                            group_position + 1,
                            worker_id,
                            work_source,
                        )
                        if split_count:
                            affinity_split_group_count += 1
                            affinity_split_contingency_count += split_count
                            break
                task_queue.task_done(work_source)
            if process.poll() is None:
                assert process.stdin is not None
                process.stdin.write('{"stop":true}\n')
                process.stdin.flush()
                process.stdin.close()
                assert process.stdout is not None
                output_lines.extend(process.stdout.readlines())
            return_code = process.wait(timeout=5.0)
            if return_code != 0:
                raise RuntimeError(
                    f"worker {worker_id} exited with status {return_code}"
                )
            return {
                "worker_id": worker_id,
                "labels": labels,
                "wall_seconds": time.perf_counter() - started,
                "screen_lane": task_queue.worker_lane(worker_id),
                "affinity_split_group_count": affinity_split_group_count,
                "affinity_split_contingency_count": (
                    affinity_split_contingency_count
                ),
            }
        except Exception:
            abort.set()
            raise
        finally:
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=3.0)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
            log_path.write_text("".join(output_lines), encoding="utf-8")

    trial_start = time.perf_counter()
    worker_records: list[dict[str, Any]] = []
    error: str | None = None
    pool = concurrent.futures.ThreadPoolExecutor(max_workers=workers)
    futures = [pool.submit(worker, worker_id) for worker_id in range(workers)]
    try:
        for future in concurrent.futures.as_completed(futures):
            worker_records.append(future.result())
    except Exception as exception:
        abort.set()
        error = str(exception)
        for future in futures:
            future.cancel()
    finally:
        pool.shutdown(wait=True, cancel_futures=True)
    wall_seconds = time.perf_counter() - trial_start
    result = {
        "workers": workers,
        "task_count": len(records),
        "completed_count": len(completed),
        "success": error is None and len(completed) == len(records),
        "error": error,
        "wall_seconds": wall_seconds,
        "throughput_per_second": len(completed) / max(wall_seconds, 1e-12),
        "max_residual": max(
            (item["max_residual"] for item in completed), default=None
        ),
        "secure_count": sum(item["secure"] for item in completed),
        "fallback_count": sum(
            item["requires_exact_fallback"] for item in completed
        ),
        "rolling_seed_selected_count": sum(
            item["rolling_corrective_seed_label"] is not None
            for item in completed
        ),
        "affinity_split_group_count": sum(
            int(item.get("affinity_split_group_count", 0))
            for item in worker_records
        ),
        "affinity_split_contingency_count": sum(
            int(item.get("affinity_split_contingency_count", 0))
            for item in worker_records
        ),
        "profiled_heavy_group_count": task_queue.initial_heavy_group_count,
        "heavy_worker_count": heavy_worker_count,
        "solver_seconds_sum": sum(
            item["solver_wall_seconds"] for item in completed
        ),
        "model_preparation_seconds_sum": sum(
            item["model_preparation_wall_seconds"] for item in completed
        ),
        "workers_detail": sorted(worker_records, key=lambda item: item["worker_id"]),
        "contingencies": sorted(completed, key=lambda item: item["label"]),
    }
    write_json(output_dir / "trial_summary.json", result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-json", type=Path, required=True)
    parser.add_argument("--base-json", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--executable", type=Path, default=DEFAULT_EXE)
    parser.add_argument("--distro", default="Ubuntu-24.04")
    parser.add_argument("--workers", type=int, nargs="+", default=[8, 12, 16])
    parser.add_argument("--task-count", type=int, default=32)
    parser.add_argument("--trial-timeout", type=float, default=180.0)
    parser.add_argument("--cooldown", type=float, default=30.0)
    parser.add_argument("--fast-power-flow-screen", action="store_true")
    parser.add_argument("--fast-only", action="store_true")
    parser.add_argument(
        "--economic-contingency-polish", action="store_true"
    )
    parser.add_argument("--fast-screen-affinity-schedule", action="store_true")
    parser.add_argument("--fast-screen-easy-first", action="store_true")
    parser.add_argument("--fast-screen-heavy-profile", type=Path)
    parser.add_argument("--fast-screen-heavy-workers", type=int, default=0)
    parser.add_argument("--additional-easy-task-count", type=int, default=0)
    parser.add_argument("--selection-offset", type=int, default=0)
    parser.add_argument("--labels", nargs="+")
    parser.add_argument("--linearized-fallback", action="store_true")
    parser.add_argument("--precomputed-fast-screen-dir", type=Path)
    parser.add_argument("--wsl-fast-screen-scratch", action="store_true")
    parser.add_argument(
        "--worker-env",
        action="append",
        default=[],
        metavar="GRAVITYX_NAME=VALUE",
        help=(
            "explicitly pass one GRAVITYX_ setting through wsl.exe to every "
            "resident worker; may be repeated"
        ),
    )
    args = parser.parse_args()

    worker_environment: dict[str, str] = {}
    for assignment in args.worker_env:
        if "=" not in assignment:
            parser.error(f"--worker-env requires NAME=VALUE: {assignment!r}")
        name, value = assignment.split("=", 1)
        if not name.startswith("GRAVITYX_"):
            parser.error(
                "--worker-env names must use the GRAVITYX_ prefix: "
                f"{name!r}"
            )
        if name in worker_environment:
            parser.error(f"duplicate --worker-env name: {name}")
        worker_environment[name] = value

    for path in (args.case_json, args.base_json, args.output_dir, args.executable):
        reject_onedrive(path)
    if args.output_dir.exists():
        raise ValueError(f"calibration output already exists: {args.output_dir}")
    if any(workers < 1 for workers in args.workers):
        raise ValueError("worker counts must be positive")
    if args.fast_only and not args.fast_power_flow_screen:
        parser.error("--fast-only requires --fast-power-flow-screen")
    if (args.economic_contingency_polish and
            not args.fast_power_flow_screen):
        parser.error(
            "--economic-contingency-polish requires --fast-power-flow-screen"
        )
    if args.fast_only and args.linearized_fallback:
        parser.error("--fast-only cannot be combined with --linearized-fallback")
    if args.wsl_fast_screen_scratch and not args.fast_only:
        parser.error("--wsl-fast-screen-scratch requires --fast-only")
    if args.precomputed_fast_screen_dir is not None:
        reject_onedrive(args.precomputed_fast_screen_dir)
        if not args.linearized_fallback:
            parser.error(
                "--precomputed-fast-screen-dir requires --linearized-fallback"
            )
    if args.fast_screen_affinity_schedule and not args.fast_only:
        parser.error("--fast-screen-affinity-schedule requires --fast-only")
    if args.fast_screen_easy_first and not args.fast_screen_affinity_schedule:
        parser.error(
            "--fast-screen-easy-first requires "
            "--fast-screen-affinity-schedule"
        )
    if ((args.fast_screen_heavy_profile is None) !=
            (args.fast_screen_heavy_workers == 0)):
        parser.error(
            "--fast-screen-heavy-profile and a positive "
            "--fast-screen-heavy-workers must be provided together"
        )
    if (args.fast_screen_heavy_profile is not None and
            not args.fast_screen_affinity_schedule):
        parser.error(
            "--fast-screen-heavy-profile requires "
            "--fast-screen-affinity-schedule"
        )
    if args.fast_screen_heavy_workers < 0:
        parser.error("--fast-screen-heavy-workers must be nonnegative")
    if any(args.fast_screen_heavy_workers > workers for workers in args.workers):
        parser.error(
            "--fast-screen-heavy-workers cannot exceed any calibrated "
            "worker count"
        )
    if args.selection_offset < 0:
        parser.error("--selection-offset must be nonnegative")
    if args.additional_easy_task_count < 0:
        parser.error("--additional-easy-task-count must be nonnegative")

    case = read_json(args.case_json)
    heavy_labels: set[str] = set()
    heavy_label_seconds: dict[str, float] = {}
    heavy_profile_metadata: dict[str, Any] | None = None
    if args.fast_screen_heavy_profile is not None:
        reject_onedrive(args.fast_screen_heavy_profile)
        (
            heavy_labels,
            heavy_label_seconds,
            heavy_profile_metadata,
        ) = load_fast_screen_heavy_profile(
            args.fast_screen_heavy_profile,
            sha256(args.case_json),
            {
                str(item["label"])
                for item in contingency_records(case)
            },
        )
    base = read_json(args.base_json)
    ordered_records = longest_first_contingencies(
        case,
        base["selected_state"],
        contingency_records(case),
    )
    if args.labels:
        records_by_label = {
            str(item["label"]): item for item in ordered_records
        }
        missing = [label for label in args.labels if label not in records_by_label]
        if missing:
            parser.error(f"unknown contingency labels: {missing}")
        records = [records_by_label[label] for label in args.labels]
    else:
        records = ordered_records[
            args.selection_offset : args.selection_offset + args.task_count
        ]
        if len(records) != args.task_count:
            parser.error(
                "selection offset and task count exceed contingency count"
            )
    if args.additional_easy_task_count:
        primary_labels = {str(item["label"]) for item in records}
        easy_records = [
            item
            for item in reversed(ordered_records)
            if str(item["label"]) not in primary_labels
        ][: args.additional_easy_task_count]
        if len(easy_records) != args.additional_easy_task_count:
            parser.error("not enough distinct easy contingencies to add")
        records.extend(easy_records)
    task_groups: list[list[dict[str, Any]]] | None = None
    if args.fast_screen_affinity_schedule:
        task_groups = fast_screen_affinity_groups(
            case,
            base["selected_state"],
            records,
            difficult_groups_first=not args.fast_screen_easy_first,
        )
        records = [item for group in task_groups for item in group]
    args.output_dir.mkdir(parents=True)
    summary: dict[str, Any] = {
        "purpose": "non-official fixed-subset resident-worker calibration",
        "task_count": len(records),
        "selection_offset": args.selection_offset,
        "fast_only": args.fast_only,
        "fast_screen_affinity_schedule": args.fast_screen_affinity_schedule,
        "fast_screen_easy_first": args.fast_screen_easy_first,
        "fast_screen_heavy_profile": heavy_profile_metadata,
        "fast_screen_heavy_workers": args.fast_screen_heavy_workers,
        "additional_easy_task_count": args.additional_easy_task_count,
        "linearized_fallback": args.linearized_fallback,
        "economic_contingency_polish": args.economic_contingency_polish,
        "wsl_fast_screen_scratch": args.wsl_fast_screen_scratch,
        "precomputed_fast_screen_dir": (
            str(args.precomputed_fast_screen_dir)
            if args.precomputed_fast_screen_dir is not None
            else None
        ),
        "worker_environment": worker_environment,
        "schedule": [item["label"] for item in records],
        "trials": [],
    }
    for position, workers in enumerate(args.workers):
        # The full experiment naturally has a warm WSL VM after Code 1. Keep
        # distro startup outside the Code-2-only calibration boundary too.
        subprocess.run(
            ["wsl", "-d", args.distro, "--", "true"],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        trial = run_trial(
            args.case_json,
            args.base_json,
            args.executable,
            args.distro,
            args.output_dir / f"workers_{workers:02d}",
            records,
            workers,
            args.trial_timeout,
            args.fast_power_flow_screen,
            args.fast_only,
            task_groups,
            args.linearized_fallback,
            args.precomputed_fast_screen_dir,
            args.wsl_fast_screen_scratch,
            heavy_labels,
            args.fast_screen_heavy_workers,
            heavy_label_seconds,
            args.economic_contingency_polish,
            worker_environment,
        )
        summary["trials"].append(trial)
        write_json(args.output_dir / "calibration_summary.json", summary)
        print(json.dumps(trial, indent=2), flush=True)
        if position + 1 < len(args.workers) and args.cooldown > 0:
            time.sleep(args.cooldown)
    successful = [item for item in summary["trials"] if item["success"]]
    if successful:
        best = min(successful, key=lambda item: item["wall_seconds"])
        summary["selected_workers"] = best["workers"]
        summary["selection_basis"] = "minimum fixed-subset wall time"
    else:
        summary["selected_workers"] = None
        summary["selection_basis"] = "no trial completed"
    write_json(args.output_dir / "calibration_summary.json", summary)
    print(json.dumps({"selected_workers": summary["selected_workers"]}, indent=2))
    return 0 if successful else 1


if __name__ == "__main__":
    raise SystemExit(main())
