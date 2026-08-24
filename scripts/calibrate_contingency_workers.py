#!/usr/bin/env python3
"""Calibrate resident contingency-worker concurrency on a fixed hard subset."""

from __future__ import annotations

import argparse
import concurrent.futures
import json
from pathlib import Path
import queue
import subprocess
import threading
import time
from typing import Any

from run_experiment import (
    DEFAULT_EXE,
    CompetitionTimeout,
    contingency_records,
    cpp_command,
    longest_first_contingencies,
    read_json,
    reject_onedrive,
    safe_label,
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
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=False)
    task_queue: queue.Queue[dict[str, Any]] = queue.Queue()
    for item in records:
        task_queue.put(item)
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
        command = cpp_command(
            executable,
            distro,
            worker_arguments,
            remaining,
        )
        output_lines: list[str] = []
        labels: list[str] = []
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
                try:
                    item = task_queue.get_nowait()
                except queue.Empty:
                    break
                label = str(item["label"])
                labels.append(label)
                result_path = output_dir / "results" / f"{safe_label(label)}.json"
                assert process.stdin is not None
                process.stdin.write(
                    json.dumps(
                        {"label": label, "output_path": to_wsl(result_path)},
                        separators=(",", ":"),
                    )
                    + "\n"
                )
                process.stdin.flush()
                acknowledgement = json.loads(
                    read_until("GRAVITYX_TASK_RESULT ")
                )
                if acknowledgement.get("label") != label:
                    raise RuntimeError(f"worker {worker_id} acknowledged wrong task")
                if not acknowledgement.get("success", False):
                    raise RuntimeError(f"calibration contingency {label} failed")
                result = read_json(result_path)
                task_queue.task_done()
                with completed_lock:
                    completed.append(
                        {
                            "label": label,
                            "worker_id": worker_id,
                            "solver_wall_seconds": result["solve"]["wall_seconds"],
                            "model_preparation_wall_seconds": result[
                                "model_preparation_wall_seconds"
                            ],
                            "iterations": result["solve"].get("iterations", -1),
                            "max_residual": result["validation"]["max_residual"],
                            "solution_method": result.get("solution_method"),
                        }
                    )
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
    args = parser.parse_args()

    for path in (args.case_json, args.base_json, args.output_dir, args.executable):
        reject_onedrive(path)
    if args.output_dir.exists():
        raise ValueError(f"calibration output already exists: {args.output_dir}")
    if any(workers < 1 for workers in args.workers):
        raise ValueError("worker counts must be positive")

    case = read_json(args.case_json)
    base = read_json(args.base_json)
    records = longest_first_contingencies(
        case,
        base["selected_state"],
        contingency_records(case),
    )[: args.task_count]
    args.output_dir.mkdir(parents=True)
    summary: dict[str, Any] = {
        "purpose": "non-official fixed-subset resident-worker calibration",
        "task_count": len(records),
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
