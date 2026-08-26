#!/usr/bin/env python3
"""Build a timing-only heavy-screen profile from a prior partial/full run."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

from run_experiment import read_json, reject_onedrive, sha256, write_json


def measured_screen_times(worker_log_dir: Path) -> dict[str, float]:
    measured: dict[str, float] = {}
    logs = sorted(worker_log_dir.glob("worker_*.log"))
    if not logs:
        raise ValueError(f"no fast-screen worker logs found in {worker_log_dir}")
    for log_path in logs:
        with log_path.open(encoding="utf-8", errors="replace") as stream:
            for line in stream:
                if not line.startswith("GRAVITYX_TASK_RESULT "):
                    continue
                acknowledgement = json.loads(line.split(" ", 1)[1])
                label = str(acknowledgement.get("label", ""))
                summary = acknowledgement.get("result_summary")
                if not label or not isinstance(summary, dict):
                    raise ValueError(
                        f"invalid task result in {log_path}: missing compact summary"
                    )
                wall_seconds = float(summary["solve"]["wall_seconds"])
                if not math.isfinite(wall_seconds) or wall_seconds < 0.0:
                    raise ValueError(
                        f"invalid solver wall time for {label}: {wall_seconds}"
                    )
                if label in measured:
                    raise ValueError(f"duplicate task result for {label}")
                measured[label] = wall_seconds
    return measured


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-json", type=Path, required=True)
    parser.add_argument("--worker-log-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threshold", type=float, default=5.0)
    args = parser.parse_args()

    for path in (args.case_json, args.worker_log_dir, args.output):
        reject_onedrive(path)
    if args.output.exists():
        raise ValueError(f"profile output already exists: {args.output}")
    if not math.isfinite(args.threshold) or args.threshold <= 0.0:
        raise ValueError("heavy-screen threshold must be positive")

    measured = measured_screen_times(args.worker_log_dir)
    selected = [
        {
            "label": label,
            "measured_solver_wall_seconds": measured[label],
        }
        for label in sorted(measured)
        if measured[label] >= args.threshold
    ]
    if not selected:
        raise ValueError("no measured screens meet the heavy threshold")

    source_run_status_path = (
        args.worker_log_dir.parent.parent / "run_status.json"
    )
    source_run_status: dict[str, Any] = (
        read_json(source_run_status_path)
        if source_run_status_path.is_file()
        else {}
    )
    profile = {
        "schema_version": 1,
        "purpose": (
            "Scheduling-only heavy-screen concurrency profile. Contains no "
            "primal, dual, commitment, network, or solver state."
        ),
        "case_sha256": sha256(args.case_json),
        "heavy_threshold_seconds": args.threshold,
        "source_measurement": {
            "git_revision": source_run_status.get("git_revision"),
            "scenario_run_name": args.worker_log_dir.parent.parent.name,
            "completed_task_result_count": len(measured),
            "worker_log_count": len(
                list(args.worker_log_dir.glob("worker_*.log"))
            ),
            "source_run_success": source_run_status.get("success"),
        },
        "contingencies": selected,
    }
    write_json(args.output, profile)
    print(
        json.dumps(
            {
                "output": str(args.output),
                "profiled_contingency_count": len(selected),
                "measured_task_count": len(measured),
                "threshold_seconds": args.threshold,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
