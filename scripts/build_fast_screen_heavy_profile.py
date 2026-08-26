#!/usr/bin/env python3
"""Build a timing-only heavy-screen profile from a prior partial/full run."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

from run_experiment import (
    contingency_records,
    read_json,
    reject_onedrive,
    sha256,
    write_json,
)


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


parse_worker_log_measurements = measured_screen_times


def merge_max_measurements(
    measurements: list[dict[str, float]],
) -> dict[str, float]:
    """Merge repeated cold measurements by conservative maximum time."""
    if not measurements:
        raise ValueError("at least one timing measurement is required")
    merged: dict[str, float] = {}
    for source in measurements:
        for label, wall_seconds in source.items():
            if not math.isfinite(wall_seconds) or wall_seconds < 0.0:
                raise ValueError(
                    f"invalid solver wall time for {label}: {wall_seconds}"
                )
            merged[label] = max(wall_seconds, merged.get(label, 0.0))
    return merged


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-json", type=Path, required=True)
    parser.add_argument(
        "--worker-log-dir",
        type=Path,
        action="append",
        required=True,
        help="Worker-log directory; repeat to merge multiple cold runs.",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threshold", type=float, default=5.0)
    parser.add_argument("--require-complete-coverage", action="store_true")
    args = parser.parse_args()

    for path in (args.case_json, *args.worker_log_dir, args.output):
        reject_onedrive(path)
    if args.output.exists():
        raise ValueError(f"profile output already exists: {args.output}")
    if not math.isfinite(args.threshold) or args.threshold <= 0.0:
        raise ValueError("heavy-screen threshold must be positive")

    measurement_sets = [
        measured_screen_times(path) for path in args.worker_log_dir
    ]
    measured = merge_max_measurements(measurement_sets)
    case = read_json(args.case_json)
    contingency_labels = {
        str(item["label"]) for item in contingency_records(case)
    }
    unknown = set(measured) - contingency_labels
    if unknown:
        raise ValueError(
            "worker logs contain unknown contingencies: "
            + ", ".join(sorted(unknown)[:10])
        )
    missing = contingency_labels - set(measured)
    if args.require_complete_coverage and missing:
        raise ValueError(
            "worker logs do not cover every contingency: "
            + ", ".join(sorted(missing)[:10])
        )
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

    sources: list[dict[str, Any]] = []
    for worker_log_dir, source_measurement in zip(
        args.worker_log_dir, measurement_sets
    ):
        source_run_status_path = (
            worker_log_dir.parent.parent / "run_status.json"
        )
        source_run_status: dict[str, Any] = (
            read_json(source_run_status_path)
            if source_run_status_path.is_file()
            else {}
        )
        sources.append(
            {
                "git_revision": source_run_status.get("git_revision"),
                "scenario_run_name": worker_log_dir.parent.parent.name,
                "completed_task_result_count": len(source_measurement),
                "worker_log_count": len(
                    list(worker_log_dir.glob("worker_*.log"))
                ),
                "source_run_success": source_run_status.get("success"),
            }
        )
    profile = {
        "schema_version": 1,
        "purpose": (
            "Scheduling-only heavy-screen concurrency profile. Contains no "
            "primal, dual, commitment, network, or solver state. Repeated "
            "cold timings are merged by maximum measured solver time."
        ),
        "case_sha256": sha256(args.case_json),
        "heavy_threshold_seconds": args.threshold,
        "source_measurement": {
            "complete_contingency_coverage": not missing,
            "merge_rule": "maximum measured solver wall time by label",
            "runs": sources,
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
                "measurement_run_count": len(measurement_sets),
                "threshold_seconds": args.threshold,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
