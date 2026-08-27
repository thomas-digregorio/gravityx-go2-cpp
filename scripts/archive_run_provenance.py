#!/usr/bin/env python3
"""Archive compact, hash-backed provenance before pruning generated runs."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


INLINE_JSON_LIMIT = 64 * 1024
TAIL_BYTE_LIMIT = 16 * 1024


def reject_onedrive(path: Path) -> None:
    if "onedrive" in str(path.resolve()).casefold():
        raise ValueError(f"OneDrive path is forbidden: {path}")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def file_record(path: Path, *, include_json: bool = False) -> dict[str, Any]:
    data = path.read_bytes()
    record: dict[str, Any] = {
        "name": path.name,
        "bytes": len(data),
        "sha256": sha256_bytes(data),
    }
    if include_json:
        value = json.loads(data.decode("utf-8"))
        record["content"] = compact_json_object(value)
    return record


def compact_json_object(value: Any) -> Any:
    if not isinstance(value, dict):
        encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
        if len(encoded) <= INLINE_JSON_LIMIT:
            return value
        return {
            "omitted": True,
            "type": type(value).__name__,
            "count": len(value) if hasattr(value, "__len__") else None,
            "json_bytes": len(encoded),
            "sha256": sha256_bytes(encoded),
        }

    compact: dict[str, Any] = {}
    for key, item in value.items():
        encoded = json.dumps(item, sort_keys=True, separators=(",", ":")).encode()
        if len(encoded) <= INLINE_JSON_LIMIT:
            compact[key] = item
        else:
            compact[key] = {
                "omitted": True,
                "type": type(item).__name__,
                "count": len(item) if hasattr(item, "__len__") else None,
                "json_bytes": len(encoded),
                "sha256": sha256_bytes(encoded),
            }
    return compact


def inventory_tree(root: Path) -> dict[str, int]:
    file_count = 0
    logical_bytes = 0
    solution_count = 0
    for dirpath, _, filenames in os.walk(root):
        for name in filenames:
            file_count += 1
            if name.startswith("solution_") and name.endswith(".txt"):
                solution_count += 1
            try:
                logical_bytes += (Path(dirpath) / name).stat().st_size
            except OSError:
                pass
    return {
        "file_count": file_count,
        "logical_bytes": logical_bytes,
        "solution_file_entries": solution_count,
    }


def archive_run(status_path: Path, runs_root: Path) -> dict[str, Any]:
    run_root = status_path.parent
    status_bytes = status_path.read_bytes()
    status = json.loads(status_bytes.decode("utf-8"))
    record: dict[str, Any] = {
        "relative_run_root": run_root.relative_to(runs_root).as_posix(),
        "run_status_bytes": len(status_bytes),
        "run_status_sha256": sha256_bytes(status_bytes),
        "run_status": compact_json_object(status),
        "inventory": inventory_tree(run_root),
        "adjacent_json": [],
        "root_log_tails": [],
    }

    for path in sorted(run_root.glob("*.json")):
        if path == status_path:
            continue
        try:
            record["adjacent_json"].append(file_record(path, include_json=True))
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            record["adjacent_json"].append(
                {"name": path.name, "archive_error": str(exc)}
            )

    for pattern in ("*.log", "*.out", "*.err"):
        for path in sorted(run_root.glob(pattern)):
            try:
                data = path.read_bytes()
            except OSError as exc:
                record["root_log_tails"].append(
                    {"name": path.name, "archive_error": str(exc)}
                )
                continue
            record["root_log_tails"].append(
                {
                    "name": path.name,
                    "bytes": len(data),
                    "sha256": sha256_bytes(data),
                    "tail": data[-TAIL_BYTE_LIMIT:].decode("utf-8", errors="replace"),
                }
            )
    return record


def archive_diagnostics(diagnostics_root: Path, pattern: str) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    if not diagnostics_root.exists():
        return records
    for path in sorted(diagnostics_root.glob(pattern)):
        if path.is_dir():
            records.append(
                {
                    "name": path.name,
                    "inventory": inventory_tree(path),
                }
            )
    return records


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--run-status-glob", default="c2fen19402_*/**/run_status.json")
    parser.add_argument("--diagnostics-root", type=Path)
    parser.add_argument("--diagnostics-glob", default="c2fen19402_*")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    runs_root = args.runs_root.resolve()
    output = args.output.resolve()
    reject_onedrive(runs_root)
    reject_onedrive(output)
    if args.diagnostics_root is not None:
        reject_onedrive(args.diagnostics_root)

    status_paths = sorted(runs_root.glob(args.run_status_glob))
    archive = {
        "schema_version": 1,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "runs_root": str(runs_root),
        "run_status_glob": args.run_status_glob,
        "run_count": len(status_paths),
        "runs": [archive_run(path, runs_root) for path in status_paths],
        "diagnostics": archive_diagnostics(
            args.diagnostics_root.resolve(), args.diagnostics_glob
        )
        if args.diagnostics_root is not None
        else [],
    }
    encoded = (json.dumps(archive, indent=2, sort_keys=True) + "\n").encode("utf-8")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    temporary.write_bytes(encoded)
    temporary.replace(output)
    print(
        json.dumps(
            {
                "output": str(output),
                "bytes": len(encoded),
                "sha256": sha256_bytes(encoded),
                "run_count": len(status_paths),
                "diagnostic_count": len(archive["diagnostics"]),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
