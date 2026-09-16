"""Read a compact checkpoint with no file handle held during parsing/display.

Use this for live monitoring instead of streaming run_status.json into a slow
console pipeline. The short byte-copy may still briefly share the file; the
writer's existing bounded transient-lock retry policy is unchanged.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from archive_run_provenance import reject_onedrive


def read_snapshot_bytes(path: Path) -> bytes:
    reject_onedrive(path)
    # Return only bytes, never a live stream: JSON decoding and output may be
    # slow, but neither should retain a Windows handle on the checkpoint.
    with path.open("rb") as stream:
        return stream.read()


def read_progress(path: Path) -> dict:
    snapshot = read_snapshot_bytes(path)
    status = json.loads(snapshot)
    if not isinstance(status, dict):
        raise ValueError("run status must be an object")
    # Do not emit the very large contingency schedule or solution vectors.
    names = ("stage", "success", "git_revision", "started_at_utc",
             "base_algorithm_wall_seconds", "completed_contingency_count",
             "screened_contingency_count", "last_completed_contingency",
             "total_wall_seconds", "code2_timed_out", "error")
    result = {name: status[name] for name in names if name in status}
    result["expected_contingency_count"] = status.get(
        "competition_timing", {}).get("contingency_count")
    if status.get("started_at_utc"):
        start = datetime.fromisoformat(status["started_at_utc"])
        if start.tzinfo is None:
            raise ValueError("start timestamp must include a timezone")
        result["seconds_since_start_at_observation"] = (
            datetime.now(timezone.utc) - start).total_seconds()
    result["snapshot_is_not_process_liveness_evidence"] = True
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("status", type=Path)
    args = parser.parse_args()
    print(json.dumps(read_progress(args.status), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
