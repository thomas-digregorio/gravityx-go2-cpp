from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "archive_run_provenance.py"


def test_archives_compact_hash_backed_run_provenance(tmp_path: Path) -> None:
    runs = tmp_path / "runs"
    run = runs / "c2fen19402_deadbee_20260826" / "scenario"
    diagnostics = runs / "diagnostics" / "c2fen19402_example"
    run.mkdir(parents=True)
    diagnostics.mkdir(parents=True)
    large_schedule = [{"label": f"CTG_{index:06d}"} for index in range(5000)]
    status = {
        "success": False,
        "stage": "code2",
        "git_revision": "deadbeef",
        "contingency_schedule": large_schedule,
    }
    (run / "run_status.json").write_text(json.dumps(status), encoding="utf-8")
    (run / "run_summary.json").write_text(
        json.dumps({"completed": 12}), encoding="utf-8"
    )
    (run / "solver.err").write_text("failure tail", encoding="utf-8")
    (run / "solution_CTG_000001.txt").write_text("generated", encoding="utf-8")
    (diagnostics / "trace.txt").write_text("trace", encoding="utf-8")
    output = tmp_path / "evidence" / "archive.json"

    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--runs-root",
            str(runs),
            "--output",
            str(output),
            "--diagnostics-root",
            str(runs / "diagnostics"),
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    reported = json.loads(completed.stdout)
    archived = json.loads(output.read_text(encoding="utf-8"))
    assert reported["run_count"] == 1
    assert reported["diagnostic_count"] == 1
    assert archived["runs"][0]["run_status"]["success"] is False
    schedule = archived["runs"][0]["run_status"]["contingency_schedule"]
    assert schedule["omitted"] is True
    assert schedule["count"] == 5000
    assert len(schedule["sha256"]) == 64
    assert archived["runs"][0]["inventory"]["solution_file_entries"] == 1
    assert archived["runs"][0]["root_log_tails"][0]["tail"] == "failure tail"
