from __future__ import annotations

import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from read_run_progress import read_snapshot_bytes, read_progress


class ProgressSnapshotTests(unittest.TestCase):
    def test_closed_reader_allows_replacement_and_retains_old_snapshot(self):
        with tempfile.TemporaryDirectory() as directory:
            status = Path(directory) / "run_status.json"
            temporary = Path(directory) / "next.tmp"
            status.write_text('{"version": 1}', encoding="utf-8")
            temporary.write_text('{"version": 2}', encoding="utf-8")
            snapshot = read_snapshot_bytes(status)
            temporary.replace(status)
            self.assertEqual(json.loads(snapshot), {"version": 1})
            self.assertEqual(json.loads(status.read_text()), {"version": 2})

    def test_read_progress_closes_checkpoint_before_decoding(self):
        with tempfile.TemporaryDirectory() as directory:
            status = Path(directory) / "run_status.json"
            temporary = Path(directory) / "next.tmp"
            status.write_text('{"stage": "code1"}', encoding="utf-8")
            temporary.write_text('{"stage": "code2"}', encoding="utf-8")
            original_loads = json.loads

            def replace_before_decoding(snapshot):
                temporary.replace(status)
                return original_loads(snapshot)

            with mock.patch("read_run_progress.json.loads", replace_before_decoding):
                self.assertEqual(read_progress(status)["stage"], "code1")
            self.assertEqual(json.loads(status.read_text())["stage"], "code2")

    @unittest.skipUnless(os.name == "nt", "Windows sharing-mode diagnostic")
    def test_ordinary_reader_can_reproduce_access_denied_on_replace(self):
        with tempfile.TemporaryDirectory() as directory:
            status = Path(directory) / "run_status.json"
            temporary = Path(directory) / "next.tmp"
            status.write_text('{"version": 1}', encoding="utf-8")
            temporary.write_text('{"version": 2}', encoding="utf-8")
            with status.open("r", encoding="utf-8"):
                with self.assertRaises(PermissionError) as failure:
                    temporary.replace(status)
                self.assertIn(failure.exception.winerror, (5, 32))
            # The same identity/path succeeds after closing the reader. This
            # demonstrates sharing contention, not an ACL change or bypass.
            temporary.replace(status)

    def test_compact_snapshot_omits_schedule_and_preserves_unknowns(self):
        with tempfile.TemporaryDirectory() as directory:
            status = Path(directory) / "run_status.json"
            status.write_text(json.dumps({
                "stage": "code2", "completed_contingency_count": 7,
                "competition_timing": {"contingency_count": 10},
                "contingency_schedule": [{"label": "do not emit"}] * 100,
            }), encoding="utf-8")
            summary = read_progress(status)
            self.assertEqual(summary["completed_contingency_count"], 7)
            self.assertEqual(summary["expected_contingency_count"], 10)
            self.assertNotIn("contingency_schedule", summary)
            self.assertNotIn("total_wall_seconds", summary)
            self.assertTrue(summary["snapshot_is_not_process_liveness_evidence"])

    def test_missing_file_is_not_hidden_or_retried(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(FileNotFoundError):
                read_progress(Path(directory) / "missing.json")

    def test_onedrive_path_is_rejected(self):
        with self.assertRaises(ValueError):
            read_progress(Path("OneDrive") / "status.json")


if __name__ == "__main__":
    unittest.main()
