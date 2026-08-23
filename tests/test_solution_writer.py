#!/usr/bin/env python3

import json
import tempfile
from pathlib import Path
import sys
import unittest
from unittest import mock


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

from run_experiment import (  # noqa: E402
    CompetitionTimeout,
    code2_time_limit,
    effective_process_timeout,
    longest_first_contingencies,
    write_json,
    write_solution,
)


class SolutionWriterTests(unittest.TestCase):
    def setUp(self):
        self.case = {
            "bus": {
                "1": {"index": 1, "present": True},
                "2": {"index": 2, "present": True},
            },
            "load": {
                "1": {"index": 1, "present": True, "source_id": ["load", 2, "L1"]},
            },
            "gen": {
                "1": {"index": 1, "present": True, "source_id": ["generator", 1, "G1"]},
            },
            "branch": {
                "1": {
                    "index": 1,
                    "present": True,
                    "transformer": False,
                    "source_id": ["branch", 1, 2, "1"],
                    "br_status": 1,
                },
                "2": {
                    "index": 2,
                    "present": True,
                    "transformer": True,
                    "source_id": ["transformer", 1, 2, 0, "T1"],
                    "br_status": 1,
                    "control_mode": 1,
                    "tm_step": -2,
                },
            },
            "shunt": {
                "1": {
                    "index": 1,
                    "present": True,
                    "dispatchable": True,
                    "source_id": ["switched shunt", 2, 0],
                    "xst": [1, 0],
                },
            },
        }
        self.state = {
            "vm": [1.0, 0.99],
            "va": [0.0, -0.01],
            "demand_factor": [1.0],
            "pg": [0.5],
            "qg": [0.1],
        }

    def render(self, contingency=None):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "solution.txt"
            write_solution(path, self.case, self.state, [1], contingency)
            return path.read_text(encoding="utf-8")

    def test_base_sections_and_transformer_step(self):
        text = self.render()
        self.assertIn("1, 1, 0", text)
        self.assertIn("2, L1, 1", text)
        self.assertIn("1, G1, 0.5, 0.10000000000000001, 1", text)
        self.assertIn("1, 2, 1, 1", text)
        self.assertIn("1, 2, T1, 1, -2", text)
        self.assertIn("2, 1, 0", text)

    def test_generator_outage_is_omitted(self):
        text = self.render({"type": "gen", "idx": 1})
        self.assertNotIn("1, G1,", text)

    def test_branch_outage_is_omitted(self):
        text = self.render({"type": "branch", "idx": 1})
        self.assertNotIn("1, 2, 1, 1", text)
        self.assertIn("1, 2, T1, 1, -2", text)


class CompetitionTimingTests(unittest.TestCase):
    def test_longest_first_schedule_uses_base_apparent_power(self):
        case = {
            "gen": {
                "1": {"index": 1},
                "2": {"index": 2},
            },
            "branch": {
                "1": {"index": 1},
                "2": {"index": 2},
            },
        }
        state = {
            "pg": [0.3, 0.8],
            "qg": [0.4, 0.6],
            "pf": [0.0, 0.6],
            "qf": [0.2, 0.8],
            "pt": [0.1, -0.5],
            "qt": [0.0, -0.5],
        }
        records = [
            {"label": "G_LOW", "type": "gen", "idx": 1},
            {"label": "B_HIGH", "type": "branch", "idx": 2},
            {"label": "G_HIGH", "type": "gen", "idx": 2},
        ]
        scheduled = longest_first_contingencies(case, state, records)
        self.assertEqual(
            [item["label"] for item in scheduled],
            ["B_HIGH", "G_HIGH", "G_LOW"],
        )
        self.assertEqual([item["schedule_rank"] for item in scheduled], [1, 2, 3])
        self.assertAlmostEqual(
            scheduled[0]["schedule_score_base_apparent_power"], 1.0
        )

    def test_code2_budget_uses_contingency_count(self):
        self.assertEqual(code2_time_limit(105, 2.0), 210.0)
        self.assertEqual(code2_time_limit(0, 2.0), 0.0)

    def test_process_timeout_is_bounded_by_stage_deadline(self):
        self.assertEqual(effective_process_timeout(300.0, 120.0, now=100.0), 20.0)
        self.assertEqual(effective_process_timeout(5.0, 120.0, now=100.0), 5.0)

    def test_expired_stage_is_rejected_before_launch(self):
        with self.assertRaises(CompetitionTimeout):
            effective_process_timeout(5.0, 100.0, now=100.0)

    def test_json_checkpoint_retries_transient_windows_lock(self):
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / "status.json"
            original_replace = Path.replace
            attempts = 0

            def flaky_replace(source, target):
                nonlocal attempts
                attempts += 1
                if attempts < 3:
                    raise PermissionError("simulated transient reader lock")
                return original_replace(source, target)

            with mock.patch.object(Path, "replace", flaky_replace):
                write_json(destination, {"success": True})

            self.assertEqual(attempts, 3)
            self.assertEqual(json.loads(destination.read_text()), {"success": True})


if __name__ == "__main__":
    unittest.main()
