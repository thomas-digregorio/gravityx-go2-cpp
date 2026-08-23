#!/usr/bin/env python3

import tempfile
from pathlib import Path
import sys
import unittest


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

from run_experiment import write_solution  # noqa: E402


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


if __name__ == "__main__":
    unittest.main()
