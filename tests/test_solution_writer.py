#!/usr/bin/env python3

import concurrent.futures
import json
import queue
import tempfile
from pathlib import Path
import sys
import threading
import time
import unittest
from unittest import mock


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

from run_experiment import (  # noqa: E402
    CompetitionTimeout,
    StreamingSerialEvaluation,
    apply_fallback_schedule_profile,
    code2_completed_within_limit,
    code2_time_limit,
    contiguous_shard_groups,
    effective_process_timeout,
    fast_screen_affinity_groups,
    longest_first_contingencies,
    progress_checkpoint_due,
    streamed_queue_get,
    validate_and_normalize_evaluation_details,
    read_contingency_blocks,
    write_contingency_subset,
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
    def test_progress_checkpoint_is_throttled(self) -> None:
        self.assertFalse(progress_checkpoint_due(10.0, 10.999))
        self.assertTrue(progress_checkpoint_due(10.0, 11.0))
        with self.assertRaises(ValueError):
            progress_checkpoint_due(10.0, 9.0)

    def test_parallel_evaluation_requires_and_normalizes_every_detail(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            internal = output / "internal"
            details = {
                "BASECASE": {"obj": {"val": 10.0}, "infeas": {"val": False}},
                "CTG_A": {"obj": {"val": 4.0}, "infeas": {"val": False}},
                "CTG_B": {"obj": {"val": 8.0}, "infeas": {"val": False}},
            }
            for label, detail in details.items():
                write_json(output / f"eval_detail_{label}.json", detail)
            write_json(
                output / "eval_summary.json",
                {
                    "num_ctg": 2,
                    "obj": 16.0,
                    "infeas": 0.0,
                    "obj_cumulative": 99.0,
                    "obj_all_cases": {"BASECASE": 8.0},
                    "infeas_cumulative": 0.0,
                    "infeas_all_cases": {"BASECASE": 0},
                },
            )
            summary = json.loads((output / "eval_summary.json").read_text())
            normalized, certificate = validate_and_normalize_evaluation_details(
                output, internal, {"CTG_A", "CTG_B"}, summary, 3
            )

            self.assertEqual(normalized["obj_cumulative"], 16.0)
            self.assertEqual(
                normalized["obj_all_cases"],
                {"BASECASE": 10.0, "CTG_A": 4.0, "CTG_B": 8.0},
            )
            self.assertEqual(len(normalized["infeas_all_cases"]), 3)
            self.assertTrue(certificate["complete_label_set"])
            self.assertTrue((internal / "eval_summary.vendor_mpi.json").exists())

    def test_parallel_evaluation_rejects_missing_detail(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            write_json(
                output / "eval_detail_BASECASE.json",
                {"obj": {"val": 10.0}, "infeas": {"val": False}},
            )
            with self.assertRaisesRegex(RuntimeError, "detail-set mismatch"):
                validate_and_normalize_evaluation_details(
                    output,
                    output / "internal",
                    {"CTG_A"},
                    {"num_ctg": 1, "obj": 10.0, "infeas": 0.0},
                    2,
                )

    def test_contingency_subset_is_exact_and_terminated(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.con"
            destination = root / "subset.con"
            source.write_text(
                "CONTINGENCY CTG_B\n"
                "OPEN BRANCH FROM BUS 1 TO BUS 2 CIRCUIT 1\n"
                "END\n"
                "CONTINGENCY CTG_A\n"
                "REMOVE UNIT 1 FROM BUS 3\n"
                "END\n"
                "END\n",
                encoding="utf-8",
            )
            write_contingency_subset(source, destination, ["CTG_A"])
            self.assertEqual(
                destination.read_text(encoding="utf-8"),
                "CONTINGENCY CTG_A\n"
                "REMOVE UNIT 1 FROM BUS 3\n"
                "END\n"
                "END\n",
            )
            self.assertEqual(
                set(read_contingency_blocks(destination)), {"CTG_A"}
            )

    def test_serial_shard_certificate_does_not_repair_mpi_fields(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            internal = output / "internal"
            for label, objective in (("BASECASE", 10.0), ("CTG_A", 4.0)):
                write_json(
                    output / f"eval_detail_{label}.json",
                    {"obj": {"val": objective}, "infeas": {"val": False}},
                )
            summary = {
                "solutions_exist": True,
                "num_ctg": 1,
                "obj": 14.0,
                "infeas": 0.0,
                "obj_cumulative": 14.0,
                "obj_all_cases": {"BASECASE": 10.0, "CTG_A": 4.0},
                "infeas_cumulative": 0.0,
                "infeas_all_cases": {"BASECASE": False, "CTG_A": False},
            }
            _, certificate = validate_and_normalize_evaluation_details(
                output,
                internal,
                {"CTG_A"},
                summary,
                4,
                "serial_shards",
            )
            self.assertEqual(certificate["parallel_mode"], "serial_shards")
            self.assertEqual(
                certificate["repaired_vendor_mpi_bookkeeping_fields"], []
            )
            self.assertFalse(
                (internal / "eval_summary.vendor_mpi.json").exists()
            )

    def test_profiled_schedule_prioritizes_measured_fallbacks(self):
        contingencies = [
            {"label": "A", "schedule_rank": 1},
            {"label": "B", "schedule_rank": 2},
            {"label": "C", "schedule_rank": 3},
        ]
        profile = {
            "B": {"worker": 1, "predicted_wall_seconds": 2.0},
            "C": {"worker": 0, "predicted_wall_seconds": 5.0},
        }
        scheduled = apply_fallback_schedule_profile(contingencies, profile)
        self.assertEqual([item["label"] for item in scheduled], ["C", "B", "A"])
        self.assertEqual([item["schedule_rank"] for item in scheduled], [1, 2, 3])
        self.assertEqual(scheduled[0]["profiled_fallback_worker"], 0)

    def test_streamed_fallback_waits_for_screening_and_then_drains(self):
        tasks = queue.Queue()
        screening_finished = threading.Event()
        abort = threading.Event()
        waiter_started = threading.Event()

        def wait_for_task():
            waiter_started.set()
            return streamed_queue_get(
                tasks, screening_finished, abort, poll_seconds=0.005
            )

        with concurrent.futures.ThreadPoolExecutor(max_workers=1) as pool:
            future = pool.submit(wait_for_task)
            self.assertTrue(waiter_started.wait(timeout=1.0))
            with self.assertRaises(concurrent.futures.TimeoutError):
                future.result(timeout=0.02)
            expected = {"label": "CTG_TEST"}
            tasks.put(expected)
            self.assertEqual(future.result(timeout=1.0), expected)

        screening_finished.set()
        self.assertIsNone(
            streamed_queue_get(
                tasks, screening_finished, abort, poll_seconds=0.005
            )
        )

        profiled = queue.PriorityQueue()
        profiled.put((-2.0, 2, "SHORT", {"label": "SHORT"}))
        profiled.put((-5.0, 1, "LONG", {"label": "LONG"}))
        shared = queue.Queue()
        shared.put({"label": "SHARED"})
        self.assertEqual(
            streamed_queue_get(
                shared,
                screening_finished,
                abort,
                poll_seconds=0.005,
                profiled_queue=profiled,
            )["label"],
            "LONG",
        )
        self.assertEqual(
            streamed_queue_get(
                shared,
                screening_finished,
                abort,
                poll_seconds=0.005,
                profiled_queue=profiled,
            )["label"],
            "SHORT",
        )
        self.assertEqual(
            streamed_queue_get(
                shared,
                screening_finished,
                abort,
                poll_seconds=0.005,
                profiled_queue=profiled,
            )["label"],
            "SHARED",
        )

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

    def test_fast_screen_affinity_groups_keep_related_outages_together(self):
        case = {
            "gen": {
                "1": {"index": 1, "gen_bus": 10},
                "2": {"index": 2, "gen_bus": 10},
                "3": {"index": 3, "gen_bus": 20},
            },
            "branch": {
                "1": {"index": 1, "f_bus": 30},
                "2": {"index": 2, "f_bus": 30},
                "3": {"index": 3, "f_bus": 40},
            },
        }
        state = {
            "pg": [0.8, 0.2, 0.6],
            "qg": [0.0, 0.0, 0.0],
            "pf": [0.7, 0.1, 0.4],
            "qf": [0.0, 0.0, 0.0],
            "pt": [-0.7, -0.1, -0.4],
            "qt": [0.0, 0.0, 0.0],
        }
        records = [
            {"label": "G1", "type": "gen", "idx": 1, "schedule_rank": 1},
            {"label": "B1", "type": "branch", "idx": 1, "schedule_rank": 2},
            {"label": "G2", "type": "gen", "idx": 2, "schedule_rank": 3},
            {"label": "B2", "type": "branch", "idx": 2, "schedule_rank": 4},
            {"label": "G3", "type": "gen", "idx": 3, "schedule_rank": 5},
            {"label": "B3", "type": "branch", "idx": 3, "schedule_rank": 6},
        ]
        groups = fast_screen_affinity_groups(case, state, records)
        labels = [[item["label"] for item in group] for group in groups]
        self.assertIn(["G2", "G1"], labels)
        self.assertIn(["B2", "B1"], labels)
        self.assertIn(["G3"], labels)
        self.assertIn(["B3"], labels)
        flattened = [item for group in groups for item in group]
        self.assertEqual(len(flattened), len(records))
        self.assertEqual(
            {item["label"] for item in flattened},
            {item["label"] for item in records},
        )
        for group in groups:
            self.assertEqual(
                [item["fast_screen_affinity_position"] for item in group],
                list(range(1, len(group) + 1)),
            )
            self.assertTrue(
                all(item["fast_screen_affinity_size"] == len(group) for item in group)
            )
        easy_first_groups = fast_screen_affinity_groups(
            case, state, records, difficult_groups_first=False
        )
        self.assertLessEqual(
            easy_first_groups[0][0]["fast_screen_affinity_score"],
            easy_first_groups[-1][0]["fast_screen_affinity_score"],
        )

    def test_contiguous_evaluation_shards_preserve_schedule(self):
        labels = [f"CTG_{index}" for index in range(7)]
        groups = contiguous_shard_groups(labels, 3)
        self.assertEqual([len(group) for group in groups], [3, 2, 2])
        self.assertEqual([label for group in groups for label in group], labels)

    def test_streaming_evaluation_runs_ready_shards_and_merges_exactly(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            case_dir = root / "case"
            output_dir = root / "output"
            internal_dir = output_dir / "internal"
            case_dir.mkdir()
            output_dir.mkdir()
            (case_dir / "case.raw").write_text("raw\n", encoding="utf-8")
            (case_dir / "case.json").write_text("{}\n", encoding="utf-8")
            labels = ["CTG_A", "CTG_B", "CTG_C", "CTG_D"]
            (case_dir / "case.con").write_text(
                "".join(
                    f"CONTINGENCY {label}\nEND\n" for label in labels
                )
                + "END\n",
                encoding="utf-8",
            )
            (output_dir / "solution_BASECASE.txt").write_text(
                "base\n", encoding="utf-8"
            )
            evaluator = root / "fake_evaluator.py"
            evaluator.write_text(
                """
import json
from pathlib import Path
import sys

case_dir = Path(sys.argv[2])
solution_dir = Path(sys.argv[3])
labels = []
for line in (case_dir / "case.con").read_text(encoding="utf-8").splitlines():
    if line.startswith("CONTINGENCY "):
        labels.append(line.split()[1])

details = {"BASECASE": {"obj": {"val": 10.0}, "infeas": {"val": False}}}
for index, label in enumerate(labels, start=1):
    details[label] = {"obj": {"val": float(index)}, "infeas": {"val": False}}
for label, detail in details.items():
    (solution_dir / f"eval_detail_{label}.json").write_text(
        json.dumps(detail), encoding="utf-8"
    )
objective = 10.0 + sum(details[label]["obj"]["val"] for label in labels) / len(labels)
(solution_dir / "eval_summary.json").write_text(
    json.dumps({
        "solutions_exist": True,
        "num_ctg": len(labels),
        "obj": objective,
        "infeas": 0.0,
    }),
    encoding="utf-8",
)
""".lstrip(),
                encoding="utf-8",
            )
            manager = StreamingSerialEvaluation(
                Path(sys.executable),
                evaluator,
                case_dir,
                output_dir,
                internal_dir,
                labels,
                2,
                1,
                time.perf_counter() + 30.0,
            )
            for label in labels:
                (output_dir / f"solution_{label}.txt").write_text(
                    f"{label}\n", encoding="utf-8"
                )
                manager.mark_completed(label)
            summary, metadata = manager.finish()

            self.assertEqual(summary["num_ctg"], 4)
            self.assertEqual(summary["infeas"], 0.0)
            self.assertEqual(metadata["mode"], "overlapped_serial_shards")
            self.assertEqual(metadata["active_shards"], 2)
            self.assertEqual(metadata["maximum_parallel_processes"], 1)
            self.assertEqual(
                {
                    path.stem.removeprefix("eval_detail_")
                    for path in output_dir.glob("eval_detail_*.json")
                },
                {"BASECASE", *labels},
            )

    def test_code2_budget_uses_contingency_count(self):
        self.assertEqual(code2_time_limit(105, 2.0), 210.0)
        self.assertEqual(code2_time_limit(0, 2.0), 0.0)

    def test_process_timeout_is_bounded_by_stage_deadline(self):
        self.assertEqual(effective_process_timeout(300.0, 120.0, now=100.0), 20.0)
        self.assertEqual(effective_process_timeout(5.0, 120.0, now=100.0), 5.0)

    def test_code2_timeout_is_not_reported_within_limit(self):
        self.assertFalse(code2_completed_within_limit(475.9, 476.0, True))
        self.assertTrue(code2_completed_within_limit(475.9, 476.0, False))

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
