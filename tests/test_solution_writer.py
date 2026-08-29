#!/usr/bin/env python3

import concurrent.futures
import json
import numpy as np
import queue
import subprocess
import tempfile
from pathlib import Path
import sys
import threading
import time
from types import SimpleNamespace
import unittest
from unittest import mock


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

from run_experiment import (  # noqa: E402
    AffinityScreenWorkQueue,
    CompetitionTimeout,
    StreamingSerialEvaluation,
    apply_fallback_schedule_profile,
    code2_completed_within_limit,
    code2_time_limit,
    completion_order_shard_groups,
    contiguous_shard_groups,
    cpp_command,
    effective_process_timeout,
    evaluator_priority_popen_options,
    evaluator_subprocess_environment,
    fast_screen_affinity_groups,
    finalization_reserve_seconds,
    idle_screen_evaluation_process_target,
    load_fast_screen_heavy_profile,
    longest_first_contingencies,
    progress_log_due,
    progress_checkpoint_due,
    persistent_evaluator_protocol_markers_present,
    require_minimum_free_space,
    stage_wsl_worker_inputs,
    streamed_queue_get,
    validate_and_normalize_evaluation_details,
    validate_in_memory_evaluation_certificate,
    read_contingency_blocks,
    write_contingency_subset,
    write_json,
    write_solution,
)
from fast_official_evaluator import (  # noqa: E402
    install_compact_output_mode,
    install_in_memory_summary_aggregation,
    install_static_case_read_cache,
    read_generated_solution,
    suppress_verbose_timing_output,
    skip_dataframe_copy,
)
from build_fast_screen_heavy_profile import (  # noqa: E402
    existing_profile_measurements,
    merge_max_measurements,
    parse_worker_log_measurements,
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
    def test_cpp_command_passes_explicit_gravityx_environment_into_wsl(self):
        command = cpp_command(
            Path("C:/repo/build/gravityx"),
            "Ubuntu-24.04",
            ["contingency-worker"],
            12.0,
            {
                "GRAVITYX_ECONOMIC_EXACT_REUSE_NUMERIC": "1",
                "GRAVITYX_ECONOMIC_EXACT_TARGET_FRACTION": "0.0125",
            },
        )
        timeout_position = command.index("timeout")
        self.assertIn(
            "GRAVITYX_ECONOMIC_EXACT_REUSE_NUMERIC=1",
            command[:timeout_position],
        )
        self.assertIn(
            "GRAVITYX_ECONOMIC_EXACT_TARGET_FRACTION=0.0125",
            command[:timeout_position],
        )

    def test_cpp_command_rejects_non_gravityx_environment(self):
        with self.assertRaisesRegex(ValueError, "GRAVITYX_ prefix"):
            cpp_command(
                Path("C:/repo/build/gravityx"),
                "Ubuntu-24.04",
                ["solve"],
                12.0,
                {"LD_LIBRARY_PATH": "unexpected"},
            )

    def test_compact_summary_uses_measured_short_finalization_reserve(self) -> None:
        self.assertEqual(finalization_reserve_seconds(False), 5.0)
        self.assertEqual(finalization_reserve_seconds(True), 0.5)
        self.assertEqual(finalization_reserve_seconds(True, 4.0), 4.0)
        with self.assertRaisesRegex(ValueError, "positive and finite"):
            finalization_reserve_seconds(True, 0.0)

    def test_progress_checkpoint_is_throttled(self) -> None:
        self.assertFalse(progress_checkpoint_due(10.0, 14.999))
        self.assertTrue(progress_checkpoint_due(10.0, 15.0))
        self.assertTrue(progress_checkpoint_due(10.0, 11.0, 1.0))
        with self.assertRaises(ValueError):
            progress_checkpoint_due(10.0, 9.0)

    def test_progress_logging_is_batched_but_includes_boundaries(self) -> None:
        self.assertTrue(progress_log_due(1, 6693))
        self.assertTrue(progress_log_due(10, 6693))
        self.assertFalse(progress_log_due(11, 6693))
        self.assertTrue(progress_log_due(100, 6693))
        self.assertTrue(progress_log_due(6693, 6693))
        with self.assertRaises(ValueError):
            progress_log_due(0, 6693)

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

    def test_affinity_screen_queue_prioritizes_failed_group_siblings(self):
        abort = threading.Event()
        first = [{"label": "A"}, {"label": "B"}, {"label": "C"}]
        second = [{"label": "D"}]
        work = AffinityScreenWorkQueue([first, second], worker_count=2)

        source, group = work.get(abort, worker_id=0)
        self.assertEqual(source, "bulk")
        self.assertEqual([item["label"] for item in group], ["A", "B", "C"])
        self.assertEqual(
            work.requeue_remaining_as_singletons(
                group,
                1,
                origin_worker_id=0,
                source=source,
            ),
            2,
        )
        work.task_done(source)

        source, group = work.get(abort, worker_id=0)
        self.assertEqual((source, group[0]["label"]), ("bulk", "D"))
        work.task_done(source)

        observed = []
        for _ in range(2):
            source, group = work.get(abort, worker_id=1)
            observed.append((source, [item["label"] for item in group]))
            work.task_done(source)
        self.assertTrue(all(source == "urgent_bulk" for source, _ in observed))
        self.assertEqual(
            {labels[0] for _, labels in observed},
            {"B", "C"},
        )
        self.assertEqual(work.remaining_group_count, 0)
        self.assertIsNone(work.get(abort, worker_id=0))

    def test_affinity_screen_queue_caps_profiled_heavy_lane(self):
        abort = threading.Event()
        groups = [
            [{"label": "H1"}],
            [{"label": "B1"}],
            [{"label": "H2"}],
            [{"label": "B2"}],
        ]
        work = AffinityScreenWorkQueue(
            groups,
            worker_count=2,
            heavy_labels={"H1", "H2"},
            heavy_worker_count=1,
        )
        self.assertEqual(work.initial_heavy_group_count, 2)
        self.assertEqual(work.worker_lane(0), "heavy")
        self.assertEqual(work.worker_lane(1), "bulk")

        heavy_source, heavy_group = work.get(abort, worker_id=0)
        bulk_source, bulk_group = work.get(abort, worker_id=1)
        self.assertEqual(
            (heavy_source, heavy_group[0]["label"]), ("heavy", "H1")
        )
        self.assertEqual(
            (bulk_source, bulk_group[0]["label"]), ("bulk", "B1")
        )
        work.task_done(heavy_source)
        work.task_done(bulk_source)

        heavy_source, heavy_group = work.get(abort, worker_id=0)
        bulk_source, bulk_group = work.get(abort, worker_id=1)
        self.assertEqual(
            (heavy_source, heavy_group[0]["label"]), ("heavy", "H2")
        )
        self.assertEqual(
            (bulk_source, bulk_group[0]["label"]), ("bulk", "B2")
        )
        work.task_done(heavy_source)
        work.task_done(bulk_source)
        self.assertEqual(work.remaining_group_count, 0)

    def test_idle_bulk_worker_spills_when_queued_bulk_lane_drains(self):
        abort = threading.Event()
        work = AffinityScreenWorkQueue(
            [
                [{"label": "H1"}],
                [{"label": "H2"}],
                [{"label": "B1"}],
            ],
            worker_count=3,
            heavy_labels={"H1", "H2"},
            heavy_worker_count=1,
        )

        heavy_source, heavy_group = work.get(abort, worker_id=0)
        bulk_source, bulk_group = work.get(abort, worker_id=1)
        self.assertEqual((heavy_source, heavy_group[0]["label"]), ("heavy", "H1"))
        self.assertEqual((bulk_source, bulk_group[0]["label"]), ("bulk", "B1"))

        # B1 and H1 both remain active. With no queued bulk work left, the
        # second bulk worker must help the heavy critical path immediately
        # instead of waiting for every active bulk group to finish.
        spill_source, spill_group = work.get(abort, worker_id=2)
        self.assertEqual((spill_source, spill_group[0]["label"]), ("heavy", "H2"))

        work.task_done(spill_source)
        work.task_done(bulk_source)
        work.task_done(heavy_source)
        self.assertEqual(work.remaining_group_count, 0)

    def test_active_heavy_group_exposes_siblings_after_bulk_queue_drains(self):
        abort = threading.Event()
        work = AffinityScreenWorkQueue(
            [
                [{"label": "H1"}, {"label": "H2"}, {"label": "H3"}],
                [{"label": "B1"}],
                [{"label": "B2"}],
            ],
            worker_count=3,
            heavy_labels={"H1"},
            heavy_worker_count=1,
        )

        heavy_source, heavy_group = work.get(abort, worker_id=0)
        bulk_source_1, bulk_group_1 = work.get(abort, worker_id=1)
        self.assertFalse(work.should_split_heavy_group(heavy_source))
        bulk_source_2, bulk_group_2 = work.get(abort, worker_id=2)
        self.assertTrue(work.should_split_heavy_group(heavy_source))

        self.assertEqual(
            work.requeue_remaining_as_singletons(
                heavy_group,
                1,
                origin_worker_id=0,
                source=heavy_source,
            ),
            2,
        )
        work.task_done(heavy_source)

        exposed = []
        for worker_id in (1, 2):
            source, group = work.get(abort, worker_id=worker_id)
            exposed.append((source, group[0]["label"]))
            work.task_done(source)
        self.assertEqual(
            set(exposed),
            {("urgent_heavy", "H2"), ("urgent_heavy", "H3")},
        )

        work.task_done(bulk_source_1)
        work.task_done(bulk_source_2)
        self.assertEqual(
            [bulk_group_1[0]["label"], bulk_group_2[0]["label"]],
            ["B1", "B2"],
        )
        self.assertEqual(work.remaining_group_count, 0)

    def test_heavy_groups_use_longest_predicted_work_first(self):
        abort = threading.Event()
        work = AffinityScreenWorkQueue(
            [
                [{"label": "H_SHORT"}],
                [{"label": "H_LONG"}, {"label": "H_LONG_SIBLING"}],
                [{"label": "B1"}],
            ],
            worker_count=2,
            heavy_labels={"H_SHORT", "H_LONG"},
            heavy_worker_count=1,
            heavy_label_seconds={"H_SHORT": 6.0, "H_LONG": 7.0},
        )

        source, group = work.get(abort, worker_id=0)
        self.assertEqual(source, "heavy")
        self.assertEqual(
            [item["label"] for item in group],
            ["H_LONG", "H_LONG_SIBLING"],
        )
        self.assertEqual(
            group[0]["fast_screen_profiled_group_wall_seconds"], 13.0
        )
        self.assertEqual(
            group[0]["fast_screen_profiled_heavy_dispatch_rank"], 1
        )
        work.task_done(source)

        source, group = work.get(abort, worker_id=0)
        self.assertEqual((source, group[0]["label"]), ("heavy", "H_SHORT"))
        self.assertEqual(
            group[0]["fast_screen_profiled_heavy_dispatch_rank"], 2
        )
        work.task_done(source)
        source, group = work.get(abort, worker_id=1)
        self.assertEqual((source, group[0]["label"]), ("bulk", "B1"))
        work.task_done(source)
        self.assertEqual(work.remaining_group_count, 0)

    def test_single_heavy_worker_can_reclaim_its_split_sibling(self):
        abort = threading.Event()
        group = [{"label": "H1"}, {"label": "H2"}]
        work = AffinityScreenWorkQueue(
            [group],
            worker_count=2,
            heavy_labels={"H1"},
            heavy_worker_count=1,
        )
        source, assigned = work.get(abort, worker_id=0)
        self.assertEqual(source, "heavy")
        self.assertEqual(
            work.requeue_remaining_as_singletons(
                assigned,
                1,
                origin_worker_id=0,
                source=source,
            ),
            1,
        )
        work.task_done(source)
        source, assigned = work.get(abort, worker_id=0)
        self.assertEqual((source, assigned[0]["label"]), ("urgent_heavy", "H2"))
        work.task_done(source)
        self.assertEqual(work.remaining_group_count, 0)

    def test_fast_screen_heavy_profile_is_timing_only_and_hash_bound(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "heavy.json"
            write_json(
                path,
                {
                    "schema_version": 1,
                    "case_sha256": "case-hash",
                    "heavy_threshold_seconds": 5.0,
                    "contingencies": [
                        {
                            "label": "A",
                            "measured_solver_wall_seconds": 7.5,
                        },
                        {
                            "label": "B",
                            "measured_solver_wall_seconds": 5.0,
                        },
                    ],
                },
            )
            labels, measured, metadata = load_fast_screen_heavy_profile(
                path, "case-hash", {"A", "B", "C"}
            )
            self.assertEqual(labels, {"A", "B"})
            self.assertEqual(measured, {"A": 7.5, "B": 5.0})
            self.assertEqual(metadata["profiled_contingency_count"], 2)
            self.assertFalse(metadata["uses_prior_solution_state"])
            with self.assertRaisesRegex(ValueError, "case hash"):
                load_fast_screen_heavy_profile(
                    path, "different-hash", {"A", "B", "C"}
                )

    def test_full_screen_profile_classifies_heavy_and_keeps_all_timings(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "full_profile.json"
            write_json(
                path,
                {
                    "schema_version": 2,
                    "case_sha256": "case-hash",
                    "heavy_threshold_seconds": 5.0,
                    "contingencies": [
                        {
                            "label": "HEAVY",
                            "measured_solver_wall_seconds": 7.0,
                        },
                        {
                            "label": "LIGHT",
                            "measured_solver_wall_seconds": 1.25,
                        },
                    ],
                },
            )
            labels, measured, metadata = load_fast_screen_heavy_profile(
                path, "case-hash", {"HEAVY", "LIGHT"}
            )
            self.assertEqual(labels, {"HEAVY"})
            self.assertEqual(measured, {"HEAVY": 7.0, "LIGHT": 1.25})
            self.assertEqual(metadata["profiled_contingency_count"], 2)
            self.assertEqual(metadata["profiled_heavy_contingency_count"], 1)

            work = AffinityScreenWorkQueue(
                [[{"label": "HEAVY"}, {"label": "LIGHT"}]],
                worker_count=1,
                heavy_labels=labels,
                heavy_worker_count=1,
                heavy_label_seconds=measured,
            )
            _, group = work.get(threading.Event(), worker_id=0)
            self.assertEqual(
                group[0]["fast_screen_profiled_group_wall_seconds"], 8.25
            )

    def test_fast_screen_profile_measurements_parse_and_merge_by_maximum(self):
        with tempfile.TemporaryDirectory() as directory:
            logs = Path(directory)
            payloads = [
                {
                    "label": "A",
                    "result_summary": {"solve": {"wall_seconds": 4.0}},
                },
                {
                    "label": "B",
                    "result_summary": {"solve": {"wall_seconds": 7.5}},
                },
            ]
            (logs / "worker_000.log").write_text(
                "noise\n"
                + "\n".join(
                    "GRAVITYX_TASK_RESULT " + json.dumps(payload)
                    for payload in payloads
                )
                + "\n",
                encoding="utf-8",
            )
            first = parse_worker_log_measurements(logs)
            self.assertEqual(first, {"A": 4.0, "B": 7.5})
            merged = merge_max_measurements(
                [first, {"A": 8.0, "B": 6.0, "C": 5.0}]
            )
            self.assertEqual(merged, {"A": 8.0, "B": 7.5, "C": 5.0})

    def test_existing_fast_screen_profile_supplies_timing_only_measurements(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "profile.json"
            write_json(
                path,
                {
                    "schema_version": 1,
                    "case_sha256": "case-hash",
                    "contingencies": [
                        {
                            "label": "A",
                            "measured_solver_wall_seconds": 7.5,
                        }
                    ],
                },
            )
            self.assertEqual(
                existing_profile_measurements(path, "case-hash"),
                {"A": 7.5},
            )
            with self.assertRaisesRegex(ValueError, "case hash"):
                existing_profile_measurements(path, "different-hash")

    def test_idle_screen_capacity_promotes_evaluation_without_oversubscription(self):
        self.assertEqual(
            idle_screen_evaluation_process_target(1, 44, 24, 24), 1
        )
        self.assertEqual(
            idle_screen_evaluation_process_target(1, 44, 24, 20), 5
        )
        self.assertEqual(
            idle_screen_evaluation_process_target(1, 44, 24, 0), 25
        )
        self.assertEqual(
            idle_screen_evaluation_process_target(4, 8, 24, 0), 8
        )

    def test_contiguous_evaluation_shards_preserve_schedule(self):
        labels = [f"CTG_{index}" for index in range(7)]
        groups = contiguous_shard_groups(labels, 3)
        self.assertEqual([len(group) for group in groups], [3, 2, 2])
        self.assertEqual([label for group in groups for label in group], labels)

    def test_completion_order_shards_can_taper_the_final_tail(self):
        labels = [f"CTG_{index}" for index in range(30)]
        groups = completion_order_shard_groups(
            labels, 6, [6, 4, 2]
        )
        self.assertEqual([len(group) for group in groups], [6, 6, 6, 6, 4, 2])
        self.assertEqual([label for group in groups for label in group], labels)
        with self.assertRaisesRegex(ValueError, "leave at least one"):
            completion_order_shard_groups(labels, 3, [6, 4, 2])

    def test_evaluator_environment_prevents_nested_thread_oversubscription(self):
        with mock.patch.dict(
            "os.environ",
            {"OPENBLAS_NUM_THREADS": "24", "OMP_NUM_THREADS": "12"},
        ):
            environment = evaluator_subprocess_environment(
                1, Path("referenced_vendor_evaluator.py")
            )
        self.assertEqual(environment["OPENBLAS_NUM_THREADS"], "1")
        self.assertEqual(environment["OMP_NUM_THREADS"], "1")
        self.assertEqual(environment["MKL_NUM_THREADS"], "1")
        self.assertTrue(
            environment["GRAVITYX_VENDOR_EVALUATOR"].endswith(
                "referenced_vendor_evaluator.py"
            )
        )

    def test_evaluator_priority_is_opt_in_and_uses_windows_priority_class(self):
        self.assertEqual(evaluator_priority_popen_options(False), {})
        if sys.platform == "win32":
            self.assertEqual(
                evaluator_priority_popen_options(True),
                {
                    "creationflags": int(
                        subprocess.BELOW_NORMAL_PRIORITY_CLASS
                    )
                },
            )
        else:
            with self.assertRaisesRegex(ValueError, "only on Windows"):
                evaluator_priority_popen_options(True)

    def test_persistent_evaluator_requires_all_protocol_markers(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            vendor = root / "vendor.py"
            vendor.write_text("print('ordinary evaluator')\n", encoding="utf-8")
            persistent = root / "persistent.py"
            persistent.write_text(
                "\n".join(
                    (
                        'MODE = "--persistent-server"',
                        'READY = "GRAVITYX_EVALUATOR_READY"',
                        'RESULT = "GRAVITYX_EVALUATION_RESULT"',
                    )
                ),
                encoding="utf-8",
            )

            self.assertFalse(
                persistent_evaluator_protocol_markers_present(vendor)
            )
            self.assertTrue(
                persistent_evaluator_protocol_markers_present(persistent)
            )
            self.assertFalse(
                persistent_evaluator_protocol_markers_present(root / "missing.py")
            )

    def test_fast_evaluator_parser_populates_vendor_arrays_from_canonical_text(self):
        with tempfile.TemporaryDirectory() as directory:
            solution = Path(directory) / "solution.txt"
            solution.write_text(
                """--bus section
i, v, theta
101, 1.01, 0.1
102, 0.99, -0.2
--load section
i, id, t
101, L1, 0.75
--generator section
i, id, p, q, x
101, G1, 1.2, -0.3, 1
--line section
iorig, idest, id, x
101, 102, 1, 1
--transformer section
iorig, idest, id, x, xst
102, 101, T1, 1, -2
--switched shunt section
i, xst1, xst2, xst3, xst4, xst5, xst6, xst7, xst8
102, 1, 2, 3
""",
                encoding="utf-8",
            )
            target = SimpleNamespace(
                num_bus_read=2,
                num_load_read=1,
                num_gen_read=1,
                num_line_read=1,
                num_xfmr_read=1,
                num_swsh_read=1,
                bus_map={101: 0, 102: 1},
                load_map={(101, "L1"): 0},
                gen_map={(101, "G1"): 0},
                line_map={(101, 102, "1"): 0},
                xfmr_map={(102, 101, "T1"): 0},
                swsh_map={102: 0},
                bus_volt_mag=np.zeros(2),
                bus_volt_ang=np.zeros(2),
                load_t=np.zeros(1),
                gen_pow_real=np.zeros(1),
                gen_pow_imag=np.zeros(1),
                gen_xon=np.zeros(1),
                line_xsw=np.zeros(1),
                xfmr_xsw=np.zeros(1),
                xfmr_xst=np.zeros(1),
                swsh_xst=np.zeros((1, 8)),
            )
            read_generated_solution(target, solution)
            skip_dataframe_copy(target)
            np.testing.assert_allclose(target.bus_volt_mag, [1.01, 0.99])
            np.testing.assert_allclose(target.bus_volt_ang, [0.1, -0.2])
            np.testing.assert_allclose(target.load_t, [0.75])
            np.testing.assert_allclose(target.gen_pow_real, [1.2])
            np.testing.assert_allclose(target.gen_pow_imag, [-0.3])
            np.testing.assert_allclose(target.gen_xon, [1.0])
            np.testing.assert_allclose(target.line_xsw, [1.0])
            np.testing.assert_allclose(target.xfmr_xsw, [1.0])
            np.testing.assert_allclose(target.xfmr_xst, [-2.0])
            np.testing.assert_allclose(
                target.swsh_xst, [[1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
            )

            numeric_solution = Path(directory) / "solution_numeric.txt"
            numeric_text = (
                solution.read_text(encoding="utf-8")
                .replace("L1", "1")
                .replace("G1", "2")
                .replace("T1", "3")
            )
            numeric_solution.write_text(numeric_text, encoding="utf-8")
            target.load_map = {(101, "1"): 0}
            target.gen_map = {(101, "2"): 0}
            target.xfmr_map = {(102, 101, "3"): 0}
            read_generated_solution(target, numeric_solution)
            self.assertEqual(
                set(target._gravityx_numeric_layout),
                {"bus", "load", "generator", "line", "transformer", "switched shunt"},
            )
            self.assertEqual(
                set(target._gravityx_exact_binary_matrix_cache),
                {"bus", "load", "generator"},
            )
            self.assertEqual(
                set(target._gravityx_binary_static_section_cache),
                {"line", "transformer"},
            )
            numeric_solution.write_text(
                numeric_text.replace("101, 102, 1, 1", "101, 102, 1, 0"),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "full line rows changed"):
                read_generated_solution(target, numeric_solution)
            numeric_solution.write_text(
                numeric_text.replace("101, 1.01, 0.1", "101, 1.02, 0.1"),
                encoding="utf-8",
            )
            read_generated_solution(target, numeric_solution)
            np.testing.assert_allclose(target.bus_volt_mag, [1.02, 0.99])
            numeric_solution.write_text(
                numeric_text.replace("1.2, -0.3, 1", "1.3, -0.3, 1"),
                encoding="utf-8",
            )
            read_generated_solution(target, numeric_solution)
            np.testing.assert_allclose(target.gen_pow_real, [1.3])

            line_outage_solution = Path(directory) / "solution_line_outage.txt"
            line_outage_solution.write_text(
                numeric_text.replace(
                    "101, 102, 1, 1\n--transformer section",
                    "--transformer section",
                ),
                encoding="utf-8",
            )
            target.num_line_read = 0
            read_generated_solution(target, line_outage_solution)
            np.testing.assert_allclose(target.line_xsw, [0.0])
            np.testing.assert_allclose(target.xfmr_xsw, [1.0])

            transformer_outage_solution = (
                Path(directory) / "solution_transformer_outage.txt"
            )
            transformer_outage_solution.write_text(
                numeric_text.replace(
                    "102, 101, 3, 1, -2\n--switched shunt section",
                    "--switched shunt section",
                ),
                encoding="utf-8",
            )
            target.num_line_read = 1
            target.num_xfmr_read = 0
            read_generated_solution(target, transformer_outage_solution)
            np.testing.assert_allclose(target.line_xsw, [1.0])
            np.testing.assert_allclose(target.xfmr_xsw, [0.0])
            np.testing.assert_allclose(target.xfmr_xst, [0.0])

    def test_fast_evaluator_reuses_exact_written_summaries_in_memory(self):
        class FakeEvaluation:
            def __init__(self):
                self.summary_all_cases = {}
                self.summary = {
                    "obj": {"val": 1.0},
                    "infeas": {"val": False},
                }
                self.ctg_label = ["CTG_A"]
                self.writes = []

            def write_final_summary_and_detail(self, path):
                del path

            def write_detail(
                self, path, case, detail_csv=False, detail_json=False
            ):
                self.writes.append((path, case, detail_csv, detail_json))
                if detail_json:
                    output = Path(path)
                    output.mkdir(parents=True, exist_ok=True)
                    (output / f"eval_detail_{case.strip()}.json").write_text(
                        json.dumps(self.summary), encoding="utf-8"
                    )

            def json_to_summary_all_cases(self, path):
                raise AssertionError(f"unexpected disk reread from {path}")

        module = SimpleNamespace(
            Evaluation=FakeEvaluation,
            clean_string=lambda value: value.strip(),
        )
        install_in_memory_summary_aggregation(module)
        with tempfile.TemporaryDirectory() as directory:
            evaluation = module.Evaluation()
            evaluation.write_detail(directory, "BASECASE", detail_json=True)
            evaluation.summary["obj"]["val"] = 2.0
            evaluation.write_detail(directory, "CTG_A", detail_json=True)
            evaluation.summary["obj"]["val"] = 3.0
            evaluation.json_to_summary_all_cases(directory)
            evaluation.write_final_summary_and_detail(directory)

            self.assertEqual(
                evaluation.writes,
                [
                    (directory, "BASECASE", False, True),
                    (directory, "CTG_A", False, True),
                ],
            )
            self.assertEqual(
                evaluation.summary_all_cases["BASECASE"]["obj"]["val"],
                1.0,
            )
            self.assertEqual(
                evaluation.summary_all_cases["CTG_A"]["obj"]["val"], 2.0
            )
            capture = json.loads(
                (
                    Path(directory)
                    / "gravityx_in_memory_detail_certificate.json"
                ).read_text(encoding="utf-8")
            )
            self.assertEqual(
                capture["objectives"], {"BASECASE": 1.0, "CTG_A": 2.0}
            )

    def test_static_case_cache_restores_pristine_independent_state(self):
        class FakeReader:
            read_count = 0

            def read(self, file_name):
                type(self).read_count += 1
                self.source = str(file_name)
                self.nested = {"values": [type(self).read_count]}

        class FakeSupplemental(FakeReader):
            read_count = 0

        module = SimpleNamespace(
            data=SimpleNamespace(Raw=FakeReader, Sup=FakeSupplemental)
        )
        statistics = install_static_case_read_cache(module)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "case.raw"
            linked = root / "linked.raw"
            source.write_text("case\n", encoding="utf-8")
            linked.hardlink_to(source)

            first = module.data.Raw()
            first.read(source)
            first.nested["values"][0] = 999
            second = module.data.Raw()
            second.read(linked)

            self.assertEqual(FakeReader.read_count, 1)
            self.assertEqual(second.nested, {"values": [1]})
            self.assertIsNot(first.nested, second.nested)
            self.assertEqual(statistics["raw_misses"], 1)
            self.assertEqual(statistics["raw_hits"], 1)

    def test_compact_evaluator_output_keeps_vendor_summary_calculation(self):
        class FakeEvaluation:
            def __init__(self):
                self.summary_all_cases = {}
                self.summary = {
                    "obj": {"val": 1.0},
                    "infeas": {"val": False},
                }
                self.ctg_label = ["CTG_A"]
                self.summary_steps = []

            def write_detail(self, path, case, **kwargs):
                if kwargs.get("detail_json"):
                    output = Path(path)
                    output.mkdir(parents=True, exist_ok=True)
                    (output / f"eval_detail_{case}.json").write_text(
                        json.dumps(self.summary), encoding="utf-8"
                    )

            def write_final_summary_and_detail(self, path):
                raise AssertionError(f"redundant final writers called for {path}")

            def json_to_summary_all_cases(self, path):
                raise AssertionError(f"unexpected disk reread from {path}")

            def summary_all_cases_to_summary(self):
                self.summary_steps.append("calculate")

            def write_summary_json(self, path):
                self.summary_steps.append("write_json")
                (Path(path) / "eval_summary.json").write_text(
                    "{}", encoding="utf-8"
                )

            def write_sol_change(self, path, case):
                raise AssertionError((path, case))

        module = SimpleNamespace(
            Evaluation=FakeEvaluation,
            clean_string=lambda value: value,
        )
        install_in_memory_summary_aggregation(module)
        install_compact_output_mode(module)
        with tempfile.TemporaryDirectory() as directory:
            evaluation = module.Evaluation()
            evaluation.write_detail(directory, "BASECASE", detail_json=True)
            evaluation.summary["obj"]["val"] = 2.0
            evaluation.write_detail(directory, "CTG_A", detail_json=True)
            evaluation.write_final_summary_and_detail(directory)
            evaluation.write_sol_change(directory, "CTG_A")

            self.assertEqual(
                evaluation.summary_steps, ["calculate", "write_json"]
            )
            capture = json.loads(
                (
                    Path(directory)
                    / "gravityx_in_memory_detail_certificate.json"
                ).read_text(encoding="utf-8")
            )
            self.assertEqual(
                capture["final_summary_mode"],
                "vendor_summary_json_without_redundant_aggregate_artifacts",
            )

    def test_in_memory_evaluator_certificate_validates_exact_artifacts(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            values = {"BASECASE": 10.0, "CTG_A": 4.0}
            detail_files = {}
            for label, objective in values.items():
                detail = output / f"eval_detail_{label}.json"
                detail.write_text(
                    json.dumps(
                        {
                            "obj": {"val": objective},
                            "infeas": {"val": False},
                        }
                    ),
                    encoding="utf-8",
                )
                detail_files[label] = {
                    "name": detail.name,
                    "size": detail.stat().st_size,
                }
            write_json(
                output / "gravityx_in_memory_detail_certificate.json",
                {
                    "schema_version": 1,
                    "capture_point": (
                        "immediately_after_unchanged_vendor_write_detail_returned"
                    ),
                    "vendor_detail_json_files_preserved": True,
                    "expected_detail_count": 2,
                    "observed_detail_count": 2,
                    "objectives": values,
                    "infeasibilities": {
                        "BASECASE": False,
                        "CTG_A": False,
                    },
                    "detail_files": detail_files,
                },
            )
            summary = {
                "solutions_exist": True,
                "num_ctg": 1,
                "obj": 14.0,
                "infeas": 0.0,
                "obj_all_cases": values,
                "infeas_all_cases": {
                    "BASECASE": False,
                    "CTG_A": False,
                },
            }
            normalized, certificate = (
                validate_in_memory_evaluation_certificate(
                    output,
                    output / "internal",
                    {"CTG_A"},
                    summary,
                )
            )
            self.assertEqual(normalized["obj"], 14.0)
            self.assertTrue(certificate["complete_label_set"])
            self.assertFalse(certificate["detail_json_files_reparsed"])
            (output / "eval_detail_CTG_A.json").write_text(
                "changed", encoding="utf-8"
            )
            with self.assertRaisesRegex(RuntimeError, "artifact changed"):
                validate_in_memory_evaluation_certificate(
                    output,
                    output / "internal2",
                    {"CTG_A"},
                    summary,
                )

    def test_fast_evaluator_quiet_mode_only_replaces_print(self):
        original = object()
        module = SimpleNamespace(print=original, marker=object())
        marker = module.marker
        suppress_verbose_timing_output(module)
        self.assertIs(module.marker, marker)
        self.assertIsNone(module.print("suppressed"))

    def test_deferred_streaming_evaluation_merges_every_shard_exactly(self):
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
                0,
                time.perf_counter() + 30.0,
                1,
                2,
            )
            for label in labels:
                (output_dir / f"solution_{label}.txt").write_text(
                    f"{label}\n", encoding="utf-8"
                )
                manager.mark_completed(label)
            self.assertEqual(manager.running_records, {})
            self.assertEqual(len(manager.ready_records), 2)
            manager.promote_maximum_processes(1, remaining_screen_groups=1)
            self.assertEqual(manager.maximum_processes, 1)
            self.assertEqual(len(manager.running_records), 1)
            self.assertEqual(len(manager.ready_records), 1)
            summary, metadata = manager.finish()

            self.assertEqual(summary["num_ctg"], 4)
            self.assertEqual(summary["infeas"], 0.0)
            self.assertEqual(metadata["mode"], "overlapped_serial_shards")
            self.assertEqual(metadata["active_shards"], 2)
            self.assertEqual(metadata["maximum_parallel_processes"], 2)
            self.assertFalse(metadata["top_level_details_materialized"])
            self.assertEqual(
                metadata["aggregate_certificate"]["observed_unique_detail_count"],
                5,
            )
            self.assertEqual(metadata["initial_maximum_parallel_processes"], 0)
            self.assertEqual(
                metadata["post_screen_maximum_parallel_processes"], 2
            )
            self.assertEqual(
                metadata["pre_tail_maximum_parallel_processes"], 1
            )
            self.assertEqual(
                metadata["screen_idle_promotion_events"],
                [{"maximum_processes": 1, "remaining_screen_groups": 1}],
            )
            self.assertTrue(metadata["incremental_shard_validation"])
            self.assertEqual(metadata["shard_validation_workers"], 2)
            self.assertEqual(list(output_dir.glob("eval_detail_*.json")), [])
            sharded_detail_labels = {
                path.stem.removeprefix("eval_detail_")
                for path in (
                    internal_dir / "streaming_serial_evaluation_shards"
                ).glob("shard_*/solutions/eval_detail_*.json")
            }
            self.assertEqual(sharded_detail_labels, {"BASECASE", *labels})

    def test_streaming_shards_can_follow_actual_solution_completion_order(self):
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
            completion_order = ["CTG_C", "CTG_A", "CTG_D", "CTG_B"]
            for label in completion_order:
                (output_dir / f"solution_{label}.txt").write_text(
                    f"{label}\n", encoding="utf-8"
                )

            manager = StreamingSerialEvaluation(
                Path(sys.executable),
                root / "unused_evaluator.py",
                case_dir,
                output_dir,
                internal_dir,
                labels,
                2,
                0,
                time.perf_counter() + 30.0,
                1,
                2,
                completion_order_groups=True,
            )
            try:
                for label in completion_order:
                    manager.mark_completed(label)

                self.assertEqual(
                    [record["labels"] for record in manager.records],
                    [["CTG_C", "CTG_A"], ["CTG_D", "CTG_B"]],
                )
                self.assertEqual(len(manager.ready_records), 2)
                self.assertEqual(manager.dynamic_group_labels, [])
                self.assertEqual(manager.dynamic_next_shard_index, 2)
                first_subset = read_contingency_blocks(
                    manager.records[0]["case_dir"] / "case.con"
                )
                self.assertEqual(set(first_subset), {"CTG_C", "CTG_A"})
            finally:
                manager.abort()

    def test_streaming_shards_reuse_persistent_evaluator_process(self):
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
                "".join(f"CONTINGENCY {label}\nEND\n" for label in labels)
                + "END\n",
                encoding="utf-8",
            )
            (output_dir / "solution_BASECASE.txt").write_text(
                "base\n", encoding="utf-8"
            )
            evaluator = root / "persistent_fake_evaluator.py"
            evaluator.write_text(
                """
import json
from pathlib import Path
import sys
import time

assert sys.argv[1] == "--persistent-server"
print("GRAVITYX_EVALUATOR_READY", flush=True)
for line in sys.stdin:
    request = json.loads(line)
    if request.get("stop"):
        break
    case_dir = Path(request["case_dir"])
    solution_dir = Path(request["solution_dir"])
    labels = [
        row.split()[1]
        for row in (case_dir / "case.con").read_text(encoding="utf-8").splitlines()
        if row.startswith("CONTINGENCY ")
    ]
    objectives = {"BASECASE": 10.0}
    infeasibilities = {"BASECASE": False}
    for label in labels:
        objectives[label] = float(ord(label[-1]) - ord("A") + 1)
        infeasibilities[label] = False
    for label, objective in objectives.items():
        (solution_dir / f"eval_detail_{label}.json").write_text(
            json.dumps({
                "obj": {"val": objective},
                "infeas": {"val": False},
            }),
            encoding="utf-8",
        )
    (solution_dir / "eval_summary.json").write_text(
        json.dumps({
            "solutions_exist": True,
            "num_ctg": len(labels),
            "obj": 10.0 + sum(objectives[label] for label in labels) / len(labels),
            "infeas": 0.0,
        }),
        encoding="utf-8",
    )
    print(
        "GRAVITYX_EVALUATION_RESULT "
        + json.dumps({
            "request_id": request["request_id"],
            "success": True,
            "wall_seconds": 0.01,
            "static_case_cache": {"raw_hits": request["request_id"]},
        }),
        flush=True,
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
                4,
                1,
                time.perf_counter() + 30.0,
                1,
                1,
                completion_order_groups=True,
                persistent_evaluator_processes=True,
            )
            for label in labels:
                (output_dir / f"solution_{label}.txt").write_text(
                    f"{label}\n", encoding="utf-8"
                )
                manager.mark_completed(label)
            summary, metadata = manager.finish()

            self.assertEqual(summary["num_ctg"], 4)
            self.assertEqual(summary["infeas"], 0.0)
            self.assertTrue(metadata["persistent_evaluator_processes"])
            self.assertEqual(metadata["persistent_evaluator_process_count"], 1)
            self.assertEqual(
                metadata["persistent_evaluator_task_counts"], {"0": 4}
            )

    def test_code2_budget_uses_contingency_count(self):
        self.assertEqual(code2_time_limit(105, 2.0), 210.0)
        self.assertEqual(code2_time_limit(0, 2.0), 0.0)

    def test_process_timeout_is_bounded_by_stage_deadline(self):
        self.assertEqual(effective_process_timeout(300.0, 120.0, now=100.0), 20.0)
        self.assertEqual(effective_process_timeout(5.0, 120.0, now=100.0), 5.0)

    def test_cold_run_free_space_guard_fails_before_writing(self):
        usage = SimpleNamespace(total=100 * 1024**3, used=99 * 1024**3, free=1 * 1024**3)
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / "future" / "run"
            with mock.patch("run_experiment.shutil.disk_usage", return_value=usage):
                with self.assertRaisesRegex(RuntimeError, "insufficient free disk"):
                    require_minimum_free_space(destination, 2.0)
                self.assertEqual(
                    require_minimum_free_space(destination, 1.0), 1.0
                )

    def test_wsl_worker_input_staging_verifies_both_hashes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            case = root / "case.json"
            base = root / "base.json"
            case.write_text('{"case":1}\n', encoding="utf-8")
            base.write_text('{"base":1}\n', encoding="utf-8")
            case_hash = __import__("hashlib").sha256(case.read_bytes()).hexdigest()
            base_hash = __import__("hashlib").sha256(base.read_bytes()).hexdigest()
            registered = []

            def fake_run(command, **kwargs):
                del kwargs
                if "sha256sum" in command:
                    return subprocess.CompletedProcess(
                        command,
                        0,
                        f"{case_hash}  {command[-2]}\n{base_hash}  {command[-1]}\n",
                    )
                return subprocess.CompletedProcess(command, 0, "")

            with mock.patch(
                "run_experiment.subprocess.run", side_effect=fake_run
            ), mock.patch(
                "run_experiment.atexit.register",
                side_effect=lambda callback: registered.append(callback),
            ):
                staged_case, staged_base, metadata = stage_wsl_worker_inputs(
                    "Test-Distro",
                    case,
                    base,
                    time.perf_counter() + 10.0,
                )
                self.assertTrue(staged_case.startswith("/dev/shm/gravityx_go2_"))
                self.assertTrue(staged_base.startswith("/dev/shm/gravityx_go2_"))
                self.assertEqual(
                    metadata["sha256"], {"case": case_hash, "base": base_hash}
                )
                self.assertEqual(len(registered), 1)
                registered[0]()

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

    def test_json_checkpoint_uses_short_temporary_name(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            filename = "official_evaluation_certificate.json"
            padding = 250 - len(str(root)) - len(filename) - 2
            self.assertGreater(padding, 0)
            destination = root / ("x" * padding) / filename
            self.assertEqual(len(str(destination)), 250)

            write_json(destination, {"success": True})

            self.assertEqual(json.loads(destination.read_text()), {"success": True})
            self.assertEqual(list(destination.parent.glob("*.tmp")), [])


if __name__ == "__main__":
    unittest.main()
