from __future__ import annotations

import copy
from pathlib import Path
import sys
import threading
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from run_experiment import (AffinityScreenWorkQueue, fast_screen_affinity_groups,
                            retain_exact_parallel_group)
from run_reliability_suite import arguments


class ExactParallelReuseTests(unittest.TestCase):
    def fixture(self):
        branch = dict(index=1, source_id=["branch", 10, 20, "1"], f_bus=10, t_bus=20,
                      br_r=0.01, br_x=0.1, b_fr=0.0, b_to=0.0, g_fr=0.0, g_to=0.0,
                      rate_a=2.0, rate_b=2.0, rate_c=2.0, tap=1.0, shift=0.0,
                      angmin=-1.0, angmax=1.0, br_status=1, present=True,
                      transformer=False, control_mode=0)
        peer = copy.deepcopy(branch)
        peer.update(index=2, source_id=["branch", 10, 20, "2"])
        different = copy.deepcopy(branch)
        different.update(index=3, t_bus=30, source_id=["branch", 10, 30, "1"])
        case = {"gen": {}, "branch": {"1": branch, "2": peer, "3": different}}
        state = dict(pg=[], qg=[], pf=[0.5] * 3, pt=[-0.5] * 3,
                     qf=[0.1] * 3, qt=[-0.1] * 3)
        records = [dict(label=f"B{i}", type="branch", idx=i, schedule_rank=i) for i in (1, 2, 3)]
        return case, state, records

    def groups(self, case, state, records):
        return fast_screen_affinity_groups(case, state, records, exact_parallel_outage_reuse=True)

    def test_exact_peer_identity_and_complete_label_coverage(self):
        case, state, records = self.fixture()
        original = copy.deepcopy((case, state, records))
        groups = self.groups(case, state, records)
        peers = [g for g in groups if retain_exact_parallel_group(g)]
        self.assertEqual(len(peers), 1)
        self.assertEqual([r["label"] for r in peers[0]], ["B1", "B2"])
        self.assertEqual(peers[0][0]["fast_screen_affinity_value"], 1)
        self.assertCountEqual([r["label"] for g in groups for r in g], ["B1", "B2", "B3"])
        self.assertEqual((case, state, records), original)

    def test_disabled_option_preserves_old_from_bus_grouping(self):
        case, state, records = self.fixture()
        groups = fast_screen_affinity_groups(case, state, records)
        self.assertEqual(len(groups), 1)
        self.assertFalse(retain_exact_parallel_group(groups[0]))
        self.assertTrue(all(r["fast_screen_affinity_type"] == "branch_from_bus" for r in groups[0]))

    def test_each_changed_source_attribute_breaks_matching(self):
        mutations = dict(br_r=0.02, br_x=0.2, rate_a=3.0, rate_b=3.0, rate_c=3.0,
                         tap=1.01, shift=0.01, br_status=0, present=False,
                         f_bus=11, t_bus=21, control_mode=1, b_fr=0.01,
                         b_to=0.01, g_fr=0.01, g_to=0.01, angmin=-0.5,
                         angmax=0.5, transformer=True, extra_source_attribute="different")
        for field, value in mutations.items():
            with self.subTest(field=field):
                case, state, records = self.fixture()
                case["branch"]["2"][field] = value
                self.assertFalse(any(retain_exact_parallel_group(g) for g in self.groups(case, state, records)))

    def test_nonfinite_fingerprint_is_not_silently_normalized(self):
        case, state, records = self.fixture()
        case["branch"]["2"]["br_x"] = float("nan")
        with self.assertRaises(ValueError):
            self.groups(case, state, records)

    def test_retention_requires_one_exact_group(self):
        case, state, records = self.fixture()
        groups = self.groups(case, state, records)
        exact = next(g for g in groups if retain_exact_parallel_group(g))
        ordinary = next(g for g in groups if not retain_exact_parallel_group(g))
        self.assertFalse(retain_exact_parallel_group([]))
        self.assertFalse(retain_exact_parallel_group(exact + ordinary))
        changed = copy.deepcopy(exact)
        changed[1]["fast_screen_affinity_value"] = 99
        self.assertFalse(retain_exact_parallel_group(changed))

    def test_successful_exact_group_stays_serial_without_lost_labels(self):
        case, state, records = self.fixture()
        queue = AffinityScreenWorkQueue(self.groups(case, state, records), worker_count=2)
        abort = threading.Event()
        source, group = queue.get(abort, worker_id=0)
        self.assertTrue(retain_exact_parallel_group(group))
        seen = [item["label"] for item in group]
        queue.task_done(source)
        source, ordinary = queue.get(abort, worker_id=1)
        seen.extend(item["label"] for item in ordinary)
        queue.task_done(source)
        self.assertCountEqual(seen, ["B1", "B2", "B3"])
        self.assertEqual(len(seen), len(set(seen)))
        self.assertEqual(queue.remaining_group_count, 0)
        self.assertIsNone(queue.get(abort, worker_id=0))

    def test_failed_group_can_still_split_and_complete_every_label(self):
        case, state, records = self.fixture()
        queue = AffinityScreenWorkQueue(self.groups(case, state, records), worker_count=2)
        abort = threading.Event()
        source, group = queue.get(abort, worker_id=0)
        seen = [group[0]["label"]]
        self.assertEqual(queue.requeue_remaining_as_singletons(group, 1, 0, source), 1)
        queue.task_done(source)
        source, ordinary = queue.get(abort, worker_id=0)
        seen.extend(item["label"] for item in ordinary)
        queue.task_done(source)
        source, peer = queue.get(abort, worker_id=1)
        seen.extend(item["label"] for item in peer)
        queue.task_done(source)
        self.assertCountEqual(seen, ["B1", "B2", "B3"])
        self.assertEqual(len(seen), len(set(seen)))
        self.assertEqual(queue.remaining_group_count, 0)
        self.assertIsNone(queue.get(abort, worker_id=0))

    def test_registered_flag_requires_economic_mode(self):
        config = dict(data_repository="local-data", source_root="local-source", python="python",
                      total_time_limit=300, minimum_free_space_gib=30, vendor_evaluator="vendor.py",
                      exact_parallel_outage_reuse=True)
        with self.assertRaises(ValueError):
            arguments(config, "19402", "095", Path("local-output"))
        config["cached_economic_contingency_polish"] = True
        command = arguments(config, "19402", "095", Path("local-output"))
        self.assertIn("--exact-parallel-outage-reuse", command)
        self.assertIn("--fast-screen-affinity-schedule", command)
        config["exact_parallel_outage_reuse"] = False
        self.assertNotIn("--exact-parallel-outage-reuse", arguments(config, "19402", "095", Path("local-output")))


if __name__ == "__main__":
    unittest.main()
