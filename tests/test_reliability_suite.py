import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from run_reliability_suite import audit_success, normalized_case, arguments
from run_experiment import load_fast_screen_heavy_profile


class ReliabilityAuditTests(unittest.TestCase):
    def setUp(self):
        self.status = {"success": True, "stage": "complete", "completed_contingency_count": 2}
        self.summary = {
            "git_revision": "frozen", "base_exact_executable_sha256": "binary",
            "fast_screen_executable_sha256": "binary", "total_wall_seconds": 10,
            "contingency_count": 2, "max_independent_contingency_residual": 1e-7,
            "base": {"base_validation": {"max_residual": 2e-7}},
            "official_evaluation": {"obj": -12, "infeas": 0.0},
            "official_evaluation_certificate": {"complete_label_set": True,
                "expected_detail_count": 3, "observed_unique_detail_count": 3,
                "reported_infeasibility": 0.0},
        }

    def check(self, summary=None, status=None):
        return audit_success(status or self.status, summary or self.summary, 2,
                             "frozen", "binary", 300)

    def test_complete_negative_objective_is_not_infeasibility(self):
        self.assertEqual(self.check(), [])

    def test_missing_or_nonfinite_fields_fail_closed(self):
        for key in ("total_wall_seconds", "max_independent_contingency_residual"):
            for value in (None, float("nan"), float("inf")):
                candidate = copy.deepcopy(self.summary)
                candidate[key] = value
                self.assertTrue(self.check(candidate))

    def test_complete_counts_do_not_replace_physical_verification(self):
        self.summary["base"]["base_validation"]["max_residual"] = 0.0557
        self.assertIn("independent residual", self.check())

    def test_missing_contingency_or_changed_executable_rejected(self):
        self.status["completed_contingency_count"] = 1
        self.assertIn("incomplete contingency count", self.check())
        self.status["completed_contingency_count"] = 2
        self.summary["fast_screen_executable_sha256"] = "changed"
        self.assertTrue(self.check())

    def test_deadline_and_certificate_required(self):
        self.summary["total_wall_seconds"] = 300.001
        self.assertIn("end-to-end deadline", self.check())
        self.summary["total_wall_seconds"] = 299
        self.summary["official_evaluation_certificate"]["complete_label_set"] = False
        self.assertIn("official certificate", self.check())

    def test_case_paths_preserve_exact_scenario_identity(self):
        self.assertEqual(normalized_case(Path("repo"), "00617", "005"),
                         Path("repo/.data/final_617_005/model.json"))
        self.assertEqual(normalized_case(Path("repo"), "19402", "095"),
                         Path("repo/.data/final_19402/scenario_095/model.json"))

    def test_unfinished_priority_is_not_a_measured_duration(self):
        profile = {"schema_version": 3, "case_sha256": "case", "heavy_threshold_seconds": 10,
                   "contingencies": [{"label": "a", "measured_solver_wall_seconds": 20}],
                   "unfinished_priority_labels": ["b"]}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "profile.json"
            path.write_text(json.dumps(profile))
            labels, ordering, meta = load_fast_screen_heavy_profile(path, "case", {"a", "b"})
            self.assertEqual(labels, {"a", "b"})
            self.assertGreater(ordering["b"], ordering["a"])
            self.assertEqual(meta["measured_solver_seconds_sum"], 20)
            self.assertEqual(meta["unfinished_priority_labels"], ["b"])
            self.assertFalse(meta["uses_prior_solution_state"])
            for bad in (["unknown"], ["b", "b"], [1]):
                profile["unfinished_priority_labels"] = bad
                path.write_text(json.dumps(profile))
                with self.assertRaises(ValueError):
                    load_fast_screen_heavy_profile(path, "case", {"a", "b"})
            profile["unfinished_priority_labels"] = ["b"]
            profile["schema_version"] = 2
            path.write_text(json.dumps(profile))
            with self.assertRaises(ValueError):
                load_fast_screen_heavy_profile(path, "case", {"a", "b"})

    def test_suite_profiles_only_change_queue_configuration(self):
        config = {"python": "python", "data_repository": "data", "source_root": "sources",
                  "vendor_evaluator": "evaluator", "total_time_limit": 300,
                  "minimum_free_space_gib": 30}
        before = arguments(config, "19402", "010", Path("run"))
        config["screen_profiles"] = {"19402/010": "config/profile.json"}
        after = arguments(config, "19402", "010", Path("run"))
        self.assertEqual(after[:len(before)], before)
        self.assertEqual(after[len(before)], "--fast-screen-heavy-profile")
        self.assertEqual(after[-2:], ["--fast-screen-heavy-workers", "4"])


if __name__ == "__main__":
    unittest.main()
