import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from run_reliability_suite import audit_success, normalized_case, arguments
from run_experiment import (load_fast_screen_heavy_profile, start_worker_if_task,
                            additional_corrective_workers, cpp_command)


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

    def test_missing_summary_is_not_reported_as_observed_hash_mismatch(self):
        failures = audit_success({"success": False, "stage": "code2",
                                  "completed_contingency_count": 1}, {}, 2,
                                 "frozen", "binary", 300)
        self.assertIn("final summary and verification evidence missing", failures)
        self.assertIn("incomplete contingency count", failures)
        self.assertFalse(any("mismatch" in message for message in failures))

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

    def test_cached_polish_is_explicit_without_external_solution_start(self):
        config = {"python": "python", "data_repository": "data", "source_root": "sources",
                  "vendor_evaluator": "evaluator", "total_time_limit": 300,
                  "minimum_free_space_gib": 30}
        baseline = arguments(config, "19402", "095", Path("run"))
        config["cached_economic_contingency_polish"] = True
        candidate = arguments(config, "19402", "095", Path("run"))
        self.assertEqual(candidate, baseline + ["--cached-economic-contingency-polish"])
        self.assertNotIn("--economic-contingency-polish", candidate)
        self.assertIn("--validated-source-base", candidate)
        self.assertNotIn("--skip-evaluation", candidate)

    def test_epigraph_changes_only_an_explicit_base_flag(self):
        config = {"python": "python", "data_repository": "data", "source_root": "sources",
                  "vendor_evaluator": "evaluator", "total_time_limit": 300,
                  "minimum_free_space_gib": 30}
        before = arguments(config, "19402", "095", Path("run"))
        config["base_pwl_epigraph"] = True
        after = arguments(config, "19402", "095", Path("run"))
        self.assertEqual(after, before + ["--base-pwl-epigraph"])
        self.assertIn("--validated-source-base", after)

    def test_early_evaluator_pool_changes_only_explicit_timed_preparation(self):
        config = {"python": "python", "data_repository": "data", "source_root": "sources",
                  "vendor_evaluator": "evaluator", "total_time_limit": 300,
                  "minimum_free_space_gib": 30}
        before = arguments(config, "19402", "095", Path("run"))
        config["prepare_evaluator_pool_early"] = True
        after = arguments(config, "19402", "095", Path("run"))
        self.assertEqual(after, before + ["--prepare-evaluator-pool-early"])
        self.assertNotIn("--base-json", after)
        self.assertIn("--validated-source-base", after)
        self.assertEqual(after[after.index("--total-time-limit") + 1], "300")

    def test_exact_hessian_requires_epigraph_and_preserves_cold_limits(self):
        config = {"python": "python", "data_repository": "data", "source_root": "sources",
                  "vendor_evaluator": "evaluator", "total_time_limit": 300,
                  "minimum_free_space_gib": 30, "base_pwl_epigraph": True}
        before = arguments(config, "19402", "095", Path("run"))
        config["base_exact_hessian"] = True
        after = arguments(config, "19402", "095", Path("run"))
        self.assertEqual(after, before + ["--base-exact-hessian"])
        self.assertEqual(after[after.index("--total-time-limit") + 1], "300")
        self.assertNotIn("--base-json", after)
        self.assertNotIn("--skip-evaluation", after)
        self.assertNotIn("--common-corrective-reference-seconds", after)
        config["base_pwl_epigraph"] = False
        with self.assertRaises(ValueError):
            arguments(config, "19402", "095", Path("run"))

    def test_common_corrective_reference_is_explicit_and_cold(self):
        config = {"python": "python", "data_repository": "data", "source_root": "sources",
                  "vendor_evaluator": "evaluator", "total_time_limit": 300,
                  "minimum_free_space_gib": 30}
        before = arguments(config, "19402", "095", Path("run"))
        config["common_corrective_reference_seconds"] = 30
        after = arguments(config, "19402", "095", Path("run"))
        self.assertEqual(after, before + ["--common-corrective-reference-seconds", "30.0"])
        self.assertIn("--validated-source-base", after)
        self.assertEqual(after[after.index("--total-time-limit") + 1], "300")
        self.assertNotIn("--skip-evaluation", after)
        self.assertNotIn("--base-json", after)
        for invalid in (-1, float("nan"), float("inf")):
            config["common_corrective_reference_seconds"] = invalid
            with self.assertRaises(ValueError):
                arguments(config, "19402", "095", Path("run"))

    def test_predictor_handoff_budget_is_explicit_hash_bound_and_not_a_measurement(self):
        profile = {"schema_version": 4, "case_sha256": "case", "heavy_threshold_seconds": 10,
                   "contingencies": [{"label": "a", "measured_solver_wall_seconds": 20}],
                   "unfinished_priority_labels": ["b"], "screen_handoff_seconds": {"b": 15}}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "profile.json"
            path.write_text(json.dumps(profile))
            labels, _, meta = load_fast_screen_heavy_profile(path, "case", {"a", "b"})
            self.assertEqual(labels, {"a", "b"})
            self.assertEqual(meta["screen_handoff_seconds"], {"b": 15})
            self.assertEqual(meta["measured_solver_seconds_sum"], 20)
            self.assertFalse(meta["uses_prior_solution_state"])
            for bad in ({"unknown": 15}, {"b": 0}, {"b": -1}, {"b": 300},
                        {"b": float("nan")}, {"b": float("inf")}, {"b": True},
                        {"b": "15"}, [15]):
                profile["screen_handoff_seconds"] = bad
                path.write_text(json.dumps(profile))
                with self.assertRaises(ValueError):
                    load_fast_screen_heavy_profile(path, "case", {"a", "b"})
            profile["screen_handoff_seconds"] = {"b": 15}
            profile["schema_version"] = 3
            path.write_text(json.dumps(profile))
            with self.assertRaises(ValueError):
                load_fast_screen_heavy_profile(path, "case", {"a", "b"})

    def test_empty_corrective_queue_never_starts_native_solver(self):
        def forbidden():
            self.fail("Empty work queue started a solver")
        self.assertEqual(start_worker_if_task(lambda: None, forbidden), (None, None))
        self.assertEqual(additional_corrective_workers(4, 8, 0), 0)

    def test_repair_logging_does_not_change_process_deadline(self):
        normal = cpp_command(Path("C:/local/solver"), "Ubuntu-24.04", ["worker"], 123)
        logged = cpp_command(Path("C:/local/solver"), "Ubuntu-24.04", ["worker"], 123,
                             repair_logging=True)
        self.assertEqual([s for s in logged if not s.startswith("GRAVITYX_")], normal)
        self.assertIn("123.000s", logged)
        self.assertIn("GRAVITYX_HIGHS_LOG=1", logged)
        self.assertIn("GRAVITYX_REPAIR_LOG=1", logged)

    def test_corrective_process_starts_only_after_claiming_real_task(self):
        order = []
        def claim():
            order.append("claim")
            return {"label": "a"}
        def launch():
            order.append("launch")
            return "native"
        self.assertEqual(start_worker_if_task(claim, launch), ({"label": "a"}, "native"))
        self.assertEqual(order, ["claim", "launch"])
        self.assertEqual(additional_corrective_workers(4, 8, 2), 2)
        self.assertEqual(additional_corrective_workers(4, 8, 100), 4)


if __name__ == "__main__":
    unittest.main()
