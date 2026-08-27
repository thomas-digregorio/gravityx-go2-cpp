import json
import math
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from analyze_retained_objective import analyze_run  # noqa: E402


def _metric(value):
    return {"val": value}


def _detail(*, benefit, bus, gen, line, xfmr=0.0, infeas=False):
    objective = benefit - bus - gen - line - xfmr
    return {
        "infeas": _metric(infeas),
        "obj": _metric(objective),
        "total_load_benefit": _metric(benefit),
        "total_bus_cost": _metric(bus),
        "total_bus_real_cost": _metric(0.75 * bus),
        "total_bus_imag_cost": _metric(0.25 * bus),
        "total_gen_cost": _metric(gen),
        "total_gen_energy_cost": _metric(gen - 1.0),
        "total_gen_on_cost": _metric(1.0),
        "total_gen_su_cost": _metric(0.0),
        "total_gen_sd_cost": _metric(0.0),
        "total_line_cost": _metric(line),
        "total_line_limit_cost": _metric(line),
        "total_line_switch_cost": _metric(0.0),
        "total_xfmr_cost": _metric(xfmr),
        "total_xfmr_limit_cost": _metric(xfmr),
        "total_xfmr_switch_cost": _metric(0.0),
    }


class RetainedObjectiveAnalysisTests(unittest.TestCase):
    def test_reconstructs_base_plus_contingency_mean_and_components(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            details = root / "internal" / "streaming_serial_evaluation_shards"
            base = _detail(benefit=100.0, bus=10.0, gen=20.0, line=1.0)
            contingencies = [
                _detail(benefit=90.0, bus=15.0, gen=18.0, line=2.0),
                _detail(benefit=80.0, bus=25.0, gen=17.0, line=3.0),
            ]
            for index, detail in enumerate(contingencies):
                shard = details / f"shard_{index:03d}" / "solutions"
                shard.mkdir(parents=True)
                (shard / "eval_detail_BASECASE.json").write_text(
                    json.dumps(base), encoding="utf-8"
                )
                (shard / f"eval_detail_CTG_{index:06d}.json").write_text(
                    json.dumps(detail), encoding="utf-8"
                )
            expected = base["obj"]["val"] + math.fsum(
                detail["obj"]["val"] for detail in contingencies
            ) / len(contingencies)
            (root / "eval_summary.json").write_text(
                json.dumps({"num_ctg": 2, "obj": expected}), encoding="utf-8"
            )
            (root / "internal" / "base.json").write_text(
                json.dumps(
                    {
                        "base_method": "fixture",
                        "base_optimization_performed": False,
                        "commitment": [1, 0],
                        "selected_state": {"startup": [0, 0], "shutdown": [0, 0]},
                    }
                ),
                encoding="utf-8",
            )

            report = analyze_run(root)

            self.assertAlmostEqual(report["reconstructed_objective"], expected)
            self.assertAlmostEqual(
                report["components"]["total_bus_cost"]["contingency_mean"],
                20.0,
            )
            self.assertEqual(report["contingency_count"], 2)
            self.assertEqual(report["base_detail_copy_count"], 2)
            self.assertEqual(report["state"]["committed_unit_count"], 1)
            self.assertAlmostEqual(
                report["components"]["total_gen_on_cost"]["aggregate"],
                2.0,
            )
            self.assertAlmostEqual(
                report["components"]["total_gen_su_cost"]["aggregate"],
                0.0,
            )
            self.assertAlmostEqual(
                report["components"]["total_gen_sd_cost"]["aggregate"],
                0.0,
            )

    def test_rejects_incomplete_contingency_evidence(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            shard = (
                root
                / "internal"
                / "streaming_serial_evaluation_shards"
                / "shard_000"
                / "solutions"
            )
            shard.mkdir(parents=True)
            base = _detail(benefit=1.0, bus=0.0, gen=0.0, line=0.0)
            (shard / "eval_detail_BASECASE.json").write_text(
                json.dumps(base), encoding="utf-8"
            )
            (root / "eval_summary.json").write_text(
                json.dumps({"num_ctg": 1, "obj": 1.0}), encoding="utf-8"
            )
            (root / "internal" / "base.json").write_text(
                json.dumps({"commitment": [], "selected_state": {}}), encoding="utf-8"
            )

            with self.assertRaisesRegex(ValueError, "detail count"):
                analyze_run(root)


if __name__ == "__main__":
    unittest.main()
