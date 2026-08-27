#!/usr/bin/env python3
"""Audit retained GO2 evaluator details without rerunning a scenario.

The official score is a base-case objective plus the arithmetic mean of all
contingency objectives.  Large campaigns retain the detail files in evaluator
shards, with one identical BASECASE detail per shard.  This module verifies
that evidence boundary, reconstructs the score, and reports every evaluator
``total_*`` component separately for the base case and contingency mean.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable


OBJECTIVE_TERMS = {
    "total_load_benefit": 1.0,
    "total_bus_cost": -1.0,
    "total_gen_cost": -1.0,
    "total_line_cost": -1.0,
    "total_xfmr_cost": -1.0,
}

CONTROL_METRICS = (
    "load_pow_real_delta_to_prior",
    "gen_pow_real_delta_to_prior",
    "gen_pow_imag_delta_to_prior",
    "xfmr_tau_delta_to_prior",
    "xfmr_phi_delta_to_prior",
    "swsh_b_delta_to_prior",
    "line_switch_up_actual",
    "line_switch_down_actual",
    "xfmr_switch_up_actual",
    "xfmr_switch_down_actual",
)

PENALTY_METRICS = (
    "max_bus_pow_real_over",
    "max_bus_pow_real_under",
    "max_bus_pow_imag_over",
    "max_bus_pow_imag_under",
    "sum_bus_pow_real_over",
    "sum_bus_pow_real_under",
    "sum_bus_pow_imag_over",
    "sum_bus_pow_imag_under",
    "max_line_viol",
    "max_xfmr_viol",
)


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as stream:
        value = json.load(stream)
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object in {path}")
    return value


def _scalar_values(detail: dict[str, Any]) -> dict[str, float]:
    values: dict[str, float] = {}
    for key, entry in detail.items():
        if not isinstance(entry, dict) or "val" not in entry:
            continue
        value = entry["val"]
        if isinstance(value, bool):
            continue
        if isinstance(value, (int, float)) and math.isfinite(float(value)):
            values[key] = float(value)
    return values


def _detail_label(path: Path) -> str:
    prefix = "eval_detail_"
    if not path.stem.startswith(prefix):
        raise ValueError(f"unexpected evaluator-detail filename: {path}")
    return path.stem[len(prefix) :]


def _close(left: float, right: float, tolerance: float = 1e-8) -> bool:
    return math.isclose(left, right, rel_tol=1e-12, abs_tol=tolerance)


def _mean(records: Iterable[dict[str, float]], key: str) -> float:
    values = [record.get(key, 0.0) for record in records]
    return math.fsum(values) / len(values) if values else 0.0


def _max(records: Iterable[dict[str, float]], key: str) -> float:
    values = [record.get(key, 0.0) for record in records]
    return max(values) if values else 0.0


def _ordered_components(
    base: dict[str, float], contingencies: list[dict[str, float]]
) -> list[str]:
    available = {
        key
        for record in [base, *contingencies]
        for key in record
        if key.startswith("total_")
    }
    nonzero = {
        key
        for record in [base, *contingencies]
        for key, value in record.items()
        if key.startswith("total_") and value != 0.0
    }
    preferred = [
        "total_load_benefit",
        "total_bus_cost",
        "total_bus_real_cost",
        "total_bus_imag_cost",
        "total_gen_cost",
        "total_gen_energy_cost",
        "total_gen_on_cost",
        "total_gen_su_cost",
        "total_gen_sd_cost",
        "total_line_cost",
        "total_line_limit_cost",
        "total_line_switch_cost",
        "total_xfmr_cost",
        "total_xfmr_limit_cost",
        "total_xfmr_switch_cost",
    ]
    return [key for key in preferred if key in available] + sorted(
        nonzero - set(preferred)
    )


def _source_state_summary(
    internal_base: dict[str, Any], model: dict[str, Any] | None
) -> dict[str, Any]:
    commitment = [int(value) for value in internal_base.get("commitment", [])]
    state = internal_base.get("selected_state", {})
    if not isinstance(state, dict):
        state = {}
    result: dict[str, Any] = {
        "base_method": internal_base.get("base_method"),
        "base_optimization_performed": internal_base.get(
            "base_optimization_performed"
        ),
        "generator_count": len(commitment),
        "committed_unit_count": sum(commitment),
        "commitment_binary": all(value in (0, 1) for value in commitment),
        "startup_count": sum(
            1 for value in state.get("startup", []) if float(value) > 0.5
        ),
        "shutdown_count": sum(
            1 for value in state.get("shutdown", []) if float(value) > 0.5
        ),
    }
    if model is None:
        return result

    generators = sorted(model.get("gen", {}).values(), key=lambda row: row["index"])
    loads = sorted(model.get("load", {}).values(), key=lambda row: row["index"])
    source_commitment = [int(row["status_prev"]) for row in generators]
    result.update(
        {
            "source_committed_unit_count": sum(source_commitment),
            "commitment_changes_from_source": sum(
                left != right
                for left, right in zip(commitment, source_commitment, strict=True)
            ),
            "source_all_generators_online": bool(source_commitment)
            and all(source_commitment),
        }
    )
    base_mva = float(model.get("baseMVA", 100.0))
    pg = [float(value) for value in state.get("pg", [])]
    qg = [float(value) for value in state.get("qg", [])]
    demand = [float(value) for value in state.get("demand_factor", [])]
    if len(pg) == len(generators):
        result["generator_dispatch_mw"] = math.fsum(pg) * base_mva
        result["source_generator_dispatch_mw"] = (
            math.fsum(float(row["pg_prev"]) for row in generators) * base_mva
        )
        result["generators_with_real_dispatch_change"] = sum(
            not _close(value, float(row["pg_prev"]), 1e-10)
            for value, row in zip(pg, generators, strict=True)
        )
    if len(qg) == len(generators):
        result["generator_reactive_dispatch_mvar"] = math.fsum(qg) * base_mva
        result["source_generator_reactive_dispatch_mvar"] = (
            math.fsum(float(row["qg_prev"]) for row in generators) * base_mva
        )
        result["generators_with_reactive_dispatch_change"] = sum(
            not _close(value, float(row["qg_prev"]), 1e-10)
            for value, row in zip(qg, generators, strict=True)
        )
    if len(demand) == len(loads):
        load_mw = [
            float(row["pd_nominal"]) * value * base_mva
            for value, row in zip(demand, loads, strict=True)
        ]
        prior_load_mw = [float(row["pd_prev"]) * base_mva for row in loads]
        prior_factor = [
            (
                float(row["pd_prev"]) / float(row["pd_nominal"])
                if abs(float(row["pd_nominal"])) > 1e-14
                else float(row["tmin"])
            )
            for row in loads
        ]
        result.update(
            {
                "load_served_mw": math.fsum(load_mw),
                "source_load_mw": math.fsum(prior_load_mw),
                "loads_with_quantity_change": sum(
                    not _close(value, prior, 1e-10)
                    for value, prior in zip(demand, prior_factor, strict=True)
                ),
                "loads_at_tmin": sum(
                    _close(value, float(row["tmin"]), 1e-10)
                    for value, row in zip(demand, loads, strict=True)
                ),
                "loads_at_tmax": sum(
                    _close(value, float(row["tmax"]), 1e-10)
                    for value, row in zip(demand, loads, strict=True)
                ),
            }
        )
    return result


def analyze_run(run_dir: Path, model_json: Path | None = None) -> dict[str, Any]:
    """Return a self-validating objective audit for one retained run."""

    run_dir = run_dir.resolve()
    summary_path = run_dir / "eval_summary.json"
    internal_base_path = run_dir / "internal" / "base.json"
    if not summary_path.is_file() or not internal_base_path.is_file():
        raise FileNotFoundError(f"retained run is incomplete: {run_dir}")

    summary = _read_json(summary_path)
    internal_base = _read_json(internal_base_path)
    paths = sorted(
        (run_dir / "internal" / "streaming_serial_evaluation_shards").glob(
            "shard_*/solutions/eval_detail_*.json"
        )
    )
    if not paths:
        raise FileNotFoundError(f"no retained evaluator details under {run_dir}")

    by_label: dict[str, list[tuple[Path, dict[str, Any]]]] = {}
    for path in paths:
        by_label.setdefault(_detail_label(path), []).append((path, _read_json(path)))
    if "BASECASE" not in by_label:
        raise ValueError("evaluator details do not contain BASECASE")

    base_records = by_label.pop("BASECASE")
    base_hashes = {
        hashlib.sha256(path.read_bytes()).hexdigest() for path, _ in base_records
    }
    if len(base_hashes) != 1:
        raise ValueError("duplicated BASECASE evaluator details are not byte-identical")
    duplicate_contingencies = {
        label: len(records) for label, records in by_label.items() if len(records) != 1
    }
    if duplicate_contingencies:
        raise ValueError(
            f"contingency evaluator-detail multiplicity error: {duplicate_contingencies}"
        )

    base = _scalar_values(base_records[0][1])
    contingency_labels = sorted(by_label)
    contingencies = [_scalar_values(by_label[label][0][1]) for label in contingency_labels]
    expected_count = int(summary["num_ctg"])
    if len(contingencies) != expected_count:
        raise ValueError(
            f"retained contingency detail count {len(contingencies)} != {expected_count}"
        )
    if any(by_label[label][0][1].get("infeas", {}).get("val") is not False for label in contingency_labels):
        raise ValueError("at least one retained contingency detail is infeasible")
    if base_records[0][1].get("infeas", {}).get("val") is not False:
        raise ValueError("retained base detail is infeasible")

    reconstructed = base["obj"] + _mean(contingencies, "obj")
    reported = float(summary["obj"])
    if not _close(reconstructed, reported, 1e-6):
        raise ValueError(
            f"objective reconstruction {reconstructed} != reported {reported}"
        )

    for label, record in [("BASECASE", base), *zip(contingency_labels, contingencies)]:
        reconstructed_case = math.fsum(
            sign * record.get(key, 0.0) for key, sign in OBJECTIVE_TERMS.items()
        )
        if not _close(reconstructed_case, record["obj"], 1e-5):
            raise ValueError(
                f"{label} objective components {reconstructed_case} != {record['obj']}"
            )

    components: dict[str, dict[str, float]] = {}
    for key in _ordered_components(base, contingencies):
        base_value = base.get(key, 0.0)
        contingency_mean = _mean(contingencies, key)
        components[key] = {
            "base": base_value,
            "contingency_mean": contingency_mean,
            "aggregate": base_value + contingency_mean,
        }

    control_metrics = {
        key: {
            "base": base.get(key, 0.0),
            "contingency_mean": _mean(contingencies, key),
            "contingency_max": _max(contingencies, key),
        }
        for key in CONTROL_METRICS
    }
    penalty_metrics = {
        key: {
            "base": base.get(key, 0.0),
            "contingency_mean": _mean(contingencies, key),
            "contingency_max": _max(contingencies, key),
        }
        for key in PENALTY_METRICS
    }

    model = _read_json(model_json.resolve()) if model_json is not None else None
    return {
        "run_dir": str(run_dir),
        "model_json": str(model_json.resolve()) if model_json is not None else None,
        "reported_objective": reported,
        "reconstructed_objective": reconstructed,
        "base_objective": base["obj"],
        "contingency_mean_objective": _mean(contingencies, "obj"),
        "contingency_count": len(contingencies),
        "base_detail_copy_count": len(base_records),
        "base_detail_sha256": next(iter(base_hashes)),
        "components": components,
        "control_metrics": control_metrics,
        "penalty_metrics": penalty_metrics,
        "state": _source_state_summary(internal_base, model),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--model-json", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = analyze_run(args.run_dir, args.model_json)
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        if "onedrive" in str(args.output).casefold():
            raise ValueError(f"refusing OneDrive output path: {args.output}")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
