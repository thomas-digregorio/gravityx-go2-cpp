"""Run a frozen cold scenario suite, retaining evidence and never retrying.

This is orchestration, not a solver.  All physical acceptance is performed by
the unchanged independent native checker and pinned official evaluator.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

from archive_run_provenance import archive_run, reject_onedrive
from run_experiment import read_contingency_blocks, write_json

REPO = Path(__file__).resolve().parents[1]


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def normalized_case(root: Path, family: str, scenario: str) -> Path:
    if family in ("17700", "19402"):
        return root / ".data" / f"final_{family}" / f"scenario_{scenario}" / "model.json"
    prefix = "617" if family == "00617" and scenario == "005" else family
    return root / ".data" / f"final_{prefix}_{scenario}" / "model.json"


def audit_success(status: dict, summary: dict, expected: int, revision: str,
                  executable_hash: str, limit: float) -> list[str]:
    failures = []
    if not status.get("success") or status.get("stage") != "complete":
        failures.append("not complete")
    if summary.get("git_revision") != revision:
        failures.append("revision mismatch")
    for field in ("base_exact_executable_sha256", "fast_screen_executable_sha256"):
        if summary.get(field) != executable_hash:
            failures.append(field + " mismatch")
    wall = summary.get("total_wall_seconds")
    if not isinstance(wall, (int, float)) or not math.isfinite(wall) or wall > limit or wall < 0:
        failures.append("end-to-end deadline")
    if (summary.get("contingency_count") != expected or
            status.get("completed_contingency_count") != expected):
        failures.append("incomplete contingency count")
    cert = summary.get("official_evaluation_certificate", {})
    if (cert.get("complete_label_set") is not True or
            cert.get("expected_detail_count") != expected + 1 or
            cert.get("observed_unique_detail_count") != expected + 1 or
            cert.get("reported_infeasibility") != 0.0):
        failures.append("official certificate")
    evaluation = summary.get("official_evaluation", {})
    obj = evaluation.get("obj")
    if (evaluation.get("infeas") != 0.0 or not isinstance(obj, (int, float))
            or not math.isfinite(obj)):
        failures.append("official objective/infeasibility")
    residuals = [summary.get("max_independent_contingency_residual"),
                 summary.get("base", {}).get("base_validation", {}).get("max_residual")]
    if any(not isinstance(x, (int, float)) or not math.isfinite(x)
           or x > 1e-5 or x < 0 for x in residuals):
        failures.append("independent residual")
    return failures


def arguments(config: dict, family: str, scenario: str, output: Path) -> list[str]:
    data = Path(config["data_repository"])
    source = Path(config["source_root"]) / f"C2FEN{family}" / f"scenario_{scenario}"
    large = family == "19402"
    command = [config["python"], "-B", str(REPO / "scripts/run_experiment.py"),
        "--case-json", str(normalized_case(data, family, scenario)),
        "--case-dir", str(source), "--output-dir", str(output),
        "--executable", str(REPO / ".build-native/gravityx_go2"),
        "--fast-screen-executable", str(REPO / ".build-native/gravityx_go2"),
        "--workers", "4", "--post-screen-workers", "8",
        "--fast-workers", "24" if large else "16", "--base-timeout", "110",
        "--contingency-timeout", "300", "--total-time-limit", str(config["total_time_limit"]),
        "--minimum-free-space-gb", str(config["minimum_free_space_gib"]),
        "--evaluation-reserve", "3", "--finalization-reserve", "2",
        "--validated-source-base", "--base-sparse-economic-refinement-seconds", "10",
        "--base-sparse-ac-economic-refinement-seconds", "40",
        "--fast-power-flow-screen", "--two-stage-contingency-screen",
        "--linearized-contingency-fallback", "--cpp-solution-writer",
        "--fast-screen-affinity-schedule", "--wsl-stage-worker-inputs",
        "--wsl-fast-screen-scratch", "--compact-final-summary",
        "--evaluator", str(REPO / "scripts/fast_official_evaluator.py"),
        "--vendor-evaluator-reference", config["vendor_evaluator"],
        "--streaming-serial-evaluation-shards", "47" if large else "12",
        "--streaming-evaluation-processes", "3" if large else "2",
        "--post-screen-streaming-evaluation-processes", "12" if large else "4",
        "--streaming-evaluation-idle-screen-ramp", "--streaming-evaluation-completion-order-shards",
        "--streaming-persistent-evaluator-processes", "--streaming-evaluation-tail-shard-sizes",
        "96,64,48,32,16" if large else "32,16,8,4"]
    profile = config.get("screen_profiles", {}).get(f"{family}/{scenario}")
    if profile:
        command += ["--fast-screen-heavy-profile", str(REPO / profile),
                    "--fast-screen-heavy-workers", "4"]
    return command


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=REPO / "config/reliability_suite.json")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--family")
    parser.add_argument("--scenario")
    opts = parser.parse_args()
    config = json.loads(opts.config.read_text())
    for path in (REPO, opts.output, Path(config["data_repository"]), Path(config["source_root"])):
        reject_onedrive(path)
    if opts.output.exists():
        raise RuntimeError("Output exists; refusing duplicate campaign")
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip()
    dirty = subprocess.check_output(["git", "diff", "HEAD", "--"], cwd=REPO, text=True)
    if dirty:
        raise RuntimeError("Commit tracked changes before running")
    executable_hash = digest(REPO / ".build-native/gravityx_go2")
    entries = []
    for network in config["families"]:
        family = network["family"]
        if opts.family and opts.family != family:
            continue
        for scenario in network["scenarios"]:
            if opts.scenario and opts.scenario != scenario:
                continue
            source = Path(config["source_root"]) / f"C2FEN{family}" / f"scenario_{scenario}"
            model = normalized_case(Path(config["data_repository"]), family, scenario)
            for file in (model, source / "case.raw", source / "case.con", source / "case.json"):
                if not file.is_file():
                    raise FileNotFoundError(file)
            expected = len(read_contingency_blocks(source / "case.con"))
            profile = config.get("screen_profiles", {}).get(f"{family}/{scenario}")
            if profile:
                reject_onedrive(REPO / profile)
            entries.append({"family": family, "buses": network["buses"], "scenario": scenario,
                            "source_contingencies": expected, "model_sha256": digest(model),
                            "source_sha256": {name: digest(source / name)
                                for name in ("case.raw", "case.json", "case.con")},
                            "screen_profile_sha256": digest(REPO / profile) if profile else None})
    if not entries:
        raise RuntimeError("No scenarios selected")
    opts.output.mkdir(parents=True)
    manifest = {"revision": revision, "executable_sha256": executable_hash,
                "config_sha256": digest(opts.config), "config": config, "scenarios": entries}
    write_json(opts.output / "manifest.json", manifest)
    results = []
    for entry in entries:
        if shutil.disk_usage(opts.output).free < config["minimum_free_space_gib"] * 1024**3:
            print("STOP: pre-run storage floor", flush=True)
            break
        if digest(REPO / ".build-native/gravityx_go2") != executable_hash:
            raise RuntimeError("Executable changed during campaign")
        source = Path(config["source_root"]) / f"C2FEN{entry['family']}" / f"scenario_{entry['scenario']}"
        model = normalized_case(Path(config["data_repository"]), entry["family"], entry["scenario"])
        if digest(model) != entry["model_sha256"] or any(
                digest(source / name) != sha for name, sha in entry["source_sha256"].items()):
            raise RuntimeError("Source input changed during campaign")
        profile = config.get("screen_profiles", {}).get(f"{entry['family']}/{entry['scenario']}")
        if profile and digest(REPO / profile) != entry["screen_profile_sha256"]:
            raise RuntimeError("Scheduling profile changed during campaign")
        name = f"C2FEN{entry['family']}_s{entry['scenario']}_cold"
        run = opts.output / name
        command = arguments(config, entry["family"], entry["scenario"], run)
        write_json(opts.output / (name + ".command.json"), command)
        print("START " + name, flush=True)
        with (opts.output / (name + ".console.log")).open("w") as log:
            result = subprocess.run(command, cwd=REPO, stdout=log, stderr=subprocess.STDOUT)
        status_path = run / "run_status.json"
        if not status_path.exists():
            raise RuntimeError("Setup failure; see " + str(log.name))
        status = json.loads(status_path.read_text())
        summary_path = run / "run_summary.json"
        summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
        failures = audit_success(status, summary, entry["source_contingencies"], revision,
                                 executable_hash, config["total_time_limit"])
        if result.returncode != 0:
            failures.append("runner nonzero exit")
        record = dict(entry, returncode=result.returncode, passed=not failures,
                      gates=failures, seconds=status.get("total_wall_seconds"),
                      objective=summary.get("official_evaluation", {}).get("obj"),
                      completed=status.get("completed_contingency_count"), error=status.get("error"))
        evidence = archive_run(status_path, opts.output)
        evidence_path = opts.output / (name + ".evidence.json")
        write_json(evidence_path, evidence)
        # Read back the archive and independently verify the primary status hash.
        reread = json.loads(evidence_path.read_text())
        if reread["run_status_sha256"] != digest(status_path):
            raise RuntimeError("Archive verification failed")
        results.append(record)
        write_json(opts.output / "results.json", results)
        print(json.dumps(record), flush=True)
    return 0 if len(results) == len(entries) and all(r["passed"] for r in results) else 1


if __name__ == "__main__":
    sys.exit(main())
