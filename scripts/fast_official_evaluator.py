#!/usr/bin/env python3
"""Run the GO2 evaluator equations with a one-pass generated-solution parser.

The official evaluator reads each solution section with a separate pandas
``read_csv`` call.  Each call starts at the beginning of the same text file,
so large generated submissions are scanned seven times per contingency.  This
launcher replaces only ``CaseSolution.read`` and its subsequent DataFrame-to-
array copy.  The official case loader, equations, tolerances, objective,
infeasibility logic, and output writers remain the referenced vendor module.

This parser intentionally accepts only the canonical section order emitted by
GravityX's C++ solution writer.  It validates every row count and source key
before assigning values to the vendor evaluator's existing NumPy arrays.
"""

from __future__ import annotations

import csv
import importlib.util
import os
from pathlib import Path
import sys
from types import ModuleType
from typing import Any, Iterable


VENDOR_EVALUATOR_ENVIRONMENT_VARIABLE = "GRAVITYX_VENDOR_EVALUATOR"


def _rows(
    lines: list[str],
    cursor: int,
    section: str,
    columns: str,
    count: int,
) -> tuple[Iterable[list[str]], int]:
    if count < 0:
        raise ValueError(f"negative row count for {section}: {count}")
    end = cursor + 2 + count
    if end > len(lines):
        raise ValueError(
            f"truncated {section}: expected {count} rows, "
            f"have at most {max(0, len(lines) - cursor - 2)}"
        )
    if lines[cursor].strip() != section:
        raise ValueError(
            f"expected {section} at line {cursor + 1}, "
            f"found {lines[cursor].strip()!r}"
        )
    if lines[cursor + 1].strip() != columns:
        raise ValueError(
            f"invalid {section} column header: {lines[cursor + 1].strip()!r}"
        )
    reader = csv.reader(
        lines[cursor + 2 : end],
        delimiter=",",
        quotechar="'",
        skipinitialspace=True,
    )
    return reader, end


def _require_fields(fields: list[str], expected: int, section: str) -> None:
    if len(fields) != expected:
        raise ValueError(
            f"{section} row has {len(fields)} fields; expected {expected}: {fields!r}"
        )


def read_generated_solution(self: Any, file_name: str | os.PathLike[str]) -> None:
    """Populate one vendor ``CaseSolution`` directly from canonical text."""
    self._gravityx_arrays_loaded = False
    with open(file_name, "r", encoding="utf-8", newline=None) as stream:
        lines = stream.readlines()

    cursor = 0

    rows, cursor = _rows(
        lines, cursor, "--bus section", "i, v, theta", self.num_bus_read
    )
    bus_indices: list[int] = []
    bus_vm: list[float] = []
    bus_va: list[float] = []
    for fields in rows:
        _require_fields(fields, 3, "bus")
        bus_indices.append(self.bus_map[int(fields[0])])
        bus_vm.append(float(fields[1]))
        bus_va.append(float(fields[2]))
    self.bus_volt_mag.fill(1.0)
    self.bus_volt_ang.fill(0.0)
    self.bus_volt_mag[bus_indices] = bus_vm
    self.bus_volt_ang[bus_indices] = bus_va

    rows, cursor = _rows(
        lines, cursor, "--load section", "i, id, t", self.num_load_read
    )
    load_indices: list[int] = []
    load_values: list[float] = []
    for fields in rows:
        _require_fields(fields, 3, "load")
        load_indices.append(self.load_map[(int(fields[0]), fields[1].strip())])
        load_values.append(float(fields[2]))
    self.load_t.fill(0.0)
    self.load_t[load_indices] = load_values

    rows, cursor = _rows(
        lines,
        cursor,
        "--generator section",
        "i, id, p, q, x",
        self.num_gen_read,
    )
    generator_indices: list[int] = []
    generator_p: list[float] = []
    generator_q: list[float] = []
    generator_x: list[float] = []
    for fields in rows:
        _require_fields(fields, 5, "generator")
        generator_indices.append(
            self.gen_map[(int(fields[0]), fields[1].strip())]
        )
        generator_p.append(float(fields[2]))
        generator_q.append(float(fields[3]))
        generator_x.append(float(fields[4]))
    self.gen_pow_real.fill(0.0)
    self.gen_pow_imag.fill(0.0)
    self.gen_xon.fill(0.0)
    self.gen_pow_real[generator_indices] = generator_p
    self.gen_pow_imag[generator_indices] = generator_q
    self.gen_xon[generator_indices] = generator_x

    rows, cursor = _rows(
        lines,
        cursor,
        "--line section",
        "iorig, idest, id, x",
        self.num_line_read,
    )
    line_indices: list[int] = []
    line_values: list[float] = []
    for fields in rows:
        _require_fields(fields, 4, "line")
        line_indices.append(
            self.line_map[
                (int(fields[0]), int(fields[1]), fields[2].strip())
            ]
        )
        line_values.append(float(fields[3]))
    self.line_xsw.fill(0.0)
    self.line_xsw[line_indices] = line_values

    rows, cursor = _rows(
        lines,
        cursor,
        "--transformer section",
        "iorig, idest, id, x, xst",
        self.num_xfmr_read,
    )
    transformer_indices: list[int] = []
    transformer_x: list[float] = []
    transformer_xst: list[float] = []
    for fields in rows:
        _require_fields(fields, 5, "transformer")
        transformer_indices.append(
            self.xfmr_map[
                (int(fields[0]), int(fields[1]), fields[2].strip())
            ]
        )
        transformer_x.append(float(fields[3]))
        transformer_xst.append(float(fields[4]))
    self.xfmr_xsw.fill(0.0)
    self.xfmr_xst.fill(0.0)
    self.xfmr_xsw[transformer_indices] = transformer_x
    self.xfmr_xst[transformer_indices] = transformer_xst

    rows, cursor = _rows(
        lines,
        cursor,
        "--switched shunt section",
        "i, xst1, xst2, xst3, xst4, xst5, xst6, xst7, xst8",
        self.num_swsh_read,
    )
    shunt_indices: list[int] = []
    shunt_values: list[list[float]] = []
    for fields in rows:
        if not 1 <= len(fields) <= 9:
            raise ValueError(
                f"switched-shunt row has {len(fields)} fields; expected 1 to 9"
            )
        shunt_indices.append(self.swsh_map[int(fields[0])])
        values = [float(value) if value.strip() else 0.0 for value in fields[1:]]
        shunt_values.append(values + [0.0] * (8 - len(values)))
    self.swsh_xst.fill(0.0)
    if shunt_indices:
        self.swsh_xst[shunt_indices, :] = shunt_values

    if any(line.strip() for line in lines[cursor:]):
        raise ValueError(
            f"unexpected nonempty content after canonical solution sections in {file_name}"
        )
    self._gravityx_arrays_loaded = True


def skip_dataframe_copy(self: Any) -> None:
    if not getattr(self, "_gravityx_arrays_loaded", False):
        raise RuntimeError("fast evaluator arrays were not loaded before copy step")


def load_vendor_evaluator(path: Path) -> ModuleType:
    if not path.is_file():
        raise FileNotFoundError(f"vendor evaluator does not exist: {path}")
    package_root = path.parent.parent
    for candidate in (package_root, path.parent):
        if str(candidate) not in sys.path:
            sys.path.insert(0, str(candidate))
    spec = importlib.util.spec_from_file_location(
        "gravityx_referenced_go2_evaluator", path
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load vendor evaluator: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def install_fast_parser(module: ModuleType) -> None:
    module.CaseSolution.read = read_generated_solution
    module.CaseSolution.set_arrays_from_dfs = skip_dataframe_copy


def main() -> int:
    raw_path = os.environ.get(VENDOR_EVALUATOR_ENVIRONMENT_VARIABLE)
    if not raw_path:
        raise RuntimeError(
            f"{VENDOR_EVALUATOR_ENVIRONMENT_VARIABLE} must name the referenced "
            "official evaluator"
        )
    module = load_vendor_evaluator(Path(raw_path).resolve())
    install_fast_parser(module)
    module.main()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
