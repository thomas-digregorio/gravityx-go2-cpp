#!/usr/bin/env python3
"""Run the GO2 evaluator equations with a one-pass generated-solution parser.

The official evaluator reads each solution section with a separate pandas
``read_csv`` call.  Each call starts at the beginning of the same text file,
so large generated submissions are scanned seven times per contingency.  This
launcher replaces ``CaseSolution.read`` and its subsequent DataFrame-to-array
copy.  It also retains each already-written per-case summary in memory so the
vendor final writer does not reread the same JSON files, and suppresses the
vendor's per-function timing prints.  The official case loader, equations,
tolerances, objective, infeasibility logic, per-case detail writer, and final
output writers remain the referenced vendor module.

This parser intentionally accepts only the canonical section order emitted by
GravityX's C++ solution writer.  It validates every row count and source key
before assigning values to the vendor evaluator's existing NumPy arrays.
"""

from __future__ import annotations

import csv
import copy
import importlib.util
import os
from pathlib import Path
import sys
from types import ModuleType
from typing import Any, Iterable
import warnings

import numpy as np


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


def _read_generated_solution_csv(
    self: Any, file_name: str | os.PathLike[str]
) -> None:
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


class _NumericParserFallback(Exception):
    """Signal that a valid canonical file needs the generic string-ID parser."""


def _numeric_matrix(
    lines: list[str],
    cursor: int,
    section: str,
    columns: str,
    count: int,
    width: int,
    cache_owner: Any | None = None,
    cache_key: str | None = None,
) -> tuple[np.ndarray, int]:
    """Parse a fixed-width numeric canonical section in NumPy's C loop."""

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
    rows = lines[cursor + 2 : end]
    if cache_owner is not None and cache_key is not None:
        cache = getattr(cache_owner, "_gravityx_exact_matrix_cache", None)
        if cache is not None:
            cached = cache.get(cache_key)
            if cached is not None and cached["rows"] == rows:
                return cached["matrix"], end
    block = "".join(rows).replace(",", " ")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        try:
            values = np.fromstring(block, dtype=np.float64, sep=" ")
        except ValueError as error:
            raise _NumericParserFallback(
                f"{section} is not an all-numeric section"
            ) from error
    expected = count * width
    if values.size != expected:
        raise _NumericParserFallback(
            f"{section} is not an all-numeric {width}-column section"
        )
    matrix = values.reshape((count, width))
    if cache_owner is not None and cache_key is not None:
        cache = getattr(cache_owner, "_gravityx_exact_matrix_cache", None)
        if cache is None:
            cache = {}
            cache_owner._gravityx_exact_matrix_cache = cache
        cache[cache_key] = {"rows": rows, "matrix": matrix}
    return matrix, end


def _integer_key_columns(values: np.ndarray, section: str) -> np.ndarray:
    keys = values.astype(np.int64)
    if not np.array_equal(values, keys):
        raise _NumericParserFallback(f"{section} source keys are not integers")
    return keys


def _cached_numeric_indices(
    self: Any,
    section: str,
    source_keys: np.ndarray,
    mapping: dict[Any, int],
    key_builder: Any,
) -> np.ndarray:
    """Validate every source key while amortizing immutable map lookups."""

    cache = getattr(self, "_gravityx_numeric_layout", None)
    if cache is None:
        cache = {}
        self._gravityx_numeric_layout = cache
    cached = cache.get(section)
    if cached is not None and np.array_equal(
        cached["source_keys"], source_keys
    ):
        return cached["indices"]
    try:
        indices = np.fromiter(
            (mapping[key_builder(row)] for row in source_keys),
            dtype=np.int64,
            count=source_keys.shape[0],
        )
    except (KeyError, TypeError, ValueError) as error:
        raise _NumericParserFallback(
            f"{section} numeric source keys do not match the vendor map"
        ) from error
    cache[section] = {
        "source_keys": source_keys.copy(),
        "indices": indices,
    }
    return indices


def _single_missing_row(
    full_rows: list[str], subset_rows: list[str], section: str
) -> int:
    if len(subset_rows) + 1 != len(full_rows):
        raise ValueError(
            f"canonical {section} row count is not full or single-outage"
        )
    missing = len(subset_rows)
    for position, (observed, expected) in enumerate(
        zip(subset_rows, full_rows)
    ):
        if observed != expected:
            missing = position
            break
    if subset_rows[missing:] != full_rows[missing + 1 :]:
        raise ValueError(
            f"canonical {section} rows are not an exact single-row subset"
        )
    return missing


def _static_numeric_matrix(
    self: Any,
    lines: list[str],
    cursor: int,
    section: str,
    columns: str,
    count: int,
    full_count: int,
    width: int,
    key_width: int,
    cache_key: str,
    mapping: dict[Any, int],
    key_builder: Any,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Reuse immutable line/transformer rows across single outages."""

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
    rows = lines[cursor + 2 : end]
    cache = getattr(self, "_gravityx_static_section_cache", None)
    cached = cache.get(cache_key) if cache is not None else None
    if cached is not None and count == full_count:
        if cached["rows"] != rows:
            raise ValueError(f"canonical full {section} rows changed")
        return cached["matrix"], cached["indices"], end
    if cached is not None and count + 1 == full_count:
        missing = _single_missing_row(cached["rows"], rows, cache_key)
        matrix = np.concatenate(
            (cached["matrix"][:missing], cached["matrix"][missing + 1 :]),
            axis=0,
        )
        indices = np.concatenate(
            (cached["indices"][:missing], cached["indices"][missing + 1 :])
        )
        return matrix, indices, end

    matrix, parsed_end = _numeric_matrix(
        lines, cursor, section, columns, count, width
    )
    source_keys = _integer_key_columns(matrix[:, :key_width], cache_key)
    indices = _cached_numeric_indices(
        self, cache_key, source_keys, mapping, key_builder
    )
    if count == full_count:
        if cache is None:
            cache = {}
            self._gravityx_static_section_cache = cache
        cache[cache_key] = {
            "rows": rows,
            "matrix": matrix,
            "indices": indices,
        }
    return matrix, indices, parsed_end


def _read_numeric_generated_solution(
    self: Any, file_name: str | os.PathLike[str]
) -> None:
    """Read the all-numeric GO2 canonical form with vectorized conversion."""

    with open(file_name, "r", encoding="utf-8", newline=None) as stream:
        lines = stream.readlines()
    cursor = 0

    bus, cursor = _numeric_matrix(
        lines,
        cursor,
        "--bus section",
        "i, v, theta",
        self.num_bus_read,
        3,
    )
    bus_keys = _integer_key_columns(bus[:, :1], "bus")
    bus_indices = _cached_numeric_indices(
        self, "bus", bus_keys, self.bus_map, lambda row: int(row[0])
    )

    load, cursor = _numeric_matrix(
        lines,
        cursor,
        "--load section",
        "i, id, t",
        self.num_load_read,
        3,
        self,
        "load",
    )
    load_keys = _integer_key_columns(load[:, :2], "load")
    load_indices = _cached_numeric_indices(
        self,
        "load",
        load_keys,
        self.load_map,
        lambda row: (int(row[0]), str(int(row[1]))),
    )

    generator, cursor = _numeric_matrix(
        lines,
        cursor,
        "--generator section",
        "i, id, p, q, x",
        self.num_gen_read,
        5,
    )
    generator_keys = _integer_key_columns(generator[:, :2], "generator")
    generator_indices = _cached_numeric_indices(
        self,
        "generator",
        generator_keys,
        self.gen_map,
        lambda row: (int(row[0]), str(int(row[1]))),
    )

    line, line_indices, cursor = _static_numeric_matrix(
        self,
        lines,
        cursor,
        "--line section",
        "iorig, idest, id, x",
        self.num_line_read,
        len(self.line_map),
        4,
        3,
        "line",
        self.line_map,
        lambda row: (int(row[0]), int(row[1]), str(int(row[2]))),
    )

    transformer, transformer_indices, cursor = _static_numeric_matrix(
        self,
        lines,
        cursor,
        "--transformer section",
        "iorig, idest, id, x, xst",
        self.num_xfmr_read,
        len(self.xfmr_map),
        5,
        3,
        "transformer",
        self.xfmr_map,
        lambda row: (int(row[0]), int(row[1]), str(int(row[2]))),
    )

    shunt_end = cursor + 2 + self.num_swsh_read
    if shunt_end > len(lines):
        raise ValueError(
            "truncated --switched shunt section: "
            f"expected {self.num_swsh_read} rows"
        )
    if lines[cursor].strip() != "--switched shunt section":
        raise ValueError(
            f"expected --switched shunt section at line {cursor + 1}"
        )
    expected_shunt_header = (
        "i, xst1, xst2, xst3, xst4, xst5, xst6, xst7, xst8"
    )
    if lines[cursor + 1].strip() != expected_shunt_header:
        raise ValueError(
            "invalid --switched shunt section column header: "
            f"{lines[cursor + 1].strip()!r}"
        )
    shunt_rows = lines[cursor + 2 : shunt_end]
    shunt_cache = getattr(self, "_gravityx_shunt_cache", None)
    if shunt_cache is not None and shunt_cache["rows"] == shunt_rows:
        shunt_values = shunt_cache["values"]
        shunt_indices = shunt_cache["indices"]
    else:
        shunt_keys = np.empty((self.num_swsh_read, 1), dtype=np.int64)
        shunt_values = np.zeros((self.num_swsh_read, 8), dtype=np.float64)
        for position, raw_line in enumerate(shunt_rows):
            fields = raw_line.split(",")
            if not 1 <= len(fields) <= 9:
                raise ValueError(
                    "switched-shunt row has "
                    f"{len(fields)} fields; expected 1 to 9"
                )
            try:
                shunt_keys[position, 0] = int(fields[0])
                if len(fields) > 1:
                    shunt_values[position, : len(fields) - 1] = [
                        float(value) if value.strip() else 0.0
                        for value in fields[1:]
                    ]
            except ValueError as error:
                raise _NumericParserFallback(
                    "switched-shunt section is not all numeric"
                ) from error
        shunt_indices = _cached_numeric_indices(
            self,
            "switched shunt",
            shunt_keys,
            self.swsh_map,
            lambda row: int(row[0]),
        )
        self._gravityx_shunt_cache = {
            "rows": shunt_rows,
            "values": shunt_values,
            "indices": shunt_indices,
        }
    cursor = shunt_end
    if any(line.strip() for line in lines[cursor:]):
        raise ValueError(
            f"unexpected nonempty content after canonical solution sections in {file_name}"
        )

    self.bus_volt_mag.fill(1.0)
    self.bus_volt_ang.fill(0.0)
    self.bus_volt_mag[bus_indices] = bus[:, 1]
    self.bus_volt_ang[bus_indices] = bus[:, 2]
    self.load_t.fill(0.0)
    self.load_t[load_indices] = load[:, 2]
    self.gen_pow_real.fill(0.0)
    self.gen_pow_imag.fill(0.0)
    self.gen_xon.fill(0.0)
    self.gen_pow_real[generator_indices] = generator[:, 2]
    self.gen_pow_imag[generator_indices] = generator[:, 3]
    self.gen_xon[generator_indices] = generator[:, 4]
    self.line_xsw.fill(0.0)
    self.line_xsw[line_indices] = line[:, 3]
    self.xfmr_xsw.fill(0.0)
    self.xfmr_xst.fill(0.0)
    self.xfmr_xsw[transformer_indices] = transformer[:, 3]
    self.xfmr_xst[transformer_indices] = transformer[:, 4]
    self.swsh_xst.fill(0.0)
    if shunt_indices.size:
        self.swsh_xst[shunt_indices, :] = shunt_values
    self._gravityx_arrays_loaded = True


def read_generated_solution(self: Any, file_name: str | os.PathLike[str]) -> None:
    """Populate vendor arrays, using a numeric fast path with exact fallback."""

    self._gravityx_arrays_loaded = False
    try:
        _read_numeric_generated_solution(self, file_name)
    except _NumericParserFallback:
        _read_generated_solution_csv(self, file_name)


def skip_dataframe_copy(self: Any) -> None:
    if not getattr(self, "_gravityx_arrays_loaded", False):
        raise RuntimeError("fast evaluator arrays were not loaded before copy step")


def install_in_memory_summary_aggregation(module: ModuleType) -> None:
    """Feed the unchanged vendor final writer from summaries it just wrote.

    The serial vendor path writes one ``eval_detail_<label>.json`` file per
    case, then rereads every one of those files solely to reconstruct
    ``summary_all_cases``.  Capture an independent deep copy at the same point
    as each detail write and replace only that redundant disk reread.  The
    individual detail files and all downstream vendor summary writers are
    unchanged.
    """

    original_write_detail = module.Evaluation.write_detail

    def write_detail_and_capture(
        self: Any,
        path: str | os.PathLike[str],
        case: str,
        detail_csv: bool = False,
        detail_json: bool = False,
    ) -> Any:
        result = original_write_detail(
            self,
            path,
            case,
            detail_csv=detail_csv,
            detail_json=detail_json,
        )
        if detail_json:
            self.summary_all_cases[module.clean_string(case)] = copy.deepcopy(
                self.summary
            )
        return result

    def summaries_already_loaded(self: Any, path: str | os.PathLike[str]) -> None:
        del path
        expected = {"BASECASE", *self.ctg_label}
        observed = set(self.summary_all_cases)
        if observed != expected:
            missing = sorted(expected - observed)
            extra = sorted(observed - expected)
            raise RuntimeError(
                "in-memory evaluator summary-set mismatch: "
                f"missing={missing}, extra={extra}"
            )

    module.Evaluation.write_detail = write_detail_and_capture
    module.Evaluation.json_to_summary_all_cases = summaries_already_loaded


def suppress_verbose_timing_output(module: ModuleType) -> None:
    """Suppress flushed diagnostic stdout without changing evaluator results."""

    module.print = lambda *args, **kwargs: None


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
    install_in_memory_summary_aggregation(module)
    suppress_verbose_timing_output(module)


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
