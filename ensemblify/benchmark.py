"""Runtime benchmarking: per-job timing records and compiling them into Parquet tables.

The workflows write one JSON record per job when run with `--config runtime=True`.
`python -m ensemblify.benchmark` (or `pixi run compile-benchmarks`) collects them into
`ensemble_benchmarks.parquet` and `scoring_benchmarks.parquet` in the output directory.
"""

from __future__ import annotations

import argparse
import json
import resource
import time
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

# Recording ---------------------------------------------------------------------


def _cpu_snapshot() -> dict[str, float]:
    # Snakemake runs each `run:` job in its own process, so these counters cover
    # exactly this job: the Python work plus all aligner/RAxML-NG child processes.
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    return {"process": time.process_time(), "children": children.ru_utime + children.ru_stime}


@contextmanager
def record_runtime(path: Path | None, **fields: Any):
    """
    Time the enclosed block and write a JSON record to `path`.

    `fields` are stored as-is (source, dataset, tool, ...). The record is written
    whether the block succeeds or raises; exceptions are re-raised. With
    `path=None` the block runs without recording.
    """
    if path is None:
        yield
        return

    started_at = datetime.now(timezone.utc)
    wall_start = time.perf_counter()
    cpu_start = _cpu_snapshot()
    record = {**fields, "status": "success", "exception_type": None, "exception_message": None}
    try:
        yield
    except Exception as exc:
        record.update(status="failed", exception_type=type(exc).__name__, exception_message=str(exc))
        raise
    finally:
        cpu_end = _cpu_snapshot()
        process = cpu_end["process"] - cpu_start["process"]
        children = cpu_end["children"] - cpu_start["children"]
        record.update(
            started_at=started_at.isoformat(),
            ended_at=datetime.now(timezone.utc).isoformat(),
            wall_seconds=time.perf_counter() - wall_start,
            cpu_seconds={
                "process_seconds": process,
                "children_seconds": children,
                "total_seconds": process + children,
            },
        )
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8") as handle:
            json.dump(record, handle, indent=2, sort_keys=True)
            handle.write("\n")


# Compiling ---------------------------------------------------------------------


@dataclass(frozen=True)
class _Kind:
    glob: str  # relative to the data root
    path_fields: tuple[str, ...]  # directory names in front of the JSON file
    output_name: str


KINDS = {
    "ensemble": _Kind("*/*/*/benchmark.json", ("source", "dataset", "tool"), "ensemble_benchmarks.parquet"),
    "scoring": _Kind("*/*/scoring_benchmark.json", ("source", "dataset"), "scoring_benchmarks.parquet"),
}

_TIMING_COLUMNS = [
    "benchmark_scope",
    "status",
    "threads",
    "started_at",
    "ended_at",
    "wall_seconds",
    "cpu_process_seconds",
    "cpu_children_seconds",
    "cpu_total_seconds",
    "exception_type",
    "exception_message",
    "benchmark_path",
]


def _flatten_record(path: Path, data_root: Path, kind: _Kind, record: dict[str, Any]) -> dict[str, Any]:
    result = dict(record)
    cpu = result.pop("cpu_seconds", {}) or {}
    if isinstance(cpu, dict):
        result["cpu_process_seconds"] = cpu.get("process_seconds")
        result["cpu_children_seconds"] = cpu.get("children_seconds")
        result["cpu_total_seconds"] = cpu.get("total_seconds")
    else:  # older records stored a single number
        result["cpu_process_seconds"] = None
        result["cpu_children_seconds"] = None
        result["cpu_total_seconds"] = cpu
    for field, value in zip(kind.path_fields, path.relative_to(data_root).parts):
        result.setdefault(field, value)
    result["benchmark_path"] = str(path)
    return result


def compile_benchmarks(data_root: Path, kind: str, output_path: Path | None = None) -> tuple[pd.DataFrame, Path]:
    """Collect all benchmark JSON files of `kind` under `data_root` into one Parquet file."""
    spec = KINDS[kind]
    records = []
    for path in sorted(Path(data_root).glob(spec.glob)):
        with path.open(encoding="utf-8") as handle:
            records.append(_flatten_record(path, Path(data_root), spec, json.load(handle)))

    df = pd.DataFrame(records)
    if df.empty:
        df = pd.DataFrame(columns=[*spec.path_fields, *_TIMING_COLUMNS])

    output_path = Path(output_path) if output_path else Path(data_root) / spec.output_name
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(output_path, index=False)
    return df, output_path


def _load_data_root(config_path: Path) -> Path:
    """`general.output` from the config. Relative paths resolve against the current
    directory, the repository root when run through `pixi run`."""
    config = yaml.safe_load(Path(config_path).read_text(encoding="utf-8")) or {}
    output = (config.get("general") or {}).get("output")
    if not output:
        raise ValueError(f"Could not find 'general.output' in {config_path}")
    return Path(output).expanduser().resolve()


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="python -m ensemblify.benchmark",
        description="Compile runtime benchmark JSON files into Parquet tables.",
    )
    parser.add_argument("--kind", choices=[*KINDS, "all"], default="all", help="Which benchmarks to compile (default: all).")
    parser.add_argument("--config", type=Path, default=Path("configs/general.yaml"), help="Config providing general.output (default: configs/general.yaml).")
    parser.add_argument("--data-root", type=Path, default=None, help="Benchmark data root. Defaults to general.output from --config.")
    parser.add_argument("--output", type=Path, default=None, help="Output Parquet path. Only with a single --kind. Defaults to DATA_ROOT/<kind>_benchmarks.parquet.")
    args = parser.parse_args(argv)

    if args.output and args.kind == "all":
        parser.error("--output needs a single --kind (ensemble or scoring)")
    data_root = args.data_root.expanduser() if args.data_root else _load_data_root(args.config)

    for kind in KINDS if args.kind == "all" else [args.kind]:
        df, path = compile_benchmarks(data_root, kind, args.output)
        print(f"Wrote {len(df)} {kind} benchmark rows to {path}")


if __name__ == "__main__":
    main()
