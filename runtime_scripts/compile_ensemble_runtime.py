"""Compile per-tool ensemble benchmark JSON files into a root-level Parquet file."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd
import yaml


DEFAULT_COLUMNS = [
    "source",
    "dataset",
    "tool",
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


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _load_data_root(config_path: Path) -> Path:
    """Read general.output from the config. Relative paths resolve against the
    repository root, matching how the workflows run (`pixi run` from the root)."""
    config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    output = (config.get("general") or {}).get("output")
    if not output:
        raise ValueError(f"Could not find 'general.output' in {config_path}")
    path = Path(output).expanduser()
    return path if path.is_absolute() else _repo_root() / path


def _flatten_record(
    path: Path, data_root: Path, record: dict[str, Any]
) -> dict[str, Any]:
    result = dict(record)
    cpu_seconds = result.pop("cpu_seconds", {}) or {}
    if isinstance(cpu_seconds, dict):
        result["cpu_process_seconds"] = cpu_seconds.get("process_seconds")
        result["cpu_children_seconds"] = cpu_seconds.get("children_seconds")
        result["cpu_total_seconds"] = cpu_seconds.get("total_seconds")
    else:
        result["cpu_process_seconds"] = None
        result["cpu_children_seconds"] = None
        result["cpu_total_seconds"] = cpu_seconds

    parts = path.relative_to(data_root).parts
    if len(parts) >= 4:
        result.setdefault("source", parts[0])
        result.setdefault("dataset", parts[1])
        result.setdefault("tool", parts[2])
    result["benchmark_path"] = str(path)
    return result


def compile_benchmarks(data_root: Path, output_path: Path) -> pd.DataFrame:
    records = []
    for benchmark_path in sorted(data_root.glob("*/*/*/benchmark.json")):
        with benchmark_path.open(encoding="utf-8") as handle:
            record = json.load(handle)
        records.append(_flatten_record(benchmark_path, data_root, record))

    df = pd.DataFrame(records)
    if df.empty:
        df = pd.DataFrame(columns=DEFAULT_COLUMNS)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(output_path, index=False)
    return df


def parse_args() -> argparse.Namespace:
    default_config = _repo_root() / "configs" / "general.yaml"
    parser = argparse.ArgumentParser(
        description="Compile OUT_DIR/source/dataset/tool/benchmark.json files."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=default_config,
        help="Path to configs/general.yaml. Used only when --data-root is omitted.",
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        default=None,
        help="Benchmark data root. Defaults to general.output from --config.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output parquet path. Defaults to DATA_ROOT/ensemble_benchmarks.parquet.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data_root = (args.data_root or _load_data_root(args.config)).expanduser()
    output_path = args.output or data_root / "ensemble_benchmarks.parquet"
    df = compile_benchmarks(data_root, output_path)
    print(f"Wrote {len(df)} ensemble benchmark rows to {output_path}")


if __name__ == "__main__":
    main()
