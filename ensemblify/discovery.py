"""Dataset discovery for the Snakemake workflows.

Kept out of the Snakefiles so that a single malformed directory (no sequences, an
unfinished ensemble) is skipped with a warning instead of breaking DAG construction
for every dataset.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path


def _warn(message: str) -> None:
    print(f"WARNING: {message}", file=sys.stderr)


def _subdirs(path: Path) -> list[str]:
    return sorted(p.name for p in path.iterdir() if p.is_dir() and not p.name.startswith("."))


def input_datasets(input_dir: Path, sources: list[str]) -> dict[str, list[str]]:
    """Datasets per source that contain `sequences.fasta`."""
    datasets = {}
    for source in sources:
        source_dir = Path(input_dir) / source
        if not source_dir.is_dir():
            _warn(f"Input source directory {source_dir} does not exist.")
            datasets[source] = []
            continue
        names = _subdirs(source_dir)
        valid = [n for n in names if (source_dir / n / "sequences.fasta").is_file()]
        if skipped := sorted(set(names) - set(valid)):
            _warn(f"Skipping {len(skipped)} dataset(s) in {source_dir} without sequences.fasta: {', '.join(skipped)}")
        datasets[source] = valid
    return datasets


def scorable_datasets(output_dir: Path, sources: list[str], tools: list[str]) -> dict[str, list[str]]:
    """Datasets per source whose ensemble is complete, i.e. every tool has a `done` flag."""
    datasets = {}
    for source in sources:
        source_dir = Path(output_dir) / source
        if not source_dir.is_dir():
            _warn(f"Output source directory {source_dir} does not exist. Run the ensemble workflow first.")
            datasets[source] = []
            continue
        valid = []
        for name in _subdirs(source_dir):
            if not (source_dir / name / "ensemble").is_dir():
                continue
            missing = [t for t in tools if not (source_dir / name / t / "done").is_file()]
            if missing:
                _warn(f"Skipping {source}/{name}: ensemble incomplete (missing {', '.join(missing)}).")
            else:
                valid.append(name)
        datasets[source] = valid
    return datasets


def wildcard_regex(names) -> str:
    """Regex for `wildcard_constraints` matching exactly the given names."""
    return "|".join(re.escape(name) for name in names)
