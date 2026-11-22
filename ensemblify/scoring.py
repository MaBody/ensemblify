"""Scoring of MSA ensembles, shared by the bench and runtime scoring workflows."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import AlignIO

from aldiscore.datastructures.alignment import Alignment
from aldiscore.datastructures.ensemble import Ensemble
from aldiscore.scoring import pairwise, pythia

from ensemblify.utils import infer_data_type

# RAxML-NG, and with it Pythia, needs at least four sequences to build a tree.
MIN_PYTHIA_SEQUENCES = 4

STAT_NAMES = ["mean", "std", "min", "p25", "p50", "p75", "max"]


def descriptive_stats(scores) -> dict[str, float]:
    """
    Summary statistics of a score array, ignoring NaN entries.

    NaN marks individual failures (e.g. one alignment RAxML-NG could not process),
    which must not discard the scores of the other alignments. `n` is the number of
    valid scores. `std` needs at least two values; all other entries need one.
    """
    scores = np.asarray(scores, dtype=float).ravel()
    valid = scores[~np.isnan(scores)]
    n = len(valid)
    if n == 0:
        stats = {name: np.nan for name in STAT_NAMES}
    else:
        stats = {
            "mean": np.mean(valid),
            "std": np.std(valid, ddof=1) if n > 1 else np.nan,
            "min": np.min(valid),
            "p25": np.percentile(valid, 25),
            "p50": np.percentile(valid, 50),
            "p75": np.percentile(valid, 75),
            "max": np.max(valid),
        }
    stats["n"] = n
    return stats


def load_ensemble(ensemble_dir: Path) -> Ensemble:
    """Load all `*.fasta` alignments in `ensemble_dir`, in sorted file name order."""
    ensemble_dir = Path(ensemble_dir)
    files = sorted(p for p in ensemble_dir.iterdir() if p.suffix == ".fasta" and p.is_file())
    if not files:
        raise FileNotFoundError(f"No .fasta alignments in {ensemble_dir}")
    return Ensemble([Alignment(AlignIO.read(f, "fasta")) for f in files])


def find_reference(*dataset_dirs: Path) -> Path | None:
    """Return the first existing `reference.fasta` in the given dataset directories."""
    for directory in dataset_dirs:
        candidate = Path(directory) / "reference.fasta"
        if candidate.is_file():
            return candidate
    return None


def score_ensemble(
    ensemble: Ensemble,
    reference_file: Path | None = None,
    raxml_path: Path | None = None,
    threads: int = 1,
) -> pd.DataFrame:
    """
    Score an ensemble. Returns one row per method (`dpos`, `pythia`, and `dpos_ref`
    if a reference alignment is given) with the columns of `descriptive_stats` plus
    `datatype`.
    """
    stats = {}
    measure = pairwise.DPosDistance(out_format="flat")
    stats["dpos"] = descriptive_stats(measure.compute(ensemble))
    stats["pythia"] = _pythia_stats(ensemble, raxml_path, threads)

    if reference_file is not None:
        reference = Alignment(AlignIO.read(reference_file, "fasta"))
        stats["dpos_ref"] = descriptive_stats(measure.compute(ensemble, reference))

    df = pd.DataFrame.from_dict(stats, orient="index")
    df.index.name = "method"
    df["datatype"] = infer_data_type(ensemble.dataset.records)
    return df


def _pythia_stats(ensemble: Ensemble, raxml_path: Path | None, threads: int) -> dict:
    n_sequences = len(ensemble.dataset.records)
    if n_sequences < MIN_PYTHIA_SEQUENCES:
        return descriptive_stats([])
    try:
        scores = pythia.compute_pythia_difficulty(ensemble, raxml_path, threads=threads)
    except Exception as e:  # noqa: BLE001 - keep scoring the remaining methods/datasets
        # aldiscore already maps per-alignment RAxML-NG failures to NaN. Anything that
        # reaches this point is unexpected, so report it loudly instead of hiding it.
        print(f"WARNING: Pythia failed ({type(e).__name__}: {e}); reporting NaN.", file=sys.stderr)
        return descriptive_stats([])
    return descriptive_stats(scores)
