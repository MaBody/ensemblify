"""End-to-end smoke test: run both Snakemake workflows on the demo data.

Uses a reduced ensemble (MAFFT, Clustal Omega, Muscle5) so it finishes quickly and
does not depend on the optional Muscle3 environment. Deselect with `-m "not workflow"`.
"""

import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent

pytestmark = pytest.mark.workflow

ENSEMBLE = {
    "MAFFT_FFT-NS-2": {
        "aligner": "MAFFT",
        "type": "default",
        "params": {"--op": [1.4025, 1.6575]},
    },
    "ClustalO": {
        "aligner": "ClustalO",
        "type": "default",
        "params": {"--max-hmm-iterations": [1, 2]},
    },
    "Muscle5": {
        "aligner": "Muscle5",
        "type": "Muscle5",
        "params": {"-replicates": 1},
    },
}


@pytest.fixture(scope="module")
def workdir(tmp_path_factory, demo_input) -> Path:
    """A Snakemake working directory with its own configs/ pointing at the demo data."""
    wd = tmp_path_factory.mktemp("workflow")
    configs = wd / "configs"
    configs.mkdir()

    aligners = yaml.safe_load((REPO_ROOT / "configs" / "aligners.yaml").read_text())
    for aligner in aligners["aligners"].values():
        aligner["threads"] = 2
    # The bioconda osx-arm64 build of Clustal Omega lacks OpenMP and rejects --threads > 1.
    aligners["aligners"]["ClustalO"]["threads"] = 1

    (configs / "aligners.yaml").write_text(yaml.safe_dump(aligners))
    (configs / "ensemble.yaml").write_text(yaml.safe_dump({"ensemble": ENSEMBLE}))
    (configs / "general.yaml").write_text(
        yaml.safe_dump(
            {
                "general": {
                    "input": str(demo_input),
                    "output": str(wd / "results"),
                    "sources": ["demo"],
                    "raxml": "raxml-ng",
                }
            }
        )
    )
    return wd


def _snakemake(workflow: str, workdir: Path) -> None:
    snakemake = shutil.which("snakemake")
    assert snakemake is not None
    proc = subprocess.run(
        [snakemake, "-s", str(REPO_ROOT / "workflow" / workflow), "-d", str(workdir), "-c", "2"],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stdout[-3000:] + proc.stderr[-3000:]


def test_ensemble_workflow(workdir):
    _snakemake("ensemble_bench.smk", workdir)
    for dataset in ["protein", "dna"]:
        ensemble_dir = workdir / "results" / "demo" / dataset / "ensemble"
        files = sorted(p.name for p in ensemble_dir.iterdir())
        assert sum(f.startswith("MAFFT_FFT-NS-2.") for f in files) == 2
        assert sum(f.startswith("ClustalO.") for f in files) == 2
        assert sum(f.startswith("Muscle5.") for f in files) >= 1


def test_scoring_workflow(workdir):
    _snakemake("scoring_bench.smk", workdir)
    stats = pd.read_parquet(workdir / "results" / "demo" / "stats.parquet")
    datasets = set(stats.index.get_level_values("dataset"))
    methods = set(stats.index.get_level_values("method"))
    assert datasets == {"protein", "dna"}
    assert {"dpos", "pythia"} <= methods
    dpos = stats.xs("dpos", level="method")
    assert dpos["mean"].between(0, 1).all()
    assert stats.xs("pythia", level="method")["mean"].notna().all()
    datatypes = stats.xs("dpos", level="method")["datatype"].droplevel("source")
    assert datatypes.to_dict() == {"protein": "AA", "dna": "DNA"}
    # demo/data/demo/protein ships a reference.fasta next to sequences.fasta, dna does not.
    ref = stats.xs("dpos_ref", level="method").droplevel("source")
    assert list(ref.index) == ["protein"]
    n_alignments = len(list((workdir / "results" / "demo" / "protein" / "ensemble").glob("*.fasta")))
    assert ref.loc["protein", "n"] == n_alignments  # one distance per alignment
    assert dpos.loc[("demo", "protein"), "n"] == n_alignments * (n_alignments - 1) // 2  # all pairs

