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


def _make_workdir(tmp_path_factory, demo_input) -> Path:
    """A Snakemake working directory with its own configs/ pointing at the demo data."""
    wd = tmp_path_factory.mktemp("workflow")
    configs = wd / "configs"
    configs.mkdir()

    aligners = yaml.safe_load((REPO_ROOT / "configs" / "aligners.yaml").read_text())
    # Single-threaded throughout: aligner output may depend on the thread count, and
    # the runtime-mode run below (always 1 thread) must reproduce this ensemble. It also
    # avoids the missing OpenMP support of Clustal Omega on osx-arm64.
    for aligner in aligners["aligners"].values():
        aligner["threads"] = 1

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
                    "seed": 13,
                }
            }
        )
    )
    return wd


@pytest.fixture(scope="module")
def workdir(tmp_path_factory, demo_input) -> Path:
    return _make_workdir(tmp_path_factory, demo_input)


@pytest.fixture(scope="module")
def runtime_workdir(tmp_path_factory, demo_input) -> Path:
    return _make_workdir(tmp_path_factory, demo_input)


def _snakemake(workflow: str, workdir: Path, *extra: str) -> None:
    snakemake = shutil.which("snakemake")
    assert snakemake is not None
    proc = subprocess.run(
        [snakemake, "-s", str(REPO_ROOT / "workflow" / workflow), "-d", str(workdir), "-c", "2", *extra],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stdout[-3000:] + proc.stderr[-3000:]


def _ensemble_contents(workdir: Path) -> dict[str, str]:
    root = workdir / "results"
    return {str(p.relative_to(root)): p.read_text() for p in sorted(root.glob("*/*/ensemble/*.fasta"))}


def test_ensemble_workflow(workdir):
    _snakemake("ensemble.smk", workdir)
    for dataset in ["protein", "dna"]:
        ensemble_dir = workdir / "results" / "demo" / dataset / "ensemble"
        files = sorted(p.name for p in ensemble_dir.iterdir())
        assert sum(f.startswith("MAFFT_FFT-NS-2.") for f in files) == 2
        assert sum(f.startswith("ClustalO.") for f in files) == 2
        assert sum(f.startswith("Muscle5.") for f in files) >= 1


def test_scoring_workflow(workdir):
    _snakemake("scoring.smk", workdir)
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



def test_runtime_mode_reproduces_ensemble_and_writes_benchmarks(workdir, runtime_workdir):
    # Second, independent run in runtime mode (single-threaded jobs) with the same seed.
    _snakemake("ensemble.smk", runtime_workdir, "--config", "runtime=True")
    _snakemake("scoring.smk", runtime_workdir, "--config", "runtime=True")

    # Same seed and config in a separate run and directory -> identical ensembles.
    assert _ensemble_contents(runtime_workdir) == _ensemble_contents(workdir)

    results = runtime_workdir / "results"
    tool_benchmarks = sorted(results.glob("demo/*/*/benchmark.json"))
    assert len(tool_benchmarks) == 2 * len(ENSEMBLE)
    assert len(list(results.glob("demo/*/scoring_benchmark.json"))) == 2
    # Normal mode writes no benchmark files.
    assert not list((workdir / "results").glob("demo/*/*/benchmark.json"))

    from ensemblify import benchmark

    benchmark.main(["--data-root", str(results)])
    df = pd.read_parquet(results / "ensemble_benchmarks.parquet")
    assert (df["status"] == "success").all() and (df["threads"] == 1).all()
