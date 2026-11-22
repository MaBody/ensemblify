"""Checks that the environment from `pixi install` is complete."""

import importlib.metadata
import re
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

import ensemblify

REPO_ROOT = Path(__file__).resolve().parent.parent
PYPROJECT = REPO_ROOT / "pyproject.toml"


def test_version_matches_pyproject():
    match = re.search(r'^version\s*=\s*"([^"]+)"', PYPROJECT.read_text(), re.MULTILINE)
    assert match is not None
    assert importlib.metadata.version("ensemblify") == match.group(1)
    assert ensemblify.__version__ == match.group(1)


def test_aldiscore_importable():
    from aldiscore.datastructures.ensemble import Ensemble  # noqa: F401
    from aldiscore.scoring import pairwise, pythia  # noqa: F401


def test_pythia_importable():
    import pypythia.predictor  # noqa: F401


@pytest.mark.parametrize("name", ["general", "aligners", "ensemble"])
def test_configs_parse(name):
    config = yaml.safe_load((REPO_ROOT / "configs" / f"{name}.yaml").read_text())
    assert name in config


@pytest.mark.parametrize("binary", ["snakemake", "mafft", "clustalo", "muscle", "raxml-ng"])
def test_binaries_on_path(binary):
    assert shutil.which(binary) is not None, f"{binary} missing; run inside `pixi run`"


def test_muscle_on_path_is_muscle5():
    out = subprocess.run(["muscle", "-version"], capture_output=True, text=True)
    assert re.search(r"muscle\s+5", out.stdout + out.stderr, re.IGNORECASE)
