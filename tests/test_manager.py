import os
import stat
import tempfile
from pathlib import Path

import pytest
from Bio import AlignIO

from ensemblify import utils
from ensemblify.manager import DefaultManager, Muscle5Manager

FASTA = ">a\nACGTACGT\n>b\nACGTTCGT\n>c\nACGAACGT\n"


def _script(path: Path, body: str) -> Path:
    path.write_text("#!/bin/sh\n" + body)
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


@pytest.fixture
def in_file(tmp_path) -> Path:
    f = tmp_path / "sequences.fasta"
    f.write_text(FASTA)
    return f


def test_run_cmd_shared_log_keeps_every_run(tmp_path):
    log = tmp_path / "run.log"
    cmd = ["sh", "-c", "echo out-$0; echo err-$0 >&2"]
    for i in range(3):
        utils.run_cmd(cmd[:2] + [cmd[2], str(i)], outfile=log, logfile=log)
    text = log.read_text()
    for i in range(3):
        assert f"out-{i}" in text and f"err-{i}" in text


def test_run_cmd_separate_out_and_log(tmp_path):
    out, log = tmp_path / "out.txt", tmp_path / "run.log"
    utils.run_cmd(["sh", "-c", "echo result; echo progress >&2"], outfile=out, logfile=log)
    assert out.read_text() == "result\n"
    assert log.read_text() == "progress\n"


def test_run_cmd_failure_message(tmp_path):
    log = tmp_path / "run.log"
    with pytest.raises(RuntimeError) as e:
        utils.run_cmd(["sh", "-c", "echo boom >&2; exit 3"], outfile=tmp_path / "o", logfile=log)
    msg = str(e.value)
    assert "exit code 3" in msg
    assert "boom" in msg
    assert str(log) in msg


def test_run_cmd_failure_without_logfile(tmp_path):
    with pytest.raises(RuntimeError, match="kaputt"):
        utils.run_cmd(["sh", "-c", "echo kaputt >&2; exit 1"], outfile=tmp_path / "o")


def _default_config(aligner: Path, params: dict) -> dict:
    return {
        "aligners": {"Fake": {"path": str(aligner), "cmd": "{aligner} {input} {output} {options}", "threads": 1}},
        "ensemble": {"FakeTool": {"aligner": "Fake", "type": "default", "params": params}},
    }


def test_default_manager_runs_all_variants_and_cleans_up(tmp_path, in_file):
    # Fake aligner: copies input to output, logs its options to stderr.
    aligner = _script(tmp_path / "fake", 'cp "$1" "$2"; shift 2; echo "run $*" >&2\n')
    out_dir = tmp_path / "FakeTool"
    out_dir.mkdir()
    log = out_dir / "ensemble.log"
    tmp_before = set(os.listdir(tempfile.gettempdir()))

    manager = DefaultManager(
        _default_config(aligner, {"--x": [1, 2, 3]}), "FakeTool", in_file, out_dir, log, 1, shuffle=True
    )
    ensemble = manager.compute()

    assert list(ensemble) == [("FakeTool", i) for i in range(3)]
    assert sorted(os.listdir(out_dir)) == ["ensemble.log", "msa.0.fasta", "msa.1.fasta", "msa.2.fasta"]
    assert set(os.listdir(tempfile.gettempdir())) - tmp_before == set()
    assert log.read_text().count("run --x") == 3


def test_save_ensemble_replaces_stale_members(tmp_path, in_file):
    ensemble_dir = tmp_path / "ensemble"
    ensemble_dir.mkdir()
    for i in range(5):
        (ensemble_dir / f"FakeTool.{i}.fasta").write_text(FASTA)
    (ensemble_dir / "FakeTool_long.0.fasta").write_text(FASTA)
    (ensemble_dir / "Other.0.fasta").write_text(FASTA)

    aligner = _script(tmp_path / "fake", 'cp "$1" "$2"\n')
    manager = DefaultManager(_default_config(aligner, {"--x": [1, 2]}), "FakeTool", in_file, tmp_path, tmp_path / "log", 1)
    manager.save_ensemble(manager.compute(), ensemble_dir)

    assert sorted(os.listdir(ensemble_dir)) == [
        "FakeTool.0.fasta",
        "FakeTool.1.fasta",
        "FakeTool_long.0.fasta",
        "Other.0.fasta",
    ]


def test_muscle5_manager_is_deterministic_and_ignores_stale_files(tmp_path, in_file):
    # Fake Muscle5: writes one file per permutation/replicate, replacing '@' in the template.
    aligner = _script(
        tmp_path / "fake-muscle5",
        'tpl="$2"\n'
        "for p in none abc acb bca; do for r in 0 1 2 3 4 5 6 7 8 9 10 11; do\n"
        '  printf ">%s\\nACGT\\n>b\\nACGT\\n>c\\nACGT\\n" "$p$r" > "${tpl%@*}$p.$r${tpl#*@}"\n'
        "done; done\n",
    )
    config = {
        "aligners": {"M5": {"path": str(aligner), "cmd": "{aligner} {input} {output} {options}", "threads": 1}},
        "ensemble": {"Muscle5": {"aligner": "M5", "type": "Muscle5", "params": {"-replicates": 12}}},
    }
    out_dir = tmp_path / "Muscle5"
    out_dir.mkdir()
    (out_dir / "msa.zzz.99.fasta").write_text(FASTA)  # left over from an earlier run

    ensemble = Muscle5Manager(config, "Muscle5", in_file, out_dir, out_dir / "ensemble.log", 1).compute()

    assert len(ensemble) == 48
    # The first record's id encodes which Muscle5 output file the alignment came from.
    sources = [ensemble[("Muscle5", i)][0].id for i in range(len(ensemble))]
    expected = [f"{p}{r}" for p in ["abc", "acb", "bca", "none"] for r in range(12)]
    assert sources == expected


def test_muscle5_manager_fails_without_output(tmp_path, in_file):
    aligner = _script(tmp_path / "fake-muscle5", "true\n")
    config = {
        "aligners": {"M5": {"path": str(aligner), "cmd": "{aligner} {input} {output}", "threads": 1}},
        "ensemble": {"Muscle5": {"aligner": "M5", "type": "Muscle5", "params": {"-replicates": 1}}},
    }
    with pytest.raises(RuntimeError, match="no alignments"):
        Muscle5Manager(config, "Muscle5", in_file, tmp_path, tmp_path / "log", 1).compute()
