import os
import stat

import pytest

from ensemblify import utils


def test_resolve_bare_name_uses_path(tmp_path, monkeypatch):
    exe = tmp_path / "fake-aligner"
    exe.write_text("#!/bin/sh\n")
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv("PATH", str(tmp_path) + os.pathsep + os.environ["PATH"])
    assert utils.resolve_executable("fake-aligner") == exe


def test_resolve_explicit_path(tmp_path):
    exe = tmp_path / "bin" / "muscle5.1.linux_intel64"
    exe.parent.mkdir()
    exe.write_text("#!/bin/sh\n")
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    assert utils.resolve_executable(str(exe)) == exe.resolve()


def test_resolve_relative_path(tmp_path, monkeypatch):
    exe = tmp_path / "tools" / "aligner"
    exe.parent.mkdir()
    exe.write_text("#!/bin/sh\n")
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    monkeypatch.chdir(tmp_path)
    assert utils.resolve_executable("tools/aligner") == exe.resolve()


@pytest.mark.parametrize("path", ["definitely-not-an-aligner", "/no/such/aligner", "", None])
def test_resolve_missing_raises(path):
    with pytest.raises(FileNotFoundError):
        utils.resolve_executable(path)


def test_resolve_non_executable_raises(tmp_path):
    f = tmp_path / "aligner"
    f.write_text("not executable")
    with pytest.raises(FileNotFoundError):
        utils.resolve_executable(str(f))


def test_check_aligners_reports_all_problems():
    config = {
        "aligners": {"A": {"path": "/no/such/a"}},
        "ensemble": {"a1": {"aligner": "A"}, "b1": {"aligner": "B"}},
    }
    with pytest.raises(FileNotFoundError) as e:
        utils.check_aligners(config)
    assert "a1" in str(e.value) and "b1" in str(e.value)


def test_option_variants_cartesian_product():
    params = {"--flag": None, "--fixed": 6, "--x": [1, 2], "--y": ["a", "b", "c"]}
    variants = utils.get_option_variants(params)
    assert len(variants) == 6
    assert all(v.startswith(utils._SPACER.join(["--flag", "--fixed", "6"])) for v in variants)


def test_format_cmd_resolves_aligner(tmp_path):
    exe = tmp_path / "aligner"
    exe.write_text("#!/bin/sh\n")
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    config = {
        "aligners": {"X": {"path": str(exe), "cmd": "{aligner} -i {input} -o {output} {options} -t {threads}"}},
        "ensemble": {"X1": {"aligner": "X"}},
    }
    options = utils.get_option_variants({"--a": 1})[0]
    cmd = utils.format_cmd(config, "X1", options, 4, tmp_path / "in.fa", tmp_path / "out.fa")
    assert cmd == [str(exe.resolve()), "-i", str(tmp_path / "in.fa"), "-o", str(tmp_path / "out.fa"), "--a", "1", "-t", "4"]
