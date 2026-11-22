import re

from ensemblify import discovery


def _touch(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("")


def test_input_datasets_require_sequences(tmp_path):
    _touch(tmp_path / "src" / "b" / "sequences.fasta")
    _touch(tmp_path / "src" / "a" / "sequences.fasta")
    (tmp_path / "src" / "empty").mkdir()
    _touch(tmp_path / "src" / "stray.txt")

    assert discovery.input_datasets(tmp_path, ["src"]) == {"src": ["a", "b"]}


def test_missing_source_dir_is_reported(tmp_path, capsys):
    assert discovery.input_datasets(tmp_path, ["nope"]) == {"nope": []}
    assert "nope" in capsys.readouterr().err


def test_scorable_datasets_need_all_tools_done(tmp_path, capsys):
    for tool in ["T1", "T2"]:
        _touch(tmp_path / "src" / "complete" / tool / "done")
    _touch(tmp_path / "src" / "complete" / "ensemble" / "T1.0.fasta")
    _touch(tmp_path / "src" / "partial" / "T1" / "done")
    _touch(tmp_path / "src" / "partial" / "ensemble" / "T1.0.fasta")

    assert discovery.scorable_datasets(tmp_path, ["src"], ["T1", "T2"]) == {"src": ["complete"]}
    err = capsys.readouterr().err
    assert "partial" in err and "T2" in err


def test_wildcard_regex_escapes_names():
    pattern = discovery.wildcard_regex(["a.b", "c+d"])
    assert re.fullmatch(pattern, "a.b") and re.fullmatch(pattern, "c+d")
    assert not re.fullmatch(pattern, "aXb")
