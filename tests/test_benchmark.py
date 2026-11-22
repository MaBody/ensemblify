import json

import pandas as pd
import pytest

from ensemblify import benchmark


def test_record_runtime_success(tmp_path):
    path = tmp_path / "a" / "benchmark.json"
    with benchmark.record_runtime(path, source="s", tool="t"):
        sum(range(1000))
    record = json.loads(path.read_text())
    assert record["status"] == "success" and record["source"] == "s"
    assert record["wall_seconds"] >= 0
    assert set(record["cpu_seconds"]) == {"process_seconds", "children_seconds", "total_seconds"}


def test_record_runtime_failure_is_recorded_and_reraised(tmp_path):
    path = tmp_path / "benchmark.json"
    with pytest.raises(ValueError):
        with benchmark.record_runtime(path):
            raise ValueError("nope")
    record = json.loads(path.read_text())
    assert record["status"] == "failed"
    assert record["exception_type"] == "ValueError" and record["exception_message"] == "nope"


def test_record_runtime_without_path_is_noop(tmp_path):
    with benchmark.record_runtime(None, source="s"):
        pass
    assert list(tmp_path.iterdir()) == []


def test_compile_cli(tmp_path):
    for tool in ["T1", "T2"]:
        with benchmark.record_runtime(tmp_path / "src" / "d1" / tool / "benchmark.json", tool=tool):
            pass
    with benchmark.record_runtime(tmp_path / "src" / "d1" / "scoring_benchmark.json"):
        pass

    benchmark.main(["--data-root", str(tmp_path)])

    ensemble = pd.read_parquet(tmp_path / "ensemble_benchmarks.parquet")
    scoring = pd.read_parquet(tmp_path / "scoring_benchmarks.parquet")
    assert sorted(ensemble["tool"]) == ["T1", "T2"]
    assert (ensemble["source"] == "src").all() and (ensemble["dataset"] == "d1").all()
    assert len(scoring) == 1 and scoring.loc[0, "dataset"] == "d1"
    assert "cpu_total_seconds" in scoring


def test_compile_empty_root_writes_schema(tmp_path):
    df, path = benchmark.compile_benchmarks(tmp_path, "scoring")
    assert df.empty and path.is_file()
    assert {"source", "dataset", "wall_seconds"} <= set(pd.read_parquet(path).columns)


def test_compile_output_requires_single_kind(tmp_path):
    with pytest.raises(SystemExit):
        benchmark.main(["--data-root", str(tmp_path), "--output", str(tmp_path / "x.parquet")])
