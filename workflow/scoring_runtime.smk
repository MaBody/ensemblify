configfile: "configs/general.yaml"
configfile: "configs/aligners.yaml"
configfile: "configs/ensemble.yaml"

import json
import pathlib
import resource
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from ensemblify import discovery, scoring
from ensemblify.utils import resolve_executable


INPUT_DIR = pathlib.Path(config["general"]["input"])
OUT_DIR = pathlib.Path(config["general"]["output"])
RAXML_PATH = resolve_executable(config["general"]["raxml"])

TOOLS = list(config["ensemble"].keys())
SOURCES = config["general"]["sources"]
DATASET_MAP = discovery.scorable_datasets(OUT_DIR, SOURCES, TOOLS)
ACTIVE_SOURCES = [source for source in SOURCES if DATASET_MAP[source]]


wildcard_constraints:
    source=discovery.wildcard_regex(SOURCES),
    dataset=discovery.wildcard_regex(d for datasets in DATASET_MAP.values() for d in datasets)


def _expanded_dataset_paths(file_name):
    return expand(
        OUT_DIR / "{source}" / "{dataset}" / file_name,
        zip,
        source=[
            source
            for source in DATASET_MAP
            for dataset in DATASET_MAP[source]
        ],
        dataset=[
            dataset
            for source in DATASET_MAP
            for dataset in DATASET_MAP[source]
        ],
    )


def _cpu_snapshot():
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    return {
        "process": time.process_time(),
        "children": children.ru_utime + children.ru_stime,
    }


def _cpu_elapsed(start, end):
    return {
        "process_seconds": end["process"] - start["process"],
        "children_seconds": end["children"] - start["children"],
        "total_seconds": (
            end["process"]
            - start["process"]
            + end["children"]
            - start["children"]
        ),
    }


def _write_benchmark(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


rule all:
    input:
        scoring_benchmark = _expanded_dataset_paths("scoring_benchmark.json"),
        scoring_stats = _expanded_dataset_paths("stats.parquet"),
        source_stats = expand(OUT_DIR / "{source}" / "stats.parquet", source=ACTIVE_SOURCES)


rule compute_scores_runtime:
    input:
        ensemble_dir = OUT_DIR / "{source}" / "{dataset}" / "ensemble"
    output:
        stats = OUT_DIR / "{source}" / "{dataset}" / "stats.parquet",
        benchmark = OUT_DIR / "{source}" / "{dataset}" / "scoring_benchmark.json"
    run:
        started_at = datetime.now(timezone.utc)
        wall_start = time.perf_counter()
        cpu_start = _cpu_snapshot()
        status = "success"
        exception_type = None
        exception_message = None

        try:
            ensemble = scoring.load_ensemble(input.ensemble_dir)
            reference = scoring.find_reference(
                INPUT_DIR / wildcards.source / wildcards.dataset,
                OUT_DIR / wildcards.source / wildcards.dataset,
            )
            stats_df = scoring.score_ensemble(ensemble, reference, RAXML_PATH, threads=threads)
            stats_df = pd.concat({(wildcards.source, wildcards.dataset): stats_df}, names=["source", "dataset"])
            stats_df.to_parquet(output.stats)
        except Exception as exc:
            status = "failed"
            exception_type = type(exc).__name__
            exception_message = str(exc)
            raise
        finally:
            ended_at = datetime.now(timezone.utc)
            wall_seconds = time.perf_counter() - wall_start
            cpu_seconds = _cpu_elapsed(cpu_start, _cpu_snapshot())
            _write_benchmark(
                output.benchmark,
                {
                    "benchmark_scope": "scoring_dataset",
                    "source": wildcards.source,
                    "dataset": wildcards.dataset,
                    "status": status,
                    "started_at": started_at.isoformat(),
                    "ended_at": ended_at.isoformat(),
                    "wall_seconds": wall_seconds,
                    "cpu_seconds": cpu_seconds,
                    "exception_type": exception_type,
                    "exception_message": exception_message,
                },
            )


rule collect:
    input:
        stats = lambda wildcards: expand(
            OUT_DIR / "{source}" / "{dataset}" / "stats.parquet",
            source=[wildcards.source],
            dataset=DATASET_MAP[wildcards.source],
        )
    output:
        stats_global = OUT_DIR / "{source}" / "stats.parquet"
    run:
        stats_global_df = pd.concat([pd.read_parquet(f) for f in input.stats], axis=0)
        stats_global_df.to_parquet(output.stats_global)
