configfile: "configs/general.yaml"
configfile: "configs/aligners.yaml"
configfile: "configs/ensemble.yaml"

import json
import os
import pathlib
import resource
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import AlignIO
from aldiscore.datastructures.alignment import Alignment
from aldiscore.datastructures.ensemble import Ensemble
from aldiscore.scoring import pairwise
from aldiscore.scoring import pythia
from ensemblify.scoring import descriptive_stats
from ensemblify.utils import infer_data_type


OUT_DIR = pathlib.Path(config["general"]["output"])
RAXML_PATH = pathlib.Path(config["general"]["raxml"])

SOURCES = config["general"]["sources"]
DATASET_MAP = {
    source: list(
        filter(
            lambda name: os.path.isdir(OUT_DIR / source / name),
            os.listdir(OUT_DIR / source),
        )
    )
    for source in SOURCES
}


wildcard_constraints:
    source="|".join(SOURCES),
    dataset="|".join([dataset for datasets in DATASET_MAP.values() for dataset in datasets])


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
        source_stats = expand(OUT_DIR / "{source}" / "stats.parquet", source=SOURCES)


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
            ensemble_dir = Path(input.ensemble_dir)
            alignments = []
            for msa_file in os.listdir(ensemble_dir):
                alignment = Alignment(AlignIO.read(ensemble_dir / msa_file, "fasta"))
                alignments.append(alignment)

            ensemble = Ensemble(alignments)
            ref_dir = ensemble_dir.parent / "reference.fasta"

            has_reference = os.path.exists(ref_dir)
            if has_reference:
                reference = Alignment(AlignIO.read(ref_dir, "fasta"))

            source = wildcards.source
            dataset = wildcards.dataset
            stats = {}
            measure = pairwise.DPosDistance(format="flat")
            scores_dpos = measure.compute(ensemble)
            stats[(source, dataset, "dpos")] = descriptive_stats(scores_dpos)

            try:  # Fewer than 4 sequences cause error in RAxML.
                scores_pythia = pythia.compute_pythia_difficulty(ensemble, RAXML_PATH)
            except Exception:
                empty = {key: np.nan for key in stats[(source, dataset, "dpos")]}
                stats[(source, dataset, "pythia")] = empty
            else:
                stats[(source, dataset, "pythia")] = descriptive_stats(scores_pythia)

            if has_reference:
                scores_dpos_ref = measure.compute(ensemble, reference)
                stats[(source, dataset, "dpos_ref")] = descriptive_stats(scores_dpos_ref)

            stats_df = pd.DataFrame(stats.values(), index=stats.keys())
            stats_df.index.names = ["source", "dataset", "method"]

            stats_df["datatype"] = infer_data_type(ensemble.dataset.records)
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
        stats_dfs = []
        for stats_file in input.stats:
            stats_df = pd.read_parquet(stats_file)
            stats_dfs.append(stats_df)

        stats_global_df = pd.concat(stats_dfs, axis=0)
        stats_global_df.to_parquet(output.stats_global)
