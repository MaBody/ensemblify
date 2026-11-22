configfile: "configs/general.yaml"
configfile: "configs/aligners.yaml"
configfile: "configs/ensemble.yaml"

# Score every complete ensemble produced by workflow/ensemble.smk.
#
#   Output: {output}/{source}/{dataset}/stats.parquet  (per dataset)
#           {output}/{source}/stats.parquet            (per source)
#
# Runtime mode (`--config runtime=True`, or `pixi run score-runtime`) runs every
# job single-threaded and writes a scoring_benchmark.json per dataset.

import pathlib

import pandas as pd

from ensemblify import discovery, scoring
from ensemblify.benchmark import record_runtime
from ensemblify.utils import resolve_executable


INPUT_DIR = pathlib.Path(config["general"]["input"])
OUT_DIR = pathlib.Path(config["general"]["output"])
RAXML_PATH = resolve_executable(config["general"]["raxml"])
RUNTIME = bool(config.get("runtime", False))

TOOLS = list(config["ensemble"].keys())
SOURCES = config["general"]["sources"]
DATASET_MAP = discovery.scorable_datasets(OUT_DIR, SOURCES, TOOLS)
ACTIVE_SOURCES = [source for source in SOURCES if DATASET_MAP[source]]

wildcard_constraints:
    source=discovery.wildcard_regex(SOURCES),
    dataset=discovery.wildcard_regex(d for datasets in DATASET_MAP.values() for d in datasets)


DATASET_DIR = OUT_DIR / "{source}" / "{dataset}"
OUTPUTS = {"stats": DATASET_DIR / "stats.parquet"}
if RUNTIME:
    OUTPUTS["benchmark"] = DATASET_DIR / "scoring_benchmark.json"


rule all:
    input:
        source_stats = expand(OUT_DIR / "{source}" / "stats.parquet", source=ACTIVE_SOURCES),
        benchmarks = [
            str(OUTPUTS["benchmark"]).format(source=source, dataset=dataset)
            for source in ACTIVE_SOURCES
            for dataset in DATASET_MAP[source]
        ] if RUNTIME else [],


rule compute_scores:
    input:
        ensemble_dir = DATASET_DIR / "ensemble"
    output:
        **OUTPUTS
    threads: 1 if RUNTIME else 4
    run:
        with record_runtime(
            output.get("benchmark"),
            benchmark_scope="scoring_dataset",
            source=wildcards.source,
            dataset=wildcards.dataset,
            threads=threads,
        ):
            ensemble = scoring.load_ensemble(input.ensemble_dir)
            reference = scoring.find_reference(
                INPUT_DIR / wildcards.source / wildcards.dataset,
                OUT_DIR / wildcards.source / wildcards.dataset,
            )
            stats_df = scoring.score_ensemble(ensemble, reference, RAXML_PATH, threads=threads)
            stats_df = pd.concat({(wildcards.source, wildcards.dataset): stats_df}, names=["source", "dataset"])
            stats_df.to_parquet(output.stats)


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
