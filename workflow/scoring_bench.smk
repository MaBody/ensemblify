configfile: "configs/general.yaml"
configfile: "configs/aligners.yaml"
configfile: "configs/ensemble.yaml"

import os
import pathlib
import pickle
import pandas as pd
import numpy as np
from Bio import AlignIO
from ensemblify.scoring import descriptive_stats
from ensemblify.utils import infer_data_type

from aldiscore.scoring import pairwise
from aldiscore.scoring import set_based
from aldiscore.scoring import pythia
from aldiscore.datastructures.alignment import Alignment
from aldiscore.datastructures.ensemble import Ensemble


OUT_DIR  = pathlib.Path(config["general"]["output"])
RAXML_PATH  = pathlib.Path(config["general"]["raxml"])

# BENCHMARKS = ["b1", "b2"]
BENCHMARKS = ["bali2dna",  "bali2dnaf",  "bali3",  "ox",  "prefab4",  "sabre"]
DATASET_MAP = {benchmark:os.listdir(OUT_DIR / benchmark) for benchmark in BENCHMARKS}


wildcard_constraints:
    benchmark="|".join(BENCHMARKS),
    dataset="|".join([dataset for datasets in DATASET_MAP.values() for dataset in datasets])


rule all:
    input:
        done = expand(OUT_DIR / "{benchmark}" / "stats.parquet", benchmark=BENCHMARKS)


rule compute_scores:
    input:
        ensemble_dir = OUT_DIR / "{benchmark}" / "{dataset}" / "ensemble"
    output:
        stats = OUT_DIR / "{benchmark}" / "{dataset}" / "stats.parquet"
    run:
        ensemble_dir = Path(input.ensemble_dir)
        alignments = []
        for msa_file in os.listdir(ensemble_dir):
            alignment = Alignment(AlignIO.read(ensemble_dir / msa_file, "fasta"))
            alignments.append(alignment)

        ensemble = Ensemble(alignments)
        ref_dir = OUT_DIR / wildcards.dataset / "reference.fasta"
        
        has_reference = os.path.exists(ref_dir)
        if has_reference:
            reference = Alignment(AlignIO.read(ref_dir))
        
        source = wildcards.benchmark
        dataset = wildcards.dataset
        stats = {}
        measure = pairwise.DPosDistance(format="flat")
        scores_dpos = measure.compute(ensemble)
        stats[(source, dataset, "dpos")] = descriptive_stats(scores_dpos)

        try: # Fewer than 4 sequences cause error in Raxml
            scores_pythia = pythia.compute_pythia_difficulty(ensemble, RAXML_PATH)
        except:
            empty = {key:np.nan for key in stats[(source, dataset, "dpos")]}
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


rule collect:
    input:
        stats = lambda wildcards: expand(
            OUT_DIR / "{benchmark}" / "{dataset}" / "stats.parquet",
            benchmark=[wildcards.benchmark],
            dataset=DATASET_MAP[wildcards.benchmark]
        )
    output:
        stats_global = OUT_DIR / "{benchmark}" / "stats.parquet"
    run:
        stats_dfs = []
        for stats_file in input.stats:
            stats_df = pd.read_parquet(stats_file)
            stats_dfs.append(stats_df)

        stats_global_df = pd.concat(stats_dfs, axis=0)
        stats_global_df.to_parquet(output.stats_global)