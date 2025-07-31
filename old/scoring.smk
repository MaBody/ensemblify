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


INPUT_DIR = pathlib.Path(config["general"]["input"])
OUT_DIR  = pathlib.Path(config["general"]["output"])
RAXML_PATH  = pathlib.Path(config["general"]["raxml"])

TOOLS = list(config["ensemble"].keys())

DATASETS = list(filter(lambda name: os.path.isdir(INPUT_DIR / name), os.listdir(INPUT_DIR)))


rule all:
    input:
        done = OUT_DIR / "stats.parquet"


rule compute_scores:
    input:
        ensemble_dir = OUT_DIR / "{dataset}" / "ensemble"
    output:
        stats = OUT_DIR / "{dataset}" / "stats.parquet"
    params:
        _dataset = lambda wildcards: wildcards.dataset
    run:
        ensemble_dir = Path(input.ensemble_dir)
        alignments = []
        for msa_file in os.listdir(ensemble_dir):
            alignment = Alignment(AlignIO.read(ensemble_dir / msa_file, "fasta"))
            alignments.append(alignment)

        ensemble = Ensemble(alignments)
        ref_dir = ensemble_dir.parent / "reference.fasta"
        
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
        stats = expand(rules.compute_scores.output.stats, dataset=DATASETS)
    output:
        stats_global = OUT_DIR / "stats.parquet"
    run:
        stats_dfs = []
        for stats_file in input.stats:
            stats_df = pd.read_parquet(stats_file)
            stats_dfs.append(stats_df)

        stats_global_df = pd.concat(stats_dfs, axis=0)
        stats_global_df.to_parquet(output.stats_global)