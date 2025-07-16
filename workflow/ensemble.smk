configfile: "configs/general.yaml"
configfile: "configs/aligners.yaml"
configfile: "configs/ensemble.yaml"

import os
import pathlib
import pickle
import pandas as pd
import numpy as np
from Bio import AlignIO
from ensemblify.manager import infer_manager

INPUT_DIR = pathlib.Path(config["general"]["input"])
OUT_DIR  = pathlib.Path(config["general"]["output"])

TOOLS = list(config["ensemble"].keys())
THREADS = {key:val["threads"] for key,val in config["aligners"].items()}

DATASETS = list(filter(lambda name: os.path.isdir(INPUT_DIR / name), os.listdir(INPUT_DIR)))


rule all:
    input:
        done = expand(OUT_DIR / "{dataset}" / "{tool}" / "done", dataset=DATASETS, tool=TOOLS)


# Ensemble generation assumes the following file structure:
#   input_dir/ 
#       dataset_1/ 
#           sequences.fasta
#       dataset_2/
#           sequences.fasta
#       ...      
rule generate_ensemble:
    threads: lambda wildcards: config["aligners"][config["ensemble"][wildcards.tool]["aligner"]]["threads"]
    input:
        in_file = OUT_DIR / "{dataset}" / "sequences.fasta"
    output:
        done = OUT_DIR / "{dataset}" / "{tool}" / "done"
    log:
        log_file = OUT_DIR / "{dataset}" / "{tool}"  / "ensemble.log"
    params:
        _tool = lambda wildcards: wildcards.tool,
        _dataset = lambda wildcards: wildcards.dataset
    run:
        in_file = input.in_file
        out_dir = OUT_DIR / params._dataset / params._tool
        log_file = pathlib.Path(log.log_file)
        manager_class = infer_manager(config["ensemble"][params._tool]["aligner"])
        manager = manager_class(config, params._tool, in_file, out_dir, log_file, threads, shuffle=True)
 
        ensemble = manager.compute()
        ensemble_dir = OUT_DIR / params._dataset / "ensemble"
        os.makedirs(ensemble_dir, exist_ok=True)
        manager.save_ensemble(ensemble, ensemble_dir)
        
        # Set output flag
        open(output.done, "a").close()


# rule collect:
#     input:
#         done = expand(rules.generate_ensemble.output.done, dataset=DATASETS, tool=TOOLS)
#     output:
#         global_done = OUT_DIR / "done"
#     run:
#         # Set output flag
#         open(output.global_done, "a").close()