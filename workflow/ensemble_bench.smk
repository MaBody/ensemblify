configfile: "configs/general.yaml"
configfile: "configs/aligners.yaml"
configfile: "configs/ensemble.yaml"

import os
import pathlib
import pickle
import pandas as pd
import numpy as np
from Bio import SeqIO, AlignIO
from ensemblify.manager import infer_manager

INPUT_DIR = pathlib.Path(config["general"]["input"])
OUT_DIR  = pathlib.Path(config["general"]["output"])

TOOLS = list(config["ensemble"].keys())
THREADS = {key:val["threads"] for key,val in config["aligners"].items()}

SOURCES = config["general"]["sources"]
DATASET_MAP = {source:list(filter(lambda name: os.path.isdir(INPUT_DIR / source / name), os.listdir(INPUT_DIR / source))) for source in SOURCES}
# print("sources:", SOURCES)
# print("#datasets:", sum(len(names) for names in DATASET_MAP.values()))

rule all:
    input:
        done = expand(OUT_DIR / "{source}" / "{dataset}" / "{tool}" / "done", 
                zip,
                source=[
                    b for b in DATASET_MAP
                        for d in DATASET_MAP[b]
                        for t in TOOLS
                ],
                dataset=[
                    d for b in DATASET_MAP
                        for d in DATASET_MAP[b]
                        for t in TOOLS
                ],
                tool=[
                    t for b in DATASET_MAP
                        for d in DATASET_MAP[b]
                        for t in TOOLS
                ]
            )

# # Convert to expected format
# rule prepare_sequences:
#     input:
#         msa = INPUT_DIR / "{source}" / "in" / "{dataset}"
#     output:
#         sequences = OUT_DIR / "{source}" / "{dataset}" / "sequences.fasta"
#     params:
#         source = lambda wildcards: wildcards.source,
#         dataset = lambda wildcards: wildcards.dataset
#     run:
#         seqs = list(SeqIO.parse(input.msa, "fasta"))
#         SeqIO.write(seqs, output.sequences, "fasta")

    

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
        # in_file = rules.prepare_sequences.output.sequences
        in_file = INPUT_DIR / "{source}" / "{dataset}" / "sequences.fasta"
    output:
        done = OUT_DIR / "{source}" / "{dataset}" / "{tool}" / "done"
    log:
        log_file = OUT_DIR / "{source}" / "{dataset}" / "{tool}"  / "ensemble.log"
    params:
        _tool = lambda wildcards: wildcards.tool,
        _dataset = lambda wildcards: wildcards.dataset
    run:
        in_file = input.in_file
        out_dir = Path(output.done).parent
        log_file = pathlib.Path(log.log_file)
        manager_type = config["ensemble"][params._tool]["type"]
        manager_class = infer_manager(manager_type)
        manager = manager_class(config, params._tool, in_file, out_dir, log_file, threads, shuffle=True)

        ensemble = manager.compute()
        ensemble_dir = out_dir.parent / "ensemble"
        os.makedirs(ensemble_dir, exist_ok=True)
        manager.save_ensemble(ensemble, ensemble_dir)
        
        # Set output flag
        open(output.done, "a").close()

