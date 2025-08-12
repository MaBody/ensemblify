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

OUT_DIR  = pathlib.Path("/hits/fast/cme/bodynems/data/paper/treebase_v1")

TOOLS = list(config["ensemble"].keys())

DATASETS = list(filter(lambda name: os.path.isdir(OUT_DIR / name), os.listdir(OUT_DIR)))


rule all:
    input:
        done = expand(OUT_DIR / "{dataset}" / "sequences.fasta", dataset=DATASETS)


rule copy_sequences:
    input:
        in_file = OUT_DIR / "{dataset}" / "ensemble" / "Muscle5.0.fasta"
    output:
        out_file = OUT_DIR / "{dataset}" / "sequences.fasta"
    run:
        msa = AlignIO.read(input.in_file, "fasta")
        with open(output.out_file, "w") as outfile:
            for sequence in msa:
                outfile.write(f">{sequence.id}\n")
                outfile.write(str(sequence.seq).replace("-", ""))
                outfile.write("\n")

