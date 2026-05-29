# ensemblify

Generate Multiple Sequence Alignment (MSA) ensembles and compute alignment difficulty metrics.

## Overview

`ensemblify` generates ensembles of MSAs using multiple alignment methods and computes dispersion-based scores to estimate alignment difficulty.

The workflow consists of two main steps:

1. **Generate an MSA ensemble**
2. **Compute pairwise positional dispersion (`d_pos`) across the ensemble**

---

## Setup

Clone this repository and navigate to the project root:

```bash
git clone https://github.com/MaBody/ensemblify.git
cd ensemblify
```

---

### a) Library + CLI (pip / setuptools build)

You can install `ensemblify` as a Python package using a standard build:

```bash
python -m build
```

This will generate a `dist/` folder containing the wheel file.

Install it with:

```bash
pip install dist/ensemblify-<version>.whl
```

Alternatively, install directly in editable mode (recommended for development):

```bash
pip install -e .
```

---

### b) Library on Conda (recommended)

If you prefer a reproducible environment, use the provided `environment.yml` file:

```bash
conda env create -f environment.yml
conda activate ensemblify
```

Then install the package itself:

```bash
pip install .
```


---

## Recommended usage

### 0. Edit Configuration files

#### `general.yaml`

Specify the location of your datasets and working directories.

The workflow expects input sequences organized as:

```text
work_dir/
└── database/
    └── dataset/
        └── sequences.fasta
```

Example:

```text
data/
└── ox/
    └── 104s10/
        └── sequences.fasta
```

so the `general.yaml` file would look like: 
```yaml
general:

  input: path/to/data/
  output: path/to/data/
  sources: ["ox"]

  raxml: "/path/to/raxml-ng"
```
---

#### `aligners.yaml`

Provide the paths to installed alignment software.

Example `aligners.yaml` for `ClustalO`:

```yaml
aligners:
    ClustalO:
        path: /path/to/miniconda3/envs/ensemblify/bin/clustalo
        cmd: "{aligner} -i {input} -o {output} -v --force {options} --threads={threads}"
        threads: 8 
```
---

#### `ensemble.yaml`

You can modify the current settings:

- Alignment methods to include
- Number of threads per aligner
- Parameters for each aligner

Example of `ClustalO` parameters in `ensemble.yaml`: 

```yaml
ensemble:
  ClustalO:
    aligner: ClustalO
    type: "default"
    params:
      "--iterations": 6
      "--max-guidetree-iterations": [1, 2, 3, 4]
      "--max-hmm-iterations": [1, 2]
```

---

### 1. Generate an MSA ensemble

Run the ensemble generation workflow with snakemake:

```bash
snakemake -s workflow/ensemble_bench.smk \
    -c 40 \
    --verbose \
    --latency-wait 60 \
    --printshellcmds \
    --keep-going \
    > stdout_ensemble.log 2>&1 & disown
```

This step generates multiple alignments from different aligners and parameterizations.

After running the ensemble workflow, the output directory will look similar to:

```text
work_dir/
└── database/
    └── dataset/
        ├── ClustalO/
        │   ├── done
        │   ├── ensemble.log
        │   └── msa.[0-7].fasta
        │
        ├── MAFFT_FFT-NS-2/
        │   └── msa.[0-7].fasta
        │
        ├── MAFFT_G-INS-i/
        │   └── msa.[0-7].fasta
        │
        ├── MAFFT_L-INS-i/
        │   └── msa.[0-7].fasta
        │
        ├── Muscle3/
        │   └── msa.[0-7].fasta
        │
        ├── Muscle5/
        │   └── msa.[0-7].fasta
        │
        ├── ensemble/
        │   ├── ClustalO.[0-7].fasta
        │   ├── MAFFT_FFT-NS-2.[0-7].fasta
        │   ├── MAFFT_G-INS-i.[0-7].fasta
        │   ├── MAFFT_L-INS-i.[0-7].fasta
        │   ├── Muscle3.[0-7].fasta
        │   └── Muscle5.[0-7].fasta
        │
        ├── sequences.fasta
        └── stats.parquet
```

---

### 2. Compute MSA difficulty scores

Run the scoring workflow with snakemake:

```bash
snakemake -s workflow/scoring_bench.smk \
    -c 20 \
    --verbose \
    --latency-wait 60 \
    --printshellcmds \
    --keep-going \
    > stdout_scoring.log 2>&1 & disown
```

This computes pairwise positional dispersion (`d_pos`) across alignments in the ensemble to estimate alignment difficulty. 
The calculated alignment difficulty per database is stored in the `stats.parquet` file

```
work_dir/database/dataset/stats.parquet
```

To convert the .parquet file into a .csv file you can

```
# Install parquet-tools
pip install parquet-tools

# convert to csv
parquet-tools csv work_dir/database/dataset/stats.parquet > work_dir/database/dataset/stats.csv
```


---

### Notes

- Make sure all external alignment tools are installed separately and configured in `aligners.yaml`
- Snakemake is required to run workflows
- Python ≥ 3.9 recommended

---
