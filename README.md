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

### c) Optional dependencies

Some utilities (e.g., improved CLI output) may require additional packages:

```bash
pip install rich_argparse
```

---

### Notes

- Make sure all external alignment tools are installed separately and configured in `aligners.yaml`
- Snakemake is required to run workflows
- Python ≥ 3.9 recommended
```

---

## Recommended usage

### 1. Generate an MSA ensemble

Run the ensemble generation workflow:

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

---

### 2. Compute MSA difficulty scores

Run the scoring workflow:

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

---

## Input structure

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

---

## Output structure

After running the ensemble workflow, the output directory will look similar to:

```text
dataset/
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

### Description

- `ClustalO/`, `MAFFT_*`, `Muscle*`
  - Individual aligner outputs and replicates.

- `ensemble/`
  - Aggregated alignments from all methods.

- `sequences.fasta`
  - Original input sequences.

- `stats.parquet`
  - Computed statistics and difficulty metrics.

---

# Configuration

## `general.yaml`

Specify the location of your datasets and working directories.

The workflow expects sequence files at:

```text
work_dir/database/dataset/sequences.fasta
```

---

## `ensemble.yaml`

Specify:

- Alignment methods to include
- Number of replicates
- Parameters for each aligner
- Ensemble generation settings

Example use cases:

- Enable/disable specific aligners
- Modify MAFFT strategy
- Adjust bootstrap or replicate settings

---

## `aligners.yaml`

Provide the paths to installed alignment software.

Example:

```yaml
ClustalO: /path/to/clustalo
MAFFT: /path/to/mafft
Muscle3: /path/to/muscle3
Muscle5: /path/to/muscle5
```

---

## Requirements

- Snakemake
- Clustal Omega
- MAFFT
- MUSCLE3
- MUSCLE5

Make sure all aligners are installed and available in `aligners.yaml`.

---

## Workflow summary

```text
sequences.fasta
        ↓
Multiple aligners + parameterizations
        ↓
MSA ensemble generation
        ↓
Combined ensemble
        ↓
Pairwise d_pos computation
        ↓
Difficulty statistics
```