# ensemblify

Generate ensembles of Multiple Sequence Alignments (MSAs) and compute dispersion-based alignment difficulty scores.

`ensemblify` runs several aligners (MAFFT, Clustal Omega, Muscle3, Muscle5) under a range of parameter settings, collects the resulting alignments into an ensemble, and scores the ensemble's dispersion with [AlDiScore](https://github.com/MaBody/aldiscore) (`d_pos`) and [Pythia](https://github.com/tschuelia/PyPythia).

The pipeline has two [Snakemake](https://snakemake.readthedocs.io) workflows:

1. **Ensemble generation**: align every dataset with every configured aligner variant.
2. **Scoring**: compute the pairwise positional dispersion (`d_pos`) and Pythia difficulty across each ensemble.

## Setup

Clone this repository and navigate to the top folder:

```shell
git clone https://github.com/MaBody/ensemblify.git
cd ensemblify
```

Install with [pixi](https://pixi.sh):

```shell
pixi install
```

This creates a conda environment with Python, Snakemake, the `ensemblify` package (editable) and AlDiScore. It also installs these prebuilt tools from conda-forge and bioconda:

| Tool          | Package    | Binary     |
| ------------- | ---------- | ---------- |
| MAFFT         | `mafft`    | `mafft`    |
| Clustal Omega | `clustalo` | `clustalo` |
| Muscle5       | `muscle>=5`| `muscle`   |
| RAxML-NG      | `raxml-ng` | `raxml-ng` |

### Muscle3

Muscle 3.8.31 and Muscle5 are the same conda package, so they cannot share one environment. Muscle3 gets its own:

```shell
pixi install -e muscle3
```

`configs/aligners.yaml` already points at `.pixi/envs/muscle3/bin/muscle`.

> **Apple Silicon (osx-arm64):**
> - bioconda has no native Muscle 3.8.31 build. pixi installs the osx-64 build instead, which runs under [Rosetta 2](https://support.apple.com/en-us/102527). Without Rosetta, remove the `Muscle3` entry from `configs/ensemble.yaml`.
> - The bioconda build of Clustal Omega has no OpenMP support. Set `threads: 1` for `ClustalO` in `configs/aligners.yaml`, or it fails with `Cannot change number of threads`.

### Manually installed aligners

Every tool path in the configs accepts either a bare executable name (looked up on `PATH`) or an explicit file path:

```yaml
aligners:
    Muscle5:
        path: /opt/tools/muscle5.1.linux_intel64   # instead of `muscle`
```

The ensemble workflow checks all configured aligners before it starts and reports any that are missing.

### Tests

```shell
pixi run test
```

The test suite checks the installation (package, AlDiScore, Pythia, all binaries) and runs both workflows on the bundled demo data.

## Quickstart

The default configuration runs on two small datasets in `demo/data/demo/` and writes to `results/`. `protein` holds 14 related protein sequences (from the AlDiScore demo), and `dna` holds 12 synthetic DNA sequences with substitutions and indels.

```shell
pixi install -e muscle3        # Muscle3 environment, see above
pixi run ensemble -c 8         # step 1: generate the MSA ensembles
pixi run score -c 8            # step 2: compute difficulty scores
```

`pixi run <task>` passes extra arguments on to Snakemake (`-c 8` sets the number of cores). For long runs, detach and log to a file:

```shell
nohup pixi run ensemble -c 40 > stdout_ensemble.log 2>&1 &
```

## Configuration

All configuration lives in `configs/`.

### `general.yaml`

Input and output locations. Input sequences must be organized as `{input}/{source}/{dataset}/sequences.fasta`:

```text
data/
└── ox/                  # source
    └── 104s10/          # dataset
        └── sequences.fasta
```

```yaml
general:
  input: path/to/data      # relative paths resolve against the repository root
  output: path/to/results  # may be the same as `input`
  sources: ["ox"]
  raxml: raxml-ng          # or an explicit path
  seed: 0                  # fixed input shuffle order; null = random
```

Before each aligner run, the input sequences are shuffled so that the ensemble does not inherit a bias from the order in `sequences.fasta`. The `seed` fixes this shuffle order. It does not make the ensemble reproducible: multithreaded aligners can return different alignments for the same input, e.g. MAFFT's iterative refinement (L-INS-i, G-INS-i).

### `aligners.yaml`

How each aligner is called: executable `path`, command template `cmd`, and `threads`:

```yaml
aligners:
    ClustalO:
        path: clustalo
        cmd: "{aligner} -i {input} -o {output} -v --force {options} --threads={threads}"
        threads: 8
```

Placeholders: `{aligner}`, `{input}`, `{output}`, `{options}`, `{threads}`. Templates without `{output}` (MAFFT) write the alignment to stdout.

### `ensemble.yaml`

Which aligner configurations make up the ensemble. A list value creates one run per value, and the runs span the Cartesian product of all list-valued parameters. `~` is a flag without a value.

```yaml
ensemble:
  ClustalO:
    aligner: ClustalO
    type: "default"
    params:
      "--iterations": 6
      "--max-guidetree-iterations": [1, 2, 3, 4]
      "--max-hmm-iterations": [1, 2]    # -> 4 x 2 = 8 alignments
```

`type: "Muscle5"` uses Muscle5's built-in ensemble mode (`-stratified`, `-replicates`) instead of parameter variants.

## Outputs

After step 1, each dataset directory contains one folder per ensemble member plus the combined `ensemble/` folder:

```text
results/
└── demo/
    └── protein/
        ├── ClustalO/
        │   ├── done
        │   ├── ensemble.log
        │   └── msa.[0-7].fasta
        ├── MAFFT_FFT-NS-2/
        ├── MAFFT_G-INS-i/
        ├── MAFFT_L-INS-i/
        ├── Muscle3/
        ├── Muscle5/
        └── ensemble/
            ├── ClustalO.[0-7].fasta
            ├── MAFFT_FFT-NS-2.[0-7].fasta
            └── ...
```

Step 2 writes per-dataset scores to `results/{source}/{dataset}/stats.parquet` and the combined scores to `results/{source}/stats.parquet`. Each row holds descriptive statistics for one method:

| Method     | Scores summarized                                                    |
| ---------- | -------------------------------------------------------------------- |
| `dpos`     | `d_pos` distance of every pair of alignments in the ensemble         |
| `pythia`   | Pythia difficulty of every alignment (NaN for fewer than 4 sequences) |
| `dpos_ref` | `d_pos` distance of every alignment to `reference.fasta`, if present |

Columns: `mean`, `std`, `min`, `p25`, `p50`, `p75`, `max`, the number of valid scores `n` (failed Pythia runs are excluded), and `datatype` (`DNA`, `AA` or `n/a`).

A `reference.fasta` (an alignment of the same sequences) is looked up next to `sequences.fasta` in the input directory first, then in the output dataset directory. The protein demo dataset ships one.

Only datasets whose ensemble is complete (every tool in `ensemble.yaml` has a `done` flag) are scored. Incomplete ones are skipped with a warning.

To convert a Parquet file to CSV:

```shell
pixi run parquet-tools csv results/demo/stats.parquet > results/demo/stats.csv
```

## Runtime benchmarks

Runtime mode runs the same workflows with every job single-threaded, so CPU times are comparable across tools. Each job also writes its wall-clock and CPU time to a `benchmark.json` (one per tool and dataset) or a `scoring_benchmark.json` (one per dataset):

```shell
pixi run ensemble-runtime -c 8     # = pixi run ensemble -c 8 --config runtime=True
pixi run score-runtime -c 8
```

To compile the JSON files into `ensemble_benchmarks.parquet` and `scoring_benchmarks.parquet` in the output directory:

```shell
pixi run compile-benchmarks                    # both
pixi run compile-benchmarks --kind scoring     # one, optionally with --output path.parquet
```

## License

GPL-3.0-or-later. See [LICENSE](LICENSE).
