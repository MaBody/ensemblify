configfile: "configs/general.yaml"
configfile: "configs/aligners.yaml"
configfile: "configs/ensemble.yaml"

# Generate one MSA ensemble per dataset.
#
#   Input:  {input}/{source}/{dataset}/sequences.fasta
#   Output: {output}/{source}/{dataset}/{tool}/msa.*.fasta  (+ done flag, ensemble.log)
#           {output}/{source}/{dataset}/ensemble/{tool}.{i}.fasta
#
# Runtime mode (`--config runtime=True`, or `pixi run ensemble-runtime`) runs every
# job single-threaded and writes a benchmark.json with wall/CPU times next to each
# done flag.

import pathlib

from ensemblify import discovery
from ensemblify.benchmark import record_runtime
from ensemblify.manager import infer_manager
from ensemblify.utils import check_aligners


INPUT_DIR = pathlib.Path(config["general"]["input"])
OUT_DIR = pathlib.Path(config["general"]["output"])
SEED = config["general"].get("seed")
RUNTIME = bool(config.get("runtime", False))

TOOLS = list(config["ensemble"].keys())
check_aligners(config)

SOURCES = config["general"]["sources"]
DATASET_MAP = discovery.input_datasets(INPUT_DIR, SOURCES)

wildcard_constraints:
    source=discovery.wildcard_regex(SOURCES),
    dataset=discovery.wildcard_regex(d for datasets in DATASET_MAP.values() for d in datasets),
    tool=discovery.wildcard_regex(TOOLS),


TOOL_DIR = OUT_DIR / "{source}" / "{dataset}" / "{tool}"
OUTPUTS = {"done": TOOL_DIR / "done"}
if RUNTIME:
    OUTPUTS["benchmark"] = TOOL_DIR / "benchmark.json"


def _aligner_threads(wildcards):
    if RUNTIME:
        return 1
    return config["aligners"][config["ensemble"][wildcards.tool]["aligner"]]["threads"]


rule all:
    input:
        [
            str(path).format(source=source, dataset=dataset, tool=tool)
            for source, datasets in DATASET_MAP.items()
            for dataset in datasets
            for tool in TOOLS
            for path in OUTPUTS.values()
        ]


rule generate_ensemble:
    threads: _aligner_threads
    input:
        in_file = INPUT_DIR / "{source}" / "{dataset}" / "sequences.fasta"
    output:
        **OUTPUTS
    log:
        log_file = TOOL_DIR / "ensemble.log"
    run:
        out_dir = pathlib.Path(output.done).parent
        manager_class = infer_manager(config["ensemble"][wildcards.tool]["type"])
        manager = manager_class(
            config, wildcards.tool, input.in_file, out_dir, pathlib.Path(log.log_file),
            threads, shuffle=True, seed=SEED,
        )
        with record_runtime(
            output.get("benchmark"),
            benchmark_scope="ensemble_tool_compute",
            source=wildcards.source,
            dataset=wildcards.dataset,
            tool=wildcards.tool,
            threads=threads,
        ):
            ensemble = manager.compute()
            ensemble_dir = out_dir.parent / "ensemble"
            ensemble_dir.mkdir(parents=True, exist_ok=True)
            manager.save_ensemble(ensemble, ensemble_dir)
            pathlib.Path(output.done).touch()
