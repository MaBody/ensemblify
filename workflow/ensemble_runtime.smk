configfile: "configs/general.yaml"
configfile: "configs/aligners.yaml"
configfile: "configs/ensemble.yaml"

import json
import os
import pathlib
import resource
import time
from datetime import datetime, timezone
from pathlib import Path

from ensemblify.manager import infer_manager


INPUT_DIR = pathlib.Path(config["general"]["input"])
OUT_DIR = pathlib.Path(config["general"]["output"])

TOOLS = list(config["ensemble"].keys())
SOURCES = config["general"]["sources"]
DATASET_MAP = {
    source: list(
        filter(
            lambda name: os.path.isdir(INPUT_DIR / source / name),
            os.listdir(INPUT_DIR / source),
        )
    )
    for source in SOURCES
}


wildcard_constraints:
    source="|".join(SOURCES),
    dataset="|".join([dataset for datasets in DATASET_MAP.values() for dataset in datasets]),
    tool="|".join(TOOLS),


def _expanded_dataset_tool_paths(file_name):
    return expand(
        OUT_DIR / "{source}" / "{dataset}" / "{tool}" / file_name,
        zip,
        source=[
            source
            for source in DATASET_MAP
            for dataset in DATASET_MAP[source]
            for tool in TOOLS
        ],
        dataset=[
            dataset
            for source in DATASET_MAP
            for dataset in DATASET_MAP[source]
            for tool in TOOLS
        ],
        tool=[
            tool
            for source in DATASET_MAP
            for dataset in DATASET_MAP[source]
            for tool in TOOLS
        ],
    )


def _cpu_snapshot():
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    return {
        "process": time.process_time(),
        "children": children.ru_utime + children.ru_stime,
    }


def _cpu_elapsed(start, end):
    return {
        "process_seconds": end["process"] - start["process"],
        "children_seconds": end["children"] - start["children"],
        "total_seconds": (
            end["process"]
            - start["process"]
            + end["children"]
            - start["children"]
        ),
    }


def _write_benchmark(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


rule all:
    input:
        ensemble_done = _expanded_dataset_tool_paths("done"),
        ensemble_benchmark = _expanded_dataset_tool_paths("benchmark.json")


rule generate_ensemble_benchmark:
    threads: lambda wildcards: config["aligners"][config["ensemble"][wildcards.tool]["aligner"]]["threads"]
    input:
        in_file = INPUT_DIR / "{source}" / "{dataset}" / "sequences.fasta"
    output:
        done = OUT_DIR / "{source}" / "{dataset}" / "{tool}" / "done",
        benchmark = OUT_DIR / "{source}" / "{dataset}" / "{tool}" / "benchmark.json"
    log:
        log_file = OUT_DIR / "{source}" / "{dataset}" / "{tool}" / "ensemble.log"
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

        started_at = datetime.now(timezone.utc)
        wall_start = time.perf_counter()
        cpu_start = _cpu_snapshot()
        status = "success"
        exception_type = None
        exception_message = None

        try:
            ensemble = manager.compute()

            ensemble_dir = out_dir.parent / "ensemble"
            os.makedirs(ensemble_dir, exist_ok=True)
            manager.save_ensemble(ensemble, ensemble_dir)
            Path(output.done).touch()
        except Exception as exc:
            status = "failed"
            exception_type = type(exc).__name__
            exception_message = str(exc)
            raise
        finally:
            ended_at = datetime.now(timezone.utc)
            wall_seconds = time.perf_counter() - wall_start
            cpu_seconds = _cpu_elapsed(cpu_start, _cpu_snapshot())
            _write_benchmark(
                output.benchmark,
                {
                    "benchmark_scope": "ensemble_tool_compute",
                    "source": wildcards.source,
                    "dataset": wildcards.dataset,
                    "tool": wildcards.tool,
                    "threads": threads,
                    "status": status,
                    "started_at": started_at.isoformat(),
                    "ended_at": ended_at.isoformat(),
                    "wall_seconds": wall_seconds,
                    "cpu_seconds": cpu_seconds,
                    "exception_type": exception_type,
                    "exception_message": exception_message,
                },
            )


