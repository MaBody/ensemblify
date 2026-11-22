import itertools
import os
import shlex
import shutil
import subprocess
from contextlib import ExitStack
import numpy as np
import pandas as pd
from typing import Literal
from pathlib import Path
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

_SPACER = "§"

GAP_CHAR = "-"
DNA_CHARS = np.array(list("ACGT"))
AA_CHARS = np.array(list("ACDEFGHIKLMNPQRSTVWY"))


def equals_none(string: str):
    if string is None:
        return True
    string = string.strip().lower()
    if (string == "null") | (string == "none") | (string == "") | (string == "~"):
        return True
    return False


def get_option_variants(params: dict):
    """Returns a list with the cartesian product of all combinations of options, formatted as strings."""
    if params is None:
        raise ValueError("Missing 'params' entry.")
        # Return empty string to avoid cmd error (necessary ?)
        # return [""]

    variants = []
    base_options = []
    product_options = {}
    for key, val in params.items():
        if equals_none(str(val)):
            # If val is none: simple flag
            base_options.append(key)
        elif isinstance(val, list):
            # If val is list: contains options
            product_options[key] = val
        else:
            # val is assumed to be a single value
            base_options.append(key)
            base_options.append(str(val))

    product_keys = list(product_options.keys())
    for product_vals in itertools.product(*product_options.values()):
        variant = base_options.copy()
        for key, val in zip(product_keys, product_vals):
            variant.extend([key, str(val)])
        variants.append(variant)
    variants = [_SPACER.join(variant) for variant in variants]
    return variants


def resolve_executable(path: str | Path) -> Path:
    """
    Resolve an executable given in a config file.

    A bare name such as ``"mafft"`` is looked up on ``PATH`` (the pixi environment
    puts the conda-installed tools there). Anything containing a path separator is
    treated as an explicit location, so manually installed binaries keep working.
    Relative locations are resolved against the current working directory, which is
    the repository root when the workflows run through ``pixi run``.

    Raises
    ------
    FileNotFoundError
        If the executable cannot be found or is not executable.
    """
    if path is None or str(path).strip() == "":
        raise FileNotFoundError("No executable path configured.")
    path_str = os.path.expanduser(str(path).strip())

    if os.sep not in path_str:
        found = shutil.which(path_str)
        if found is None:
            raise FileNotFoundError(
                f"Executable '{path_str}' not found on PATH. Run inside the pixi "
                f"environment (`pixi run ...`) or set an explicit path in the config."
            )
        return Path(found)

    resolved = Path(path_str).resolve()
    if not resolved.is_file() or not os.access(resolved, os.X_OK):
        raise FileNotFoundError(
            f"Executable '{path_str}' (resolved to '{resolved}') does not exist or "
            f"is not executable. Check the path in the config."
        )
    return resolved


def check_aligners(config: dict) -> None:
    """Fail early if any aligner referenced in the ensemble config cannot be resolved."""
    errors = []
    for tool_name, tool_cfg in config["ensemble"].items():
        aligner = tool_cfg["aligner"]
        if aligner not in config["aligners"]:
            errors.append(f"{tool_name}: aligner '{aligner}' missing from aligners config")
            continue
        try:
            resolve_executable(config["aligners"][aligner]["path"])
        except FileNotFoundError as e:
            errors.append(f"{tool_name} ({aligner}): {e}")
    if errors:
        raise FileNotFoundError(
            "Aligner configuration problems:\n  " + "\n  ".join(errors)
        )


def format_cmd(
    config: dict,
    tool_name: str,
    options: str,
    threads: int,
    in_file: Path,
    out_file: Path,
):

    aligner_name = config["ensemble"][tool_name]["aligner"]
    kwargs = dict(
        aligner=str(resolve_executable(config["aligners"][aligner_name]["path"])),
        input=str(in_file),
        output=str(out_file),
        options=options,
        threads=threads,
    )
    cmd_template = config["aligners"][aligner_name]["cmd"]
    cmd_template = _SPACER.join(map(str.strip, cmd_template.split(" ")))
    temp = cmd_template.format(**kwargs).split(_SPACER)
    commands = []
    for command in temp:
        section = command.strip(" ")
        if section:
            commands.append(section)
    return commands


def run_cmd(
    cmd: list[str],
    outfile: Path,
    logfile: Path = None,
    log_write_mode: Literal["a", "w"] = "a",
    out_write_mode: Literal["a", "w"] = "w",
) -> None:
    """
    Run a command, writing its stdout to `outfile` and its stderr to `logfile`.

    Parameters
    ----------
    cmd : list[str]
        Command to run as list of strings.
    outfile : pathlib.Path
        File that receives the command's stdout.
    logfile : pathlib.Path, optional
        File that receives the command's stderr. If it is the same file as `outfile`,
        stdout and stderr share one handle opened with `log_write_mode`, so repeated
        calls append to the log instead of truncating it.
    log_write_mode, out_write_mode : {"a", "w"}
        Append to or overwrite the log and output files.

    Raises
    ------
    RuntimeError
        If the command exits with a non-zero code. The message contains the command
        and the last lines of its error output.
    """
    for name, mode in [("log_write_mode", log_write_mode), ("out_write_mode", out_write_mode)]:
        if mode not in ("a", "w"):
            raise ValueError(f"Invalid {name} '{mode}'. Allowed: 'a' (append), 'w' (write).")

    outfile = Path(outfile)
    logfile = Path(logfile) if logfile else None
    shared = logfile is not None and outfile.resolve() == logfile.resolve()

    with ExitStack() as stack:
        if shared:
            out = stack.enter_context(outfile.open(log_write_mode))
            err = subprocess.STDOUT
        else:
            out = stack.enter_context(outfile.open(out_write_mode))
            err = stack.enter_context(logfile.open(log_write_mode)) if logfile else subprocess.PIPE
        proc = subprocess.run(cmd, stdout=out, stderr=err, encoding="utf-8")

    if proc.returncode != 0:
        if logfile is None:
            details = (proc.stderr or "").strip()
        else:
            details = f"Log file: {logfile.absolute()}\nLast lines:\n{_tail(logfile)}"
        raise RuntimeError(
            f"Command failed with exit code {proc.returncode}: {shlex.join(map(str, cmd))}\n{details}"
        )


def _tail(path: Path, n_lines: int = 20) -> str:
    try:
        lines = path.read_text(errors="replace").splitlines()
    except OSError:
        return "<log file not readable>"
    return "\n".join(lines[-n_lines:])


def infer_data_type(
    records: list[SeqRecord] | MultipleSeqAlignment, max_residues: int = 10_000
) -> Literal["DNA", "AA", "n/a"]:
    """
    Infer the biological data type (DNA or protein) from a collection of sequence records.

    Pools the ungapped residues of the records (up to `max_residues`) and decides by
    character composition: DNA if all four bases occur and make up at least 90% of
    the residues, protein if the 20 amino acids make up at least 90%, else "n/a".

    Parameters
    ----------
    records : list[SeqRecord] | MultipleSeqAlignment
        Collection of sequence records.
    max_residues : int
        Upper bound on the number of residues inspected.

    Returns
    -------
    data_type : Literal['DNA', 'AA', 'n/a']
        Inferred data type. "n/a" for empty input or an unknown alphabet.
    """
    chunks = []
    total = 0
    for record in records:
        residues = str(record.seq).upper().replace(GAP_CHAR, "").replace(".", "")
        chunks.append(residues[: max_residues - total])
        total += len(chunks[-1])
        if total >= max_residues:
            break
    residues = "".join(chunks)
    if not residues:
        return "n/a"

    counts = pd.Series(list(residues)).value_counts()

    def share(chars: np.ndarray) -> float:
        return counts.reindex(chars, fill_value=0).sum() / len(residues)

    all_bases_present = all(counts.get(base, 0) > 0 for base in DNA_CHARS)
    if all_bases_present and share(DNA_CHARS) >= 0.9:
        return "DNA"
    if share(AA_CHARS) >= 0.9:
        return "AA"
    return "n/a"
