import itertools
import tempfile
import subprocess
import numpy as np
import pandas as pd
from collections import defaultdict
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
        aligner=str(config["aligners"][aligner_name]["path"]),
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
    Runs a given command using the subprocess Python module. The results of the run are stored in the given `outfile`
    and all logging and errors are written to `logfile` using the given write mode.

    Parameters
    ----------
    cmd : list[str]
        Command to run as list of strings.
    outfile : pathlib.Path
        File where to store the output of the given command execution in.
    logfile : pathlib.Path
        File where to store all logging/errors in.
    log_write_mode : str
        Write mode to open the logfile with.
        Allowed options are `'a'` (append logs to existing logs) and `'w'` (overwrite all existing logs).

    """
    if log_write_mode not in ["a", "w"]:
        raise ValueError(
            f"Invalid write mode for logfile given: {log_write_mode}. Allowed options are 'a' (append) and 'w' (write)."
        )
    if not logfile:
        logfile = Path(tempfile.NamedTemporaryFile("w").name)
    with logfile.open(log_write_mode) as log, outfile.open(out_write_mode) as out:
        try:
            subprocess.run(cmd, stdout=out, stderr=log, check=True, encoding="utf-8")
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"Error running CMD: {' '.join(cmd)}"
                f"Check the logfile {logfile.absolute()} for details."
                f"Logfile: {logfile.open().read()}"
            )


def infer_data_type(records: list[SeqRecord] | MultipleSeqAlignment):
    """
    Infer the biological data type (DNA or protein) from a collection of sequence records.

    Examines the first sufficiently long ungapped sequence and determines the type
    based on character composition.

    Parameters
    ----------
    records : list[SeqRecord] | MultipleSeqAlignment
        Collection of sequence records.

    Returns
    -------
    data_type : Literal['DNA', 'AA', 'n/a']
        Inferred data type.

    Raises TODO: Remove!
    ------
    AssertionError
        If the alphabet cannot be determined.
    """
    ungapped = None
    for seq in records:
        ungapped = str(seq.seq).replace(GAP_CHAR, "")[:1000].upper()
        if len(ungapped) >= 100:
            break
    char_counts = pd.Series(list(ungapped)).value_counts()
    unique_chars = char_counts.index.to_series()

    dna_included = np.isin(DNA_CHARS, unique_chars)
    protein_included = np.isin(AA_CHARS, unique_chars)

    protein_missing = AA_CHARS[~protein_included]
    missing_counts = pd.Series(0, index=protein_missing)
    char_counts = pd.concat((char_counts, missing_counts))

    # If all DNA chars are included and their cumulative share is above 90%: DNA!
    if dna_included.all() & (char_counts[DNA_CHARS].sum() / len(ungapped) >= 0.9):
        return "DNA"
    # If the cumulative share of protein chars is above 90%: Protein!
    elif char_counts[AA_CHARS].sum() / len(ungapped) >= 0.9:
        return "AA"
    else:
        return "n/a"
        # raise AssertionError(
        #     "Unknown alphabet detected for residues '{}'".format(ungapped)
        # )
