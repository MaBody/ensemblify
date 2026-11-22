from abc import ABC
from ensemblify import utils
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from pathlib import Path
import hashlib
import random
import os
import re


class Manager(ABC):
    def __init__(
        self,
        config: dict,
        tool_name: str,
        in_file: Path,
        out_dir: Path,
        log_file: Path,
        threads: int,
        shuffle: bool = False,
        seed: int | None = None,
    ):
        """
        Parameters
        ----------
        shuffle : bool
            Shuffle the input sequence order before each aligner run, so the
            ensemble does not inherit a bias from the order in the input file.
        seed : int, optional
            Makes the shuffling reproducible. `None` shuffles differently on every run.
        """
        self.config = config
        self.tool_name = tool_name
        self.in_file = in_file
        self.out_dir = Path(out_dir)
        self.log_file = log_file
        self.threads = threads
        self.shuffle = shuffle
        self.seed = seed
        self.variants = utils.get_option_variants(
            self.config["ensemble"][tool_name]["params"]
        )

    def _rng(self, variant: int) -> random.Random:
        """Random generator for one aligner run.

        Derived from the seed, the input's content, the tool and the variant index:
        reproducible across machines and input paths, and independent between the
        parallel jobs of one workflow run. (`random.Random` hashes string seeds with
        SHA-512, so this does not depend on PYTHONHASHSEED.)
        """
        if self.seed is None:
            return random.Random()
        digest = hashlib.sha256(Path(self.in_file).read_bytes()).hexdigest()
        return random.Random(f"{self.seed}:{digest}:{self.tool_name}:{variant}")

    def save_ensemble(
        self, ensemble_dict: dict[tuple, MultipleSeqAlignment], out_dir: Path
    ):
        """Write this tool's alignments to `out_dir` as `{tool}.{i}.fasta`.

        Earlier members of the same tool are removed first. Otherwise a rerun with
        fewer variants would leave stale alignments behind, and the scoring workflow
        reads every file in the ensemble directory.
        """
        out_dir = Path(out_dir)
        stale = re.compile(rf"^{re.escape(self.tool_name)}\.[0-9]+\.fasta$")
        for name in os.listdir(out_dir):
            if stale.match(name):
                (out_dir / name).unlink()
        for keys, msa in ensemble_dict.items():
            AlignIO.write(msa, out_dir / (".".join(map(str, keys)) + ".fasta"), "fasta")

    def compute(self) -> dict[tuple, MultipleSeqAlignment]:
        # Compute alignments associated with given config
        # Output should always be a collection of .fasta files in the correct output directory
        raise NotImplementedError()


class DefaultManager(Manager):
    def compute(self) -> dict[tuple, MultipleSeqAlignment]:
        self.alignments = {}
        for i, options in enumerate(self.variants):
            out_file = self.out_dir / f"msa.{i}.fasta"

            aligner = self.config["ensemble"][self.tool_name]["aligner"]
            # Mafft does not have an -output flag, we have to use stdout as intended
            if not "{output}" in self.config["aligners"][aligner]["cmd"]:
                stdout = out_file
            else:
                # All other aligners have dedicated -output flags, using log file as stdout to avoid bugs
                stdout = self.log_file

            # The shuffled copy lives next to the outputs (not in /tmp, which can be
            # small on cluster nodes) and is removed even if the aligner fails.
            shuffled = None
            if self.shuffle:
                seq_records = list(SeqIO.parse(self.in_file, "fasta"))
                self._rng(i).shuffle(seq_records)
                shuffled = self.out_dir / f".input.{i}.fasta"
                SeqIO.write(seq_records, shuffled, format="fasta")
            in_file = shuffled or self.in_file

            try:
                cmd = utils.format_cmd(
                    self.config, self.tool_name, options, self.threads, in_file, out_file
                )
                utils.run_cmd(
                    cmd, outfile=stdout, logfile=self.log_file, log_write_mode="a"
                )
            finally:
                if shuffled is not None:
                    shuffled.unlink(missing_ok=True)

            msa = AlignIO.read(out_file, "fasta")
            self.alignments[(self.tool_name, i)] = msa

        return self.alignments


class Muscle5Manager(Manager):
    # Muscle5 replaces '@' in the output template with "<permutation>.<replicate>",
    # e.g. msa.abc.3.fasta. Permutations: none, abc, acb, bca.
    _REGEX = re.compile(r"^msa\.([a-z]{3,4})\.([0-9]+)\.fasta$")

    def compute(self) -> dict[tuple, MultipleSeqAlignment]:
        self.alignments = {}

        in_file = self.in_file
        out_template = self.out_dir / "msa.@.fasta"

        # Remove outputs of earlier runs so they cannot leak into this ensemble.
        for name in os.listdir(self.out_dir):
            if self._REGEX.match(name):
                (self.out_dir / name).unlink()

        assert (
            len(self.variants) == 1
        ), "Explicit Muscle5 variants are not supported (internal variation)."
        options = self.variants[0]

        cmd = utils.format_cmd(
            self.config, self.tool_name, options, self.threads, in_file, out_template
        )
        utils.run_cmd(
            cmd, outfile=self.log_file, logfile=self.log_file, log_write_mode="a"
        )

        # Sort by (permutation, replicate) so indices are stable across runs and machines.
        matches = [self._REGEX.match(name) for name in os.listdir(self.out_dir)]
        matches = sorted(
            (m for m in matches if m), key=lambda m: (m.group(1), int(m.group(2)))
        )
        if not matches:
            raise RuntimeError(
                f"Muscle5 produced no alignments in {self.out_dir}. "
                f"Check the log file {self.log_file}."
            )
        for i, match in enumerate(matches):
            msa = AlignIO.read(self.out_dir / match.group(0), "fasta")
            self.alignments[(self.tool_name, i)] = msa

        return self.alignments


def infer_manager(aligner: str) -> Manager:
    match (aligner):
        case "Muscle5":
            return Muscle5Manager
        case "default":
            return DefaultManager
        case _:
            raise NotImplementedError(f"Aligner type '{aligner}'")
        # TODO: Add more manager classes for other alignment programs
