from abc import ABC
from ensemblify import utils
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from pathlib import Path
from tempfile import NamedTemporaryFile
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
    ):
        self.config = config
        self.tool_name = tool_name
        self.in_file = in_file
        self.out_dir = out_dir
        self.log_file = log_file
        self.threads = threads
        self.shuffle = shuffle
        self.variants = utils.get_option_variants(
            self.config["ensemble"][tool_name]["params"]
        )

    @classmethod
    def save_ensemble(
        cls, ensemble_dict: dict[tuple, MultipleSeqAlignment], out_dir: Path
    ):
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

            if self.shuffle:
                seq_records = list(SeqIO.parse(self.in_file, "fasta"))
                random.shuffle(seq_records)
                temp_file = NamedTemporaryFile("w", delete=False)
                in_file = Path(temp_file.name)
                SeqIO.write(seq_records, in_file, format="fasta")
            else:
                in_file = self.in_file

            cmd = utils.format_cmd(
                self.config, self.tool_name, options, self.threads, in_file, out_file
            )
            utils.run_cmd(
                cmd, outfile=stdout, logfile=self.log_file, log_write_mode="a"
            )

            msa = AlignIO.read(out_file, "fasta")
            self.alignments[(self.tool_name, i)] = msa

            if self.shuffle:
                temp_file.close()

        return self.alignments


class Muscle5Manager(Manager):
    _REGEX = re.compile(r"^[a-z]{3,4}\.[0-9]+\.fasta$")

    def compute(self) -> dict[tuple, MultipleSeqAlignment]:
        self.alignments = {}

        in_file = self.in_file
        out_template = self.out_dir / f"msa.@.fasta"

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

        for i, file_name in enumerate(
            filter(lambda name: re.match(self._REGEX, name), os.listdir(self.out_dir))
        ):
            msa = AlignIO.read(self.out_dir / file_name, "fasta")
            self.alignments[(self.tool_name, i)] = msa

        return self.alignments


def infer_manager(aligner: str) -> Manager:
    match (aligner):
        case "Muscle5":
            return Muscle5Manager
        case "default":
            return DefaultManager
        case _:
            raise NotImplementedError()
        # TODO: Add more manager classes for other alignment programs
