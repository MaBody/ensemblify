from abc import ABC
from ensemblify import utils
from Bio import SeqRecord
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from pathlib import Path
from tempfile import NamedTemporaryFile
import random
import os


class Manager(ABC):
    def __init__(self, config:dict, tool_name:str, in_file:Path, out_dir:Path, log_file:Path, threads:int, shuffle:bool=False):
        self.config = config
        self.tool_name = tool_name
        self.in_file = in_file
        self.out_dir = out_dir
        self.log_file = log_file
        self.threads = threads
        self.shuffle = shuffle
        self.variants = utils.get_option_variants(self.config["ensemble"][tool_name]["params"])


    @classmethod
    def save_ensemble(cls, ensemble_dict: dict[tuple, MultipleSeqAlignment], out_dir:Path):
        for keys, msa in ensemble_dict.items():
            AlignIO.write(msa, out_dir / ("-".join(map(str, keys)) + ".fasta"), "fasta")


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

            cmd = utils.format_cmd(self.config, self.tool_name, options, self.threads, in_file, out_file)
            utils.run_cmd(cmd, outfile=stdout, logfile=self.log_file, log_write_mode="a")

            msa = AlignIO.read(out_file, "fasta")
            self.alignments[(self.tool_name, i)] = msa
            
            if self.shuffle:
                temp_file.close()
    
        return self.alignments


class Muscle5Manager(Manager):
    pass



def infer_manager(aligner:str) -> Manager:
    match(aligner):
        case "Muscle5":
            return Muscle5Manager
        case _: 
            return DefaultManager
        # TODO: Add more manager classes for other alignment programs


   

    # def _compute_muscle5(
    #     self,
    #     muscle5: ToolEnum,
    #     out_file: pathlib.Path,
    #     log_file: pathlib.Path,
    #     threads: int,
    # ):
    #     options = ToolEnum.get_option_variants(muscle5)
    #     if len(options) > 1:
    #         raise NotImplementedError(
    #             "Combinatorial options not implemented for Muscle5"
    #         )
    #     options = options[0]
    #     # Using -stratified for 16 alignments at once
    #     cmd = self._format_cmd(muscle5, out_file, options, threads)
    #     if self._benchmarking:
    #         start = perf_counter()
    #     features.utils.run_cmd(
    #         cmd, outfile=log_file, logfile=log_file, log_write_mode="a"
    #     )
    #     if self._benchmarking:
    #         stop = perf_counter()
    #         self._perf_dict[str(muscle5)].append(stop - start)
    #     ensemble = Ensemble.from_efa(AlignerEnum.MUSCLE5.path, out_file)
    #     return ensemble

    # def _format_cmd(
    #     self,
    #     tool_enum: ToolEnum,
    #     out_file: pathlib.Path,
    #     options: str,
    #     threads: int,
    #     index: int = None,
    # ):
    #     shuffle_seqs = index is not None
    #     if shuffle_seqs:
    #         sequence_records = self._dataset.sequences
    #         random.shuffle(sequence_records)
    #         seq_file = out_file.parent / f"sequences.{index}.fasta"
    #         SeqIO.write(sequence_records, seq_file, format="fasta")
    #     else:
    #         seq_file = self._dataset.file_path

    #     kwargs = dict(
    #         aligner=str(tool_enum.aligner.path),
    #         input=str(seq_file),
    #         output=str(out_file),
    #         options=options,
    #         threads=threads,
    #     )
    #     return ToolEnum.format_cmd(tool_enum, kwargs)
