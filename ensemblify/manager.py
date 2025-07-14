from abc import ABC
from ensemblify import utils
from Bio.SeqRecord import SeqRecord
from pathlib import Path

class BaseManager(ABC):
    # Job: execute cmd 
    pass

    def compute(self):
        # Compute alignments associated with given config
        # Output should always be a collection of .fasta files in the correct output directory 
        pass

class DefaultManager(BaseManager):
    def __init__(self, config:dict, tool_name:str, in_file:Path, out_dir:Path, log_file:Path, threads:int):
        self.config = config
        self.tool_name = tool_name
        self.in_file = in_file
        self.out_dir = out_dir
        self.log_file = log_file
        self.threads = threads
        self.variants = utils.get_option_variants(self.config)

    def compute(self):
        for i, options in enumerate(self.variants):
            out_file = self.output / f"msa.{i}.fasta"
            
            # Mafft does not have an -output flag, we have to use stdout as intended
            if not "{output}" in self.config["aligner"]["cmd"]:
                stdout = out_file
            else:
                # All other aligners have dedicated -output flags, using log file as stdout to avoid bugs
                stdout = self.log
            
            in_file = 
            cmd = utils.format_cmd(self.config, self.tool_name, options, self.threads, in_file, out_file,)
            utils.run_cmd(
                cmd, outfile=stdout, logfile=self.log, log_write_mode="a"
            )
            if self._benchmarking:
                stop = perf_counter()
                self._perf_dict[str(tool_enum)].append(stop - start)
            alignment = Alignment(out_file)
            alignments.append(alignment)

        ensemble = Ensemble(alignments)
        ensemble.to_efa(out_file)
        return ensemble


class Muscle5Manager(BaseManager):
    pass



def infer_manager(aligner:str):
    match(aligner):
        case "Muscle5":
            return Muscle5Manager
        case _: 
            return DefaultManager
        # TODO: Add more manager classes for other alignment programs


class Aligner:
    def __init__(
        self,
        dataset: list[SeqRecord],
        seed: int = 0,
        # benchmarking: bool = False,
    ):
        self._dataset = dataset
        self._seed = seed
        # self._benchmarking = benchmarking
        # if benchmarking:
        #     self._perf_dict = defaultdict(list)

    def compute(
        self,
        config: dict,
        out_file: Path,
        log_file: Path,
        threads: int,
    ):
        if config.aligner == "default":
            return self._compute_muscle5(tool_enum, out_file, log_file, threads)
        else:
            return self._compute_singletons(tool_enum, out_file, log_file, threads)

    def _compute_singletons(
        self,
        tool_enum: ToolEnum,
        out_file: pathlib.Path,
        log_file: pathlib.Path,
        threads: int,
    ):
        alignments = []
        # get number of option combinations:
        for i, options in enumerate(ToolEnum.get_option_variants(tool_enum)):
            file_name = options.replace(SPACE_CHAR, "_").strip(" ")
            if not file_name:
                file_name = "_"
            msa_out_file = out_file.parent / f"msa.{i}.{file_name}.fasta"
            # Mafft does not have an -output flag, we have to use stdout as intended
            if tool_enum.aligner == AlignerEnum.MAFFT:
                stdout = msa_out_file
            else:
                # All other aligners have dedicated -output flags, using log file as stdout to avoid bugs
                stdout = log_file

            cmd = self._format_cmd(tool_enum, msa_out_file, options, threads, index=i)
            if self._benchmarking:
                start = perf_counter()
            features.utils.run_cmd(
                cmd, outfile=stdout, logfile=log_file, log_write_mode="a"
            )
            if self._benchmarking:
                stop = perf_counter()
                self._perf_dict[str(tool_enum)].append(stop - start)
            alignment = Alignment(msa_out_file)
            alignments.append(alignment)

        ensemble = Ensemble(alignments)
        ensemble.to_efa(out_file)
        return ensemble

    def _compute_muscle5(
        self,
        muscle5: ToolEnum,
        out_file: pathlib.Path,
        log_file: pathlib.Path,
        threads: int,
    ):
        options = ToolEnum.get_option_variants(muscle5)
        if len(options) > 1:
            raise NotImplementedError(
                "Combinatorial options not implemented for Muscle5"
            )
        options = options[0]
        # Using -stratified for 16 alignments at once
        cmd = self._format_cmd(muscle5, out_file, options, threads)
        if self._benchmarking:
            start = perf_counter()
        features.utils.run_cmd(
            cmd, outfile=log_file, logfile=log_file, log_write_mode="a"
        )
        if self._benchmarking:
            stop = perf_counter()
            self._perf_dict[str(muscle5)].append(stop - start)
        ensemble = Ensemble.from_efa(AlignerEnum.MUSCLE5.path, out_file)
        return ensemble

    def _format_cmd(
        self,
        tool_enum: ToolEnum,
        out_file: pathlib.Path,
        options: str,
        threads: int,
        index: int = None,
    ):
        shuffle_seqs = index is not None
        if shuffle_seqs:
            sequence_records = self._dataset.sequences
            random.shuffle(sequence_records)
            seq_file = out_file.parent / f"sequences.{index}.fasta"
            SeqIO.write(sequence_records, seq_file, format="fasta")
        else:
            seq_file = self._dataset.file_path

        kwargs = dict(
            aligner=str(tool_enum.aligner.path),
            input=str(seq_file),
            output=str(out_file),
            options=options,
            threads=threads,
        )
        return ToolEnum.format_cmd(tool_enum, kwargs)
