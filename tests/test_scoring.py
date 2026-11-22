import numpy as np
import pytest
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ensemblify import scoring
from ensemblify.utils import infer_data_type


def _records(*seqs):
    return [SeqRecord(Seq(s), id=f"s{i}") for i, s in enumerate(seqs)]


def _write_msa(path, rows):
    path.write_text("".join(f">s{i}\n{row}\n" for i, row in enumerate(rows)))


# descriptive_stats -----------------------------------------------------------


def test_stats_basic():
    stats = scoring.descriptive_stats(np.array([0.0, 0.5, 1.0]))
    assert stats["mean"] == pytest.approx(0.5)
    assert stats["min"] == 0.0 and stats["max"] == 1.0
    assert stats["n"] == 3


def test_stats_ignore_nan():
    # One failed RAxML-NG run must not wipe out the scores of all other alignments.
    stats = scoring.descriptive_stats(np.array([0.2, np.nan, 0.4]))
    assert stats["mean"] == pytest.approx(0.3)
    assert stats["n"] == 2


@pytest.mark.parametrize("scores", [np.array([]), np.array([np.nan, np.nan])])
def test_stats_empty(scores):
    stats = scoring.descriptive_stats(scores)
    assert stats["n"] == 0
    assert all(np.isnan(v) for k, v in stats.items() if k != "n")


def test_stats_single_value_has_nan_std():
    stats = scoring.descriptive_stats(np.array([0.7]))
    assert stats["mean"] == 0.7 and np.isnan(stats["std"])


# infer_data_type -------------------------------------------------------------


def test_infer_dna_pools_short_sequences():
    # No single sequence reaches 100 residues, and none contains all four bases.
    seqs = ["ACAC-ACAC" * 5, "GTGT.GTGT" * 5, "ACGG" * 10]
    assert infer_data_type(_records(*seqs)) == "DNA"


def test_infer_protein():
    assert infer_data_type(_records("MKVLAAGIVGLLLAQWERTYIPASDFGHKLCNM" * 4)) == "AA"


def test_infer_empty_is_na():
    assert infer_data_type([]) == "n/a"
    assert infer_data_type(_records("----", "..")) == "n/a"


def test_infer_accepts_msa():
    msa = MultipleSeqAlignment(_records("ACGT-ACGT" * 20, "ACGTTACGT" * 20))
    assert infer_data_type(msa) == "DNA"


# load_ensemble / score_ensemble ----------------------------------------------


def test_load_ensemble_ignores_non_fasta_and_sorts(tmp_path):
    _write_msa(tmp_path / "B.0.fasta", ["AC-GT", "ACGGT", "A-CGT"])
    _write_msa(tmp_path / "A.0.fasta", ["ACG-T", "ACGGT", "AC-GT"])
    (tmp_path / ".DS_Store").write_bytes(b"\x00\x01")
    (tmp_path / "notes.txt").write_text("hello")

    ensemble = scoring.load_ensemble(tmp_path)
    assert len(ensemble.alignments) == 2


def test_load_ensemble_empty_dir_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        scoring.load_ensemble(tmp_path)


def test_score_ensemble_small_dataset_skips_pythia(tmp_path):
    # Three sequences: too few for RAxML-NG, so Pythia is NaN without calling it.
    _write_msa(tmp_path / "A.0.fasta", ["ACG-TACGTA", "ACGGTACGTA", "AC-GTACGTA"])
    _write_msa(tmp_path / "A.1.fasta", ["AC-GTACGTA", "ACGGTACGTA", "A-CGTACGTA"])
    _write_msa(tmp_path / "A.2.fasta", ["ACG-TACGTA", "ACGGTACGTA", "AC-GTACGTA"])
    ref = tmp_path.parent / "reference.fasta"
    _write_msa(ref, ["ACG-TACGTA", "ACGGTACGTA", "AC-GTACGTA"])

    ensemble = scoring.load_ensemble(tmp_path)
    stats = scoring.score_ensemble(ensemble, reference_file=ref, raxml_path=None)

    assert list(stats.index) == ["dpos", "pythia", "dpos_ref"]
    assert stats.loc["dpos", "n"] == 3  # 3 alignments -> 3 pairs
    assert stats.loc["dpos", "min"] == 0.0  # A.0 and A.2 are identical
    assert stats.loc["pythia", "n"] == 0
    assert stats.loc["dpos_ref", "n"] == 3
    assert (stats["datatype"] == "DNA").all()
