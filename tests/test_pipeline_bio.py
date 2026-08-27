"""Tier 2: the top-level integration tests, calling real ANARCI + AbNumber
through the full `run_anarci_toolz` / CLI entry points. Requires
ANARCI/AbNumber/HMMER actually installed — skipped cleanly otherwise so
Tier 1 stays runnable without these deps.
"""

import pandas as pd
import pytest
from Bio.Seq import Seq

from conftest import require_real_module

require_real_module("anarci")
require_real_module("abnumber")

from anarci_toolz import pipeline
from anarci_toolz.pipeline import run_anarci_toolz

pytestmark = pytest.mark.bio

# One arbitrary codon per amino acid, sufficient to build a DNA sequence that
# translates cleanly back via Bio.Seq (not meant to reflect real codon usage).
_AA_TO_CODON = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT",
    "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA", "L": "CTT",
    "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT",
    "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT",
}


def _aa_to_dna(aa_seq: str) -> str:
    return "".join(_AA_TO_CODON[aa] for aa in aa_seq)


class TestRunAnarciToolzEndToEnd:
    def test_full_fixture_matches_readme_worked_example(self, therasabdab_full_df):
        df = therasabdab_full_df.copy()

        result = run_anarci_toolz(
            df,
            scheme="imgt",
            allowed_species=["human"],
            seq_aa_header="sequence_aa",
        )

        assert len(result) == len(therasabdab_full_df)
        assert "seq_id" not in result.columns

        for col in (
            "passed_abnumber",
            "v_gene",
            "j_gene",
            "cdr3_aa",
            "passed_anarci",
            "e_value",
            "bitscore",
            "bias",
        ):
            assert col in result.columns

        timi_row = result[result["Therapeutic"] == "Timigutuzumab"].iloc[0]
        assert timi_row["v_gene"] == "IGHV3-66*01"
        assert timi_row["j_gene"] == "IGHJ4*01"
        assert timi_row["cdr3_aa"] == "SRWGGDGFYAMDY"

        inebi_row = result[result["Therapeutic"] == "Inebilizumab"].iloc[0]
        assert inebi_row["v_gene"] == "IGHV3-66*01"
        assert inebi_row["j_gene"] == "IGHJ4*01"
        assert inebi_row["cdr3_aa"] == "ARSGFITTVRDFDY"

        # e_value/bitscore/bias can drift slightly across HMMER versions;
        # only assert they're present and numerically sane, per CLAUDE.md.
        assert (result["e_value"] >= 0).all()
        assert result["bitscore"].notna().all()

    def test_skip_run_base_anarci_omits_anarci_columns(self, therasabdab_subset_df):
        df = therasabdab_subset_df.copy()

        result = run_anarci_toolz(
            df,
            scheme="imgt",
            allowed_species=["human"],
            seq_aa_header="sequence_aa",
            skip_run_base_anarci=True,
        )

        for col in ("passed_anarci", "e_value", "bitscore", "bias"):
            assert col not in result.columns
        assert "v_gene" in result.columns  # AbNumber-derived columns still present

    def test_dna_mode_slices_regions_consistently_with_aa(self, therasabdab_subset_df):
        df = therasabdab_subset_df.head(3).copy()
        df["sequence_dna"] = df["sequence_aa"].apply(_aa_to_dna)

        # sanity-check our fixture-building helper before trusting the result
        for aa, dna in zip(df["sequence_aa"], df["sequence_dna"]):
            assert str(Seq(dna).translate()) == aa

        result = run_anarci_toolz(
            df,
            scheme="imgt",
            allowed_species=["human"],
            seq_aa_header="sequence_aa",
            seq_dna_header="sequence_dna",
        )

        for region in ("cdr1", "cdr2", "cdr3", "fr1", "fr2", "fr3", "fr4"):
            dna_col, aa_col = f"{region}_dna", f"{region}_aa"
            assert dna_col in result.columns
            for dna_seq, aa_seq in zip(result[dna_col], result[aa_col]):
                assert len(dna_seq) == len(aa_seq) * 3
                assert str(Seq(dna_seq).translate()) == aa_seq


class TestCliEntryPoint:
    def test_main_writes_tagged_output_file(self, tmp_path, monkeypatch, therasabdab_subset_df):
        monkeypatch.chdir(tmp_path)

        input_csv = tmp_path / "mini_sample.csv"
        therasabdab_subset_df.head(2)[["Therapeutic", "sequence_aa"]].to_csv(
            input_csv, index=False
        )

        pipeline.main(
            input=str(input_csv),
            scheme="imgt",
            allowed_species=["human"],
            seq_aa_header="sequence_aa",
            seq_dna_header=None,
            retain_indices=False,
            display_residue_view=False,
            display_unformatted_residue_view=False,
            skip_run_base_anarci=True,
            num_cpu=1,
            tag="_customtag",
        )

        output_file = tmp_path / "anarci_annot" / "mini_sample_customtag.csv"
        assert output_file.exists()

        result = pd.read_csv(output_file)
        assert len(result) == 2
        assert "v_gene" in result.columns
