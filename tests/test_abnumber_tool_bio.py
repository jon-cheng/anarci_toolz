"""Tier 2: bio-dependent tests that call real AbNumber. Requires
ANARCI/AbNumber/HMMER actually installed — skipped cleanly otherwise so
Tier 1 stays runnable without these deps.
"""

import pytest

from conftest import require_real_module

require_real_module("abnumber")

from anarci_toolz.abnumber_tool import (
    ABNUMBER_TOOL_SCHEMA,
    create_chain,
    get_region_seqs,
    run_parallel_abnumber,
)
from anarci_toolz.utils import create_row_id

pytestmark = pytest.mark.bio

TIMIGUTUZUMAB_SEQ = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTIS"
    "ADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)


class TestCreateChainReal:
    def test_known_good_sequence_resolves(self):
        chain = create_chain(
            TIMIGUTUZUMAB_SEQ, "Timigutuzumab", "imgt", allowed_species=["human"]
        )

        assert chain.chain_type == "H"
        assert chain.species == "human"
        assert chain.v_gene == "IGHV3-66*01"
        assert chain.j_gene == "IGHJ4*01"


class TestGetRegionSeqsReal:
    def test_real_sequences_populate_all_regions(self, therasabdab_subset_df):
        for _, row in therasabdab_subset_df.head(3).iterrows():
            result = get_region_seqs(
                row["sequence_aa"],
                row["Therapeutic"],
                "imgt",
                ["human"],
                "fr1_seq",
                "cdr1_seq",
                "fr2_seq",
                "cdr2_seq",
                "fr3_seq",
                "cdr3_seq",
                "fr4_seq",
            )
            region_info = result[0][row["Therapeutic"]]
            assert region_info["passed_abnumber"] is True

            for region in ("fr1_seq", "cdr1_seq", "fr2_seq", "cdr2_seq", "fr3_seq", "cdr3_seq", "fr4_seq"):
                seq_data = region_info[region]
                assert seq_data["sequence"]
                assert seq_data["start_idx"] < seq_data["end_idx"]

    @pytest.mark.parametrize("garbage_seq", ["AAAAAAAAAAAAAAAAAAAA", ""])
    def test_garbage_sequence_fails_gracefully(self, garbage_seq):
        result = get_region_seqs(
            garbage_seq, "garbage", "imgt", ["human"], "fr1_seq", "cdr1_seq"
        )
        region_info = result[0]["garbage"]

        assert region_info["passed_abnumber"] is False
        for region in ("fr1_seq", "cdr1_seq"):
            assert region_info[region]["sequence"] is None
            assert region_info[region]["start_idx"] is None
            assert region_info[region]["end_idx"] is None


class TestRunParallelAbnumberReal:
    def test_schema_and_readme_spot_check(self, therasabdab_subset_df):
        df, seqs = create_row_id(therasabdab_subset_df.copy(), "sequence_aa")

        result = run_parallel_abnumber(
            df,
            seqs,
            scheme="imgt",
            seq_aa_header="sequence_aa",
            seq_dna_header=None,
            allowed_species=["human"],
            retain_indices=False,
            display_residue_view=False,
            display_unformatted_residue_view=False,
            num_cpu=1,
        )

        for col in ABNUMBER_TOOL_SCHEMA:
            assert col in result.columns

        # README worked example, Timigutuzumab row (seq_id "1" = first row of
        # the fixture, since create_row_id enumerates from 1)
        timi_row = result[result["seq_id"] == "1"].iloc[0]
        assert timi_row["v_gene"] == "IGHV3-66*01"
        assert timi_row["j_gene"] == "IGHJ4*01"
        assert timi_row["cdr3_aa"] == "SRWGGDGFYAMDY"
        assert timi_row["chain_type"] == "H"
