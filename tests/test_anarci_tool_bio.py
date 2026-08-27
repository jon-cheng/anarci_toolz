"""Tier 2: bio-dependent tests that call real ANARCI. Requires ANARCI/HMMER
actually installed (see CLAUDE.md's Docker/environment notes) — skipped
cleanly otherwise so Tier 1 stays runnable without these deps.
"""

import pytest

from conftest import require_real_module

require_real_module("anarci")

from anarci_toolz.anarci_tool import (
    ANARCI_RESULTS_SCHEMA,
    get_anarci_alignment,
    run_parallel_anarci,
)
from anarci_toolz.utils import create_row_id

pytestmark = pytest.mark.bio


class TestGetAnarciAlignmentReal:
    def test_real_sequences_pass(self, therasabdab_subset_df):
        for _, row in therasabdab_subset_df.head(3).iterrows():
            result = get_anarci_alignment(
                row["sequence_aa"], row["Therapeutic"], "imgt", allowed_species=["human"]
            )
            _, passed_anarci, start, end, chain_type, v_call, j_call, e_value, _, _ = result

            assert passed_anarci is True
            assert chain_type in {"H", "K", "L"}
            assert isinstance(v_call, str) and v_call.startswith("IG")
            assert isinstance(j_call, str) and j_call.startswith("IG")
            assert isinstance(e_value, float) and e_value >= 0
            assert start < end

    @pytest.mark.parametrize("garbage_seq", ["AAAAAAAAAAAAAAAAAAAA", ""])
    def test_garbage_sequence_fails_gracefully(self, garbage_seq):
        result = get_anarci_alignment(
            garbage_seq, "garbage", "imgt", allowed_species=["human"]
        )
        name, passed_anarci = result[0], result[1]

        assert name == "garbage"
        assert passed_anarci is False
        assert all(v is None for v in result[2:])


class TestRunParallelAnarciReal:
    def test_schema_and_real_sequences_pass(self, therasabdab_subset_df):
        df, seqs = create_row_id(therasabdab_subset_df.copy(), "sequence_aa")

        result = run_parallel_anarci(
            df, seqs, scheme="imgt", allowed_species=["human"], seq_aa_header="sequence_aa"
        )

        assert list(result.columns) == ["seq_id"] + ANARCI_RESULTS_SCHEMA
        assert len(result) == len(df)
        assert result["passed_anarci"].all()
        assert result["e_value"].notna().all()
