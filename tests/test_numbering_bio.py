"""Tier 2: bio-dependent tests for the wide residue-view formatting, which
needs genuine AbNumber `Chain.to_dataframe` output — Tier 1 skipped these
since they assume real AbNumber output shapes.
"""

import pytest

from conftest import require_real_module

require_real_module("abnumber")

from anarci_toolz.abnumber_tool import REGIONS, get_region_seqs
from anarci_toolz.numbering import get_seq_view_dataframe

pytestmark = pytest.mark.bio


def _build_result(df, scheme):
    """Mirrors what `parallel_get_region_seqs` produces: a list of
    `get_region_seqs` outputs, one per sequence."""
    result = []
    for _, row in df.iterrows():
        result.append(
            get_region_seqs(
                row["sequence_aa"], row["Therapeutic"], scheme, ["human"], *REGIONS
            )
        )
    return result


class TestGetSeqViewDataframeReal:
    def test_imgt_scheme_formats_headers_and_orders_columns(self, therasabdab_subset_df):
        df = therasabdab_subset_df.head(2)
        result = _build_result(df, "imgt")

        seq_view_df = get_seq_view_dataframe(
            result, scheme="imgt", display_unformatted_residue_view=False
        )

        assert "seq_id" in seq_view_df.columns
        pos_cols = [c for c in seq_view_df.columns if c != "seq_id"]
        assert pos_cols, "expected at least one residue-position column"
        assert all(c.startswith("imgt_pos_") for c in pos_cols)
        assert len(seq_view_df) == len(df)

    def test_kabat_scheme_formats_headers(self, therasabdab_subset_df):
        df = therasabdab_subset_df.head(2)
        result = _build_result(df, "kabat")

        seq_view_df = get_seq_view_dataframe(
            result, scheme="kabat", display_unformatted_residue_view=False
        )

        pos_cols = [c for c in seq_view_df.columns if c != "seq_id"]
        assert pos_cols
        assert all(c.startswith("kabat_pos_") for c in pos_cols)

    def test_unformatted_residue_view_skips_header_formatting(self, therasabdab_subset_df):
        df = therasabdab_subset_df.head(2)
        result = _build_result(df, "kabat")

        seq_view_df = get_seq_view_dataframe(
            result, scheme="kabat", display_unformatted_residue_view=True
        )

        pos_cols = [c for c in seq_view_df.columns if c != "seq_id"]
        assert pos_cols
        assert not any(c.startswith("kabat_pos_") for c in pos_cols)
