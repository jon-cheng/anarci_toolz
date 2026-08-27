import pandas as pd

from anarci_toolz.abnumber_tool import (
    ABNUMBER_TOOL_SCHEMA,
    ABNUMBER_TOOL_SCHEMA_DNA,
    drop_start_end,
    generate_dna_indices,
    reorder_abnumber_tool_result,
    unpack_data,
)


class TestGenerateDnaIndices:
    def test_renames_and_triples_aa_start_end_columns(self):
        df_ab = pd.DataFrame({"cdr1_aa_start": [1, 2], "cdr1_aa_end": [4, 5], "other": ["x", "y"]})
        result = generate_dna_indices(df_ab)

        assert result["cdr1_dna_start"].tolist() == [3, 6]
        assert result["cdr1_dna_end"].tolist() == [12, 15]
        # original aa columns are untouched
        assert result["cdr1_aa_start"].tolist() == [1, 2]

    def test_no_start_end_columns_is_noop(self):
        df_ab = pd.DataFrame({"other": ["x", "y"]})
        result = generate_dna_indices(df_ab)
        assert list(result.columns) == ["other"]


class TestDropStartEnd:
    def test_drops_start_and_end_columns(self):
        df = pd.DataFrame(
            {
                "cdr1_aa_start": [1],
                "cdr1_aa_end": [4],
                "cdr1_aa": ["ABC"],
                "seq_id": ["1"],
            }
        )
        result = drop_start_end(df)
        assert list(result.columns) == ["cdr1_aa", "seq_id"]

    def test_no_start_end_columns_keeps_all(self):
        df = pd.DataFrame({"seq_id": ["1"], "cdr1_aa": ["ABC"]})
        result = drop_start_end(df)
        assert list(result.columns) == ["seq_id", "cdr1_aa"]


def _full_schema_df(include_dna=False):
    cols = {c: ["x"] for c in ABNUMBER_TOOL_SCHEMA}
    if include_dna:
        cols.update({c: ["x"] for c in ABNUMBER_TOOL_SCHEMA_DNA})
    cols["seq_id"] = ["1"]
    return pd.DataFrame(cols)


class TestReorderAbnumberToolResult:
    def test_without_dna_header(self):
        df_ab = _full_schema_df(include_dna=False)
        result = reorder_abnumber_tool_result(df_ab, seq_dna_header=None)

        assert list(result.columns) == ["seq_id"] + ABNUMBER_TOOL_SCHEMA
        for dna_col in ABNUMBER_TOOL_SCHEMA_DNA:
            assert dna_col not in result.columns

    def test_with_dna_header(self):
        df_ab = _full_schema_df(include_dna=True)
        result = reorder_abnumber_tool_result(df_ab, seq_dna_header="sequence_dna")

        assert list(result.columns) == ["seq_id"] + ABNUMBER_TOOL_SCHEMA + ABNUMBER_TOOL_SCHEMA_DNA
        for dna_col in ABNUMBER_TOOL_SCHEMA_DNA:
            assert dna_col in result.columns


class TestUnpackData:
    def test_flattens_one_row_per_region(self):
        results = [
            [
                {
                    "seq1": {
                        "sequence_alignment_aa": "EVQL",
                        "passed_abnumber": True,
                        "scheme": "imgt",
                        "species": "human",
                        "chain_type": "H",
                        "v_gene": "IGHV3-66*01",
                        "j_gene": "IGHJ4*01",
                        "cdr1_seq": {"sequence": "GFT", "start_idx": 25, "end_idx": 33},
                        "cdr2_seq": {"sequence": "ISG", "start_idx": 50, "end_idx": 58},
                    }
                }
            ]
        ]
        df = unpack_data(results)

        assert len(df) == 2
        assert set(df["region"]) == {"cdr1_seq", "cdr2_seq"}
        assert df.loc[df["region"] == "cdr1_seq", "sequence"].iloc[0] == "GFT"
        assert df.loc[df["region"] == "cdr1_seq", "start_idx"].iloc[0] == 25
        assert (df["seq_id"] == "seq1").all()
        assert (df["v_gene"] == "IGHV3-66*01").all()

    def test_empty_results(self):
        df = unpack_data([])
        assert df.empty
        assert list(df.columns) == [
            "seq_id",
            "sequence_alignment_aa",
            "passed_abnumber",
            "scheme",
            "species",
            "chain_type",
            "v_gene",
            "j_gene",
            "region",
            "sequence",
            "start_idx",
            "end_idx",
        ]
