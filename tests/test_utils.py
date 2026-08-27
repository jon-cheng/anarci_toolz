from multiprocessing import cpu_count

import pandas as pd
import pytest

from anarci_toolz.utils import create_row_id, get_num_cpu, slice_seqs_from_indices


class TestCreateRowId:
    def test_assigns_sequential_string_ids(self):
        df = pd.DataFrame({"sequence_aa": ["AAA", "BBB", "CCC"]})
        result_df, id_to_seq = create_row_id(df, "sequence_aa")

        assert result_df["seq_id"].tolist() == ["1", "2", "3"]
        assert result_df.columns[0] == "seq_id"
        assert id_to_seq == {"1": "AAA", "2": "BBB", "3": "CCC"}

    def test_empty_dataframe(self):
        df = pd.DataFrame({"sequence_aa": []})
        result_df, id_to_seq = create_row_id(df, "sequence_aa")

        assert result_df.columns[0] == "seq_id"
        assert result_df.empty
        assert id_to_seq == {}


class TestSliceSeqsFromIndices:
    def test_normal_slicing(self):
        df = pd.DataFrame(
            {
                "sequence_aa": ["ABCDEFGH", "IJKLMNOP"],
                "cdr1_aa_start": [1, 2],
                "cdr1_aa_end": [4, 5],
            }
        )
        result = slice_seqs_from_indices(
            df, region="cdr1", seq_dna_header=None, seq_aa_header="sequence_aa", mode="aa"
        )
        assert result["cdr1_aa"].tolist() == ["BCD", "KLM"]

    def test_nan_start_end_produces_empty_string(self):
        df = pd.DataFrame(
            {
                "sequence_aa": ["ABCDEFGH", "IJKLMNOP"],
                "cdr1_aa_start": pd.array([1, None], dtype="Int64"),
                "cdr1_aa_end": pd.array([4, None], dtype="Int64"),
            }
        )
        result = slice_seqs_from_indices(
            df, region="cdr1", seq_dna_header=None, seq_aa_header="sequence_aa", mode="aa"
        )
        assert result["cdr1_aa"].tolist() == ["BCD", ""]


class TestGetNumCpu:
    def test_returns_given_value(self):
        assert get_num_cpu(4) == 4

    def test_defaults_to_cpu_count(self):
        assert get_num_cpu(None) == cpu_count()
        assert get_num_cpu() == cpu_count()
