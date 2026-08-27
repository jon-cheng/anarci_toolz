import pandas as pd
import pytest

from anarci_toolz.numbering import (
    format_header_style,
    get_zero_pad_pos,
    sort_numbers,
    zero_pad_int_3d,
)


class TestSortNumbers:
    def test_kabat_sorts_numerically_then_alphabetically(self):
        result = sort_numbers(["10", "2", "2A", "1"], scheme="kabat")
        assert result == ["1", "2", "2A", "10"]

    def test_imgt_requires_order_dict(self):
        with pytest.raises(ValueError):
            sort_numbers(["1", "2"], scheme="imgt")

    def test_imgt_forward_order(self):
        result = sort_numbers(["2A", "2"], scheme="imgt", order_dict={2: "F"})
        assert result == ["2", "2A"]

    def test_imgt_reverse_order(self):
        result = sort_numbers(["2", "2A"], scheme="imgt", order_dict={2: "R"})
        assert result == ["2A", "2"]

    def test_unsupported_scheme_raises(self):
        with pytest.raises(ValueError):
            sort_numbers(["1", "2"], scheme="chothia")


class TestZeroPadInt3d:
    def test_pads_to_three_digits(self):
        assert zero_pad_int_3d(5) == "005"

    def test_no_padding_needed(self):
        assert zero_pad_int_3d(123) == "123"


class TestGetZeroPadPos:
    def test_with_suffix(self):
        assert get_zero_pad_pos("5A") == "005A"

    def test_without_suffix(self):
        assert get_zero_pad_pos("12") == "012"


class TestFormatHeaderStyle:
    def test_renames_position_columns_and_keeps_seq_id(self):
        df = pd.DataFrame({"seq_id": ["1", "2"], "5": ["A", "B"], "12A": ["C", "D"]})
        result = format_header_style(df, scheme="imgt")

        assert "seq_id" in result.columns
        assert "imgt_pos_005" in result.columns
        assert "imgt_pos_012A" in result.columns
        assert result["seq_id"].tolist() == ["1", "2"]
        assert result["imgt_pos_005"].tolist() == ["A", "B"]
