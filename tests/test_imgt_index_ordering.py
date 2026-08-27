import pandas as pd
import pytest

from anarci_toolz.imgt_index_ordering import (
    MergeConflictError,
    assign_order,
    extract_leading_integer,
    extract_trailing_alpha,
    get_abnumber_called_indices,
    get_comprehensive_dicts_list,
    get_comprehensive_order,
    get_per_seq_order_lookup,
    is_empty_string_before_A,
    merge_dicts_with_conflict_check,
)

# `sort_numbers` is actually defined in `numbering.py` (not this module,
# despite superficially sitting next to IMGT-ordering logic) — see
# TestSortNumbers in test_numbering.py for its coverage.


class TestExtractLeadingInteger:
    def test_extracts_leading_integer(self):
        assert extract_leading_integer("123A") == 123

    def test_no_match_returns_none(self):
        assert extract_leading_integer("A") is None


class TestExtractTrailingAlpha:
    def test_extracts_trailing_alpha(self):
        assert extract_trailing_alpha("123A") == "A"

    def test_no_match_returns_empty_string(self):
        assert extract_trailing_alpha("123") == ""


class TestIsEmptyStringBeforeA:
    def test_true_when_empty_before_a(self):
        assert is_empty_string_before_A(["", "A"]) is True

    def test_false_when_a_before_empty(self):
        assert is_empty_string_before_A(["A", ""]) is False

    def test_raises_when_neither_present(self):
        with pytest.raises(ValueError):
            is_empty_string_before_A(["B", "C"])


class TestAssignOrder:
    def test_true_maps_to_f(self):
        assert assign_order(True) == "F"

    def test_false_maps_to_r(self):
        assert assign_order(False) == "R"


class TestGetPerSeqOrderLookup:
    def test_forward_order(self):
        result = get_per_seq_order_lookup(["1", "2", "2A", "3"])
        assert result == {2: "F"}

    def test_reverse_order(self):
        result = get_per_seq_order_lookup(["2A", "2"])
        assert result == {2: "R"}

    def test_no_suffixed_indices_returns_empty_dict(self):
        assert get_per_seq_order_lookup(["1", "2", "3"]) == {}


class TestMergeDictsWithConflictCheck:
    def test_no_conflict_merge(self):
        result = merge_dicts_with_conflict_check([{1: "F"}, {2: "R"}])
        assert result == {1: "F", 2: "R"}

    def test_conflict_raises(self):
        with pytest.raises(MergeConflictError) as exc_info:
            merge_dicts_with_conflict_check([{1: "F"}, {1: "R"}])
        assert exc_info.value.conflicts == {1: [1]}


class TestGetComprehensive:
    def test_get_comprehensive_dicts_list(self):
        lists = [["1", "2", "2A"], ["3A", "3"]]
        result = get_comprehensive_dicts_list(lists)
        assert result == [{2: "F"}, {3: "R"}]

    def test_get_comprehensive_order(self):
        lists = [["1", "2", "2A"], ["3A", "3"]]
        result = get_comprehensive_order(lists)
        assert result == {2: "F", 3: "R"}

    def test_get_comprehensive_order_conflict_raises(self):
        lists = [["2", "2A"], ["2A", "2"]]
        with pytest.raises(MergeConflictError):
            get_comprehensive_order(lists)


class TestGetAbnumberCalledIndices:
    def test_drops_metadata_columns_and_skips_none(self):
        df1 = pd.DataFrame(
            {"species": ["human"], "chain_type": ["H"], "seq_id": ["1"], "5": ["A"]}
        )
        df2 = None
        df3 = pd.DataFrame({"seq_id": ["2"], "6": ["B"], "6A": ["C"]})
        result = get_abnumber_called_indices([df1, df2, df3])
        assert result == [["5"], ["6", "6A"]]
