import pandas as pd
import pytest

from anarci_toolz.validation import (
    Validator,
    disallow_existing_columns,
    disallow_existing_residue_view_columns,
)
from anarci_toolz.abnumber_tool import ABNUMBER_TOOL_SCHEMA
from anarci_toolz.anarci_tool import ANARCI_RESULTS_SCHEMA


class TestDisallowExistingColumns:
    def test_passes_on_clean_columns(self):
        df = pd.DataFrame({"sequence_aa": ["AAA"], "name": ["ab1"]})
        assert disallow_existing_columns(df) is True

    def test_raises_on_schema_column(self):
        df = pd.DataFrame({"sequence_aa": ["AAA"], ABNUMBER_TOOL_SCHEMA[0]: ["x"]})
        with pytest.raises(ValueError):
            disallow_existing_columns(df)


class TestDisallowExistingResidueViewColumns:
    @pytest.mark.parametrize("col", ["123", "123A", "imgt_pos_5", "imgt_pos_5A"])
    def test_raises_on_matching_columns(self, col):
        df = pd.DataFrame({"sequence_aa": ["AAA"], col: ["x"]})
        with pytest.raises(ValueError):
            disallow_existing_residue_view_columns(df)

    @pytest.mark.parametrize("col", ["sequence_aa", "name", "abc", "pos_5"])
    def test_passes_on_non_matching_columns(self, col):
        df = pd.DataFrame({"sequence_aa": ["AAA"], col: ["x"]})
        assert disallow_existing_residue_view_columns(df) is True


class TestValidator:
    def test_passes_on_clean_dataframe(self):
        df = pd.DataFrame({"sequence_aa": ["AAA", "BBB"]})
        validator = Validator(seq_aa_header="sequence_aa", display_residue_view=False)
        assert validator.validate(df) is True

    @pytest.mark.parametrize(
        "bad_col",
        [ABNUMBER_TOOL_SCHEMA[0], ANARCI_RESULTS_SCHEMA[0]],
    )
    def test_raises_on_schema_column(self, bad_col):
        df = pd.DataFrame({"sequence_aa": ["AAA"], bad_col: ["x"]})
        validator = Validator(seq_aa_header="sequence_aa", display_residue_view=False)
        with pytest.raises(Exception):
            validator.validate(df)

    def test_raises_on_residue_view_column_when_enabled(self):
        df = pd.DataFrame({"sequence_aa": ["AAA"], "123A": ["x"]})
        validator = Validator(seq_aa_header="sequence_aa", display_residue_view=True)
        with pytest.raises(Exception):
            validator.validate(df)

    def test_ignores_residue_view_column_when_disabled(self):
        df = pd.DataFrame({"sequence_aa": ["AAA"], "123A": ["x"]})
        validator = Validator(seq_aa_header="sequence_aa", display_residue_view=False)
        assert validator.validate(df) is True
