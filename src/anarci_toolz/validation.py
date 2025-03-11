import re
import pandas as pd
import pandera as pa
from pandera import Column, DataFrameSchema, Check
from loguru import logger
from typing import List
from .abnumber_tool import ABNUMBER_TOOL_SCHEMA, ABNUMBER_TOOL_SCHEMA_DNA
from .anarci_tool import ANARCI_RESULTS_SCHEMA


# Define the schema
def disallow_existing_columns(df: pd.DataFrame) -> bool:
    detected_disallowed_columns = [
        col
        for col in df.columns
        if col
        in ABNUMBER_TOOL_SCHEMA + ABNUMBER_TOOL_SCHEMA_DNA + ANARCI_RESULTS_SCHEMA
    ]
    if detected_disallowed_columns:
        detected_columns_str = ", ".join(detected_disallowed_columns)
        raise ValueError(
            f"Columns: {detected_columns_str} already exist in the DataFrame. This is not allowed because they are part of the output schema. Consider dropping these columns before running the pipeline."
        )
    return True


def disallow_existing_residue_view_columns(df: pd.DataFrame) -> bool:
    # Define the patterns to disallow
    patterns = [
        r"^\d+$",  # <int>
        r"^\d+[A-Za-z]$",  # <int><alphabet>
        r"^[A-Za-z]+_pos_\d+$",  # <str>_pos_<int>
        r"^[A-Za-z]+_pos_\d+[A-Za-z]$",  # <str>_pos_<int><alphabet>
    ]

    # Check if any column name matches any of the patterns
    for col in df.columns:
        if any(re.match(pattern, col) for pattern in patterns):
            raise ValueError(
                f"Columns that denote residue positions exist in the DataFrame. This is not allowed because it is part of the output schema. Consider dropping these columns before running the pipeline."
            )
    return True


class Validator:
    def __init__(self, seq_aa_header: str, display_residue_view: bool):
        self.seq_aa_header = seq_aa_header
        self.display_residue_view = display_residue_view

    def validate(self, df: pd.DataFrame) -> bool:
        schema = DataFrameSchema(
            {self.seq_aa_header: Column(pa.String)},
            checks=[
                Check(
                    disallow_existing_columns,
                    error="Disallowed columns already exist in the DataFrame.",
                ),
                Check(
                    disallow_existing_residue_view_columns,
                    error="Columns that denote residue positions exist in the DataFrame.",
                ),
            ],
        )

        if not self.display_residue_view:
            # Only run the disallow_existing_columns check
            schema = DataFrameSchema(
                {self.seq_aa_header: Column(pa.String)},
                checks=[
                    Check(
                        disallow_existing_columns,
                        error="Disallowed columns already exist in the DataFrame.",
                    )
                ],
            )

        # Validate the DataFrame
        schema.validate(df)
        return True
