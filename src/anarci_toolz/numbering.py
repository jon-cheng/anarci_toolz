import re
import json
import pandas as pd
from typing import List, Dict, Optional
from collections import OrderedDict
from typing import List, Optional, Dict
from .imgt_index_ordering import (
    extract_leading_integer,
    extract_trailing_alpha,
    get_abnumber_called_indices,
    get_comprehensive_order,
)

# Helper functions for numbering schemes and sequence dataframes


def sort_numbers(
    numbers: List[str], scheme: str, order_dict: Optional[Dict[int, str]] = None
) -> List[str]:
    """Sorts numbering strings into order based on a custom order dictionary and scheme.

    Args:
        numbers (List[str]): List of numbering strings.
        scheme (str): The scheme to use for sorting ('kabat' or 'imgt').
        order_dict (Dict[int, str], optional): Dictionary specifying 'forward' or 'reverse' order for each integer position.

    Returns:
        List[str]: Sorted list of numbering strings.
    """

    def kabat_key(s):
        # Split the string into numeric and alphabetic parts
        match = re.match(r"(\d+)([A-Za-z]?)", s)
        if match:
            number = int(match.group(1))
            suffix = match.group(2)
            return (number, suffix)
        return (float("inf"), "")  # Handle unexpected format

    def imgt_key(x):
        if order_dict is not None:
            number, suffix = kabat_key(x)
            if order_dict.get(number) == "R":
                return (number, -ord(suffix) if suffix else float("inf"))
            # For forward ordering, force empty strings to precede alphabetic characters
            return (number, ord(suffix) if suffix else -1)
        else:
            raise ValueError(
                "Order dictionary must be provided for IMGT numbering scheme."
            )

    if scheme == "kabat":
        return sorted(numbers, key=kabat_key)

    elif scheme == "imgt":
        return sorted(numbers, key=imgt_key)

    else:
        raise ValueError(
            "Unsupported scheme for column ordering. Please specify either 'kabat' or 'imgt' scheme."
        )


def convert_abnumber_positions_to_json_string(chain):
    return json.dumps(OrderedDict((str(k), v) for k, v in chain.positions.items()))


def get_seq_alignment_dataframes(result):
    seq_ids = []
    sequence_alignment_aa_dataframes = []

    for r in result:
        sequence_dict = r[0]
        for seq_id, sequence_info in sequence_dict.items():
            seq_ids.append(seq_id)
            sequence_alignment_aa_dataframes.append(
                sequence_info["sequence_alignment_aa_dataframe"]
            )
    return seq_ids, sequence_alignment_aa_dataframes


def concatenate_seq_dataframes(seq_ids, sequence_alignment_aa_dataframes):

    concat_dfs = []
    for seq_id, df in zip(seq_ids, sequence_alignment_aa_dataframes):
        if df is not None:
            # Add seq_id as a column to the DataFrame
            df["seq_id"] = seq_id
            concat_dfs.append(df)
        else:
            # Create a DataFrame with seq_id only
            concat_dfs.append(pd.DataFrame({"seq_id": [seq_id]}))

    # Concatenate all DataFrames
    result_df = pd.concat(concat_dfs, ignore_index=True)
    return result_df


def clean_seq_dataframe(
    result_df,
    scheme,
    display_unformatted_residue_view,
    comprehensive_order,
):
    result_df = result_df.drop(["species", "chain_type"], axis=1)
    # Reorder columns to retain the order of appearance, use custom sorting order derived from the dataset
    df_f = result_df[
        sort_numbers(result_df.columns.tolist(), scheme, comprehensive_order)
    ]

    # Format headers to Numbering Standard Style
    if not display_unformatted_residue_view:
        df_f = format_header_style(df_f, scheme)
    return df_f


def zero_pad_int_3d(number: int) -> str:
    return f"{number:03d}"


def get_zero_pad_pos(pattern: str) -> str:
    return zero_pad_int_3d(extract_leading_integer(pattern)) + extract_trailing_alpha(
        pattern
    )


def format_header_style(df: pd.DataFrame, scheme: str) -> pd.DataFrame:
    pos_cols = [i for i in df.columns.tolist() if i != "seq_id"]
    new_pos_cols = [(f"{scheme}" + "_pos_" + get_zero_pad_pos(i)) for i in pos_cols]
    df = df.rename(columns=dict(zip(pos_cols, new_pos_cols)))
    return df


def get_seq_view_dataframe(result, scheme, display_unformatted_residue_view):
    seq_ids, sequence_alignment_aa_dataframes = get_seq_alignment_dataframes(result)

    # IMGT numbering logic
    if scheme == "imgt":
        sequence_indices = get_abnumber_called_indices(sequence_alignment_aa_dataframes)
        comprehensive_order = get_comprehensive_order(sequence_indices)
    else:
        comprehensive_order = None

    result_df = concatenate_seq_dataframes(seq_ids, sequence_alignment_aa_dataframes)
    result_df = clean_seq_dataframe(
        result_df,
        scheme,
        display_unformatted_residue_view,
        comprehensive_order=comprehensive_order,
    )
    return result_df
