import pandas as pd
from typing import Tuple, Dict
from multiprocessing import cpu_count


def create_row_id(
    df: pd.DataFrame, seq_aa_header: str
) -> Tuple[pd.DataFrame, Dict[str, str]]:
    """Assign an internal ID column to input dataframe, enumerating all rows, not necessarily all unique seqs"""
    df = df.reset_index(drop=True)
    df["seq_id"] = df.index + 1
    df["seq_id"] = df["seq_id"].astype(str)
    id_to_seq = dict(zip(df["seq_id"], df[seq_aa_header]))

    # rearrange "seq_id" to be the first column
    df = df.reindex(columns=["seq_id"] + [col for col in df.columns if col != "seq_id"])

    return df, id_to_seq


def slice_seqs_from_indices(
    df: pd.DataFrame,
    region: str,
    seq_dna_header: str,
    seq_aa_header: str,
    mode: str = "dna",
):
    if mode == "dna":
        header = seq_dna_header
    else:
        header = seq_aa_header

    df[f"{region}_{mode}"] = df.apply(
        lambda row: (
            row[header][row[f"{region}_{mode}_start"] : row[f"{region}_{mode}_end"]]
            if pd.notna(row[f"{region}_{mode}_start"])
            and pd.notna(row[f"{region}_{mode}_end"])
            else ""
        ),
        axis=1,
    )

    return df


def get_num_cpu(num_cpu: int = None) -> int:
    """Get number of CPUs to use for parallel processing

    :param num_cpu: Number of CPUs to use, defaults to None
    :type num_cpu: int, optional
    :return: Number of CPUs to use
    :rtype: int
    """
    processes = num_cpu if num_cpu is not None else cpu_count()
    return processes
