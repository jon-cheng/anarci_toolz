import re
import pandas as pd
from typing import List, Tuple, Dict


def extract_leading_integer(pattern: str) -> int | None:
    match = re.match(r"(\d+)", pattern)
    if match:
        return int(match.group(1))
    return None


def extract_trailing_alpha(pattern: str) -> str:
    match = re.search(r"([A-Z])$", pattern)
    if match:
        return match.group(1)
    else:
        return ""


def extract_index_tuple(pattern: str) -> List[Tuple[int, str]]:
    return (extract_leading_integer(pattern), extract_trailing_alpha(pattern))


def get_all_indices_with_and_without_suffixes(
    sequence_indices: List[str],
) -> List[Tuple[int, str]]:
    """_summary_

    Args:
        sequence_indices (List[str]): Raw sequence indices derived from AbNumber - list of strings of only numbers or numbers with suffixes representing the sequence indices of the scheme (IMGT) convention

    Returns:
        List[Tuple[int, str]]: _description_
    """
    all_indices_with_and_without_suffixes = [
        extract_index_tuple(item) for item in sequence_indices
    ]
    return all_indices_with_and_without_suffixes


def filter_for_suffixed_indices(
    all_indices_with_and_without_suffixes: List[Tuple[int, str]]
) -> List[Tuple[int, str]]:
    return [i for i in all_indices_with_and_without_suffixes if i[1] != ""]


def get_unique_integers(tuples_list: List[Tuple[int, str]]) -> List[int]:
    unique_integers = {t[0] for t in tuples_list}
    return sorted(list(unique_integers))


def is_empty_string_before_A(lst: List[str]) -> bool:
    try:
        index_empty = lst.index("")
        index_A = lst.index("A")
        return index_empty < index_A
    except ValueError:
        raise ValueError("Either '' or 'A' is not present in the position suffix list")


def assign_order(called_order: bool) -> str:
    if called_order:
        return "F"
    else:
        return "R"


def get_per_seq_ordered_lookup(
    unique_integers: List[int],
    all_indices_with_and_without_suffixes: List[Tuple[int, str]],
) -> Dict[int, str]:
    """Returns a dict of the order of the suffixes for each unique integer in the sequence indices"""
    per_seq_order_lookup = {}
    for u in unique_integers:
        per_pos_indices = [
            i[1] for i in all_indices_with_and_without_suffixes if i[0] == u
        ]
        per_seq_order_lookup[u] = assign_order(
            is_empty_string_before_A(per_pos_indices)
        )
    return per_seq_order_lookup


def backfill_per_seq_order_lookup(
    all_indices_with_and_without_suffixes: List[int], per_seq_order_lookup: dict
) -> dict:
    sorted_all_indices = sorted(
        list(set([i[0] for i in all_indices_with_and_without_suffixes]))
    )
    for i in sorted_all_indices:
        if i not in per_seq_order_lookup:
            per_seq_order_lookup[i] = "N"
    return per_seq_order_lookup


def get_per_seq_order_lookup(sequence_indices: List[str]) -> Dict[int, str]:
    all_indices_with_and_without_suffixes = get_all_indices_with_and_without_suffixes(
        sequence_indices
    )
    suffixed_indices = filter_for_suffixed_indices(
        all_indices_with_and_without_suffixes
    )
    unique_suffixed_integers = get_unique_integers(suffixed_indices)
    per_seq_order_lookup = get_per_seq_ordered_lookup(
        unique_suffixed_integers, all_indices_with_and_without_suffixes
    )
    # per_seq_order_lookup = backfill_per_seq_order_lookup(all_indices_with_and_without_suffixes, per_seq_order_lookup)
    return per_seq_order_lookup


class MergeConflictError(Exception):
    def __init__(self, message, conflicts):
        super().__init__(message)
        self.conflicts = conflicts


def merge_dicts_with_conflict_check(dicts_list: List[dict]) -> dict:
    merged_dict = {}
    conflicts = {}

    for idx, du in enumerate(dicts_list):
        for key, value in du.items():
            if key in merged_dict and merged_dict[key] != value:
                if key not in conflicts:
                    conflicts[key] = []
                conflicts[key].append(idx)
            else:
                merged_dict[key] = value

    if conflicts:
        conflict_keys = [str(i) for i in list(conflicts.keys())]
        conflict_message = f"Conflicts found in positions: {', '.join(conflict_keys)}"
        raise MergeConflictError(conflict_message, conflicts)

    return merged_dict


def get_comprehensive_dicts_list(lists: List[List[str]]) -> List[Dict[int, str]]:
    ls = []
    for sequence_indices in lists:
        ls.append(get_per_seq_order_lookup(sequence_indices))
    return ls


def get_comprehensive_order(lists: List[List[str]]) -> Dict[int, str]:
    return merge_dicts_with_conflict_check(get_comprehensive_dicts_list(lists))


def get_abnumber_called_indices(
    sequence_alignment_aa_dataframes: List[pd.DataFrame],
) -> List[List[str]]:
    """Generates a list of lists of index positions from the AbNumber sequence dataframes

    Args:
        sequence_alignment_aa_dataframes (List[pd.DataFrame]): _description_

    Returns:
        List[List[str]]: _description_
    """
    return [
        df.drop(
            columns=[
                "species",
                "chain_type",
                "seq_id",
            ],
            errors="ignore",
        ).columns.tolist()
        for df in sequence_alignment_aa_dataframes
        if df is not None
    ]
