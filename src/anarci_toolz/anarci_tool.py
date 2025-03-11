import pandas as pd
from anarci import anarci
from typing import List, Dict
from tqdm import tqdm
from loguru import logger
from multiprocessing import Pool, cpu_count


ANARCI_RESULTS_COLS = [
    "passed_anarci",
    "variable_region_start_index",
    "variable_region_end_index",
    "chain_type",
    "v_call",
    "j_call",
    "e_value",
    "bitscore",
    "bias",
]

ANARCI_RESULTS_SCHEMA = [
    "passed_anarci",
    "variable_region_start_index",
    "variable_region_end_index",
    "e_value",
    "bitscore",
    "bias",
]


def get_alignment_details(
    seq: list, name: str, scheme: list, allowed_species: List[str]
) -> list:

    numbering, alignment_details, hit_tables = anarci(
        [(name, seq)],
        scheme=scheme,
        assign_germline=True,
        allowed_species=allowed_species,
        output=False,
    )
    return alignment_details


def extract_useful_info(alignment_details: list):
    variable_region_start_index = alignment_details[0][0]["query_start"]
    variable_region_end_index = alignment_details[0][0]["query_end"]
    chain_type = alignment_details[0][0]["chain_type"]
    v_call = alignment_details[0][0]["germlines"]["v_gene"][0][1]
    j_call = alignment_details[0][0]["germlines"]["j_gene"][0][1]
    e_value = alignment_details[0][0]["evalue"]
    bitscore = alignment_details[0][0]["bitscore"]
    bias = alignment_details[0][0]["bias"]
    return (
        variable_region_start_index,
        variable_region_end_index,
        chain_type,
        v_call,
        j_call,
        e_value,
        bitscore,
        bias,
    )


def get_anarci_alignment(
    seq: str, name: str, scheme: str, allowed_species: List[str]
) -> tuple:
    try:
        alignment_details = get_alignment_details(seq, name, scheme, allowed_species)
        results_tuple = extract_useful_info(alignment_details)
        passed_anarci = True
        results_tuple = (name,) + (passed_anarci,) + results_tuple

    except TypeError as e:
        results_tuple = (None,) * (len(ANARCI_RESULTS_COLS) - 1)
        passed_anarci = False
        results_tuple = (name,) + (passed_anarci,) + results_tuple

    return results_tuple


def run_parallel_anarci(
    df: pd.DataFrame,
    seqs: Dict[str, str],
    scheme: str,
    allowed_species: List[str],
    seq_aa_header: str,
) -> pd.DataFrame:

    logger.info("Starting ANARCI processing ...")

    assert (
        seq_aa_header in df.columns
    ), f"{seq_aa_header} is not a column in the DataFrame, please rename your amino acid input sequence column to {seq_aa_header}"

    tasks = [(seq, seq_id, scheme, allowed_species) for seq_id, seq in seqs.items()]

    with Pool(processes=cpu_count()) as p:
        results = p.starmap(
            get_anarci_alignment,
            tqdm(tasks, total=len(tasks), unit="seq", desc="Processing sequences"),
        )
        logger.success("ANARCI parallel processing complete")

    df_anarci_result = pd.DataFrame(
        data=results, columns=["seq_id"] + ANARCI_RESULTS_COLS
    )

    df_anarci_result["variable_region_start_index"] = df_anarci_result[
        "variable_region_start_index"
    ].astype("Int64")
    df_anarci_result["variable_region_end_index"] = df_anarci_result[
        "variable_region_end_index"
    ].astype("Int64")

    df_anarci_result = df_anarci_result[["seq_id"] + ANARCI_RESULTS_SCHEMA]

    return df_anarci_result
