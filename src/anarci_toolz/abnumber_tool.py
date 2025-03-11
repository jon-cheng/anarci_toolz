from abnumber import Chain, Position
from .utils import slice_seqs_from_indices, create_row_id, get_num_cpu
from .numbering import get_seq_view_dataframe
from typing import List, Tuple, Dict, Union, Optional
from tqdm import tqdm
from loguru import logger
import pandas as pd
from multiprocessing import Pool


ABNUMBER_TOOL_SCHEMA = [
    "scheme",  # input scheme
    "passed_abnumber",  # whether AbNumber call was successful
    "sequence_alignment_aa",  # extracted variable region sequence determined by AbNumber
    "species",  # species determined by AbNumber
    "chain_type",  # one of H, K, L
    "v_gene",  # v-gene call
    "j_gene",  # j-gene call
    "cdr1_aa",
    "cdr2_aa",
    "cdr3_aa",
    "fr1_aa",
    "fr2_aa",
    "fr3_aa",
    "fr4_aa",
]

ABNUMBER_TOOL_SCHEMA_DNA = [
    "cdr1_dna",
    "cdr2_dna",
    "cdr3_dna",
    "fr1_dna",
    "fr2_dna",
    "fr3_dna",
    "fr4_dna",
]

REGIONS = (
    "fr1_seq",
    "cdr1_seq",
    "fr2_seq",
    "cdr2_seq",
    "fr3_seq",
    "cdr3_seq",
    "fr4_seq",
)


def run_parallel_abnumber(
    df: pd.DataFrame,
    seqs: Dict[str, str],
    scheme: str,
    seq_aa_header: str,
    seq_dna_header: str,
    allowed_species: List[str],
    retain_indices: bool,
    display_residue_view: bool,
    display_unformatted_residue_view: bool,
    num_cpu: int,
) -> pd.DataFrame:

    logger.info("Starting AbNumber processing ...")

    # parser for allowed_species input; AbNumber can handle either str or list for allowed_species
    if (allowed_species is not None) and (len(allowed_species) == 1):
        allowed_species = allowed_species[0]

    assert (
        seq_aa_header in df.columns
    ), f"{seq_aa_header} is not a column in the DataFrame, please rename your amino acid input sequence column to {seq_aa_header}"

    tasks = [
        (seq, seq_id, scheme, allowed_species, *REGIONS) for seq_id, seq in seqs.items()
    ]
    processes = get_num_cpu(num_cpu)
    result = parallel_get_region_seqs(tasks, processes)
    df_ab_long = unpack_data(result)

    df_ab_long["region"] = df_ab_long["region"].apply(lambda x: x.split("_seq")[0])
    df_ab_long = df_ab_long.rename(
        columns={"sequence": "aa", "start_idx": "aa_start", "end_idx": "aa_end"}
    )
    df_ab_long["aa_start"] = df_ab_long["aa_start"].astype("Int64")
    df_ab_long["aa_end"] = df_ab_long["aa_end"].astype("Int64")

    df_ab = df_ab_long.pivot_table(
        index="seq_id",
        columns="region",
        values=["aa", "aa_start", "aa_end"],
        aggfunc="first",
        dropna=False,
    )

    # collapse multiindex column names
    df_ab.columns = df_ab.columns.map(lambda x: "_".join(x[::-1]))
    df_ab = df_ab.reset_index()

    # get other data on per-sequence level
    df_ab_vseq = df_ab_long.pivot_table(
        index="seq_id",
        values=[
            "sequence_alignment_aa",
            "passed_abnumber",
            "scheme",
            "species",
            "chain_type",
            "v_gene",
            "j_gene",
        ],
        aggfunc="first",
        dropna=False,
    ).reset_index(drop=False)

    df_ab = df_ab.merge(df_ab_vseq, on="seq_id", how="left")

    # handle DNA slicing
    df_ab = process_data(df, df_ab, seq_dna_header, seq_aa_header, retain_indices)
    df_ab = reorder_abnumber_tool_result(df_ab, seq_dna_header)

    # Get wide-display residue visualization
    if display_residue_view:
        df_seqs_wide = get_seq_view_dataframe(
            result, scheme, display_unformatted_residue_view
        )
        df_ab = df_ab.merge(df_seqs_wide, on="seq_id", how="left")

    return df_ab


def parallel_get_region_seqs(tasks: List[Tuple], processes: int) -> List:
    """Parallelizes AbNumber calls for multiple sequences"""

    with Pool(processes=processes) as p:
        results = p.starmap(
            get_region_seqs,
            tqdm(tasks, total=len(tasks), unit="seq", desc="Processing sequences"),
        )
        logger.success("AbNumber parallel processing complete")
    return results


def get_region_seqs(
    seq: str,
    name: str,
    scheme: str,
    allowed_species: Union[str, List[str]],
    *regions,
) -> List[dict]:
    """Get AbNumber Chain sub-sequence attributes for multiple regions, such as CDR and FWR sequences, with optional indices.

    Args:
        seq (str): The sequence to search within.
        scheme (str): The numbering scheme to use: "imgt", "kabat"
        *regions (str): Variable length list of regions, e.g., 'cdr1_seq', 'fr2_seq'.

    Returns:
        List of dictionaries with region sequence and, optionally, its start and end indices. Returns None for each part if not found.
    """
    results = []
    try:
        chain = create_chain(seq, name, scheme, allowed_species=allowed_species)
        pos_index_dict = get_pos_index_dict(chain)
        v_gene, j_gene = get_gene_calls(chain)

        # result for one sequence
        result_dict = {
            name: {
                "sequence_alignment_aa": str(chain),
                "sequence_alignment_aa_dataframe": chain.to_dataframe([chain]),
                "passed_abnumber": True,
                "scheme": scheme,
                "species": chain.species,
                "chain_type": chain.chain_type,
                "v_gene": v_gene,
                "j_gene": j_gene,
            }
        }
        for region in regions:
            region_seq = getattr(
                chain, region, None
            )  # Use getattr to dynamically access the attribute, with a fallback to None
            start_idx, end_idx = get_region_indices(chain, pos_index_dict, region)
            result_dict[name][region] = {
                "sequence": region_seq,
                "start_idx": start_idx,
                "end_idx": end_idx,
            }
        results.append(result_dict)

    except Exception as e:
        result_dict = {
            name: {
                "sequence_alignment_aa": None,
                "sequence_alignment_aa_dataframe": None,
                "passed_abnumber": False,
                "scheme": scheme,
                "species": None,
                "chain_type": None,
                "v_gene": None,
                "j_gene": None,
            }
        }
        for region in regions:
            result_dict[name][region] = {
                "sequence": None,
                "start_idx": None,
                "end_idx": None,
            }
        results.append(result_dict)

    return results


def create_chain(
    seq: str,
    name: str,
    scheme: str,
    allowed_species: Union[str, List[str]],
    assign_germline: bool = True,
):
    """Create a Chain object with the given sequence, scheme, and species.

    Args:
        seq (str): The sequence to search within.
        scheme (str): The numbering scheme to use: "imgt", "kabat"
        species (str): The species to use.

    Returns:
        Chain: The created Chain object.
    """
    chain = Chain(
        seq,
        name=name,
        scheme=scheme,
        allowed_species=allowed_species,
        assign_germline=assign_germline,
    )
    return chain


def get_gene_calls(chain: Chain) -> Tuple[str, str]:
    """Note: ANARCI does not support D gene calls, only V and J gene calls"""
    return chain.v_gene, chain.j_gene


def get_pos_index_dict(chain: Chain) -> Dict[Position, int]:
    """Generate mapping of antibody, e.g. kabat, imgt, etc. numbering position to linear index"""
    pos_ls = list(chain.positions.keys())

    pos_index_dict = {}
    for idx, k in enumerate(pos_ls):
        pos_index_dict[k] = idx

    return pos_index_dict


def get_region_indices(
    chain: Chain, pos_index_dict: Dict[Position, int], region: str
) -> Tuple[int, int]:
    """Generate start and end linear indices for a given region"""
    region = region.split("_seq")[0].upper()
    region_positions = chain.regions[region]

    # If not an empty ordered dict:
    if region_positions:
        first_key = next(iter(region_positions))
        last_key = next(reversed(region_positions))
        start_idx = pos_index_dict[first_key]
        end_idx = pos_index_dict[last_key]

        # end_idx + 1 to include the last index, as ANARCI convention is to include the last index
        return start_idx, end_idx + 1

    else:
        return None, None


def unpack_data(results: List) -> pd.DataFrame:

    unpacked_data = []

    for result in results:
        sequence_dict = result[0]
        for seq_id, sequence_info in sequence_dict.items():
            sequence_alignment_aa = sequence_info["sequence_alignment_aa"]
            passed_abnumber = sequence_info["passed_abnumber"]
            scheme = sequence_info["scheme"]
            species = sequence_info["species"]
            chain_type = sequence_info["chain_type"]
            v_gene = sequence_info["v_gene"]
            j_gene = sequence_info["j_gene"]
            for sequence_type, sequence_details in sequence_info.items():
                if isinstance(sequence_details, dict):
                    row = [
                        seq_id,
                        sequence_alignment_aa,
                        passed_abnumber,
                        scheme,
                        species,
                        chain_type,
                        v_gene,
                        j_gene,
                        sequence_type,
                        sequence_details.get("sequence", None),
                        sequence_details.get("start_idx", None),
                        sequence_details.get("end_idx", None),
                    ]
                    unpacked_data.append(row)

    df = pd.DataFrame(
        unpacked_data,
        columns=[
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
        ],
    )
    return df


def generate_dna_indices(df_ab):
    """Given AA indices, generate DNA indices by multiplying by 3
    :param df_ab: _description_
    :type df_ab: _type_
    :return: _description_
    :rtype: _type_
    """
    rename_dict = {
        col: col.replace("_aa", "_dna")
        for col in df_ab.columns
        if col.endswith("start") or col.endswith("end")
    }

    # Create new columns with the renamed column names
    for old_name, new_name in rename_dict.items():
        df_ab[new_name] = df_ab[old_name] * 3

    return df_ab


def get_all_region_slices(df, seq_aa_header, seq_dna_header):

    for region in REGIONS:
        region = region.split("_seq")[0]
        df_slice_result = slice_seqs_from_indices(
            df=df,
            region=region,
            seq_dna_header=seq_dna_header,
            seq_aa_header=seq_aa_header,
            mode="dna",
        )
    return df_slice_result


def drop_start_end(df):
    return df.drop(
        [col for col in df.columns if col.endswith("start") or col.endswith("end")],
        axis=1,
    )


def merge_aa_and_dna(df, df_ab, seq_aa_header, seq_dna_header):
    df_ab = generate_dna_indices(df_ab)
    df_merge = pd.merge(df, df_ab, on="seq_id", how="left")
    df_merge = get_all_region_slices(df_merge, seq_aa_header, seq_dna_header)
    return df_merge


def merge_aa_only(df, df_ab):
    return pd.merge(df, df_ab, on="seq_id", how="left")


def process_data(df, df_ab, seq_dna_header, seq_aa_header, retain_indices):
    """Processing logic after AbNumbering, merges to input DataFrame, slices sequences, drops indices"""
    if seq_dna_header is not None:
        df_merge = merge_aa_and_dna(df, df_ab, seq_aa_header, seq_dna_header)
    else:
        df_merge = merge_aa_only(df, df_ab)

    if not retain_indices:
        df_merge = drop_start_end(df_merge)
    return df_merge


def reorder_abnumber_tool_result(df_ab, seq_dna_header):

    if seq_dna_header is not None:
        ordered_cols = ABNUMBER_TOOL_SCHEMA + ABNUMBER_TOOL_SCHEMA_DNA
    else:
        ordered_cols = ABNUMBER_TOOL_SCHEMA

    # Get the columns in df but not in ordered_cols
    other_cols = [col for col in df_ab.columns if col not in ordered_cols]

    # Append ordered_cols to other_cols
    final_cols = other_cols + ordered_cols

    # Reorder the DataFrame
    df_ab = df_ab[final_cols]

    return df_ab
