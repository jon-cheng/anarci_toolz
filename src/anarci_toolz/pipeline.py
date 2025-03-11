import os
import time
import pandas as pd
from loguru import logger
from Bio.Seq import Seq
from . import anarci_tool as an
from . import abnumber_tool as ab
from .utils import create_row_id
from typing import List, Optional
from .validation import Validator


def file_or_dir_input(Input_Path_List):
    Cleaned_IPL = []
    for IP in Input_Path_List:
        print(IP)
        if os.path.isdir(IP):
            Input_Files = [
                os.path.join(IP, filename)
                for filename in os.listdir(IP)
                if filename.endswith(".csv")
            ]
            Cleaned_IPL.extend(Input_Files)
        else:
            Cleaned_IPL.append(IP)
    print(Cleaned_IPL)
    return Cleaned_IPL


def run_anarci_toolz(
    df: pd.DataFrame,
    scheme: str,
    allowed_species: List[str],
    seq_aa_header: str,
    seq_dna_header: str = None,
    retain_indices: bool = False,
    display_residue_view: bool = False,
    display_unformatted_residue_view: bool = False,
    skip_run_base_anarci: bool = False,
    num_cpu: Optional[int] = None,
):
    # Perform data validation:
    validator = Validator(
        seq_aa_header=seq_aa_header, display_residue_view=display_residue_view
    )

    validation_result = validator.validate(df)

    # Perform ANARCI / AbNumber annotation on one dataframe
    df, seqs = create_row_id(df, seq_aa_header)

    assert (
        seq_aa_header in df.columns
    ), f"{seq_aa_header} is not a column in the DataFrame, please rename your amino acid input sequence column to {seq_aa_header}"

    if seq_dna_header is not None:
        df = validate_aa_translation(df, seq_aa_header, seq_dna_header)

    df_ab = ab.run_parallel_abnumber(
        df,
        seqs,
        scheme=scheme,
        allowed_species=allowed_species,
        seq_aa_header=seq_aa_header,
        seq_dna_header=seq_dna_header,
        retain_indices=retain_indices,
        display_residue_view=display_residue_view,
        display_unformatted_residue_view=display_unformatted_residue_view,
        num_cpu=num_cpu,
    )

    if not skip_run_base_anarci:
        df_an = an.run_parallel_anarci(
            df,
            seqs,
            scheme=scheme,
            allowed_species=allowed_species,
            seq_aa_header=seq_aa_header,
        )
        df_result_pd = df_ab.merge(df_an, on="seq_id")
    else:
        df_result_pd = df_ab

    df_result_pd.drop(["seq_id"], axis=1, inplace=True)

    logger.success(f"anarci-toolz run complete.")

    return df_result_pd


def main(
    input,
    scheme,
    allowed_species,
    seq_aa_header,
    seq_dna_header,
    retain_indices,
    display_residue_view,
    display_unformatted_residue_view,
    skip_run_base_anarci,
    num_cpu,
    tag,
):
    start_time = time.time()

    input_files = file_or_dir_input(input.split(","))

    out_dir = "anarci_annot"
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)  # Will create the folder in the current working directory

    for input_file in input_files:

        base_name = os.path.basename(input_file)  # get the base name
        file_name = os.path.splitext(base_name)[0]  # remove the file suffix

        logger.info(f"Processing {file_name} ...")

        df = pd.read_csv(input_file)

        df_result_pd = run_anarci_toolz(
            df=df,
            scheme=scheme,
            allowed_species=allowed_species,
            seq_aa_header=seq_aa_header,
            seq_dna_header=seq_dna_header,
            retain_indices=retain_indices,
            display_residue_view=display_residue_view,
            display_unformatted_residue_view=display_unformatted_residue_view,
            skip_run_base_anarci=skip_run_base_anarci,
            num_cpu=num_cpu,
        )

        output_filename = file_name + tag + ".csv"
        logger.info(f"Writing {output_filename} ...")
        df_result_pd.to_csv(os.path.join(out_dir, output_filename), index=False)

    end_time = time.time()
    elapsed_time = end_time - start_time
    logger.info(f"Total runtime: {round(elapsed_time,2)} seconds")


def validate_aa_translation(df, seq_aa_header, seq_dna_header):
    """Validates the amino acid sequence is direct translation of input DNA sequence. Generates a Boolean column "is_correct_translation" to be used for where to apply downstream functions, e.g. ANARCI-based DNA slicing
    :param df: DataFrame containing the amino acid sequence
    :type df: pd.DataFrame
    :param seq_aa_header: Header of the amino acid sequence column, defaults to 'sequence_aa'
    :type seq_aa_header: str, optional
    :return: DataFrame with the amino acid sequence and the validation status
    :rtype: pd.DataFrame
    """
    df[f"{seq_aa_header}_translation"] = df[seq_dna_header].apply(
        lambda x: str(Seq(x).translate())
    )
    df["is_correct_translation"] = (
        df[f"{seq_aa_header}_translation"] == df[seq_aa_header]
    )
    df = df.drop([f"{seq_aa_header}_translation"], axis=1)
    return df
