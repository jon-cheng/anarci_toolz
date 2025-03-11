from Bio.Seq import Seq


def validate_aa_translation(df, seq_aa_header, dna_header):
    """Validates the amino acid sequence is direct translation of input DNA sequence. Generates a Boolean column "is_correct_translation" to be used for where to apply downstream functions, e.g. ANARCI-based DNA slicing

    :param df: DataFrame containing the amino acid sequence
    :type df: pd.DataFrame
    :param seq_aa_header: Header of the amino acid sequence column, defaults to 'sequence_aa'
    :type seq_aa_header: str, optional
    :return: DataFrame with the amino acid sequence and the validation status
    :rtype: pd.DataFrame
    """
    df[f"{seq_aa_header}_translation"] = df[dna_header].apply(
        lambda x: str(Seq(x).translate())
    )
    df["is_correct_translation"] = (
        df[f"{seq_aa_header}_translation"] == df[seq_aa_header]
    )
    df = df.drop([f"{seq_aa_header}_translation"], axis=1)
    return df
