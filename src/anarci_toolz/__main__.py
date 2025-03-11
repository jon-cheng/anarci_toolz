import argparse
from .pipeline import main as pipeline_main


def main():
    parser = argparse.ArgumentParser(
        description="Generates ANARCI and AbNumber annotations on pandas Dataframes of antibody sequences"
    )

    parser.add_argument(
        "--input", type=str, required=True, help="The root directory of input data"
    )
    parser.add_argument(
        "--scheme",
        type=str,
        required=True,
        choices=["imgt", "chothia", "kabat", "aho"],
        help="Numbering scheme: One of 'imgt', 'chothia', 'kabat', 'aho'",
    )
    parser.add_argument(
        "--seq_aa_header",
        type=str,
        required=True,
        help="The name of the column containing the amino acid sequences",
    )
    parser.add_argument(
        "--seq_dna_header",
        type=str,
        required=False,
        default=None,
        help="The name of the column containing the DNA sequences",
    )
    parser.add_argument(
        "--allowed_species",
        type=str,
        required=True,
        nargs="+",
        choices=["human", "mouse", "rat", "rabbit", "rhesus", "pig", "alpaca", "all"],
        help="Allowed species for germline assignment. Use 'all' to allow all species or one or multiple of 'human', 'mouse', 'rat', 'rabbit', 'rhesus', 'pig', 'alpaca'",
    )
    parser.add_argument(
        "--retain_indices",
        action="store_true",
        help="Whether to retain linear index columns in sequence annotation; defaults to False",
    )
    parser.add_argument(
        "--display_residue_view",
        action="store_true",
        help="Whether to display wideform view showing all residue positions as columns; defaults to False",
    )
    parser.add_argument(
        "--display_unformatted_residue_view",
        action="store_true",
        help="Override default formatted residue view that includes scheme prefix, only displaying the scheme residue positions as columns; defaults to False",
    )

    parser.add_argument(
        "--skip_run_base_anarci",
        action="store_true",
        help="Whether to skip running base ANARCI for detailed alignment quality metrics like e-value; reduces run-time; defaults to False",
    )
    parser.add_argument(
        "--num_cpu",
        type=int,
        default=None,
        help="The number of CPUs to use for parallel processing; defaults to max available CPUs on the machine",
    )
    parser.add_argument(
        "--tag",
        type=str,
        required=False,
        default="_anarci_annot",
        help="suffix added to the end of the output filenames",
    )

    args = parser.parse_args()

    if args.allowed_species == ["all"]:
        args.allowed_species = None
    pipeline_main(
        input=args.input,
        scheme=args.scheme,
        allowed_species=args.allowed_species,
        seq_aa_header=args.seq_aa_header,
        seq_dna_header=args.seq_dna_header,
        retain_indices=args.retain_indices,
        display_residue_view=args.display_residue_view,
        display_unformatted_residue_view=args.display_unformatted_residue_view,
        skip_run_base_anarci=args.skip_run_base_anarci,
        num_cpu=args.num_cpu,
        tag=args.tag,
    )


if __name__ == "__main__":
    main()
