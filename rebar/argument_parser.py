#!/usr/bin/env python3

import argparse
from .wrappers import dataset, run
from . import version


def add_alignment_param(parser, required=False):
    text = "Path to alignment fasta."
    parser.add_argument("--alignment", required=required, type=str, help=text)


def add_barcodes_param(parser, required=False):
    text = "Input barcodes csv, from `barcodes` subcommand."
    parser.add_argument("--barcodes", required=required, type=str, help=text)


def add_dataset_param(parser, required=False):
    text = "Path to dataset directory, output of `dataset` subcommand ."
    parser.add_argument("--dataset", required=required, type=str, help=text)


def add_dataset_name_param(parser, required=False):
    text = "Dataset name (Default: sars-cov-2)."
    parser.add_argument(
        "--name", required=required, type=str, choices=["sars-cov-2"], help=text
    )


def add_dataset_tag_param(parser, required=False):
    text = "Dataset tag (Default: latest)."
    parser.add_argument(
        "--tag", required=required, type=str, choices=["latest"], help=text
    )


def add_debug_param(parser, required=False):
    text = "Enable debugging mode."
    parser.add_argument("--debug", required=required, action="store_true", help=text)


def add_no_edge_cases_param(parser, required=False):
    text = "Disable sensitive edge case handling for recombinants with few mutations."
    parser.add_argument(
        "--no-edge-cases", required=required, action="store_true", help=text
    )


def add_exclude_non_recomb_param(parser, required=False):
    text = "Exclude non-recombinant samples from output files"
    parser.add_argument(
        "--exclude-non-recomb", required=required, action="store_true", help=text
    )


def add_log_param(parser, required=False):
    text = "Log file path."
    parser.add_argument("--log", required=required, type=str, help=text)


def add_lineages_param(parser, required=False):
    text = "Comma-separated list of lineages to test (ex. 'XBF,XBB.1.5')."
    parser.add_argument("--lineages", required=required, type=str, help=text)


def add_mask_param(parser, required=False):
    text = "Number of bases to mask at 5' and 3' end of genome (Default: 200)."
    parser.add_argument("--mask", required=required, default=200, type=int, help=text)


def add_max_depth_param(parser, required=False):
    text = "Maximum search depth to look for parents (Default: 3)."
    parser.add_argument(
        "--max-depth", required=required, default=3, type=int, help=text
    )


def add_min_consecutive_param(parser, required=False):
    text = (
        "Minimum number of consecutive barcode positions from each parent (Default: 3)."
    )
    parser.add_argument(
        "--min-consecutive", required=required, default=3, type=int, help=text
    )


def add_min_length_param(parser, required=False):
    text = "Minimum length of regions contributed by each parent (Default: 500)."
    parser.add_argument(
        "--min-length", required=required, default=500, type=int, help=text
    )


def add_min_subs_param(parser, required=False):
    text = "Minimum number of lineage-determining substitutions from each parent (Default: 1)."
    parser.add_argument("--min-subs", required=required, default=1, type=int, help=text)


def add_outdir_param(parser, required=False, default="output"):
    text = "Output directory for results (Default: {}).".format(default)
    parser.add_argument(
        "--outdir", required=required, default=default, type=str, help=text
    )


def add_output_all_param(parser, required=False):
    text = "Produce all possible output files."
    parser.add_argument(
        "--output-all", required=required, action="store_true", help=text
    )


def add_output_fasta_param(parser, required=False):
    text = "Output FASTA results (alignments with non-barcode positions masked)."
    parser.add_argument(
        "--output-fasta", required=required, action="store_true", help=text
    )


def add_output_plot_param(parser, required=False):
    text = "Output snipit plots."
    parser.add_argument(
        "--output-plot", required=required, action="store_true", help=text
    )


def add_output_barcode_param(parser, required=False):
    text = "Output barcode mutations as table."
    parser.add_argument(
        "--output-barcode", required=required, action="store_true", help=text
    )


def add_output_tsv_param(parser, required=False):
    text = "Output TSV summary."
    parser.add_argument(
        "--output-tsv", required=required, action="store_true", help=text
    )


def add_output_yaml_param(parser, required=False):
    text = "Output YAML summary."
    parser.add_argument(
        "--output-yaml", required=required, action="store_true", help=text
    )


def add_reference_param(parser, required=False):
    text = "Path to reference fasta."
    parser.add_argument("--reference", required=required, type=str, help=text)


def add_shared_param(parser, required=False):
    text = "Include mutations shared by all parents when exporting."
    parser.add_argument("--shared", required=required, action="store_true", help=text)


def add_snipit_format_param(parser, required=False):
    text = "snipit format extension for figures"
    parser.add_argument(
        "--snipit-format",
        required=required,
        default="png",
        choices=["pdf", "png"],
        type=str,
        help=text,
    )


def add_subs_param(parser, required=False):
    text = "Input subs tsv, from `subs` subcommand."
    parser.add_argument("--subs", required=required, type=str, help=text)


def add_threads_param(parser, required=False):
    text = "Number of threads to use (Default: 1)."
    parser.add_argument("--threads", required=required, type=int, default=1, help=text)


def add_tree_param(parser, required=False):
    text = "Input newick tree, from `tree` subcommand."
    parser.add_argument("--tree", required=required, type=str, help=text)


def add_validate_param(parser, required=False):
    text = "Validate lineages against expected values."
    parser.add_argument("--validate", required=required, action="store_true", help=text)


def add_params(parser, subcommand=None):

    parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")

    if subcommand == "dataset":

        # Mandatory
        add_dataset_name_param(parser=required, required=True)
        add_dataset_tag_param(parser=required, required=True)
        add_outdir_param(parser=optional, required=False)

        # Optional General
        add_debug_param(parser=optional, required=False)
        add_log_param(parser=optional, required=False)
        add_threads_param(parser=optional, required=False)

    elif subcommand == "run":

        # Mandatory
        add_dataset_param(parser=required, required=True)

        # Mutually exclusive arguments
        add_lineages_param(parser=optional, required=False)
        add_alignment_param(parser=optional, required=False)

        # Optional General
        add_debug_param(parser=optional, required=False)
        add_log_param(parser=optional, required=False)
        add_outdir_param(parser=optional, required=False)
        add_threads_param(parser=optional, required=False)

        # Optional Specific
        add_no_edge_cases_param(parser=optional, required=False)
        add_mask_param(parser=optional, required=False)
        add_exclude_non_recomb_param(parser=optional, required=False)
        add_max_depth_param(parser=optional, required=False)
        add_min_length_param(parser=optional, required=False)
        add_min_consecutive_param(parser=optional, required=False)
        add_min_subs_param(parser=optional, required=False)
        add_shared_param(parser=optional, required=False)
        add_snipit_format_param(parser=optional, required=False)
        add_validate_param(parser=optional, required=False)

        # Optional Output
        add_output_all_param(parser=optional, required=False)
        add_output_fasta_param(parser=optional, required=False)
        add_output_plot_param(parser=optional, required=False)
        add_output_barcode_param(parser=optional, required=False)
        add_output_tsv_param(parser=optional, required=False)
        add_output_yaml_param(parser=optional, required=False)


# -----------------------------------------------------------------------------
def make_parser():

    # descriptions
    rebar_desc = "rebar: REcombination BARcode detection"
    dataset_desc = "Download and create the rebar data model."
    run_desc = "Run rebar on an alignment or user-specified lineages."
    version_desc = "Print version."
    help_desc = "Print subcommands."
    parser = argparse.ArgumentParser(description="", usage=rebar_desc)

    rebar_subcommand_desc = (
        "\n\n"
        + "\tdataset\t\t"
        + dataset_desc
        + "\n"
        + "\trun\t\t"
        + run_desc
        + "\n"
        + "\thelp\t\t"
        + help_desc
        + "\n"
        + "\tversion\t\t"
        + version_desc
        + "\n"
    )
    parser.set_defaults(func=lambda x: print(rebar_desc + rebar_subcommand_desc))

    subparsers = parser.add_subparsers()

    # help
    help_parser = subparsers.add_parser("help", description=help_desc)
    help_parser.set_defaults(func=lambda x: print(rebar_desc + rebar_subcommand_desc))

    # version
    version_parser = subparsers.add_parser("version", description=version_desc)
    version_parser.set_defaults(func=lambda x: print("rebar v" + version))

    # dataset
    dataset_parser = subparsers.add_parser("dataset", description=dataset_desc)
    add_params(dataset_parser, subcommand="dataset")
    dataset_parser.set_defaults(func=dataset)

    # run
    run_parser = subparsers.add_parser("run", description=run_desc)
    add_params(run_parser, subcommand="run")
    run_parser.set_defaults(func=run)

    return parser
