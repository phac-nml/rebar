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
    text = "Dataset name."
    parser.add_argument(
        "--name", required=required, type=str, choices=["sars-cov-2"], help=text
    )


def add_debug_param(parser, required=False):
    text = "Log debugging output."
    parser.add_argument("--debug", required=required, action="store_true", help=text)


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
    text = "Maximum search depth to look for parents (Default: 20)."
    parser.add_argument(
        "--max-depth", required=required, default=20, type=int, help=text
    )


def add_min_consecutive_param(parser, required=False):
    text = (
        "Minimum number of consecutive barcode positions from each parent (Default: 3)."
    )
    parser.add_argument(
        "--min-consecutive", required=required, default=3, type=int, help=text
    )


def add_min_length_param(parser, required=False):
    text = "Minimum length of regions contributed by each parent (Default: 1000)."
    parser.add_argument(
        "--min-length", required=required, default=1000, type=int, help=text
    )


def add_min_subs_param(parser, required=False):
    text = "Minimum number of lineage-determining substitutions from each parent (Default: 1)."
    parser.add_argument("--min-subs", required=required, default=1, type=int, help=text)


def add_outdir_param(parser, required=False, default="output"):
    text = "Output directory for results (Default: {}).".format(default)
    parser.add_argument(
        "--outdir", required=required, default=default, type=str, help=text
    )


def add_plot_param(parser, required=False):
    text = "Create snipit plots."
    parser.add_argument("--plot", required=required, action="store_true", help=text)


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
        default="pdf",
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


def add_params(parser, subcommand=None):

    required = parser.add_argument_group("required arguments")

    # Universal to all subcommands, optional
    add_debug_param(parser, required=False)
    add_log_param(parser, required=False)
    add_outdir_param(parser, required=False)
    add_threads_param(parser, required=False)

    if subcommand == "dataset":
        add_dataset_name_param(parser=required, required=True)

    # Note: tree has no specific subcommands

    if subcommand == "run":
        add_dataset_param(parser=required, required=True)

        # One of these must be specified
        add_lineages_param(parser=parser, required=False)
        add_alignment_param(parser=parser, required=False)

        add_mask_param(parser=parser, required=False)
        add_exclude_non_recomb_param(parser=parser, required=False)
        add_max_depth_param(parser=parser, required=False)
        add_min_length_param(parser=parser, required=False)
        add_min_consecutive_param(parser=parser, required=False)
        add_min_subs_param(parser=parser, required=False)
        add_plot_param(parser=parser, required=False)
        add_shared_param(parser=parser, required=False)
        add_snipit_format_param(parser=parser, required=False)

    # # subs subcommand
    # if subcommand == "subs":
    #     add_alignment_param(parser=required, required=True)
    #     add_reference_param(parser=required, required=True)

    #     add_mask_param(parser=parser, required=False)

    # if subcommand == "recombination":
    #     add_barcodes_param(parser=required, required=True)
    #     add_subs_param(parser=required, required=True)
    #     add_tree_param(parser=required, required=True)

    #     add_exclude_non_recomb_param(parser=parser, required=False)
    #     add_max_depth_param(parser=parser, required=False)
    #     add_min_length_param(parser=parser, required=False)
    #     add_min_consecutive_param(parser=parser, required=False)
    #     add_min_subs_param(parser=parser, required=False)
    #     add_shared_param(parser=parser, required=False)
    #     add_snipit_format_param(parser=parser, required=False)

    # # for run, alignment is mandatory
    # if subcommand == "run":
    #     add_alignment_param(parser=required, required=True)
    # # for validate, alignment is optional
    # if subcommand == "validate":
    #     add_alignment_param(parser=required, required=False)

    # # run, validate, and test need everything
    # if subcommand in ["run", "validate", "test"]:

    #     add_reference_param(parser=required, required=True)

    #     add_barcodes_param(parser=required, required=False)
    #     add_subs_param(parser=required, required=False)
    #     add_tree_param(parser=required, required=False)

    #     add_exclude_non_recomb_param(parser=parser, required=False)
    #     add_mask_param(parser=parser, required=False)
    #     add_max_depth_param(parser=parser, required=False)
    #     add_min_length_param(parser=parser, required=False)
    #     add_min_consecutive_param(parser=parser, required=False)
    #     add_min_subs_param(parser=parser, required=False)
    #     add_shared_param(parser=parser, required=False)
    #     add_snipit_format_param(parser=parser, required=False)


# -----------------------------------------------------------------------------
def make_parser():

    # descriptions
    rebar_desc = "rebar: REcombination BARcode detection"
    dataset_desc = "Download and create the rebar data model."
    run_desc = "Run rebar on alignment or user-specified lineages."
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

    # # subs
    # subs_desc = "Create nextclade-like TSV output from alignment."
    # subs_parser = subparsers.add_parser("subs", description=subs_desc)
    # add_params(subs_parser, subcommand="subs")
    # subs_parser.set_defaults(func=analyze_subs)

    # # recombination
    # rec_desc = "Detect recombination."
    # recombination_parser = subparsers.add_parser("recombination", description=rec_desc)
    # add_params(recombination_parser, subcommand="recombination")
    # recombination_parser.set_defaults(func=detect_recombination)

    # # test
    # test_desc = "Test rebar on select lineages.\n\n"
    # test_parser = subparsers.add_parser("test", description=test_desc)
    # add_params(test_parser, subcommand="test")
    # test_parser.set_defaults(func=run_test)

    # # validate
    # validate_desc = "Validate rebar on all designated lineages.\n\n"
    # validate_parser = subparsers.add_parser("validate", description=validate_desc)
    # add_params(validate_parser, subcommand="validate")
    # validate_parser.set_defaults(func=validate)

    # # run
    # run_desc = "Run the full rebar pipeline."
    # run_parser = subparsers.add_parser("run", description=run_desc)
    # add_params(run_parser, subcommand="run")
    # run_parser.set_defaults(func=run)

    return parser
