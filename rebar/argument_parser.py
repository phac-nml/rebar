#!/usr/bin/env python3

import argparse
from .wrappers import (
    run,
    run_test,
    validate,
    lineage_tree,
    download_barcodes,
    detect_recombination,
    analyze_subs,
)
from . import version


def add_alignment_param(parser, required=False):
    text = "Path to alignment fasta."
    parser.add_argument("--alignment", required=required, type=str, help=text)


def add_barcodes_param(parser, required=False):
    text = "Input barcodes csv, from `barcodes` subcommand."
    parser.add_argument("--barcodes", required=required, type=str, help=text)


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


def add_min_length_param(parser, required=False):
    text = (
        "Minimum number of consecutive barcode positions from each parent (Default: 3)."
    )
    parser.add_argument(
        "--min-length", required=required, default=3, type=int, help=text
    )


def add_min_subs_param(parser, required=False):
    text = "Minimum number of lineage-determining substitutions from each parent (Default: 1)."
    parser.add_argument("--min-subs", required=required, default=1, type=int, help=text)


def add_outdir_param(parser, required=False):
    text = "Output directory for results (Default: output)."
    parser.add_argument(
        "--outdir", required=required, default="output", type=str, help=text
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

    # Note: tree has no specific subcommands

    if subcommand == "test":
        add_lineages_param(parser=required, required=True)

    # subs subcommand
    if subcommand == "subs":
        add_alignment_param(parser=required, required=True)
        add_reference_param(parser=required, required=True)

        add_mask_param(parser=parser, required=False)

    if subcommand == "recombination":
        add_barcodes_param(parser=required, required=True)
        add_subs_param(parser=required, required=True)
        add_tree_param(parser=required, required=True)

        add_exclude_non_recomb_param(parser=parser, required=False)
        add_max_depth_param(parser=parser, required=False)
        add_min_length_param(parser=parser, required=False)
        add_min_subs_param(parser=parser, required=False)
        add_shared_param(parser=parser, required=False)
        add_snipit_format_param(parser=parser, required=False)

    # for run, alignment is mandatory
    if subcommand == "run":
        add_alignment_param(parser=required, required=True)
    # for validate, alignment is optional
    if subcommand == "validate":
        add_alignment_param(parser=required, required=False)

    # run, validate, and test need everything
    if subcommand in ["run", "validate", "test"]:

        add_reference_param(parser=required, required=True)

        add_barcodes_param(parser=required, required=False)
        add_subs_param(parser=required, required=False)
        add_tree_param(parser=required, required=False)

        add_exclude_non_recomb_param(parser=parser, required=False)
        add_mask_param(parser=parser, required=False)
        add_max_depth_param(parser=parser, required=False)
        add_min_length_param(parser=parser, required=False)
        add_min_subs_param(parser=parser, required=False)
        add_shared_param(parser=parser, required=False)
        add_snipit_format_param(parser=parser, required=False)


# -----------------------------------------------------------------------------
def make_parser():

    rebar_desc = "rebar: REcombination BARcode detection\n\n"
    parser = argparse.ArgumentParser(description="", usage=rebar_desc)

    rebar_subcommand_desc = (
        "rebar implements the following sub-commands:\n\n"
        "\trun\t\tRun all subcommands.\n"
        "\tvalidate\tValidate rebar on all designated lineages.\n"
        "\ttree\t\tCreate nomenclature tree of designated pango lineages.\n"
        "\tbarcodes\tDownload lineage barcodes from freyja-data.\n"
        "\tsubs\t\tSummarize substitutions observed in alignment.\n"
        "\trecombination\tDetect recombination from nextclade output.\n"
        "\thelp\t\tPrint rebar subcommands.\n"
        "\tversion\t\tPrint rebar version.\n"
        "\n"
    )
    parser.set_defaults(func=lambda x: print(rebar_desc + rebar_subcommand_desc))

    subparsers = parser.add_subparsers()

    # help
    help_parser = subparsers.add_parser("help", description="Print rebar subcommands.")
    help_parser.set_defaults(func=lambda x: print(rebar_desc + rebar_subcommand_desc))

    # version
    version_parser = subparsers.add_parser(
        "version", description="Print rebar version."
    )
    version_parser.set_defaults(func=lambda x: print("rebar v" + version))

    # barcodes
    barcodes_desc = (
        "Download barcodes csv of mutations observed in each pango-lineage."
        " Uses the usher_barcodes.csv from https://github.com/andersen-lab/Freyja-data\n\n"
    )
    barcodes_parser = subparsers.add_parser("barcodes", description=barcodes_desc)
    add_params(barcodes_parser, subcommand="barcodes")
    barcodes_parser.set_defaults(func=download_barcodes)

    # tree
    tree_desc = (
        "Creates a nomenclature tree of designated pango lineages."
        " Uses the lineages_notes.txt from"
        " https://github.com/cov-lineages/pango-designation"
        " and the pango_aliasor library to resolve descendants and aliases.\n\n"
    )
    tree_parser = subparsers.add_parser("tree", description=tree_desc)
    add_params(tree_parser, subcommand="tree")
    tree_parser.set_defaults(func=lineage_tree)

    # subs
    subs_desc = "Create nextclade-like TSV output from alignment."
    subs_parser = subparsers.add_parser("subs", description=subs_desc)
    add_params(subs_parser, subcommand="subs")
    subs_parser.set_defaults(func=analyze_subs)

    # recombination
    rec_desc = "Detect recombination."
    recombination_parser = subparsers.add_parser("recombination", description=rec_desc)
    add_params(recombination_parser, subcommand="recombination")
    recombination_parser.set_defaults(func=detect_recombination)

    # test
    test_desc = "Test rebar on select lineages.\n\n"
    test_parser = subparsers.add_parser("test", description=test_desc)
    add_params(test_parser, subcommand="test")
    test_parser.set_defaults(func=run_test)

    # validate
    validate_desc = "Validate rebar on all designated lineages.\n\n"
    validate_parser = subparsers.add_parser("validate", description=validate_desc)
    add_params(validate_parser, subcommand="validate")
    validate_parser.set_defaults(func=validate)

    # run
    run_desc = "Run the full rebar pipeline."
    run_parser = subparsers.add_parser("run", description=run_desc)
    add_params(run_parser, subcommand="run")
    run_parser.set_defaults(func=run)

    return parser
