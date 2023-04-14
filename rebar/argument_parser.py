#!/usr/bin/env python3

import argparse
from .wrappers import (
    run,
    validate,
    lineage_tree,
    download_barcodes,
    detect_recombination,
    analyze_subs,
)
from . import version


def add_params(parser, subcommand=None):

    required = parser.add_argument_group("required arguments")

    # Universal to all subcommands, optional
    parser.add_argument(
        "--threads",
        required=False,
        type=int,
        default=1,
        help="Number of threads to use.",
    )
    parser.add_argument("--log", required=False, type=str, help="Log file path.")
    parser.add_argument(
        "--outdir",
        required=False,
        type=str,
        default=".",
        help="Output directory for results.",
    )

    if subcommand in ["subs", "run", "validate"]:
        required.add_argument(
            "--reference", required=True, type=str, help="Path to reference fasta."
        )
        parser.add_argument(
            "--mask",
            required=False,
            default=200,
            type=int,
            help="Number of bases to mask at 5' and 3' end of genome.",
        )

    if subcommand in ["recombination"]:
        required.add_argument(
            "--tree",
            required=True,
            type=str,
            help="Input newick tree, from `tree` subcommand.",
        )
        required.add_argument(
            "--barcodes",
            required=True,
            type=str,
            help="Input barcodes csv, from `barcodes` subcommand.",
        )
        required.add_argument(
            "--subs",
            required=True,
            type=str,
            help="Input subs tsv, from `subs` subcommand.",
        )

    if subcommand in ["recombination", "run", "validate"]:
        parser.add_argument(
            "--include-shared",
            required=False,
            action="store_true",
            help="Include mutations shared by all parents when exporting.",
        )
        parser.add_argument(
            "--max-depth",
            required=False,
            default=5,
            type=int,
            help="The maximum search depth to look for parents through the top lineage matches.",
        )
        parser.add_argument(
            "--min-subs",
            required=False,
            default=2,
            type=int,
            help="The minimum number of consecutive barcode positions from each parent.",
        )
        parser.add_argument(
            "--min-length",
            required=False,
            default=1,
            type=int,
            help="The minimum genomic region length each parent must contribute.",
        )
        parser.add_argument(
            "--exclude-non-recomb",
            required=False,
            action="store_true",
            help="Exclude non-recombinant samples from output files",
        )
        parser.add_argument(
            "--snipit-format",
            required=False,
            default="pdf",
            choices=["pdf", "png"],
            type=str,
            help="snipit format extension for figures",
        )

    if subcommand in ["subs", "recombination", "run"]:
        required.add_argument(
            "--alignment", required=True, type=str, help="Path to alignment fasta."
        )

    if subcommand in ["validate"]:
        parser.add_argument(
            "--alignment", required=False, type=str, help="Path to alignment fasta."
        )
        required.add_argument(
            "--tree",
            required=False,
            type=str,
            help="Input newick tree, from `tree` subcommand.",
        )
        required.add_argument(
            "--barcodes",
            required=False,
            type=str,
            help="Input barcodes csv, from `barcodes` subcommand.",
        )
        required.add_argument(
            "--subs",
            required=False,
            type=str,
            help="Input subs tsv, from `subs` subcommand.",
        )


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
