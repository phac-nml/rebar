#!/usr/bin/env python3

import argparse
from .wrappers import (
    lineage_tree,
    download_barcodes,
    detect_recombination,
    analyze_subs,
)
from . import version

rebar_description = "rebar: REcombination BARcode detection\n\n"
rebar_subcommand_description = (
    "rebar implements the following sub-commands:\n\n"
    "\ttree\t\tCreate nomenclature tree of designated pango lineages.\n"
    "\tbarcodes\tDownload lineage barcodes from freyja-data.\n"
    "\tsubs\t\tSummarize substitutions observed in alignment.\n"
    "\trecombination\tDetect recombination from nextclade output.\n"
    "\thelp\t\tPrint rebar subcommands.\n"
    "\tversion\t\tPrint rebar version.\n"
    "\n"
)
# rebar_usage_description = ""

log_description = "Output log file path.\n\n"
threads_description = "Number of threads to use."

# -----------------------------------------------------------------------------
tree_description = (
    "Creates a nomenclature tree of designated pango lineages."
    " Uses the lineages_notes.txt from"
    " https://github.com/cov-lineages/pango-designation"
    " and the pango_aliasor library to resolve descendants and aliases.\n\n"
)

tree_output_description = "Output file path for the newick tree.\n\n"

# -----------------------------------------------------------------------------
barcodes_description = (
    "Download barcodes csv of mutations observed in each pango-lineage."
    " Uses the usher_barcodes.csv from https://github.com/andersen-lab/Freyja-data\n\n"
)
barcodes_output_description = "Output file path for the barcodes csv.\n\n"

# -----------------------------------------------------------------------------
subs_description = "Create nextclade-like TSV output from alignment."
subs_alignment_description = "Input file path for the alignment fasta.\n\n"
subs_reference_description = "Input file path for the reference fasta.\n\n"
subs_output_description = "Output file path for the subs tsv.\n\n"

# -----------------------------------------------------------------------------
recombination_description = "Detect recombination.\n\n"
recombination_tree_description = "Input newick path for the lineage tree.\n\n"
recombination_barcodes_description = "Input csv path for the lineage barcodes.\n\n"
recombination_nextclade_description = "Input tsv path for the nextclade qc.\n\n"
recombination_outdir_description = "Output directory for detection results.\n\n"
recombination_exclude_shared_description = (
    "Exclude mutations shared by all parents when exporting.\n\n"
)


def add_barcodes_group(parser):
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "--output", required=True, type=str, help=barcodes_output_description
    )
    parser.add_argument("--log", required=False, type=str, help=log_description)


def add_recombination_group(parser):
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "--tree", required=True, type=str, help=recombination_tree_description
    )
    required.add_argument(
        "--barcodes", required=True, type=str, help=recombination_barcodes_description
    )
    required.add_argument(
        "--nextclade", required=True, type=str, help=recombination_nextclade_description
    )
    required.add_argument(
        "--outdir", required=True, type=str, help=recombination_outdir_description
    )
    parser.add_argument(
        "--exclude-shared",
        required=False,
        action="store_true",
        help=recombination_exclude_shared_description,
    )
    parser.add_argument("--log", required=False, type=str, help=log_description)


def add_subs_group(parser):
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "--output", required=True, type=str, help=subs_output_description
    )
    required.add_argument(
        "--alignment", required=True, type=str, help=subs_alignment_description
    )
    required.add_argument(
        "--reference", required=True, type=str, help=subs_reference_description
    )
    required.add_argument(
        "--threads", required=True, type=int, help=threads_description
    )
    parser.add_argument("--log", required=False, type=str, help=log_description)


def add_tree_group(parser):
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "--output", required=True, type=str, help=tree_output_description
    )
    parser.add_argument("--log", required=False, type=str, help=log_description)


def make_parser():
    parser = argparse.ArgumentParser(description="", usage=rebar_description)
    subparsers = parser.add_subparsers()

    # help
    help_parser = subparsers.add_parser("help", description="Print rebar subcommands.")
    help_parser.set_defaults(
        func=lambda x: print(rebar_description + rebar_subcommand_description)
    )

    # version
    version_parser = subparsers.add_parser(
        "version", description="Print rebar version."
    )
    version_parser.set_defaults(func=lambda x: print("rebar v" + version))

    # tree
    tree_parser = subparsers.add_parser("tree", description=tree_description)
    add_tree_group(tree_parser)
    tree_parser.set_defaults(func=lineage_tree)

    # barcodes
    barcodes_parser = subparsers.add_parser(
        "barcodes", description=barcodes_description
    )
    add_barcodes_group(barcodes_parser)
    barcodes_parser.set_defaults(func=download_barcodes)

    # subs
    subs_parser = subparsers.add_parser("subs", description=subs_description)
    add_subs_group(subs_parser)
    subs_parser.set_defaults(func=analyze_subs)

    # recombination
    recombination_parser = subparsers.add_parser(
        "recombination", description=recombination_description
    )
    add_recombination_group(recombination_parser)
    recombination_parser.set_defaults(func=detect_recombination)

    return parser
