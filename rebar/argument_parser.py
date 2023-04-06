#!/usr/bin/env python3

import argparse
from .wrappers import lineage_tree, download_barcodes
from . import version

rebar_description = "rebar: REcombination BARcode detection\n\n"
rebar_subcommand_description = (
    "rebar implements the following sub-commands:\n\n"
    "\ttree\t\tCreate nomenclature tree of designated pango lineages.\n"
    "\tbarcodes\tDownload lineage barcodes from freyja-data.\n"
    "\thelp\t\tPrint rebar subcommands.\n"
    "\tversion\t\tPrint rebar version.\n"
    "\n"
)
# rebar_usage_description = ""

log_description = "Output log file path.\n\n"

tree_description = (
    "Creates a nomenclature tree of designated pango lineages."
    " Uses the lineages_notes.txt from"
    " https://github.com/cov-lineages/pango-designation"
    " and the pango_aliasor library to resolve descendants and aliases.\n\n"
)

tree_output_description = "Output file path for the newick tree.\n\n"

barcodes_description = (
    "Download barcodes csv of mutations observed in each pango-lineage."
    " Uses the usher_barcodes.csv from https://github.com/andersen-lab/Freyja-data\n\n"
)
barcodes_output_description = "Output file path for the barcodes csv.\n\n"


def add_barcodes_group(parser):
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "--output", required=True, type=str, help=barcodes_output_description
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

    return parser
