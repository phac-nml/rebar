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

rebar_desc = "rebar: REcombination BARcode detection\n\n"
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
# rebar_usage_desc = ""

log_desc = "Output log file path.\n\n"
threads_desc = "Number of threads to use."
alignment_desc = "Input file path for the alignment fasta.\n\n"
reference_desc = "Input file path for the reference fasta.\n\n"
genome_length_desc = (
    "Length of the reference genome, if genome_len column is not provided in subs.\n\n"
)
tree_input_desc = "Input newick path for the lineage tree.\n\n"
outdir_desc = "Output directory for results.\n\n"
mask_desc = "Number of bases to mask at 5' and 3' end of genome."
email_desc = "Email for authenticating NCBI reference genome download."
# -----------------------------------------------------------------------------
run_desc = "Run the full rebar pipeline."

# -----------------------------------------------------------------------------
tree_desc = (
    "Creates a nomenclature tree of designated pango lineages."
    " Uses the lineages_notes.txt from"
    " https://github.com/cov-lineages/pango-designation"
    " and the pango_aliasor library to resolve descendants and aliases.\n\n"
)

tree_output_desc = "Output file path for the newick tree.\n\n"

# -----------------------------------------------------------------------------
barcodes_desc = (
    "Download barcodes csv of mutations observed in each pango-lineage."
    " Uses the usher_barcodes.csv from https://github.com/andersen-lab/Freyja-data\n\n"
)
barcodes_output_desc = "Output file path for the barcodes csv.\n\n"

# -----------------------------------------------------------------------------
subs_desc = "Create nextclade-like TSV output from alignment."
subs_output_desc = "Output file path for the subs tsv.\n\n"

# -----------------------------------------------------------------------------
rec_desc = "Detect recombination.\n\n"
rec_barcodes_desc = "Input csv path for the lineage barcodes.\n\n"
rec_subs_desc = "Input tsv path from the subs subcommand output (or nextclade tsv).\n\n"
rec_max_depth_desc = (
    "The maximum search depth to look for parents through the top lineage matches.\n\n"
)
rec_min_subs_desc = (
    "The minimum number of consecutive barcode positions from each parent.\n\n"
)
rec_min_len_desc = "The minimum genomic region length each parent must contribute.\n\n"
rec_include_shared_desc = "Include mutations shared by all parents when exporting.\n\n"
rec_exclude_non_desc = "Include non-recombinant samples in output files.\n\n"
rec_snipit_format_desc = "snipit format extension for figures.\n\n"

# -----------------------------------------------------------------------------
validate_desc = "Validate rebar on all designated lineages.\n\n"

# -----------------------------------------------------------------------------
# run subcommand parameters


def add_run_group(parser):
    required = parser.add_argument_group("required arguments")

    required.add_argument("--alignment", required=True, type=str, help=alignment_desc)
    required.add_argument("--reference", required=True, type=str, help=reference_desc)

    parser.add_argument("--threads", required=False, type=int, help=threads_desc)
    parser.add_argument("--log", required=False, type=str, help=log_desc)
    parser.add_argument(
        "--outdir", required=False, type=str, default=".", help=outdir_desc
    )


# -----------------------------------------------------------------------------
# barcodes subcommand parameters
def add_barcodes_group(parser):
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "--output",
        required=False,
        type=str,
        default="barcodes.csv",
        help=barcodes_output_desc,
    )
    parser.add_argument("--log", required=False, type=str, help=log_desc)


# -----------------------------------------------------------------------------
# tree subcommand parameters
def add_tree_group(parser):
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "--output", required=False, type=str, default="tree.nwk", help=tree_output_desc
    )
    parser.add_argument("--log", required=False, type=str, help=log_desc)


# -----------------------------------------------------------------------------
# subs subcommand parameters
def add_subs_group(parser):

    required = parser.add_argument_group("required arguments")

    required.add_argument("--alignment", required=True, type=str, help=alignment_desc)
    required.add_argument("--reference", required=True, type=str, help=reference_desc)

    parser.add_argument("--mask", required=False, default=200, type=int, help=mask_desc)
    parser.add_argument(
        "--output", required=True, type=str, default="subs.tsv", help=subs_output_desc
    )
    parser.add_argument(
        "--threads", required=False, type=int, default=1, help=threads_desc
    )
    parser.add_argument("--log", required=False, type=str, help=log_desc)


# -----------------------------------------------------------------------------
# recombination subcommand parameters
def add_recombination_group(parser):

    required = parser.add_argument_group("required arguments")

    required.add_argument("--tree", required=True, type=str, help=tree_input_desc)
    required.add_argument("--barcodes", required=True, type=str, help=rec_barcodes_desc)
    required.add_argument("--subs", required=True, type=str, help=rec_subs_desc)

    parser.add_argument(
        "--outdir", required=False, type=str, default=".", help=outdir_desc
    )
    parser.add_argument(
        "--include-shared",
        required=False,
        action="store_true",
        help=rec_include_shared_desc,
    )
    parser.add_argument(
        "--threads", required=False, type=int, default=1, help=threads_desc
    )
    parser.add_argument("--log", required=False, type=str, help=log_desc)
    parser.add_argument(
        "--genome-length", required=False, type=int, help=genome_length_desc
    )
    parser.add_argument(
        "--max-depth", required=False, default=5, type=int, help=rec_max_depth_desc
    )
    parser.add_argument(
        "--min-subs", required=False, default=3, type=int, help=rec_min_subs_desc
    )
    parser.add_argument(
        "--min-length", required=False, default=1, type=int, help=rec_min_len_desc
    )
    parser.add_argument(
        "--exclude-non-recomb",
        required=False,
        action="store_true",
        help=rec_exclude_non_desc,
    )
    parser.add_argument(
        "--snipit-format",
        required=False,
        default="pdf",
        type=str,
        help=rec_snipit_format_desc,
    )


# -----------------------------------------------------------------------------
# validate subcommand parameters
def add_validate_group(parser):

    # validate specific params
    parser.add_argument("--email", required=False, type=str, help=email_desc)
    parser.add_argument("--reference", required=False, type=str, help=reference_desc)
    parser.add_argument("--alignment", required=False, type=str, help=alignment_desc)

    parser.add_argument("--subs", required=False, type=str, help=rec_subs_desc)
    parser.add_argument("--barcodes", required=False, type=str, help=rec_barcodes_desc)
    parser.add_argument("--tree", required=False, type=str, help=tree_input_desc)
    parser.add_argument("--mask", required=False, default=200, type=int, help=mask_desc)

    # recombination params
    parser.add_argument(
        "--include-shared",
        required=False,
        action="store_true",
        help=rec_include_shared_desc,
    )
    parser.add_argument(
        "--genome-length", required=False, type=int, help=genome_length_desc
    )
    parser.add_argument(
        "--max-depth", required=False, default=5, type=int, help=rec_max_depth_desc
    )
    parser.add_argument(
        "--min-subs", required=False, default=3, type=int, help=rec_min_subs_desc
    )
    parser.add_argument(
        "--min-length", required=False, default=1, type=int, help=rec_min_len_desc
    )
    parser.add_argument(
        "--exclude-non-recomb",
        required=False,
        action="store_true",
        help=rec_exclude_non_desc,
    )
    parser.add_argument(
        "--snipit-format",
        required=False,
        default="pdf",
        type=str,
        help=rec_snipit_format_desc,
    )

    # General
    parser.add_argument(
        "--threads", required=False, type=int, default=1, help=threads_desc
    )
    parser.add_argument("--log", required=False, type=str, help=log_desc)
    parser.add_argument(
        "--outdir", required=False, type=str, default=".", help=outdir_desc
    )


# -----------------------------------------------------------------------------
def make_parser():
    parser = argparse.ArgumentParser(description="", usage=rebar_desc)
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

    # run
    run_parser = subparsers.add_parser("run", description=run_desc)
    add_run_group(run_parser)
    run_parser.set_defaults(func=run)

    # tree
    tree_parser = subparsers.add_parser("tree", description=tree_desc)
    add_tree_group(tree_parser)
    tree_parser.set_defaults(func=lineage_tree)

    # barcodes
    barcodes_parser = subparsers.add_parser("barcodes", description=barcodes_desc)
    add_barcodes_group(barcodes_parser)
    barcodes_parser.set_defaults(func=download_barcodes)

    # subs
    subs_parser = subparsers.add_parser("subs", description=subs_desc)
    add_subs_group(subs_parser)
    subs_parser.set_defaults(func=analyze_subs)

    # recombination
    recombination_parser = subparsers.add_parser("recombination", description=rec_desc)
    add_recombination_group(recombination_parser)
    recombination_parser.set_defaults(func=detect_recombination)

    # validate
    validate_parser = subparsers.add_parser("validate", description=validate_desc)
    add_validate_group(validate_parser)
    validate_parser.set_defaults(func=validate)

    return parser
