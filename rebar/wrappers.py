# Standard libraries
import os
from datetime import datetime
import logging
import yaml

# PyPI libraries
from Bio import SeqIO, Phylo


# rebar objects
from rebar.utils import (
    download_consensus_sequences,
    download_reference_sequence,
    parse_alignment,
    create_annotations,
    create_tree,
    create_barcodes,
    detect_recombination,
)
from rebar import RebarError

# Quiet URL fetch requests messages
logging.getLogger("requests").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)

# -----------------------------------------------------------------------------
# Dataset Subcommand
# -----------------------------------------------------------------------------


def dataset(params):
    """
    Download and create the rebar data model.

    Parameters
    ----------
        output : str
            file path for output barcodes csv.
        log : str
            file path for output log.
    """
    start_time = datetime.now()

    logger = params.logger
    logger.info(str(datetime.now()) + "\t" + "-" * 40)
    logger.info(str(datetime.now()) + "\tBEGINNING SUBCOMMAND: dataset.")

    info = {"name": params.name, "tag": params.tag}

    # Initialize tag
    if params.name == "sars-cov-2":

        if params.tag == "latest":

            create_annotations(params)
            accession = "MN908947.3"
            info["reference"] = download_reference_sequence(params, accession)
            info["sequences"] = download_consensus_sequences(params)
            info["tree"] = create_tree(params)
            info["barcodes"] = create_barcodes(params)

            # Summarize dataset in yaml file
            info_yaml = yaml.dump(info, sort_keys=False, indent=2)
            info_path = os.path.join(params.outdir, "dataset.yaml")
            logger.info(
                str(datetime.now()) + "\tExporting dataset summary: " + info_path
            )
            with open(info_path, "w") as outfile:
                outfile.write(info_yaml + "\n")

    # Runtime
    end_time = datetime.now()
    runtime = end_time - start_time
    # more than an hour
    if runtime.seconds > 3600:
        runtime = round(runtime.seconds / 3600, 1)
        units = "hours"
    # more than a minute
    elif runtime.seconds > 60:
        runtime = round(runtime.seconds / 60, 1)
        units = "minutes"
    else:
        runtime = round(runtime.seconds, 1)
        units = "seconds"

    # Finish
    logger.info(str(datetime.now()) + "\t" + "-" * 40)
    logger.info(str(datetime.now()) + "\tRUNTIME: " + str(runtime) + " " + units)
    logger.info(str(datetime.now()) + "\tFINISHED SUBCOMMAND: dataset.")


# -----------------------------------------------------------------------------
# Test Subcommand
# -----------------------------------------------------------------------------


def run(params):
    """
    Run rebar on an alignment or user-specified lineages.
    """
    start_time = datetime.now()

    logger = params.logger
    logger.info(str(datetime.now()) + "\t" + "-" * 40)
    logger.info(str(datetime.now()) + "\tBEGINNING SUBCOMMAND: run.")

    # Check for either lineages or alignment
    if not params.lineages and not params.alignment:
        raise SystemExit(
            RebarError(
                "RebarError: Either --lineages or --alignment must be specified."
            )
        )
    if params.lineages and params.alignment:
        raise SystemExit(
            RebarError(
                "RebarError: Both --lineages and --alignment cannot be specified."
            )
        )

    reference_path = os.path.join(params.dataset, "reference.fasta")
    consensus_path = os.path.join(params.dataset, "alignment.fasta")
    tree_path = os.path.join(params.dataset, "tree.nwk")
    barcodes_path = os.path.join(params.dataset, "barcodes.tsv")
    clade_path = os.path.join(params.dataset, "lineage_to_clade.tsv")
    required_files = [reference_path, consensus_path, tree_path, barcodes_path]

    # Check all dataset files exist
    for file_path in required_files:
        if not os.path.exists(file_path):
            raise SystemExit(
                RebarError("RebarError: Can't find dataset file " + file_path)
            )

    # Search consensus sequences for target lineages
    if params.lineages:
        tree = Phylo.read(tree_path, "newick")
        lineage_list = params.lineages.split(",")

        # Check for "*" character to imply descendants
        for lineage in lineage_list:

            if lineage.endswith("*"):
                lineage_tree = next(tree.find_clades(lineage[:-1]))
                lineage_descendants = [c.name for c in lineage_tree.find_clades()]
                lineage_list.remove(lineage)
                lineage_list += lineage_descendants

        params.lineages = ",".join(lineage_list)

        logger.info(
            str(datetime.now())
            + "\tSearching consensus sequences for: "
            + params.lineages
        )
        fasta_lines = []
        records = SeqIO.parse(consensus_path, "fasta")
        for record in records:
            if record.id in lineage_list:
                fasta_lines.append(">" + str(record.id))
                fasta_lines.append(str(record.seq))

        if len(fasta_lines) == 0:
            raise SystemExit(
                RebarError("RebarError: No sequences were found for input lineages.")
            )

        alignment_path = os.path.join(params.outdir, "alignment.fasta")
        logger.info(
            str(datetime.now()) + "\tExporting lineage alignment to: " + alignment_path
        )
        with open(alignment_path, "w") as outfile:
            outfile.write("\n".join(fasta_lines) + "\n")

        params.alignment = alignment_path

    # Run pipeline
    params.reference = reference_path
    params.barcodes = barcodes_path
    params.tree = tree_path
    params.lineage_to_clade = clade_path

    params.subs = os.path.join(params.outdir, "subs.tsv")
    parse_alignment(params)

    params.edge_cases = True
    detect_recombination(params)

    # Runtime
    end_time = datetime.now()
    runtime = end_time - start_time
    # more than an hour
    if runtime.seconds > 3600:
        runtime = round(runtime.seconds / 3600, 1)
        units = "hours"
    # more than a minute
    elif runtime.seconds > 60:
        runtime = round(runtime.seconds / 60, 1)
        units = "minutes"
    else:
        runtime = round(runtime.seconds, 1)
        units = "seconds"

    # Finish
    logger.info(str(datetime.now()) + "\t" + "-" * 40)
    logger.info(str(datetime.now()) + "\tRUNTIME: " + str(runtime) + " " + units)
    logger.info(str(datetime.now()) + "\tFINISHED SUBCOMMAND: run.")
