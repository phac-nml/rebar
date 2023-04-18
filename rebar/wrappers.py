# Standard libraries
import os
from datetime import datetime
import functools
import requests
import logging
import re
from io import StringIO

# PyPI libraries
import pandas as pd
from tqdm import tqdm
from multiprocess import Pool  # Note that we are importing "multiprocess", no "ing"!
from pango_aliasor.aliasor import Aliasor
from Bio import Phylo, SeqIO
from Bio.Phylo.BaseTree import Clade

# rebar objects
from .utils import (
    download_consensus_sequences,
    NO_DATA_CHAR,
    LINEAGES_URL,
    BARCODES_USHER_URL,
    BARCODES_NEXTCLADE_URL,
    PROBLEMATIC_LINEAGES,
)
from .genome import genome_mp
from .export import Export
from . import RebarError
from .substitution import Substitution

# Quiet URL fetch requests messages
logging.getLogger("requests").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)


def run(params):
    """
    Run the full rebar pipeline.

    Parameters
    ----------
        output : str
            file path for output barcodes csv.
        log : str
            file path for output log.
    """
    #  run barcodes subcommand
    if not params.barcodes:
        params.output = os.path.join(params.outdir, "barcodes.csv")
        download_barcodes(params)
        params.barcodes = params.output

    # run tree subcommand
    if not params.tree:
        params.output = os.path.join(params.outdir, "tree.nwk")
        lineage_tree(params)
        params.tree = params.output

    # run subs subcommand
    if not params.subs:
        params.output = os.path.join(params.outdir, "subs.tsv")
        analyze_subs(params)
        params.subs = params.output

    # detect recombination
    detect_recombination(params)

    return 0


def run_test(params):
    """
    Test rebar on select lineages.
    """
    logger = params.logger
    logger.info(str(datetime.now()) + "\t" + "-" * 40)
    logger.info(str(datetime.now()) + "\tBEGINNING SUBCOMMAND: test.")

    # Download latest lineage consensus sequences
    logger.info(str(datetime.now()) + "\tDownloading lineage consensus sequences.")
    consensus_path = download_consensus_sequences(outdir=params.outdir)

    # Search consensus sequences for target lineages
    lineage_list = params.lineages.split(",")
    logger.info(
        str(datetime.now()) + "\tSearching consensus sequences for: " + params.lineages
    )
    fasta_lines = []
    records = SeqIO.parse(consensus_path, "fasta")
    for record in records:
        if record.id in lineage_list:
            fasta_lines.append(">" + str(record.id))
            fasta_lines.append(str(record.seq))

    alignment_path = os.path.join(params.outdir, "alignment.fasta")
    logger.info(str(datetime.now()) + "\tExporting alignment to: " + alignment_path)
    with open(alignment_path, "w") as outfile:
        outfile.write("\n".join(fasta_lines) + "\n")

    # Delete the consensus sequences file
    os.remove(consensus_path)

    # Run pipeline
    params.alignment = alignment_path
    run(params)

    # Finish
    logger.info(str(datetime.now()) + "\t" + "-" * 40)
    logger.info(str(datetime.now()) + "\tFINISHED SUBCOMMAND: test.")


def validate(params):
    """
    Validate rebar on all designated lineages.
    """
    logger = params.logger
    logger.info(str(datetime.now()) + "\t" + "-" * 40)
    logger.info(str(datetime.now()) + "\tBEGINNING SUBCOMMAND: validate.")

    # Download latest lineage consensus sequences
    logger.info(str(datetime.now()) + "\tDownloading lineage consensus sequences.")
    consensus_path = download_consensus_sequences(outdir=params.outdir)

    params.alignment = consensus_path

    run(params)

    logger.info(str(datetime.now()) + "\t" + "-" * 40)
    logger.info(str(datetime.now()) + "\tFINISHED SUBCOMMAND: validate.")


def download_barcodes(params):
    """
    Download csv of lineage barcodes from freyja data.

    Parameters
    ----------
        output : str
            file path for output barcodes csv.
        log : str
            file path for output log.
    """
    logger = params.logger
    logger.info(str(datetime.now()) + "\t" + "-" * 40)
    logger.info(str(datetime.now()) + "\tBEGINNING SUBCOMMAND: barcodes.")

    # -------------------------------------------------------------------------
    # Nextclade barcodes

    logger.info(str(datetime.now()) + "\tDownloading Nextclade barcodes.")
    r = requests.get(BARCODES_NEXTCLADE_URL)
    barcodes_data = r.json()
    barcodes_dict = {
        lineage: barcodes_data[lineage]["nucSubstitutions"] for lineage in barcodes_data
    }

    # -------------------------------------------------------------------------
    # UShER Barcodes

    logger.info(str(datetime.now()) + "\tDownloading UShER barcodes.")
    r = requests.get(BARCODES_USHER_URL)
    barcodes_text = r.text
    barcodes_usher_df = pd.read_csv(StringIO(barcodes_text), sep=",")
    # Rename the empty first column that should be lineage
    barcodes_usher_df.rename(columns={"Unnamed: 0": "lineage"}, inplace=True)
    # # Convert to dict
    # for rec in barcodes_usher_df.iterrows():
    #     lineage = rec[1]["lineage"]
    #     print(lineage)
    #     print(type(rec[1]))
    #     print(dir(rec[1]))
    #     x = rec[1].apply(lambda col: col.sum() == 1)
    #     print(x)
    #     break

    # Convert to dataframe
    logger.info(str(datetime.now()) + "\tConverting barcodes to dataframe.")

    lineages = list(barcodes_dict.keys())
    subs = [item for sublist in barcodes_dict.values() for item in sublist]
    subs = [s for s in set(subs) if s != ""]
    subs_detections = {s: [0] * len(barcodes_dict) for s in subs}
    for s in subs:
        for i, lineage in enumerate(lineages):
            if s in barcodes_dict[lineage]:
                subs_detections[s][i] = 1

    barcodes_nextclade_df = pd.DataFrame(subs_detections)
    barcodes_nextclade_df.insert(loc=0, column="lineage", value=lineages)

    # -------------------------------------------------------------------------
    # UShER Barcodes

    logger.info(str(datetime.now()) + "\tDownloading UShER barcodes.")
    r = requests.get(BARCODES_USHER_URL)
    barcodes_text = r.text
    barcodes_usher_df = pd.read_csv(StringIO(barcodes_text), sep=",")
    # Rename the empty first column that should be lineage
    barcodes_usher_df.rename(columns={"Unnamed: 0": "lineage"}, inplace=True)

    logger.info(str(datetime.now()) + "\tSupplementing missing Nextclade lineages.")

    nextclade_lineages = list(barcodes_nextclade_df["lineage"])
    usher_lineages = list(barcodes_usher_df["lineage"])

    usher_uniq = [l for l in usher_lineages if l not in nextclade_lineages]

    for lineage in usher_uniq:
        logger.info(str(datetime.now()) + "\t\tAdding UShER lineage " + lineage + ".")
        lineage_row = pd.DataFrame({s: [0] for s in barcodes_nextclade_df.columns})
        lineage_row["lineage"] = lineage

        barcodes_nextclade_df = pd.concat(
            [barcodes_nextclade_df, lineage_row], ignore_index=True
        )
        lineage_i = len(barcodes_nextclade_df) - 1

        df = barcodes_usher_df[barcodes_usher_df["lineage"] == lineage]

        detections = list(df.columns[df.apply(lambda col: col.sum() == 1)])
        for sub in detections:
            if sub not in barcodes_nextclade_df:
                barcodes_nextclade_df[sub] = 0
            barcodes_nextclade_df.at[lineage_i, sub] = 1

    # Sort columns by genomic position
    subs_order = sorted([Substitution(s) for s in barcodes_nextclade_df.columns[1:]])
    subs_order_str = [str(s) for s in subs_order]
    cols_order = ["lineage"] + subs_order_str
    barcodes_df = barcodes_nextclade_df[cols_order]

    logger.info(
        str(datetime.now())
        + "\tRemoving problematic lineages: "
        + ",".join(PROBLEMATIC_LINEAGES)
    )
    barcodes_df = barcodes_df[~barcodes_df["lineage"].isin(PROBLEMATIC_LINEAGES)]

    # Export
    barcodes_path = os.path.join(params.outdir, "barcodes.csv")
    logger.info(str(datetime.now()) + "\tExporting barcodes: " + barcodes_path)
    # Export to csv, to be consistent with usher barcodes format
    barcodes_df.to_csv(barcodes_path, sep=",", index=False)

    # Finish
    logger.info(str(datetime.now()) + "\tFINISHED SUBCOMMAND: barcodes.")
    return 0


def lineage_tree(params):
    """
    Create nomenclature tree of designated pango lineages.

    Parameters
    ----------
        output : str
            file path for output newick tree.
        log : str
            file path for output log.
    """

    logger = params.logger
    logger.info(str(datetime.now()) + "\t" + "-" * 40)
    logger.info(str(datetime.now()) + "\tBEGINNING SUBCOMMAND: tree.")

    # -------------------------------------------------------------------------
    # Download latest designated lineages from pango-designation

    logger.info(str(datetime.now()) + "\tDownloading designated lineage summaries.")
    r = requests.get(LINEAGES_URL)
    lineage_text = r.text

    # Attempt semi-structure text parsing for parents
    # This is most useful for identifying recursive recombinants
    recombinant_parents = {}

    # Convert the text table to list
    lineages = []

    # This could be parallelized, but it's already extremely fast
    for line in lineage_text.split("\n"):
        if "Withdrawn" in line or line.startswith("Lineage"):
            continue

        lineage = line.split("\t")[0]
        if lineage == "":
            continue
        lineages.append(lineage)

        # Check for recursive recombination
        # Recombinant lineage of .... Delta and BA.1
        if "Recombinant lineage" in line:
            line_split = line.split("Recombinant lineage")[1]
            line_words = line_split.split(" ")
            # Check for words that match lineage format
            # Bit risky, also catches "UK" if not comma after
            patterns = "^([A-Z]{2,}$|[A-Z]+\\.[0-9]+)"
            for word in line_words:
                lineage_matches = re.findall(patterns, word)
                if len(lineage_matches) > 0:
                    parent = word.replace(",", "").replace("*", "")
                    if lineage not in recombinant_parents:
                        recombinant_parents[lineage] = []
                    recombinant_parents[lineage].append(parent)

    # Initialize the aliasor, which will download the latest aliases
    logger.info(str(datetime.now()) + "\tInitialising aliases.")
    aliasor = Aliasor()

    # -------------------------------------------------------------------------
    # Construct Tree

    logger.info(str(datetime.now()) + "\tConstructing lineage tree.")

    # Create a tree with a root node "MRCA"
    tree = Clade(name="MRCA", clades=[], branch_length=1)
    # Add an "X" parent for recombinants
    clade = Clade(name="X", clades=[], branch_length=1)
    tree.clades.append(clade)

    # This can't be parallelized, sequential processing is required
    for lineage in lineages:

        # Identify the parent
        lineage_uncompress = aliasor.uncompress(lineage)
        parent_uncompress = ".".join(lineage_uncompress.split(".")[0:-1])
        parent = aliasor.compress(parent_uncompress)

        # Manual parents setting for A and B
        if lineage == "A":
            parent = "MRCA"

        elif lineage == "B":
            parent = "A"

        # Special handling for recombinants
        elif lineage.startswith("X") and parent == "":
            parent = "X"

        # Check for recursive recombinant
        if lineage in recombinant_parents:
            recursive_parents = [
                p for p in recombinant_parents[lineage] if p.startswith("X")
            ]
            if len(recursive_parents) > 0:
                parent = recursive_parents[0]

        parent_clade = [c for c in tree.find_clades(parent)]
        # If we found a parent, as long as the input list is formatted correctly
        # this should always be true
        if len(parent_clade) == 1:
            parent_clade = parent_clade[0]
            clade = Clade(name=lineage, clades=[], branch_length=1)
            parent_clade.clades.append(clade)

    # -------------------------------------------------------------------------
    # Export

    tree_path = os.path.join(params.outdir, "tree.nwk")
    logger.info(str(datetime.now()) + "\tExporting newick tree: " + tree_path)
    Phylo.write(tree, tree_path, "newick")

    # Finish
    logger.info(str(datetime.now()) + "\tFINISHED SUBCOMMAND: tree.")
    return 0


def analyze_subs(params):
    """
    Create nextclade-like TSV output from alignment.

    Parameters
    ----------
        output : str
            file path for output newick tree.
        log : str
            file path for output log.
    """

    logger = params.logger
    logger.info(str(datetime.now()) + "\t" + "-" * 40)
    logger.info(str(datetime.now()) + "\tBeginning subcommand: subs.")

    # Import reference
    logger.info(str(datetime.now()) + "\tImporting reference: " + params.reference)
    records = SeqIO.parse(params.reference, "fasta")
    ref_rec = next(records)

    # Import alignment
    logger.info(str(datetime.now()) + "\tImporting alignment: " + params.alignment)
    num_records = len(list(SeqIO.parse(params.alignment, "fasta")))
    records = SeqIO.parse(params.alignment, "fasta")

    # Parse substitutions
    logger.info(str(datetime.now()) + "\tParsing substitutions from alignment.")

    # Process genomes in parallel
    pool = Pool(params.threads)
    iterator = records
    task = functools.partial(
        genome_mp,
        reference=ref_rec,
        mask=params.mask,
        debug=params.debug,
        logger=params.logger,
    )
    total = num_records
    genomes = list(tqdm(pool.imap_unordered(task, iterator), total=total))

    # Pool memory management, don't accept anymore new tasks and wait
    pool.close()
    pool.join()

    # Export
    subs_path = os.path.join(params.outdir, "subs.tsv")
    logger.info(str(datetime.now()) + "\tExporting results to: " + subs_path)
    dfs = [genome.to_dataframe() for genome in genomes]
    df = pd.concat(dfs)
    df.to_csv(subs_path, sep="\t", index=False)

    # Finish
    logger.info(str(datetime.now()) + "\tFinished subcommand: subs.")
    return 0


def detect_recombination(params):
    """
    Detect recombination using lineage barcodes.

    Parameters
    ----------
        tree : str
            file path of input tree newick.
        barcodes : str
            file path of input barcodes csv.
        subs : str
            file path of input subs tsv.
        genome_length : int
            length of genome, only required if using nextclade TSV.
        logger : str
            output logging object (created in __main__.py).
    """

    logger = params.logger
    logger.info(str(datetime.now()) + "\t" + "-" * 40)
    logger.info(str(datetime.now()) + "\tBEGINNING SUBCOMMAND: recombination.")

    # -------------------------------------------------------------------------
    # PHASE 1: SETUP
    # Objectives:
    #   1. Import the substitutions dataframe of each sample.
    #   2. Import the lineage-determining mutation barcodes dataframe.
    #   3. Import the lineage nomenclature tree, to identify designated recombinants.

    # Import the dataframe from the `subs` module, or alternatively from nextclade
    logger.info(str(datetime.now()) + "\tImporting substitutions: " + params.subs)
    subs_df = pd.read_csv(params.subs, sep="\t").fillna(NO_DATA_CHAR)

    # We're going to need to know the genome length later, raise an error
    # if it wasn't provided.
    if "genome_length" not in subs_df.columns and not params.genome_length:
        raise RebarError(
            "If 'genome_length' is not a column in --subs TSV,"
            " --genome-length must be provided."
        )
    # If the user provided nextclade TSV as input instead of the rebar subs output
    # the ID column will be called "seqName" instead of strain
    if "seqName" in subs_df.columns:
        subs_df.rename({"seqName": "strain"}, inplace=True)

    subs_df.set_index("strain", inplace=True)
    subs_df["strain"] = subs_df.index

    # Import lineage barcodes
    logger.info(str(datetime.now()) + "\tImporting barcodes: " + params.barcodes)
    barcodes_df = pd.read_csv(params.barcodes, sep=",")

    # Import tree
    logger.info(str(datetime.now()) + "\tImporting tree: " + params.tree)
    tree = Phylo.read(params.tree, "newick")

    # Identify which lineages are known recombinants
    recombinant_tree = [c for c in tree.find_clades("X")][0]
    recombinant_lineages = [c.name for c in recombinant_tree.find_clades()]

    # -------------------------------------------------------------------------
    # PHASE 2: DETECT RECOMBINATION
    # Objectives:
    #   1. Store substitutions, deletions, and missing data as class attributes.
    #     - Ex. `genome.substitutions`, `genome.deletions`
    #   2. Summarize matches to lineage barcodes.
    #     - Ex. `genome.barcode_summary`
    #   3. Assign a lineage to each genome based on the highest barcode matches.
    #   4. Flag any samples that are definitively not recombinants.
    #     - Based on a perfect match to a non-recombinant lineage barcode.
    #     - These will be excluded from downstream processing to save time.
    #   5. Flag any samples that are recursive recombinants.
    #     - A perfect match to XBB.1.5 (BA.2.75 x BA.2.10) is not recursive.
    #     - A perfect match to XBL (XBB.1 x BA.2.75) is a recursive recombinant.
    #     - All other cases will have recursive set to `None`, for uncertainty
    #     - This will determine whether we allow recombinants to be the original
    #       parents in subsequent phases.

    logger.info(str(datetime.now()) + "\tDetecting recombination.")

    # Process genomes in parallel
    pool = Pool(params.threads)

    iterator = subs_df.iterrows()
    total = len(subs_df)

    # Debugging
    iterator = [rec for rec in subs_df.iterrows() if rec[1]["strain"].startswith("X")]
    # iterator = [rec for rec in subs_df.iterrows() if rec[1]["strain"] == "XD"]
    iterator = iterator[0:10]
    total = len(iterator)

    task = functools.partial(
        genome_mp,
        debug=params.debug,
        logger=params.logger,
        barcodes=barcodes_df,
        tree=tree,
        recombinant_tree=recombinant_tree,
        recombinant_lineages=recombinant_lineages,
        max_depth=params.max_depth,
        min_subs=params.min_subs,
        min_consecutive=params.min_consecutive,
        min_length=params.min_length,
    )

    genomes = list(tqdm(pool.imap(task, iterator), total=total))
    # Pool memory management, don't accept new tasks and wait
    pool.close()
    pool.join()

    # -------------------------------------------------------------------------
    # Pass 3: Export

    logger.info(str(datetime.now()) + "\tPreparing to export.")

    # If requested, exclude non-recombinants from output
    if params.exclude_non_recomb:
        genomes = [g for g in genomes if len(g.recombination.breakpoints) > 0]
    export = Export(genomes=genomes, outdir=params.outdir, include_shared=params.shared)

    # YAML
    logger.info(str(datetime.now()) + "\tExporting YAML.")
    export.to_yaml()

    logger.info(str(datetime.now()) + "\tExporting dataframe.")
    export.to_dataframe()

    logger.info(str(datetime.now()) + "\tExporting dataframes by parent.")
    export.to_dataframes()

    # Dataframe
    logger.info(str(datetime.now()) + "\tExporting alignments by parent.")
    export.to_alignments()

    # Snipit figures
    logger.info(str(datetime.now()) + "\tExporting snipit by parent.")
    export.to_snipit(ext=params.snipit_format)

    # Finish
    logger.info(str(datetime.now()) + "\tFINISHED SUBCOMMAND: recombination")
    return 0
