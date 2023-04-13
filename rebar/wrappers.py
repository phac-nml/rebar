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
from .utils import create_logger, NO_DATA_CHAR, LINEAGES_URL, BARCODES_URL
from .genome import genome_mp
from .export import Export
from . import RebarError

# Quiet URL fetch requests messages
logging.getLogger("requests").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)


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
    logger.info(str(datetime.now()) + "\tBeginning module: barcodes.")

    # Download latest lineage barcodes
    logger.info(str(datetime.now()) + "\tDownloading lineage barcodes.")
    r = requests.get(BARCODES_URL)
    barcodes_text = r.text

    # Validate barcodes structure
    logger.info(str(datetime.now()) + "\tValidating barcode structure.")
    try:
        barcodes_df = pd.read_csv(StringIO(barcodes_text), sep=",")
        logger.info(str(datetime.now()) + "\tValidation successful.")
    except Exception as e:
        logger.info(str(datetime.now()) + "\tValidation failed.")
        logger.info(str(datetime.now()) + "\t\t" + e.message + " " + e.args)
        return 1

    # Export
    logger.info(str(datetime.now()) + "\tExporting barcodes: " + params.output)
    barcodes_df.to_csv(params.output, sep=",", index=False)

    # Finish
    logger.info(str(datetime.now()) + "\tFinished module: barcodes.")
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

    logger = create_logger(params.log)
    logger.info(str(datetime.now()) + "\tBeginning module: tree.")

    # Create output directory if it doesn't exist
    outdir = os.path.dirname(params.output)
    if not os.path.exists(outdir) and outdir != "":
        logger.info(str(datetime.now()) + "\tCreating output directory: " + outdir)
        os.mkdir(outdir)

    # -------------------------------------------------------------------------
    # Download latest designated lineages from pango-designation

    logger.info(str(datetime.now()) + "\tDownloading list of designated lineages.")
    r = requests.get(LINEAGES_URL)
    lineage_text = r.text

    # Attempt semi-structure text parsing for parents
    # This is most useful for identifying recursive recombinants
    recombinant_parents = {}

    # Convert the text table to list
    lineages = []
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

    logger.info(str(datetime.now()) + "\tExporting newick tree: " + params.output)
    Phylo.write(tree, params.output, "newick")

    logger.info(str(datetime.now()) + "\tFinished module: tree.")

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
    logger.info(str(datetime.now()) + "\tBeginning module: subs.")

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

    p = Pool(params.threads)
    genomes = list(
        tqdm(
            p.imap(functools.partial(genome_mp, reference=ref_rec), records),
            total=num_records,
        )
    )
    # Pool memory management, don't accept anymore new tasks and wait for them
    p.close()
    p.join()

    # Export
    logger.info(str(datetime.now()) + "\tExporting results to: " + params.output)
    dfs = [genome.to_dataframe() for genome in genomes]
    df = pd.concat(dfs)
    df.to_csv(params.output, sep="\t", index=False)

    # Finish
    logger.info(str(datetime.now()) + "\tFinished module: subs.")
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
    logger.info(str(datetime.now()) + "\tBeginning module: recombination.")

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
    # Rename the empty first column that should be lineage
    barcodes_df.rename(columns={"Unnamed: 0": "lineage"}, inplace=True)

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

    # The objects we're going to process in parallel
    # We specifically use iterrows, otherwise it iterates over the df cols
    iterator = subs_df.iterrows()
    # The function we're going to apply to the iterator in parallel
    # Note the use of functools.partial to later pass named arguments to `imap`
    # `genome_mp` is a multiprocessing wrapper for the Genome class
    task = functools.partial(
        genome_mp,
        barcodes=barcodes_df,
        tree=tree,
        recombinant_tree=recombinant_tree,
        recombinant_lineages=recombinant_lineages,
        max_depth=params.max_depth,
        min_subs=params.min_subs,
        min_length=params.min_length,
    )
    # The total number of objects, for the progress bar
    total = len(subs_df)
    # Create pool for multiprocessing tasks
    pool = Pool(params.threads)
    # Reminder: list() is required, to run the tasks and update progress bar
    genomes = list(tqdm(pool.imap_unordered(task, iterator), total=total))

    # Pool memory management, don't accept new tasks and wait for them to finish
    pool.close()
    pool.join()

    # -------------------------------------------------------------------------
    # Pass 3: Export

    logger.info(str(datetime.now()) + "\tPreparing to export.")

    # Unless requested, exclude non-recombinants from output
    if not params.export_non_rec:
        genomes = [g for g in genomes if len(g.recombination.breakpoints) > 0]
    export = Export(
        genomes=genomes, outdir=params.outdir, include_shared=params.include_shared
    )

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
    logger.info(str(datetime.now()) + "\tFinished module: recombination")
    return 0
