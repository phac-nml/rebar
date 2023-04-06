import os
from datetime import datetime
import requests
import logging
from io import StringIO
import pandas as pd
from .utils import create_logger
from pango_aliasor.aliasor import Aliasor
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade

logging.getLogger("requests").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)

# URL of pango-designation lineages for lineage_tree
LINEAGES_URL = (
    "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/"
    "lineage_notes.txt"
)

BARCODES_URL = "https://github.com/andersen-lab/Freyja-data/raw/main/usher_barcodes.csv"


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
    logger = create_logger(params.log)
    logger.info(str(datetime.now()) + "\tBeginning module: barcodes.")

    # Create output directory if it doesn't exist
    outdir = os.path.dirname(params.output)
    if not os.path.exists(outdir) and outdir != "":
        logger.info(str(datetime.now()) + "\tCreating output directory: " + outdir)
        os.mkdir(outdir)

    # -------------------------------------------------------------------------
    # Download latest lineage barcodes

    logger.info(str(datetime.now()) + "\tDownloading lineage barcodes.")
    r = requests.get(BARCODES_URL)
    barcodes_text = r.text

    logger.info(str(datetime.now()) + "\tValidating barcode structure.")
    try:
        df = pd.read_csv(StringIO(barcodes_text), sep=",")
        logger.info(str(datetime.now()) + "\tValidation successful.")
    except Exception as e:
        logger.info(str(datetime.now()) + "\tValidation failed.")
        logger.info(str(datetime.now()) + "\t\t" + e.message + " " + e.args)
        return 1

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

    # Convert the text table to list
    lineages = []
    for line in lineage_text.split("\n"):
        if "Withdrawn" in line or line.startswith("Lineage"):
            continue

        lineage = line.split("\t")[0]
        if lineage == "":
            continue
        lineages.append(lineage)

    # Initialize the aliasor, which will download the latest aliases
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
