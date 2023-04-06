import os
from datetime import datetime
import requests
import logging
from io import StringIO
import pandas as pd
from .utils import create_logger, NO_DATA_CHAR, GENOME_LEN
from .genome import Genome
from pango_aliasor.aliasor import Aliasor
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade

logging.getLogger("requests").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)

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


def detect_recombination(params):
    """
    Detect recombination using lineage barcodes.

    Parameters
    ----------
        tree : str
            file path of input lineage nomenclature tree.
        barcodes : str
            file path of input barcodes csv.
        nextclade : str
            file path of input nextclade tsv.
        log : str
            file path for output log.
    """
    logger = create_logger(params.log)
    logger.info(str(datetime.now()) + "\tBeginning module: recombination")

    # Create output directory if it doesn't exist
    outdir = os.path.dirname(params.outdir)
    if not os.path.exists(outdir) and outdir != "":
        logger.info(str(datetime.now()) + "\tCreating output directory: " + outdir)
        os.mkdir(outdir)

    # -------------------------------------------------------------------------
    # Inputs

    # barcodes
    logger.info(str(datetime.now()) + "\tImporting barcodes: " + params.barcodes)
    barcodes_df = pd.read_csv(params.barcodes, sep=",")
    # Rename the empty column that should be lineage
    barcodes_df.rename(columns={"Unnamed: 0": "lineage"}, inplace=True)

    # tree
    logger.info(str(datetime.now()) + "\tImporting tree: " + params.tree)
    tree = Phylo.read(params.tree, "newick")
    # Identify which lineages are known recombinants
    recombinant_tree = [c for c in tree.find_clades("X")][0]
    recombinant_lineages = [c.name for c in recombinant_tree.find_clades()]

    # nextclade
    logger.info(str(datetime.now()) + "\tImporting nextclade: " + params.nextclade)
    nextclade_df = pd.read_csv(params.nextclade, sep="\t").fillna(NO_DATA_CHAR)

    # parse subs into genomic record
    logger.info(str(datetime.now()) + "\tCreating genome records.")
    strains = list(nextclade_df["seqName"])
    genomes = []
    for strain in strains:
        nexclade_row = nextclade_df.query("seqName == @strain")
        genome = Genome(
            id=strain,
            genome_len=GENOME_LEN,
            nextclade_row=nexclade_row,
            barcode=barcodes_df,
        )
        genomes.append(genome)

    # -------------------------------------------------------------------------
    # PASS 0: Lineage assignment
    # Identify the lineage (or set of lineages) whose barcodes most closely
    # match each genomes observed substitutions.

    logger.info(str(datetime.now()) + "\tAssigning lineage.")
    for genome in genomes:

        backbone = genome.identify_backbone(
            barcode_summary=genome.barcode_summary,
            barcodes_df=barcodes_df,
            tree=tree,
            recombinant_lineages=recombinant_lineages,
            include_recombinants=True,
        )
        genome.lineage = backbone["lineage"]

    # -------------------------------------------------------------------------
    # PASS 1: Backbone
    # Identify the lineage (or set of lineages) whose barcodes most closely
    # match each genomes observed substitutions.

    logger.info(str(datetime.now()) + "\tIdentifying backbone lineage.")
    for genome in genomes:

        genome.backbone = genome.identify_backbone(
            barcode_summary=genome.barcode_summary,
            barcodes_df=barcodes_df,
            tree=tree,
            recombinant_lineages=recombinant_lineages,
            include_recombinants=False,
        )

        # If no subs from the barcode is missing, this is not a recombinant
        num_missing_subs = len(genome.backbone["missing_subs"])
        if num_missing_subs == 0:
            genome.is_recombinant = False

    # -------------------------------------------------------------------------
    # Pass 2: Secondary Parent

    logger.info(str(datetime.now()) + "\tIdentifying secondary parent.")
    parents_dict = {}

    for genome in genomes:

        if genome.is_recombinant == False:
            continue

        max_lineages = genome.backbone["max_lineages"]
        barcode_summary_filter = genome.barcode_summary[
            ~genome.barcode_summary["lineage"].isin(max_lineages)
        ]
        alt_backbone = genome.identify_backbone(
            barcode_summary=barcode_summary_filter,
            barcodes_df=barcodes_df,
            tree=tree,
            recombinant_lineages=recombinant_lineages,
            include_recombinants=False,
        )

        result = genome.identify_breakpoints(alt_backbone)
        genome.recombination_subs = result["subs_df"]
        genome.recombination_breakpoints = result["breakpoints"]
        genome.recombination_regions = result["regions"]
        # parents = [r["parent"] for r in genome.recombination_regions.values()]
        # parents = ",".join(np.unique(sorted(parents)))
        # if parents not in parents_dict:
        #     parents_dict[parents] = []

        # parents_dict[parents].append(genome)

    # -------------------------------------------------------------------------
    # Pass 3: Export

    cols = ["strain", "lineage", "parents", "breakpoints", "regions", "subs"]

    export_dict = {col: [] for col in cols}

    for genome in genomes:
        if genome.is_recombinant == False:
            for col in cols:
                if col == "strain":
                    export_dict["strain"].append(genome.id)
                elif col == "lineage":
                    export_dict["lineage"].append(genome.lineage)
                else:
                    export_dict[col].append(NO_DATA_CHAR)
            continue

        regions = [
            "{}:{}|{}".format(r["start"], r["end"], r["parent"])
            for r in genome.recombination_regions.values()
        ]
        subs = [
            "{}|{}".format(",".join([str(s) for s in r["subs"]]), r["parent"])
            for r in genome.recombination_regions.values()
        ]
        parents = [r["parent"] for r in genome.recombination_regions.values()]

        export_dict["strain"].append(genome.id)
        export_dict["lineage"].append(genome.lineage)
        export_dict["parents"].append(",".join(parents))
        export_dict["regions"].append(",".join(regions))
        export_dict["breakpoints"].append(",".join(genome.recombination_breakpoints))
        export_dict["subs"].append(";".join(subs))

    export_df = pd.DataFrame(export_dict)
    file_name = "recombination.tsv"
    file_path = os.path.join(outdir, file_name)
    export_df.to_csv(file_path, sep="\t", index=False)
    # print(export_df)

    # for parent in parents_dict:
    #     print(parent)
    #     print(parents_dict[parent])
    quit()
    logger.info(str(datetime.now()) + "\tExporting results to: " + params.outdir)
    recombinants = [g for g in genomes]
    for genome in recombinants:

        file_name = "test.tsv"
        file_path = os.path.join(params.outdir, file_name)
        genome.recombination_subs.to_csv(file_path, sep="\t", index=False)

        if params.exclude_shared:
            genome.recombination_subs = genome.recombination_subs.query(
                "parent != 'shared'"
            )

        # Create the alignment of Ns + informative positions
        genome.write_snipit(exclude_shared=params.exclude_shared, outdir=params.outdir)
        break

    logger.info(str(datetime.now()) + "\tFinished module: recombination")
