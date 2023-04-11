import os
from datetime import datetime
from multiprocess import Pool  # Note that we are importing "multiprocess", no "ing"!
import functools
import requests
import logging
import re
from io import StringIO
import pandas as pd
from tqdm import tqdm
from .utils import create_logger, NO_DATA_CHAR, GENOME_LEN
from .genome import Genome, GenomeAlt
from .backbone import Backbone
from pango_aliasor.aliasor import Aliasor
from Bio import Phylo, SeqIO
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
    logger.info(str(datetime.now()) + "\tBeginning module: tree.")

    # Create output directory if it doesn't exist
    outdir = os.path.dirname(params.output)
    if not os.path.exists(outdir) and outdir != "":
        logger.info(str(datetime.now()) + "\tCreating output directory: " + outdir)
        os.mkdir(outdir)

    # reference
    logger.info(str(datetime.now()) + "\tImporting reference: " + params.reference)
    records = SeqIO.parse(params.reference, "fasta")
    ref_rec = next(records)

    # alignment
    logger.info(str(datetime.now()) + "\tImporting alignment: " + params.alignment)
    num_records = len(list(SeqIO.parse(params.alignment, "fasta")))
    records = SeqIO.parse(params.alignment, "fasta")

    logger.info(str(datetime.now()) + "\tParsing sequences.")
    p = Pool(params.threads)

    # Use starmap for positional parameters
    # genomes = p.starmap(GenomeAlt, zip(records, repeat(ref_rec)))

    # Use imap + functools.partial for named parameters
    # genomes = p.imap(functools.partial(GenomeAlt, reference=ref_rec), records)

    # Reminder: An enclosing list() statement waits for the iterator to end.
    genomes = list(
        tqdm(
            p.imap(functools.partial(GenomeAlt, reference=ref_rec), records),
            total=num_records,
        )
    )

    logger.info(str(datetime.now()) + "\tExporting results to: " + params.output)

    cols = ["strain", "substitutions", "deletions", "missing"]
    export_data = {col: [] for col in cols}

    dfs = [genome.to_dataframe() for genome in genomes]
    df = pd.concat(dfs)
    df.to_csv(params.output, sep="\t", index=False)

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
    # PHASE 1: IDENTIFY NON-RECOMBINANTS
    # Identify the lineage (or set of lineages) whose barcodes most closely
    # match each genomes observed substitutions.

    logger.info(str(datetime.now()) + "\tIdentifying non-recombinants.")
    for genome in genomes:

        genome.backbone.search(
            genome=genome,
            barcode_summary=genome.barcode_summary,
            barcodes_df=barcodes_df,
            tree=tree,
            recombinant_lineages=recombinant_lineages,
            include_recombinants=True,
        )
        genome.lineage = genome.backbone.lineage

        # Option 1: Top backbone lineage is recombinant
        #           Investigate further in next phase
        if genome.backbone.lineage in recombinant_lineages:

            # Identify the generic recombinant type
            recombinant_path = recombinant_tree.get_path(genome.backbone.lineage)
            # Move backwards up the path, until we find a parental lineage that
            # starts with "X", because it might be an alias ("EK")

            for c in recombinant_path[::-1]:
                if c.name.startswith("X"):
                    genome.recombinant = c.name.split(".")[0]
                    break

            # Check if this is a recursive recombinant
            # In the get_path method, the root is excluded
            node_path = recombinant_tree.get_path(genome.recombinant)
            # So if node_path just has more than one clade, it's recursive
            if len(node_path) > 1:
                genome.recursive = True

        # Option 2: Perfect match to non-recombinant
        #           Stop further investigation
        elif len(genome.backbone.conflict_subs_ref) == 0:
            genome.recombinant = False

        # else. Has conflicts with non-recombinant lineage
        #           Investigate further in next phase

    quit()

    # -------------------------------------------------------------------------
    # PHASE 2: BACKBONE SEARCH
    # Identify the lineage (or set of lineages) whose barcodes most closely
    # match each genomes observed substitutions, EXCLUDING RECOMBINANTS

    logger.info(str(datetime.now()) + "\tIdentifying recombinant backbones")
    for genome in genomes:

        # Skip clear non-recombinants
        if genome.recombinant == False:
            continue

        # Redo the backbone search, excluding select lineages
        barcode_summary = genome.barcode_summary

        # Option 1. Not Recursive, exclude recombinant lineages from search
        if not genome.recursive:
            barcode_summary = barcode_summary[
                ~barcode_summary["lineage"].isin(recombinant_lineages)
            ]
        # Uhhh, TBD
        else:
            ...

        genome.backbone = Backbone()
        genome.backbone.search(
            genome=genome,
            barcode_summary=barcode_summary,
            barcodes_df=barcodes_df,
            tree=tree,
            recombinant_lineages=recombinant_lineages,
            include_recombinants=False,
        )

    # -------------------------------------------------------------------------
    # Pass 2: Secondary Parent

    logger.info(str(datetime.now()) + "\tIdentifying secondary parent.")
    parents_dict = {}
    max_depth = 3

    for genome in genomes:

        if genome.recombinant == False:
            continue

        print(genome)
        barcode_summary_filter = genome.barcode_summary[
            ~genome.barcode_summary["lineage"].isin(genome.backbone.max_lineages)
        ]
        alt_backbone = Backbone()
        alt_backbone.search(
            genome=genome,
            barcode_summary=barcode_summary_filter,
            barcodes_df=barcodes_df,
            tree=tree,
            recombinant_lineages=recombinant_lineages,
            include_recombinants=False,
        )

        result = genome.identify_breakpoints(alt_backbone)
        print(result)
        break
        genome.recombination_subs = result["subs_df"]
        genome.recombination_breakpoints = result["breakpoints"]
        genome.recombination_regions = result["regions"]
        # parents = [r["parent"] for r in genome.recombination_regions.values()]
        # parents = ",".join(np.unique(sorted(parents)))
        # if parents not in parents_dict:
        #     parents_dict[parents] = []

        # parents_dict[parents].append(genome)

    quit()

    # -------------------------------------------------------------------------
    # Pass 3: Export

    cols = [
        "strain",
        "is_recombinant",
        "lineage",
        "parents",
        "breakpoints",
        "regions",
        "subs",
    ]

    export_dict = {col: [] for col in cols}

    for genome in genomes:
        if genome.is_recombinant == False:
            for col in cols:
                if col == "strain":
                    export_dict["strain"].append(genome.id)
                elif col == "lineage":
                    export_dict["lineage"].append(genome.backbone.lineage)
                elif col == "is_recombinant":
                    export_dict["is_recombinant"].append(genome.is_recombinant)
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
