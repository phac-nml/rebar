import sys
import logging
import os
from datetime import datetime
from io import StringIO
import re
import requests
import urllib
import functools

import zstandard as zstd
import pandas as pd
from pango_aliasor.aliasor import Aliasor
from Bio import Phylo, Entrez, SeqIO
from Bio.Phylo.BaseTree import Clade
from multiprocess import Pool  # Note that we are importing "multiprocess", no "ing"!
from tqdm import tqdm

from .genome import genome_mp
from .substitution import Substitution
from .constants import (
    NO_DATA_CHAR,
    PANGO_SEQUENCES_URL,
    BARCODES_NEXTCLADE_URL,
    BARCODES_USHER_URL,
    BARCODE_MANUAL_EDITS,
    PROBLEMATIC_LINEAGES,
    LINEAGE_SUMMARY_URL,
)
from .export import Export

# -----------------------------------------------------------------------------
# Classes
# -----------------------------------------------------------------------------


class Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------


def create_logger(logfile=None):
    """
    Create logging object for help messages.

    Parameters
    ----------
        logfile : str
            file path to write log to.
    """
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    # create file handler which logs even debug messages
    if logfile:
        handler = logging.FileHandler(logfile, "w+")
    else:
        handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    return logger


def url_header_info(url):

    info = {}

    # Download the file  and parse info from header
    url_handle = urllib.request.urlopen(url)
    headers = url_handle.info()

    # Date format: Wed, 19 Apr 2023 16:19:59 GMT
    file_date_str = " ".join(headers["date"].split(" ")[1:4])
    file_date = datetime.strptime(file_date_str, "%d %b %Y").date()

    info["url"] = url
    info["date"] = str(file_date)
    info["etag"] = headers["etag"].replace('"', "")

    return info


def download_reference_sequence(params, accession):
    """
    Download reference sequence from genbank.

    Parameters
    ----------
        accession : str
            Genbank nucleotide accession of reference genome.
        params.logger : logging.RootLogger
            logging object for messages.
        params.outdir : str
            output directory to write fasta sequence to.
    """
    logger = params.logger

    # object to hold information about downloaded files
    info = {}

    logger.info(str(datetime.now()) + "\t" + "-" * 40)
    logger.info(str(datetime.now()) + "\tDownloading reference.")
    Entrez.email = "Your.Name.Here@example.org"
    handle = Entrez.efetch(
        db="nucleotide", id=accession, rettype="fasta", retmode="text"
    )
    record = SeqIO.read(handle, "fasta")

    # Export
    file_name = "reference.fasta"
    file_path = os.path.join(params.outdir, file_name)
    logger.info(str(datetime.now()) + "\tExporting reference: " + file_path)
    SeqIO.write(record, file_path, "fasta")

    # Update info
    info["accession"] = accession
    info["file"] = file_path

    # Finish
    logger.info(str(datetime.now()) + "\tFinished downloading reference.")

    return info


def download_consensus_sequences(params):
    """
    Download consensus sequences of designated sars-cov-2 lineages.
    Sources:
      - github.com/corneliusroemer/pango-sequences

    Parameters
    ----------
        logger : logging.RootLogger
            logging object for messages
        outdir : str
            output directory to write fasta sequences to.
    """
    logger = params.logger

    # object to hold information about downloaded files
    info = {}

    # The sequences are a .fasta.zst file, remove the .zst as we're decompressing
    fasta_name = "sequences.fasta"
    fasta_path = os.path.join(params.outdir, fasta_name)

    logger.info(str(datetime.now()) + "\t" + "-" * 40)
    logger.info(str(datetime.now()) + "\tDownloading lineage sequences.")

    # Summarize url file info
    info = url_header_info(PANGO_SEQUENCES_URL)
    info["file"] = fasta_name

    response = requests.get(PANGO_SEQUENCES_URL, stream=True)

    # Decompress the zstd format
    decomp = zstd.ZstdDecompressor()
    # Write decompressed contents to file (tmp)
    with open(fasta_path, "wb") as outfile:
        decomp.copy_stream(response.raw, outfile)

    records = SeqIO.parse(fasta_path, "fasta")
    logger.info(str(datetime.now()) + "\tApplying edge-case curation.")
    fasta_lines = []

    for record in records:
        for lineage in BARCODE_MANUAL_EDITS:
            if record.id != lineage:
                continue
            for sub in BARCODE_MANUAL_EDITS[lineage]:
                logger.info(
                    str(datetime.now())
                    + "\t\tAdding "
                    + lineage
                    + " barcode "
                    + str(sub)
                )
                sub = Substitution(sub)
                # genome coordinates are 1 based
                sub_i = sub.coord - 1
                record.seq = record.seq[:sub_i] + sub.alt + record.seq[sub_i + 1 :]
        fasta_lines.append(">" + str(record.id))
        fasta_lines.append(str(record.seq))

    logger.info(str(datetime.now()) + "\tExported lineage sequences: " + fasta_path)
    with open(fasta_path, "w") as outfile:
        outfile.write("\n".join(fasta_lines) + "\n")

    # Finish
    logger.info(str(datetime.now()) + "\tFinished downloading lineage sequences.")

    return info


def create_barcodes(params):
    """
    Create csv of lineage barcodes from nextclade and usher.
    Sources:
      - github.com/corneliusroemer/pango-sequences
      - github.com/andersen-lab/Freyja-data

    Parameters
    ----------
        logger : logging.RootLogger
            logging object for messages
        output : str
            file path for output barcodes csv.
    """

    info = {}

    logger = params.logger
    logger.info(str(datetime.now()) + "\t" + "-" * 40)
    logger.info(str(datetime.now()) + "\tCreating barcodes.")

    file_name = "barcodes.tsv"
    # -------------------------------------------------------------------------
    # Nextclade barcodes

    logger.info(str(datetime.now()) + "\tDownloading Nextclade barcodes.")

    # Summarize url file info
    info = url_header_info(BARCODES_NEXTCLADE_URL)
    info["url_nextclade"] = info["url"]

    r = requests.get(BARCODES_NEXTCLADE_URL)
    barcodes_data = r.json()
    barcodes_dict = {
        lineage: barcodes_data[lineage]["nucSubstitutions"] for lineage in barcodes_data
    }
    # Mapping of lineage to clade information
    lineage_to_clade = {
        lineage: barcodes_data[lineage]["nextstrainClade"] for lineage in barcodes_data
    }
    lineage_to_clade_df = pd.DataFrame(
        {
            "lineage": list(lineage_to_clade.keys()),
            "nextstrainClade": list(lineage_to_clade.values()),
        }
    )

    # -------------------------------------------------------------------------
    # UShER Barcodes

    logger.info(str(datetime.now()) + "\tDownloading UShER barcodes.")

    # Summarize url file info
    info = url_header_info(BARCODES_USHER_URL)
    info["url_usher"] = info["url"]
    info["file"] = file_name

    r = requests.get(BARCODES_USHER_URL)
    barcodes_text = r.text
    barcodes_usher_df = pd.read_csv(StringIO(barcodes_text), sep=",")
    # Rename the empty first column that should be lineage
    barcodes_usher_df.rename(columns={"Unnamed: 0": "lineage"}, inplace=True)

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

    logger.info(str(datetime.now()) + "\tApplying edge-case curation.")

    for lineage in BARCODE_MANUAL_EDITS:
        lineage_i = barcodes_nextclade_df[
            barcodes_nextclade_df["lineage"] == lineage
        ].index.values[0]

        for sub, value in BARCODE_MANUAL_EDITS[lineage].items():
            logger.info(
                str(datetime.now())
                + "\t\tAdding "
                + lineage
                + " barcode "
                + sub
                + "="
                + str(value)
            )
            if sub not in barcodes_nextclade_df.columns:
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
    barcodes_path = os.path.join(params.outdir, file_name)
    logger.info(str(datetime.now()) + "\tExporting barcodes: " + barcodes_path)
    # Export to csv, to be consistent with usher barcodes format
    barcodes_df.to_csv(barcodes_path, sep="\t", index=False)

    clade_path = os.path.join(params.outdir, "lineage_to_clade.tsv")
    logger.info(
        str(datetime.now()) + "\tExporting lineage to clade mapping: " + clade_path
    )
    lineage_to_clade_df.to_csv(clade_path, sep="\t", index=False)

    # Finish
    logger.info(str(datetime.now()) + "\tFinished creating barcodes.")
    return info


def create_tree(params):
    """
    Create nomenclature tree of designated pango lineages.

    Parameters
    ----------
        logger : logging.RootLogger
            logging object for messages
        output : str
            file path for output newick tree.
    """

    info = {}

    logger = params.logger
    logger.info(str(datetime.now()) + "\t" + "-" * 40)
    logger.info(str(datetime.now()) + "\tCreating tree.")
    file_name = "tree.nwk"

    # -------------------------------------------------------------------------
    # Download latest designated lineages from pango-designation

    logger.info(str(datetime.now()) + "\tDownloading designated lineage summaries.")

    # Summarize url file info
    info = url_header_info(LINEAGE_SUMMARY_URL)
    info["file"] = file_name

    r = requests.get(LINEAGE_SUMMARY_URL)
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

    tree_path = os.path.join(params.outdir, file_name)
    logger.info(str(datetime.now()) + "\tExporting newick tree: " + tree_path)
    Phylo.write(tree, tree_path, "newick")

    # Finish
    logger.info(str(datetime.now()) + "\tFinished creating tree.")
    return info


def parse_alignment(params):
    """
    Parse alignment for substitutions, deletions, and missing data.

    Parameters
    ----------
        reference : str
            file path to reference genome.
        alignment : str
            file_path to alignment.
        mask : int
            number of bases to mask at 5' and 3' end.
        logger : logging.RootLogger
            logging object for messages
        threads : int
            number of CPUs to use.
        outdir : str
            directory path for output files.
    """

    logger = params.logger
    logger.info(str(datetime.now()) + "\t" + "-" * 40)
    logger.info(str(datetime.now()) + "\tParsing substitutions from alignment.")

    # Import reference
    logger.info(str(datetime.now()) + "\tImporting reference: " + params.reference)
    records = SeqIO.parse(params.reference, "fasta")
    ref_rec = next(records)

    # Import alignment
    logger.info(str(datetime.now()) + "\tImporting alignment: " + params.alignment)
    num_records = len(list(SeqIO.parse(params.alignment, "fasta")))
    records = SeqIO.parse(params.alignment, "fasta")

    # Parse substitutions
    # Process genomes in parallel, `genome_mp` is a multiprocessing wrapper
    # function for the `Genome` class.
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
    task_progress = tqdm(pool.imap(task, iterator), total=total)
    task_description = (
        str(datetime.now()) + "      Parsing substitutions from alignment"
    )
    task_progress.set_description(task_description, refresh=True)
    genomes = list(task_progress)

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
    logger.info(str(datetime.now()) + "\tFinished parsing substitutions.")
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
        lineage_to_clade : str
            file path of mapping lineages to clades.
        logger : logging.RootLogger
            logging object for messages.
        threads : int
            number of CPUs to use.
        outdir : str
            directory path for output files.
    """

    logger = params.logger
    logger.info(str(datetime.now()) + "\t" + "-" * 40)
    logger.info(str(datetime.now()) + "\tDetecting recombination.")

    # Import the dataframe from the `subs` module, or alternatively from nextclade
    logger.info(str(datetime.now()) + "\tImporting substitutions: " + params.subs)
    subs_df = pd.read_csv(params.subs, sep="\t").fillna(NO_DATA_CHAR)
    subs_df.set_index("strain", inplace=True)
    subs_df["strain"] = subs_df.index

    # Import lineage barcodes
    logger.info(str(datetime.now()) + "\tImporting barcodes: " + params.barcodes)
    barcodes_df = pd.read_csv(params.barcodes, sep="\t")

    # Import tree
    logger.info(str(datetime.now()) + "\tImporting tree: " + params.tree)
    tree = Phylo.read(params.tree, "newick")

    # Import mapping of lineages to clades
    logger.info(
        str(datetime.now())
        + "\tImporting lineage to clade mapping: "
        + params.lineage_to_clade
    )
    lineage_to_clade = pd.read_csv(params.lineage_to_clade, sep="\t")

    # Identify which lineages are known recombinants
    # ie. descended from the "X" recombinant MRCA node
    recombinant_tree = [c for c in tree.find_clades("X")][0]
    recombinant_lineages = [c.name for c in recombinant_tree.find_clades()]

    # Detect recombination in samples.
    # Process genomes in parallel, `genome_mp` is a multiprocessing wrapper
    # function for the `Genome` class.
    pool = Pool(params.threads)
    iterator = subs_df.iterrows()
    total = len(subs_df)

    # Debugging
    # iterator = [rec for rec in subs_df.iterrows() if rec[1]["strain"].startswith("X")]
    # total = len(iterator)

    task = functools.partial(
        genome_mp,
        debug=params.debug,
        logger=params.logger,
        barcodes=barcodes_df,
        tree=tree,
        recombinant_tree=recombinant_tree,
        recombinant_lineages=recombinant_lineages,
        lineage_to_clade=lineage_to_clade,
        max_depth=params.max_depth,
        min_subs=params.min_subs,
        min_consecutive=params.min_consecutive,
        min_length=params.min_length,
        edge_cases=params.edge_cases,
    )

    task_progress = tqdm(pool.imap(task, iterator), total=total)
    task_description = str(datetime.now()) + "      Detecting recombination"
    task_progress.set_description(task_description, refresh=True)
    genomes = list(task_progress)

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
    if params.output_all or params.output_yaml:
        outpath = os.path.join(params.outdir, "summary.yaml")
        logger.info(str(datetime.now()) + "\tExporting YAML: " + outpath)
        export.to_yaml()

    if params.output_all or params.output_tsv:
        outpath = os.path.join(params.outdir, "summary.tsv")
        logger.info(str(datetime.now()) + "\tExporting TSV: " + outpath)
        export.to_dataframe()

    if params.output_all or params.output_barcode:
        outpath = os.path.join(params.outdir, "barcode_<parent_1>_<parent_2>.tsv")
        logger.info(str(datetime.now()) + "\tExporting barcode mutations: " + outpath)
        export.to_barcodes()

    # Dataframe
    if params.output_all or params.output_fasta:
        outpath = os.path.join(params.outdir, "alignment_<parent_1>_<parent_2>.fasta")
        logger.info(str(datetime.now()) + "\tExporting FASTA alignments: " + outpath)
        export.to_alignments()

    # Snipit plots
    if params.output_all or params.output_plot:
        outpath = os.path.join(
            params.outdir, "snipit_<parent_1>_<parent_2>." + params.snipit_format
        )
        logger.info(str(datetime.now()) + "\tExporting snipit plots: " + outpath)
        export.to_snipit(ext=params.snipit_format)

    # Finish
    logger.info(str(datetime.now()) + "\tFinished detecting recombination.")
    return 0
