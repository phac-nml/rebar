import sys
import logging
import os
import requests
import zstandard as zstd

NO_DATA_CHAR = "NA"

LINEAGES_URL = (
    "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/"
    "lineage_notes.txt"
)

BARCODES_USHER_URL = (
    "https://github.com/andersen-lab/Freyja-data/raw/main/usher_barcodes.csv"
)
BARCODES_NEXTCLADE_URL = (
    "https://raw.githubusercontent.com/corneliusroemer/pango-sequences/"
    "main/data/pango-consensus-sequences_summary.json"
)

PANGO_SEQUENCES_URL = (
    "https://raw.githubusercontent.com/corneliusroemer/pango-sequences/"
    "main/data/pango-consensus-sequences_genome-nuc.fasta.zst"
)

PROBLEMATIC_LINEAGES = ["BA.2.85"]


def create_logger(logfile=None):
    # create logger
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


def download_consensus_sequences(outdir):
    file_name = PANGO_SEQUENCES_URL.split("/")[-1]
    # The sequences are a .fasta.zst file, remove the .zst as we're decompressing
    fasta_path = os.path.join(outdir, os.path.splitext(file_name)[0])
    response = requests.get(PANGO_SEQUENCES_URL, stream=True)

    decomp = zstd.ZstdDecompressor()
    with open(fasta_path, "wb") as outfile:
        decomp.copy_stream(response.raw, outfile)

    return fasta_path
