import sys
import logging

NO_DATA_CHAR = "NA"

LINEAGES_URL = (
    "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/"
    "lineage_notes.txt"
)

BARCODES_URL = "https://github.com/andersen-lab/Freyja-data/raw/main/usher_barcodes.csv"


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
