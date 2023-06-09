# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------

NO_DATA_CHAR = "NA"

LINEAGE_SUMMARY_URL = (
    "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/"
    "lineage_notes.txt"
)

ALIAS_KEY_URL = (
    "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/"
    "pango_designation/alias_key.json"
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

# Lineages that produce undesireable results across all recombinants.
# PROBLEMATIC_LINEAGES = ["BA.2.85"]
PROBLEMATIC_LINEAGES = []

BARCODE_MANUAL_EDITS = {
    # This substitution is found in the UShER tree, but missing from Nextclade
    "XAE": {"A26530G": 1}
}

MASK = 200
MAX_DEPTH = 3
MIN_LENGTH = 500
MIN_SUBS = 1
MIN_CONSECUTIVE = 3
MAX_BREAKPOINTS = 10
PLOT_EXT = "png"

# These are known edge case recombinants, which generally required
# more relaxed parameters.
EDGE_CASE_RECOMBINANTS = [
    "XB",
    "XP",
    "XR",
    "XAD",
    "XAE",
    "XAJ",
    "XAS",
    "XAV",
    "XAY",
    "XAZ",
    "XBC",
    "XBK",
    "XBQ",
    "XBZ",
]
