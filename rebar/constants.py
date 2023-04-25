# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------

NO_DATA_CHAR = "NA"

LINEAGE_SUMMARY_URL = (
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

# Lineages that produce undesireable results across all recombinants.
PROBLEMATIC_LINEAGES = ["BA.2.85"]

BARCODE_MANUAL_EDITS = {
    # This substitution is found in the UShER tree, but missing from Nextclade
    "XAE": {"A26530G": 1}
}
