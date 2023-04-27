import pytest
from rebar import utils
from .constants import (
    MASK,
    MAX_DEPTH,
    MIN_LENGTH,
    MIN_SUBS,
    MIN_CONSECUTIVE,
    MAX_BREAKPOINTS,
    PLOT_FORMAT,
)


@pytest.fixture(scope="module")
def params():
    """Return a NameSpace of generic testing params."""

    params = utils.Namespace(
        outdir="",
        logger=utils.create_logger(),
        threads=1,
        debug=True,
        reference="test/tmp/download_reference_sequence/reference.fasta",
        alignment="test/tmp/download_consensus_sequences/alignment.fasta",
        subs="test/tmp/parse_alignment/subs.tsv",
        tree="test/tmp/create_tree/tree.nwk",
        barcodes="test/tmp/create_barcodes/barcodes.tsv",
        lineage_to_clade="test/tmp/create_barcodes/lineage_to_clade.tsv",
        mask=MASK,
        max_depth=MAX_DEPTH,
        min_length=MIN_LENGTH,
        min_subs=MIN_SUBS,
        min_consecutive=MIN_CONSECUTIVE,
        max_breakpoints=MAX_BREAKPOINTS,
        plot_format=PLOT_FORMAT,
        edge_cases=True,
        output_all=True,
        output_fasta=True,
        output_tsv=True,
        output_barcode=True,
        output_plot=True,
        output_yaml=True,
        exclude_non_recomb=False,
        shared=False,
        # dataset
        name="sars-cov-2",
        tag="latest",
        # run
        lineages="AY.4,BA.5.2,XD,XBB.1.5.1,XBL",
        dataset="test/tmp/dataset/sars-cov-2-latest",
    )
    return params
