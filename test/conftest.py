import pytest
from rebar import utils


@pytest.fixture(scope="module")
def params():
    """Return a NameSpace of generic testing params."""

    params = utils.Namespace(
        outdir="",
        logger=utils.create_logger(),
        threads=1,
        debug=True,
        reference="test/tmp/download_reference_sequence/reference.fasta",
        alignment="test/tmp/download_consensus_sequences/sequences.fasta",
        subs="test/tmp/parse_alignment/subs.tsv",
        tree="test/tmp/create_tree/tree.nwk",
        barcodes="test/tmp/create_barcodes/barcodes.tsv",
        lineage_to_clade="test/tmp/create_barcodes/lineage_to_clade.tsv",
        mask=1,
        max_depth=3,
        min_length=500,
        min_subs=1,
        min_consecutive=3,
        edge_cases=True,
        output_all=True,
        output_fasta=True,
        output_tsv=True,
        output_barcode=True,
        output_plot=True,
        output_yaml=True,
        exclude_non_recomb=False,
        shared=False,
        snipit_format="png",
        # dataset
        name="sars-cov-2",
        tag="latest",
        # run
        lineages="AY.4,BA.5.2,XD,XBB.1.5.1,XBL",
        dataset="test/tmp/dataset/sars-cov-2-latest",
    )
    return params
