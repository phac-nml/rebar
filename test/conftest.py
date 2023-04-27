import os
import pytest

from rebar import utils
from rebar.constants import (
    MASK,
    MAX_DEPTH,
    MIN_LENGTH,
    MIN_SUBS,
    MIN_CONSECUTIVE,
    MAX_BREAKPOINTS,
    PLOT_EXT,
)


@pytest.fixture(scope="module")
def params():
    """Return a NameSpace of generic testing params."""

    dataset = "test/tmp/dataset/sars-cov-2-latest"
    name = "sars-cov-2"
    tag = "latest"

    params = utils.Namespace(
        outdir="",
        logger=utils.create_logger(),
        threads=1,
        debug=True,
        dataset=dataset,
        reference=os.path.join(dataset, "reference.fasta"),
        alignment=os.path.join(dataset, "alignment.fasta"),
        tree=os.path.join(dataset, "tree.nwk"),
        barcodes=os.path.join(dataset, "barcodes.tsv"),
        lineage_to_clade=os.path.join(dataset, "lineage_to_clade.tsv"),
        subs="test/tmp/parse_alignment/subs.tsv",
        mask=MASK,
        max_depth=MAX_DEPTH,
        min_length=MIN_LENGTH,
        min_subs=MIN_SUBS,
        min_consecutive=MIN_CONSECUTIVE,
        max_breakpoints=MAX_BREAKPOINTS,
        plot_ext=PLOT_EXT,
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
        name=name,
        tag=tag,
        # run
        lineages="AY.4,BA.5.2,XD,XBB.1.5.1,XBL",
        validate=True,
    )
    return params
