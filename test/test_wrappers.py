from rebar import utils, wrappers
import os


def test_dataset_sarscov2_latest(params):
    """Test function wrappers.dataset."""
    params.outdir = "test/tmp/wrappers/dataset/sars-cov-2-latest"
    if not os.path.exists(params.outdir):
        os.makedirs(params.outdir)
    params.logger = utils.create_logger(os.path.join(params.outdir, "test.log"))

    wrappers.dataset(params)


def test_run_sarscov2_latest(params):
    """Test function wrappers.run."""
    params.outdir = "test/tmp/wrappers/run"
    if not os.path.exists(params.outdir):
        os.makedirs(params.outdir)
    # Disable alignment, we'll use lineage list from conftest
    params.alignment = None
    params.dataset = "test/tmp/wrappers/dataset/sars-cov-2-latest"
    params.logger = utils.create_logger(os.path.join(params.outdir, "test.log"))

    wrappers.run(params)


# def test_no_lineages_or_alignment(params):
#     """Test function wrappers.run with missing lineages and alignment."""
#     params.outdir = "test/tmp/wrappers/no_lineages_or_alignment"
#     if not os.path.exists(params.outdir):
#         os.makedirs(params.outdir)
#     params.logger = utils.create_logger(os.path.join(params.outdir, "test.log"))

#     params.alignment = None
#     params.lineages = None
#     wrappers.run(params)
