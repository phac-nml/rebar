from rebar import edge_cases, utils, wrappers
import os


def test_handle_edge_cases(params):
    """Test function edge_cases.handle_edge_cases."""

    # Parse edge_cases file as plain text to get edge case lineages
    file_path = edge_cases.__file__
    with open(file_path) as infile:
        file_lines = infile.readlines()

    # Search for the str patter, example:
    #   'genome.lineage.recombinant in ["XBK", "XBQ"]'
    edge_case_lineages = []

    for line in file_lines:
        if "recombinant in [" in line:
            split_line = line.split('["')[1]
            split_line = split_line.split('"]')[0]

            clean_line = split_line.replace('"', "").replace(" ", "")
            for lineage in clean_line.split(","):
                edge_case_lineages.append(lineage)

    edge_case_lineages = ",".join(edge_case_lineages)

    # Now, run with the specified edge case lineages
    params.outdir = "test/tmp/edge_cases"
    if not os.path.exists(params.outdir):
        os.makedirs(params.outdir)

    params.alignment = None
    params.lineages = edge_case_lineages
    # Use the dataset created from test_wrappers
    params.dataset = "test/tmp/wrappers/dataset/sars-cov-2-latest"
    params.logger = utils.create_logger(os.path.join(params.outdir, "test.log"))

    wrappers.run(params)
