from rebar import wrappers, genome, edge_cases, utils
import os
import pandas as pd
from Bio import Phylo


def test_prep_edge_cases(params):
    """Prepare input for function edge_cases.handle_edge_cases."""

    params.outdir = "test/tmp/edge_cases"
    if not os.path.exists(params.outdir):
        os.makedirs(params.outdir)

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
    params.alignment = None
    params.lineages = edge_case_lineages
    # Use the dataset created from test_wrappers
    params.dataset = "test/tmp/wrappers/dataset/sars-cov-2-latest"
    params.logger = utils.create_logger(os.path.join(params.outdir, "test.log"))

    wrappers.run(params)


def test_handle_edge_cases(params, edge_cases_expected):
    """Test function edge_cases.handle_edge_cases."""

    subs_path = "test/tmp/edge_cases/subs.tsv"
    subs_df = pd.read_csv(subs_path, sep="\t")

    barcodes_path = "test/tmp/wrappers/dataset/sars-cov-2-latest/barcodes.tsv"
    barcodes = pd.read_csv(barcodes_path, sep="\t")

    tree_path = "test/tmp/wrappers/dataset/sars-cov-2-latest/tree.nwk"
    tree = Phylo.read(tree_path, "newick")
    recombinant_tree = [c for c in tree.find_clades("X")][0]
    recombinant_lineages = [c.name for c in recombinant_tree.find_clades()]

    lineage_to_clade_path = (
        "test/tmp/wrappers/dataset/sars-cov-2-latest/lineage_to_clade.tsv"
    )
    lineage_to_clade = pd.read_csv(lineage_to_clade_path, sep="\t")

    edge_case_lineages = list(subs_df["strain"])

    conflict = False

    for lineage in edge_case_lineages:
        # Convert type from dataframe to series for func
        subs_row = subs_df[subs_df["strain"] == lineage].squeeze()

        lineage_genome = genome.Genome(
            subs_row=subs_row,
            barcodes=barcodes,
            lineage_to_clade=lineage_to_clade,
            tree=tree,
            recombinant_tree=recombinant_tree,
            recombinant_lineages=recombinant_lineages,
        )
        result = edge_cases.handle_edge_cases(
            genome=lineage_genome,
            barcode_summary=lineage_genome.barcode_summary,
            tree=tree,
            min_subs=params.min_subs,
            min_length=params.min_length,
            min_consecutive=params.min_consecutive,
        )

        if lineage not in edge_cases_expected:
            conflict = True
            print("lineage: {} has no expected values.".format(lineage))

        else:
            for param in edge_cases_expected[lineage]:
                expected = edge_cases_expected[lineage][param]
                actual = result[param]

                if expected != actual:
                    conflict = True
                    print(
                        "lineage: {}, param: {}, expected: {}, actual: {}".format(
                            lineage, param, expected, actual
                        )
                    )

    for lineage in edge_cases_expected:
        if lineage not in edge_case_lineages:
            print("lineage: {} has no actual values.".format(lineage))
            conflict = True

    assert conflict is False
