# from rebar import wrappers, genome, edge_cases, utils, constants
# import os
# import pandas as pd
# from Bio import Phylo


# def test_prep_edge_cases(params):
#     """Prepare input for function edge_cases.handle_edge_cases."""

#     params.outdir = "test/tmp/edge_cases"
#     if not os.path.exists(params.outdir):
#         os.makedirs(params.outdir)

#     edge_case_recombinants = ",".join(constants.EDGE_CASE_RECOMBINANTS)
#     params.alignment = None
#     params.lineages = edge_case_recombinants
#     # Use the dataset created from test_wrappers
#     params.dataset = "test/tmp/wrappers/dataset/sars-cov-2-latest"
#     params.logger = utils.create_logger(os.path.join(params.outdir, "test.log"))

#     wrappers.run(params)


# # Note: I Think this needs to wait for a test_genome suite to be developed first.
# def test_handle_edge_cases(params, edge_cases_expected):
#     """Test function edge_cases.handle_edge_cases."""

#     subs_path = "test/tmp/edge_cases/subs.tsv"
#     subs_df = pd.read_csv(subs_path, sep="\t")

#     barcodes_path = "test/tmp/wrappers/dataset/sars-cov-2-latest/barcodes.tsv"
#     barcodes = pd.read_csv(barcodes_path, sep="\t")

#     diagnostic_path = "test/tmp/wrappers/dataset/sars-cov-2-latest/diagnostic.tsv"
#     diagnostic = pd.read_csv(diagnostic_path, sep="\t")

#     tree_path = "test/tmp/wrappers/dataset/sars-cov-2-latest/tree.nwk"
#     tree = Phylo.read(tree_path, "newick")
#     recombinant_tree = [c for c in tree.find_clades("X")][0]
#     recombinant_lineages = [c.name for c in recombinant_tree.find_clades()]

#     lineage_to_clade_path = (
#         "test/tmp/wrappers/dataset/sars-cov-2-latest/lineage_to_clade.tsv"
#     )
#     lineage_to_clade = pd.read_csv(lineage_to_clade_path, sep="\t")

#     edge_case_lineages = list(subs_df["strain"])
#     edge_case_recombinants = []

#     # Debugging
#     #edge_case_lineages = ["XP"]

#     conflict = False

#     for lineage in edge_case_lineages:
#         # Convert type from dataframe to series for func
#         subs_row = subs_df[subs_df["strain"] == lineage].squeeze()

#         lineage_genome = genome.Genome(
#             subs_row=subs_row,
#             barcodes=barcodes,
#             diagnostic=diagnostic,
#             lineage_to_clade=lineage_to_clade,
#             tree=tree,
#             recombinant_tree=recombinant_tree,
#             recombinant_lineages=recombinant_lineages,
#         )
#         result = edge_cases.handle_edge_cases(
#             genome=lineage_genome,
#             barcode_summary=lineage_genome.barcode_summary,
#             tree=tree,
#             min_subs=params.min_subs,
#             min_length=params.min_length,
#             min_consecutive=params.min_consecutive,
#         )

#         recombinant = lineage_genome.lineage.recombinant

#         if recombinant is False:
#             conflict = True
#             print("lineage: {}, no recombination detected.".format(lineage))
#             continue

#         if recombinant not in edge_case_recombinants:
#             edge_case_recombinants.append(recombinant)

#         if recombinant not in edge_cases_expected:
#             conflict = True
#             print(
#                 "lineage: {}, recombinant: {} has no expected values.".format(
#                     lineage, recombinant
#                 )
#             )

#         else:
#             for param in edge_cases_expected[recombinant]:
#                 expected = edge_cases_expected[recombinant][param]
#                 actual = result[param]

#                 if expected != actual:
#                     conflict = True
#                     print(
#                         "lineage: {}, param: {}, expected: {}, actual: {}".format(
#                             lineage, param, expected, actual
#                         )
#                     )

#     for recombinant in edge_cases_expected:
#         if recombinant not in edge_case_recombinants:
#             print("recombinant: {} has no actual values.".format(recombinant))
#             conflict = True

#     assert conflict is False
