# The function `handle_edge_cases` is called in `genome.py` at Line 610.
# The currently implemented params that can be customized for a sample are:
#   1. min_consecutive
#   2. min_subs
#   3. min_length
#   4. barcode_summary: summary of barcode matches to find parent_1
# More parameters could be added as needed, with the correponding updates
# added to `genome.py`


def handle_edge_cases(
    genome, barcode_summary, tree, min_subs, min_length, min_consecutive
):
    """
    Edge case handling for recombinants with few mutations.
    """

    # These are the parameters we might need to adjust in an edge case
    result = {
        "min_consecutive": min_consecutive,
        "min_subs": min_subs,
        "min_length": min_length,
        "barcode_summary": barcode_summary,
    }

    # ---------------------------------------------------------------------
    # XB: top_lineages are tied exactly B.1.631 and B.1.634
    #     force the first parent to be B.1.631
    #     there is maybe a small second breakpoint (~100 nuc)
    if genome.lineage.recombinant in ["XB"]:
        genome.lineage.edge_case = True
        include_tree = next(tree.find_clades("B.1.631"))
        include_descendants = [c.name for c in include_tree.find_clades()]
        result["barcode_summary"] = barcode_summary[
            barcode_summary["lineage"].isin(include_descendants)
        ]
        result["min_length"] = 100

    # ---------------------------------------------------------------------
    # XP: second parent (BA.2) comes from only one barcode position: A29510C
    #     force first parent to be BA.2 and relax region/subs lengths
    elif genome.lineage.recombinant in ["XP"]:
        genome.lineage.edge_case = True
        include_tree = next(tree.find_clades("BA.2"))
        include_descendants = [c.name for c in include_tree.find_clades()]
        result["barcode_summary"] = barcode_summary[
            barcode_summary["lineage"].isin(include_descendants)
        ]
        result["min_consecutive"] = 1
        result["min_length"] = 1

    # ---------------------------------------------------------------------
    # XR: no diagnostic subs from second parent, only 2 consecutive barcodes
    #     relax region/subs lengths
    elif genome.lineage.recombinant in ["XR"]:
        genome.lineage.edge_case = True
        result["min_subs"] = 0
        result["min_consecutive"] = 2

    # ---------------------------------------------------------------------
    # XAD, XAE: second_parent only has 1 conflict sub
    #     force the first parent to be the minor parent (BA.1)
    elif genome.lineage.recombinant in ["XAD", "XAE"]:
        genome.lineage.edge_case = True
        include_tree = next(tree.find_clades("BA.1"))
        include_descendants = [c.name for c in include_tree.find_clades()]
        result["barcode_summary"] = barcode_summary[
            barcode_summary["lineage"].isin(include_descendants)
        ]
        result["min_consecutive"] = 5

    # ---------------------------------------------------------------------
    # XAJ: Parent 2 (BA.4) is tiny (120 nuc), force parent 1 to be BA.4
    #      and relax min_length.
    elif genome.lineage.recombinant in ["XAJ"]:
        genome.lineage.edge_case = True
        include_tree = next(tree.find_clades("BA.4"))
        include_descendants = [c.name for c in include_tree.find_clades()]
        result["barcode_summary"] = barcode_summary[
            barcode_summary["lineage"].isin(include_descendants)
        ]
        result["min_length"] = 4

    # ---------------------------------------------------------------------
    # XAS: The pango designation required deletions to resolve the first parent
    elif genome.lineage.recombinant in ["XAS"]:
        genome.lineage.edge_case = True
        include_tree = next(tree.find_clades("BA.4"))
        include_descendants = [c.name for c in include_tree.find_clades()]
        result["barcode_summary"] = barcode_summary[
            barcode_summary["lineage"].isin(include_descendants)
        ]

    # ---------------------------------------------------------------------
    # XAV: no diagnostic subs from second parent, only 2 consecutive barcodes
    #      BA.5.1.24 interferes
    elif genome.lineage.recombinant in ["XAV"]:
        genome.lineage.edge_case = True
        exclude_tree = next(tree.find_clades("BA.5.1.24"))
        exclude_descendants = [c.name for c in exclude_tree.find_clades()]
        result["min_subs"] = 0
        result["min_consecutive"] = 2
        result["barcode_summary"] = barcode_summary[
            ~barcode_summary["lineage"].isin(exclude_descendants)
        ]

    # ---------------------------------------------------------------------
    # XAY: Many breakpoints, some very small (200 nuc)
    elif genome.lineage.recombinant in ["XAY"]:
        genome.lineage.edge_case = True
        result["min_length"] = 200

    # ---------------------------------------------------------------------
    # XAZ: no diagnostic subs from BA.2, only 1 consecutive barcode from BA.2 parent
    #     force the minor parent (BA.2) to be the first parent
    #     this improves the search for the major parent (BA.5)
    elif genome.lineage.recombinant in ["XAZ"]:
        genome.lineage.edge_case = True
        include_tree = next(tree.find_clades("BA.2"))
        include_descendants = [c.name for c in include_tree.find_clades()]
        result["barcode_summary"] = barcode_summary[
            barcode_summary["lineage"].isin(include_descendants)
        ]
        result["min_subs"] = 0
        result["min_consecutive"] = 1
        result["min_length"] = 1

    # ---------------------------------------------------------------------
    # XBC: only 2 consecutive barcodes for first breakpoint
    elif genome.lineage.recombinant in ["XBC"]:
        result["min_consecutive"] = 2

    # ---------------------------------------------------------------------
    # XBK, XBQ: only 2 consecutive barcodes
    elif genome.lineage.recombinant in ["XBK", "XBQ"]:
        genome.lineage.edge_case = True
        include_tree = next(tree.find_clades("BA.2"))
        include_descendants = [c.name for c in include_tree.find_clades()]
        result["barcode_summary"] = barcode_summary[
            barcode_summary["lineage"].isin(include_descendants)
        ]
        result["min_consecutive"] = 2

    # ---------------------------------------------------------------------
    # XBZ: only 2 consecutive barcodes, extremely short parent 2 length
    elif genome.lineage.recombinant in ["XBZ"]:
        genome.lineage.edge_case = True
        result["min_consecutive"] = 2
        result["min_length"] = 300

    return result
