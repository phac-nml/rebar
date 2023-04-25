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
        # This is the summary of lineages to search for parent_1
        "barcode_summary": barcode_summary,
    }

    # ---------------------------------------------------------------------
    # XB: top_lineages are tied exactly B.1.631 and B.1.634
    #     force the first parent to be B.1.631
    if genome.lineage.recombinant in ["XB"]:
        genome.lineage.edge_case = True
        include_tree = next(tree.find_clades("B.1.631"))
        include_descendants = [c.name for c in include_tree.find_clades()]
        result["barcode_summary"] = barcode_summary[
            barcode_summary["lineage"].isin(include_descendants)
        ]

    # ---------------------------------------------------------------------
    # XP: second parent comes from one barcode position: A29510C
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
    elif genome.lineage.recombinant in ["XR"]:
        genome.lineage.edge_case = True
        result["min_subs"] = 0
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

    # ---------------------------------------------------------------------
    # XAS: The pango designation required deletions to resolve the first parent
    #     force the first parent to be BA.2
    elif genome.lineage.recombinant in ["XAS"]:
        genome.lineage.edge_case = True
        include_tree = next(tree.find_clades("BA.2"))
        include_descendants = [c.name for c in include_tree.find_clades()]
        result["barcode_summary"] = barcode_summary[
            barcode_summary["lineage"].isin(include_descendants)
        ]

    # ---------------------------------------------------------------------
    # XAE: second_parent only has 1 conflict sub
    #     force the first parent to be the minor parent (BA.1)
    elif genome.lineage.recombinant in ["XAE"]:
        genome.lineage.edge_case = True
        include_tree = next(tree.find_clades("BA.1"))
        include_descendants = [c.name for c in include_tree.find_clades()]
        result["barcode_summary"] = barcode_summary[
            barcode_summary["lineage"].isin(include_descendants)
        ]
        result["min_consecutive"] = 5

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

    return result
