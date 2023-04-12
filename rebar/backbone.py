# import pandas as pd
from .substitution import Substitution
from . import RebarError


class Backbone:
    def __init__(
        self,
        genome=None,
        barcode_summary=None,
        barcodes=None,
        tree=None,
    ):
        # Initialize attributes
        self.lineage = None
        self.top_lineages = None
        self.barcode = None
        self.support = None
        self.missing = None
        self.conflict_ref = None
        self.conflict_alt = None

        # Run search
        if genome and barcode_summary and barcodes and tree:
            self.search(
                genome,
                barcode_summary,
                barcodes,
                tree,
            )

    def __repr__(self):
        text = (
            "lineage:      "
            + str(self.lineage)
            + "\ntop_lineages: "
            + str(self.top_lineages)
            + "\nbarcode:      "
            + str(self.barcode)
            + "\nsupport:      "
            + str(self.support)
            + "\nmissing:      "
            + str(self.missing)
            + "\nconflict_ref: "
            + str(self.conflict_ref)
            + "\nconflict_alt: "
            + str(self.conflict_alt)
            + "\n"
        )
        return text

    def search(self, genome, barcode_summary, barcodes, tree):

        # Identify the lineages with the largest number of barcode matches
        max_barcodes = barcode_summary["total"].max()
        top_lineages = list(
            barcode_summary[barcode_summary["total"] == max_barcodes]["lineage"]
        )

        # Get MRCA of all top_lineages
        lineage = tree.common_ancestor(top_lineages).name

        # Query the levels of support/conflict in this barcode
        lineage_row = barcodes.query("lineage == @lineage")
        lineage_subs = sorted(
            [
                Substitution(s)
                for s in lineage_row.columns[1:]
                if list(lineage_row[s])[0] == 1
            ]
        )

        # Get the barcode subs that were observed
        support = sorted([s for s in lineage_subs if s in genome.substitutions])
        # Get the barcodes subs that were missing data
        missing = sorted(
            [
                s
                for s in lineage_subs
                if s not in genome.substitutions and s.coord in genome.missing
            ]
        )
        # Get the barcode subs that were ref instead
        conflict_ref = sorted(
            [
                s
                for s in lineage_subs
                if s not in genome.substitutions and s.coord not in genome.missing
            ]
        )
        # Get non-barcode subs, that were alt and unexpected
        # TBD: Unlabeled private sub coords excluded?
        # TBD: deletions excluded?
        # privates_unlabeled_coord = [p.coord for p in genome.privates_unlabeled]
        conflict_alt = sorted(
            [
                s
                for s in genome.substitutions
                if s not in lineage_subs  # and s.coord not in privates_unlabeled_coord
            ]
        )
        conflict_subs = sorted(conflict_ref + conflict_alt)

        # Check if any of the top_lineages contain the conflict subs
        if len(top_lineages) > 1:
            for max_lin in top_lineages:
                # Filter the barcodes to just this lineage
                barcode_summary_filter = barcode_summary[
                    barcode_summary["lineage"] == max_lin
                ]
                max_lin_backbone = Backbone()
                max_lin_backbone.search(
                    genome=genome,
                    barcode_summary=barcode_summary_filter,
                    barcodes=barcodes,
                    tree=tree,
                )
                conflict_subs = [
                    s for s in conflict_subs if s not in max_lin_backbone.support
                ]

        conflict_ref = [s for s in conflict_ref if s in conflict_subs]
        conflict_alt = [s for s in conflict_alt if s in conflict_subs]

        self.lineage = lineage
        self.top_lineages = top_lineages
        self.barcode = lineage_subs
        self.support = support
        self.missing = missing
        self.conflict_ref = conflict_ref
        self.conflict_alt = conflict_alt


def backbone_search_lineage_mp(
    iterator, barcodes, tree, recombinant_lineages, recombinant_tree
):
    """
    Run a backbone search using multiprocess.
    This is for lineage assignment, not recombination detection.
    Used to control the named parameters that are passed.
    """
    genome = iterator

    genome.backbone.search(
        genome=genome,
        barcode_summary=genome.barcode_summary,
        barcodes=barcodes,
        tree=tree,
    )

    genome.lineage = genome.backbone.lineage

    # Option 1: Top backbone lineage is recombinant
    #           Investigate further in next phase
    if genome.backbone.lineage in recombinant_lineages:

        # Identify the generic recombinant type
        recombinant_path = recombinant_tree.get_path(genome.backbone.lineage)
        # Move backwards up the path, until we find a parental lineage that
        # starts with "X", because it might be an alias ("EK")

        for c in recombinant_path[::-1]:
            if c.name.startswith("X"):
                genome.recombinant = c.name.split(".")[0]
                break

        # Check if this is a recursive recombinant
        # In the get_path method, the root is excluded
        node_path = recombinant_tree.get_path(genome.recombinant)
        # So if node_path just has more than one clade, it's recursive
        if len(node_path) > 1:
            genome.recursive = True

    # Option 2: Perfect match to non-recombinant
    #           Stop further investigation
    elif len(genome.backbone.conflict_ref) == 0:
        genome.recombinant = False

    # else. Has conflicts with non-recombinant lineage
    #           Investigate further in next phase
    return genome


def backbone_search_recombination_mp(
    iterator, parent, barcodes, tree, recombinant_lineages, recombinant_tree
):
    """
    Run a backbone search using multiprocess.
    This is for the first recombination parent detection.
    Used to control the named parameters that are passed.
    """
    genome = iterator

    # Skip clear non-recombinants
    if genome.recombinant == False:
        return genome

    # Save a copy of the barcode summary, before we modify it
    barcode_summary = genome.barcode_summary

    # Option 1. Definitely not a recursive recombinant.
    #           Exclude all recombinant lineages from new search.
    #           Ex. XBB.1.5 is not a recursive recombinant (BA.2.10* and BA.2.75*)
    #           If we remove all recombinant lineages from it's barcode summary
    #           the top lineage will become BJ.1.1 (BA.2.10*)
    if not genome.recursive:
        barcode_summary = barcode_summary[
            ~barcode_summary["lineage"].isin(recombinant_lineages)
        ]
    # Option 2. Potentially recursive recombinant
    #           Exclude only original backbone lineages from new search.
    #           Ex. XBL is a recursive recombinant (XBB.1* and BA.2.75*)
    else:
        barcode_summary = barcode_summary[
            ~barcode_summary["lineage"].isin(genome.backbone.top_lineages)
        ]

    if parent == "parent_1":
        genome.parent_1.search(
            genome=genome, barcode_summary=barcode_summary, barcodes=barcodes, tree=tree
        )
    elif parent == "parent_2":
        # Exclude the backbone and parent_1 lineages
        exclude_lineages = genome.backbone.top_lineages + genome.parent_1.top_lineages
        barcode_summary = barcode_summary[
            ~barcode_summary["lineage"].isin(exclude_lineages)
        ]
        print("-" * 80)
        print(genome)
        alt_backbone = Backbone()
        alt_backbone.search(
            genome=genome, barcode_summary=barcode_summary, barcodes=barcodes, tree=tree
        )
        print(alt_backbone)
        # Check for breakpoints
        result = genome.breakpoints(alt_backbone)
        print(result)
    else:
        raise RebarError("Unknown parent supplied to backbone_search_recombination_mp.")

    return genome
