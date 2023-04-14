import pandas as pd
from .substitution import Substitution


class Barcode:
    def __init__(
        self,
        genome=None,
        barcode_summary=None,
        barcodes=None,
        tree=None,
        recombinant_lineages=None,
        recombinant_tree=None,
    ):
        # Initialize attributes
        self.lineage = None
        self.top_lineages = []
        self.recombinant = None
        self.recursive = None
        self.barcode = []
        self.support = []
        self.missing = []
        self.conflict_ref = []
        self.conflict_alt = []

        # Run search
        if (
            genome
            and tree
            and type(barcode_summary) == pd.core.frame.DataFrame
            and type(barcodes) == pd.core.frame.DataFrame
        ):
            self.search(
                genome,
                barcode_summary,
                barcodes,
                tree,
            )

        # Set recombinant status (self.recombinant and self.recursive)
        if recombinant_lineages and recombinant_tree:
            self.set_recombinant_status(
                recombinant_lineages=recombinant_lineages,
                recombinant_tree=recombinant_tree,
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

    def to_dict(self):
        barcode_dict = {
            "lineage": self.lineage,
            "top_lineages": ",".join(self.top_lineages),
            "barcode": ",".join([str(s) for s in self.barcode]),
            "support": ",".join([str(s) for s in self.support]),
            "missing": ",".join([str(s) for s in self.missing]),
            "conflict_ref": ",".join([str(s) for s in self.conflict_ref]),
            "conflict_alt": ",".join([str(s) for s in self.conflict_alt]),
        }
        return barcode_dict

    def search(self, genome, barcode_summary, barcodes, tree, max_lineages=10):

        # Identify the lineages with the largest number of barcode matches
        max_barcodes = barcode_summary["total"].max()
        top_lineages = list(
            barcode_summary[barcode_summary["total"] == max_barcodes]["lineage"]
        )
        top_lineages = top_lineages[0:max_lineages]

        # Get MRCA of all top_lineages
        lineage = tree.common_ancestor(top_lineages).name

        # Query the levels of support/conflict in this barcode
        lineage_row = barcodes.query("lineage == @lineage")
        print(lineage)
        print(lineage_row)
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
        conflict_ref = [s for s in conflict_ref if s in conflict_subs]
        conflict_alt = [s for s in conflict_alt if s in conflict_subs]

        self.lineage = lineage
        self.top_lineages = top_lineages
        self.barcode = lineage_subs
        self.support = support
        self.missing = missing
        self.conflict_ref = conflict_ref
        self.conflict_alt = conflict_alt

    def set_recombinant_status(self, recombinant_lineages, recombinant_tree):

        # Option 1: Top backbone lineage is s recombinant
        if self.lineage in recombinant_lineages:

            # Identify the generic recombinant type (XBB.1.5 = XBB)
            recombinant_path = recombinant_tree.get_path(self.lineage)
            # Move backwards up the path, until we find a parental lineage that
            # starts with "X", because it might be an alias ("EK")
            for c in recombinant_path[::-1]:
                if c.name.startswith("X"):
                    self.recombinant = c.name.split(".")[0]
                    break

            # Check if this is a recursive recombinant
            # Note: In the get_path method, the root is excluded
            node_path = recombinant_tree.get_path(self.recombinant)
            # So if node_path just has more than one clade (ex. [XBB, XBL]),
            # it's recursive
            if len(node_path) > 1:
                self.recursive = True

        # Option 2: Perfect match to non-recombinant
        elif len(self.conflict_ref) == 0:
            self.recombinant = False

        return 0
