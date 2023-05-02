import yaml
import statistics
import random
from datetime import datetime

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
        lineage_to_clade=None,
        top_n=1,
    ):
        # Initialize attributes
        self.name = None
        self.clade = None
        self.clade_lineage = None
        self.top_lineages = []
        self.top_lineages_subsample = []
        self.outlier_lineages = []
        self.recombinant = None
        self.recursive = None
        self.edge_case = False
        self.barcode = []
        self.support = []
        self.missing = []
        self.conflict_ref = []
        self.conflict_alt = []
        self.definition = None
        self.definition_aa = None

        # Run search
        if (
            genome
            and tree
            and type(barcode_summary) == pd.core.frame.DataFrame
            and type(barcodes) == pd.core.frame.DataFrame
            and type(lineage_to_clade) == pd.core.frame.DataFrame
        ):
            self.search(
                genome,
                barcode_summary,
                barcodes,
                tree,
                lineage_to_clade,
                top_n=top_n,
            )

        # Set recombinant status (self.recombinant and self.recursive)
        if genome and recombinant_lineages and recombinant_tree:
            self.set_recombinant_status(
                genome=genome,
                recombinant_lineages=recombinant_lineages,
                recombinant_tree=recombinant_tree,
            )

    def __repr__(self):
        text = (
            "lineage:      "
            + str(self.name)
            + "definition:      "
            + str(self.definition)
            + "clade:      "
            + str(self.clade)
            + "clade_lineage:      "
            + str(self.clade_lineage)
            + "\ntop_lineages: "
            + str(self.top_lineages)
            + "\ntop_lineages_subsample: "
            + str(self.top_lineages_subsample)
            + "\noutlier_lineages: "
            + str(self.outlier_lineages)
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
            + "\nrecombinant: "
            + str(self.recombinant)
            + "\nrecursive: "
            + str(self.recursive)
            + "\nedge_case: "
            + str(self.edge_case)
            + "\n"
        )
        return text

    def to_dict(self):
        barcode_dict = {
            "lineage": self.name,
            "definition": self.definition,
            "clade": self.clade,
            "clade_lineage": self.clade_lineage,
            "top_lineages": ",".join(self.top_lineages),
            "top_lineages_subsample": ",".join(self.top_lineages_subsample),
            "outlier_lineages": ",".join(self.outlier_lineages),
            "barcode": ",".join([str(s) for s in self.barcode]),
            "support": ",".join([str(s) for s in self.support]),
            "missing": ",".join([str(s) for s in self.missing]),
            "conflict_ref": ",".join([str(s) for s in self.conflict_ref]),
            "conflict_alt": ",".join([str(s) for s in self.conflict_alt]),
            "recombinant": str(self.recombinant),
            "recursive": str(self.recursive),
            "edge_case": str(self.edge_case),
        }
        return barcode_dict

    def to_yaml(self, indent=2):
        """
        Convert Barcode object to yaml.

        Returns
        -------
        genome_yaml : yaml
            YAML representation of Barcode.
        """

        barcode_yaml = (
            yaml.dump(self.to_dict(), sort_keys=False, indent=indent)
            .replace("null", "")
            .replace("''", "")
            + "\n"
        )
        return barcode_yaml

    def search(
        self,
        genome,
        barcode_summary,
        barcodes,
        tree,
        lineage_to_clade,
        subsample_threshold=10,
        top_n=1,
    ):

        # No barcode matches
        if len(barcode_summary) == 0:
            return 0

        # Identify the lineages with the largest number of barcode matches
        largest_totals = sorted(list(set(barcode_summary["total"])))
        largest_totals.reverse()
        max_barcodes = largest_totals[0:top_n]

        # Subs must be within top_n of largest total
        max_total = largest_totals[0]
        max_barcodes = [t for t in largest_totals[0:top_n] if t > (max_total - top_n)]

        # max_barcodes = barcode_summary["total"].max()
        top_lineages = list(
            barcode_summary[barcode_summary["total"].isin(max_barcodes)]["lineage"]
        )
        top_lineages_subsample = top_lineages

        # Get MRCA of all top_lineages
        lineage = tree.common_ancestor(top_lineages).name
        distances = {}
        # Check for top lineage outliers
        for l in top_lineages:
            distances[l] = tree.distance(lineage, l)

        outlier_lineages = []

        # ---------------------------------------------------------------------
        # Outlier Detection

        # if there's only one or two top lineages, there's no outliers
        if len(top_lineages) < 2:
            outlier_lineages = []

        else:

            # Fix the seed, so that results are the same every rerun
            random.seed(123456)
            # If there are a large number of top_lineages, subsample down for speed
            if len(top_lineages) > subsample_threshold:
                top_lineages_subsample = random.choices(
                    top_lineages, k=subsample_threshold
                )

            # 1ST OUTLIER STAGE: phylogenetic distance to mrca
            #   Ex. needed for XAT
            distances_summary = {}
            for l1 in top_lineages_subsample:
                distances = []
                for l2 in top_lineages_subsample:
                    if l1 == l2:
                        continue
                    distances.append(tree.distance(l1, l2))
                distances_summary[l1] = statistics.mean(distances)

            distances_mode = statistics.mode(distances_summary.values())

            # Identify all the non-outliers in the subsampled lineages
            keep_lineages = [
                l for l, d in distances_summary.items() if d <= distances_mode
            ]
            # Grab the MRCA+descendants of the non-outliers subsample
            lineage_tree = tree.common_ancestor(keep_lineages)
            lineage = lineage_tree.name
            lineage_descendants = [c.name for c in lineage_tree.find_clades()]
            # Filter out top_lineages if outside the MRCA clade
            outlier_lineages += [
                l for l in top_lineages if l not in lineage_descendants
            ]

            print("outlier_lineages_1:", outlier_lineages)

            # 2ND OUTLIER STAGE: minimum conflicts, aka parsimony
            top_lineages_subsample_filtered = [
                l for l in top_lineages_subsample if l not in outlier_lineages
            ]
            conflict_count = {lin: 0 for lin in top_lineages_subsample_filtered}
            for lin in top_lineages_subsample_filtered:
                row = barcodes.query("lineage == @lin")
                subs = sorted(
                    [Substitution(s) for s in row.columns[1:] if list(row[s])[0] == 1]
                )
                conflict_alt = [s for s in genome.substitutions if s not in subs]
                conflict_ref = [s for s in subs if s not in genome.substitutions]

                conflict_total = len(conflict_alt) + len(conflict_ref)
                conflict_count[lin] = conflict_total

            min_conflict_count = min(conflict_count.values())

            # Identify all the non-outliers in the subsampled lineages
            keep_lineages = [
                lin
                for lin, count in conflict_count.items()
                if count == min_conflict_count
            ]

            # Grab the MRCA+descendants of the non-outliers subsample
            lineage_tree = tree.common_ancestor(keep_lineages)
            lineage = lineage_tree.name
            lineage_descendants = [c.name for c in lineage_tree.find_clades()]
            # Filter out top_lineages if outside the MRCA clade
            outlier_lineages += [
                l for l in top_lineages if l not in lineage_descendants
            ]

            print("outlier_lineages_2:", outlier_lineages)

        # ---------------------------------------------------------------------
        # Lineage to Clade

        # Get clade of lineage
        if lineage in list(lineage_to_clade["lineage"]):
            clade = lineage_to_clade[lineage_to_clade["lineage"] == lineage][
                "nextstrainClade"
            ].values[0]
            clade_lineage = lineage_to_clade[lineage_to_clade["lineage"] == lineage][
                "nextstrainClade_lineage"
            ].values[0]
        elif lineage in ["MRCA", "X"]:
            clade = lineage
            clade_lineage = lineage
        else:
            clade = None
            clade_lineage = None
            if genome.debug:
                genome.logger.info(
                    str(datetime.now())
                    + "\t\t\tWARNING: unknown clade for lineage "
                    + str(lineage)
                )

        # ---------------------------------------------------------------------
        # Substitution Conflicts

        # There might be a case where a sub conflicts with the mrca of top_lineages
        # but is still found in all of the top_lineages. Don't consider these subs
        # to be conflicts.
        # Ex. XAJ. T15009C is a conflict for MRCA BA.2.12.1, but all top_lineages
        #          (BG*) have that sub.
        # Identify the subs that all the top lineages shared
        top_lineages_subs = []
        for lin in top_lineages_subsample:
            if lin in outlier_lineages:
                continue
            row = barcodes.query("lineage == @lin")
            subs = sorted(
                [Substitution(s) for s in row.columns[1:] if list(row[s])[0] == 1]
            )
            top_lineages_subs += subs
        top_lineages_subs_shared = sorted(
            [
                s
                for s in set(top_lineages_subs)
                if top_lineages_subs.count(s) == len(top_lineages)
            ]
        )

        # Query the levels of support/conflict in this barcode
        lineage_row = barcodes.query("lineage == @lineage")

        # No lineages matching, possibly because it's MRCA:
        if len(lineage_row) == 0:
            lineage_subs = []
        else:
            lineage_subs = [
                Substitution(s)
                for s in lineage_row.columns[1:]
                if list(lineage_row[s])[0] == 1
            ]
            # Add in the top_lineages shared subs
            # This interferes with plotting of XBJ. Where is it a benefit?
            # lineage_subs = sorted(list(set(lineage_subs + top_lineages_subs_shared)))

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
        # TBD: deletions excluded?
        conflict_alt = sorted(
            [
                s
                for s in genome.substitutions
                if s not in lineage_subs and s not in top_lineages_subs_shared
            ]
        )
        conflict_subs = sorted(conflict_ref + conflict_alt)
        conflict_ref = [s for s in conflict_ref if s in conflict_subs]
        conflict_alt = [s for s in conflict_alt if s in conflict_subs]

        # ---------------------------------------------------------------------
        # Update Attributes

        self.name = lineage
        self.definition = lineage
        self.clade = clade
        self.clade_lineage = clade_lineage
        self.top_lineages = top_lineages
        self.top_lineages_subsample = top_lineages_subsample
        self.outlier_lineages = outlier_lineages
        self.barcode = lineage_subs
        self.support = support
        self.missing = missing
        self.conflict_ref = conflict_ref
        self.conflict_alt = conflict_alt

        return 0

    def set_recombinant_status(self, genome, recombinant_lineages, recombinant_tree):

        if self.name == "X":
            self.recombinant = "X"
            self.recursive = False

        # Option 1: Top backbone lineage is s recombinant
        elif self.name in recombinant_lineages:

            # Identify the generic recombinant type (XBB.1.5 = XBB)
            recombinant_path = recombinant_tree.get_path(self.name)
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

    def set_definition(self):
        self.definition = self.name
        if len(self.conflict_alt) > 0:
            self.definition += "+" + ",".join([str(s) for s in self.conflict_alt])
