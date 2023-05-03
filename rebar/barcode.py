import yaml
import statistics

# import random
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
        diagnostic=None,
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
        self.diagnostic = []
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
            and recombinant_lineages
            and type(barcode_summary) == pd.core.frame.DataFrame
            and type(barcodes) == pd.core.frame.DataFrame
            and type(diagnostic) == pd.core.frame.DataFrame
            and type(lineage_to_clade) == pd.core.frame.DataFrame
        ):
            self.search(
                genome=genome,
                barcode_summary=barcode_summary,
                barcodes=barcodes,
                tree=tree,
                recombinant_lineages=recombinant_lineages,
                lineage_to_clade=lineage_to_clade,
                top_n=top_n,
                diagnostic=diagnostic,
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
            + "\ndiagnostic:      "
            + str(self.diagnostic)
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
        recombinant_lineages,
        tree,
        lineage_to_clade,
        diagnostic,
        max_lineages=10,
        top_n=1,
    ):
        # No barcode matches, stop the search
        if len(barcode_summary) == 0:
            return 0

        # ---------------------------------------------------------------------
        # Iterative Searches

        # print(barcode_summary)
        # Search Method #1 : Candidate Top Lineage Matches
        top_lineages = self.search_candidate_matches(barcode_summary, top_n)
        # print("search_1:", top_lineages)

        # Assume all lineages are outliers to begin with, we will remove
        # lineages from this list that survive all search methods.
        outlier_lineages = top_lineages

        # Search Method #2: Lineage-diagnostic mutations
        top_lineages = self.search_diagnostic_mutations(
            genome, top_lineages, diagnostic, tree
        )
        # print("search_2:", top_lineages)

        # If our top_lineages list is too long ( > subsample_threshold ), subsample
        top_lineages_subsample = top_lineages
        if len(top_lineages) > max_lineages:
            # Option 1: Take top subsample_threshold samples
            top_lineages_subsample = top_lineages[0:max_lineages]
            # Option 2: Randomly select samples, this is not working well
            # Fix the seed, so that results are the same every rerun
            # random.seed(123456)
            # top_lineages_subsample = random.choices(
            #     top_lineages, k=subsample_threshold
            # )

        # print("search_2:", top_lineages_subsample)
        # Search Method #3: Pairwise Distance
        #   If our top_lineages are a mix of recombinants and non-recombinants
        #   don't use this filter, because distances will be uninformative
        #   since recombinants have pseudo-tree placements under a fake "X" node
        top_lineages_rec = [l for l in top_lineages if l in recombinant_lineages]
        top_lineages_non_rec = [l for l in top_lineages if l not in top_lineages_rec]
        if len(top_lineages_rec) == 0 or len(top_lineages_non_rec) == 0:
            top_lineages = self.search_pairwise_distance(top_lineages_subsample, tree)
        else:
            top_lineages = top_lineages_subsample
        # print("search_3:", top_lineages)

        # Search Method #4: Maximum Parsimony
        top_lineages = self.search_maximum_parsimony(genome, top_lineages, barcodes)
        # print("search_4:", top_lineages)

        # ---------------------------------------------------------------------
        # Summarize Search

        outlier_lineages = [l for l in outlier_lineages if l not in top_lineages]
        lineage = tree.common_ancestor(top_lineages).name
        clade, clade_lineage = self.convert_lineage_to_clade(
            genome, lineage, lineage_to_clade
        )
        (
            lineage_barcode,
            support,
            conflict_alt,
            conflict_ref,
            missing,
        ) = self.summarise_top_lineages(genome, lineage, top_lineages, barcodes)

        # ---------------------------------------------------------------------
        # Update Attributes

        self.name = lineage
        self.definition = lineage
        self.clade = clade
        self.clade_lineage = clade_lineage
        self.top_lineages = top_lineages
        self.top_lineages_subsample = top_lineages_subsample
        self.outlier_lineages = outlier_lineages
        self.barcode = lineage_barcode
        self.support = support
        self.missing = missing
        self.conflict_ref = conflict_ref
        self.conflict_alt = conflict_alt

        return 0

    def search_candidate_matches(self, barcode_summary, top_n):
        # Identify the lineage(s) with the largest number of barcode matches
        # taking the top_n matches. Lineages such as XV require relaxation of
        # the top_n parameter (top_n=3) because XV is NOT the lineage with
        # the hightest number of matches.
        #   Example: XV
        #   top_n=3
        #   XJ [61, yes], XV [60, yes],  XY [60, yes], XAF [59, yes], XE [58, no]
        largest_totals = sorted(list(set(barcode_summary["total"])))
        largest_totals.reverse()
        max_barcodes = largest_totals[0:top_n]

        # Restrict to lineages with subs within top_n of largest total
        #   Lineage: B.1.634
        #   largest_total=27
        #   min_subs=25
        #   B.1.634 [27, yes], XB [15, no], XAY.3 [10, no], XAY.2.3 [10, no]
        max_total = largest_totals[0]
        min_subs = max_total - top_n - 1
        max_barcodes = [t for t in largest_totals[0:top_n] if t > min_subs]

        # Identify our top lineages based on the above criteria
        top_lineages = list(
            barcode_summary[barcode_summary["total"].isin(max_barcodes)]["lineage"]
        )

        return top_lineages

    def search_diagnostic_mutations(self, genome, top_lineages, diagnostic, tree):
        # ---------------------------------------------------------------------
        # Outlier Detection #1: Lineage Diagnostic Mutations
        #   If a sub in genome.substitutions is diagnostic for a particular
        #   lineage or its descendants, we will retain only those lineages.
        #
        #   Example: XBK
        #     XBK.1 [88], XBK [88], XBQ [87], CJ.1 [86], CJ [86], ...
        #     There are so many CJ.1 close matches, diagnostic mutations
        #     helps us resolve that XBK* is actually the best match.

        # if there's only one top lineage, just return that
        if len(top_lineages) <= 1:
            return top_lineages

        keep_lineages = []

        # Search for diagnostic mutations in the genome subs
        for s in genome.substitutions:

            s_row = diagnostic[diagnostic["mutation"] == str(s)]

            if len(s_row) == 0:
                continue

            s_lin = s_row["lineage"].values[0]
            s_include_desc = s_row["include_descendants"].values[0]

            s_top_lin = []
            # Complex, descendant match
            if s_include_desc:
                s_desc = [c.name for c in next(tree.find_clades(s_lin)).find_clades()]
                s_top_lin += [l for l in s_desc if l in top_lineages]

            # Simple, exact match
            elif s_lin in top_lineages:
                s_top_lin.append(s_lin)

            keep_lineages += s_top_lin

        # Remove duplicates
        keep_lineages = list(set(keep_lineages))

        # If we found keepers, return those
        if len(keep_lineages) > 0:
            return keep_lineages
        # otherwise, just return original top_lineages
        else:
            return top_lineages

    def search_pairwise_distance(self, top_lineages, tree):
        # ---------------------------------------------------------------------
        # Outlier Detection #2: Pairwise-Phylogenetic Distance
        #  If a lineage is too far away from the other candidate linages
        #  (ie. phylogenetic outlier), we will remove it. The value of this
        #  method is mainly when no diagnostic mutations are observed for
        #  detection method #1.
        #
        #  Example: XAT

        # if there's 2 or less top lineage, can't use this method
        if len(top_lineages) <= 2:
            return top_lineages

        else:
            distances_summary = {}

            # Calculate all pairwise distances between lineages
            # this is why we subsample, otherwise extraordinarily slow
            for l1 in top_lineages:
                distances = []
                for l2 in top_lineages:
                    if l1 == l2:
                        continue
                    distances.append(tree.distance(l1, l2))

                # Summarize the pairwise distances for this lineage by `mean`
                distances_summary[l1] = statistics.mean(distances)

            # The mode of all mean distances (confusing, I know) is how
            # we'll find the threshold for outliers
            distances_mode = statistics.mode(distances_summary.values())

            # keeper lineages are ones where there mean pairwise distance
            # was less than or equal to the mode (most frequently observed distance)
            keep_lineages = [
                l for l, d in distances_summary.items() if d <= distances_mode
            ]

        # If we found keepers, return those
        if len(keep_lineages) > 0:
            return keep_lineages
        # otherwise, just return original top_lineages
        else:
            return top_lineages

    def search_maximum_parsimony(self, genome, top_lineages, barcodes):
        # ---------------------------------------------------------------------
        # Outlier Detection #3: Maximum Parsimony (ie. minimum conflicts)
        # If lineages are tied for best match at this point, we will prefer the
        # lineage with the least sub conflicts with the genome.substitutions.
        # This is deliberately put after the pairwise-distance method, for reasons
        # I need to document.

        # Subsampling was already done in detection #2. But we need to exclude
        # any additional outliers from it. We don't directly modify the variable
        # top_lineages_subsample, because we want to display it in the end for debug.

        # if there's only one top lineage, there are no outliers
        if len(top_lineages) <= 1:
            return top_lineages

        else:
            parsimony_summary = {lin: 0 for lin in top_lineages}

            for lin in top_lineages:
                row = barcodes.query("lineage == @lin")
                subs = sorted(
                    [Substitution(s) for s in row.columns[1:] if list(row[s])[0] == 1]
                )

                # support: sub in genome also in candidate's barcode
                support = [s for s in genome.substitutions if s in subs]
                # conflict_alt: sub in genome that is not in candidate's barcode
                #               ie. unexpected ALT base.
                conflict_alt = [s for s in genome.substitutions if s not in subs]
                # conflict_ref: sub in candidate's barcode that is not in genome
                #               ie. unexpected REF base.
                conflict_ref = [
                    s
                    for s in subs
                    if s not in genome.substitutions
                    and s.coord not in genome.missing
                    and s.coord not in genome.deletions
                ]

                # our parsimony score is support - conflict
                parsimony_score = len(support) - (len(conflict_alt) + len(conflict_ref))
                parsimony_summary[lin] = parsimony_score

            # Identify maximum parsimony lineage
            max_parsimony_count = max(parsimony_summary.values())

            # Identify all the non-outliers in the subsampled lineages
            keep_lineages = [
                l for l, c in parsimony_summary.items() if c == max_parsimony_count
            ]

        # If we found keepers, return those
        if len(keep_lineages) > 0:
            return keep_lineages
        # otherwise, just return original top_lineages
        else:
            return top_lineages

    def convert_lineage_to_clade(self, genome, lineage, lineage_to_clade):
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
        return clade, clade_lineage

    def summarise_top_lineages(self, genome, lineage, top_lineages, barcodes):
        # There might be a case where a sub conflicts with the mrca of top_lineages
        # but is still found in all of the top_lineages. Don't consider these subs
        # to be conflicts.
        # Ex. XAJ. T15009C is a conflict for MRCA BA.2.12.1, but all top_lineages
        #          (BG*) have that sub.

        # Get the full barcode for the final lineage
        lineage_row = barcodes[barcodes["lineage"] == lineage]
        # No lineages may match if it was MRCA:
        lineage_barcode = []
        if len(lineage_row) > 0:
            lineage_barcode = [
                Substitution(s)
                for s in lineage_row.columns[1:]
                if list(lineage_row[s])[0] == 1
            ]

        # Identify subs for each top lineage
        top_lineages_subs = []
        for lin in top_lineages:
            row = barcodes[barcodes["lineage"] == lin]
            subs = sorted(
                [Substitution(s) for s in row.columns[1:] if list(row[s])[0] == 1]
            )
            top_lineages_subs += subs

        # Identify the subs that are shared among all
        top_lineages_subs_shared = sorted(
            [
                s
                for s in set(top_lineages_subs)
                if top_lineages_subs.count(s) == len(top_lineages)
            ]
        )

        # Get the barcode subs that were observed
        support = sorted([s for s in lineage_barcode if s in genome.substitutions])
        # Get the barcodes subs that were missing data
        missing = sorted(
            [
                s
                for s in lineage_barcode
                if s not in genome.substitutions and s.coord in genome.missing
            ]
        )
        # Get the barcode subs that were ref instead
        conflict_ref = sorted(
            [
                s
                for s in lineage_barcode
                if s not in genome.substitutions and s.coord not in genome.missing
            ]
        )
        # Get non-barcode subs, that were alt and unexpected
        # TBD: deletions excluded?
        conflict_alt = sorted(
            [
                s
                for s in genome.substitutions
                if s not in lineage_barcode and s not in top_lineages_subs_shared
            ]
        )
        conflict_subs = sorted(conflict_ref + conflict_alt)
        conflict_ref = [s for s in conflict_ref if s in conflict_subs]
        conflict_alt = [s for s in conflict_alt if s in conflict_subs]

        return lineage_barcode, support, conflict_alt, conflict_ref, missing

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
