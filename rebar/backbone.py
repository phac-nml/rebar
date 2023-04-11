from .substitution import Substitution


class Backbone:
    def __init__(
        self,
        genome=None,
        barcode_summary=None,
        barcodes_df=None,
        tree=None,
        recombinant_lineages=None,
        include_recombinants=False,
    ):
        # Initialize attributes
        self.lineage = None
        self.max_lineages = None
        self.support_subs = None
        self.missing_subs = None
        self.conflict_subs_ref = None
        self.conflict_subs_alt = None

        # Run search
        if genome and barcode_summary and barcodes_df and tree:
            self.search(
                genome,
                barcode_summary,
                barcodes_df,
                tree,
                recombinant_lineages,
                include_recombinants,
            )

    def __repr__(self):
        text = (
            "lineage:           "
            + str(self.lineage)
            + "\nmax_lineages:      "
            + str(self.max_lineages)
            + "\nsupport_subs:      "
            + str(self.support_subs)
            + "\nmissing_subs:      "
            + str(self.missing_subs)
            + "\nconflict_subs_ref: "
            + str(self.conflict_subs_ref)
            + "\nconflict_subs_alt: "
            + str(self.conflict_subs_alt)
        )
        return text

    def search(
        self,
        genome,
        barcode_summary,
        barcodes_df,
        tree,
        recombinant_lineages=None,
        include_recombinants=False,
    ):

        # Optionally filter the barcodes to remove recombinant lineages
        if not include_recombinants and recombinant_lineages != False:
            barcode_summary = barcode_summary[
                ~barcode_summary["lineage"].isin(recombinant_lineages)
            ]

        # Identify the lineages with the largest number of barcode matches
        max_barcodes = barcode_summary["total"].max()
        max_lineages = list(
            barcode_summary[barcode_summary["total"] == max_barcodes]["lineage"]
        )

        # Get MRCA of all max_lineages
        lineage = tree.common_ancestor(max_lineages).name

        # Query the levels of support/conflict in this barcode
        lineage_row = barcodes_df.query("lineage == @lineage")
        lineage_subs = sorted(
            [
                Substitution(s)
                for s in lineage_row.columns[1:]
                if list(lineage_row[s])[0] == 1
            ]
        )

        # Get the barcode subs that were observed
        support_subs = sorted([s for s in lineage_subs if s in genome.subs])
        # Get the barcodes subs that were missing data
        missing_subs = sorted(
            [
                s
                for s in lineage_subs
                if s not in genome.subs and s.coord in genome.missing
            ]
        )
        # Get the barcode subs that were ref instead
        conflict_subs_ref = sorted(
            [
                s
                for s in lineage_subs
                if s not in genome.subs and s.coord not in genome.missing
            ]
        )
        # Get non-barcode subs, that were alt and unexpected
        # TBD: Unlabeled private sub coords excluded?
        # TBD: deletions excluded?
        # privates_unlabeled_coord = [p.coord for p in genome.privates_unlabeled]
        conflict_subs_alt = sorted(
            [
                s
                for s in genome.subs
                if s not in lineage_subs  # and s.coord not in privates_unlabeled_coord
            ]
        )
        conflict_subs = sorted(conflict_subs_ref + conflict_subs_alt)

        # Check if any of the max_lineages contain the conflict subs
        if len(max_lineages) > 1:
            for max_lin in max_lineages:
                # Filter the barcodes to just this lineage
                barcode_summary_filter = barcode_summary[
                    barcode_summary["lineage"] == max_lin
                ]
                max_lin_backbone = Backbone()
                max_lin_backbone.search(
                    genome=genome,
                    barcode_summary=barcode_summary_filter,
                    barcodes_df=barcodes_df,
                    tree=tree,
                    recombinant_lineages=recombinant_lineages,
                    include_recombinants=include_recombinants,
                )
                conflict_subs = [
                    s for s in conflict_subs if s not in max_lin_backbone.support_subs
                ]

        conflict_subs_ref = [s for s in conflict_subs_ref if s in conflict_subs]
        conflict_subs_alt = [s for s in conflict_subs_alt if s in conflict_subs]

        self.lineage = lineage
        self.max_lineages = max_lineages
        self.support_subs = support_subs
        self.missing_subs = missing_subs
        self.conflict_subs_ref = conflict_subs_ref
        self.conflict_subs_alt = conflict_subs_alt
