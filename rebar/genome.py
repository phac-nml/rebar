import copy
import os
import subprocess
import pandas as pd
import numpy as np
from .utils import NO_DATA_CHAR
from .substitution import Substitution


class Genome:
    def __init__(self, id, genome_len, nextclade_row, barcode):
        self.id = id

        self.genome_len = genome_len
        self.subs = self.parse_subs(nextclade_row, col="substitutions")
        self.privates_reversion = self.parse_subs(
            nextclade_row, col="privateNucMutations.reversionSubstitutions"
        )
        self.privates_labeled = self.parse_subs(
            nextclade_row, col="privateNucMutations.labeledSubstitutions"
        )
        self.privates_unlabeled = self.parse_subs(
            nextclade_row, col="privateNucMutations.unlabeledSubstitutions"
        )
        self.missing = self.parse_intervals(nextclade_row, col="missing")
        self.deletions = self.parse_intervals(nextclade_row, col="deletions")

        # Check which substitutions are "barcodes" (lineage-defining in combination)
        self.barcode_summary = self.identify_barcode(barcode)

        # Attributes that will be set later
        self.lineage = None
        self.is_recombinant = None
        self.backbone = None

    def parse_intervals(self, nextclade_row, col):
        col_str = nextclade_row[col].values[0].split(",")
        col_coords = []

        for coords in col_str:
            if coords == NO_DATA_CHAR:
                continue
            coords_split = [int(m) for m in coords.split("-")]
            if len(coords_split) == 1:
                coords = [coords_split[0]]
            else:
                coords = list(range(coords_split[0], coords_split[1] + 1))
            col_coords += coords

        return col_coords

    def identify_barcode(self, barcode):

        barcode_subs = [str(s) for s in self.subs if str(s) in barcode.columns]
        # Count up barcode mutations by lineage
        df = copy.copy(barcode[["lineage"] + barcode_subs])
        df["total"] = df[barcode_subs].sum(axis=1)
        total_df = (
            df[["lineage", "total"]]
            .query("total > 0")
            .sort_values(by="total", ascending=False)
        )
        return total_df

    def parse_subs(self, nextclade_row, col):
        """Parse substitutions from Nextclade records."""
        subs_str = nextclade_row[col].values[0].split(",")

        # labeled mutation is special format
        if col == "privateNucMutations.labeledSubstitutions":
            subs_str = [s.split("|")[0] for s in subs_str if "|" in s]

        subs = sorted([Substitution(s) for s in set(subs_str) if s != NO_DATA_CHAR])
        return subs

    def identify_backbone(
        self,
        barcode_summary,
        barcodes_df,
        tree,
        recombinant_lineages=None,
        include_recombinants=False,
    ):

        backbone = {}

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
        support_subs = sorted([s for s in lineage_subs if s in self.subs])
        # Get the barcode subs that were not observed
        missing_subs = sorted([s for s in lineage_subs if s not in self.subs])
        # Get the sample subs that are not in the barcode
        # Unlabeled private sub coords should not be included in conflict
        privates_unlabeled_coord = [p.coord for p in self.privates_unlabeled]
        conflict_subs = sorted(
            [
                s
                for s in self.subs
                if s not in lineage_subs
                if s.coord not in privates_unlabeled_coord
            ]
        )

        # Check if any of the max_lineages contain the conflict subs
        if len(max_lineages) > 1:
            for max_lin in max_lineages:
                # Filter the barcodes to just this lineage
                barcode_summary_filter = barcode_summary[
                    barcode_summary["lineage"] == max_lin
                ]
                max_lin_backbone = self.identify_backbone(
                    barcode_summary=barcode_summary_filter,
                    barcodes_df=barcodes_df,
                    tree=tree,
                    recombinant_lineages=recombinant_lineages,
                    include_recombinants=include_recombinants,
                )
                conflict_subs = [
                    s
                    for s in conflict_subs
                    if s not in max_lin_backbone["support_subs"]
                ]

        backbone["lineage"] = lineage
        backbone["max_lineages"] = max_lineages
        backbone["support_subs"] = support_subs
        backbone["missing_subs"] = missing_subs
        backbone["conflict_subs"] = conflict_subs

        return backbone

    def identify_breakpoints(self, alt_backbone, region_min_subs=1, region_min_len=1):

        result = {}

        parent_1 = self.backbone["lineage"]
        parent_2 = alt_backbone["lineage"]

        # Identify the substitutions that are uniq to each parent
        parent_1_subs = [
            s
            for s in self.backbone["support_subs"]
            if s not in alt_backbone["support_subs"]
        ]
        parent_2_subs = [
            s
            for s in alt_backbone["support_subs"]
            if s not in self.backbone["support_subs"]
        ]
        # Identify the substitutions that are shared by both parents
        shared_subs = [
            s
            for s in self.backbone["support_subs"]
            if s in alt_backbone["support_subs"]
        ]

        # Organize into a dataframe for easy parsing
        subs_df = pd.DataFrame(
            {
                "substitution": shared_subs + parent_1_subs + parent_2_subs,
                "parent": ["shared"] * len(shared_subs)
                + [parent_1] * len(parent_1_subs)
                + [parent_2] * len(parent_2_subs),
            }
        ).sort_values(by="substitution")
        subs_df["strain"] = [self.id] * len(subs_df)

        # Search for genomic blocks for each parent
        # Just look at the subs that are uniq to one parent
        subs_uniq_df = subs_df.query("parent != 'shared'")

        regions = {}
        p_prev = None
        start = 0
        end = 0

        for rec in subs_uniq_df.iterrows():
            p_curr = rec[1]["parent"]
            sub = rec[1]["substitution"]
            # First region
            if not p_prev:
                start = sub.coord
                end = sub.coord
                regions[start] = {
                    "start": start,
                    "end": end,
                    "parent": p_curr,
                    "subs": [sub],
                }
            # Same parent, region continues
            elif p_curr == p_prev:
                regions[start]["end"] = sub.coord
                regions[start]["subs"].append(sub)
            # Parent change, start of new region
            elif p_curr != p_prev:
                start = sub.coord
                end = sub.coord
                regions[start] = {
                    "start": start,
                    "end": end,
                    "parent": p_curr,
                    "subs": [sub],
                }

            end = sub.coord
            p_prev = p_curr

        # Filter recombination regions
        #   - At least X consecutive subs per region
        #   - At least X bases per region
        regions_filter = {}
        for start in regions:
            num_subs = len(regions[start]["subs"])
            region_len = regions[start]["end"] - start
            if num_subs < region_min_subs:
                continue
            if region_len < region_min_len:
                continue
            regions_filter[start] = regions[start]

        # If we're left with one filtered region, the alt parent fails
        if len(regions_filter) < 2:
            return 0

        # Identify breakpoints
        breakpoints = []
        prev_start_coord = None
        prev_end_coord = None

        for start_coord in regions_filter:

            end_coord = regions_filter[start_coord]["end"]

            # Skip the first record for breakpoints
            if prev_start_coord:
                breakpoint_start = prev_end_coord + 1
                breakpoint_end = start_coord - 1
                breakpoint = "{}:{}".format(breakpoint_start, breakpoint_end)
                breakpoints.append(breakpoint)

            prev_start_coord = start_coord
            prev_end_coord = end_coord

        # Check if the missing and conflict subs from parent 1 (self.backbone) can
        # be explained by recombination from parent 2 (alt_backbone)
        # parent_1_missing = self.cleanup_conflicts(
        #     self.backbone, regions_filter, "missing_subs"
        # )
        # parent_2_missing = self.cleanup_conflicts(
        #     alt_backbone, regions_filter, "missing_subs"
        # )

        parent_1_conflict = self.cleanup_conflicts(
            self.backbone, regions_filter, "conflict_subs"
        )
        parent_2_conflict = self.cleanup_conflicts(
            alt_backbone, regions_filter, "conflict_subs"
        )

        # Missing subs, hmm

        # Any remaining conflicts will be private mutations
        private_subs = parent_1_conflict + parent_2_conflict
        subs_df = pd.concat(
            [
                subs_df,
                pd.DataFrame(
                    {
                        "substitution": private_subs,
                        "parent": ["private"] * len(private_subs),
                    }
                ),
            ]
        ).sort_values(by="substitution")

        result["subs_df"] = subs_df
        result["regions"] = regions_filter
        result["breakpoints"] = breakpoints

        return result

    def cleanup_conflicts(self, backbone, regions_filter, col):
        subs = copy.copy(backbone[col])
        lineage = backbone["lineage"]
        for s in backbone[col]:
            region_found = False
            for r in regions_filter.values():
                start = r["start"]
                end = r["end"]
                parent = r["parent"]

                # Does it fall within another parent's region?
                if s.coord >= start and s.coord <= end:
                    region_found = True
                    if lineage != parent:
                        subs.remove(s)

            # Does it fall within no known parental region?
            if not region_found:
                subs.remove(s)

        return subs

    def write_snipit(self, outdir, exclude_shared=False, ext="pdf"):
        if exclude_shared:
            self.recombination_subs = self.recombination_subs.query(
                "parent != 'shared'"
            )

        parents = list(np.unique(self.recombination_subs["parent"]))
        if "shared" in parents:
            parents.remove("shared")
        if "private" in parents:
            parents.remove("private")

        # This will be our crude object to hold fasta lines
        output_lines = []
        blank_seq = ["N"] * self.genome_len

        for strain in ["Reference"] + parents + [self.id]:

            header = ">" + strain
            # Start off with all Ns, we'll patch in subs as we go
            seq = blank_seq

            for rec in self.recombination_subs.iterrows():
                sub = rec[1]["substitution"]
                sub_parent = rec[1]["parent"]

                # Decide whether we want to write the ref or alt
                # mutation, depending on the currently sample
                # and the mutation's parental origin

                # Reference strain, always use ref
                if strain == "Reference":
                    seq[sub.coord - 1] = sub.ref

                # genomic sample, always use alt
                elif strain == self.id:
                    seq[sub.coord - 1] = sub.alt

                # Shared by all parents, always use alt
                elif sub_parent == "shared":
                    seq[sub.coord - 1] = sub.alt

                # Sub originates from the current strain, use alt
                elif sub_parent == strain:
                    seq[sub.coord - 1] = sub.alt

                # Private mutation, and not a genomic sample
                elif strain != self.id and sub_parent == "private":
                    seq[sub.coord - 1] = sub.ref

                # otherwise, use ref
                else:
                    seq[sub.coord - 1] = sub.ref

            output_lines.append(header)
            output_lines.append("".join(seq))

        # Write the alignment
        fasta_name = "snipit_{}_{}.fasta".format(parents[0], parents[1])
        fasta_path = os.path.join(outdir, fasta_name)
        with open(fasta_path, "w") as outfile:
            outfile.write("\n".join(output_lines) + "\n")

        # Write the figure
        fig_name = "snipit_{}_{}".format(parents[0], parents[1])
        fig_path = os.path.join(outdir, fig_name)
        cmd_str = (
            "snipit {fasta_path}"
            " --output-file {fig_path}"
            " --recombi-mode"
            " --recombi-references {parent_1},{parent_2}"
            " --flip-vertical"
            " --format {ext}"
        ).format(
            fasta_path=fasta_path,
            fig_path=fig_path,
            parent_1=parents[0],
            parent_2=parents[1],
            ext=ext,
        )
        subprocess.run(cmd_str, shell=True)

    def __repr__(self):
        return self.id
