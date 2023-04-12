import copy

# import os
# import subprocess
import pandas as pd

# import numpy as np
from Bio.SeqRecord import SeqRecord
from .utils import NO_DATA_CHAR
from .substitution import Substitution
from .backbone import Backbone
from . import RebarError


class Genome:
    def __init__(self, record=None, reference=None, subs_row=None, barcodes=None):

        # Generic genomic features
        self.id = None
        self.seq = None
        self.substitutions = []
        self.deletions = []
        self.missing = []
        self.genome_length = None

        # Lineage features
        self.barcode_summary = None
        self.lineage = None

        # Recombination features
        self.recombinant = None
        self.recursive = None
        self.backbone = Backbone()
        self.parent_1 = Backbone()
        self.parent_2 = Backbone()

        # Entry point #1, from fasta alignment
        if record:
            self.id = record.id
            self.seq = record.seq

        if reference and self.seq:
            self.parse_sequence(reference)

        # Entry point #2, from subs dataframe
        if type(subs_row) == pd.core.series.Series:
            self.id = subs_row["strain"]
            self.genome_length = subs_row["genome_length"]
            self.substitutions = self.parse_subs(subs_row, "substitutions")
            self.deletions = self.ranges_to_coords(values=subs_row["deletions"])
            self.missing = self.ranges_to_coords(values=subs_row["missing"])

        # Check which substitutions are "barcodes" (lineage-defining in combination)
        if type(barcodes) == pd.core.frame.DataFrame:
            self.barcode_summary = self.parse_barcodes(barcodes)

    def __repr__(self):
        text = (
            self.id
            + ":"
            + "\n  substitutions: "
            + str(self.substitutions)
            + "\n  deletions: "
            + str(self.deletions)
            + "\n  missing: "
            + str(self.missing)
            + "\n  lineage: "
            + str(self.lineage)
            + "\n  recombinant: "
            + str(self.recombinant)
            + "\n  recursive: "
            + str(self.recursive)
            + "\n  backbone:    "
            + "\n    "
            + str(self.backbone).replace("\n", "\n    ")
            + "\n  parent_1:    "
            + "\n    "
            + str(self.parent_1).replace("\n", "\n    ")
            + "\n  parent_2:    "
            + "\n    "
            + str(self.parent_2).replace("\n", "\n    ")
        )
        return text

    def parse_sequence(self, reference):
        """Parse fasta alignment into TSV records."""
        # Dashes (-) at the 5' and 3' end should be treated as
        # missing data, not true deletions
        terminal_5p = True
        terminal_3p = False
        terminal_char = ["N", "-", "N-"]

        for i, bases in enumerate(zip(reference.seq, self.seq)):
            r = bases[0]
            s = bases[1]
            # Genomic coordinates are 1 based
            coord = i + 1

            # Collapse all downstream and upstream bases
            upstream_bases = "".join(set(self.seq[:i]))
            downstream_bases = "".join(set(self.seq[i + 1 :]))

            if (
                len(upstream_bases) > 0
                and terminal_5p
                and upstream_bases not in terminal_char
            ):
                terminal_5p = False
            if (
                len(downstream_bases) > 0
                and not terminal_3p
                and downstream_bases in terminal_char
            ):
                terminal_3p = True

            # Missing Data
            if s == "N":
                self.missing.append(coord)
            elif (s == "-" and terminal_5p) or (s == "-" and terminal_3p):
                self.missing.append(coord)

            # Deletions
            elif s == "-":
                self.deletions.append(coord)

            # Substitutions
            elif r == "N":
                continue

            elif r != s:
                sub = "{ref}{coord}{alt}".format(ref=r, coord=coord, alt=s)
                self.substitutions.append(Substitution(sub))

        self.genome_length = coord
        return 0

    def parse_subs(self, subs_row, col):
        """Parse substitutions from dataframe row."""
        subs_str = subs_row[col].split(",")
        subs = sorted([Substitution(s) for s in set(subs_str) if s != NO_DATA_CHAR])
        return subs

    def parse_barcodes(self, barcodes):

        barcodes_subs = [
            str(s) for s in self.substitutions if str(s) in barcodes.columns
        ]
        # Count up barcode mutations by lineage
        df = copy.copy(barcodes[["lineage"] + barcodes_subs])
        df["total"] = df[barcodes_subs].sum(axis=1)
        summary_df = (
            df[["lineage", "total"]]
            .query("total > 0")
            .sort_values(by="total", ascending=False)
        )
        return summary_df

    def coords_to_ranges(self, attr):
        # Author: @cs95
        # Source: https://stackoverflow.com/a/52302366
        values = getattr(self, attr)
        if len(values) == 0:
            return values

        coords = pd.Series([str(c) for c in values])
        diffs = coords.astype(int).diff().bfill().ne(1).cumsum()
        ranges = (
            coords.groupby(diffs)
            .apply(lambda x: "-".join(x.values[[0, -1]]) if len(x) > 1 else x.item())
            .tolist()
        )
        return ranges

    def ranges_to_coords(self, attr=None, values=None):
        if attr:
            values = getattr(self, attr)
        elif not values:
            raise RebarError("ranges_to_coords: attr or values must be specified.")
        values_split = values.split(",")
        coords = []

        for c in values_split:
            if c == NO_DATA_CHAR:
                continue
            c_split = [int(m) for m in c.split("-")]
            if len(c_split) == 1:
                c = [c_split[0]]
            else:
                c = list(range(c_split[0], c_split[1] + 1))
            coords += c

        return coords

    def to_dataframe(self):
        df = pd.DataFrame(
            {
                "strain": [self.id],
                "substitutions": [",".join([str(s) for s in self.substitutions])],
                "deletions": [",".join(self.coords_to_ranges("deletions"))],
                "missing": [",".join(self.coords_to_ranges("missing"))],
                "genome_length": self.genome_length,
            }
        )
        return df

    def cleanup_conflicts(self, backbone, regions_filter, col):
        subs = copy.copy(getattr(backbone, col))
        lineage = backbone.lineage
        for s in getattr(backbone, col):
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

    def breakpoints(self, parent_2, region_min_subs=1, region_min_len=1):

        result = {}

        parent_1 = self.parent_1

        # Identify barcode subs which uniq to each parent
        all_subs = sorted(list(set(parent_1.barcode + parent_2.barcode)))
        parent_1_subs = [s for s in parent_1.barcode if s not in parent_2.barcode]
        parent_2_subs = [s for s in parent_2.barcode if s not in parent_1.barcode]
        # Store the origin of each mutation
        all_subs_origin = [
            parent_1.lineage
            if s in parent_1_subs
            else parent_2.lineage
            if s in parent_2_subs
            else "shared"
            for s in all_subs
        ]

        # Organize into a dataframe where rows are substitutions
        subs_df = pd.DataFrame(
            {
                "substitution": all_subs,
                "parent": all_subs_origin,
            }
        ).sort_values(by="substitution")

        # Search for genomic blocks from each parent
        # Just look at the subs that are uniq to one parent and detected in sample
        subs_uniq_df = subs_df[
            (subs_df["parent"] != "shared")
            & (subs_df["substitution"].isin(self.substitutions))
        ]
        print(subs_uniq_df)

        # The subs_df has to have at minimum these uniq subs from each parent
        if len(subs_uniq_df) < (region_min_subs * 2):
            return None

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

        print(regions_filter)
        # If we're left with one filtered region, the alt parent fails
        if len(regions_filter) < 2:
            return None

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

        # parent_1_conflict = self.cleanup_conflicts(
        #     parent_1, regions_filter, "conflict_ref"
        # )
        # parent_2_conflict = self.cleanup_conflicts(
        #     parent_2, regions_filter, "conflict_ref"
        # )

        # # Missing subs, hmm

        # # Any remaining conflicts will be private mutations
        # private_subs = parent_1_conflict + parent_2_conflict
        # subs_df = pd.concat(
        #     [
        #         subs_df,
        #         pd.DataFrame(
        #             {
        #                 "substitution": private_subs,
        #                 "parent": ["private"] * len(private_subs),
        #             }
        #         ),
        #     ]
        # ).sort_values(by="substitution")

        # Add bases for reference, parents, and genome
        subs_df["Reference"] = [s.ref for s in all_subs]
        subs_df[parent_1.lineage] = (
            [s.alt if s in parent_1_subs else s.ref for s in all_subs],
        )
        subs_df[parent_2.lineage] = (
            [s.alt if s in parent_2_subs else s.ref for s in all_subs],
        )
        subs_df[self.id] = [
            "N"
            if s.coord in self.missing
            else "-"
            if s.coord in self.deletions
            else s.alt
            if s in self.substitutions
            else s.ref
            for s in all_subs
        ]
        print(subs_df)

        result["subs_df"] = subs_df
        result["regions"] = regions_filter
        result["breakpoints"] = breakpoints

        return result


# class Genome:
#     def __init__(self, id, genome_len, nextclade_row, barcode):
#         self.id = id

#         self.genome_len = genome_len
#         self.subs = self.parse_subs(nextclade_row, col="substitutions")
#         self.privates_reversion = self.parse_subs(
#             nextclade_row, col="privateNucMutations.reversionSubstitutions"
#         )
#         self.privates_labeled = self.parse_subs(
#             nextclade_row, col="privateNucMutations.labeledSubstitutions"
#         )
#         self.privates_unlabeled = self.parse_subs(
#             nextclade_row, col="privateNucMutations.unlabeledSubstitutions"
#         )
#         self.missing = self.parse_intervals(nextclade_row, col="missing")
#         self.deletions = self.parse_intervals(nextclade_row, col="deletions")

#         # Check which substitutions are "barcodes" (lineage-defining in combination)
#         self.barcode_summary = self.identify_barcode(barcode)

#         # Attributes that will be set later
#         self.lineage = None
#         self.recombinant = None
#         self.recursive = False
#         self.backbone = Backbone()

#     def parse_intervals(self, nextclade_row, col):
#         col_str = nextclade_row[col].values[0].split(",")
#         col_coords = []

#         for coords in col_str:
#             if coords == NO_DATA_CHAR:
#                 continue
#             coords_split = [int(m) for m in coords.split("-")]
#             if len(coords_split) == 1:
#                 coords = [coords_split[0]]
#             else:
#                 coords = list(range(coords_split[0], coords_split[1] + 1))
#             col_coords += coords

#         return col_coords

#     def identify_barcode(self, barcode):

#         barcode_subs = [str(s) for s in self.subs if str(s) in barcode.columns]
#         # Count up barcode mutations by lineage
#         df = copy.copy(barcode[["lineage"] + barcode_subs])
#         df["total"] = df[barcode_subs].sum(axis=1)
#         total_df = (
#             df[["lineage", "total"]]
#             .query("total > 0")
#             .sort_values(by="total", ascending=False)
#         )
#         return total_df

#     def parse_subs(self, nextclade_row, col):
#         """Parse substitutions from Nextclade records."""
#         subs_str = nextclade_row[col].values[0].split(",")

#         # labeled mutation is special format
#         if col == "privateNucMutations.labeledSubstitutions":
#             subs_str = [s.split("|")[0] for s in subs_str if "|" in s]

#         subs = sorted([Substitution(s) for s in set(subs_str) if s != NO_DATA_CHAR])
#         return subs

#     def identify_breakpoints(self, alt_backbone, region_min_subs=1, region_min_len=1):

#         result = {}

#         parent_1 = self.backbone.lineage
#         parent_2 = alt_backbone.lineage

#         # Identify the substitutions that are uniq to each parent
#         parent_1_subs = [
#             s for s in self.backbone.support_subs
#             if s not in alt_backbone.support_subs
#         ]
#         parent_2_subs = [
#             s for s in alt_backbone.support_subs
#             if s not in self.backbone.support_subs
#         ]
#         # Identify the substitutions that are shared by both parents
#         shared_subs = [
#             s for s in self.backbone.support_subs if s in alt_backbone.support_subs
#         ]

#         # Organize into a dataframe for easy parsing
#         subs_df = pd.DataFrame(
#             {
#                 "substitution": shared_subs + parent_1_subs + parent_2_subs,
#                 "parent": ["shared"] * len(shared_subs)
#                 + [parent_1] * len(parent_1_subs)
#                 + [parent_2] * len(parent_2_subs),
#             }
#         ).sort_values(by="substitution")
#         subs_df["strain"] = [self.id] * len(subs_df)

#         # Search for genomic blocks for each parent
#         # Just look at the subs that are uniq to one parent
#         subs_uniq_df = subs_df.query("parent != 'shared'")

#         regions = {}
#         p_prev = None
#         start = 0
#         end = 0

#         for rec in subs_uniq_df.iterrows():
#             p_curr = rec[1]["parent"]
#             sub = rec[1]["substitution"]
#             # First region
#             if not p_prev:
#                 start = sub.coord
#                 end = sub.coord
#                 regions[start] = {
#                     "start": start,
#                     "end": end,
#                     "parent": p_curr,
#                     "subs": [sub],
#                 }
#             # Same parent, region continues
#             elif p_curr == p_prev:
#                 regions[start]["end"] = sub.coord
#                 regions[start]["subs"].append(sub)
#             # Parent change, start of new region
#             elif p_curr != p_prev:
#                 start = sub.coord
#                 end = sub.coord
#                 regions[start] = {
#                     "start": start,
#                     "end": end,
#                     "parent": p_curr,
#                     "subs": [sub],
#                 }

#             end = sub.coord
#             p_prev = p_curr

#         # Filter recombination regions
#         #   - At least X consecutive subs per region
#         #   - At least X bases per region
#         regions_filter = {}
#         for start in regions:
#             num_subs = len(regions[start]["subs"])
#             region_len = regions[start]["end"] - start
#             if num_subs < region_min_subs:
#                 continue
#             if region_len < region_min_len:
#                 continue
#             regions_filter[start] = regions[start]

#         # If we're left with one filtered region, the alt parent fails
#         if len(regions_filter) < 2:
#             return 0

#         # Identify breakpoints
#         breakpoints = []
#         prev_start_coord = None
#         prev_end_coord = None

#         for start_coord in regions_filter:

#             end_coord = regions_filter[start_coord]["end"]

#             # Skip the first record for breakpoints
#             if prev_start_coord:
#                 breakpoint_start = prev_end_coord + 1
#                 breakpoint_end = start_coord - 1
#                 breakpoint = "{}:{}".format(breakpoint_start, breakpoint_end)
#                 breakpoints.append(breakpoint)

#             prev_start_coord = start_coord
#             prev_end_coord = end_coord

#         # Check if the missing and conflict subs from parent 1 (self.backbone) can
#         # be explained by recombination from parent 2 (alt_backbone)
#         # parent_1_missing = self.cleanup_conflicts(
#         #     self.backbone, regions_filter, "missing_subs"
#         # )
#         # parent_2_missing = self.cleanup_conflicts(
#         #     alt_backbone, regions_filter, "missing_subs"
#         # )

#         parent_1_conflict = self.cleanup_conflicts(
#             self.backbone, regions_filter, "conflict_subs_ref"
#         )
#         parent_2_conflict = self.cleanup_conflicts(
#             alt_backbone, regions_filter, "conflict_subs_ref"
#         )

#         # Missing subs, hmm

#         # Any remaining conflicts will be private mutations
#         private_subs = parent_1_conflict + parent_2_conflict
#         subs_df = pd.concat(
#             [
#                 subs_df,
#                 pd.DataFrame(
#                     {
#                         "substitution": private_subs,
#                         "parent": ["private"] * len(private_subs),
#                     }
#                 ),
#             ]
#         ).sort_values(by="substitution")

#         result["subs_df"] = subs_df
#         result["regions"] = regions_filter
#         result["breakpoints"] = breakpoints

#         return result

#     def cleanup_conflicts(self, backbone, regions_filter, col):
#         subs = copy.copy(getattr(backbone, col))
#         lineage = backbone.lineage
#         for s in getattr(backbone, col):
#             region_found = False
#             for r in regions_filter.values():
#                 start = r["start"]
#                 end = r["end"]
#                 parent = r["parent"]

#                 # Does it fall within another parent's region?
#                 if s.coord >= start and s.coord <= end:
#                     region_found = True
#                     if lineage != parent:
#                         subs.remove(s)

#             # Does it fall within no known parental region?
#             if not region_found:
#                 subs.remove(s)

#         return subs

#     def write_snipit(self, outdir, exclude_shared=False, ext="pdf"):
#         if exclude_shared:
#             self.recombination_subs = self.recombination_subs.query(
#                 "parent != 'shared'"
#             )

#         parents = list(np.unique(self.recombination_subs["parent"]))
#         if "shared" in parents:
#             parents.remove("shared")
#         if "private" in parents:
#             parents.remove("private")

#         # This will be our crude object to hold fasta lines
#         output_lines = []
#         blank_seq = ["N"] * self.genome_len

#         for strain in ["Reference"] + parents + [self.id]:

#             header = ">" + strain
#             # Start off with all Ns, we'll patch in subs as we go
#             seq = blank_seq

#             for rec in self.recombination_subs.iterrows():
#                 sub = rec[1]["substitution"]
#                 sub_parent = rec[1]["parent"]

#                 # Decide whether we want to write the ref or alt
#                 # mutation, depending on the currently sample
#                 # and the mutation's parental origin

#                 # Reference strain, always use ref
#                 if strain == "Reference":
#                     seq[sub.coord - 1] = sub.ref

#                 # genomic sample, always use alt
#                 elif strain == self.id:
#                     seq[sub.coord - 1] = sub.alt

#                 # Shared by all parents, always use alt
#                 elif sub_parent == "shared":
#                     seq[sub.coord - 1] = sub.alt

#                 # Sub originates from the current strain, use alt
#                 elif sub_parent == strain:
#                     seq[sub.coord - 1] = sub.alt

#                 # Private mutation, and not a genomic sample
#                 elif strain != self.id and sub_parent == "private":
#                     seq[sub.coord - 1] = sub.ref

#                 # otherwise, use ref
#                 else:
#                     seq[sub.coord - 1] = sub.ref

#             output_lines.append(header)
#             output_lines.append("".join(seq))

#         # Write the alignment
#         fasta_name = "snipit_{}_{}.fasta".format(parents[0], parents[1])
#         fasta_path = os.path.join(outdir, fasta_name)
#         with open(fasta_path, "w") as outfile:
#             outfile.write("\n".join(output_lines) + "\n")

#         # Write the figure
#         fig_name = "snipit_{}_{}".format(parents[0], parents[1])
#         fig_path = os.path.join(outdir, fig_name)
#         cmd_str = (
#             "snipit {fasta_path}"
#             " --output-file {fig_path}"
#             " --recombi-mode"
#             " --recombi-references {parent_1},{parent_2}"
#             " --flip-vertical"
#             " --format {ext}"
#         ).format(
#             fasta_path=fasta_path,
#             fig_path=fig_path,
#             parent_1=parents[0],
#             parent_2=parents[1],
#             ext=ext,
#         )
#         subprocess.run(cmd_str, shell=True)

#     def reverse_iter_collapse(
#         regions,
#         min_len,
#         max_breakpoint_len,
#         start_coord,
#         end_coord,
#         parent,
#     ):
#         """Collapse adjacent regions from the same parent into one region."""

#         coord_list = list(regions.keys())
#         coord_list.reverse()

#         for coord in coord_list:
#             prev_start_coord = coord
#             prev_end_coord = regions[prev_start_coord]["end"]
#             prev_region_len = (prev_end_coord - prev_start_coord) + 1
#             prev_parent = regions[coord]["parent"]
#             breakpoint_len = start_coord - prev_end_coord

#             # If the previous region was too short AND from a different parent
#             # Delete that previous region, it's an intermission
#             if prev_region_len < min_len and parent != prev_parent:
#                 del regions[prev_start_coord]

#             # If the previous breakpoint was too long AND from a different parent
#             # Don't add the current region
#             elif (
#                 start_coord != prev_start_coord
#                 and parent != prev_parent
#                 and (max_breakpoint_len != -1 and breakpoint_len > max_breakpoint_len)
#             ):
#                 break

#             # Collapse the current region into the previous one
#             elif parent == prev_parent:
#                 regions[prev_start_coord]["end"] = end_coord
#                 break

#             # Otherwise, parents differ and this is the start of a new region
#             else:
#                 regions[start_coord] = {"parent": parent, "end": end_coord}
#                 break

#         # Check if the reveres iter collapse wound up deleting all the regions
#         if len(regions) == 0:
#             regions[start_coord] = {"parent": parent, "end": end_coord}

#     def __repr__(self):
#         text = (
#             self.id
#             + "\n\trecombinant: "
#             + str(self.recombinant)
#             + "\n\tbackbone:    "
#             + "\n\t\t"
#             + str(self.backbone).replace("\n", "\n\t\t")
#         )
#         return text


def genome_mp(iterator, **kwargs):
    """
    Create genome with multiprocess.
    Used to control the named parameters that are passed.
    """

    # When creating from FASTA, iterator is a SeqRecord
    if type(iterator) == SeqRecord:
        kwargs["record"] = iterator

    # When creating from SUBS df, iterator is a tuple
    elif type(iterator) == tuple:
        # First value is index, second value is series
        kwargs["subs_row"] = iterator[1]

    else:
        raise RebarError("Unknown iterator was passed to genome_mp.")

    genome = Genome(**kwargs)
    return genome
