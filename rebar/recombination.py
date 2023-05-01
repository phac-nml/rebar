import yaml
import pandas as pd
from datetime import datetime

from .barcode import Barcode
from .substitution import Substitution


class Recombination:
    def __init__(
        self,
        genome=None,
        parent_1=None,
        parent_2=None,
        max_breakpoints=1,
        min_subs=1,
        min_length=1,
        min_consecutive=1,
    ):

        self.parent_1 = Barcode()
        self.parent_2 = Barcode()
        self.breakpoints = []
        self.regions = {}
        self.dataframe = None
        self.depth = 0

        if parent_1:
            self.parent_1 = parent_1
        if parent_2:
            self.parent_2 = parent_2

        if genome and parent_1 and parent_2:
            self.search(
                genome,
                parent_1,
                parent_2,
                max_breakpoints,
                min_subs,
                min_length,
                min_consecutive,
            )

    def __repr__(self):
        """
        Printable representation of a Recombination object.

        Returns
        -------
        text : str
            String representation.
        """
        return self.to_yaml()

    def to_dict(self):
        recombination_dict = {
            "breakpoints": ",".join(self.breakpoints),
            "regions": ",".join(
                [
                    "{}-{}|{}".format(r["start"], r["end"], r["parent"])
                    for r in self.regions.values()
                ]
            ),
            "parent_1": self.parent_1.to_dict(),
            "parent_2": self.parent_2.to_dict(),
            "depth": self.depth,
        }
        return recombination_dict

    def to_yaml(self, indent=2):
        """
        Convert Recombination object to yaml.

        Returns
        -------
        genome_yaml : yaml
            YAML representation.
        """

        recombination_yaml = (
            yaml.dump(self.to_dict(), sort_keys=False, indent=indent)
            .replace("null", "")
            .replace("''", "")
            + "\n"
        )
        return recombination_yaml

    def search(
        self,
        genome,
        parent_1,
        parent_2,
        max_breakpoints=1,
        min_subs=1,
        min_length=1,
        min_consecutive=1,
    ):

        # ---------------------------------------------------------------------
        # Initialize Barcode Dataframe
        # Create a dataframe where rows are coordinates and columns are
        #   coord, parent, Reference, <parent_1>, <parent_2>, <genome>

        # Identify which subs are non-bi-allelic, these will wind up being
        # duplicate rows, which we'll need to reconcile and collapse
        all_subs = sorted(
            list(set(parent_1.barcode + parent_2.barcode + genome.substitutions))
        )
        all_coords = [s.coord for s in all_subs]
        dup_coords = set([c for c in all_coords if all_coords.count(c) > 1])

        # Re-do all subs list just with parents
        all_subs = sorted(list(set(parent_1.barcode + parent_2.barcode)))

        parent_1_subs = [s for s in parent_1.barcode if s not in parent_2.barcode]
        parent_2_subs = [s for s in parent_2.barcode if s not in parent_1.barcode]

        parent_1_coords = [s.coord for s in parent_1_subs]
        parent_2_coords = [s.coord for s in parent_2_subs]
        genome_coords = [s.coord for s in genome.substitutions]

        # Create the barcode dataframe as described.
        subs_df = pd.DataFrame(
            {
                "coord": [s.coord for s in all_subs],
                "Reference": [s.ref for s in all_subs],
                parent_1.name: [
                    s.alt if s in parent_1_subs else s.ref for s in all_subs
                ],
                parent_2.name: [
                    s.alt if s in parent_2_subs else s.ref for s in all_subs
                ],
                genome.id: [
                    "N"
                    if s.coord in genome.missing
                    else "-"
                    if s.coord in genome.deletions
                    else s.alt
                    if s in genome.substitutions
                    else s.ref
                    for s in all_subs
                ],
            }
        ).sort_values(by="coord")

        # ---------------------------------------------------------------------
        # Collapse duplicate rows from non bi-allelic sites

        for coord in dup_coords:

            # Get base of reference
            ref_base = [s.ref for s in all_subs if s.coord == coord][0]

            # Get base of parent 1
            parent_1_base = ref_base
            if coord in parent_1_coords:
                parent_1_base = [s.alt for s in parent_1_subs if s.coord == coord][0]

            # Get base of parent 2
            parent_2_base = ref_base
            if coord in parent_2_coords:
                parent_2_base = [s.alt for s in parent_2_subs if s.coord == coord][0]

            # If alt's of parent1 and parent2 are same, just exclude, not helpful
            if parent_1_base == parent_2_base:
                # Remove the old duplicate rows
                subs_df = subs_df[subs_df["coord"] != coord]
                continue
            # Otherwise, we'll tidy up and collapse the duplicates
            else:
                # Get base of genomic sample
                genome_base = ref_base
                if coord in genome_coords:
                    genome_base = [
                        s.alt for s in genome.substitutions if s.coord == coord
                    ][0]
                elif coord in genome.missing:
                    genome_base = "N"
                elif coord in genome.deletions:
                    genome_base = "-"

                data = {
                    "coord": [coord],
                    "Reference": [ref_base],
                    parent_1.name: [parent_1_base],
                    parent_2.name: [parent_2_base],
                    genome.id: [genome_base],
                }
                row = pd.DataFrame(data)

                # Remove the old duplicate rows
                subs_df = subs_df[subs_df["coord"] != coord]
                # Add new deduplicated row
                subs_df = pd.concat([subs_df, row]).sort_values(by="coord")

        # Identify private genome substitutions and exclude these
        private_sub_coords = list(
            subs_df[
                (subs_df[genome.id] != subs_df[parent_1.name])
                & (subs_df[genome.id] != subs_df[parent_2.name])
                & (subs_df[genome.id] != subs_df["Reference"])
            ]["coord"]
        )
        subs_df = subs_df[~subs_df["coord"].isin(private_sub_coords)]

        # ---------------------------------------------------------------------
        # Annotate dataframe with parental origin

        # Identify genome sub origins by parent, this is not an efficient method
        genome_subs_origin = []
        for rec in subs_df.iterrows():
            genome_base = rec[1][genome.id]
            parent_1_base = rec[1][parent_1.name]
            parent_2_base = rec[1][parent_2.name]
            if genome_base == parent_1_base and genome_base == parent_2_base:
                origin = "shared"
            elif genome_base == parent_1_base:
                origin = parent_1.name
            elif genome_base == parent_2_base:
                origin = parent_2.name

            genome_subs_origin.append(origin)

        subs_df.insert(loc=1, column="parent", value=genome_subs_origin)

        # ---------------------------------------------------------------------
        # Remove non-discriminating sites

        # Search for genomic blocks from each parent
        # Just look at the subs/barcodes that are uniq to one parent and in sample
        subs_df = subs_df[(subs_df["parent"] != "shared")]

        if genome.debug:
            genome.logger.info(str(datetime.now()) + "\t\t\tBARCODE DISCRIMINATING:")
            subs_md = subs_df.to_markdown(index=False)
            subs_str = subs_md.replace("\n", "\n" + "\t" * 7)
            genome.logger.info(str(datetime.now()) + "\t\t\t\t" + subs_str)

        # Each parent must have at least x min_subs that are lineage-determining
        parent_1_uniq = subs_df[
            (subs_df["parent"] == parent_1.name)
            & (subs_df[parent_1.name] != subs_df["Reference"])
        ]
        parent_1_num_uniq = len(parent_1_uniq)
        parent_2_uniq = subs_df[
            (subs_df["parent"] == parent_2.name)
            & (subs_df[parent_2.name] != subs_df["Reference"])
        ]
        parent_2_num_uniq = len(parent_2_uniq)

        if parent_1_num_uniq < min_subs:
            if genome.debug:
                genome.logger.info(
                    str(datetime.now())
                    + "\t\t\tInsufficient unique substitutions from parent_1: "
                    + parent_1.name
                )
            return None
        if parent_2_num_uniq < min_subs:
            if genome.debug:
                genome.logger.info(
                    str(datetime.now())
                    + "\t\t\tInsufficient unique substitutions from parent_2: "
                    + parent_2.name
                )
            return None

        # ---------------------------------------------------------------------
        # Identify and filter parental regions

        # First: 5' -> 3'
        regions_5p = self.identify_regions(subs_df, genome)
        regions_5p = self.filter_regions_5p(regions_5p, min_consecutive, 0)
        regions_5p = self.filter_regions_5p(regions_5p, 0, min_length)

        if genome.debug:
            genome.logger.info(
                str(datetime.now()) + "\t\t\tREGIONS 5': " + str(regions_5p)
            )

        # Second: 3' to 5'
        regions_3p = self.identify_regions(subs_df, genome)
        regions_3p = dict(reversed(regions_3p.items()))
        regions_3p = self.filter_regions_3p(regions_3p, min_consecutive, 0)
        regions_3p = self.filter_regions_3p(regions_3p, 0, min_length)
        regions_3p = dict(reversed(regions_3p.items()))

        if genome.debug:
            genome.logger.info(
                str(datetime.now()) + "\t\t\tREGIONS 3': " + str(regions_3p)
            )

        # Reconcile 5' vs. 3' differences, by increasing uncertainty
        regions_intersect = self.intersect_regions(regions_5p, regions_3p)

        if genome.debug:
            genome.logger.info(
                str(datetime.now())
                + "\t\t\tREGIONS INTERSECT: "
                + str(regions_intersect)
            )

        # If we're left with one filtered parental region, no recombination
        if len(regions_intersect) < 2:
            if genome.debug:
                genome.logger.info(
                    str(datetime.now()) + "\t\t\t" + "No breakpoints detected."
                )
            return None

        # ---------------------------------------------------------------------
        # Identify breakpoints

        breakpoints = self.identify_breakpoints(regions_intersect)

        if genome.debug:
            genome.logger.info(
                str(datetime.now()) + "\t\t\tBREAKPOINTS: " + str(breakpoints)
            )

        if len(breakpoints) > max_breakpoints:
            if genome.debug:
                genome.logger.info(
                    str(datetime.now()) + "\t\t\tNumber of breakpoints exceeds maximum."
                )
            return None

        # Finish, update class attributes
        self.dataframe = subs_df
        self.regions = regions_intersect
        self.breakpoints = breakpoints

        return 0

    def identify_regions(self, df, genome):
        # Identifying parental regions
        regions = {}
        p_prev = None
        start = 0
        end = 0

        for rec in df.iterrows():
            p_curr = rec[1]["parent"]
            sub = Substitution(
                "{}{}{}".format(rec[1]["Reference"], rec[1]["coord"], rec[1][genome.id])
            )
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

        return regions

    def filter_regions_5p(self, regions, min_consecutive, min_length):
        regions_filter = {}
        prev_start = None
        prev_parent = None

        for start in regions:
            parent = regions[start]["parent"]
            end = regions[start]["end"]
            subs = regions[start]["subs"]
            num_consecutive = len(subs)
            region_length = (end - start) + 1

            # Option 1. First filtered region, or a different parent
            if not prev_parent or prev_parent != parent:
                # Is the new parental region long enough?
                if num_consecutive >= min_consecutive and region_length >= min_length:
                    regions_filter[start] = regions[start]
                    prev_parent = parent
                    prev_start = start
                # Otherwise, continue to next region
                else:
                    continue

            # Option 2. Prev parent continuation
            elif prev_parent == parent:
                # Update end coordinates and subs
                regions_filter[prev_start]["end"] = end
                regions_filter[prev_start]["subs"] += subs
                continue

        return regions_filter

    def filter_regions_3p(self, regions, min_consecutive, min_length):
        regions_filter = {}
        prev_start = None
        prev_parent = None

        for start in regions:
            parent = regions[start]["parent"]
            end = regions[start]["end"]
            subs = regions[start]["subs"]
            num_consecutive = len(subs)
            region_length = (end - start) + 1

            # First filtered region, or different parent from previous region
            if not prev_parent or prev_parent != parent:
                # Is the new parental region long enough?
                if num_consecutive >= min_consecutive and region_length >= min_length:
                    regions_filter[start] = regions[start]
                    prev_parent = parent
                    prev_start = start
                else:
                    continue

            # A region that continues the previous parent
            # intermissions from other parents were skipped over
            elif prev_parent == parent:
                # Update the previous regions coordinates
                regions_filter[start] = regions[prev_start]
                regions_filter[start]["subs"] = sorted(
                    regions_filter[prev_start]["subs"] + subs
                )
                regions_filter[start]["start"] = regions_filter[start]["subs"][0].coord
                regions_filter[start]["end"] = regions_filter[start]["subs"][-1].coord
                regions_filter.pop(prev_start)
                prev_start = start
                continue

        return regions_filter

    def intersect_regions(self, regions_1, regions_2):
        regions_intersect = {}
        for r1 in regions_1.values():

            r1_parent = r1["parent"]
            r1_subs = set(r1["subs"])

            for r2 in regions_2.values():

                r2_parent = r2["parent"]
                if r1_parent != r2_parent:
                    continue

                r2_subs = set(r2["subs"])
                subs_intersect = r1_subs.intersection(r2_subs)
                if len(subs_intersect) == 0:
                    continue

                start = min(subs_intersect).coord
                end = max(subs_intersect).coord
                regions_intersect[start] = {
                    "start": start,
                    "end": end,
                    "parent": r1_parent,
                    "subs": sorted(subs_intersect),
                }

        return regions_intersect

    def identify_breakpoints(self, regions):
        breakpoints = []
        prev_start_coord = None
        prev_end_coord = None

        for start_coord in regions:

            end_coord = regions[start_coord]["end"]

            # Skip the first record for breakpoints
            if prev_start_coord:
                breakpoint_start = prev_end_coord + 1
                breakpoint_end = start_coord - 1
                breakpoint = "{}:{}".format(breakpoint_start, breakpoint_end)
                breakpoints.append(breakpoint)

            prev_start_coord = start_coord
            prev_end_coord = end_coord

        return breakpoints
