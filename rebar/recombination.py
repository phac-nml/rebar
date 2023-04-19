import yaml
import pandas as pd
from .barcode import Barcode
from datetime import datetime
from .substitution import Substitution


class Recombination:
    def __init__(
        self,
        genome=None,
        parent_1=None,
        parent_2=None,
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
                genome, parent_1, parent_2, min_subs, min_length, min_consecutive
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
        self, genome, parent_1, parent_2, min_subs=1, min_length=1, min_consecutive=1
    ):

        # Identify barcode subs which uniq to each parent
        all_subs = sorted(list(set(parent_1.barcode + parent_2.barcode)))
        parent_1_subs = [s for s in parent_1.barcode if s not in parent_2.barcode]
        parent_2_subs = [s for s in parent_2.barcode if s not in parent_1.barcode]

        # Organize into a dataframe where rows are substitutions, and columns
        # will indicate strain and parental origin
        subs_df = pd.DataFrame(
            {
                # "substitution": all_subs,
                "coord": [s.coord for s in all_subs],
                "Reference": [s.ref for s in all_subs],
                parent_1.lineage: [
                    s.alt if s in parent_1_subs else s.ref for s in all_subs
                ],
                parent_2.lineage: [
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

        # Identify private genome substitutions and exclude these
        private_sub_coords = list(
            subs_df[
                (subs_df[genome.id] != subs_df[parent_1.lineage])
                & (subs_df[genome.id] != subs_df[parent_2.lineage])
                & (subs_df[genome.id] != subs_df["Reference"])
            ]["coord"]
        )
        subs_df = subs_df[~subs_df["coord"].isin(private_sub_coords)]

        # Identify genome sub origins by parent, this is not an efficient method
        genome_subs_origin = []
        for rec in subs_df.iterrows():
            genome_base = rec[1][genome.id]
            parent_1_base = rec[1][parent_1.lineage]
            parent_2_base = rec[1][parent_2.lineage]
            if genome_base == parent_1_base and genome_base == parent_2_base:
                origin = "shared"
            elif genome_base == parent_1_base:
                origin = parent_1.lineage
            elif genome_base == parent_2_base:
                origin = parent_2.lineage

            genome_subs_origin.append(origin)

        subs_df.insert(loc=1, column="parent", value=genome_subs_origin)

        # Search for genomic blocks from each parent
        # Just look at the subs that are uniq to one parent and detected in sample
        subs_uniq_df = subs_df[(subs_df["parent"] != "shared")]

        if genome.debug:
            genome.logger.info(str(datetime.now()) + "\t\t\tBARCODE ORIGINS:")
            subs_md = subs_uniq_df.to_markdown(index=False)
            subs_str = subs_md.replace("\n", "\n" + "\t" * 7)
            genome.logger.info(str(datetime.now()) + "\t\t\t\t" + subs_str)

        # Each parent must have at least x min_subs that are lineage-determining
        parent_1_uniq = subs_df[
            (subs_df["parent"] == parent_1.lineage)
            & (subs_df[parent_1.lineage] != subs_df["Reference"])
            & (subs_df[parent_1.lineage] != subs_df[parent_2.lineage])
        ]
        parent_1_num_uniq = len(parent_1_uniq)
        parent_2_uniq = subs_df[
            (subs_df["parent"] == parent_2.lineage)
            & (subs_df[parent_2.lineage] != subs_df["Reference"])
            & (subs_df[parent_2.lineage] != subs_df[parent_1.lineage])
        ]
        parent_2_num_uniq = len(parent_2_uniq)

        if (parent_1_num_uniq < min_subs) or (parent_2_num_uniq < min_subs):
            if genome.debug:
                genome.logger.info(
                    str(datetime.now()) + "\t\t\tInsufficient unique substitutions."
                )
            return None

        # Identifying breakpoint regions
        regions = {}
        p_prev = None
        start = 0
        end = 0

        for rec in subs_uniq_df.iterrows():
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

        # Filter recombination regions
        #   - At least X consecutive subs per region
        #   - At least X bases per region
        regions_filter = {}
        prev_start = None
        prev_parent = None

        for start in regions:
            parent = regions[start]["parent"]
            end = regions[start]["end"]
            num_consecutive = len(regions[start]["subs"])
            region_length = (end - start) + 1

            if region_length < min_length or num_consecutive < min_consecutive:
                continue

            # First filtered region
            if not prev_parent:
                regions_filter[start] = regions[start]

            # A region that continues the previous parent
            # intermissions from other parents were skipped over
            elif prev_parent == parent:
                # Update the end coordinates
                regions_filter[prev_start]["end"] = end
                continue
            elif prev_parent != parent:
                # start new region
                regions_filter[start] = regions[start]

            prev_parent = parent
            prev_start = start

        if genome.debug:
            genome.logger.info(
                str(datetime.now()) + "\t\t\tREGIONS: " + str(regions_filter)
            )

        # If we're left with one filtered parental region, no recombination
        if len(regions_filter) < 2:
            if genome.debug:
                genome.logger.info(
                    str(datetime.now()) + "\t\t\t" + "No breakpoints detected."
                )
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

        if genome.debug:
            genome.logger.info(
                str(datetime.now()) + "\t\t\tBREAKPOINTS: " + str(breakpoints)
            )

        self.dataframe = subs_df
        self.regions = regions_filter
        self.breakpoints = breakpoints

        return 0
