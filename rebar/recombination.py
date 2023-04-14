import yaml
import pandas as pd
from .barcode import Barcode


class Recombination:
    def __init__(
        self, genome=None, parent_1=None, parent_2=None, min_subs=1, min_length=1
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
            self.search(genome, parent_1, parent_2, min_subs, min_length)

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

    def search(self, genome, parent_1, parent_2, min_subs=1, min_length=1):

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

        # print(list(subs_df["substitution"]))
        # print(list(subs_df["parent"]))
        # Search for genomic blocks from each parent
        # Just look at the subs that are uniq to one parent and detected in sample
        subs_uniq_df = subs_df[
            (subs_df["parent"] != "shared")
            & (subs_df["substitution"].isin(genome.substitutions))
        ]
        # print(subs_uniq_df)

        # The subs_df has to have at minimum these uniq subs from each parent
        if len(subs_uniq_df) < (min_subs * 2):
            return None

        # Identifying breakpoint regions
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
            region_len = (regions[start]["end"] - start) + 1
            if num_subs < min_subs:
                continue
            if region_len < min_length:
                continue
            regions_filter[start] = regions[start]

        # If we're left with one filtered parental region, no recombination
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

        # Add bases for reference, parents, and genome
        subs_df["Reference"] = [s.ref for s in all_subs]
        subs_df[parent_1.lineage] = [
            s.alt if s in parent_1_subs else s.ref for s in all_subs
        ]
        subs_df[parent_2.lineage] = [
            s.alt if s in parent_2_subs else s.ref for s in all_subs
        ]
        subs_df[genome.id] = [
            "N"
            if s.coord in genome.missing
            else "-"
            if s.coord in genome.deletions
            else s.alt
            if s in genome.substitutions
            else s.ref
            for s in all_subs
        ]

        self.dataframe = subs_df
        self.regions = regions_filter
        self.breakpoints = breakpoints

        return 0
