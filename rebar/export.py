import os
import pandas as pd

from .plot import plot


class Export:
    def __init__(
        self,
        genomes,
        dataset,
        outdir,
    ):
        self.genomes = genomes
        self.dataset = dataset
        self.outdir = outdir
        self.recombinants = self.collect_recombinants()
        self.dataframe = None
        self.barcodes = {}
        self.alignments = {}
        self.annotations = self.load_annotations(dataset)

    def load_annotations(self, dataset):
        """
        Load annotations dataframe from dataset path.
        """
        annot_path = os.path.join(self.dataset, "annotations.tsv")
        annot_df = pd.read_csv(annot_path, sep="\t")
        return annot_df

    def collect_recombinants(self):
        """
        Collate genomes by recombinant and parents.
        """
        recombinants = list(
            set([g.lineage.recombinant for g in self.genomes if g.lineage.recombinant])
        )
        result = {r: {} for r in recombinants}

        for recombinant in recombinants:

            for genome in self.genomes:
                if genome.lineage.recombinant != recombinant:
                    continue

                parents = [
                    genome.recombination.parent_1.name,
                    genome.recombination.parent_2.name,
                ]
                # No parent 2, no recombation
                if not parents[1]:
                    continue

                parents = "{}_{}".format(parents[0], parents[1])
                if parents not in result[recombinant]:
                    result[recombinant][parents] = {}

                breakpoints_str = "_".join(genome.recombination.breakpoints)

                if breakpoints_str not in result[recombinant][parents]:
                    result[recombinant][parents][breakpoints_str] = []
                result[recombinant][parents][breakpoints_str].append(genome)

        return result

    def to_yaml(self):
        yaml_data = "\n".join([genome.to_yaml() for genome in self.genomes])
        file_path = os.path.join(self.outdir, "summary.yaml")
        with open(file_path, "w") as outfile:
            outfile.write(yaml_data + "\n")

    def to_dataframe(self):
        dataframe = pd.DataFrame()
        for genome in self.genomes:
            genome_dataframe = genome.to_dataframe(df_type="full")
            if len(dataframe) == 0:
                dataframe = genome_dataframe
            else:
                dataframe = pd.concat([dataframe, genome_dataframe], ignore_index=True)

        file_path = os.path.join(self.outdir, "summary.tsv")
        dataframe.to_csv(file_path, sep="\t", index=False)
        self.dataframe = dataframe

        return dataframe

    def to_barcodes(self):

        # Create output dir
        outdir_barcodes = os.path.join(self.outdir, "barcodes")
        if not os.path.exists(outdir_barcodes):
            os.makedirs(outdir_barcodes)

        # Process recombinant groups
        for recombinant in self.recombinants:
            self.barcodes[recombinant] = {}

            # Process parent groups within recombinant
            for parents in self.recombinants[recombinant]:
                self.barcodes[recombinant][parents] = {}

                parents_data = self.recombinants[recombinant][parents]
                parent_1 = parents.split("_")[0]
                parent_2 = parents.split("_")[1]

                # No parent 2 means no recombination
                if not parent_2:
                    continue

                # Process breakpoint groups within parents
                for breakpoints in parents_data:

                    # First pass, identify all ref and parent subs
                    parent_dict = {
                        "coord": [],
                        "Reference": [],
                        parent_1: [],
                        parent_2: [],
                    }

                    breakpoints_data = parents_data[breakpoints]
                    for genome in breakpoints_data:
                        df = genome.recombination.dataframe
                        # If we know this is a recombinant, but couldn't detect
                        # breakpoints, skip over
                        if type(df) != pd.core.frame.DataFrame:
                            continue
                        for rec in df.iterrows():
                            coord = rec[1]["coord"]
                            ref = rec[1]["Reference"]
                            p1 = rec[1][parent_1]
                            p2 = rec[1][parent_2]

                            if coord not in parent_dict["coord"]:
                                parent_dict["coord"].append(coord)
                                parent_dict["Reference"].append(ref)
                                parent_dict[parent_1].append(p1)
                                parent_dict[parent_2].append(p2)

                    parent_df = pd.DataFrame(parent_dict).sort_values(by=["coord"])

                    # Second pass, add sample subs
                    for genome in breakpoints_data:
                        df = genome.recombination.dataframe
                        bases = []

                        for rec in parent_df.iterrows():
                            coord = rec[1]["coord"]
                            ref = rec[1]["Reference"]
                            # If this coord is missing, set to ref
                            base = ref
                            coord_row = df[df["coord"] == coord]
                            if len(coord_row) == 0:
                                base = ref
                            else:
                                base = coord_row[genome.id].values[0]

                            bases.append(base)

                        parent_df[genome.id] = bases

                    self.barcodes[recombinant][parents][breakpoints] = parent_df
                    file_path = os.path.join(
                        outdir_barcodes,
                        "{}_{}_{}.tsv".format(recombinant, parents, breakpoints),
                    )
                    parent_df.to_csv(file_path, sep="\t", index=False)

        return self.barcodes

    def to_plot(self, ext):

        # Create output dir
        outdir_plot = os.path.join(self.outdir, "plots")
        if not os.path.exists(outdir_plot):
            os.makedirs(outdir_plot)

        # Create summary dataframe first
        if type(self.dataframe) != pd.core.series.Series:
            self.to_dataframe()

        # Create individual barcodes
        if len(self.barcodes) == 0:
            self.to_barcodes()

        # Process recombinant groups
        for recombinant in self.recombinants:
            # print(recombinant)
            self.alignments[recombinant] = {}

            # Process parent groups within recombinant
            for parents in self.recombinants[recombinant]:
                # print("\t", parents)
                parents_data = self.recombinants[recombinant][parents]

                # In the summary table, parents are seperated by comma
                parents_csv = parents.replace("_", ",")

                for breakpoints in parents_data:
                    # Get barcodes and summary for this recombinant
                    barcodes_df = self.barcodes[recombinant][parents][breakpoints]

                    # If we know this is a recombinant, but couldn't detect
                    # breakpoints, skip over
                    if type(barcodes_df) != pd.core.frame.DataFrame:
                        continue

                    summary_df = self.dataframe[
                        (self.dataframe["recombinant"] == recombinant)
                        & (self.dataframe["parents_lineage"] == parents_csv)
                        * (
                            self.dataframe["breakpoints"]
                            == breakpoints.replace("_", ",")
                        )
                    ]
                    output_path = os.path.join(
                        self.outdir,
                        "plots",
                        "{}_{}_{}.{}".format(recombinant, parents, breakpoints, ext),
                    )
                    plot(
                        barcodes_df=barcodes_df,
                        summary_df=summary_df,
                        annot_df=self.annotations,
                        output=output_path,
                    )

        return 0
