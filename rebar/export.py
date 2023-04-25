import os
import subprocess
import pandas as pd

from .plot import plot


class Export:
    def __init__(self, genomes, outdir, include_shared):
        self.genomes = genomes
        self.outdir = outdir
        self.recombinants = self.collect_recombinants()
        self.include_shared = include_shared
        self.dataframe = None
        self.dataframes = {}
        self.alignments = {}

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

                parents = sorted(
                    [
                        genome.recombination.parent_1.name,
                        genome.recombination.parent_2.name,
                    ]
                )
                parents = "{}_{}".format(parents[0], parents[1])
                if parents not in result[recombinant]:
                    result[recombinant][parents] = []
                result[recombinant][parents].append(genome)

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

        return dataframe

    def to_barcodes(self):

        # Create output dir
        outdir_barcodes = os.path.join(self.outdir, "barcodes")
        if not os.path.exists(outdir_barcodes):
            os.makedirs(outdir_barcodes)

        # Process recombinant groups
        for recombinant in self.recombinants:
            self.dataframes[recombinant] = {}

            # Process parent groups within recombinant
            for parents in self.recombinants[recombinant]:

                parents_data = self.recombinants[recombinant][parents]
                parent_1 = parents.split("_")[0]
                parent_2 = parents.split("_")[1]

                # First pass, identify all ref and parent subs
                parent_dict = {
                    "coord": [],
                    "Reference": [],
                    parent_1: [],
                    parent_2: [],
                }
                for genome in parents_data:
                    df = genome.recombination.dataframe
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
                for genome in parents_data:
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

                self.dataframes[recombinant][parents] = parent_df
                file_path = os.path.join(
                    outdir_barcodes, "{}_{}.tsv".format(recombinant, parents)
                )
                parent_df.to_csv(file_path, sep="\t", index=False)

    def to_alignments(self):

        # Create output dir
        outdir_alignments = os.path.join(self.outdir, "alignments")
        if not os.path.exists(outdir_alignments):
            os.makedirs(outdir_alignments)

        # Create dataframes first
        if len(self.dataframes) == 0:
            self.to_barcodes()

        # Process recombinant groups
        for recombinant in self.recombinants:
            self.alignments[recombinant] = {}

            # Process parent groups within recombinant
            for parents in self.recombinants[recombinant]:

                parents_data = self.recombinants[recombinant][parents]

                # Get genome length from first genome
                genome_length = parents_data[0].genome_length
                df = self.dataframes[recombinant][parents]

                coords = [coord for coord in list(df["coord"])]
                aln_lines = []
                # Construct a pseudo-alignment where every position is 'N'
                # unless it's a barcode substitution in the dataframe
                for header in df.columns[1:]:
                    seq = ["N"] * genome_length
                    bases = list(df[header])
                    for coord, base in zip(coords, bases):
                        # genomic coordinates are 1-based, adjust for 0-based index
                        seq[coord - 1] = base
                        # if coord == 25810:
                        #    print(base)

                    aln_lines.append(">" + header)
                    aln_lines.append("".join(seq))

                file_path = os.path.join(
                    outdir_alignments, "{}_{}.fasta".format(recombinant, parents)
                )
                with open(file_path, "w") as outfile:
                    outfile.write("\n".join(aln_lines) + "\n")

                self.alignments[recombinant][parents] = file_path

    def to_snipit(self, ext):

        # Create output dir
        outdir_snipit = os.path.join(self.outdir, "snipit")
        if not os.path.exists(outdir_snipit):
            os.makedirs(outdir_snipit)

        # Create alignments first
        if len(self.alignments) == 0:
            self.to_alignments()

        # Process recombinant groups
        for recombinant in self.recombinants:

            # Process parent groups within recombinant
            for parents in self.recombinants[recombinant]:

                # Ex. "BJ.1_CJ.1.1"
                parent_1 = parents.split("_")[0]
                parent_2 = parents.split("_")[1]
                fasta_path = self.alignments[recombinant][parents]
                fig_prefix = os.path.join(
                    outdir_snipit, "{}_{}".format(recombinant, parents)
                )

                cmd_str = (
                    "snipit {fasta_path}"
                    " --output-file {fig_prefix}"
                    " --recombi-mode"
                    " --recombi-references {parent_1},{parent_2}"
                    " --reference Reference"
                    " --flip-vertical"
                    " --format {ext}"
                    " --solid-background"
                ).format(
                    fasta_path=fasta_path,
                    fig_prefix=fig_prefix,
                    parent_1=parent_1,
                    parent_2=parent_2,
                    ext=ext,
                )
                subprocess.run(cmd_str, shell=True)

    def to_plot(self):

        # Create output dir
        outdir_plot = os.path.join(self.outdir, "plot")
        if not os.path.exists(outdir_plot):
            os.makedirs(outdir_plot)

        # Create summary dataframe first
        if type(self.dataframe) != pd.core.series.Series:
            self.to_dataframe()

        # Create individual dataframes
        if len(self.dataframes) == 0:
            self.to_barcodes()

        # Process recombinant groups
        for recombinant in self.recombinants:
            self.alignments[recombinant] = {}

            # Process parent groups within recombinant
            for parents in self.recombinants[recombinant]:
                barcodes_path = os.path.join(
                    self.outdir, "barcodes", "{}_{}.tsv".format(recombinant, parents)
                )
                summary_path = os.path.join(self.outdir, "summary.tsv")
                output_path = os.path.join(
                    self.outdir, "plots", "{}_{}.png".format(recombinant, parents)
                )
                plot(barcodes=barcodes_path, summary=summary_path, output=output_path)

        return 0
