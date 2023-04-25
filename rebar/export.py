import os
import subprocess
import pandas as pd


class Export:
    def __init__(self, genomes, outdir, include_shared):
        self.genomes = genomes
        self.outdir = outdir
        self.parents = self.collect_parents()
        self.include_shared = include_shared
        self.dataframes = {}
        self.alignments = {}

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

    def collect_parents(self):
        result = {}

        for genome in self.genomes:
            # No recombination detected, skip
            if len(genome.recombination.breakpoints) == 0:
                continue
            # Recombination detected, proceed
            parents = sorted(
                [
                    genome.recombination.parent_1.name,
                    genome.recombination.parent_2.name,
                ]
            )
            parents = "{}_{}".format(parents[0], parents[1])
            if parents not in result:
                result[parents] = []
            result[parents].append(genome)
        return result

    def to_barcodes(self):

        for parents in self.parents:

            parent_1 = parents.split("_")[0]
            parent_2 = parents.split("_")[1]

            # First pass, identify all ref and parent subs
            parent_dict = {
                "coord": [],
                "Reference": [],
                parent_1: [],
                parent_2: [],
            }
            for genome in self.parents[parents]:
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
            for genome in self.parents[parents]:
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

            self.dataframes[parents] = parent_df
            file_path = os.path.join(self.outdir, "barcode_{}.tsv".format(parents))
            parent_df.to_csv(file_path, sep="\t", index=False)

    def to_alignments(self):

        # Create dataframes first
        if len(self.dataframes) == 0:
            self.to_barcodes()

        for parents in self.parents:

            # Get genome length from first genome
            genome_length = self.parents[parents][0].genome_length
            df = self.dataframes[parents]

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

            file_path = os.path.join(self.outdir, "alignment_{}.fasta".format(parents))
            with open(file_path, "w") as outfile:
                outfile.write("\n".join(aln_lines) + "\n")

            self.alignments[parents] = file_path

    def to_snipit(self, ext):
        # Create alignments first
        if len(self.alignments) == 0:
            self.to_alignments()

        for parents in self.parents:
            # Ex. "BJ.1_CJ.1.1"
            parent_1 = parents.split("_")[0]
            parent_2 = parents.split("_")[1]
            fasta_path = self.alignments[parents]
            fig_prefix = os.path.join(self.outdir, "snipit_{}".format(parents))

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
