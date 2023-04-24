# Standard Libraries
import yaml
from datetime import datetime
from copy import copy

# PyPI libraries
import pandas as pd
from Bio.SeqRecord import SeqRecord
from pango_aliasor.aliasor import Aliasor

# rebar custom
from .constants import NO_DATA_CHAR
from .substitution import Substitution
from .barcode import Barcode
from .recombination import Recombination
from .edge_cases import handle_edge_cases
from . import RebarError


class Genome:
    """
    Genomes defines a genomic sample object with interface methods for parsing
    features from the sequence (substitutions, deletions, missing data), summarizing
    matches to lineage barcodes, and identifying recombinant parents and breakpoints.
    """

    def __init__(
        self,
        record=None,
        reference=None,
        subs_row=None,
        barcodes=None,
        tree=None,
        recombinant_tree=None,
        recombinant_lineages=None,
        lineage_to_clade=None,
        max_depth=1,
        min_subs=1,
        min_length=1,
        min_consecutive=1,
        mask=0,
        debug=False,
        logger=None,
        edge_cases=False,
        validate=None,
        dataset_info=None,
    ):
        """
        Genome constructor. Parses genomic features from a sequence records or
        substitutions table.

        Parameters
        ----------
        record : Bio.SeqRecord.SeqRecord
            Sequence record of single sample from multiple alignment.

        reference : Bio.SeqRecord.SeqRecord
            Sequence record of reference genome.

        subs_row : pd.core.series.Series
            Row from substitutions dataframe, either from `rebar subs` subcommand
            or Nextclade TSV.

        barcodes : pd.core.frame.DataFrame
            Dataframe of lineage barcodes, from `rebar barcodes` subcommand
        """
        # Debugging option
        self.debug = debug
        self.logger = logger

        # Generic genomic features
        self.id = None
        self.seq = None
        self.substitutions = []
        self.deletions = []
        self.missing = []
        self.genome_length = None

        # Dataset information
        self.dataset = dataset_info
        # Lineage features
        self.barcode_summary = None
        self.lineage = Barcode()

        # Recombination features
        self.recombination = Recombination()
        self.validate = None

        # Entry point #1, from fasta alignment
        if record:
            self.id = record.id
            # This debugging is very noisy, and not useful?
            # if self.debug:
            #     self.logger.info(
            #         str(datetime.now()) + "\t\tParsing sample: " + record.id
            #     )
            self.seq = str(record.seq)

            # Mask genome sequence
            self.seq = "".join(
                (["N"] * mask)  # 5' masking
                + [self.seq[mask:-mask]]  # in between, unmasked bases
                + ["N"] * mask  # 3' masking
            )

        if reference and self.seq:
            # Mask genome sequence
            reference.seq = str(reference.seq)
            reference.seq = "".join(
                (["N"] * mask)  # 5' masking
                + [reference.seq[mask:-mask]]  # in between, unmasked bases
                + ["N"] * mask  # 3' masking
            )
            self.parse_sequence(reference)

        # Entry point #2, from subs dataframe
        if type(subs_row) == pd.core.series.Series:
            self.id = subs_row["strain"]
            if self.debug:
                self.logger.info(str(datetime.now()) + "\tParsing sample: " + self.id)
            self.genome_length = subs_row["genome_length"]
            self.substitutions = self.parse_substitutions(subs_row=subs_row)
            self.deletions = self.ranges_to_coords(values=subs_row["deletions"])
            self.missing = self.ranges_to_coords(values=subs_row["missing"])

        # Check which substitutions are "barcodes" (lineage-defining in combination)
        if type(barcodes) == pd.core.frame.DataFrame:
            self.barcode_summary = self.summarise_barcodes(barcodes)

        # Perform lineage and parent assignment
        if (
            type(self.barcode_summary) == pd.core.frame.DataFrame
            and type(lineage_to_clade) == pd.core.frame.DataFrame
            and tree
            and recombinant_lineages
            and recombinant_tree
        ):
            if self.debug:
                self.logger.info(str(datetime.now()) + "\t\t" + "LINEAGE ASSIGNMENT:")
            self.lineage = self.lineage_assignment(
                barcode_summary=self.barcode_summary,
                barcodes=barcodes,
                tree=tree,
                recombinant_lineages=recombinant_lineages,
                recombinant_tree=recombinant_tree,
                lineage_to_clade=lineage_to_clade,
            )

            self.parent_assignment(
                barcodes=barcodes,
                tree=tree,
                recombinant_lineages=recombinant_lineages,
                recombinant_tree=recombinant_tree,
                lineage_to_clade=lineage_to_clade,
                max_depth=max_depth,
                min_subs=min_subs,
                min_length=min_length,
                min_consecutive=min_consecutive,
                edge_cases=edge_cases,
            )

        # Validate
        if self.validate and tree and recombinant_lineages:
            self.validate = self.validate_recombination(tree, recombinant_lineages)

    def __repr__(self):
        """
        Printable representation of a Genome object.

        Returns
        -------
        text : str
            String representation.
        """
        # return(self.to_yaml())
        return self.id

    def parse_sequence(self, reference):
        """
        Parse genomic features from sequence.

        Parameters
        ----------
        reference : Bio.SeqRecord.SeqRecord
            Sequence record of reference genome.

        Attributes Modified
        -------
        self.substitutions : list
        self.deletions     : list
        self.missing       : list
        self.genome_length : int
        """

        coord = 0

        for i, bases in enumerate(zip(reference.seq, self.seq)):
            r = bases[0]
            s = bases[1]
            # Genomic coordinates are 1 based
            coord = i + 1

            # print(i, r, s)

            # Missing Data
            if s == "N":
                self.missing.append(coord)

            # Deletions
            elif s == "-":
                self.deletions.append(coord)

            # Substitution, missing ref data
            elif r == "N":
                continue

            # Substitution, true
            elif r != s:
                sub = "{ref}{coord}{alt}".format(ref=r, coord=coord, alt=s)
                self.substitutions.append(Substitution(sub))
            next

        self.genome_length = coord
        return 0

    def parse_substitutions(self, subs_row):
        """
        Parse substitutions column from the subs dataframe row.

        Parameters
        ----------
        subs_row : pd.core.series.Series
            Row from substitutions dataframe, either from `rebar subs` subcommand
            or Nextclade TSV.

        Returns
        -------
        features : list
            List of Substitution objects
        """
        subs_str = subs_row["substitutions"].split(",")
        subs = sorted([Substitution(s) for s in set(subs_str) if s != NO_DATA_CHAR])
        return subs

    def summarise_barcodes(self, barcodes, barcodes_subs=None):
        """
        Summarise detected barcode substitutions.

        Parameters
        ----------
        barcodes : pd.core.frame.DataFrame
            Dataframe of lineage barcodes, from `rebar barcodes` subcommand

        Returns
        -------
        summary_df : pd.core.frame.DataFrame
            Dataframe with columns 'lineage' and 'total' to summarize barcode
            substitution detections.
        """
        if not barcodes_subs:
            barcodes_subs = [
                str(s) for s in self.substitutions if str(s) in barcodes.columns
            ]
        else:
            barcodes_subs = [
                str(s) for s in barcodes_subs if str(s) in barcodes.columns
            ]

        # Count up barcode mutations by lineage
        cols = ["lineage"] + barcodes_subs
        df = copy(barcodes[cols])
        df["total"] = df[barcodes_subs].sum(axis=1)
        summary_df = (
            df[["lineage", "total"]]
            .query("total > 0")
            .sort_values(by=["total", "lineage"], ascending=False)
        )

        return summary_df

    def coords_to_ranges(self, attr=None, values=None):
        """
        Convert list of coordinates to ranges.

        Parameters
        ----------
        attr : str
            Genome attribute name to convert (ex. deletions, missing)

        Examples
        --------
        self.deletions = [1,2,3,10,12]
        coords_to_ranges(attr="deletions")
        ['1-3', '10', '12']

        Returns
        -------
        ranges : list
            List of ranges in string representation.
        """
        # Author: @cs95
        # Source: https://stackoverflow.com/a/52302366
        if attr:
            values = getattr(self, attr)
        elif not values:
            raise SystemExit(
                RebarError(
                    "RebarError: coords_to_ranges: attr or values must be specified."
                )
            )

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
        """
        Convert string representation of ranges to coordinates.

        Parameters
        ----------
        attr : str
            Genome attribute name to convert (ex. deletions, missing)
        values : list
            List of ranges in string representation (ex. ['1-3', '10',])

        Examples
        --------
        self.deletions = ['1-3', '10', '12']
        ranges_to_coords(attr="deletions")
        [1,2,3,10,12]

        Returns
        -------
        coords : list
            List of coordinates as integers.
        """
        if attr:
            values = getattr(self, attr)
        elif not values:
            raise SystemExit(
                RebarError(
                    "RebarError: ranges_to_coords: attr or values must be specified."
                )
            )
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

    def simplify_parents_lineage(self):
        """
        Simplify parental lin

        Examples
        --------
        "BA.5.2,BE.4.1" --> "BA.5.2,BA.5.3"
        """
        parent_1 = self.recombination.parent_1.lineage
        parent_2 = self.recombination.parent_2.lineage

        aliasor = Aliasor()
        parent_1_uncompress = aliasor.uncompress(parent_1)
        parent_2_uncompress = aliasor.uncompress(parent_2)

        # Compress the parts that are shared
        parent_1_split = parent_1_uncompress.split(".")
        parent_2_split = parent_2_uncompress.split(".")
        if len(parent_1_split) < len(parent_2_split):
            min_level = len(parent_1_split)
        else:
            min_level = len(parent_2_split)

        # Find the level in which they differ
        for level in range(0, min_level):
            if parent_1_split[level] != parent_2_split[level]:
                break
        up_to_level = level // 3
        level_remainder = level % 3

        parent_1_compress = aliasor.partial_compress(
            parent_1_uncompress, up_to=up_to_level
        )
        parent_2_compress = aliasor.partial_compress(
            parent_2_uncompress, up_to=up_to_level
        )

        parent_1_simplify = ".".join(
            parent_1_compress.split(".")[0 : (level_remainder + 1)]
        )
        parent_2_simplify = ".".join(
            parent_2_compress.split(".")[0 : (level_remainder + 1)]
        )

        return parent_1_simplify, parent_2_simplify

    def to_dataframe(self, df_type="subs"):
        """
        Convert Genome object to dataframe.

        Returns
        -------
        genome_dataframe : pd.core.frame.DataFrame
            Dataframe representation of genome.
        """
        if df_type == "subs":
            genome_dataframe = pd.DataFrame(
                {
                    "strain": [self.id],
                    "substitutions": [",".join([str(s) for s in self.substitutions])],
                    "deletions": [",".join(self.coords_to_ranges("deletions"))],
                    "missing": [",".join(self.coords_to_ranges("missing"))],
                    "genome_length": self.genome_length,
                }
            )
        else:
            recombination_dict = self.recombination.to_dict()

            # only write parents if recombination detected:
            if not self.recombination.parent_2.lineage:
                parents_lineage = ""
                parents_lineage_simplify = ""
                parents_clade = ""
            else:
                parents_lineage = "{},{}".format(
                    self.recombination.parent_1.lineage,
                    self.recombination.parent_2.lineage,
                )
                parent_1_simplify, parent_2_simplify = self.simplify_parents_lineage()
                parents_lineage_simplify = "{},{}".format(
                    parent_1_simplify,
                    parent_2_simplify,
                )
                parents_clade = "{},{}".format(
                    self.recombination.parent_1.clade,
                    self.recombination.parent_2.clade,
                )

            genome_dataframe = pd.DataFrame(
                {
                    "strain": [self.id],
                    "lineage": [self.lineage.lineage],
                    "clade": [self.lineage.clade],
                    "recombinant": [self.lineage.recombinant],
                    "definition": [self.lineage.definition],
                    "validate": [self.validate],
                    "parents_lineage": parents_lineage,
                    "parents_lineage_simplify": parents_lineage_simplify,
                    "parents_clade": parents_clade,
                    "breakpoints": recombination_dict["breakpoints"],
                    "regions": recombination_dict["regions"],
                    "dataset_name": self.dataset["name"],
                    "dataset_tag": self.dataset["tag"],
                    "barcodes_date": self.dataset["barcodes"]["date"],
                    "barcodes_tag": self.dataset["barcodes"]["tag"],
                    "tree_date": self.dataset["tree"]["date"],
                    "tree_tag": self.dataset["tree"]["tag"],
                    "sequences_date": self.dataset["sequences"]["date"],
                    "sequences_tag": self.dataset["sequences"]["tag"],
                }
            )

        return genome_dataframe

    def to_dict(self):
        """
        Convert Genome object to dict.

        Returns
        -------
        genome_dict : dict
            Dictionary representation of genome.
        """
        genome_dict = {
            self.id: {
                "substitutions": ",".join([str(s) for s in self.substitutions]),
                "deletions": ",".join(self.coords_to_ranges("deletions")),
                "missing": ",".join(self.coords_to_ranges("missing")),
                "lineage": self.lineage.to_dict(),
                "recombination": self.recombination.to_dict(),
            }
        }
        return genome_dict

    def to_yaml(self, indent=2):
        """
        Convert Genome object to yaml.

        Returns
        -------
        genome_yaml : yaml
            YAML representation of genome.
        """

        genome_yaml = (
            yaml.dump(self.to_dict(), sort_keys=False, indent=indent)
            .replace("null", "")
            .replace("''", "")
            + "\n"
        )
        return genome_yaml

    def lineage_assignment(
        self,
        barcode_summary,
        barcodes,
        tree,
        recombinant_lineages,
        recombinant_tree,
        lineage_to_clade,
    ):
        """
        Assign genome to a lineage based on the top barcode matches.

        Parameters
        ----------
        barcode_summary : pd.core.frame.DataFrame
            Dataframe of barcode counts from Barcode.search().
        tree : Bio.Phylo.Tree
            Phylogenetic tree of lineage nomenclature.
        recombinant_lineages: list
            List of recombinant lineages.
        recombinant_tree: Bio.Phylo.Tree
            Phylogenetic tree of the 'X' clade (recombinant MRCA)

        Returns
        -------
        barcode : Barcode
            Summary of barcode detections, supports, and conflicts.
        """

        barcode = Barcode(
            genome=self,
            barcode_summary=barcode_summary,
            barcodes=barcodes,
            tree=tree,
            recombinant_lineages=recombinant_lineages,
            recombinant_tree=recombinant_tree,
            lineage_to_clade=lineage_to_clade,
        )

        if self.debug:
            lineage_str = barcode.to_yaml().replace("\n", "\n" + "\t" * 6)
            self.logger.info(str(datetime.now()) + "\t\t\t" + lineage_str)

        return barcode

    def parent_assignment(
        self,
        barcodes,
        tree,
        recombinant_lineages,
        recombinant_tree,
        lineage_to_clade,
        max_depth=1,
        min_subs=1,
        min_length=1,
        min_consecutive=1,
        edge_cases=False,
    ):
        """
        Assign genome to a lineage based on the top barcode matches.

        Parameters
        ----------
        barcodes : pd.core.frame.DataFrame
            Dataframe of lineage barcodes, from `rebar barcodes` subcommand
        tree : Bio.Phylo.Tree
            Phylogenetic tree of lineage nomenclature.
        recombinant_lineages: list
            List of recombinant lineages.
        recombinant_tree: Bio.Phylo.Tree
            Phylogenetic tree of the 'X' clade (recombinant MRCA)
        max_depth : int
            Maximum search depth of the top lineages.
        min_subs : int
            Minimum number of consecutive barcode subs contributed by a parent.
        min_consecutive : int
            Minimum number of consecutive bases contributed by a parent.

        Attributes Modified
        -------
        self.recombination : Recombination
        """

        # Skip clear non-recombinants
        if self.lineage.recombinant == False:
            return 0

        # Save a copy of the barcode summary, before we modify it
        barcode_summary = copy(self.barcode_summary)
        # Keep a list to exclude from parent search
        exclude_lineages = []

        # Option 1. Definitely a recursive recombinant.
        #           Exclude recombinant lineages that are not the known parent
        if self.lineage.recursive:
            exclude_lineages += self.lineage.top_lineages
            lineage_path = recombinant_tree.get_path(self.lineage.recombinant)
            lineage_parent = lineage_path[-2].name
            exclude_lineages += [l for l in recombinant_lineages if l != lineage_parent]
        # Option 2. Definitely not a recursive recombinant.
        #           Exclude all recombinant lineages from new search.
        #           Ex. XBB.1.5 is not a recursive recombinant (BA.2.10* and BA.2.75*)
        #           If we remove all recombinant lineages from it's barcode summary
        #           the top lineage will become BJ.1.1 (BA.2.10*)
        elif not self.lineage.recursive:
            exclude_lineages += recombinant_lineages
        # Option 3. Potentially recursive recombinant
        #           Exclude only original backbone lineages from new search.
        #           Ex. XBL is a recursive recombinant (XBB.1* and BA.2.75*)
        else:
            exclude_lineages += self.lineage.top_lineages

        barcode_summary = barcode_summary[
            ~barcode_summary["lineage"].isin(exclude_lineages)
        ].sort_values(by=["total", "lineage"])

        # ---------------------------------------------------------------------
        # EDGE CASES
        # This section is for legacy detection of SARS-CoV-2 lineages

        if edge_cases:
            # `handle_edge_cases` will adjust these global parameters, just
            #   for this genome if it's an edge case.
            result = handle_edge_cases(
                self, barcode_summary, tree, min_subs, min_length, min_consecutive
            )
            min_subs = result["min_subs"]
            min_length = result["min_length"]
            min_consecutive = result["min_consecutive"]
            barcode_summary = result["barcode_summary"]

        # ---------------------------------------------------------------------
        # Assign parent_1

        if self.debug:
            self.logger.info(str(datetime.now()) + "\t\tPARENT 1:")

        self.recombination.parent_1 = self.lineage_assignment(
            barcode_summary=barcode_summary,
            barcodes=barcodes,
            tree=tree,
            recombinant_lineages=recombinant_lineages,
            recombinant_tree=recombinant_tree,
            lineage_to_clade=lineage_to_clade,
        )

        # If parent_1 has no conflict_refs, don't search for more parents
        # i.e. it's a perfect match, no evidence of recombination
        if self.debug and len(self.recombination.parent_1.conflict_ref) == 0:
            self.logger.info(
                str(datetime.now())
                + "\t\t"
                + self.recombination.parent_1.lineage
                + " is a perfect match, halting recombinant search."
            )
            self.lineage.recombinant = False
            return 0

        # ---------------------------------------------------------------------
        # Assign parent_2

        # First, exclude all descendants of parent_1 from the search
        parent_1_tree = next(tree.find_clades(self.recombination.parent_1.lineage))
        parent_1_descendants = [c.name for c in parent_1_tree.find_clades()]
        exclude_lineages += parent_1_descendants

        # Next, restrict barcodes to only lineages with the
        # conflict_alt (subs that are not in parent_1's barcode)
        # keep lineages that have ANY number of these substitutions, which means
        # the final retained lineages will be very permissive/sensitive.
        conflict_alt_summary = self.summarise_barcodes(
            barcodes=barcodes, barcodes_subs=self.recombination.parent_1.conflict_alt
        )
        # if self.debug:
        #     df_md = conflict_alt_summary.to_markdown(index=False).replace(
        #         "\n", "\n" + "\t" * 7
        #     )
        #     self.logger.info(
        #         str(datetime.now())
        #         + "\t\t\tCONFLICT ALT (INCLUDE):\n"
        #         + ("\t") * 7 + df_md
        #     )

        # Remove lineages with the conflict_ref (ref bases
        # where parent_1 has a mutation)
        conflict_ref_summary = self.summarise_barcodes(
            barcodes=barcodes, barcodes_subs=self.recombination.parent_1.conflict_ref
        )
        # exclude lineages that have ALL ref bases, which means the final
        # retained lineages are very permissive/sensitive.
        conflict_ref_summary = conflict_ref_summary[
            conflict_ref_summary["total"]
            == len(self.recombination.parent_1.conflict_alt)
        ]
        exclude_lineages += list(conflict_ref_summary["lineage"])
        # if self.debug:
        #     df_md = conflict_ref_summary.to_markdown(index=False).replace(
        #         "\n", "\n" + "\t" * 7
        #     )
        #     self.logger.info(
        #         str(datetime.now())
        #         + "\t\t\tCONFLICT REF (EXCLUDE):\n"
        #         + ("\t") * 7 + df_md
        #     )

        # The new barcode_summary is just lineages that will help
        # us resolve these conflicts
        barcode_summary = conflict_alt_summary[
            ~conflict_alt_summary["lineage"].isin(exclude_lineages)
        ]

        # Now, we search through the barcodes
        recombination_detected = False
        depth = 0

        # Search through the top lineages for a suitable parent 2
        # Keep searching unless we max out the depth counter or find recombination
        while depth < max_depth and not recombination_detected:
            depth += 1

            # Exclude the previous loops lineages
            barcode_summary = barcode_summary[
                ~barcode_summary["lineage"].isin(exclude_lineages)
            ]

            # If we've run out of barcodes, no recombination!
            if len(barcode_summary) == 0:
                if self.debug:
                    self.logger.info(
                        str(datetime.now()) + "\t\tNo more barcodes to parse."
                    )
                break

            # Summarize the barcode support for the next top lineages
            if self.debug:
                self.logger.info(
                    str(datetime.now())
                    + "\t\tPARENT 2 | DEPTH: {} / {}".format(depth, max_depth)
                )

            # Summarize the barcode support for the next top lineages
            parent_2 = self.lineage_assignment(
                barcode_summary=barcode_summary,
                barcodes=barcodes,
                tree=tree,
                recombinant_lineages=recombinant_lineages,
                recombinant_tree=recombinant_tree,
                lineage_to_clade=lineage_to_clade,
            )

            # Detect recombination
            recombination = Recombination(
                genome=self,
                parent_1=self.recombination.parent_1,
                parent_2=parent_2,
                min_subs=min_subs,
                min_length=min_length,
                min_consecutive=min_consecutive,
            )
            recombination.depth = depth

            # If recombination was detected, break free of search loop!
            if len(recombination.breakpoints) > 0:
                recombination_detected = True
                self.recombination = recombination

            # Otherwise, update our exclude lineages for the next search
            else:
                # exclude_lineages += parent_2.top_lineages
                exclude_lineages += [
                    l
                    for l in parent_2.top_lineages
                    if l not in parent_2.outlier_lineages
                ]

        if not recombination_detected and not self.lineage.recombinant:
            self.lineage.recombinant = False

        return 0

    def validate_recombination(self, tree, recombinant_lineages):
        # Identify which lineages are known recombinants
        # ie. descended from the "X" recombinant MRCA node
        lineages = [c.name for c in tree.find_clades()]
        status = None
        warn = False

        if self.id in lineages:
            if self.id in recombinant_lineages:
                expected = "positive"
            else:
                expected = "negative"
            # Correct positive
            if self.lineage.recombinant and expected == "positive":
                status = "positive"
            # Correct negative
            elif not self.lineage.recombinant and expected == "negative":
                status = "negative"
            # False positive
            elif self.lineage.recombinant and expected == "negative":
                status = "false_positive"
                warn = True
            # False negative
            elif not self.lineage.recombinant and expected == "positive":
                status = "false_negative"
                warn = True

            msg = (
                "Validation fail for {}".format(self.id)
                + ", expected='{}'".format(expected)
                + ", actual='{}'".format(status)
            )
            if warn:
                self.logger.info(str(datetime.now()) + "\t\tWARNING: " + msg)
                # # Full error raise
                # # has multiprocess hang complications
                # if self.validate_fail:
                #     raise SystemExit(RebarError("RebarError: " + msg))
                # Just a warning
                # else:
                #    self.logger.info(str(datetime.now()) + "\t\tWARNING: " + msg)

        return status


def genome_mp(iterator, **kwargs):
    """
    Create Genome with multiprocess.
    Used to control the named parameters that are passed.
    """

    # When creating from FASTA, iterator is a SeqRecord
    if type(iterator) == SeqRecord:
        kwargs["record"] = iterator

    # When creating from SUBS df, iterator is a tuple
    elif type(iterator) == tuple:
        # First value is index, second value is series
        kwargs["subs_row"] = iterator[1]

    # else:
    #    raise RebarError("Unknown iterator was passed to genome_mp.")

    genome = Genome(**kwargs)
    return genome
