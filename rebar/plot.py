#!/usr/bin/env/python3

import os
import logging
from itertools import cycle

import matplotlib.pyplot as plt
from matplotlib import patches, colors

logging.getLogger("matplotlib.font_manager").disabled = True


def plot(barcodes_df, summary_df, annot_df, output):

    # Debugging
    # barcodes_df = barcodes_df.head(35)

    # Create output directory if it doesn't exist
    outdir = os.path.dirname(output)
    if not os.path.exists(outdir) and outdir != "." and outdir != "":
        os.makedirs(outdir)

    genome_length = int(summary_df["genome_length"].values[0])

    # Identify the first and second parent lineage
    parent_1 = barcodes_df.columns[2]
    parent_2 = barcodes_df.columns[3]

    # Identify clade of each parent
    parents = summary_df["parents_clade_lineage"].values[0]
    parent_1_clade_lineage = parents.split(",")[0]
    parent_2_clade_lineage = parents.split(",")[1]

    # Set up palette to distinguish parents and mutation source
    cmap = plt.get_cmap("tab20", 20)
    palette = [colors.rgb2hex(cmap(i)) for i in range(cmap.N)]
    parent_1_mut_color = palette[0]  # Blue Dark
    parent_1_ref_color = palette[1]  # Blue Light
    parent_2_mut_color = palette[2]  # Orange Dark
    parent_2_ref_color = palette[3]  # Orange Light

    ref_color = "#c4c4c4"

    records = barcodes_df.columns[1:]
    num_records = len(records)
    num_subs = len(barcodes_df)

    # What is going to be the width and height of each sub box?
    x_inc = genome_length / num_subs
    # Make y size the same, so it produces squares
    y_inc = x_inc

    # This is the y position position to add extra space
    y_break = num_records - 4
    y_break_gap = y_inc * 0.5

    # -------------------------------------------------------------------------
    # Plot Setup

    # # Fontsize of 9 looks good for 20 subs
    fontsize = 9
    linewidth = 0.5
    guide_tick_intervals = 5000

    target_height = 20
    target_width = 20

    # True matplotlib maximum
    max_height = 50
    max_width = 50

    # Default width and height
    width = 10
    height = 10

    if num_subs > 10:

        width = num_subs
        if num_subs > target_width:
            width = num_subs * 0.25

        if width > max_width:
            width = max_width

    if num_records > 5:
        height = num_records
        if num_records > target_height:
            height = max_height
        if height > max_height:
            height = max_height

    if width < 10:
        fontsize = 7
    else:
        guide_tick_intervals = 2500

    fig, ax = plt.subplots(1, 1, dpi=200, figsize=(width, height))

    # Current height of the section for plotting
    current_y = 0

    # -----------------------------------------------------------------------------
    # PLOT SAMPLES AND SNPS

    # This will be x position in units of genomic coordinates, increment
    x_coord = x_inc

    # Sub box sizes
    box_w = x_inc * 0.8
    box_h = y_inc * 0.8

    # Iterate through subs, which are columns in plot
    for rec_i, rec in enumerate(barcodes_df.iterrows()):

        # Identify base origins
        genome_coord = rec[1]["coord"]
        ref_base = rec[1]["Reference"]
        parent_1_base = rec[1][parent_1]
        parent_2_base = rec[1][parent_2]

        box_x = x_coord - (box_w / 2)

        # This will be y position in units of genomic coordinates, increment
        y_coord = (y_inc) * 0.5

        # Iterate through samples, which are rows in plot
        for i, label in enumerate(reversed(records)):
            base = rec[1][label]

            # Draw box
            box_y = y_coord - (box_h / 2)

            if label == "Reference":
                box_c = ref_color
            elif base == parent_1_base:
                if base == ref_base:
                    box_c = parent_1_ref_color
                else:
                    box_c = parent_1_mut_color
            elif base == parent_2_base:
                if base == ref_base:
                    box_c = parent_2_ref_color
                else:
                    box_c = parent_2_mut_color
            else:
                box_c = "white"

            rect = patches.Rectangle(
                (box_x, box_y),
                box_w,
                box_h,
                alpha=0.90,
                fill=True,
                edgecolor="none",
                facecolor=box_c,
            )
            ax.add_patch(rect)

            # On the first time parsing sub, write sample label
            if rec_i == 0:
                # If parent, also include clade
                if label == parent_1:
                    label = "{} ({})".format(parent_1, parent_1_clade_lineage)
                elif label == parent_2:
                    label = "{} ({})".format(parent_2, parent_2_clade_lineage)
                ax.text(
                    -50,
                    y_coord,
                    label,
                    size=fontsize,
                    ha="right",
                    va="center",
                )

            # Draw sub text bases
            ax.text(x_coord, y_coord, base, size=fontsize, ha="center", va="center")

            # Add an extra gap after recombinant parents
            if i == y_break:
                y_coord += y_inc + (y_break_gap)
            elif i == (num_records) - 1:
                y_coord += y_inc / 2
            else:
                y_coord += y_inc

        x_coord += x_inc

    top_sample_y = y_coord
    current_y = y_coord

    # -----------------------------------------------------------------------------
    # Breakpoint Labels

    # Parse breakpoints from summary
    breakpoints = list(summary_df["breakpoints"])
    # Check if there are multiple different breakpoint versions
    # Use whichever one is most common
    if len(set(breakpoints)) > 1:
        max_count = 0
        max_breakpoints = None
        for bp in breakpoints:
            count = breakpoints.count(bp)
            if count > max_count:
                max_count = count
                max_breakpoints = bp

        breakpoints = max_breakpoints
    else:
        breakpoints = breakpoints[0]

    breakpoints_split = breakpoints.split(",")

    # Align bottom of breakpoints dividor with sub box bottom
    bp_y1 = y_inc * 0.1
    # This is the y coordinate for the text, bottom aligned
    bp_y2 = current_y + y_break + (y_inc * 1.5)

    # Write breakpoints label of left
    ax.text(
        -50,
        bp_y2,
        "Breakpoints",
        size=fontsize,
        ha="right",
        va="center",
    )

    for i, breakpoint in enumerate(breakpoints_split):
        # get region start and end

        start = int(breakpoint.split(":")[0]) - 1
        end = int(breakpoint.split(":")[1]) + 1

        start_match = barcodes_df[barcodes_df["coord"] == start]
        end_match = barcodes_df[barcodes_df["coord"] == end]

        if len(start_match) == 0:
            continue

        start_i = start_match.index.values[0]
        # End coordinates are not inclusive, move back one
        end_i = end_match.index.values[0] - 1

        # If start and end are adjacent coordinates, plot dashed line
        if (end_i - start_i) == 0:
            x_coord = (x_inc) * (start_i + 1.5)

            ax.plot(
                [x_coord, x_coord],
                [bp_y1, bp_y2],
                ls="--",
                lw=linewidth,
                color="black",
                clip_on=False,
            )
        # Plot greyed, out dashed rectangle
        else:

            box_x = (x_inc) * (start_i + 1.5)
            box_x2 = (x_inc) * (end_i + 1.5)
            box_y = bp_y1
            box_w = box_x2 - box_x
            box_h = (y_inc) * num_records + (y_inc * 0.8)

            # Greyed out rectangle
            rect = patches.Rectangle(
                (box_x, box_y),
                box_w,
                box_h,
                alpha=0.4,
                fill=True,
                facecolor="black",
                edgecolor="none",
                lw=linewidth,
                ls="--",
            )
            ax.add_patch(rect)

            # Dashed rectangle outline
            rect = patches.Rectangle(
                (box_x, box_y),
                box_w,
                box_h,
                alpha=1,
                fill=False,
                facecolor="none",
                edgecolor="black",
                lw=linewidth,
                ls="--",
            )
            ax.add_patch(rect)

            # Bonus little dashed tick in center
            line_x = box_x + (box_w / 2)
            line_y1 = box_y + box_h
            line_y2 = bp_y2

            ax.plot(
                [line_x, line_x],
                [line_y1, line_y2],
                ls="--",
                lw=linewidth,
                color="black",
                clip_on=False,
            )

            x_coord = box_x + (box_x2 - box_x) / 2

        ax.text(
            x_coord,
            bp_y2,
            "Breakpoint #" + str(i + 1),
            size=fontsize,
            ha="center",
            va="center",
            bbox=dict(
                facecolor="white", edgecolor="black", lw=linewidth, boxstyle="round"
            ),
        )

    current_y = bp_y2

    # -----------------------------------------------------------------------------
    # Genome Coordinates Guide

    guide_x = x_inc * 0.5
    guide_y1 = current_y + y_break + (y_inc * 2)
    guide_y2 = guide_y1 + y_inc
    guide_h = guide_y2 - guide_y1

    # Write breakpoints label of left
    ax.text(
        -50,
        guide_y1 + (guide_y2 - guide_y1) / 2,
        "Guide",
        size=fontsize,
        ha="right",
        va="center",
    )

    # Write grey square as background
    rect = patches.Rectangle(
        (guide_x, guide_y1),
        genome_length,
        guide_h,
        alpha=0.2,
        fill=True,
        edgecolor="none",
        facecolor="dimgrey",
    )
    ax.add_patch(rect)

    # Add ticks
    for x in list(range(0, genome_length, guide_tick_intervals)) + [genome_length]:
        tick_text = str(x)
        x += guide_x
        tick_y2 = guide_y1
        tick_y1 = tick_y2 - (y_inc * 0.25)
        text_y = tick_y1 - (y_inc * 0.25)

        ax.plot(
            [x, x],
            [tick_y1, tick_y2],
            lw=linewidth,
            color="black",
        )
        ax.text(
            x,
            text_y,
            str(tick_text),
            size=fontsize - 2,
            ha="center",
            va="top",
            color="black",
        )

    # Iterate through subs to show on guide
    guide_sub_x = x_inc
    # Sub box sizes
    box_w = x_inc * 0.8

    for rec in barcodes_df.iterrows():

        # Identify base origins
        genome_coord = rec[1]["coord"] + (x_inc * 0.5)
        ref_base = rec[1]["Reference"]
        parent_1_base = rec[1][parent_1]
        parent_2_base = rec[1][parent_2]

        box_x = guide_sub_x - (box_w / 2)

        guide_color = "black"

        # Plot a line on the guide
        ax.plot(
            [genome_coord, genome_coord],
            [guide_y1, guide_y2],
            lw=linewidth,
            color=guide_color,
            clip_on=False,
        )

        # Draw a polygon from guide to top row of subs
        # P1: Top Guide, P2: Top left of SNP box, P3: Top right of SNP box
        x1 = box_x
        x2 = box_x + box_w
        x3 = genome_coord
        y1 = top_sample_y - (y_inc * 0.1)
        y2 = guide_y1

        poly_coords = [[x1, y1], [x2, y1], [x3, y2]]
        poly = patches.Polygon(
            poly_coords,
            closed=True,
            alpha=0.2,
            fill=True,
            edgecolor="none",
            facecolor="dimgrey",
        )
        ax.add_patch(poly)

        guide_sub_x += x_inc

    current_y = guide_y2

    # -----------------------------------------------------------------------------
    # Parental Regions

    regions_x = x_inc * 0.5
    regions_y1 = current_y + y_break + (y_inc * 1)
    regions_y2 = regions_y1 + y_inc
    regions_h = regions_y2 - regions_y1

    # Write label to left
    ax.text(
        -50,
        regions_y1 + (regions_y2 - regions_y1) / 2,
        "Regions",
        size=fontsize,
        ha="right",
        va="center",
    )

    # Parse regions from summary
    # Use whichever region matches the previously identified max breakpoints
    regions = summary_df[summary_df["breakpoints"] == breakpoints]["regions"].values[0]

    regions_split = regions.split(",")
    prev_end = None

    # Example: 3796-22599|B.1.617.2
    for region in regions_split:
        coords = region.split("|")[0]
        parent = region.split("|")[1]
        # First region, plot from beginning
        if region == regions_split[0]:
            start = 0
        else:
            start = int(coords.split("-")[0])
        # Last region, plot to end
        if region == regions_split[-1]:
            end = genome_length
        else:
            end = int(coords.split("-")[1])

        # Adjust x coord
        start += regions_x
        end += regions_x

        if parent == parent_1:
            color = parent_1_mut_color
        else:
            color = parent_2_mut_color

        # Parental region box
        rect = patches.Rectangle(
            (start, regions_y1),
            end - start,
            regions_h,
            fill=True,
            edgecolor="black",
            lw=linewidth,
            facecolor=color,
        )
        ax.add_patch(rect)

        text_x = start + (end - start) / 2
        # Offset height of every other text label
        text_y = regions_y1 + (regions_h * 1.5)
        ax.plot([text_x, text_x], [regions_y2, text_y], lw=linewidth, c="black")
        ax.text(text_x, text_y, parent, size=fontsize, ha="center", va="bottom")

        # Breakpoints box
        if prev_end:
            rect = patches.Rectangle(
                (prev_end, regions_y1),
                start - prev_end,
                regions_h,
                fill=True,
                edgecolor="black",
                lw=linewidth,
                facecolor=ref_color,
            )
            ax.add_patch(rect)

        prev_end = end

    current_y = regions_y2

    # -----------------------------------------------------------------------------
    # Genome Annotations
    # Above genome guide

    annot_x = x_inc * 0.5
    annot_y1 = current_y + (y_inc * 1.75)
    annot_y2 = annot_y1 + y_inc
    annot_h = annot_y2 - annot_y1
    # Exclude the sub colors from the palette
    annot_palette = cycle(palette[4:])

    # Write label to left
    ax.text(
        -50,
        annot_y1 + (annot_y2 - annot_y1) / 2,
        "Genes",
        size=fontsize,
        ha="right",
        va="center",
    )

    for i, rec in enumerate(annot_df.iterrows()):
        start = rec[1]["start"]
        # End coordinates non-inclusive
        end = rec[1]["end"] - 1

        annot_x1 = rec[1]["start"] + annot_x
        annot_x2 = end = (rec[1]["end"] - 1) + annot_x
        annot_w = annot_x2 - annot_x1

        annot_c = next(annot_palette)

        # Plot annotation rectangle
        rect = patches.Rectangle(
            (annot_x1, annot_y1),
            annot_w,
            annot_h,
            alpha=1,
            fill=True,
            edgecolor="black",
            linewidth=linewidth,
            facecolor=annot_c,
        )
        ax.add_patch(rect)

        gene = rec[1]["abbreviation"]

        text_x = annot_x1 + (annot_w * 0.5)
        # Offset height of every other text label
        if i % 2 == 0:
            text_y = annot_y1 + (annot_h * 1.5)
        else:
            text_y = annot_y1 + (annot_h * 2.25)
        ax.plot(
            [text_x, text_x],
            [annot_y2, text_y],
            ls="-",
            lw=linewidth,
            color="black",
            clip_on=False,
        )
        ax.text(
            text_x, text_y, gene, size=fontsize, ha="center", va="bottom", color="black"
        )

    # -----------------------------------------------------------------------------
    # PLOT WRAPUP

    # Hide the axes lines
    for spine in ax.spines:
        ax.spines[spine].set_visible(False)

    # Hide the ticks on all axes, we've manually done it
    xticks = [(x_inc * i) for i in range(1, (len(barcodes_df) + 1))]
    xtick_labels = [str(c) for c in list(barcodes_df["coord"])]

    plt.xticks(xticks)
    ax.set_xticklabels(xtick_labels, rotation=90, fontsize=fontsize)
    # ax.set_xticks([])
    ax.set_yticks([])

    # # Set the Y axis limits to the genome length
    # Extra y_inc for:
    #  1. Breakpoints, 2. Genome Guide,
    ax.set_ylim(0, (num_records * y_inc) + y_break_gap + (y_inc * 10))
    ax.set_xlim(0, genome_length + x_inc)
    ax.set_xlabel("Genomic Coordinate", fontsize=fontsize)

    ax.axes.set_aspect("equal")
    # Export
    # plt.tight_layout()
    plt.savefig(output)

    # Close for memory management
    plt.close()


# # Testing code
# import pandas as pd
# annot_df = pd.read_csv("dataset/sars-cov-2-latest/annotations.tsv", sep="\t")

# barcodes_df = pd.read_csv("output/XBC/barcodes/XBC_CJ.1_B.1.617.2.tsv", sep="\t")
# summary_df = pd.read_csv("output/XBC/summary.tsv", sep="\t")
# plot(
#     barcodes_df=barcodes_df,
#     summary_df=summary_df,
#     annot_df=annot_df,
#     output="output/XBC/plots/XBC_CJ.1_B.1.617.2.png",
# )

# barcodes_df = pd.read_csv("output/XBJ/barcodes/XBJ_BA.2.3.20_BA.5.2.tsv", sep="\t")
# summary_df = pd.read_csv("output/XBJ/summary.tsv", sep="\t")
# plot(
#     barcodes_df=barcodes_df,
#     summary_df=summary_df,
#     annot_df=annot_df,
#     output="output/XBJ/plots/XBJ_BA.2.3.20_BA.5.2.png",
# )

# barcodes_df = pd.read_csv("output/XBB.1.16/barcodes/XBB_BJ.1_CJ.1.tsv", sep="\t")
# summary_df = pd.read_csv("output/XBB.1.16/summary.tsv", sep="\t")
# plot(
#     barcodes_df=barcodes_df,
#     summary_df=summary_df,
#     annot_df=annot_df,
#     output="output/XBB.1.16/plots/XBB_BJ.1_CJ.1.png",
# )
