#!/usr/bin/env/python3

import os
import logging

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

    # Parse data from summary
    breakpoints = summary_df["breakpoints"].values[0]
    breakpoints_split = breakpoints.split(",")
    genome_length = int(summary_df["genome_length"].values[0])

    # Identify the first and second parent
    parent_1 = barcodes_df.columns[2]
    parent_2 = barcodes_df.columns[3]

    # Set up palette to distinguish parents and mutation source
    cmap = plt.get_cmap("tab20", 20)
    palette = [colors.rgb2hex(cmap(i)) for i in range(cmap.N)]
    parent_1_mut_color = palette[0]  # Blue Dark
    parent_1_ref_color = palette[1]  # Blue Light
    parent_2_mut_color = palette[6]  # Red Dark
    parent_2_ref_color = palette[7]  # Red Light
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
    # y_break_gap = 0
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

    if num_subs >= 10:

        width = num_subs
        if num_subs > target_width:
            width = num_subs * 0.25

        if width > max_width:
            width = max_width

    if num_records >= 10:
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

    # -----------------------------------------------------------------------------
    # Genome Coordinates Guide

    guide_x = x_inc * 0.5
    guide_y1 = (y_inc * num_records) + y_break + (y_inc * 3)
    guide_y2 = guide_y1 + y_inc
    guide_h = guide_y2 - guide_y1
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

    # -----------------------------------------------------------------------------
    # Genome Annotations

    annot_y1 = guide_y1 + (y_inc * 1.5)
    annot_y2 = annot_y1 + y_inc
    annot_h = annot_y2 - annot_y1

    for i, rec in enumerate(annot_df.iterrows()):
        start = rec[1]["start"]
        # End coordinates non-inclusive
        end = rec[1]["end"] - 1

        annot_x1 = rec[1]["start"] + guide_x
        annot_x2 = end = (rec[1]["end"] - 1) + guide_x
        annot_w = annot_x2 - annot_x1

        annot_c = palette[i]

        rect = patches.Rectangle(
            (annot_x1, annot_y1),
            annot_w,
            annot_h,
            alpha=1,
            fill=True,
            edgecolor="none",
            facecolor=annot_c,
        )
        ax.add_patch(rect)

        gene = rec[1]["abbreviation"]

        text_x = annot_x1 + (annot_w * 0.5)
        # text_y = annot_y1 + (annot_h * 0.5)
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
    # PLOT SAMPLES AND SNPS

    # This will be x position in units of genomic coordinates, increment
    x_coord = x_inc

    # Snp box sizes
    box_w = x_inc * 0.8
    box_h = y_inc * 0.8
    top_sample_y = (y_inc * num_records) + y_break + (y_inc * 0.4)

    # Iterate through subs, which are columns in plot
    for rec_i, rec in enumerate(barcodes_df.iterrows()):

        # Identify base origins
        genome_coord = rec[1]["coord"]
        ref_base = rec[1]["Reference"]
        parent_1_base = rec[1][parent_1]
        parent_2_base = rec[1][parent_2]

        box_x = x_coord - (box_w / 2)

        # Plot this SNP on the genome guide
        if parent_1_base != ref_base:
            guide_color = parent_1_mut_color
        elif parent_2_base != ref_base:
            guide_color = parent_2_mut_color
        else:
            guide_color = "black"

        # Plot a colored line on the guide
        ax.plot(
            [genome_coord, genome_coord],
            [guide_y1, guide_y2],
            lw=1,
            color=guide_color,
            clip_on=False,
        )

        # Draw a polygon from guide to top row of subs
        # P1: Top Guide, P2: Top left of SNP box, P3: Top right of SNP box
        x1 = box_x
        x2 = box_x + box_w
        x3 = genome_coord
        y1 = top_sample_y
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

            if rec_i == 0:
                ax.text(
                    -50,
                    y_coord,
                    label,
                    size=fontsize,
                    ha="right",
                    va="center",
                    fontweight="bold",
                )

            # Draw sub text bases
            ax.text(x_coord, y_coord, base, size=fontsize, ha="center", va="center")

            # Add an extra gap after recombinant parents
            if i == y_break:
                y_coord += y_inc + (y_break_gap)
            else:
                y_coord += y_inc

        x_coord += x_inc

    # -----------------------------------------------------------------------------
    # PLOT BREAKPOINTS

    y_min = y_inc * 0.1
    y_max = (y_inc * num_records) + (y_inc)
    for i, breakpoint in enumerate(breakpoints_split):
        # get region start and end

        start = int(breakpoint.split(":")[0]) - 1
        end = int(breakpoint.split(":")[1]) + 1

        start_match = barcodes_df[barcodes_df["coord"] == start]
        end_match = barcodes_df[barcodes_df["coord"] == end]

        if len(start_match) == 0:
            continue

        start_i = start_match.index.values[0]
        end_i = end_match.index.values[0]

        # If start and end are adjacent coordinates, plot dashed line
        if (end_i - start_i) == 1:
            x_coord = (x_inc) * (start_i + 1.5)

            ax.plot(
                [x_coord, x_coord],
                [y_min, y_max],
                ls="--",
                lw=linewidth,
                color="black",
                clip_on=False,
            )
        else:

            box_x = (x_inc) * (start_i + 1.5)
            box_x2 = (x_inc) * (end_i + 1.5)
            box_y = y_min
            box_w = box_x2 - box_x
            box_h = (y_inc) * num_records + (y_inc * 0.8)

            rect = patches.Rectangle(
                (box_x, box_y),
                box_w,
                box_h,
                alpha=1,
                fill=False,
                edgecolor="black",
                lw=linewidth,
                ls="--",
            )
            ax.add_patch(rect)

            x_coord = box_x + (box_x2 - box_x) / 2

        ax.text(
            x_coord,
            y_max,
            "Breakpoint #" + str(i + 1),
            size=fontsize,
            ha="center",
            va="bottom",
            bbox=dict(
                facecolor="white", edgecolor="black", lw=linewidth, boxstyle="round"
            ),
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
    ax.set_xlabel("Genomic Coordinate", fontsize=fontsize, fontweight="bold")

    ax.axes.set_aspect("equal")
    # Export
    # plt.tight_layout()
    plt.savefig(output)

    # Close for memory management
    plt.close()
