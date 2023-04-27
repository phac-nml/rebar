#!/usr/bin/env/python3

import os
import logging

import matplotlib.pyplot as plt
from matplotlib import patches, colors
import pandas as pd

logging.getLogger("matplotlib.font_manager").disabled = True


def plot(barcodes_df, summary_df, output):

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

    if num_subs < 10:
        width = 10
    else:
        width = num_subs * 0.25

    if num_records < 10:
        height = 5
    else:
        height = num_records

    fig, ax = plt.subplots(1, 1, dpi=200, figsize=(width, height))

    # -----------------------------------------------------------------------------
    # PLOT GENOME GUIDE

    guide_x = x_inc * 0.5
    guide_y = (y_inc * num_records) + y_break + (y_inc * 2)
    rect = patches.Rectangle(
        (guide_x, guide_y),
        genome_length,
        y_inc,
        alpha=0.2,
        fill=True,
        edgecolor="none",
        facecolor="dimgrey",
    )
    ax.add_patch(rect)

    for coord in list(barcodes_df["coord"]):
        ax.plot(
            [coord, coord],
            [guide_y, guide_y + y_inc],
            lw=linewidth,
            color="black",
            clip_on=False,
        )
        break

    # -----------------------------------------------------------------------------
    # PLOT SAMPLES

    # This will be x position in units of genomic coordinates, increment
    x_coord = 0

    # Iterate through subs, which are columns in plot
    for rec_i, rec in enumerate(barcodes_df.iterrows()):

        # Identify base origins
        ref_base = rec[1]["Reference"]
        parent_1_base = rec[1][parent_1]
        parent_2_base = rec[1][parent_2]

        x_coord += x_inc
        # This will be y position in units of genomic coordinates, increment
        y_coord = (y_inc) * 0.5

        # Iterate through samples, which are rows in plot
        for i, label in enumerate(reversed(records)):
            base = rec[1][label]

            # Draw box
            box_w = x_inc * 0.8
            box_h = y_inc * 0.8
            box_x = x_coord - (box_w / 2)
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

    # -----------------------------------------------------------------------------
    # PLOT BREAKPOINTS

    y_min = y_inc * 0.1
    y_max = (y_inc * num_records) + (y_inc * 0.5)
    for i, breakpoint in enumerate(breakpoints_split):
        # get region start and end

        start = int(breakpoint.split(":")[0]) - 1
        end = int(breakpoint.split(":")[1]) + 1

        start_match = barcodes_df[barcodes_df["coord"] == start]
        end_match = barcodes_df[barcodes_df["coord"] == end]

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
            box_h = (y_inc) * num_records + (y_inc * 0.4)

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
    ax.set_ylim(0, (num_records * y_inc) + y_break_gap + (y_inc * 3))
    ax.set_xlim(0, genome_length + x_inc)
    ax.set_xlabel("Genomic Coordinate", fontsize=fontsize, fontweight="bold")

    ax.axes.set_aspect("equal")
    # Export
    # plt.tight_layout()
    plt.savefig(output)


barcodes_df = pd.read_csv("output/XBC/barcodes/XBC_CJ.1_B.1.617.2.tsv", sep="\t")
summary_df = pd.read_csv("output/XBC/summary.tsv", sep="\t")
plot(barcodes_df=barcodes_df, summary_df=summary_df, output="test.png")
