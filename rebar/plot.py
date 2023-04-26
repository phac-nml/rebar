#!/usr/bin/env/python3

import os
import logging

import matplotlib.pyplot as plt
from matplotlib import patches, colors
import pandas as pd

logging.getLogger("matplotlib.font_manager").disabled = True


def plot(barcodes, summary, output):

    # Create output directory if it doesn't exist
    outdir = os.path.dirname(output)
    if not os.path.exists(outdir) and outdir != ".":
        os.makedirs(outdir)

    # Import dataframes
    barcode_df = pd.read_csv(barcodes, sep="\t")
    summary_df = pd.read_csv(summary, sep="\t")

    # Parse data from summary
    breakpoints = summary_df["breakpoints"].values[0]
    breakpoints_split = breakpoints.split(",")
    genome_length = int(summary_df["genome_length"].values[0])

    # Identify the first and second parent
    parent_1 = barcode_df.columns[2]
    parent_2 = barcode_df.columns[3]

    # Set up palette to distinguish parents and mutation source
    cmap = plt.get_cmap("tab20", 20)
    palette = [colors.rgb2hex(cmap(i)) for i in range(cmap.N)]
    parent_1_mut_color = palette[0]
    parent_1_ref_color = palette[1]
    parent_2_mut_color = palette[2]
    parent_2_ref_color = palette[3]
    ref_color = "#c4c4c4"

    num_records = len(barcode_df.columns[1:])
    num_subs = len(barcode_df)

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
    # PLOT

    # This will be x position in units of genomic coordinates, increment
    x_coord = 0

    for rec_i, rec in enumerate(barcode_df.iterrows()):

        # Identify base origins
        ref_base = rec[1]["Reference"]
        parent_1_base = rec[1][parent_1]
        parent_2_base = rec[1][parent_2]

        x_coord += x_inc
        # This will be y position in units of genomic coordinates, increment
        y_coord = (y_inc) * 0.5

        for i, label in enumerate(reversed(barcode_df.columns[1:])):
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

    # Breakpoints
    y_min = y_inc * 0.1
    y_max = (y_inc * num_records) + (y_inc)
    for i, breakpoint in enumerate(breakpoints_split):
        # get region start and end

        start = int(breakpoint.split(":")[0]) - 1
        end = int(breakpoint.split(":")[1]) + 1

        start_match = barcode_df[barcode_df["coord"] == start]
        end_match = barcode_df[barcode_df["coord"] == end]

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
    xticks = [(x_inc * i) for i in range(1, (len(barcode_df) + 1))]
    xtick_labels = [str(c) for c in list(barcode_df["coord"])]

    plt.xticks(xticks)
    ax.set_xticklabels(xtick_labels, rotation=90, fontsize=fontsize)
    # ax.set_xticks([])
    ax.set_yticks([])

    # # Set the Y axis limits to the genome length
    ax.set_ylim(0, (num_records * y_inc) + y_break_gap)
    ax.set_xlim(0, genome_length + x_inc)
    ax.set_xlabel("Genomic Coordinate", fontsize=fontsize, fontweight="bold")

    ax.axes.set_aspect("equal")
    # Export
    plt.tight_layout()
    plt.savefig(output)


plot(
    barcodes="output/debug/barcodes/XBJ_BA.2.3.20_BA.5.2.tsv",
    summary="output/debug/summary.tsv",
    output="output/debug/plots/XBJ_BA.2.3.20_BA.5.2.png",
)
plot(
    barcodes="output/debug/barcodes/XBJ_BA.2.3.20_CK.2.1.1.tsv",
    summary="output/debug/summary.tsv",
    output="output/debug/plots/XBJ_BA.2.3.20_CK.2.1.1.png",
)
plot(
    barcodes="output/XBC/barcodes/XBC_B.1.617.2_CJ.1.tsv",
    summary="output/XBC/summary.tsv",
    output="output/XBC/plots/XBC_B.1.617.2_CJ.1.png",
)
