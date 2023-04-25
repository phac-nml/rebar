#!/usr/bin/env/python3

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import patches, colors
import os
import logging

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

    # The maximum y coordinate will depend on the number of
    # records/strains in barcodes
    y_coord = len(barcode_df.columns[1:])
    # Add an extra half coordinate to add y_space between parents and samples
    y_coord += 0.5

    # What is going to be the width and height of each sub box?
    x_inc = genome_length / (len(barcode_df) + 1)
    y_inc = (x_inc * 0.8 * y_coord) / genome_length

    # -----------------------------------------------------------------------------
    # PLOT SETUP

    # Width
    if len(barcode_df) < 10:
        width = 10
    else:
        width = 0.25 * len(barcode_df)

    # Height
    if y_coord < 5:
        height = 5
    else:
        height = y_inc * 3 * y_coord + y_inc * 2

    fig, ax = plt.subplots(1, 1, figsize=(width, height), dpi=300)

    # -----------------------------------------------------------------------------
    # PLOT

    x_coord = 0

    for rec_i, rec in enumerate(barcode_df.iterrows()):

        # Identify base origins
        ref_base = rec[1]["Reference"]
        parent_1_base = rec[1][parent_1]
        parent_2_base = rec[1][parent_2]

        # This will be the true x position where the sub should be plotted
        x_coord += x_inc
        # This is the y position, increment as we process samples
        y_break = len(barcode_df.columns) - 4
        y_coord = 0

        for i, label in enumerate(reversed(barcode_df.columns[1:])):
            if i == y_break:
                y_coord += y_inc * 1.5
            else:
                y_coord += y_inc
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
                    size=9,
                    ha="right",
                    va="center",
                    fontweight="bold",
                )

            # Draw sub text bases
            ax.text(x_coord, y_coord, base, size=7, ha="center", va="center")

    # Breakpoints
    for i, breakpoint in enumerate(breakpoints_split):
        # -1 because sub in region
        start = int(breakpoint.split(":")[0]) - 1
        x_ind = barcode_df[barcode_df["coord"] == start].index.values[0]
        x_coord = (x_inc) * (x_ind + 1.5)
        y_min = y_inc / 2
        y_max = y_inc * len(barcode_df.columns)

        ax.plot(
            [x_coord, x_coord],
            [y_min, y_max],
            ls="--",
            lw=0.5,
            color="black",
            clip_on=False,
        )
        ax.text(
            x_coord,
            y_max,
            "Breakpoint #" + str(i + 1),
            size=7,
            ha="center",
            va="bottom",
        )

    # -----------------------------------------------------------------------------
    # PLOT WRAPUP

    # Hide the axes lines
    for spine in ax.spines:
        ax.spines[spine].set_visible(False)

    # Hide the ticks on all axes, we've manually done it
    xticks = [i + (x_inc * i) for i in range(1, (len(barcode_df) + 1))]
    xtick_labels = [str(c) for c in list(barcode_df["coord"])]

    plt.xticks(xticks)
    ax.set_xticklabels(xtick_labels, rotation=90, fontsize=9)
    ax.set_yticks([])

    # Set the X axis limits to the genome length
    ax.set_xlim(0, genome_length)
    # ax.tick_params(axis='x', labelsize=8)
    ax.set_xlabel("Genomic Coordinate", fontsize=11, fontweight="bold")

    # Export
    plt.tight_layout()
    plt.savefig(output)
