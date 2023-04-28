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
    # The first three are reference, parent_1, parent_2
    y_break = 3

    # -------------------------------------------------------------------------
    # Plot Setup

    # # Fontsize of 9 looks good for 20 subs
    fs = 9
    linewidth = 0.5
    plt.rcParams["hatch.linewidth"] = linewidth
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
        fs = 7
    else:
        guide_tick_intervals = 2500

    # Current height of the section for plotting
    current_y = 0

    # section_x = x_inc * 0.5
    section_x = 0
    section_label_x = -(x_inc * 0.5)
    # How weight/tall substitution boxes should be
    sub_box_w = x_inc * 0.8
    sub_box_h = sub_box_w

    fig, ax = plt.subplots(1, 1, dpi=200, figsize=(width, height))

    # -----------------------------------------------------------------------------
    # SECTION: REGIONS

    # General section dimensions
    section_y2 = current_y
    section_y1 = section_y2 - y_inc
    section_h = section_y2 - section_y1
    section_label_y = section_y1 + (section_y2 - section_y1) / 2
    section_label = "Parents"

    # Write section label
    ax.text(
        section_label_x, section_label_y, "Parents", size=fs, ha="right", va="center"
    )

    # Parse regions from summary
    regions = list(summary_df["regions"])
    # Check if there are multiple different breakpoint versions
    # Use whichever one is most common
    if len(set(regions)) > 1:
        max_count = 0
        max_regions = None
        for r in regions:
            count = regions.count(r)
            if count > max_count:
                max_count = count
                max_regions = r

        regions = max_regions
    else:
        regions = regions[0]

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
        start += section_x
        end += section_x

        if parent == parent_1:
            color = parent_1_mut_color
        else:
            color = parent_2_mut_color

        # Parental region box
        rect = patches.Rectangle(
            (start, section_y1),
            end - start,
            section_h,
            fill=True,
            edgecolor="black",
            lw=linewidth,
            facecolor=color,
        )
        ax.add_patch(rect)

        text_x = start + (end - start) / 2
        # Offset height of every other text label
        text_y = section_y1 + (section_h * 1.5)
        ax.plot([text_x, text_x], [section_y2, text_y], lw=linewidth, c="black")
        ax.text(text_x, text_y, parent, size=fs, ha="center", va="bottom")

        # Breakpoints box
        if prev_end:
            rect = patches.Rectangle(
                (prev_end, section_y1),
                start - prev_end,
                section_h,
                fill=True,
                edgecolor="black",
                lw=linewidth,
                facecolor=ref_color,
                hatch="////",
            )
            ax.add_patch(rect)

        prev_end = end

    current_y = section_y1

    # -----------------------------------------------------------------------------
    # SECTION: GENOME ANNOTATIONS AND SUB MARKERS

    # General section dimensions
    section_y2 = current_y - (y_inc * 2)
    section_y1 = section_y2 - y_inc
    section_h = section_y2 - section_y1
    section_label_y = section_y1 + (section_y2 - section_y1) / 2
    section_label = "Genome"

    # Write section label
    ax.text(
        section_label_x,
        section_label_y,
        section_label,
        size=fs,
        ha="right",
        va="center",
    )

    # Write grey rectangle as background
    rect = patches.Rectangle(
        (section_x, section_y1),
        genome_length,
        section_h,
        alpha=0.2,
        fill=True,
        edgecolor="none",
        facecolor="dimgrey",
    )
    ax.add_patch(rect)

    # Genome Annotations on top of guide rect

    # Exclude the sub colors from the palette
    annot_palette = cycle(palette[4:])

    for i, rec in enumerate(annot_df.iterrows()):

        start = rec[1]["start"]
        # End coordinates non-inclusive
        end = rec[1]["end"] - 1

        annot_x1 = rec[1]["start"] + section_x
        annot_x2 = end = (rec[1]["end"] - 1) + section_x
        annot_w = annot_x2 - annot_x1

        annot_c = next(annot_palette)

        # Plot annotation rectangle
        rect = patches.Rectangle(
            (annot_x1, section_y1),
            annot_w,
            section_h,
            alpha=1,
            fill=True,
            edgecolor="none",
            linewidth=linewidth,
            facecolor=annot_c,
        )
        ax.add_patch(rect)

        gene = rec[1]["abbreviation"]

        text_x = annot_x1 + (annot_w * 0.5)

        # Offset height of every other text label
        if i % 2 == 0:
            text_y = section_y1 + (section_h * 1.25)
        else:
            text_y = section_y1 + (section_h * 2)
        ax.plot(
            [text_x, text_x],
            [section_y2, text_y],
            ls="-",
            lw=linewidth,
            color="black",
            clip_on=False,
        )
        ax.text(text_x, text_y, gene, size=fs, ha="center", va="bottom", color="black")

    # Add subsitution ticks
    for x in list(range(0, genome_length, guide_tick_intervals)) + [genome_length]:
        tick_text = str(x)
        x += section_x
        tick_y2 = section_y1
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
            size=fs - 2,
            ha="center",
            va="top",
            color="black",
        )

    # Iterate through subs to show on guide
    for rec in barcodes_df.iterrows():

        genome_coord = rec[1]["coord"] + section_x

        # Identify base origins
        ref_base = rec[1]["Reference"]
        parent_1_base = rec[1][parent_1]
        parent_2_base = rec[1][parent_2]

        guide_color = "black"

        # Plot a line on the guide
        ax.plot(
            [genome_coord, genome_coord],
            [section_y1, section_y2],
            lw=linewidth,
            color=guide_color,
            clip_on=False,
        )

    current_y = section_y1

    # # -----------------------------------------------------------------------------
    # # SECTION: GUIDE TO SUB POLYGONS

    # section_y2 = current_y
    # # Add a little bit extra based on box_size
    # section_y1 = current_y - (y_inc * 3) - (x_inc - sub_box_h)

    # # Starting x point for sub boxes
    # sub_box_x = section_x - x_inc + (x_inc - sub_box_w) / 2

    # # Iterate through subs, which are columns in plot
    # for rec in barcodes_df.iterrows():

    #     sub_box_x += x_inc

    #     genome_coord = rec[1]["coord"] + section_x

    #     # Draw a polygon from guide to top row of subs
    #     # P1: Top Guide, P2: Top left of SNP box, P3: Top right of SNP box
    #     x1 = sub_box_x
    #     x2 = sub_box_x + sub_box_w
    #     x3 = genome_coord
    #     y1 = section_y1
    #     y2 = section_y2

    #     poly_coords = [[x1, y1], [x2, y1], [x3, y2]]
    #     poly = patches.Polygon(
    #         poly_coords,
    #         closed=True,
    #         alpha=0.2,
    #         fill=True,
    #         edgecolor="none",
    #         facecolor="dimgrey",
    #     )
    #     ax.add_patch(poly)

    # current_y = section_y1

    # # -----------------------------------------------------------------------------
    # # SECTION: BREAKPOINT LINES AND LABELS

    # # General section dimensions
    # section_y2 = section_y2
    # section_y1 = section_y1
    # section_label_y = section_y1 + (section_y2 - section_y1) / 2.25
    # section_label = "Breakpoints"

    # # Next section
    # next_section_y1 = section_y1 - (num_records * y_inc) - (y_inc)

    # # Write section label
    # ax.text(section_label_x, section_label_y, section_label, size=fs, ha="right")

    # # Parse breakpoints from summary based on previous regions
    breakpoints = summary_df[summary_df["regions"] == regions]["breakpoints"].values[0]
    # breakpoints_split = breakpoints.split(",")

    # # Align bottom of breakpoints dividor with next section sub box bottom
    # bp_y1 = next_section_y1 + (y_inc - sub_box_h) / 2
    # bp_y2 = section_label_y

    # for i, breakpoint in enumerate(breakpoints_split):

    #     # get region start and end
    #     start = int(breakpoint.split(":")[0]) - 1
    #     end = int(breakpoint.split(":")[1]) + 1

    #     start_match = barcodes_df[barcodes_df["coord"] == start]
    #     end_match = barcodes_df[barcodes_df["coord"] == end]

    #     if len(start_match) == 0:
    #         continue

    #     start_i = start_match.index.values[0]
    #     # End coordinates are not inclusive, move back one
    #     end_i = end_match.index.values[0] - 1

    #     # If start and end are adjacent coordinates, plot dashed line
    #     if (end_i - start_i) == 0:
    #         bp_x = section_x + (start_i * x_inc)

    #         ax.plot(
    #             [bp_x, bp_x],
    #             [bp_y1, bp_y2],
    #             ls="--",
    #             lw=linewidth,
    #             color="black",
    #             clip_on=False,
    #         )
    #     # # Plot greyed, out dashed rectangle
    #     # else:

    #     #     bp_box_x = (x_inc) * (start_i + 1.5)
    #     #     bp_box_x2 = (x_inc) * (end_i + 1.5)
    #     #     bp_box_y = bp_y1
    #     #     bp_box_w = bp_box_x2 - bp_box_x
    #     #     bp_box_h = (y_inc) * num_records + (y_inc * 0.8)

    #     #     # Greyed out rectangle
    #     #     rect = patches.Rectangle(
    #     #         (box_x, box_y),
    #     #         box_w,
    #     #         box_h,
    #     #         alpha=0.4,
    #     #         fill=True,
    #     #         facecolor="black",
    #     #         edgecolor="none",
    #     #         lw=linewidth,
    #     #         ls="--",
    #     #     )
    #     #     ax.add_patch(rect)

    #     #     # Dashed rectangle outline
    #     #     rect = patches.Rectangle(
    #     #         (box_x, box_y),
    #     #         box_w,
    #     #         box_h,
    #     #         alpha=1,
    #     #         fill=False,
    #     #         facecolor="none",
    #     #         edgecolor="black",
    #     #         lw=linewidth,
    #     #         ls="--",
    #     #     )
    #     #     ax.add_patch(rect)

    #     #     # Bonus little dashed tick in center
    #     #     line_x = box_x + (box_w / 2)
    #     #     line_y1 = box_y + box_h
    #     #     line_y2 = bp_y2

    #     #     ax.plot(
    #     #         [line_x, line_x],
    #     #         [line_y1, line_y2],
    #     #         ls="--",
    #     #         lw=linewidth,
    #     #         color="black",
    #     #         clip_on=False,
    #     #     )

    #     #     x_coord = box_x + (box_x2 - box_x) / 2

    #     # ax.text(
    #     #     x_coord,
    #     #     bp_y2,
    #     #     "Breakpoint #" + str(i + 1),
    #     #     size=fs,
    #     #     ha="center",
    #     #     va="center",
    #     #     bbox=dict(
    #     #         facecolor="white", edgecolor="black", lw=linewidth, boxstyle="round"
    #     #     ),
    #     # )

    # current_y = section_y1

    # # -----------------------------------------------------------------------------
    # # SECTION: SAMPLES AND SUBSTITUTION BOXES

    # # General section dimensions
    # section_y2 = current_y
    # section_y1 = section_y2 - (num_records * y_inc) - (y_inc)

    # # Starting x point for sub boxes, we start with iter
    # sub_x = section_x

    # # Iterate through subs, which are columns in plot
    # for rec_i, rec in enumerate(barcodes_df.iterrows()):

    #     # Adjust box coord based on width
    #     sub_box_x = sub_x + (x_inc - sub_box_w) / 2

    #     # Identify base origins
    #     ref_base = rec[1]["Reference"]
    #     parent_1_base = rec[1][parent_1]
    #     parent_2_base = rec[1][parent_2]

    #     # y coordinates for substitutions box
    #     sub_y = section_y2

    #     # Iterate through samples, which are rows in plot
    #     for i, label in enumerate(records):

    #         # Adjust box coord based on height
    #         sub_box_y = sub_y - (y_inc - sub_box_h) / 2

    #         # Add extra gap after recombinant parents
    #         if i == y_break:
    #             sub_box_y -= y_inc

    #         base = rec[1][label]

    #         if label == "Reference":
    #             box_c = ref_color
    #         elif base == parent_1_base:
    #             if base == ref_base:
    #                 box_c = parent_1_ref_color
    #             else:
    #                 box_c = parent_1_mut_color
    #         elif base == parent_2_base:
    #             if base == ref_base:
    #                 box_c = parent_2_ref_color
    #             else:
    #                 box_c = parent_2_mut_color
    #         else:
    #             box_c = "white"

    #         rect = patches.Rectangle(
    #             (sub_box_x, sub_box_y),
    #             sub_box_w,
    #             sub_box_h,
    #             alpha=0.90,
    #             fill=True,
    #             edgecolor="none",
    #             facecolor=box_c,
    #         )
    #         ax.add_patch(rect)

    #         # Debug
    #         rect = patches.Rectangle(
    #             (sub_x, sub_y - y_inc),
    #             x_inc,
    #             y_inc,
    #             alpha=1,
    #             fill=True,
    #             edgecolor="black",
    #             facecolor="none",
    #         )
    #         ax.add_patch(rect)

    #         text_y = sub_box_y + (sub_box_h / 2)

    #         # On the first time parsing sub, write sample label
    #         if rec_i == 0:
    #             # If parent, also include clade
    #             if label == parent_1:
    #                 label = "{} ({})".format(parent_1, parent_1_clade_lineage)
    #             elif label == parent_2:
    #                 label = "{} ({})".format(parent_2, parent_2_clade_lineage)
    #             ax.text(
    #                 section_label_x,
    #                 text_y,
    #                 label,
    #                 size=fs,
    #                 ha="right",
    #                 va="center",
    #             )
    #         text_x = sub_box_x + (sub_box_w / 2)

    #         # Draw sub text bases
    #         ax.text(text_x, text_y, base, size=fs, ha="center", va="center")

    #     sub_x += x_inc

    # current_y = section_y1

    # # -----------------------------------------------------------------------------
    # # SECTION: SUBSTITUTION X AXIS TICKS

    # section_y2 = current_y - (y_inc * 0.25)
    # section_y1 = section_y2 - y_inc

    # tick_y1 = section_y2
    # tick_y2 = section_y2 - (y_inc * 0.25)

    # tick_x = section_x - x_inc + (x_inc - sub_box_w) / 2

    # tick_text_y = tick_y2

    # # Iterate through subs, which are columns in plot
    # for rec_i, rec in enumerate(barcodes_df.iterrows()):

    #     tick_x += x_inc

    #     genome_coord = rec[1]["coord"]

    #     ax.plot([tick_x, tick_x], [tick_y1, tick_y2], lw=linewidth, c="black")
    #     ax.text(
    #         tick_x,
    #         tick_y2,
    #         str(genome_coord) + " ",
    #         size=fs,
    #         ha="center",
    #         va="top",
    #         rotation=90,
    #     )

    # current_y = section_y1

    # # -----------------------------------------------------------------------------
    # # Legend

    # # Build legend colors, labels from bottom up
    legend_colors = [
        parent_2_ref_color,
        parent_2_mut_color,
        parent_1_ref_color,
        parent_1_mut_color,
        ref_color,
    ]
    legend_labels = [
        parent_2 + " Reference",
        parent_2 + " Mutation",
        parent_1 + " Reference",
        parent_1 + " Mutation",
        "Reference",
    ]

    # # General section dimensions
    # legend_x = x_inc * 0.5
    # legend_y1 = current_y
    # legend_y2 = legend_y1 + (y_inc * (len(legend_labels) + 1))
    # legend_w = legend_x + (x_inc * 8) # 8 inc is an estimate
    # legend_h = legend_y2 - legend_y1

    # # Write label to left
    # ax.text(
    #     -50,
    #     legend_y1 + (legend_y2 - legend_y1) / 2,
    #     "Legend",
    #     size=fontsize,
    #     ha="right",
    #     va="center",
    # )

    # # Draw Legend Frame
    # legend_frame = patches.Rectangle(
    #     (legend_x, legend_y1),
    #     legend_w,
    #     legend_h,
    #     edgecolor="black",
    #     lw=linewidth,
    #     facecolor="none",
    # )
    # ax.add_patch(legend_frame)

    # # Coordinates for section elements
    # box_x = legend_x + x_inc
    # box_y = legend_y1 + (y_inc * 0.5)
    # text_x = box_x + x_inc
    # text_y =box_y + (sub_box_h / 2)

    # for color,label in zip(legend_colors,legend_labels):

    #     box = patches.Rectangle(
    #     (box_x, box_y), sub_box_w, sub_box_h, ec="black", lw=linewidth, fc=color)
    #     ax.add_patch(box)
    #     ax.text(text_x, text_y, label, size=fontsize, ha="left", va="center")

    #     box_y += y_inc
    #     text_y += y_inc

    # current_y = legend_y2

    # -----------------------------------------------------------------------------
    # PLOT WRAPUP

    # Flip y axis
    # plt.gca().invert_yaxis()

    # Hide the axes lines
    for spine in ax.spines:
        ax.spines[spine].set_visible(False)

    # Remove axis labels
    ax.set_xlabel("")
    ax.set_ylabel("")

    # Remove tick labels
    ax.set_xticks([])
    ax.set_yticks([])

    # # Set the Y axis limits to the genome length
    # Extra y_inc for additional row modules
    # Inverted y axis
    x_min = 0
    x_max = genome_length + x_inc
    ax.set_xlim(x_min, x_max)

    y_min = current_y
    y_max = y_inc
    # y_max = (num_records * y_inc)
    ax.set_ylim(y_min, y_max)

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
