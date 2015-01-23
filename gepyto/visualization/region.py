
# Utilities to handle genomic region plotting.

# This file is part of gepyto.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = ("Copyright 2014 Marc-Andre Legault and Louis-Philippe "
                 "Lemieux Perreault. All rights reserved.")
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


import re
import logging
from collections import defaultdict

import pandas as pd

from ..utils.genes import ensembl_genes_in_region


def plot_genes_in_region(figure, axe, region, build):
    """Plot genes in a region.

    :param figure: The figure that contains the axe.
    :type region: :py:class:`matplotlib.figure.Figure`

    :param axe: The axe that will contain the genes in the region.
    :type axe: :py:class:`matplotlib.axes.Axes`

    :param region: The region
    :type region: str

    :param build: The build
    :type build: str

    """
    # Fetching the gene in the region
    logging.info("Fetching genes in region {} ({})".format(region, build))
    genes = ensembl_genes_in_region(region, bare=True, build=build)
    logging.debug("Found {:,d} genes".format(len(genes)))

    # Getting the start and end positions
    r = re.search(r"(chr)?(\w+):(\d+)-(\d+)", region)
    start = r.group(3)
    end = r.group(4)
    try:
        start = int(start)
        end = int(end)
    except ValueError:
        raise ValueError("{}: malformed region".format(region))

    # Creating a Pandas data frame
    genes = [(g.symbol, g.start, g.end, g.strand, g.biotype) for g in genes]
    genes = pd.DataFrame(genes, columns=["symbol", "start", "end", "strand",
                                         "biotype"])
    genes = genes.sort(columns=["start", "end"])

    # Adding the X limits
    axe.set_xlim(start, end)

    # Getting the renderer
    renderer = figure.canvas.get_renderer()

    # Adding genes one by one
    last_t_obj = {}
    last_end = defaultdict(int)
    for i in range(genes.shape[0]):
        gene_start = genes.iloc[i, :].start
        gene_end = genes.iloc[i, :].end
        gene_name = genes.iloc[i, :].symbol

        # Checking the starting position of the gene
        if gene_start < start:
            gene_start = start

        # Checking the ending position of the gene
        if gene_end > end:
            gene_end = end

        # Updating the gene label
        gene_label = None
        if genes.iloc[i, :].strand == 1:
            gene_label = gene_name + r"$\rightarrow$"
        else:
            gene_label = r"$\leftarrow$" + gene_name

        # We find the first j where we can put the line
        j = 0
        while True:
            if last_end[j] < gene_start:
                break
            j -= 1

        # Trying to put the label there
        t = axe.text((gene_start + gene_end) / 2, j - 0.15, gene_label,
                     fontsize=5, ha="center", va="top")

        # Is there a bbox in this location?
        if j in last_t_obj:
            # Getting the bbox
            bb = t.get_window_extent(renderer=renderer)
            last_bb = last_t_obj[j].get_window_extent(renderer=renderer)

            while last_bb.overlaps(bb):
                # BBoxes overlap
                logging.info("{} overlaps".format(gene_name))
                j -= 1
                t.set_y(j - 0.15)

                # Last j?
                if j not in last_t_obj:
                    break

                # Need to update both bboxes
                bb = t.get_window_extent(renderer=renderer)
                last_bb = last_t_obj[j].get_window_extent(renderer=renderer)

        # Plotting the line
        logging.info("Putting {} at position {}".format(gene_name, j))
        marker = "-"
        other_param = {}
        if (gene_end - gene_start) < 3e-3:
            # Too small
            marker = "s"
            other_param["ms"] = 1.8
        axe.plot([gene_start, gene_end], [j, j], marker, lw=2, color="#000000",
                 clip_on=False, **other_param)

        # Saving the last position (last end and bbox)
        last_end[j] = gene_end + 3e-3
        last_t_obj[j] = t
