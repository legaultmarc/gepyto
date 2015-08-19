
# Utilities to handle chromosome ideogram plotting.

# This file is part of gepyto.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import operator
from collections import namedtuple

import matplotlib.pyplot as plt
import matplotlib.patches as patches

from ..db import ucsc


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = ("Copyright 2014 Marc-Andre Legault and Louis-Philippe "
                 "Lemieux Perreault. All rights reserved.")
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


# The default cytoband color
_ucsc_cytoband_color = {
    "acen":    "#963232",
    "gneg":    "#E3E3E3",
    "gpos100": "#000000",
    "gpos25":  "#8E8E8E",
    "gpos50":  "#555555",
    "gpos75":  "#393939",
    "gvar":    "#000000",
    "stalk":   "#7F7F7F",
}


def plot_chr_ideogram(ax, chrom, loc="bottom", cyto_color=_ucsc_cytoband_color,
                      proportion=0.05):
    """Plot the chromosome ideogram at the bottom of the axe.

    :param ax: the axe that will contain the genes in the region.
    :type ax: :py:class:`matplotlib.axes.Axes`

    :param chrom: the chromosome for which to plot the ideogram
    :type chrom: str

    :param loc: the location of the ideogram (bottom, top, left, right)
    :type loc: str

    :param cyto_color: color of the cytobands
    :type cyto_color: dict

    :param proportion: the proportion of the axe for the ideogram
    :type proportion: float

    .. note::

        The ``cyto_color`` dictionary can contain only a subset of band/color
        values (default UCSC colors will be taken for missing values).

    .. note::

        The size of the ideogram is determined using the proportion. If the
        location is top or bottom, the height is determined as
        :math:`(y_{max} - y_min) * proportion`. If the location is left or
        right, the width is determined as
        :math:`(x_{max} - x_{min}) * proportion`.

    """
    # Checking the location value
    assert loc in {"bottom", "top", "left", "right"}, "invalid location"

    # Getting the axe limits (so that it doesn't change by mistake)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    # Getting the chromosome's bands
    for band in _get_ucsc_cytobands(chrom=chrom):
        # Getting the band color
        col = cyto_color.get(
            band.gie_stain,
            _ucsc_cytoband_color.get(band.gie_stain, "#000000"),
        )

        # Getting the positions and the dimension of the rectangle
        if loc == "bottom":
            x = band.start + 1
            y = min(ylim)
            width = band.end - band.start
            height = operator.sub(*ylim[::-1]) * proportion

        elif loc == "top":
            x = band.start + 1
            y = max(ylim) - (operator.sub(*ylim[::-1]) * proportion)
            width = band.end - band.start
            height = operator.sub(*ylim[::-1]) * proportion

        elif loc == "left":
            x = min(xlim)
            y = band.start + 1
            width = operator.sub(*xlim[::-1]) * proportion
            height = band.end - band.start

        elif loc == "right":
            x = max(xlim) - (operator.sub(*xlim[::-1]) * proportion)
            y = band.start + 1
            width = operator.sub(*xlim[::-1]) * proportion
            height = band.end - band.start

        # Adding the band
        ax.add_patch(patches.Rectangle((x, y), width, height, lw=0, ec=col,
                                       fc=col, zorder=0))

    # Setting the axe limits as were (just to be sure)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)


def _get_ucsc_cytobands(chrom):
    """Gets the cytobands from UCSC's MySQL database.

    :param chrom: the chromosome
    :type chrom: str

    .. note::

        This function is a generator.

    """
    # Encoding the chromosome
    chrom = str(chrom)
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom

    res = None
    with ucsc.UCSC() as ucsc_connection:
        res = ucsc_connection.raw_sql(
            ("SELECT chrom, chromStart, chromEnd, name, gieStain "
             "FROM cytoBandIdeo "
             "WHERE chrom=%s"),
            (chrom, ),
        )

    # Creating a named tuple
    CytoBand = namedtuple("CytoBand",
                          ["chrom", "start", "end", "name", "gie_stain"])
    for values in res:
        yield CytoBand(*values)
