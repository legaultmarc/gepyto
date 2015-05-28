#
# Utility functions for parsers.
#
# This file is part of gepyto.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.

from __future__ import division, print_function

__author__ = "Marc-Andre Legault"
__copyright__ = ("Copyright 2014 Marc-Andre Legault and Louis-Philippe "
                 "Lemieux Perreault. All rights reserved.")
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


import gzip
import functools

def get_opener(fn):
    """Infers the filetype and returns an opened file object.

    :param fn: The filename.
    :type fn: str

    :returns: An opened file object.
    :rtype: file

    The current implementation can only detect gzip if the first three bytes
    are b"\x1F\x8B\x08". Alternatively, it will use the stdlib open.

    """

    with open(fn, "rb") as f:
        magic = f.read(3)

    if magic == b"\x1F\x8B\x08":
        return gzip.open

    return functools.partial(open, mode="r")
