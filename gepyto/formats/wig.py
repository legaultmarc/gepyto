"""
Parser for Wiggle Track Format files.
"""

# This file is part of gepyto.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.

__author__ = "Marc-Andre Legault"
__copyright__ = ("Copyright 2014 Marc-Andre Legault and Louis-Philippe "
                 "Lemieux Perreault. All rights reserved.")
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


import os

import pandas as pd
import six

from ..structures.region import _Segment, Region



class WiggleFile(object):
    """Parser for WIG files.
    
    This returns a pandas dataframe with all the necessary information. In the
    process, all the inherent compactness of the Wiggle format is lost in
    exchange for an easier to manage representation. This means that more
    efficient parsers should be used for large chunks of data.

    """
    def __init__(self, stream):
        self.stream = stream
        if isinstance(stream, six.string_types):
            if os.path.isfile(stream):
                self.stream = open(stream, "r")
            else:
                raise IOError("Can't find file '{}'.".format(stream))

        mode, first_header = self._parse_header(next(self.stream))
        if mode == "fixedStep":
            self.data = self._parse_fixed_step(header=first_header)
        else:
            raise NotImplementedError("fixedStep is the only implemented mode "
                                      "for now.")

        # Use categories for chrom to save space.
        self.data["chrom"] = self.data["chrom"].astype("category")

        # Check if regions or only 1 bases
        # If so use pos instead of start, end.
        if (self.data["start"] == self.data["end"]).all():
            self.data = self.data.drop("end", axis=1)
            self.data.columns = ("chrom", "pos", "value")

    def __enter__(self):
        return self

    def __exit__(self, *params):
        self.close()

    def close(self):
        # This will close the file if it's a file.
        try:
            self.stream.close()
        except AttributeError:
            pass

    def as_dataframe(self):
        return self.data


    def _parse_fixed_step(self, header=None):
        data = []
        for line in self.stream:
            if self._is_header(line):
                mode, header = self._parse_header(line)
                assert (
                    mode == "fixedStep"
                ), "Can't change mode after parsing started."

            else:
                data.append((
                    header["chrom"],
                    header["pos"],
                    header["pos"] + header["span"] - 1,
                    float(line.rstrip())
                ))
                header["pos"] += header["step"]

        return pd.DataFrame(
            data,
            columns=("chrom", "start", "end", "value")
        )


    @staticmethod
    def _parse_header(line):
        line = line.rstrip().split()
        mode = line[0]

        line = line[1:]
        header = dict([field.split("=") for field in line])
        header["start"] = int(header["start"])
        header["step"] = int(header["step"])
        header["span"] = int(header.get("span", 1))

        header["pos"] = header["start"]

        return mode, header

    @staticmethod
    def _is_header(line):
        return (
            line.startswith("variableStep") or line.startswith("fixedStep")
        )
