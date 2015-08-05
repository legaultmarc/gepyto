#
# Parser for Variant Call Format (VCF) files.
# As of now, focus is on the v4.1 version of the specification available at:
# http://samtools.github.io/hts-specs/VCFv4.1.pdf
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

import re

from .utils import get_opener


# This pretty regex captures and validates INFO fields from the VCF spec.
_INFO_REGEX = re.compile(
    r"^INFO=<ID=(?P<id>.+),Number=(?P<number>[0-9]+|[GA.]),"
     "Type=(?P<type>Integer|Float|Flag|Character|String),"
     "Description=\"(?P<description>.+)\">$"
)


_FILTER_REGEX = re.compile(
    r"^FILTER=<ID=(?P<id>.+),Description=\"(?P<description>.+)\">$"
)


class VCFFile(object):
    """Class representing a VCF file."""

    def __init__(self, fn):
        self._filename = fn
        self._file = get_opener(fn)(fn)  # get_opener returns a function.

        self.info = {}
        self.filters = {}

        self._parse_headers()

    def _parse_headers(self):
        """Parse the VCF headers.

        Recognized meta-information lines are:

            - fileformat (maps to 'version' attribute)
            - INFO

        """
        cur = self._file.tell()
        line = next(self._file)

        while line.startswith("#"):
            # We are reading a header line.
            line = line.lstrip("#")
            line = line.rstrip()

            if line.startswith("fileformat"):
                self.version = line.split("=")[1]

            elif line.startswith("INFO"):
                info = self._parse_info(line)
                info_id = info.pop("id")
                self.info[info_id] = info

            elif line.startswith("FILTER"):
                applied_filter = self._parse_filter(line)
                filter_id = applied_filter.pop("id")
                self.filters[filter_id] = applied_filter

            cur = self._file.tell()
            line = next(self._file)

        self._file.seek(cur)

    def _parse_info(self, s):
        """Parse and validate an INFO field into a standardized data structure.

        INFO fields have an ID, a number representing their length (A for
        allele and G for genotype are special lengths), a type (Integer, Float,
        Flag, Character or String) and a description.

        """
        # Validate the structure of the line.
        matcher = _INFO_REGEX.match(s)
        if not matcher:
            raise ValueError("The following INFO line did not validate: \n\n"
                             "{}".format(s))

        # Parse the different fields.
        return {
            "id": matcher.group("id"),
            "number": matcher.group("number"),
            "type": matcher.group("type"),
            "description": matcher.group("description"),
        }

    def _parse_filter(self, s):
        """Parse and validate a FILTER field.

        Expected format is ##FILTER<ID=ID,Description="description">

        """
        # Validate.
        matcher = _FILTER_REGEX.match(s)
        if not matcher:
            raise ValueError("The following FILTER line did not validate: \n\n"
                             "{}".format(s))

        # Parse the fields.
        return {
            "id": matcher.group("id"),
            "description": matcher.group("description"),
        }
