#
# Parser for GTF/GFF files.
#
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


import gzip
import functools
import datetime
from collections import namedtuple

from ..structures.sequences import Sequence

class GTFFile(object):
    """Class representing a GTF file."""

    Line = namedtuple("_Line", ["seqname", "source", "features", "start",
                                "end", "score", "strand", "frame",
                                "attributes"])

    def __init__(self, fn):
        self._filename = fn

        if fn.endswith(".gz"):
            opener = gzip.open
        else:
            opener = functools.partial(open, mode="r")

        self._file = opener(fn)

        self._read_headers()

    def __next__(self):
        line = next(self._file)
        if line is None:
            raise StopIteration()

        return parse_line(line)

    next = __next__

    def readline(self):
        return self.next()

    def __iter__(self):
        return self

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        self._file.close()

    @staticmethod
    def parse_line(line):
        assert not line.startswith("#")
        if line.startswith("chr"):
            line = line[3:]

        line = line.strip().split("\t")
        seqname, source, feature, start, end, score, strand, frame = line[:8]

    def _read_headers(self):
        """Skip generic headers and parse paraseable headers of the GTF file.

        Recognized headers are `described here
        <http://www.sanger.ac.uk/resources/software/gff/spec.html#t_2>`_ :

        They are:

        - gff-version
        - source-version
        - date
        - Type
        - DNA
        - RNA
        - Protein

        If available, all of these will be parsed as attributes to the GTFFile
        object.

        """
        line = next(self._file)
        while line.startswith("#"):

            if not line.startswith("##"):
                line = next(self._file)
                continue

            line = line.lstrip("#")
            line = line.split()

            if line[0] == "gff-version":
                self.gff_version = int(line[1])

            elif line[0] == "source-version":
                if not hasattr(self, "source_version"):
                    self.source_version = {}
                self.source_version[line[1]] = line[2]

            elif line[0] == "date":
                self.date = datetime.datetime.strptime(line[1], "%Y-%m-%d")

            elif line[0] == "Type":
                if line[1] not in ("DNA", "Protein", "RNA"):
                    logging.warning("Unknown Type in GTF: '{}'. Known "
                                    "types are 'DNA', 'Protein' and 'RNA'."
                                    "".format(line[1]))
                if len(line) == 2:
                    self.type = line[1]
                elif len(line == 3):
                    if not hasattr(self, "type"):
                        self.type = {}
                    self.type[line[2]] = line[1]

            elif line[0] in ("DNA", "RNA", "Protein"):
                seq_type = line[0]
                name = line[1]
                line = next(self._file).lstrip("#").strip()
                seq = ""

                while line != "end-{}".format(seq_type):
                    seq += line.lstrip("#").strip()
                    try:
                        line = next(self._file).lstrip("#").strip()
                    except StopIteration:
                        raise Exception("'{}' sequence header found, but "
                                        "not closed.".format(name))

                if not hasattr(self, "sequences"):
                    self.sequences = {}

                if seq_type in ("DNA", "RNA"):
                    self.sequences[name] = Sequence(name, seq, seq_type)
                else:
                    self.sequences[name] = Sequence(name, seq, "AA")

            elif line[0] == "sequence-region":
                if not hasattr(self, "sequence_regions"):
                    self.sequence_regions = {}
                self.sequence_regions[line[1]] = (
                    int(line[2]),
                    int(line[3])
                )

            line = next(self._file)

