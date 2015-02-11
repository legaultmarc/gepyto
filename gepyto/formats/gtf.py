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


import re
import gzip
import functools
import datetime
import logging
from collections import namedtuple

from ..structures.sequences import Sequence


class InvalidGTF(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


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

        self._line_accumulator = None

        self._read_headers()

    def __next__(self):
        if self._line_accumulator is None:
            line = next(self._file)
            if line is None:
                raise StopIteration()
        else:
            line = self._line_accumulator
            self._line_accumulator = None

        return GTFFile.parse_line(line)

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
        if len(line) < 8:
            raise InvalidGTF("Mandatory fields are missing.")

        seqname, source, feature, start, end, score, strand, frame = line[:8]
        attributes = None

        # Check the format and parse the attributes field (optional).
        if len(line) >= 9:
            # Merge attribute fields if tabs are in them.
            line[8] = " ".join(line[8:])

            # Extra attributes are available.
            attributes = line[8].strip().split(";")
            parsed_attributes = {}
            for attr in attributes:
                if attr == "":  # If lines finish with a ';'
                    continue
                attr = attr.strip()
                tag = re.match(r"^([A-Za-z][A-Za-z0-9_]*)\s", attr)
                if tag is None:
                    raise InvalidGTF("Invalid tag in attributes field \"{}\"."
                                     "".format(attr))
                tag = tag.group(1)
                value = attr[len(tag):].strip()
                value = value.replace('"', '')
                parsed_attributes[tag] = value

            attributes = parsed_attributes

        # Type checks.
        if seqname.startswith("chr"):
            seqname = seqname[3:]

        try:
            start = int(start)
            end = int(end)
            assert start <= end
            assert start >= 1

            score = float(score) if score != "." else None

            strand = strand if strand != "." else None
            assert strand in ("+", "-") or strand is None

            frame = frame if frame != "." else None
            assert frame in ("0", "1", "2") or frame is None

        except AssertionError as e:
            message = ("Some fields of the GTF were invalid. The parsed "
                       "parameters are as follows:\n"
                       "\tsequname: {seqname}\n"
                       "\tsource: {source}\n"
                       "\tfeature: {feature}\n"
                       "\tstart: {start}\n"
                       "\tend: {end}\n"
                       "\tscore: {score}\n"
                       "\tstrand: {strand}\n"
                       "\tframe: {frame}\n")

            message = message.format(seqname=seqname, source=source,
                                     feature=feature, start=start, end=end,
                                     score=score, strand=strand, frame=frame)

            logging.critical(message)
            raise e

        return GTFFile.Line(seqname, source, feature, start, end, score,
                            strand, frame, attributes)


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

        self._line_accumulator = line  # This is not a header line.

