
# Structure to handle genomic regions or ranges.
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

__all__ = []


import re

from .. import settings
from ..db import ucsc
from . import sequences


class _Segment(object):
    def __init__(self, chrom, start, end):
        chrom = str(chrom)
        assert re.match(settings.CHROM_REGEX, chrom)
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        assert self.start < self.end

    def overlaps_with(self, segment):
        return (self.chrom == segment.chrom and
                self.start <= segment.end and
                self.end >= segment.start)

    def __eq__(self, seg):
        return (self.chrom == seg.chrom and self.start == seg.start and
                self.end == seg.end)

    def __ne__(self, seg):
        return not self.__eq__(seg)

    def to_sequence(self):
        return sequences.Sequence.from_reference(self.chrom, self.start,
                                                 self.end)

    @staticmethod
    def merge_segments(li):
        """Merge overlapping segments in a sorted list."""

        for i in range(len(li) - 1):
            cur = li[i]
            nxt = li[i + 1]
            if nxt.start < cur.start:
                raise Exception("Only sorted lists of segments can be "
                                "merged. Sort using the position first.")

            if cur.chrom != nxt.chrom:
                raise Exception("Can't merge segments from different "
                                "chromosomes.")

        merged_segments = []
        i = 0
        while i < len(li) - 1:
            j = i

            # Walk as long as the segments are overlapping.
            cur = li[i]
            nxt = li[i + 1]
            block = [cur.start, cur.end]
            if nxt.start <= cur.end:
                block[1] = max(block[1], nxt.end)

            while nxt.start <= block[1] and j + 1 < len(li) - 1:
                block[1] = max(block[1], nxt.end)
                j += 1
                cur = li[j]
                nxt = li[j + 1]

            merged_segments.append(
                _Segment(cur.chrom, block[0], block[1])
            )
            i = j + 1

        if li[-1].start > li[-2].end:
            merged_segments.append(li[-1])

        return merged_segments

    def __repr__(self):
        return "<_Segment object chr{}:{}-{}>".format(self.chrom, self.start,
                                                      self.end)


class Region(object):
    def __init__(self, chrom, start, end):
        self.segments = []
        self.segments.append(_Segment(chrom, start, end))

    def union(self, region):
        segments = self.segments + region.segments
        segments = sorted(segments, key=lambda x: x.start)
        segments = _Segment.merge_segments(segments)
        return Region._from_segments(segments)

    def overlaps_with(self, region):
        for seg1 in self.segments:
            for seg2 in region.segments:
                if seg1.overlaps_with(seg2):
                    return True

        return False

    @property
    def is_contiguous(self):
        return len(self.segments) == 1

    @property
    def sequence(self):
        """Builds a Sequence object representing the region.

        If the region is non contiguous, a tuple of sequences is returned.
        """
        sequences = []
        for seg in self.segments:
            sequences.append(seg.to_sequence())

        return sequences[0] if len(sequences) == 1 else tuple(sequences)

    @staticmethod
    def _from_segments(segments):
        region = Region(1, 1, 2)
        region.segments = segments
        return region

    def __repr__(self):
        return "<{}Region: {}>".format(
            "Contiguous" if self.is_contiguous else "NonContiguous",
            self.segments
        )

def get_telomere(chromosome):
    """Returns a Noncontiguous region representing the telomeres of a
    chromosome.
    
    :param chromosome: The chromosome, _e.g._ "3"
    :type chromosome: str

    :returns: A region corresponding to the telomeres.
    :rtype: :py:class:`Region`

    This is done by connecting to the UCSC MySQL server.

    """
    if chromosome.startswith("chr"):
        pass
    else:
       chromosome = "chr" + chromosome

    with ucsc.UCSC() as ucsc_connection:
        telomeres = ucsc_connection.raw_sql(
            ("SELECT chromStart + 1, chromEnd + 1 "
             "FROM gap "
             "WHERE chrom=%s AND type='telomere'"),
            chromosome,
        )
        assert len(telomeres) == 2, "UCSC did not return two telomeres."

    # Create a region for both telomeres and use a union to return the full
    # region.
    telo1 = Region(chromosome.lstrip("chr"), telomeres[0][0], telomeres[0][1])
    telo2 = Region(chromosome.lstrip("chr"), telomeres[1][0], telomeres[1][1])
    return telo1.union(telo2)
