#
# Utilities to interact with the UCSC database.
#
# This file is part of gepyto.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


from __future__ import division


__author__ = "Marc-Andre Legault"
__copyright__ = ("Copyright 2014 Marc-Andre Legault and Louis-Philippe "
                 "Lemieux Perreault. All rights reserved.")
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["UCSC", ]


import collections

import numpy as np

from .. import settings
from ..structures.region import Region


class UCSC(object):
    """Provides raw access to the UCSC MySQL database.

    The database will be set to the `db` parameter or to the `BUILD` as defined
    in `gepyto`'s settings.

    Later versions could implement features like throttling, but for now this
    is a very simple interface.

    """
    def __init__(self, db=None):
        import pymysql

        if db is None:
            db = settings.BUILD

        if db.lower() == "grch37":
            db = "hg19"
        elif db.lower() == "grch38":
            db = "hg38"
        else:
            raise Exception("Invalid genome reference '{}'".format(db))

        self.con = pymysql.connect(user="genome",
                                   host="genome-mysql.cse.ucsc.edu",
                                   database=db)

        self.cur = self.con.cursor()

    def raw_sql(self, sql, params):
        """Execute a raw SQL query."""
        self.cur.execute(sql, params)
        return self.cur.fetchall()

    def query_gap_table(self, chromosome, ucsc_type):
        """Get either the "telomere" or "centromere" of a given chromosome.

        """
        if not chromosome.startswith("chr"):
            chromosome = "chr" + chromosome

        valid_types = ("telomere", "centromere")
        if ucsc_type not in valid_types:
            msg = "'{}' is not a valid type: use {}.".format(
                ucsc_type,
                valid_types
            )
            raise TypeError(msg)

        return self.raw_sql(
            ("SELECT chromStart + 1, chromEnd + 1 "
             "FROM gap "
             "WHERE chrom=%s AND type=%s"),
            (chromosome, ucsc_type),
        )

    def close(self):
        self.con.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


def get_telomere(chromosome):
    """Returns a Noncontiguous region representing the telomeres of a
    chromosome.

    :param chromosome: The chromosome, _e.g._ "3"
    :type chromosome: str

    :returns: A region corresponding to the telomeres.
    :rtype: :py:class:`Region`

    This is done by connecting to the UCSC MySQL server.

    """
    chromosome = str(chromosome)
    with UCSC() as ucsc:
        telomeres = ucsc.query_gap_table(chromosome, "telomere")

    assert len(telomeres) == 2, ("UCSC did not return two telomeres (chrom={}"
                                 ").".format(chromosome))

    # Create a region for both telomeres and use a union to return the full
    # region.
    telo1 = Region(chromosome.lstrip("chr"), telomeres[0][0], telomeres[0][1])
    telo2 = Region(chromosome.lstrip("chr"), telomeres[1][0], telomeres[1][1])
    return telo1.union(telo2)


def get_centromere(chromosome):
    """Returns a contiguous region representing the centromere of a chromosome.

    :param chromosome: The chromosome, _e.g._ "3"
    :type chromosome: str

    :returns: A region corresponding to the centromere.
    :rtype: :py:class:`Region`

    This is done by connecting to the UCSC MySQL server.

    """
    chromosome = str(chromosome)
    with UCSC() as ucsc:
        centromere = ucsc.query_gap_table(chromosome, "centromere")

    assert len(centromere) == 1, "UCSC returned {} centromere(s).".format(
        len(centromere)
    )
    centromere = centromere[0]

    return Region(chromosome.lstrip("chr"), centromere[0], centromere[1])


def get_phylop_100_way(region):
    """Get a vector of phyloP conservation scores for a given region.

    Scores represent the -log(p-value) under a H0 of neutral evolution.
    Positive values represent conservation and negative values represent
    fast-evolving bases.

    The UCSC MySQL database only contains aggregate scores for chunks of
    1024 bases. We return the results for the subset of the required region
    that is fully contained in the UCSC bins.

    Because UCSC uses 0-based indexing, we adjust the gepyto region before
    querying. This means that the user should use 1-based indexing, as
    usual when creating the Region object.

    """
    with UCSC() as ucsc:
        sql = (
            "SELECT * FROM phyloP100wayAll WHERE "
            "   chrom=%s AND "
            "   chromStart>=%s AND "
            "   chromEnd<=%s "
        )

        chrom = region.chrom
        if not chrom.startswith("chr"):
            chrom = "chr{}".format(chrom)

        ucsc.cur.execute(sql, (chrom, region.start - 1, region.end - 1))

        n = ucsc.cur.rowcount

        if not n:
            return  # No results.

        results = iter(ucsc.cur)

    phylop = np.empty(n)

    Row = collections.namedtuple(
        "Row",
        ("bin", "chrom", "chromStart", "chromEnd", "name", "span", "count",
         "offset", "file", "lowerLimit", "dataRange", "validCount",
         "sumData", "sumSquares")
    )

    start = end = None
    for i, row in enumerate(results):
        row = Row(*row)  # Bind column names.

        # Adjust types so we can do integer operations.
        row_start = int(row.chromStart)
        row_end = int(row.chromEnd)
        sum_data = int(row.sumData)
        valid_count = int(row.validCount)

        end = row_end

        if start is None:
            start = row_start + 1

        phylop[i] = sum_data / valid_count

    chrom = chrom[3:]
    return {
        "phylop_scores": phylop,
        "n_bins": n,
        "actual_region": Region(chrom, start, end)
    }
