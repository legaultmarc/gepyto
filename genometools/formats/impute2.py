#
# Implementation of the IMPUTE2 format into Python objects.
#
# This file is part of genometools.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.

from __future__ import division, print_function

__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = ("Copyright 2014 Marc-Andre Legault and Louis-Philippe "
                 "Lemieux Perreault. All rights reserved.")
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"

from collections import Counter, namedtuple
from multiprocessing import Process, Queue
import functools

import numpy as np
import pandas as pd

import gzip

_Line = namedtuple(
    "Line",
    ["name", "chrom", "pos", "a1", "a2", "probabilities"]
)

# Two possible ways of iterating over Impute2File
DOSAGE = "dosage"
LINE = "line"


class Impute2File(object):
    """Class representing an Impute2File.

    This is used to either generate a dosage matrix where columns represent
    variants and rows represent samples or to read the file line by line
    using the generator syntax.

    This also implements the context manager interface.

    Usage: ::

        # Read as probabilities (Line tuples).
        with open(Impute2File(fn)) as f:
            for line in f:
                # line has name, chrom, pos, a1, a2, probabilities
                print(line)

        # Read as dosage.
        with open(Impute2File(fn), "dosage") as f:
            for dosage_vector, info in f:
                pass

        # Read as a matrix.
        with open(Impute2File(fn)) as f:
            # 1 row per sample and 1 column per variant. Values between 0 and 2
            m = f.as_matrix()

    If you use the ``dosage`` mode, you can also add additional arguments:

        - prob_threshold: Genotype probability cutoff for no call values (NaN).
        - is_chr12: Not implemented yet, but dosage is computed differently
                    for sexual chromosomes for men (hemizygote).
        - sex_vector: Not implemented yet, but this is a vector representing
                      the gender of every sample (for dosage computation on
                      sexual chromosomes).

    .. warning::

        Be careful with the :py:func:`Impute2File.as_matrix()` function as it
        will  try to load the WHOLE Impute2 file in memory.

    .. todo::

        Add thorough testing for this class.

    """

    def __init__(self, fn, mode=LINE, **kwargs):
        self._filename = fn

        if fn.endswith(".gz"):
            opener = gzip.open
        else:
            opener = functools.partial(open, mode="r")

        self._file = opener(fn)

        assert mode in (DOSAGE, LINE)
        self._mode = mode

        self.dosage_arguments = {}
        if self._mode is DOSAGE:
            # Parse kwargs that can be passed to the _compute_dosage function.
            kw = ("prob_threshold", "is_chr23", "sex_vector")
            for keyword in kwargs:
                if keyword in kw:
                    self.dosage_arguments[keyword] = kwargs[keyword]
                else:
                    self._file.close()
                    raise TypeError("__init__() got an unexpected keyword "
                        "argument '{}'".format(keyword))
        else:
            if len(kwargs) > 0:
                self._file.close()
                raise TypeError("__init__() got an unexpected keyword "
                        "argument '{}'".format(kwargs.keys()[0]))

    def as_matrix(self):
        """Creates a numpy dosage matrix from this file.

        :returns: A numpy matrix where columns represent variant dosage
                  between 0 and 2 and a dataframe describing the variants
                  (major, minor, maf).
        :type: tuple

        .. warning::
        
            This will attempt to load the whole file in memory.

        """
        prev_pos = self._file.tell()
        self._file.seek(0)
        
        prev_mode = self._mode
        self._mode = DOSAGE

        snp_vector_list = []
        snp_info_list = []
        for v, info in self:
            snp_vector_list.append(v)
            snp_info_list.append(
                (info["name"], info["major"], info["minor"], info["maf"])
            )

        m = np.array(snp_vector_list) # snp x sample
        m = m.T # We transpose to get sample x snp matrix (standard for stats)

        # Make the information df.
        df = pd.DataFrame(
            snp_info_list, columns=["name", "major", "minor", "maf"]
        )

        # Put the file like it was.
        self._file.seek(prev_pos)
        self._mode = prev_mode

        return m, df

    def __next__(self):
        line = next(self._file)
        if line is None:
            # Done with the file.
            raise StopIteration()

        if self._mode is DOSAGE:
            return _compute_dosage(
                _read_impute2_line(line), 
                **self.dosage_arguments
            )
        elif self._mode is LINE:
            return _read_impute2_line(line)

    next = __next__

    def readline(self):
        """Read a single line from the Impute2File.

        This will return either a ``Line`` including the genotype probabilities
        or a dosage vector. This depends on the `mode` (the second argument
        given to the file when it was opened).

        Available modes are ``dosage`` and ``line``.

        """
        return self.next()

    def __iter__(self):
        return self

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        
    def close(self):
        self._file.close()


def _compute_dosage(line, prob_threshold=0, is_chr23=False,
                    sex_vector=None):
    """Computes dosage from probabilities (IMPUTE2)."""
    # Type check
    type(line) is line

    # Constructing the Pandas data frame
    data = pd.DataFrame({
        "d1": line.probabilities[:, 0],
        "d2": line.probabilities[:, 1],
        "d3": line.probabilities[:, 2],
    })

    # Getting the maximal probabilities for each sample
    data["dmax"] = data.max(axis=1)

    # Computing the A1 frequency
    if is_chr23:
        # TODO: Implement this, please
        raise NotImplementedError("dosage for chromosome 23 is not yet "
                                  "supported (because dosage is computed "
                                  "differently for males on chromosome 23)")

    else:
        sub_data = data[data.dmax >= prob_threshold][["d1", "d2", "d3"]]
        count = Counter(sub_data.idxmax(axis=1))
        a1_freq = ((count["d1"] * 2) + count["d2"]) / (sum(count.values()) * 2)

    # Computing the dosage
    minor_allele = "d1"
    minor = line.a1
    major = line.a2
    maf = a1_freq
    if a1_freq >= 0.5:
        minor_allele = "d3"
        minor = line.a2
        major = line.a1
        maf = 1 - a1_freq
    data["dosage"] = (data[minor_allele] + (data.d2 / 2)) * 2

    # Setting values to NaN when max prob is lower than threshold
    data.loc[data.dmax < prob_threshold, :] = np.nan

    return (data.dosage.values, {
        "major": major,
        "minor": minor,
        "maf": maf,
        "name": line.name,
    })


def _read_impute2_line(line):
    """Parse an IMPUTE2 line (a single marker).

    :param line: a line from an impute 2 file.
    :type line: str

    :returns: A namedtuple with the name of the variant, chromosome, position,
              allele1, allele2 and probability matrix.

    The matrix of shape ``samples x 3`` represents the probability of the 
    the following genotypes ``(aa, ab, bb)``.

    """
    # Splitting the line
    row = line.rstrip("\n").split(" ")

    # Getting marker information
    name = row[1]
    chrom = row[0]
    pos = int(row[2])
    a1 = row[3]
    a2 = row[4]

    # Constructing the genotype
    prob = np.array(row[5:], dtype=float)
    prob.shape = (prob.shape[0] // 3, 3)

    return _Line(name, chrom, pos, a1, a2, prob)
