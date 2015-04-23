#
# Implementation of the IMPUTE2 format into Python objects.
#
# This file is part of gepyto.
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
HARD_CALL = "hard_call"


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
        - is_chr23: Not implemented yet, but dosage is computed differently
                    for sexual chromosomes for men (hemizygote).
        - sex_vector: Not implemented yet, but this is a vector representing
                      the gender of every sample (for dosage computation on
                      sexual chromosomes).

    .. warning::

        Be careful with the :py:func:`Impute2File.as_matrix()` function as it
        will  try to load the WHOLE Impute2 file in memory.

    """

    def __init__(self, fn, mode=LINE, **kwargs):
        self._filename = fn

        if fn.endswith(".gz"):
            opener = gzip.open
        else:
            opener = functools.partial(open, mode="r")

        self._file = opener(fn)

        assert mode in (DOSAGE, LINE, HARD_CALL)
        self._mode = mode

        # The special function arguments
        self.dosage_arguments = {}
        self.hard_calls_arguments = {}

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

        elif self._mode is HARD_CALL:
            # Parse kwargs that can be passed to the _compute_hard_calls
            # function.
            kw = ("prob_threshold",)
            for keyword in kwargs:
                if keyword in kw:
                    self.hard_calls_arguments[keyword] = kwargs[keyword]
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
            information_fields = info.keys()
            snp_vector_list.append(v)
            snp_info_list.append([info[k] for k in information_fields])

        m = np.array(snp_vector_list)  # snp x sample
        m = m.T  # We transpose to get sample x snp matrix (standard for stats)

        # Make the information df.
        df = pd.DataFrame(snp_info_list, columns=information_fields)

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
            return _compute_dosage(_read_impute2_line(line),
                                   **self.dosage_arguments)

        elif self._mode is HARD_CALL:
            return _compute_hard_calls(_read_impute2_line(line),
                                       **self.hard_calls_arguments)

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
    assert type(line) is _Line

    # We set the dosage as being the expected number of "a2" alleles for
    # all samples.
    # Given that the rows represent p_aa, p_ab, p_bb, we have that
    # E(#b) = 2 * p_bb + p_ab, the dosage.
    dosage = 2 * line.probabilities[:, 2] + line.probabilities[:, 1]

    # Probability filtering. We mask samples with low probabilities.
    if prob_threshold > 0:
        dosage[~np.any(line.probabilities > prob_threshold, axis=1)] = np.nan

    # Compute the maf.
    mac = np.nansum(dosage)
    maf = mac / (2 * np.sum(~np.isnan(dosage)))

    # If maf > 0.5, we need to flip.
    major, minor = (line.a1, line.a2)
    if maf > 0.5:
        maf = 1 - maf
        dosage = 2 - dosage  # 0 -> 2, 1 -> 1, 2 -> 0.
        major, minor = minor, major

    if is_chr23:
        # TODO: Implement this, please
        raise NotImplementedError("dosage for chromosome 23 is not yet "
                                  "supported (because dosage is computed "
                                  "differently for males on chromosome 23)")

    return (dosage, {
        "major": major,
        "minor": minor,
        "maf": maf,
        "minor_allele_count": mac,
        "name": line.name,
        "chrom": line.chrom,
        "pos": line.pos,
    })


def _compute_hard_calls(line, prob_threshold=0):
    """Computes hard calls from probabilities (IMPUTE2)."""
    # Getting the possible genotypes
    possible_geno = np.array([" ".join([line.a1] * 2),
                              " ".join([line.a1, line.a2]),
                              " ".join([line.a2] * 2)])

    # The final genotype
    final_geno = possible_geno[np.argmax(line.probabilities, axis=1)]

    # The threshold
    low_quality = np.max(line.probabilities, axis=1) < prob_threshold
    final_geno[low_quality] = "0 0"

    return (
        final_geno,
        {"name": line.name,
         "chrom": line.chrom,
         "pos": line.pos},
    )


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
