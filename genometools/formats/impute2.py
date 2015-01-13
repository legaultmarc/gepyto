#
# Implementation of the IMPUTE2 format into Python objects.
#
# This file is part of genometools.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = ("Copyright 2014 Marc-Andre Legault and Louis-Philippe "
                 "Lemieux Perreault. All rights reserved.")
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


from collections import Counter

import numpy as np
import pandas as pd


def _compute_dosage(impute2_data, prob_threshold=0, is_chr23=False,
                    sex_vector=None):
    """Computes dosage from probabilities (IMPUTE2)."""
    # Some assertions
    assert "probabilities" in impute2_data
    assert "a1" in impute2_data
    assert "a2" in impute2_data

    # Constructing the Pandas data frame
    data = pd.DataFrame({
        "d1": impute2_data["probabilities"][:, 0],
        "d2": impute2_data["probabilities"][:, 1],
        "d3": impute2_data["probabilities"][:, 2],
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
    minor = impute2_data["a1"]
    major = impute2_data["a2"]
    maf = a1_freq
    if a1_freq >= 0.5:
        minor_allele = "d3"
        minor = impute2_data["a2"]
        major = impute2_data["a1"]
        maf = 1 - a1_freq
    data["dosage"] = (data[minor_allele] + (data.d2 / 2)) * 2

    # Setting values to NaN when max prob is lower than threshold
    data.loc[data.dmax < prob_threshold, :] = np.nan

    return {
        "major": major,
        "minor": minor,
        "maf": maf,
        "dosage": data.dosage.values,
    }


def _read_impute2_line(line):
    """Parse an IMPUTE2 line (a single marker)."""
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

    return {
        "name": name,
        "chrom": chrom,
        "pos": pos,
        "a1": a1,
        "a2": a2,
        "probabilities": geno,
    }
