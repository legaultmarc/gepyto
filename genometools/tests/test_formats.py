
# This file is part of genometools.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


__author__ = "Marc-Andre Legault"
__copyright__ = ("Copyright 2014 Marc-Andre Legault and Louis-Philippe Lemieux "
                 "Perreault. All rights reserved.")
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


import unittest
import tempfile

import numpy as np

from .. import formats as fmts


def compare_lines(l1, l2):
    """Compare line tuples. """
    # Compare the first 4 elements (snpid, chrom, pos, a1, a2).
    for i in range(5):
        if l1[i] != l2[i]:
            return False

    # Compare the actual matrices.
    return (l1[5] == l2[5]).all()


def compare_dosages(self, d1, d2):
    v1, info1 = d1
    v2, info2 = d2

    for k in info1:
        if k != "maf":
            self.assertEqual(info1[k], info2[k])
        else:
            self.assertAlmostEqual(info1[k], info2[k])
    
    return (v1 == v2).all()


class TestImpute2Class(unittest.TestCase):
    """Tests the formats.impute2.Impute2File class.

    We basically compare probability matrices and dosage vectors with hard
    coded values.

    """

    def setUp(self):
        # So we have 2 SNPs and 3 samples.
        self.f = tempfile.NamedTemporaryFile("w")
        self.f.write("""
1 rs12345 1231415 A G 1 0 0 0.988 0.002 0 0 0.997 0.003
1 rs23456 3214569 T C 0.869 0.130 0 0.903 0.095 0.002 1 0 0
""".strip())
        self.f.seek(0)

        self.prob_snp1 = (
            "rs12345",
            "1",
            1231415,
            "A",
            "G",
            np.array([[1, 0, 0], [0.988, 0.002, 0], [0, 0.997, 0.003]])
            # AA, AA, AG
        )

        self.prob_snp2 = (
            "rs23456",
            "1",
            3214569,
            "T",
            "C",
            np.array([[0.869, 0.130, 0], [0.903, 0.095, 0.002], [1, 0, 0]])
            # TT, TT, TT
        )

        self.dosage_snp1 = (
            np.array([0., 0.002, 1.003]), 
            {"alt": "G", "ref": "A", "maf": 1 / 6.0}
        )

        self.dosage_snp2 = (
            np.array([0.130, 0.099, 0]), 
            {"alt": "C", "ref": "T", "maf": 0}
        )

        self.dosage_snp2_thresh = (
            np.array([np.nan, 0.099, 0]), 
            {"alt": "C", "ref": "T", "maf": 0}
        )

    def tearDown(self):
        # Closing the temporary file
        self.f.close()

    def test_syntax(self):
        """This is mostly to test the syntax for object initialization. """
        f = fmts.impute2.Impute2File(self.f.name)
        lines = [f.readline(), ]
        lines.append(f.readline())
        f.close()

        # Test the context manager and iterator compliance.
        with fmts.impute2.Impute2File(self.f.name) as f:
            for i, l in enumerate(f):
                self.assertTrue(compare_lines(l, lines[i]))

        # Make sure using invalid arguments raises error.
        self.assertRaises(
            TypeError,
            fmts.impute2.Impute2File,
            self.f.name,
            prob_threshold=0.8
        )

        # But if we're in dosage mode, this should work.
        f = fmts.impute2.Impute2File(
            self.f.name, "dosage", prob_threshold=0.8
        )
        f.close()

    def test_proba_reader(self):
        """This is to read the probabilities matrix. """

        with fmts.impute2.Impute2File(self.f.name) as f:
            for i, line in enumerate(f):
                if i == 0:
                    self.assertTrue(compare_lines(line, self.prob_snp1))
                elif i == 1:
                    self.assertTrue(compare_lines(line, self.prob_snp2))
                else:
                    raise Exception()

    def test_imputation(self):
        """This is to read as dosage vectors. """

        with fmts.impute2.Impute2File(self.f.name, "dosage") as f:
            for i, line in enumerate(f):
                if i == 0:
                    self.assertTrue(compare_dosages(
                        self,
                        line,
                        self.dosage_snp1
                    ))

                elif i == 1:
                    self.assertTrue(compare_dosages(
                        self,
                        line,
                        self.dosage_snp2
                    ))
                else:
                    raise Exception()
