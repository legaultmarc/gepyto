
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
    
    # Checking the dosage values (nan != nan in numpy)
    v1_nan = np.isnan(v1)
    v2_nan = np.isnan(v2)
    return ((v1 == v2) | (v1_nan & v2_nan)).all()


class TestImpute2Class(unittest.TestCase):
    """Tests the formats.impute2.Impute2File class.

    We basically compare probability matrices and dosage vectors with hard
    coded values.

    """

    def setUp(self):
        # So we have 2 SNPs and 3 samples. (first IMPUTE2 file)
        self.f = tempfile.NamedTemporaryFile("w")
        self.f.write("""
1 rs12345 1231415 A G 1 0 0 0.988 0.002 0 0 0.997 0.003
1 rs23456 3214569 T C 0.869 0.130 0 0.903 0.095 0.002 1 0 0
""".strip())
        self.f.seek(0)

        # Second IMPUTE2 file
        self.f2 = tempfile.NamedTemporaryFile("w")
        self.f2.write("""
1 rs1234567 1234567 A T 1 0 0 0.1 0.3 0.6 0.1 0.35 0.55 0 1 0
""".strip())
        self.f2.seek(0)

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

        self.prob_snp3 = (
            "rs1234567",
            "1",
            1234567,
            "A",
            "T",
            np.array([[1, 0, 0], [0.1, 0.3, 0.6], [0.1, 0.35, 0.55],
                      [0, 1, 0]]),
        )

        self.dosage_snp1 = (
            np.array([0., 0.002, 1.003]), 
            {"minor": "G", "major": "A", "maf": 1 / 6.0, "name": "rs12345",
             "chrom": "1", "pos": 1231415, }
        )

        self.dosage_snp2 = (
            np.array([0.130, 0.099, 0]), 
            {"minor": "C", "major": "T", "maf": 0, "name": "rs23456", 
             "chrom": "1", "pos": 3214569}
        )

        self.dosage_snp3_thresh_0 = (
            np.array([2, 0.5, 0.55, 1]),
            {"minor": "A", "major": "T", "maf": 3 / 8.0, "name": "rs1234567",
             "chrom": "1", "pos": 1234567},
        )

        self.dosage_snp3_thresh_9 = (
            np.array([0, np.nan, np.nan, 1]),
            {"minor": "T", "major": "A", "maf": 1 / 4.0, "name": "rs1234567",
             "chrom": "1", "pos": 1234567},
        )

    def tearDown(self):
        # Closing the temporary file
        self.f.close()
        self.f2.close()

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

    def test_dosage(self):
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

    def test_dosage_maf(self):
        """Test the maf computation with different prob threshold."""
        # Checking for probability threshold of 0
        with fmts.impute2.Impute2File(self.f2.name, "dosage") as f:
            self.assertTrue(compare_dosages(
                self,
                f.readline(),
                self.dosage_snp3_thresh_0,
            ))

        # Checking for probability threshold of 0.9
        with fmts.impute2.Impute2File(self.f2.name, "dosage",
                                      prob_threshold=0.9) as f:
            self.assertTrue(compare_dosages(
                self,
                f.readline(),
                self.dosage_snp3_thresh_9,
            ))
