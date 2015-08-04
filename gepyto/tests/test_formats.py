
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


import os
import datetime
import unittest
import tempfile

import numpy as np

from ..formats import impute2
from .. import formats as fmts
from ..structures.sequences import Sequence


def compare_vectors(v1, v2):
    try:
        np.testing.assert_array_almost_equal(v1, v2)
    except AssertionError:
        return False
    return True


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
        if k == "maf":
            self.assertAlmostEqual(info1[k], info2[k])
        elif k == "minor_allele_count":
            self.assertEqual(round(info1[k]), round(info2[k]))
        else:
            self.assertEqual(info1[k], info2[k])

    if v1.dtype is np.dtype(float):
        # Checking the dosage values (nan != nan in numpy)
        return compare_vectors(v1, v2)

    if v1.dtype.char == "U" or v1.dtype.char == "S":
        return (v1 == v2).all()


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
1 rs23456 3214569 T C 0.869 0.130 0 0.903 0.095 0.002 0 0 1
1 rs23457 3214570 T TC 0.869 0.130 0 0 1 0 0 0 1
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
            np.array([[0.869, 0.130, 0], [0.903, 0.095, 0.002], [0, 0, 1]])
            # TT, TT, CC
        )

        self.prob_indel = (
            "rs23457",
            "1",
            3214570,
            "T",
            "TC",
            np.array([[0.869, 0.130, 0], [0, 1, 0], [0, 0, 1]]),
            # T/T, T/TC, TC/TC
        )

        self.prob_snp3 = (
            "rs1234567",
            "1",
            1234567,
            "A",
            "T",
            np.array([[1, 0, 0], [0.1, 0.3, 0.6], [0.1, 0.35, 0.55],
                      [0, 1, 0]]),
            # AA, TT, TT, AT
            # 0, 1.5, 1.45, 1
        )

        self.dosage_snp1 = (
            np.array([0., 0.002, 1.003]),
            {"minor": "G", "major": "A", "maf": 1.005 / 6.0, "name": "rs12345",
             "chrom": "1", "pos": 1231415, "minor_allele_count": 1.005}
        )

        self.dosage_snp2 = (
            np.array([0.130, 0.099, 2]),
            {"minor": "C", "major": "T", "maf": 2.229 / 6.0, "name": "rs23456",
             "chrom": "1", "pos": 3214569, "minor_allele_count": 2.229}
        )

        self.dosage_indel = (
            np.array([1.87, 1, 0]),
            {"minor": "T", "major": "TC", "maf": 2.87 / 6.0, "name": "rs23457",
             "chrom": "1", "pos": 3214570, "minor_allele_count": 2.87}
        )

        self.dosage_snp3_thresh_0 = (
            np.array([0, 1.5, 1.45, 1]),
            {"minor": "T", "major": "A", "maf": 3.95 / 8.0,
             "name": "rs1234567", "chrom": "1", "pos": 1234567,
             "minor_allele_count": 3.95},
        )

        self.dosage_snp3_thresh_9 = (
            np.array([0, np.nan, np.nan, 1]),
            {"minor": "T", "major": "A", "maf": 1 / 4.0, "name": "rs1234567",
             "chrom": "1", "pos": 1234567, "minor_allele_count": 1},
        )

        self.hard_call_snp1 = (
            np.array(["A A", "A A", "A G"]),
            {"name": "rs12345", "chrom": "1", "pos": 1231415},
        )

        self.hard_call_snp2 = (
            np.array(["T T", "T T", "C C"]),
            {"name": "rs23456", "chrom": "1", "pos": 3214569},
        )

        self.hard_call_indel = (
            np.array(["T T", "T TC", "TC TC"]),
            {"name": "rs23457", "chrom": "1", "pos": 3214570},
        )

        self.hard_call_snp2_thresh_9 = (
            np.array(["0 0", "T T", "C C"]),
            {"name": "rs23456", "chrom": "1", "pos": 3214569},
        )

        self.hard_call_indel_thresh_9 = (
            np.array(["0 0", "T TC", "TC TC"]),
            {"name": "rs23457", "chrom": "1", "pos": 3214570},
        )

    def tearDown(self):
        # Closing the temporary file
        self.f.close()
        self.f2.close()

    def test_syntax(self):
        """This is mostly to test the syntax for object initialization. """
        # Reading the three lines
        f = impute2.Impute2File(self.f.name)
        lines = [f.readline(), ]
        lines.append(f.readline())
        lines.append(f.readline())
        f.close()

        # Test the context manager and iterator compliance.
        with impute2.Impute2File(self.f.name) as f:
            for i, l in enumerate(f):
                self.assertTrue(compare_lines(l, lines[i]))

        # Make sure using invalid arguments raises error.
        self.assertRaises(
            TypeError,
            impute2.Impute2File,
            self.f.name,
            prob_threshold=0.8
        )

        # But if we're in dosage mode, this should work.
        f = impute2.Impute2File(
            self.f.name, "dosage", prob_threshold=0.8
        )
        f.close()

    def test_proba_reader(self):
        """This is to read the probabilities matrix. """
        with impute2.Impute2File(self.f.name) as f:
            for i, line in enumerate(f):
                if i == 0:
                    self.assertTrue(compare_lines(line, self.prob_snp1))
                elif i == 1:
                    self.assertTrue(compare_lines(line, self.prob_snp2))
                elif i == 2:
                    self.assertTrue(compare_lines(line, self.prob_indel))
                else:
                    raise Exception()

    def test_dosage(self):
        """This is to read as dosage vectors. """
        with impute2.Impute2File(self.f.name, "dosage") as f:
            for i, line in enumerate(f):
                if i == 0:
                    self.assertTrue(compare_dosages(
                        self,
                        line,
                        self.dosage_snp1,
                    ))

                elif i == 1:
                    self.assertTrue(compare_dosages(
                        self,
                        line,
                        self.dosage_snp2,
                    ))

                elif i == 2:
                    self.assertTrue(compare_dosages(
                        self,
                        line,
                        self.dosage_indel,
                    ))

                else:
                    raise Exception()

    def test_dosage_maf(self):
        """Test the maf computation with different prob threshold."""
        # Checking for probability threshold of 0
        with impute2.Impute2File(self.f2.name, "dosage") as f:
            self.assertTrue(compare_dosages(
                self,
                f.readline(),
                self.dosage_snp3_thresh_0,
            ))

        # Checking for probability threshold of 0.9
        with impute2.Impute2File(self.f2.name, "dosage",
                                 prob_threshold=0.9) as f:
            self.assertTrue(compare_dosages(
                self,
                f.readline(),
                self.dosage_snp3_thresh_9,
            ))

    def test_hard_call(self):
        """Test the hard calling of imputed markers."""
        with impute2.Impute2File(self.f.name, "hard_call") as f:
            for i, line in enumerate(f):
                if i == 0:
                    self.assertTrue(compare_dosages(
                        self,
                        line,
                        self.hard_call_snp1,
                    ))

                elif i == 1:
                    self.assertTrue(compare_dosages(
                        self,
                        line,
                        self.hard_call_snp2,
                    ))

                elif i == 2:
                    self.assertTrue(compare_dosages(
                        self,
                        line,
                        self.hard_call_indel,
                    ))

                else:
                    raise Exception()

        with impute2.Impute2File(self.f.name, "hard_call",
                                 prob_threshold=0.9) as f:
            for i, line in enumerate(f):
                if i == 0:
                    self.assertTrue(compare_dosages(
                        self,
                        line,
                        self.hard_call_snp1
                    ))

                elif i == 1:
                    self.assertTrue(compare_dosages(
                        self,
                        line,
                        self.hard_call_snp2_thresh_9,
                    ))

                elif i == 2:
                    self.assertTrue(compare_dosages(
                        self,
                        line,
                        self.hard_call_indel_thresh_9,
                    ))

                else:
                    raise Exception()


class TestGTF(unittest.TestCase):
    """Test the GTF file parser."""
    def setUp(self):
        self.cls = fmts.gtf.GTFFile

    def test_real_life_body(self):
        ans = [
            ("O60503", "UniProtKB", "Chain", 1, 1353, None, None, None,
             {"ID": "PRO_0000195708", "Note": "Adenylate cyclase type 9"}),
            ("O60503", "UniProtKB", "Topological domain", 1, 117, None, None,
             None, {"Note": "Cytoplasmic", "evidence": "ECO:0000255"}),
            ("O60503", "UniProtKB", "Transmembrane", 118, 138, None, None,
             None, {"Note": "Helical", "evidence": "ECO:0000255"}),
        ]

        # Check line parsing.
        with BasicGTF() as f:
            gtf = self.cls(f.name)
            for i, line in enumerate(gtf):
                for j in range(len(line)):
                    self.assertEqual(ans[i][j], line[j])
            gtf.close()

    def test_real_life_meta(self):
        with BasicGTF() as f:
            gtf = self.cls(f.name)
            self.assertEqual(gtf.gff_version, 3)
            self.assertEqual(gtf.sequence_regions["O60503"], (1, 1353))
        gtf.close()

    def test_readline(self):
        lines = []
        with BasicGTF() as f:
            gtf = self.cls(f.name)
            lines = list(gtf)
            gtf.close()

        with BasicGTF() as f:
            gtf = self.cls(f.name)
            for i in range(len(lines)):
                self.assertEqual(lines[i], gtf.readline())
            gtf.close()

    def test_context_mgr(self):
        with BasicGTF() as f:
            with self.cls(f.name) as gtf:
                for i, line in enumerate(gtf):
                    pass
            self.assertEqual(i, 2)

    def test_headers(self):
        f = tempfile.NamedTemporaryFile("w")
        f.write("""
#Nonparsed comment here.
##gff-version 3
#Nonparsed comment here.
##source-version myProgram v1.2
##date 2015-12-29
##Type DNA my_dna_sequence
##Type RNA my_rna_sequence
##Type Protein my_protein_sequence
#A random comment
##DNA myseq
##acggctcggattggcgctggatgatagatcagacgac
##tccccgcaaactcgggcaggttgatttattagaa
##end-DNA
##Protein myprot
##MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSF
##end-Protein
#A comment: done with that.
##sequence-region my_dna_sequence 2341 14510
""".lstrip())
        f.seek(0)

        seq1 = ("acggctcggattggcgctggatgatagatcagacgactccccgcaaactcgggcaggttg"
                "atttattagaa")

        seq2 = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSF"

        with self.cls(f.name) as gtf:
            self.assertEqual(gtf.gff_version, 3)
            self.assertEqual(gtf.source_version["myProgram"], "v1.2")
            self.assertEqual(
                gtf.date,
                datetime.datetime.strptime("2015-12-29", "%Y-%m-%d")
            )
            self.assertEqual(gtf.type["my_dna_sequence"], "DNA")
            self.assertEqual(gtf.type["my_rna_sequence"], "RNA")
            self.assertEqual(gtf.type["my_protein_sequence"], "Protein")
            self.assertEqual(
                gtf.sequences["myseq"],
                Sequence("myseq", seq1, "DNA")
            )
            self.assertEqual(
                gtf.sequences["myprot"],
                Sequence("myprot", seq2, "AA")
            )

        gtf.close()
        f.close()


class BasicGTF(object):
    def __init__(self):
        self.f = tempfile.NamedTemporaryFile("w", delete=False)
        self.f.write("""
##gff-version 3
##sequence-region O60503 1 1353
O60503	UniProtKB	Chain	1	1353	.	.	.	ID=PRO_0000195708;Note=Adenylate cyclase type 9
O60503	UniProtKB	Topological domain	1	117	.	.	.	Note=Cytoplasmic;evidence=ECO:0000255
O60503	UniProtKB	Transmembrane	118	138	.	.	.	Note=Helical;evidence=ECO:0000255
""".lstrip())
        self.f.seek(0)

    def __enter__(self):
        return self.f

    def __exit__(self, *args):
        self.f.close()
        os.remove(self.f.name)


class TestGFF(TestGTF):
    """Test the GFF binding."""
    def setUp(self):
        self.cls = fmts.gff.GFFFile
