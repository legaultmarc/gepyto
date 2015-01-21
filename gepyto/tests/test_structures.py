
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

import unittest

from .. import structures as struct

class TestVariant(unittest.TestCase):
    """Tests the struct.variants module.

    We compare the initialized values with what is expected from manual 
    queries in the database for the ensembl init. function.

    If this changes, or if Ensembl drops GRCh37 support. These will have to
    be rewritten.

    We also test the other small methods.

    """

    def setUp(self):
        self.snp_rs = "rs111715315"
        self.indel_rs = "rs72301544"

        self.snp = struct.variants.SNP("19", 55663495, self.snp_rs, "C", "T")
        self.indel  = struct.variants.Indel("19", 55663540, 55663542, 
            self.indel_rs, "TC", "-")

    def test_snp(self):
        snp = struct.variants.Variant.from_ensembl_api(self.snp_rs)
        self.assertEqual(snp, self.snp)

    def test_indel(self):
        indel = struct.variants.Variant.from_ensembl_api(self.indel_rs)
        self.assertEqual(indel, self.indel)

    def test_snp_init_from_str(self):
        snp = struct.variants.SNP.from_str("chr19:55663495_C/T")
        snp.rs = self.snp_rs

        self.assertEqual(snp, self.snp)

    def test_snp_get_position(self):
        self.assertEqual(self.snp.get_position(), "chr19:55663495")

    def test_indel_get_position(self):
        self.assertEqual(self.indel.get_position(), "chr19:55663540-55663542")

    def test_variant_in(self):
        snp_g_in = struct.genes.Gene(build="GRCh37", chrom="19", 
            start=55653495, end=55673495, xrefs={}, strand=1, transcripts=[])

        self.assertTrue(self.snp.in_gene(snp_g_in))
        self.assertTrue(self.indel.in_gene(snp_g_in))

        snp_g_out = struct.genes.Gene(build="GRCh37", chrom="19", 
            start=15653495, end=15673495, xrefs={}, strand=1, transcripts=[])

        self.assertFalse(self.snp.in_gene(snp_g_out))
        self.assertFalse(self.indel.in_gene(snp_g_out))

        snp_g_out2 = struct.genes.Gene(build="GRCh37", chrom="18", 
            start=55653495, end=55673495, xrefs={}, strand=1, transcripts=[])

        self.assertFalse(self.snp.in_gene(snp_g_out2))
        self.assertFalse(self.indel.in_gene(snp_g_out2))


class TestGene(unittest.TestCase):
    """Tests for the struct.genes module.

    """

    def setUp(self):
        # Gene is BRCA2
        self.gene_37 = struct.genes.Gene.factory_ensembl_id("ENSG00000139618")
        self.gene_38 = struct.genes.Gene.factory_ensembl_id(
            "ENSG00000139618", 
            build="GRCh38"
        )

    def test_factory_ensembl_id(self):


        # Chrom
        self.assertEqual(self.gene_37.chrom, "13")
        self.assertEqual(self.gene_38.chrom, "13")

        # Region
        self.assertEqual(self.gene_37.start, 32889611)
        self.assertEqual(self.gene_37.end, 32973805)

        self.assertEqual(self.gene_38.start, 32315474)
        self.assertEqual(self.gene_38.end, 32400266)

        # Has some transcripts...
        self.assertTrue(len(self.gene_37.transcripts) > 1)
        self.assertTrue(len(self.gene_38.transcripts) > 1)

        # Symbol
        self.assertEqual(self.gene_37.symbol, "BRCA2")
        self.assertEqual(self.gene_38.symbol, "BRCA2")

        # Has some exons.
        self.assertTrue(len(self.gene_37.exons) > 1)
        self.assertTrue(len(self.gene_38.exons) > 1)


class TestSequence(unittest.TestCase):
    """Tests for the struct.sequences module.

    """

    def setUp(self):
        self.dna = ("atggccaaggaaagatgcctgaaaaagtcctttcaagatagtcttgaagacataaag"
                    "aagcgaatgaaagagaaaaggaataaaaacttggcagagattggcaaacgcaggtct"
                    "tttatagctgcaccatgccaaataatcaccaacacttctacactgctgaaaaattac"
                    "caagacaacaacaaaatgttagttttagctttggaaaatgaaaaatccaaagtgaaa"
                    "gaagcccaagatatcatcctacagctgagaaaagaatgttactatctcacatgtcag"
                    "ctatatgcattgaaaggaaaacttacatcacaacaaacagtagaacctgctcagaac"
                    "caggaaatatgttcctctggaatggaccccaatagtgatgacagctccagaaattta"
                    "tttgtgaaggatttaccgcaaattcctcttgaagaaactgaacttccaggacaagga"
                    "gaatcatttcaaatagaagatcagatacctactattcctcaagacacactgggagtt"
                    "gattttgattcaggtgaagctaagtctactgataatgtcttacctagaactgtatct"
                    "gttcgtagcagtttaaagaaacattgttaa")

        self.rna = self.dna.replace("t", "u")

        self.protein = ("MAKERCLKKSFQDSLEDIKKRMKEKRNKNLAEIGKRRSFIAAPCQIITNTSTL"
                        "LKNYQDNNKMLVLALENEKSKVKEAQDIILQLRKECYYLTCQLYALKGKLTSQ"
                        "QTVEPAQNQEICSSGMDPNSDDSSRNLFVKDLPQIPLEETELPGQGESFQIED"
                        "QIPTIPQDTLGVDFDSGEAKSTDNVLPRTVSVRSSLKKHC").lower()

    def test_translation(self):
        dna_seq = struct.sequences.Sequence("test_dna", self.dna, "DNA")
        rna_seq = struct.sequences.Sequence("test_rna", self.rna, "RNA")
        pro_seq = struct.sequences.Sequence("test_pro", self.protein, "AA")

        self.assertEqual(dna_seq.translate().seq, pro_seq.seq)
        self.assertEqual(rna_seq.translate().seq, pro_seq.seq)

