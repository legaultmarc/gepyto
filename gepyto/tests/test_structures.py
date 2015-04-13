
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

import numpy as np
from ..structures import variants, genes, sequences


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

        # Simple variations
        self.snp = variants.SNP(chrom="19", pos=55663495,
                                       rs=self.snp_rs, ref="C", alt="T")
        self.indel = variants.Indel(chrom="19", pos=55663539,
                                           rs=self.indel_rs,
                                           ref="TTC", alt="T")

        # Insertions
        self.insertion_rs = "rs11273285"
        self.insertion = variants.Indel(chrom="2", pos=221032,
                                               rs=self.insertion_rs,
                                               ref="T", alt="TCAGGCACGTGG")

        # Deletions
        self.deletion_rs = "rs376881461"
        self.deletion = variants.Indel(chrom="1", pos=873775,
                                              rs=self.deletion_rs,
                                              ref="CAGAGCCT", alt="C")

    def test_format(self):
        test_format = "{:chr%c:%p-%p}"
        self.assertEqual("chr19:55663495-55663495",
                         test_format.format(self.snp))
        self.assertEqual("chr19:55663539-55663539",
                         test_format.format(self.indel))

        test_format = "{:%c,%p,1,%r/%a}"
        self.assertEqual("19,55663495,1,C/T", test_format.format(self.snp))
        self.assertEqual("19,55663539,1,TTC/T", test_format.format(self.indel))

        test_format = "{:%c%z}"
        self.assertRaises(KeyError, test_format.format, self.snp)
        self.assertRaises(KeyError, test_format.format, self.indel)

        test_format = "{:%c%C}"
        self.assertRaises(KeyError, test_format.format, self.snp)
        self.assertRaises(KeyError, test_format.format, self.indel)

        test_format = "{:%cc\t%i}"
        self.assertEqual("2c\trs11273285", test_format.format(self.insertion))
        self.assertEqual("1c\trs376881461", test_format.format(self.deletion))

    def test_snp(self):
        snp = variants.SNP.from_ensembl_api(self.snp_rs)
        self.assertEqual(snp, [self.snp])

    def test_indel(self):
        indel = variants.Indel.from_ensembl_api(self.indel_rs)
        self.assertEqual(indel, [self.indel])

    def test_insertion(self):
        indel = variants.Indel.from_ensembl_api(self.insertion_rs)
        self.assertEqual(indel, [self.insertion])

    def test_deletion(self):
        indel = variants.ShortVariant.from_ensembl_api(self.deletion_rs)
        self.assertEqual(indel, [self.deletion])

    def test_snp_init_from_str(self):
        snp = variants.SNP.from_str("chr19:55663495_C/T")
        snp[0].rs = self.snp_rs

        self.assertEqual(snp, [self.snp])

    def test_variant_in(self):
        snp_g_in = genes.Gene(build="GRCh37", chrom="19",
                                     start=55653495, end=55673495, xrefs={},
                                     strand=1, transcripts=[])

        self.assertTrue(self.snp.in_gene(snp_g_in))
        self.assertTrue(self.indel.in_gene(snp_g_in))

        snp_g_out = genes.Gene(build="GRCh37", chrom="19",
                                      start=15653495, end=15673495, xrefs={},
                                      strand=1, transcripts=[])

        self.assertFalse(self.snp.in_gene(snp_g_out))
        self.assertFalse(self.indel.in_gene(snp_g_out))

        snp_g_out2 = genes.Gene(build="GRCh37", chrom="18",
                                       start=55653495, end=55673495, xrefs={},
                                       strand=1, transcripts=[])

        self.assertFalse(self.snp.in_gene(snp_g_out2))
        self.assertFalse(self.indel.in_gene(snp_g_out2))

    def test_variant_comparison(self):
        snp1 = variants.SNP("22", 25855459, "rs12345", "g", "a")
        snp2 = variants.SNP("22", 25855459, None, "g", "a")  # Same as 1
        snp3 = variants.SNP("13", 32942179, "rs9567605", "t", "a")
        indel = variants.Indel("13", 32940014, "rs11571729", "g", "gt")

        self.assertEqual(snp1, snp2)
        self.assertFalse(snp1 == snp3)
        self.assertFalse(snp1 == indel)

        self.assertEqual(set([snp1, snp2]), set([snp1]))
        self.assertTrue(len(set([snp1, snp2, snp3])) == 2)
        self.assertTrue(len(set([snp1, snp2, snp3, indel])) == 3)

    def test_df(self):
        snp1 = variants.SNP("22", 25855459, "rs12345", "g", "a")
        snp2 = variants.SNP("13", 32942179, None, "t", "a")
        indel = variants.Indel("13", 32940014, "rs11571729", "g", "gt")

        # Create a dataframe with those variants.
        df = variants.variant_list_to_dataframe([snp1, snp2, indel])
        self.assertEqual(list(df["chrom"].values), ["22", "13", "13"])
        self.assertEqual(
            list(df["pos"].values),
            [25855459, 32942179, 32940014]
        )
        self.assertEqual(
            list(df["rs"].values),
            ["rs12345", np.nan, "rs11571729"]
        )
        self.assertEqual(list(df["ref"].values), ["G", "T", "G"])
        self.assertEqual(list(df["alt"].values), ["A", "A", "GT"])

        # Now try with an extra field.
        snp1._info = {"test": -1}
        snp2._info = {"test": 0}
        indel._info = {"test": 1}
        df = variants.variant_list_to_dataframe([snp1, snp2, indel])
        self.assertEqual(list(df["test"].values), [-1, 0, 1])

        # Now make sure there is a problem on inconsistent fields.
        snp2._info["test2"] = "error"
        self.assertRaises(
            AssertionError,
            variants.variant_list_to_dataframe,
            [snp1, snp2, indel]
        )

class TestGene(unittest.TestCase):
    """Tests for the struct.genes module.

    """

    def setUp(self):
        # Gene is BRCA2
        self.gene_37 = genes.Gene.factory_ensembl_id("ENSG00000139618")
        self.gene_38 = genes.Gene.factory_ensembl_id(
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
        self.dna = ("ATGGCCAAGGAAAGATGCCTGAAAAAGTCCTTTCAAGATAGTCTTGAAGACATAAAG"
                    "AAGCGAATGAAAGAGAAAAGGAATAAAAACTTGGCAGAGATTGGCAAACGCAGGTCT"
                    "TTTATAGCTGCACCATGCCAAATAATCACCAACACTTCTACACTGCTGAAAAATTAC"
                    "CAAGACAACAACAAAATGTTAGTTTTAGCTTTGGAAAATGAAAAATCCAAAGTGAAA"
                    "GAAGCCCAAGATATCATCCTACAGCTGAGAAAAGAATGTTACTATCTCACATGTCAG"
                    "CTATATGCATTGAAAGGAAAACTTACATCACAACAAACAGTAGAACCTGCTCAGAAC"
                    "CAGGAAATATGTTCCTCTGGAATGGACCCCAATAGTGATGACAGCTCCAGAAATTTA"
                    "TTTGTGAAGGATTTACCGCAAATTCCTCTTGAAGAAACTGAACTTCCAGGACAAGGA"
                    "GAATCATTTCAAATAGAAGATCAGATACCTACTATTCCTCAAGACACACTGGGAGTT"
                    "GATTTTGATTCAGGTGAAGCTAAGTCTACTGATAATGTCTTACCTAGAACTGTATCT"
                    "GTTCGTAGCAGTTTAAAGAAACATTGTTAA")

        self.rna = self.dna.replace("T", "U")

        self.protein = ("MAKERCLKKSFQDSLEDIKKRMKEKRNKNLAEIGKRRSFIAAPCQIITNTSTL"
                        "LKNYQDNNKMLVLALENEKSKVKEAQDIILQLRKECYYLTCQLYALKGKLTSQ"
                        "QTVEPAQNQEICSSGMDPNSDDSSRNLFVKDLPQIPLEETELPGQGESFQIED"
                        "QIPTIPQDTLGVDFDSGEAKSTDNVLPRTVSVRSSLKKHC")

    def test_translation(self):
        dna_seq = sequences.Sequence("test_dna", self.dna, "DNA")
        rna_seq = sequences.Sequence("test_rna", self.rna, "RNA")
        pro_seq = sequences.Sequence("test_pro", self.protein, "AA")

        self.assertEqual(dna_seq.translate().seq, pro_seq.seq)
        self.assertEqual(rna_seq.translate().seq, pro_seq.seq)

    def test_gc_content(self):
        seq = sequences.Sequence("test", "ATGC", "DNA")
        self.assertEqual(seq.gc_content(), 0.5)

        seq = sequences.Sequence("test", "TAGTTACTAT", "DNA")
        self.assertEqual(seq.gc_content(), 0.2)

    def test_reverse_complement(self):
        seq = sequences.Sequence("test", "TAGTVTAMCTATK", "DNA")
        expected = "MATAGKTABACTA"
        self.assertEqual(seq.reverse_complement().seq, expected)
