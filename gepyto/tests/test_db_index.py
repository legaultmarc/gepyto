
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
import string
import random
import os

from .. import db


class TestIndex(unittest.TestCase):
    """Tests the db.index module.

    We make sure that file indexing works.

    """

    @classmethod
    def setUpClass(cls):
        cls.fn = ".test_index_gepyto.txt"
        cls.f = open(cls.fn, "w")

        # The third column is 1 for the first occurence of the locus and
        # gibberish for the other columns.
        cls.positions = [
            ("chr1", 1, 1), (1, 11, 1), (1, 123, 1), (2, 1, 1), (2, 11, 1), 
            (2, 11, "gibberish"), (2, 11, "blablatidoo"), (2, 11, "potatoes"),
            (2, 21, 1), (2, 34, 1), (3, 1, 1), (3, 2, 1), (3, 3, 1), (3, 4, 1),
            (3, 4, "banana"), (3, 5, 1), ("chr3", 6, 1), (3, 7, 1), (3, 8, 1),
            (3, 9, 1), (4, 1, 1), (5, 2, 1), (5, 2, "duplicate"), ("X", 3, 1),
            ("Y", 2, 1), ("Y", 2, 131), ("Y", 2, "duplicate_again"),
        ]
        for fields in cls.positions:
            cls.f.write("\t".join([str(i) for i in fields]) + "\n")

        cls.f.close()
        cls.f = open(cls.fn, "r")

    @classmethod
    def tearDownClass(cls):
        cls.f.close()
        os.remove(cls.fn)
        os.remove(cls.fn + ".gtidx")

    def test_sparse_queries(self):
        db.index.build_index(TestIndex.fn, 0, 1, index_rate=0.9)
        idx = db.index.get_index(TestIndex.fn)

        loci = TestIndex.positions
        random.shuffle(loci)
        for chrom, pos, _ in loci:
            # We only test indexed chromosomes...
            if chrom in idx[0]["chrom_codes"]:
                # Search the entry.
                self.assertTrue(
                    db.index.goto(TestIndex.f, idx, chrom, pos)
                )

        # Negative examples
        for chrom, pos in [(1, 0), ("X", 1), (3, 10), ("Y", 3)]:
            try:
                self.assertFalse(
                    db.index.goto(TestIndex.f, idx, chrom, pos)
                )
            except db.index.ChromosomeNotIndexed:
                # If chromosomes are not indexed, we expect gepyto to tell
                # you without returning True or False.
                pass

    def test_full_queries(self):
        db.index.build_index(TestIndex.fn, 0, 1, index_rate=1)
        idx = db.index.get_index(TestIndex.fn)

        loci = TestIndex.positions
        random.shuffle(loci)
        for chrom, pos, _ in loci:
            # Search the entry.
            self.assertTrue(
                db.index.goto(TestIndex.f, idx, chrom, pos)
            )

        # Negative examples
        for chrom, pos in [(1, 0), ("X", 1), (3, 10), ("Y", 3)]:
            self.assertFalse(
                db.index.goto(TestIndex.f, idx, chrom, pos)
            )

    def test_repeated_regions(self):
        db.index.build_index(TestIndex.fn, 0, 1, index_rate=0.8)
        idx = db.index.get_index(TestIndex.fn)

        loci = TestIndex.positions
        random.shuffle(loci)
        for chrom, pos, _ in loci:
            try:
                db.index.goto(TestIndex.f, idx, chrom, pos)
                # The column should always be one.
                line = TestIndex.f.readline().rstrip().split("\t")
                self.assertEqual(line[0], str(chrom))
                self.assertEqual(line[1], str(pos))
                self.assertEqual(line[2], "1")

            except db.index.ChromosomeNotIndexed:
                pass
