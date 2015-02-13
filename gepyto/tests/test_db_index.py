
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

        cls.positions = [
            ("chr1", 1), (1, 11), (1, 123), (2, 1), (2, 11), (2, 11), (2, 11),
            (2, 11), (2, 21), (2, 34), (3, 1), (3, 2), (3, 3), (3, 4), (3, 4),
            (3, 5), ("chr3", 6), (3, 7), (3, 8), (3, 9), (4, 1), (5, 2),
            ("X", 3), ("Y", 2),
        ]
        for chrom, pos in cls.positions:
            # Generate a random third column.
            n = random.randint(100, 125)
            s = "".join(
                [random.choice(string.ascii_lowercase) for i in range(n)]
            )
            cls.f.write("\t".join((str(chrom), str(pos), s)) + "\n")

        cls.f.close()
        cls.f = open(cls.fn, "r")

    @classmethod
    def tearDownClass(cls):
        cls.f.close()
        os.remove(cls.fn)
        os.remove(cls.fn + ".gtidx")

    def setUp(self):
        self.idx_info = db.index.build_index(TestIndex.fn, 0, 1, index_rate=1)

    def test_query(self):
        idx = db.index.get_index(TestIndex.fn)

        loci = TestIndex.positions
        random.shuffle(loci)
        for chrom, pos in loci:
            if type(chrom) is str:
                chrom = chrom.lstrip("chr")
            # Search the entry.
            self.assertTrue(
                db.index.goto(TestIndex.f, idx, chrom, pos)
            )

        # Negative examples
        for chrom, pos in [(1, 0), (5, 1), (3, 10)]:
            self.assertFalse(
                db.index.goto(TestIndex.f, idx, chrom, pos)
            )
