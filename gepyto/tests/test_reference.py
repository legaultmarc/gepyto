
# This file is part of gepyto.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = ("Copyright 2014 Marc-Andre Legault and Louis-Philippe "
                 "Lemieux Perreault. All rights reserved.")
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


import unittest

from ..reference import Reference, InvalidMapping


class TestReference(unittest.TestCase):

    def setUp(self):
        """Sets up the test."""
        self.remote_ref = Reference(remote=True)

    def test_check_variant_reference(self):
        """Tests the 'check_variant_reference' function."""
        pass

    def test_get_nucleotide(self):
        """Tests the 'test_get_nucleotide' function."""
        self.assertEqual(self.remote_ref.get_nucleotide(22, 25855459), "G")
        self.assertEqual(self.remote_ref.get_nucleotide(22, 25855460), "A")

    def test_get_sequence(self):
        """Tests the 'get_sequence' function."""
        self.assertEqual(
            self.remote_ref.get_sequence(22, 25855459, length=5),
            "GACTT"
        )

        self.assertEqual(
            self.remote_ref.get_sequence(22, 25855459, 25855463),
            "GACTT"
        )

        self.assertRaises(
            InvalidMapping,
            self.remote_ref.get_sequence,
            30, 1234,
            length=3
        )
