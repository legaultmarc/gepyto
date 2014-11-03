
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

from .. import utils
from ..structures.variants import SNP, Indel


class TestVariantUtils(unittest.TestCase):

    def setUp(self):
        self.region = "chr19:55663495-55663541"
        self.snp = SNP("19", 55663495, "rs111715315", "C", "T")
        self.indel  = Indel("19", 55663540, 55663542, "rs72301544", "TC", "-")
        self.indel2 = Indel("19", 55663541, 55663542, "rs56007758", "CT", "-")

    def test_ensembl_variant_in_region(self):
        """Test the Ensembl region query for variants.

        This should include:

        - The SNP:  rs111715315 (chr19:55663495_C/T)
        - The Indel: rs72301544 (chr19:55663540_TC/-)
        - The Indel: rs56007758 (chr19:55663541_CT/-)

        Note that this could change if variants are added in this small region.
        Also, this test only covers the build GRCh37.

        """
        
        variants = utils.variants.ensembl_variants_in_region(self.region)

        self.assertEqual(variants[0], self.snp)
        self.assertEqual(variants[1], self.indel)
        self.assertEqual(variants[2], self.indel2)

if __name__ == '__main__':
    unittest.main()

