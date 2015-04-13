
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

from ..utils import variants
from ..structures.variants import SNP, Indel


class TestVariantUtils(unittest.TestCase):

    def setUp(self):
        self.region = "chr19:55663495-55663541"
        self.snp = SNP("19", 55663495, "rs111715315", "C", "T")
        self.indel1 = Indel("19", 55663539, "rs72301544", "TTC", "T")
        self.indel2 = Indel("19", 55663540, "rs56007758", "TCT", "T")

    def test_ensembl_variant_in_region(self):
        """Test the Ensembl region query for variants.

        This should include:

        - The SNP:  rs111715315 (chr19:55663495_C/T)
        - The Indel: rs72301544 (chr19:55663540_TC/-)
        - The Indel: rs56007758 (chr19:55663541_CT/-)

        Note that this could change if variants are added in this small region.
        Also, this test only covers the build GRCh37.

        """

        region_vars = variants.ensembl_variants_in_region(self.region)

        # Testing that some variants are in there is more robust to changes
        # on the Ensembl side (i.e. different builds).
        self.assertTrue(self.snp in region_vars)
        self.assertTrue(self.indel1 in region_vars)
        self.assertTrue(self.indel2 in region_vars)

if __name__ == '__main__':
    unittest.main()
