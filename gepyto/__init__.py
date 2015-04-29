
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
__credits__ = ["Marc-Andre Legault", "Louis-Philippe Lemieux Perreault"]
__maintainer__ = "Marc-Andre Legault"
__email__ = "legaultmarc@gmail.com"
__status__ = "Development"


# Loading the version
try:
    from .version import gepyto_version as __version__

except ImportError:
    __version__ = None


def test(verbosity=1):
    import logging
    import unittest
    from .tests import test_suite

    logging.disable(logging.CRITICAL)

    unittest.TextTestRunner(verbosity=verbosity).run(test_suite)

    logging.disable(logging.NOTSET)
