# Data provider abstract classes and interfaces.
#
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

"""
Interfaces and abstract classes for data formats supporting similar
functionality.
"""

class AbstractSequenceDictionary(object):
    """Abstract data structure representing a collection of named sequences.

    An example of an implementation for this class could be a fasta file
    parser.

    """
    def __init__(self):
        raise NotImplementedError()

    def __getitem__(self, key):
        """Get sequences or substrings.

        Implementations should support double indexing of the form:  ::

            sequence_dict["chr1"][pos]
            sequence_dict["chr1"][start:end]

        """
        raise NotImplementedError()

    def get(self, key):
        """Get a given sequence from the container.

        As opposed to regular indexing, this should return None if the sequence
        can't be found.
        """
        raise NotImplementedError()

    def keys(self):
        """Get a list of the sequence IDs in the dictionary."""
        raise NotImplementedError()
