# Fasta file format parser and utilities.
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


import pyfaidx
from .utils import AbstractSequenceDictionary


class FastaFile(AbstractSequenceDictionary):
    def __init__(self, filename):
        if filename.endswith(".gz"):
            raise ValueError("Support for compressed fasta files is not "
                             "yet implemented.")

        self.f = pyfaidx.Fasta(filename)

    def __getitem__(self, key):
        return self.f[key]

    def keys(self):
        return self.f.keys()

    def get(self, key):
        seq = self.f.get(key)
        if seq:
            return seq.seq
        return None

    def close(self):
        self.f.close()
