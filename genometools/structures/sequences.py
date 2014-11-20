
# Structures to handle biological sequences.
#
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

__all__ = ["Sequence", ]

import textwrap

class Sequence(object):
    """Object to represent biological sequences.

    :param uid: The identifier for this sequence.
    :type uid: str

    :param s: The actual sequence.
    :type s: str

    :param seq_type: The sequence type (DNA, RNA or AA).
    :type seq_type: str

    :param info: A python dict of extra parameters (optional).
    :type info: dict

    Common examples for the info attributes:

    * ``species``: Homo sapiens
    * ``species_ncbi_tax_id``: 9606
    * ``description``: dystroglycan 1
    * ``db_name``: RefSeq
    * ``db_acc``: NM_004393

    """

    types = set(["DNA", "RNA", "AA"])

    def __init__(self, uid, s, seq_type, info=None):
        if seq_type not in Sequence.types:
            raise ValueError("Invalid sequence type {}. Allowed types are: "
                "{}".format(seq_type, list(Sequence.types)))

        self.uid = uid
        self.seq = s
        self.seq_type = seq_type
        self.info = info

    def __repr__(self):
        return "<Sequence: {}>".format(self.uid)

    def to_fasta(self, line_len=80, full_header=False):
        """Converts the sequence to a valid fasta string.

        :param line_len: The maximum line length for the sequence.
        :type line_len: int

        :param full_header: Add the contents of the info field to the header.
                            (default: False).
        :type full_header: bool

        :returns: A fasta string.
        :rtype: str

        """

        s = "> {}".format(self.uid)
        if full_header:
            s += " - "
            s += ", ".join(
                ["'{}'='{}'".format(k, v) for k, v in self.info.items()]
            )
        s = [s, ] + textwrap.wrap(self.seq, line_len)
        return "\n".join(s) + "\n"
