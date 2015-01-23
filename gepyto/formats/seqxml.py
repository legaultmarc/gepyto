#
# Implementation of the SeqXML format into Python objects with
# some elementary operations.
# See http://orthoxml.org/xml/Main.html for more information on the format.
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


import collections
import gzip
import xml.etree.ElementTree as etree

from ..structures.sequences import Sequence


class SeqXML(object):
    """Parses the SeqXML format representing sequence data.

    :param fn: The filename of the SeqXML file. The format description is
               available at
               `orthoxml.org <http://seqxml.org/0.4/seqxml_doc_v0.4.html>`_
               (visited Nov. 2014).
    :type fn: str

    The returned object will have a list of entries which are
    :py:class:`Sequence` objects.

    """

    seq_xml_seqtypes = {
        "DNAseq": "DNA",
        "RNAseq": "RNA",
        "AAseq": "AA",
    }

    def __init__(self, fn):

        opener = gzip.open if fn.endswith(".gz") else open
        with opener(fn) as f:
            tree = etree.parse(f)

        self.root = tree.getroot()

        self.entries = []
        self.id_index = {}

        for entry in self.root:
            # Mandatory fields
            uid = entry.attrib.get("id")
            seq = None
            seq_type = None

            # Additional information
            info = {}

            # Parse the sequence entry.
            seq_types = set(SeqXML.seq_xml_seqtypes.keys())
            for elem in entry:
                # This is the biological sequence.
                if elem.tag in seq_types:
                    seq_type = SeqXML.seq_xml_seqtypes[elem.tag]
                    seq = elem.text
                # Those are all "info" fields.
                elif elem.tag == "property":
                    info[elem.attrib["name"]] = elem.attrib.get("value", 1)
                elif elem.tag == "species":
                    info["species"] = elem.attrib["name"]
                    info["species_ncbi_tax_id"] = elem.attrib["ncbiTaxID"]
                elif elem.tag == "description":
                    info["description"] = elem.text
                elif elem.tag == "DBRef":
                    info["db_name"] = elem.attrib["source"]
                    info["db_acc"] = elem.attrib["id"]

            # Create the Sequence object.
            seq = Sequence(uid, seq, seq_type, info)
            self.entries.append(seq)

        # Build the id_index that allows fast Sequence lookup by id.
        for entry in self.entries:
            self.id_index[entry.uid] = entry

    def get_seq(self, uid):
        """Get a sequence from it's unique identifier.

        :param uid: The sequence id.
        :type uid: str

        """

        if uid not in self.id_index:
            return None
        else:
            return self.id_index[uid]
