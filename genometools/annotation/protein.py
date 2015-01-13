
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


import logging
import contextlib
import xml.etree.ElementTree as etree

from .. import structures as struct

NS = "{http://uniprot.org/uniprot}"

try:
    # Python 2 support
    from urllib2 import urlopen, HTTPError
except ImportError:
    # Python 3 support
    from urllib.request import urlopen, HTTPError


def annotate_protein_sequence(s, uniprot_id=None):
    """Query uniprot to get protein information.

    :param s: A protein sequence.
    :type s: :py:class:`genometools.structures.sequences.Sequence`

    :param uniprot_id: Optional paremeter specifying a UniProt accession. If
                       no value is provided, the sequence's ``uid`` field is
                       used.
    :type uniprot_id: str

    """

    assert s.seq_type == "AA"

    if uniprot_id is None:
        uniprot_id = s.uid

    url = "http://www.uniprot.org/uniprot/{}.xml".format(uniprot_id)

    with contextlib.closing(urlopen(url)) as stream:
        tree = etree.parse(stream)
        root = tree.getroot()
        entry = root.find(NS + "entry")

    # We have an entry, we will parse the feature tags.
    for feature in entry.findall(NS + "feature"):
        anno_feat = _parse_uniprot_feature(feature)
        if anno_feat is not None:
            anno_feat.parent = s
            s._annotations.append(anno_feat)

    return s._annotations


def _parse_uniprot_feature(node):
    """Takes a feature node from UniProt and returns a corresponding SequenceAnnotation.

    :param node: A node representing a ``feature``
    :type node: xml.etree.ElementTree.Element

    :returns: A feature annotation that will be bound to the Sequence object.
    :rtype: TODO

    """

    # SEQUENCE_CONFLICT is a uniprot specific annotation.
    if not hasattr(struct.sequences.SequenceAnnotation, "SEQUENCE_CONFLICT"):
        struct.sequences.SequenceAnnotation.register_type(
            "SEQUENCE_CONFLICT"
        )

    type_map = {
        "modified residue": "MODIFIED_RESIDUE",
        "binding site": "BINDING_SITE",
        "chain": "CHAIN",
        "short sequence motif": "MOTIF",
        "helix": "HELIX",
        "strand": "STRAND",
        "sequence variant": "SEQUENCE_VARIANT",
        "sequence conflict": "SEQUENCE_CONFLICT",
        "turn": "TURN",
        "mutagenesis site": "MUTAGENESIS_SITE",
        "mutagenesis site": "MUTAGENESIS_SITE",
        "repeat": "REPEAT",
        "region of interest": "REGION_OF_INTEREST",
    }

    # Parse annotation type
    anno_type = node.attrib["type"]
    if anno_type not in type_map:
        logging.warning("Ignoring unknown features '{}'".format(anno_type))
        return None

    anno_type = type_map[anno_type]

    # Parse description
    desc = node.attrib.get("description")

    # Parse location
    loc = node.find(NS + "location")
    if loc.find(NS + "position") is not None:
        start = end = int(loc.find(NS + "position").attrib["position"])
    else:
        start = int(loc.find(NS + "begin").attrib["position"])
        end = int(loc.find(NS + "end").attrib["position"])

    # Build the object.
    return struct.sequences.SequenceAnnotation(
        None, # Not our responsibility to bind the parent.
        anno_type,
        start,
        end, 
        desc,
    )

