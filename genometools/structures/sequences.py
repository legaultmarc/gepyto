
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
import re

import numpy as np

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
        self.seq = s.lower()
        self.seq_type = seq_type
        self.info = info
        self._annotations = []

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

    def bbc(self, k=10, alphabet=None):
        """Shortcut to base_base_correlation.

        """
        return self.base_base_correlation(k, alphabet)

    def base_base_correlation(self, k=10, alphabet=None):
        """Compute the base base correlation (BBC) for the sequence.

        :param k: k is a parameter of the BBC. Intuitively, it represents
                  the maximum distance to observe correlation between bases.
        :type k: int

        :param alphabet: List of possible characters. This can be used to avoid
                         autodetection of the alphabet in the case where
                         sequences with missing letters are to be compared.
        :type alphabet: iterable

        :returns: A 16 dimensional vector representing the BBC.
        :rtype: :py:class:`np.ndarray`

        A description of the method can be found here:
        http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4272582

        Liu, Zhi-Hua, et al. "Base-Base Correlation a Novel Sequence Feature
        and its Applications." Bioinformatics and Biomedical Engineering, 2007.
        ICBBE 2007. The 1st International Conference on. IEEE, 2007.

        This implementation is generalized for any sequence type.

        """

        s = self.seq
        if k > len(s) - 2:
            raise Exception("Sequence too short to compute BBC with "
                "k={}".format(k))

        if alphabet is None:
            alphabet = set(s)
        alphabet = sorted(list(alphabet))
        alphabet = dict(zip(alphabet, range(len(alphabet))))
        L = len(alphabet)

        # Compute the base probabilities for every character.
        p = np.zeros(L)
        for c in s:
            p[alphabet[c]] += 1
        p /= np.sum(p)
        p.shape = (1, L)

        # Now we need to compute 
        bbc = np.zeros((L, L))
        for l in range(1, k + 2):
            # We need to compute $p_{ij}(l)$ representing the probability of
            # observing the bases i and j separated by l "gaps".
            # We will compute it for all 16 combinations of alleles.
            l_dist_correlations = np.zeros((L, L))
            for i in range(len(s) - l):
                nuc1 = alphabet[s[i]]
                nuc2 = alphabet[s[i + l]]
                l_dist_correlations[nuc1][nuc2] += 1
            l_dist_correlations /= np.sum(l_dist_correlations)

            # We can now compute the D_{ij}(l) which is the deviation from
            # statistical independance.
            # $D_{ij}(l) = p_{ij}(l) - p_i p_j$
            D = l_dist_correlations - np.dot(p.T, p)

            bbc += D + (D ** 2 / 2 * np.dot(p.T ** 2, p ** 2)) + D ** 3

        # We can now flatten the bbc into a 16 feature vector.
        bbc.shape = (1, L * L)

        return bbc


class SequenceAnnotation(object):
    """Annotation of a 1D sequence of characters. 

    Examples of sequence annotations include DNA elements or protein domains.
    These annotations will be supported by the ``visualisation`` module and
    can also be used to generate automatic reports describing the consequence
    of mutations.

    """ 

    _new_anno_int = -1
    types = {}

    def __init__(self, parent, anno_type, start, end, description=""):
        # Initialize the default types if needed.
        cls = self.__class__
        if cls._new_anno_int == -1:
            cls._init_types()

        self.parent = parent
        self.anno_type = anno_type 

        self.start = int(start)
        self.end = int(end)
        if parent is not None:
            assert start < parent.start and end < parent.end
        assert start <= end

        self.description = description

    @classmethod
    def _init_types(cls):
        """Initialize the default types. 
        
        Internally, the annotation types are stored as integers. This binds
        the default types to the class.

        """

        types = ("MODIFIED_RESIDUE", "BINDING_SITE", "CHAIN", "MOTIF",
                 "HELIX", "STRAND")
        for i, t in enumerate(types):
            setattr(cls, t, i)
            cls.types[i] = t

        cls._new_anno_int = i + 1

    @classmethod
    def register_type(cls, type_name):
        """Register a new annotation type.

        :param type_name: The parameter is a name describing the new annotation
                          name to use in your code. The naming convention is
                          that you use CAPITALS_AND_UNDERSCORES.
        :type type_name: str

        After using this method, you can use SequenceAnnotation.MY_TYPE when
        initializing new SequenceAnnotation objects.

        """

        if not re.match(r"^[A-Z0-9_]+$", type_name):
            raise ValueError("Annotation type names need to follow the Python "
                "global constant naming convention (e.g. BINDING_DOMAIN).")

        if hasattr(cls, type_name):
            raise ValueError("Can't register type '{}' as it already "
                "exists.".format(type_name))

        
        setattr(cls, type_name, cls._new_anno_int)
        cls.types[cls._new_anno_int] = type_name
        cls._new_anno_int += 1

    @classmethod
    def unregister_type(cls, type_name):
        """Remove a registered type.

        :param type_name: The type to remove.
        :type type_name: str

        """

        delattr(cls, type_name)

    def __getattribute__(self, key):
        """We override the attribute lookups to translate the types into their human readable form. 
        
        :param key: The field to lookup.
        :type key: str

        """

        parentf = super(SequenceAnnotation, self).__getattribute__

        if key == "anno_type":
            return parentf("types")[parentf("anno_type")]
        else:
            return parentf(key)

    def __setattr__(self, a, v):
        """We override the attribute setting to be consistent with the getting.

        :param a: The attribute to set.
        :type a: str
        :param b: The attribute value.
        :type b: any

        """

        parentf = super(SequenceAnnotation, self).__setattr__

        if a == "anno_type":
            # Let getattr ici change en string.
            try:
                parentf("anno_type", getattr(self, v))
            except AttributeError:
                raise AttributeError("Unknown annotation type ('{}'). Register "
                                     "it first using the `register_type` "
                                     "method.".format(v))
        else:
            parentf(a, v)

    def __repr__(self):
        return "<{} ({}) {}-{}>".format(
            self.__class__.__name__,
            self.anno_type,
            self.start,
            self.end
        )


