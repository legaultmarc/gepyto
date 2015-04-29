
# Interface to the human genome reference. This is useful to quickly extract
# sequences or to make sure the reference alleles are ok for variants.

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


import os
import functools
import re
from collections import namedtuple

from pyfaidx import Fasta

from . import settings
from .db.ensembl import query_ensembl, get_url_prefix


class _RemoteChromosome(object):
    def __init__(self, url):
        self.url = url

    def __getitem__(self, key):
        # This also makes sure that we used 1-based indexing.
        if type(key) is int:
            key += 1
            return self.do_query(key, key + 1)
        elif type(key) is slice:
            return self.do_query(key.start + 1, key.stop + 1)

    def do_query(self, start, end):
        start, end = sorted([start, end])
        res = query_ensembl(
            self.url.format(start=start, end=end - 1)
        )

        if res is None:
            raise InvalidMapping("Invalid remote region query.")

        # Note the "comp" field is ignored (see pyfaidx.Sequence).
        seq_obj = namedtuple("Sequence", ["name", "seq", "start", "end"])
        return seq_obj(res["id"], res["seq"], start, end)


class _RemoteReference(object):
    """Imitates the pyfaidx.Fasta object to query the human genome reference.

    This works by using the Ensembl REST API to do sequence region queries.
    It should be able to transparently replace the :py:class:`pyfaidx.Fasta`
    class.

    .. note::

        This is a fairly low level class. Users should be able to fill most
        of their needs using the :py:class:`Reference` object.

    """
    def __init__(self, ref):
        self.url = "{}sequence/region/homo_sapiens/".format(
            get_url_prefix(ref)
        )
        self.url += "{chrom}"

    def __getitem__(self, key):
        """This is to implement the ref[] behaviour.

        It assumes that the argument will be a chromosome.

        """
        if not re.match(settings.CHROM_REGEX, key):
            raise KeyError("Invalid chromosome '{}'.".format(key))

        url = (self.url.format(chrom=key) +
               ":{start}-{end}?content-type=application/json")

        return _RemoteChromosome(url)

    def get(self, key):
        try:
            return self.__getitem__(key)
        except KeyError:
            return None

    def close(self):
        pass


class Reference(object):
    """Interface to the human genome reference file.

    This class uses ``pyfaidx`` to parse the genome reference file referenced
    by ``settings.REFERENCE_PATH``.

    This can only be a single plain fasta file.

    Also note that if the path is not in the ``~/.gtconfig/gtrc.ini`` file,
    gepyto will look for an environment variable named ``REFERENCE_PATH``.

    If the genome file can't be found, this class fallbacks to the Ensembl
    remote API to get the sequences.

    This behaviour can also be forced by using the ``remote=True`` argument.

    """
    def __init__(self, remote=False):
        if not remote:
            try:
                self.ref = Fasta(settings.REFERENCE_PATH)
            except IOError:
                self.ref = _RemoteReference(settings.BUILD)
        else:
            self.ref = _RemoteReference(settings.BUILD)

        # Add a get method. This will not be sensitive to "chr" prefixes.
        def get(fasta, chrom):
            chr_prefix = chrom.startswith("chr")
            try:
                return fasta[chrom]
            except KeyError:
                pass
            try:
                # If there was a prefix, we try without.
                if chr_prefix:
                    return fasta[chrom[3:]]
                # If there was no prefix, we try with.
                else:
                    return fasta["chr{}".format(chrom)]
            except KeyError:
                # If it is a true mismatch, we return None.
                return None

        self.ref.get = functools.partial(get, self.ref)

    def check_variant_reference(self, variant, flip=False):
        """Given a variant, makes sure that the 'ref' allele is consistent
        with the human genome reference.

        :param variant: The variant to verify.
        :type variant: :py:class:`gepyto.structures.variants.Variant` subclass

        :param flip: If ``True`` incorrect ``(ref, alt)`` pairs will be
                     flipped (Default: False).
        :type flip: bool

        :returns: If flip is True, it returns the correct variant or raises
                  a ``ValueError`` in case it is not salvageable. If flip
                  is False, a bool is simply returned.

        """

        type_message = ("Unsupported argument to check_variant_reference. "
                        "A SNP or Indel object has to be provided.")

        if not (hasattr(variant, "chrom") and
                hasattr(variant, "pos") and
                hasattr(variant, "ref") and
                hasattr(variant, "alt")):
            raise TypeError(type_message)

        if (len(variant.ref) == len(variant.alt) == 1 and
            "-" not in (variant.ref + variant.alt)):
            return check_snp_reference(variant, self, flip)
        else:
            return check_indel_reference(variant, self, flip)

    def get_nucleotide(self, chrom, pos):
        """Get the nucleotide at the given genomic position. """
        return self.get_sequence(str(chrom), pos, length=1)

    def get_sequence(self, chrom, start, end=None, length=None):
        """Get the nucleotide sequence at the given genomic locus.

        :param chrom: The chromosome.
        :type chrom: str

        :param start: The start position of the locus.
        :type start: int

        :param end: The end position.
        :type end: int

        :param length: The length of the sequence to fetch.
        :type length: int

        Either an ``end`` or a ``length`` parameter has to be provided.

        The ranges are incluse, this means that (start, end) positions will
        both be included in the sequence.

        """
        if (end is None and length is None) or (end and length):
            raise TypeError("get_sequence needs either an 'end' OR 'length' "
                            "argument.")

        if length:
            end = start + length - 1

        try:
            seq = self.ref[str(chrom)][start - 1: end]
        except KeyError:
            seq = None

        if seq is None:
            error_message = "chr{}:{}-{} is an invalid genomic mapping"
            error_message = error_message.format(chrom, start, end)
            raise InvalidMapping(error_message)

        return str(seq.seq).upper()

    def close(self):
        self.ref.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


class InvalidMapping(Exception):
    """Exception representing an invalid mapping that we can't fix
    automatically.

    This can happen if the provided allele is incorrect for non-SNP variants.
    In this case we don't know if the locus is bad or if the sequence is bad.
    Since this is ambiguous, we raise this exception for the user to fix.

    """
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def check_snp_reference(snp, ref, flip):
    """Utility function to check if a snp has the correct reference allele.

    :param snp: The :py:class:`gepyto.structures.variants.SNP` object.
    :type snp: :py:class:`gepyto.structures.variants.SNP`

    :param ref: The :py:class:`Reference` reference object.
    :type ref: :py:class:`Reference`

    :param flip: A flag. If True, the return value is a variant with alleles
                 flipped if necessary. If False, a bool is returned: True
                 if the alleles are correct.
    :type flip: bool

    :returns: Either a :py:class:`gepyto.structures.variant.SNP` with
              flipped alleles or a bool.

    This is used internally by :py:class:`Reference`, but it is also
    available to users, but you need to provide a pyfaidx Fasta object.

    """
    # Get the reference.
    ref_allele = ref.get_nucleotide(snp.chrom, snp.pos)
    if ref_allele is None:
        s = "Can't find locus 'chr{}:{}' in the reference genome.".format(
            snp.chrom, snp.pos
        )
        raise ValueError(s)

    # Check if it is the correct allele.
    if ref_allele == snp.ref:
        # All is good in Russia.
        return snp if flip else True

    # Check if alt is the correct reference (then flip).
    if ref_allele == snp.alt:
        if flip:
            snp.ref, snp.alt = snp.alt, snp.ref
            return snp
        else:
            return False

    # The correct allele is nowhere to be found.
    raise ValueError("Invalid alleles for variant {}. Reference is {}".format(
        snp, ref_allele
    ))


def check_indel_reference(indel, ref, fix):
    """Check and/or fix alleles for Indels.

    :param ref: A reference object.
    :type ref: :py:class:`Reference`

    In fix mode, this function will try to standardise the alleles for the
    given indel. This means that the VCF format will be enforced. No "-"
    alleles will be authorized.

    _e.g._ ref: 'TC', alt: '-' will become ref: 'CTC', alt: 'C' given that
    the previous nucleotide in the reference is a 'C'. The position will be
    adjusted accordingly.

    In the regular mode, the only test will be that the `ref` allele is
    consistent with the reference. That is, the sequence given as the `ref`
    allele equals the one on the same length starting at `pos` in the genome.

    """
    ref_ok = False
    if fix:
        # Insertions:
        # We will use the VCF format, so if the ref is '-', we will make change
        # the start (-1) and add the preceding nucleotide.
        if indel.ref == "-":
            ref_ok = True
            indel.pos -= 1
            indel.ref = ref.get_nucleotide(indel.chrom, indel.pos)
            # We also need to prepend this nucleotide in the alt.
            indel.alt = indel.ref + indel.alt

        # Deletions:
        elif indel.alt == "-":
            indel.pos -= 1
            indel.alt = ref.get_nucleotide(indel.chrom, indel.pos)
            indel.ref = indel.alt + indel.ref

    # Verify the reference allele.
    if not ref_ok:
        seq = ref.get_sequence(chrom=indel.chrom, start=indel.pos,
                               length=len(indel.ref))
        if seq != indel.ref:
            # The provided sequence is not valid: we can't fix it because it
            # could be a bad mapping...
            err = ("chr{}:{} is an invalid locus. Verify that indel '{}' is "
                   "correct.").format(indel.chrom, indel.pos, indel.ref)
            if fix:
                raise InvalidMapping(err)
            else:
                return False

    return True if not fix else indel
