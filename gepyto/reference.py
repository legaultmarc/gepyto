
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

from . import settings

from pyfaidx import Fasta

class Reference(object):
    """Interface to the human genome reference file.
    
    This class uses ``pyfaidx`` to parse the genome reference file referenced
    by ``settings.REFERENCE_PATH``.

    This can only be a single plain fasta file.

    Also note that if the path is not in the ``~/.gtconfig/gtrc.ini`` file,
    gepyto will look for an environment variable named ``REFERENCE_PATH``.

    .. todo::

        This should transparently fallback to a remote version (e.g. Ensembl
        REST API) to query the reference.

    """
    def __init__(self):
        self.ref = Fasta(settings.REFERENCE_PATH)

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
                        "A SNP object has to be provided.")

        if not (hasattr(variant, "chrom") and 
                hasattr(variant, "pos") and
                hasattr(variant, "ref") and
                hasattr(variant, "alt")):
            raise TypeError(type_message)

        if len(variant.ref) == len(variant.alt) == 1:
            return check_snp_reference(variant, self.ref, flip)
        else:
            # return check_indel_reference(variant, self.ref, flip)
            raise TypeError(type_message)

    def get_nucleotide(self, chrom, pos):
        """Get the nucleotide at the given genomic position. """
        return str(self.ref.get(chrom)[pos - 1].seq)

    def close(self):
        self.ref.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


def check_snp_reference(snp, ref, flip):
    """Utility function to check if a snp has the correct reference allele.

    :param snp: The :py:class:`gepyto.structures.variants.SNP` object.
    :type snp: :py:class:`gepyto.structures.variants.SNP`

    :param ref: The ``pyfaidx`` reference object.
    :type ref: :py:class:`pyfaidx.Fasta`

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
    ref_chrom = ref.get(snp.chrom)
    if ref_chrom is None:
        s = "Can't find chromosome '{}' in the reference genome.".format(
            snp.chrom
        )
        raise ValueError(s)

    # Now check the allele at this position. We assome that the SNP positions
    # are 1-based and that the reference object is 0-based (default with
    # pyfaidx.
    ref_allele = ref_chrom[snp.pos - 1].seq.upper()

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


def check_indel_reference(indel, ref, flip):
    """TODO """
    # Get the reference.
    ref_chrom = ref.get(indel.chrom)
    if ref_chrom is None:
        s = "Can't find chromosome '{}' in the reference genome.".format(
            indel.chrom
        )
        raise ValueError(s)

    # Insertions:
    # We will use the VCF format, so if the ref is '-', we will make change
    # the start (-1) and add the preceding nucleotide.
    if indel.ref == "-":
        indel.start -= 1
        # Remember that we are using 0-based indexing in pyfaidx.
        prev_nuc = ref_chrom[indel.start - 1]
        indel.ref = prev_nuc
        # We also need to prepend this nucleotide in the alt.
        indel.alt = prev_nuc + indel.alt

    # For insertions, we also need to check that the end - start makes sense.

    # Deletions:
    # The ref should be ok for this case.

    return
