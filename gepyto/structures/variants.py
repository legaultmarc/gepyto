
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


import re
import contextlib
import logging
import sqlite3

import sys
import traceback

from .. import settings
from ..reference import Reference, check_indel_reference
from ..db import query_ensembl
from . import genes


__all__ = ["SNP", "Indel", "Variant"]


class Variant(object):
    """Base class for genome Variants.

    """

    def __init__(self, *args, **kwargs):
        """Super constructor for Variant objects.

        It will initalize the object by assuming one of the following
        scenarios:

        1) The position arguments (*args) correspond to the elements of the
           _PARAMETERS argument which is overriden by the subclasses when the
           parent constructor is called.

        2) All of the _PARAMETERS elements should be specified as keyword
           arguments (**kwargs).

        A sanity check will be done after initalization to make sure that all
        of the parameters are defined.

        """

        if "_PARAMETERS" in kwargs:
            _PARAMETERS = kwargs.pop("_PARAMETERS")
        else:
            _PARAMETERS = ["chrom", "start", "end"]

        if len(args) == len(_PARAMETERS):
            for i, p in enumerate(_PARAMETERS):
                setattr(self, p, args[i])
        else:
            keys = list(kwargs.keys())
            for k in keys:
                if k not in _PARAMETERS:
                    raise Exception("Unknown parameter {}.".format(k))

                setattr(self, k, kwargs.pop(k))

        # Make sure everything was set.
        for p in _PARAMETERS:
            if getattr(self, p, "-1") == "-1":
                raise Exception("Missing parameter {}.".format(p))

    def vcf_header(self):
        """Returns a valid VCF header line. """
        return "\t".join([
            "#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
        ])

    def in_gene(self, gene):
        """Method to test if the variant is in the gene.

        :param: The gene object to test against.
                By duck typing, this could be anything with a ``chrom``,
                ``start`` and ``end`` (or a ``pos``).

        :returns: True is the variant is in the gene (region comparison).
        :rtype: bool

        """

        if hasattr(self, "pos"):
            start = self.pos
        elif hasattr(self, "start"):
            start = self.start
        else:
            raise NotImplementedError()

        if hasattr(self, "end"):
            end = self.end
        elif hasattr(self, "pos"):
            end = self.pos
        else:
            raise NotImplementedError()

        chrom = self.chrom

        ret = (str(chrom) == str(gene.chrom) and
               int(start) > int(gene.start) and
               int(end) < int(gene.end))

        return ret

    def get_position(self):
        """Abstract method.

        Returns the genomic position of the variant using the standard format:
        chrXX:00000-11111.

        """
        raise NotImplementedError()

    @staticmethod
    def from_ensembl_api(rs, build=settings.BUILD):
        """Builds the correct Variant subclass for the specified rs number.

        :param rs: The rs number describing the variant of interest.
        :type rs: str

        :param build: The genome build (either GRCh37 or GRCh38, default is
                      defined in the settings.
        :type build: str

        This method will fetch the data from the Ensembl rest api and try
        to guess the right Variant subclass (e.g. SNP or Indel).

        """

        if build == "GRCh37":
            url = ("http://grch37.rest.ensembl.org/variation/homo_sapiens/"
                   "{snp}?content-type=application/json")
        elif build == "GRCh38":
            url = ("http://rest.ensembl.org/variation/homo_sapiens/{snp}"
                   "?content-type=application/json")
        else:
            raise Exception("Unknown build '{}'.".format(build))

        url = url.format(snp=rs)
        response = query_ensembl(url)

        if response is None:
            return response

        for mapping in response["mappings"]:
            if mapping["assembly_name"] == build:
                # This is a SNP
                if mapping["start"] == mapping["end"]:
                    pos = "chr{}_{}".format(
                        mapping["location"].split("-")[0],
                        mapping["allele_string"]
                    )

                    # Note that for multi (>2) allelic loci, this will
                    # return a list of variants.
                    return SNP.from_str(pos, rs=rs)

                # Otherwise an Indel
                else:
                    return Indel._parse_ensembl_indel(rs, mapping)


class Indel(Variant):
    """Class representing short insertions/deletions (Indels).

    Either initialize with the parameters corresponding to:
    ``chrom, pos, rs, ref, alt`` or by using the corresponding
    named parameters.

    The notation we are using is consistent with the VCF format, but it is
    different from some API (e.g. Ensembl). We will try to standardize
    everything to comply with the VCF format.

    The latter represents deletions by using the preceding nucleotide for the
    reference.

    _e.g._ If the genomic sequence is AAGAA -> AAAA (deletion of the G)
    the indel will be represented as follows:

        - Start: 2
        - Ref: `AG`
        - Alt: `A`
        - length: 1 (difference between allele lengths)

    _e.g._ If the genomic sequence is AAGAA -> AAGCAA (insertion of the C)
    the indel will be represented as follows:

        - Start: 3
        - Ref: `G`
        - Alt: `GC`
        - length: 1

    """

    def __init__(self, *args, **kwargs):
        if "skip_allele_check" in kwargs:
            skip_allele_check = kwargs.pop("skip_allele_check")
        else:
            skip_allele_check = False

        _PARAMETERS = [
            "chrom",
            "pos",
            "rs",
            "ref",
            "alt",
        ]
        kwargs["_PARAMETERS"] = _PARAMETERS
        super(Indel, self).__init__(*args, **kwargs)

        self.pos = int(self.pos)
        self.ref = self.ref.lower()
        self.alt = self.alt.lower()
        try:
            assert re.match(settings.CHROM_REGEX, str(self.chrom))
            assert self.rs is None or re.match(r"^rs[0-9]+$", self.rs)
            assert type(self.pos) is int
            assert type(self.ref) is str
            assert type(self.alt) is str
            if not skip_allele_check:
                assert "-" not in self.ref
                assert "-" not in self.alt
        except AssertionError as e:
            logging.critical(
                "Assertion failed constructing the Indel object. \n"
                "Parameters were: \n"
                "\tchrom: {chrom} ({chrom_type})\n"
                "\tpos: {pos} ({pos_type})\n"
                "\trs: {rs} ({rs_type})\n"
                "\tref: {ref} ({ref_type})\n"
                "\talt: {alt} ({alt_type})\n".format(
                    chrom=self.chrom, chrom_type=type(self.chrom),
                    pos=self.pos, pos_type=type(self.pos),
                    rs=self.rs, rs_type=type(self.rs),
                    ref=self.ref, ref_type=type(self.ref),
                    alt=self.alt, alt_type=type(self.alt),
                )
            )
            traceback.print_tb(sys.exc_info()[2])
            raise e

    @property
    def length(self):
        """Computes the size of the indel.

        This is the difference between the length of the alleles.

        """
        return abs(len(self.ref) - len(self.alt))

    @classmethod
    def _parse_ensembl_indel(cls, rs, mapping):
        chrom = str(mapping["seq_region_name"])
        pos = mapping["start"]
        alleles = str(mapping["allele_string"]).split("/")
        ref = alleles[0]
        alts = alleles[1:]

        indels = []
        with Reference() as reference:
            for alt in list(alts):
                indel = cls(chrom, pos, rs, ref, alt, skip_allele_check=True)
                indel = check_indel_reference(indel, reference, True)
                indels.append(indel)

        return indels

    def get_position(self, zero_based=False):
        """Returns an indel in the standard chrXX:POS notation.

        :param zero_based: Indicates if zero based coordinates should be used.
                           This is False by default (1 based coordinates).

        """
        pos = self.pos - 1 if zero_based else self.pos
        return "chr{}:{}".format(self.chrom, pos)

    def vcf_line(self):
        """Returns a line describing the current variant as expected by the VCF
           format.

        """
        return "\t".join(str(i) for i in [
            self.chrom,
            self.pos,
            self.rs if self.rs else ".",
            self.ref,
            self.alt,
            ".",  # No Quality
            "PASS",  # PASS for filter
            ".",  # Info
        ])

    def __repr__(self):
        return "chr{}:{}_{}/{}".format(
            self.chrom, self.pos, self.ref, self.alt
        )

    def __eq__(self, other):
        if type(self) is not type(other):
            return False

        return (
            self.chrom == other.chrom and
            self.pos == other.pos and
            self.ref == other.ref and
            self.alt == other.alt
        )

    @classmethod
    def from_ensembl_api(cls, rs, build=settings.BUILD):
        """Gets the information from the Ensembl REST API.

        :param rs: The rs number for the variant of interest.
        :param build: The genome build (e.g. GRCh37 or GRCh38).

        """
        variants = Variant.from_ensembl_api(rs, build)
        return [v for v in variants if v.__class__ is cls]


class SNP(Variant):
    """Class representing a Single Nucleotide Polymorphism (SNP).

    Instances can be created in two ways: either by providing ordered fields:
    ``chrom, pos, rs, ref, alt`` or by using named parameters.

    """
    def __init__(self, *args, **kwargs):
        _PARAMETERS = [
            "chrom",
            "pos",
            "rs",
            "ref",
            "alt",
        ]
        kwargs["_PARAMETERS"] = _PARAMETERS
        super(SNP, self).__init__(*args, **kwargs)

        # Make sure everthing has the right type.
        self.pos = int(self.pos)
        self.ref = self.ref.lower()
        self.alt = self.alt.lower()
        try:
            assert re.match(settings.CHROM_REGEX, str(self.chrom))
            assert self.rs is None or re.match(r"^rs[0-9]+$", self.rs)
            assert type(self.ref) is str and len(self.ref) == 1
            assert type(self.alt) is str and len(self.alt) == 1
        except AssertionError as e:
            logging.critical(
                "Assertion failed constructing the SNP object. \n"
                "Parameters were: \n"
                "\tchrom: {chrom} ({chrom_type})\n"
                "\tpos: {pos} ({pos_type})\n"
                "\trs: {rs} ({rs_type})\n"
                "\tref: {ref} ({ref_type})\n"
                "\talt: {alt} ({alt_type})\n".format(
                    chrom=self.chrom, chrom_type=type(self.chrom),
                    pos=self.pos, pos_type=type(self.pos),
                    rs=self.rs, rs_type=type(self.rs),
                    ref=self.ref, ref_type=type(self.ref),
                    alt=self.alt, alt_type=type(self.alt),
                )
            )
            traceback.print_tb(sys.exc_info()[2])
            raise e

    def get_position(self, zero_based=False):
        """Returns a variant in the standard chrXX:pos notation.

        :param zero_based: Indicates if zero based coordinates should be used.
                           This is False by default (1 based coordinates).

        """
        pos = self.pos - 1 if zero_based else self.pos
        return "chr{}:{}".format(self.chrom, pos)

    def vcf_line(self):
        """Returns a line describing the current variant as expected by the VCF
        format.

        """
        return "\t".join(str(i) for i in [
            self.chrom,
            self.pos,
            self.rs if self.rs else ".",
            self.ref,
            self.alt,
            ".",  # No Quality
            "PASS",  # PASS for filter
            ".",  # Info
        ])

    @classmethod
    def from_ensembl_api(cls, rs, build=settings.BUILD):
        """Gets the information from the Ensembl REST API.

        :param rs: The rs number for the variant of interest.
        :param build: The genome build (e.g. GRCh37 or GRCh38).

        """
        variants = Variant.from_ensembl_api(rs, build)
        return [v for v in variants if v.__class__ is cls]

    @classmethod
    def from_str(cls, s, rs=None):
        """Parses a variant object from a str of the form chrXX:YYYY_R/A.

        :param s: The string to parse the SNP from (Format: chrXX:YYY_R/A).
        :param rs: An optional parameter specifying the rs number.

        :returns: A list of SNP objects.

        If it is a multi-allelic loci, the list will contain one SNP per
        alternative allele.

        If it is not, this will be a list of length 1...

        """
        s = s.lstrip("chr")
        chrom, s = s.split(":")
        pos, s = s.split("_")
        alleles = s.split("/")
        ref = alleles[0]
        alts = list(alleles[1:])

        var_list = []
        for alt in alts:
            v = cls(chrom, pos, rs, ref, alt)
            var_list.append(v)

        return var_list

    def __eq__(self, other):
        if type(self) is not type(other):
            return False

        return (
            self.chrom == other.chrom and
            self.pos == other.pos and
            self.ref == other.ref and
            self.alt == other.alt
        )

    def __repr__(self):
        return "chr{}:{}_{}/{}".format(self.chrom, self.pos, self.ref,
                                       self.alt)
