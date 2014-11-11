
# Utilities to handle variant data.

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


import re
import contextlib
import logging

from ..settings import BUILD
from ..structures.variants import SNP, Indel
from ..db.ensembl import query_ensembl


def ensembl_variants_in_region(region, build=BUILD):
    """Queries a genome region of the form chr3:123-456 for variants using Ensembl API.

    :param region: The region to query.
    :type region: str

    :param build: The genome build to use (GRCh37 or GRCh38).
    :type build: str

    :returns: A list of :py:class:`genometools.structures.variants.Variant`
              or a subclass.
    :rtype: list

    """

    url = ("rest.ensembl.org/overlap/region/homo_sapiens/{region}"
           "?feature=variation"
           "&content-type=application/json")

    if build == "GRCh37":
        url = "http://grch37." + url
    elif build == "GRCh38":
        url = "http://" + url
    else:
        raise Exception("Unknown build '{}'.".format(build))

    if region.startswith("chr"):
        region = region.lstrip("chr")

    url = url.format(region=region)

    res = query_ensembl(url)

    variants = []
    for variant in res:
        # Check some stuff.
        assert variant["feature_type"] == "variation"
        assert variant["assembly_name"] == build

        # Build the variant.
        chrom = str(variant["seq_region_name"])
        start = variant["start"]
        end = variant["end"]
        rs = str(variant["id"])
        if type(rs) is str and not rs.startswith("rs"):
            # We ignore the id if it's not from dbSNP.
            rs = None

        if len(variant["alt_alleles"]) < 2:
            # Weirdly, we have less than two alleles.
            logging.warning("RS: {} has only one allele (ignored).".format(rs))
            continue

        ref = str(variant["alt_alleles"][0])

        is_snp = True
        for allele in variant["alt_alleles"]:
            if allele == "-" or len(allele) > 1:
                is_snp = False

        for alt in variant["alt_alleles"][1:]:
            if is_snp:
                variant_obj = SNP(chrom, start, rs, ref, str(alt))
            else:
                variant_obj = Indel(chrom, start, end, rs, ref, str(alt))

            variants.append(variant_obj)

    if len(variants) == 0:
        logging.warning("No SNP detected in region {}.".format(region))

    return variants

