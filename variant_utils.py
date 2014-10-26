# Utilities to handle variant data.

import re
import json
import urllib2
import contextlib
import logging

from settings import BUILD
from data_structures import SNP

def ensembl_snp_in_region(region, build=BUILD):
    """Queries a genome region of the form chr3:123-456 for variants using Ensembl API.

    :param region: The region to query.
    :type region: str

    :param build: The genome build to use (GRCh37 or GRCh38).
    :type build: str

    :returns: A list of :class:`data_structures.SNP`.
    :rtype: list

    Note: For now, as the name implies, this function is limited to SNPs and
    not to any kind of short variation.

    """

    url = ("rest.ensembl.org/overlap/region/homo_sapiens/{region}"
           "?feature=variation"
           "&content-type=application/json")

    if build == "GRCh37":
        url = "http://grch37." + url
    elif build != "GRCh38":
        raise Exception("Unknown build '{}'.".format(build))

    if region.startswith("chr"):
        region = region.lstrip("chr")

    url = url.format(region=region)

    with contextlib.closing(urllib2.urlopen(url)) as stream:
        res = json.load(stream)            

    variants = []
    for variant in res:
        if variant["start"] == variant["end"]:
            # Variant is a SNP.
            # Check some stuff.
            assert variant["feature_type"] == "variation"
            assert variant["assembly_name"] == build

            # Build the variant.
            chrom = str(variant["seq_region_name"])
            pos = variant["start"]
            rs = str(variant["id"])
            ref = str(variant["alt_alleles"][0])
            alt = str(variant["alt_alleles"][1])
            if len(variant["alt_alleles"]) > 2:
                logging.warning("{} is not biallelic (only one allele "
                    "considered).".format(rs))

            if type(rs) is str and not rs.startswith("rs"):
                # We ignore the id if it's not from dbSNP.
                rs = None

            variant_obj = SNP(chrom, pos, rs, ref, alt)
            variants.append(variant_obj)

    return variants

