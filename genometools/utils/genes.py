
# Utilities to handle gene data.

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


# TODO: CLEAN UP IMPORT

import logging

from ..settings import BUILD
from ..structures.genes import Gene
from ..db.ensembl import query_ensembl


def ensembl_genes_in_region(region, build=BUILD):
    """Queries a genome region of the form chr3:123-456 for genes using Ensembl API.

    :param region: The region to query.
    :type region: str

    :param build: The genome build to use (GRCh37 or GRCh38).
    :type build: str

    :returns: A list of :class:`data_structures.Gene`.
    :rtype: list

    """

    url = ("rest.ensembl.org/overlap/region/homo_sapiens/{region}"
           "?feature=gene"
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

    genes = []
    for gene in res:
        # Check some stuff
        assert gene["feature_type"] == "gene"
        assert gene["assembly_name"] == build

        # Building the gene
        g_obj = Gene.factory_symbol(gene["external_name"], build=build)

        genes.append(g_obj)

    if len(genes) == 0:
        logging.warning("No gene detected in region {}.".format(region))

    return genes

