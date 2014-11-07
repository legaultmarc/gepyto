
# Utilities to interact with RegulomeDB.
# This module contains code strictly to interact with RegulomeDB and not to 
# interpret the response or the underlying biology.
# RegulomeDB: http://regulome.stanford.edu

# This file is part of genometools.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = ("Copyright 2014 Marc-Andre Legault and Louis-Philippe Lemieux "
                 "Perreault. All rights reserved.")
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


import re
import time
import logging

try:
    # Python 2 support
    from urllib import urlencode
    from urllib2 import Request, urlopen, HTTPError
except ImportError:
    # Python 3 support
    from urllib.error import HTTPError
    from urllib.parse import urlencode
    from urllib.request import Request, urlopen

from ..structures.variants import SNP, Indel


__all__ = ["query_regulomedb", ]


def query_regulomedb(query, build="GRCh37"):
    """Query RegulomeDB for a particular region.

    :param query: The genomic region to query (hg19).
    :type query: str or :py:class:`genometools.structures.variants.Variant`

    :param build: The genomic build (only GRCh37 is accepted by RegulomeDB).
    :type build: str

    :returns: A list of annotation for the region or variant.
    :rtype: list

    The ``query`` might be either a genomic region (0-based) or a variant object
    (:py:class:`genometools.structures.variants.SNP` or
    :py:class:`genometools.structures.variants.Indel`).

    """

    # This tool only works with hg19 for now
    assert(build == "GRCh37")

    # Generating the data to send
    data = "{chrom}:{start}-{end}"
    if isinstance(query, SNP):
        # The query is a SNP
        data = data.format(chrom=query.chrom, start=query.pos-1, end=query.pos)

    elif isinstance(query, Indel):
        data = data.format(chrom=query.chrom, start=query.start-1,
                           end=query.end)

    elif isinstance(query, str):
        # The query is a string
        # TODO: Check the query string
        data = query

    else:
        raise ValueError("{}: not a valid query".format(query))

    # Submitting the data
    url = "http://regulome.stanford.edu/results"
    values = {"data": data}
    req = Request(url, urlencode(values).encode("utf-8"))

    try:
        response = urlopen(req).read().decode()
    except HTTPError:
        logging.warning("Request failed for query '{}'.".format(data)) 
        return []

    # Getting the sid value
    sid = None
    r = re.search("name=['\"]sid['\"] value=['\"](\\w+)['\"]", response)
    if r is not None:
        sid = r.group(1)

    # TODO: Better check if sid is still None...
    assert(sid is not None)

    # Fetching the result
    url = "http://regulome.stanford.edu/download"
    values = {"format": "full",
              "sid": sid,
              "download_token_value_id": str(int(time.time()*100))}
    req = Request(url, urlencode(values).encode("utf-8"))

    return _parse_regulomedb(urlopen(req).read().decode())


def _parse_regulomedb(regulome_output):
    """Takes a file object containing the output from CADD and parses it.

    :param regulome_output: An stream representing the output from RegulomeDB.
    :type regulome_output: str

    :returns: A list representation of the RegulomeDB output.
    :rtype: list

    """
    # Is there a result?
    if regulome_output == "":
        return []

    # Splitting the stream by lines
    regulome_output = regulome_output.splitlines()

    if len(regulome_output) == 1:
        # There is only the header
        return []

    # The header
    header = {name: i for i, name in enumerate(regulome_output[0].split("\t"))}
    for name in ("#chromosome", "coordinate", "rsid", "hits", "score"):
        if name not in header:
            logging.warnning("Invalid RegulomeDB output")
            return []

    # The data
    result = []
    for line in regulome_output[1:]:
        row = line.split("\t")
        result.append({
            "chrom": re.sub("^chr", "", row[header["#chromosome"]]),
            "pos": int(row[header["coordinate"]]),
            "rs": row[header["rsid"]],
            "hits": [re.split(r"\|", i)
                                for i in re.split(", ", row[header["hits"]])],
            "score": row[header["score"]],
        })

    return result

