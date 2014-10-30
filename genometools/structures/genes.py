# Structures to handle genes and their associated transcripts and their
# proteins.

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

import urllib2
import contextlib
import json
import logging
import re

from .. import settings
from ..db import query_ensembl

__all__ = ["Gene", ]

class Gene(object):
    def __init__(self, **kwargs):
        """Python object representing a gene.

        Store the following information:

        Required

            - build: The genome build.
            - chrom: The chromosome.
            - start and end: The genomic positions for the gene.
            - xrefs: A dict of id mappings to multiple databases.
            - transcripts: A list of Transcript objects.

        Optional 

            - symbol: An HGNC symbol.
            - desc: A short description.

        You can only pass kwargs to build the genes. This makes for more
        eloquent code and avoids mistakes.

        """
        _PARAMETERS = {
            "build": str,
            "chrom": str,
            "start": int,
            "end": int,
            "xrefs": dict,
            "transcripts": list,
        }

        _OPTIONAL_PARAMETERS = {
            "symbol": str,
            "desc": str,
        }

        _ALL_PARAMS = dict(_PARAMETERS.items() + _OPTIONAL_PARAMETERS.items())

        # Store the passed parameters.
        for arg, val in kwargs.iteritems():
            if arg not in _ALL_PARAMS:
                raise Exception("Unknown parameter {}.".format(arg))
            setattr(self, arg, _ALL_PARAMS[arg](val))

        # Check that all required arguments were passed.
        for p in _PARAMETERS:
            if getattr(self, p, "-1") == "-1":
                raise Exception("Missing parameter {}.".format(p))

        # Some basic assertions.
        assert self.build in ("GRCh37", "GRCh38")
        assert re.match(r"([0-9]{1,2}|MT|X|Y)", self.chrom)
        assert self.start < self.end

    @classmethod
    def factory_symbol(cls, symbol, build=settings.BUILD):
        """Builds a gene object from it's HGNC symbol.

        :param symbol: The HGNC symbol.
        :type symbol: str

        :returns: The Gene object.
        :rtype: :py:class`Gene`

        """
        xrefs = Gene.get_xrefs_from_symbol(symbol)
        ensembl_id = xrefs.get("ensembl_id")

        if ensembl_id is None:
            raise Exception("Could not initialize gene from symbol {} because "
                "the Ensembl ID (ENSG) could not be found.".format(
                    symbol 
                ))

        url = "http://grch37." if build == "GRCh37" else "http://"
        url += ("rest.ensembl.org/overlap/id/{}"
                "?content-type=application/json"
                "&feature=gene"
                "&feature=transcript"
                "&feature=exon")

        response = query_ensembl(url.format(ensembl_id))

        # The response contains both gene information, the underlying
        # transcripts and the exons. We need to parse all that.
        for elem in response:
            if elem["feature_type"] == "gene":
                gene_info = _parse_gene()
            elif elem["feature_type"] == "transcript":
                transcripts.append(_parse_transcript())
            elif elem["feature_type"] == "exon":
                exons.append(_parse_exon())

        return Gene(
            **gene_info,
            xrefs=xrefs,
            symbol=symbol,
        )


    @classmethod
    def get_xrefs_from_symbol(cls, symbol):
        """Fetches the HGNC (HUGO Gene Nomenclature Commitee) service to get a gene ID for other databases.

        :param symbol: The gene symbol to query.
        :type symbol: str

        :returns: A dict representing information on the gene.
        :rtype: dict

        If no gene with this symbol can be found, `None` is returned.

        """

        url = "http://rest.genenames.org/search/symbol:{}"

        headers = {
            "Accept": "application/json",    
        }

        req = urllib2.Request(url.format(symbol), headers=headers)
        with contextlib.closing(urllib2.urlopen(req)) as stream:
            res = json.load(stream)

        # We take the top search hit and run a fetch.
        if res["response"]["numFound"] > 0:
            doc = res["response"]["docs"][0]
            symbol = doc["symbol"]
            assert doc["score"] == res["response"]["maxScore"]
        else:
            logging.warning("No gene with HGNC symbol {} found.".format(symbol))
            return None

        # Use the HGNC Fetch.
        url = "http://rest.genenames.org/fetch/symbol/{}"
        req = urllib2.Request(url.format(symbol), headers=headers)
        with contextlib.closing(urllib2.urlopen(req)) as stream:
            res = json.load(stream)       

        # Parse the cross references.
        if res["response"]["numFound"] > 0:
            doc = res["response"]["docs"][0]
            id_dict = {
                "name": doc.get("name"),
                "ncbi_id": doc.get("entrez_id"),
                "cosmic_id": doc.get("cosmic"),
                "refseq_ids": doc.get("refseq_accession"),
                "ensembl_id": doc.get("ensembl_gene_id"),
                "omim_ids": doc.get("omim_id"),
                "uniprot_ids": doc.get("uniprot_ids"),
                "ucsc_id": doc.get("ucsc_id"),
            }
            id_dict = {k: str(v) for (k, v) in id_dict.iteritems()}

            return id_dict
        else:
            raise Exception("No gene returned by HGNC fetch on "
                "symbol {}.".format(symbol))


class Transcript(object):
    pass

def _parse_gene(o):
    """Parse gene information from an Ensembl `overlap` query.

    :param o: A Python dict representing an entry from the response.
    :type o: dict

    :returns: A standardised dict of gene information.
    :rtype: dict

    """
    assert o["feature_type"] == "gene"

    d = {
        "desc": o.get("description"),
        "build": o.get("assembly_name"),
        "chrom": o.get("seq_region_name"),
        "start": o.get("start"),
        "end": o.get("end"),
    }

    return d

def _parse_transcript(o):
    """Parse transcript information from an Ensembl `overlap` query.

    :param o: A Python dict representing an entry from the response.
    :type o: dict

    :returns: A Transcript object representing the transcript.
    :rtype: :py:class`Transcript`

    """
    pass

def _parse_exon(o):
    """Parse exon information from an Ensembl `overlap` query.

    :param o: A Python dict representing an entry from the response.
    :type o: dict

    :returns: A tuple of (start, end) positions representing the exon.
    :rtype: tuple

    """
    assert o["feature_type"] == "exon"
    return (o["start"], o["end"])

