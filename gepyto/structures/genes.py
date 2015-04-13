
# Structures to handle genes and their associated transcripts and their
# proteins.

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


import sys
import contextlib
import json
import logging
import re

try:
    # Python 2 support
    from urllib2 import urlopen, Request
except ImportError:
    # Python 3 support
    from urllib.request import urlopen, Request

from .. import settings
from ..db.ensembl import query_ensembl
from ..db.appris import get_category_for_transcript
from . import sequences


__all__ = ["Gene", "Transcript"]


class Gene(object):
    """Python object representing a gene.

    Store the following information:

    Required

        - build: The genome build.
        - chrom: The chromosome.
        - start and end: The genomic positions for the gene.
        - strand: The strand (either 1 or -1).
        - xrefs: A dict of id mappings to multiple databases.
        - transcripts: A list of Transcript objects.

    Optional

        - symbol: An HGNC symbol.
        - desc: A short description.
        - exons: A list of pairs of positions for exons.

    You can only pass kwargs to build the genes. This makes for more
    eloquent code and avoids mistakes.

    """

    def __init__(self, **kwargs):

        _PARAMETERS = {
            "build": str,
            "chrom": str,
            "start": int,
            "end": int,
            "strand": int,
            "xrefs": dict,
        }

        _OPTIONAL_PARAMETERS = {
            "symbol": str,
            "desc": str,
            "transcripts": list,
            "exons": list,
            "biotype": str,
        }

        _ALL_PARAMS = dict(
            list(_PARAMETERS.items()) + list(_OPTIONAL_PARAMETERS.items())
        )

        # Store the passed parameters.
        for arg, val in kwargs.items():
            if arg not in _ALL_PARAMS:
                raise Exception("Unknown parameter {}.".format(arg))
            try:
                setattr(self, arg, _ALL_PARAMS[arg](val))
            except TypeError:
                raise TypeError("Invalid type for argument {}. Expected "
                                "a {}.".format(arg, _ALL_PARAMS[arg]))

        # Check that all required arguments were passed.
        for p in _PARAMETERS:
            if getattr(self, p, "-1") == "-1":
                raise Exception("Missing parameter {}.".format(p))

        # Some basic assertions.
        assert self.build in ("GRCh37", "GRCh38")
        assert re.match(r"([0-9]{1,2}|MT|X|Y)", self.chrom)
        assert self.start < self.end

    def __repr__(self):
        # Getting the best name
        name_repr = None
        if hasattr(self, "symbol"):
            name_repr = self.symbol
        else:
            for ref in ("ensembl_gene_id", "entrez_id", "ucsc_id"):
                if ref in self.xrefs:
                    name_repr = "{} ({})".format(self.xrefs[ref], ref)
                    break

        if name_repr is None:
            name_repr = "Unknown"

        return "<Gene:{} (chr{}:{}-{})>".format(
            name_repr,
            self.chrom,
            self.start,
            self.end
        )

    def __contains__(self, o):
        """Test if an object is in the gene.

        This compares the chrom, start and end (or pos) attributes of the
        given object with the gene.

        You can use the syntax: `rs123 in gene`.
        Where the `rs123` is a SNP object (or other) and `gene` is a `Gene`
        object.

        """
        if hasattr(o, "pos"):
            start = end = int(o.pos)
        elif hasattr(o, "start") and hasattr(o, "end"):
            start, end = [int(i) for i in (o.start, o.end)]

        if not hasattr(o, "chrom"):
            raise TypeError("Testing overlap with gene requires a `chrom` "
                            "attribute.")

        ret = self.start < start and self.end > end
        ret = ret and str(self.chrom) == str(o.chrom)
        return ret

    def get_ortholog_sequences(self):
        """Queries Ensembl to get Sequence objects representing orthologs.

        :returns: A list of :py:class:`gepyto.structures.sequences.Sequence`
        :rtype: list

        """
        return self._homology("orthologs")

    def get_paralog_sequences(self):
        """Queries Ensembl to get Sequence objects representing paralogs.

        :returns: A list of :py:class:`gepyto.structures.sequences.Sequence`
        :rtype: list

        """
        return self._homology("paralogs")

    def _homology(self, homo_type="orthologs"):
        """Retrieve homology information from Ensembl.

        :param homo_type: The type of homologuous sequence to query. Can be
                          either 'orthologs' or 'paralogs'.
        :type homo_type: str

        """

        homo_types = {
            "orthologs": "orthologues",
            "paralogs": "paralogues"
        }
        if homo_type not in homo_types.keys():
            raise Exception("Invalid homology type. Valid parameters are: "
                            "{}".foramt(homo_types.keys()))

        if "ensembl_id" not in self.xrefs:
            raise Exception("Can't retrieve homology information from Ensembl "
                            "without an 'ensembl_id' in the cross references "
                            "(xrefs).")

        url = ("http://rest.ensembl.org/homology/id/{}?"
               "content-type=application/json&"
               "type={}")
        url = url.format(self.xrefs["ensembl_id"], homo_types[homo_type])

        res = query_ensembl(url)

        homolog_sequences = []
        for hit in res["data"][0]["homologies"]:
            hit = hit["target"]

            seq_id = hit["protein_id"]
            seq = hit["align_seq"]

            seq_obj = sequences.Sequence(
                seq_id,
                "".join([c for c in seq if c != '-']),
                "AA",
                {
                    "species": hit["species"],
                    "species_ncbi_tax_id": hit["taxon_id"],
                    "perc_id": hit["perc_id"]
                }
            )
            homolog_sequences.append(seq_obj)

        return homolog_sequences

    @classmethod
    def factory_symbol(cls, symbol, build=settings.BUILD):
        """Builds a gene object from it's HGNC symbol.

        :param symbol: The HGNC symbol.
        :type symbol: str

        :returns: The Gene object.
        :rtype: :py:class:`Gene`

        """
        xrefs = Gene.get_xrefs_from_symbol(symbol, build=build)
        ensembl_id = xrefs.get("ensembl_id")

        if ensembl_id is None:
            raise Exception("Could not initialize gene from symbol {} because "
                            "the Ensembl ID (ENSG) could not be "
                            "found.".format(symbol))

        return Gene.factory_ensembl_id(ensembl_id, xrefs=xrefs, build=build)

    @classmethod
    def factory_ensembl_id(cls, ensembl_id, xrefs=None, build=settings.BUILD):
        """Builds a gene object from it's Ensembl ID.

        :param ensembl_id: The Ensembl ID.
        :type ensembl_id: str

        :returns: The Gene object.
        :rtype: :py:class:`Gene`

        """
        url = "http://grch37." if build == "GRCh37" else "http://"
        url += ("rest.ensembl.org/overlap/id/{}"
                "?content-type=application/json"
                "&feature=gene"
                "&feature=transcript"
                "&feature=exon")

        response = query_ensembl(url.format(ensembl_id))

        # The response contains both gene information, the underlying
        # transcripts and the exons. We need to parse all that.
        transcripts = []
        exons = []
        gene_info = None
        for elem in response:
            # We skip irrelevent entries (not on this build).
            if elem.get("assembly_name") and elem["assembly_name"] != build:
                continue

            if elem["feature_type"] == "gene" and elem["id"] == ensembl_id:
                gene_info = _parse_gene(elem)
            elif elem["feature_type"] == "transcript":
                tr = _parse_transcript(elem)
                # Sometimes, ensembl returns odd things that are unrelated,
                # we filter here for those.
                if tr.parent == ensembl_id:
                    transcripts.append(tr)
            elif elem["feature_type"] == "exon":
                exons.append(_parse_exon(elem))

        if gene_info is None:
            raise Exception("Could not find gene {} in Ensembl.".format(
                ensembl_id,
            ))

        if xrefs is None:
            xrefs = Gene.get_xrefs_from_ensembl_id(ensembl_id, build=build)

        if "symbol" in xrefs:
            gene_info["symbol"] = xrefs.pop("symbol")
        gene_info["xrefs"] = xrefs
        gene_info["transcripts"] = transcripts
        gene_info["exons"] = exons

        g = Gene(**gene_info)
        for tr in transcripts:
            tr.parent = g

        return g

    @classmethod
    def get_xrefs_from_ensembl_id(cls, ensembl_id, build=settings.BUILD):
        """Fetches the HGNC (HUGO Gene Nomenclature Commitee) service to get a
           gene ID for other databases.

        :param ensembl_id: The gene Ensembl ID to query.
        :type ensembl_id: str

        :returns: A dict representing information on the gene.
        :rtype: dict

        If no gene with this Ensembl ID can be found, `None` is returned.

        """
        return Gene.get_xrefs("ensembl_gene_id", ensembl_id)

    @classmethod
    def get_xrefs_from_symbol(cls, symbol, build=settings.BUILD):
        """Fetches the HGNC (HUGO Gene Nomenclature Commitee) service to get a
           gene ID for other databases.

        :param symbol: The gene symbol to query.
        :type symbol: str

        :returns: A dict representing information on the gene.
        :rtype: dict

        If no gene with this symbol can be found, `None` is returned.

        """
        return Gene.get_xrefs("symbol", symbol)

    @classmethod
    def get_xrefs(cls, field, query, build=settings.BUILD):
        """Fetches the HGNC (HUGO Gene Nomenclature Commitee) service to get a
           gene ID for other databases.

        :param field: A searchable fields.
        :type field: str

        :param query: The query.
        :type query: str

        :returns: A dict representing information on the gene.
        :rtype: dict

        If no gene with this symbol can be found, `None` is returned.

        """
        # Checks if the field is valid
        assert field in {"ccds_id", "ensembl_gene_id", "entrez_id", "hgnc_id",
                         "name", "symbol", "ucsc_id", "uniprot_ids", "vega_id"}

        # The url and header
        url = "http://rest.genenames.org/search/{field}:{query}"
        headers = {"Accept": "application/json"}

        req = Request(url.format(field=field, query=query), headers=headers)
        with contextlib.closing(urlopen(req)) as stream:
            res = json.loads(stream.read().decode())

        # We take the top search hit and run a fetch.
        if res["response"]["numFound"] > 0:
            doc = res["response"]["docs"][0]
            assert doc["score"] == res["response"]["maxScore"]

            # Use the HGNC Fetch.
            url = "http://rest.genenames.org/fetch/{field}/{query}"
            req = Request(url.format(field=field, query=query),
                          headers=headers)
            with contextlib.closing(urlopen(req)) as stream:
                res = json.loads(stream.read().decode())

            # Parse the cross references.
            if res["response"]["numFound"] > 0:
                doc = res["response"]["docs"][0]

                # Check the right symbol was found
                assert doc.get(field) == query

                id_dict = {
                    "name": doc.get("name"),
                    "symbol": doc.get("symbol"),
                    "ncbi_id": doc.get("entrez_id"),
                    "cosmic_id": doc.get("cosmic"),
                    "refseq_ids": doc.get("refseq_accession"),
                    "ensembl_id": doc.get("ensembl_gene_id"),
                    "omim_ids": doc.get("omim_id"),
                    "uniprot_ids": doc.get("uniprot_ids"),
                    "ucsc_id": doc.get("ucsc_id"),
                }

                return id_dict
            else:
                raise Exception("No gene returned by HGNC fetch on "
                                "{field} {query}.".format(field=field,
                                                          query=query))

        elif field == "ensembl_gene_id":
            # Get from Ensembl
            url = "http://grch37." if build == "GRCh37" else "http://"
            url += ("rest.ensembl.org/xrefs/id/{}"
                    "?content-type=application/json")
            url = url.format(query)

            response = query_ensembl(url)

            id_dict = {"ensembl_gene_id": query}
            for db_info in response:
                # Entrez
                if db_info["dbname"] == "EntrezGene":
                    id_dict["entrez_id"] = db_info["primary_id"]
                    continue

                # UniPROT
                if db_info["dbname"] == "Uniprot_gn":
                    id_dict["uniprot_ids"] = db_info["primary_id"]

            return id_dict

        else:
            return {}


class Transcript(object):
    """Python object representing a transcript.

    Store the following information:

    Required

        - build: The genome build.
        - chrom: The chromosome.
        - start and end: The genomic positions for the gene.
        - enst: The corresponding Ensembl transcript id.

    Optional

        - appris_cat: The APPRIS category.
        - parent: The corresponding Gene object.
        - biotype: The biotype as given by Ensembl.

    """

    def __init__(self, **kwargs):
        dummy = lambda x: x

        _PARAMETERS = {
            "build": str,
            "chrom": str,
            "start": int,
            "end": int,
            "enst": str,
        }

        _OPTIONAL_PARAMETERS = {
            "appris_cat": str,
            "parent": dummy,
            "biotype": str,
        }

        _ALL_PARAMS = dict(
            list(_PARAMETERS.items()) + list(_OPTIONAL_PARAMETERS.items())
        )

        # Store the passed parameters.
        for arg, val in kwargs.items():
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

    def get_sequence(self, seq_type="genomic"):
        """Build a Sequence object representing the transcript.

        :param seq_type: This can be either genomic, cds, cdna or protein.
        :type seq_type: str

        :returns: A Sequence object representing the feature.
        :rtype: :py:class:`gepyto.structures.sequences.Sequence`

        """

        # Get the sequence from Ensembl.
        seq_types = ("genomic", "cds", "cdna", "protein")
        if seq_type not in seq_types:
            raise Exception("Invalid sequence type ({}). Known types are: "
                            "{}".format(seq_type, ", ".join(seq_types)))

        url = ("http://rest.ensembl.org/sequence/id/{}?"
               "content-type=application/json&"
               "type={}")
        url = url.format(self.enst, seq_type)

        res = query_ensembl(url)

        # Build the Sequence.
        seq = sequences.Sequence(
            res["id"],
            res["seq"],
            "AA" if seq_type == "protein" else "DNA",
        )

        return seq

    @classmethod
    def factory_position(cls, region, build=settings.BUILD):
        """Gets a list of transcripts overlapping with the given position.

        :param pos: A genomic position of the form `chr2:12345-12347`.
        :type pos: str

        :returns: A list of :py:class:`Transcript`
        :rtype: list

        This method uses the Ensembl API.

        """

        assert re.match(r"^chr([0-9]{1,2}|MT|X|Y):[0-9]+-[0-9]+$", region)
        region = region.lstrip("chr")
        region = region.replace("-", "..")

        url = "http://grch37." if build == "GRCh37" else "http://"
        url += ("rest.ensembl.org/overlap/region/homo_sapiens/{}"
                "?feature=transcript"
                "&content-type=application/json")
        url = url.format(region)

        response = query_ensembl(url)
        transcripts = []
        for tr in response:
            transcripts.append(_parse_transcript(tr))

        return transcripts

    def __repr__(self):
        return "<Transcript:{} (chr{}:{}-{})>".format(
            self.enst,
            self.chrom,
            self.start,
            self.end,
        )


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
        "strand": o.get("strand"),
        "biotype": o.get("biotype"),
    }

    return d


def _parse_transcript(o):
    """Parse transcript information from an Ensembl `overlap` query.

    :param o: A Python dict representing an entry from the response.
    :type o: dict

    :returns: A Transcript object representing the transcript.
    :rtype: :py:class:`Transcript`

    """
    assert o["feature_type"] == "transcript"

    d = {
        "build": o.get("assembly_name"),
        "chrom": o.get("seq_region_name"),
        "start": o.get("start"),
        "end": o.get("end"),
        "enst": o.get("id"),
        "biotype": o.get("biotype"),
        "parent": o.get("Parent"),
    }

    # Get the APPRIS annotation.
    try:
        d["appris_cat"] = get_category_for_transcript(d["enst"])
    except Exception:
        pass

    return Transcript(**d)


def _parse_exon(o):
    """Parse exon information from an Ensembl `overlap` query.

    :param o: A Python dict representing an entry from the response.
    :type o: dict

    :returns: A tuple of (start, end) positions representing the exon.
    :rtype: tuple

    """
    assert o["feature_type"] == "exon"
    return (o["start"], o["end"])
