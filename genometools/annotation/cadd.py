
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
import time
import uuid
import os
import gzip

from .. import structures as struct


__all__ = ["cadd_score", ]

def cadd_score(variants):
    """Annotate the variants using CADD (cadd.gs.washington.edu).

    :param vcf: A list of Variant (or subclass) objects.
    :type vcf: :py:class:`genometools.structures.variants.Variant`

    :returns: A list of annotations.
    :rtype: list

    The format for the annotations is a list of tuples of the form
    (``Transript``, ``Variant``, ``C score``, ``info``). If the annotation is
    for an Ensembl regulatory feature, the ID (ENSR) replaces the 
    :py:class:`genometools.structures.genes.Transcript` object.

    The info dictionary contains extra annotations from the CADD output.

    """

    # We import here not to break on optional modules.
    import requests # Install http://docs.python-requests.org/en/latest/

    # We need to send a form with the following information:
    # file: a vcf file.
    # inclAnno=Yes
    # to http://cadd.gs.washington.edu/upload
    # enctype="multipart/form-data"
    # with POST

    # Then we parse the html to find a link to a anno_.+\.tsv\.gz
    # We check every 20s to see if we can download the file.

    # TODO: We need to parse the html instead of this hacky regex.

    form = {
        "inclAnno": "Yes",
    }

    # Create a VCF file from the list of variants.
    fn = str(uuid.uuid4()) + ".vcf"
    with open(fn, "wb") as f:
        f.write("##fileformat=VCFv4.1\n")
        f.write(variants[0].vcf_header() + "\n")
        for v in variants:
            f.write(v.vcf_line() + "\n")

    with open(fn, "rb") as f:
        file_upload = {
            "file": f,
        }

        r = requests.post(
            "http://cadd.gs.washington.edu/upload",
            data=form,
            files=file_upload
        )
    os.remove(fn)

    if r.status_code == requests.codes.ok:
        feedback = r.text
        # Get the link
        m = re.search(r"href=\".+/(.+\.tsv\.gz)\"", feedback)
        if m:
            dld_fn = m.group(1)

            # Retry to download every 5 minutes.
            download_url = "http://cadd.gs.washington.edu/static/finished/{}"
            download_url = download_url.format(dld_fn)
            success = False
            for i in xrange(24): # Try, wait for 2 minutes, retry.
                time.sleep(2 * 60)

                r = requests.get(download_url, stream=True)
                if r.status_code == 200:
                    # Results are ready, YAY!
                    with open(fn, "wb") as f:
                        for chunk in r.iter_content():
                            f.write(chunk)

                    success = True
                    break

            if not success:
                raise Exception("Timeout while fetching CADD Results.")

            # Read the file and add the C Scores to a copy of the input file.
            with gzip.open(fn) as f:
                annotations = _parse_annotation(f)
            os.remove(fn)

            return annotations

        else:
            print(m)
            print(feedback)
            raise Exception("Could not find download link in CADD Response.")

    else:
        raise Exception("Request to CADD Failed (HTTP {}).".format(
            r.status_code 
        ))


def _parse_annotation(cadd_output):
    """Takes a file object containing the output from CADD and parses it.

    :param cadd_output: An open file object representing the output from the 
                        CADD website.
    :type cadd_output: Mosty likely :py:class:`gzip.GzipFile`

    :returns: A list of annotation tuples of the form (``Transript``, 
              ``Variant``, ``C score``, ``info``). If the annotation is for 
              an Ensembl regulatory feature, the ID (ENSR) replaces the
              ``Transcript`` object.
    :rtype: list

    """

    # We will build a list of tuples as documented.
    annotations = []

    # We keep objects in memory to avoid redundancy.
    variants = {}
    transcripts = {}

    # Parse the CADD annotation file.
    header = None
    for line in cadd_output:
        if line.startswith("##"):
            # Skip header lines that are not column labels.
            continue

        line = line.rstrip()
        line = [i.upper() for i in line.split("\t")]

        if line[0] == "#CHROM":
            line[0] = line[0].lstrip("#")
            header = {col: idx for idx, col in enumerate(line)}
            continue

        # Check if we already have this variant.
        chrom = line[header["CHROM"]]
        pos = int(line[header["POS"]])
        ref = line[header["REF"]]
        alt = line[header["ALT"]]
        var_type = line[header["TYPE"]]

        tu = (chrom, pos, ref, alt, var_type)
        if tu in variants:
            v = variants[tu]
        else:
            # We need to create the variant.
            if var_type == "SNV" and len(ref) == 1 and len(alt) == 1:
                # Create a SNP.
                v = struct.variants.SNP(
                    chrom=chrom,
                    pos=pos,
                    ref=ref,
                    alt=alt,
                    rs=None,
                )
            else:
                # Create an indel.
                v = struct.variants.Indel(
                    chrom=chrom,
                    start=pos,
                    end=pos + int(line[header["LENGTH"]]),
                    ref=ref,
                    alt=alt,
                    rs=None,
                )

            # We save it for future use.
            variants[tu] = v

        # Check if we already have this transcript.
        feature_id = line[header["FEATUREID"]]
        if feature_id in transcripts:
            feature = transcripts[feature_id]
        elif feature_id.startswith("ENSR"):
            feature = feature_id
        else:
            # We need to create the transcript.
            ov_transcripts = struct.genes.Transcript.factory_position(
                "chr{}:{}-{}".format(chrom, pos, pos)
            )
            for tr in ov_transcripts:
                if tr.enst == feature_id:
                    feature = tr
                    transcripts[feature_id] = tr
                    break

        c = float(line[header["CSCORE"]])
        info = {
            "poly": line[header["POLYPHENCAT"]],
            "sift": line[header["SIFTCAT"]],
            "gene_id": line[header["GENEID"]],
            "vert_phylop": line[header["VERPHYLOP"]],
            "annotation_type": line[header["ANNOTYPE"]],
            "consequence": line[header["CONSEQUENCE"]],
        }
        annotations.append(
            (feature, v, c, info)
        )

    return annotations

