
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
import urllib2
import time
import uuid
import os
import gzip

import requests # Install http://docs.python-requests.org/en/latest/

__all__ = ["cadd", ]

def cadd(variants):
    """Annotate the variants using CADD (cadd.gs.washington.edu).

    :param vcf: A list of Variant (or subclass) objects.
    :type vcf: :py:class`genometools.structures.variants.Variant`

    :returns: A list of the variants with an added "c_score" attribute.
    :rtype: list

    """

    # We need to send a form with the following information:
    # file: a vcf file.
    # inclAnno=Yes
    # to http://cadd.gs.washington.edu/upload
    # enctype="multipart/form-data"
    # with POST

    # Then we parse the html to find a link to a anno_.+\.tsv\.gz
    # We check every 20s to see if we can download the file.

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
            fn = m.group(1)

            # Retry to download every 5 minutes.
            download_url = "http://cadd.gs.washington.edu/static/finished/{}"
            download_url = download_url.format(fn)
            success = False
            for i in xrange(24): # Try, wait for 5 minutes, retry.
                time.sleep(5 * 60)

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
                variants = annotate_variants(f, variants)
            os.remove(fn)

            return variants

        else:
            print m
            print feedback
            raise Exception("Could not find download link in CADD Response.")

    else:
        raise Exception("Request to CADD Failed (HTTP {}).".format(
            r.status_code 
        ))


def annotate_variants(cadd_output, variants):
    """Takes a file object containing the output from CADD and parses it to annotate the variants in the list.

    """

    for line in cadd_output:
        print line

