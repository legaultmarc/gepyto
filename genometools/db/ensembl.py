
# Utilities to interact with the Ensembl database.
# This module contains code strictly to interact with Ensembl and not to 
# interpret the response or the underlying biology. 

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


import contextlib
import json
import logging
import time

try:
    # Python 2 support
    from urllib2 import urlopen, HTTPError
except ImportError:
    # Python 3 support
    from urllib.request import urlopen, HTTPError


__all__ = ["query_ensembl", ]

LAST_QUERY = 0

def query_ensembl(url):
    """Query the given (Ensembl rest api) url and get a json reponse.

    :param url: The API url to query.
    :type url: str

    :returns: A python object loaded from the JSON response from the server.

    """
    global LAST_QUERY

    this_query_time = time.time()
    delta_t = this_query_time - LAST_QUERY
    LAST_QUERY = this_query_time

    try:
        with contextlib.closing(urlopen(url)) as stream:
            response = json.loads(stream.read().decode())
            response_info = stream.info()

            limit = int(response_info["X-RateLimit-Limit"]) # Allowed / h
            reset = int(response_info["X-RateLimit-Reset"]) # Time to reset
            period = int(response_info["X-RateLimit-Period"])
            remaining = int(response_info["X-RateLimit-Remaining"])

            # Max time for request (s / request) to not exceed quota:
            max_t = 1.0 * reset / remaining
            if delta_t < max_t:
                time.sleep(max_t - delta_t + 0.5) # We add a buffer of 0.5s.

    except HTTPError as e:
        logging.warning("Request '{}' failed.".format(url)) 
        logging.warning("[{}] {}".format(e.code, e.reason)) 
        # If we busted we wait what they ask us to wait.
        if e.code == 429:
            sleep_time = float(e.info().getheader("Retry-After"))
            logging.warning("Waiting {}s before next Ensembl request (at "
                         "the server's request).".format(sleep_time))
            time.sleep(sleep_time)

            return query_ensembl(url)

        return None

    return response


def mysql_connect(ensembl_version, build):
    import MySQLdb

    if build == "GRCh37":
        build = "37"
    elif build == "GRCh38":
        build = "38"

    core_db_name = "homo_sapiens_core_{ensver}_{build}".format(
        ensver=ensembl_version,
        build=build,
    )


