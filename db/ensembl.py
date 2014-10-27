
# Utilities to interact with the Ensembl database.
# This module contains code strictly to interact with Ensembl and not to 
# interpret the response or the underlying biology. 

import contextlib
import json
import urllib2
import logging
import time

__all__ = ["query_ensembl", ]

LAST_QUERY = 0

def query_ensembl(url):
    global LAST_QUERY

    this_query_time = time.time()
    delta_t = this_query_time - LAST_QUERY
    LAST_QUERY = this_query_time

    try:
        with contextlib.closing(urllib2.urlopen(url)) as stream:
            response = json.load(stream)
            response_info = stream.info()

            limit = int(response_info.getheader("X-RateLimit-Limit")) # Allowed / h
            reset = int(response_info.getheader("X-RateLimit-Reset")) # Time to reset
            period = int(response_info.getheader("X-RateLimit-Period"))
            remaining = int(response_info.getheader("X-RateLimit-Remaining"))

            # Max rate (requests/s) to not exceed quota:
            max_rate = 1.0 * remaining / reset
            cur_rate = 2.0 / delta_t
            if cur_rate > max_rate:
                sleep(0.5)

    except urllib2.HTTPError as e:
        logging.warning("Request '{}' failed.".format(url)) 
        logging.warning("[{}] {}".format(e.code, e.reason)) 
        return None

    return response

