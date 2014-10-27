
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
        with contextlib.closing(urllib2.urlopen(url)) as stream:
            response = json.load(stream)
            response_info = stream.info()

            limit = int(response_info.getheader("X-RateLimit-Limit")) # Allowed / h
            reset = int(response_info.getheader("X-RateLimit-Reset")) # Time to reset
            period = int(response_info.getheader("X-RateLimit-Period"))
            remaining = int(response_info.getheader("X-RateLimit-Remaining"))

            # Max time for request (s / request) to not exceed quota:
            max_t = 1.0 * reset / remaining
            if delta_t < max_t:
                time.sleep(max_t - delta_t + 0.5) # We add a buffer of 0.5s.

    except urllib2.HTTPError as e:
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

