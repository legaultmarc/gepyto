# Module to create MongoDB databases containing Gepyto objects.
# This can be useful to get better query performance or for data persistence.
#
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


import functools
import json

try:
    import pymongo
    _pymongo_installed = True
except ImportError:
    _pymongo_installed = False


class GepytoDB(object):
    """Main class wrapping operations on the Gepyto Database."""
    def __init__(self, db_name, host="localhost", port=27017):
        if not _pymongo_installed:
            raise ImportError("Install pymongo to use this module.")

        self.client = pymongo.MongoClient(host, port)
        self.db = self.client[db_name]

    def insert(self, collection, obj):
        """Insert the gepyto object into the database."""
        return self.bulk_insert(collection, [obj])

    def bulk_insert(self, collection, li):
        """Insert a list of objects to the database."""
        if all([hasattr(i, "as_document") for i in li]):
            docs = [i.as_document() for i in li]

        elif all([type(i) is dict for i in li]):
            docs = li

        elif all([type(i) is str for i in li]):
            docs = [json.loads(i) for i in li]

        else:
            docs = [i.__dict__ for i in li]

        self.db[collection].insert(docs)

    def find_one(self, collection, *args):
        return self.db[collection].find_one(*args)

    def find(self, collection, *args):
        return self.db[collection].find(*args)
