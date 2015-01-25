#
# Utilities to interact with the UCSC database.
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


__all__ = ["UCSC", ]


from .. import settings


class UCSC(object):
    """Provides raw access to the UCSC MySQL database.

    The database will be set to the `db` parameter or to the `BUILD` as defined
    in `gepyto`'s settings.

    Later versions could implement features like throttling, but for now this
    is a very simple interface.

    """
    def __init__(self, db=None):
        import pymysql

        if db is None:
            db = settings.BUILD

        if db.lower() == "grch37":
            db = "hg19"
        elif db.lower() == "grch38":
            db = "hg38"
        else:
            raise Exception("Invalid genome reference '{}'".format(db))

        self.con = pymysql.connect(user="genome",
                                   host="genome-mysql.cse.ucsc.edu",
                                   database=db)

        self.cur = self.con.cursor()

    def raw_sql(self, sql, params):
        """Execute a raw SQL query."""
        self.cur.execute(sql, params)
        return self.cur.fetchall()

    def close(self):
        self.con.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
