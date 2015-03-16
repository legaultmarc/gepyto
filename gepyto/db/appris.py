# Module to get principal transcripts from the appris database.
# appris annotates transcripts and provides a main isoform for proteins.
# http://appris.bioinfo.cnio.es/
# The data included in this module was fetched from the APPRIS website on
# 2014-10-30.
# URL: http://appris.bioinfo.cnio.es/download/data/homo_sapiens/appris_data.principal.txt

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


import gzip
import collections
import sqlite3

from pkg_resources import resource_filename


APPRIS_CUR = None


def _load_appris():
    """Load the APPRIS database in memory as a sqlite3 database.

    :returns: The database cursor.

    """
    con = sqlite3.connect(":memory:")
    cur = con.cursor()

    cur.execute(
        "CREATE TABLE appris ("
        "  symbol TEXT,"
        "  ensembl_gene TEXT,"
        "  ensembl_transcript TEXT,"
        "  ccds TEXT,"
        "  category TEXT"
        ")"
    )

    fn = resource_filename(__name__, "data/appris_data.principal.txt.gz")
    db = collections.defaultdict(list)
    with gzip.open(fn) as f:
        for line in f:
            line = line.rstrip().decode("UTF-8")
            tu = tuple(line.split("\t"))
            cur.execute("INSERT INTO appris VALUES (?, ?, ?, ?, ?)", tu)

    con.commit()
    return cur


def init_db():
    """This is an initialization method for the database.

    We use this to load the database only if a function is called.

    """
    global APPRIS_CUR

    if APPRIS_CUR is None:
        APPRIS_CUR = _load_appris()


def get_transcripts_for_gene(ensg):
    """Fetches the transcripts and their annotation for a given gene (ENSG).

    :param ensg: The Ensembl gene id.
    :type ensg: str

    :returns: A list of transcript IDs and their categories (tuples).
    :rtype: tuple

    """

    init_db()
    APPRIS_CUR.execute(
        "SELECT ensembl_transcript, category FROM appris WHERE ensembl_gene=?",
        (ensg, )
    )
    return APPRIS_CUR.fetchall()


def get_category_for_transcript(enst):
    """Fetches the annotation for a transcript (ENST).

    :param enst: The Ensembl transcript id.
    :type enst: str

    :returns: The APPRIS category for this transcript.
    :rtype: str

    """

    init_db()
    APPRIS_CUR.execute(
        "SELECT category FROM appris WHERE ensembl_transcript=?",
        (enst, )
    )
    return APPRIS_CUR.fetchone()[0]


def get_main_transcripts(ensg):
    """Gets the main Ensembl transcript id for the provided gene based on the
       APPRIS annotation.

    :param ensg: The Ensembl gene number (ENSG000000).
    :param ensg: str

    :returns: The "main" transcrit (ENST). If there is an `appris_principal`
              annotation, this will be returned. If it is not the case, the
              order of priority is the following:
              `appris_candidate_longest_ccds`, `appris_candidate_ccds`,
              `appris_candidate_longest_seq`, `appris_candidate`.
    :rtype: str

    """

    init_db()
    APPRIS_CUR.execute(
        "SELECT ensembl_transcript, category FROM appris WHERE ensembl_gene=?",
        (ensg, )
    )
    li = APPRIS_CUR.fetchall()

    top_category = None
    categories = set([tu[1] for tu in li])

    ordered_cats = (
        "appris_principal", "appris_candidate_longest_ccds",
        "appris_candidate_ccds", "appris_candidate_longest_seq",
        "appris_candidate"
    )

    for cat in ordered_cats:
        if cat in categories:
            top_category = cat
            break

    top_transcripts = []
    for enst, cat in li:
        if cat == top_category:
            top_transcripts.append(enst)

    return enst
