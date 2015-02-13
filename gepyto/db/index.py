# Methods build generic indices for text files containing genomic coordinates.
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

import atexit
import sqlite3
import logging
import re
import os
import bisect

import numpy as np


_registered_connexions = []
_current_ckey = 1


def _create_index_db(fn):
    if os.path.isfile(fn):
        os.remove(fn)

    con = sqlite3.connect(fn)
    cur = con.cursor()

    # The meta table contains information on the file.
    cur.execute("CREATE TABLE meta (delim TEXT, chrom INTEGER, pos INTEGER);")

    # The key index contains large values aimed at encoding chromosomes.
    cur.execute("CREATE TABLE key ("
                "  chrom TEXT NOT NULL,"
                "  ckey INTEGER NOT NULL,"
                "  PRIMARY KEY (chrom),"
                "  UNIQUE (ckey)"
                ");")

    # The index table contains the actual index.
    cur.execute("CREATE TABLE idx ("
                "  code INTEGER NOT NULL,"
                "  seek INTEGER NOT NULL,"
                "  PRIMARY KEY (code)"
                ");")

    con.commit()

    return con, cur


def _db_insert(cur, chrom, pos, tell):
    global _current_ckey

    # Get the code for the chromosome.
    cur.execute("SELECT ckey FROM key WHERE chrom=?", (chrom, ))

    try:
        code = cur.fetchone()[0]
    except TypeError:
        # We need to insert the code for this chromosome.
        code = _current_ckey * 10 ** 10
        cur.execute(
            "INSERT INTO key VALUES (?, ?)",
            (chrom, code)
        )

        # Increment the chromosome key for the next one.
        _current_ckey += 1

    # Now do the insert.
    code += pos
    cur.execute("INSERT INTO idx VALUES (?, ?)", (code, tell))


def build_index(fn, chrom_col, pos_col, delimiter='\t', skip_lines=0,
                index_rate=0.2, ignore_startswith=None):
    """Build a index for the given file.

    :param fn: The filename
    :type fn: str

    :param chrom_col: The column representing the chromosome (0 based).
    :type chrom_col: int

    :param pos_col: The column for the position on the chromosome (0 based).
    :type pos_col: int

    :param delimiter: The delimiter for the columns (default tab).
    :type delimiter: str

    :param skip_lines: Number of header lines to skip.
    :type skip_lines: int

    :param index_rate: The approximate rate of line indexing. As an example,
                       a file with 1000 lines and the default index_rate of
                       0.2 will have an index with ~200 entries.
    :type index_rate: float

    :param ignore_startswith: Ignore lines that start with a given string.
                              This can be used to skip headers, but will not
                              be used to parse the rest of the file.
    :type ignore_startswith: str

    :returns: The index filename.
    :rtype: str

    """

    idx_fn = _get_index_fn(fn)

    con, cur = _create_index_db(idx_fn)

    # Add information about the indexed file in the db.
    cur.execute(
        "INSERT INTO meta VALUES (?, ?, ?)",
        (delimiter, chrom_col, pos_col)
    )
    con.commit()

    size = os.path.getsize(fn)  # Total filesize
    start = 0  # Byte position of the meat of the file (data).
    with open(fn, "r") as f:
        if skip_lines > 0:
            for i in range(skip_lines):
                # Skip header lines if needed.
                _ = f.readline()

        current_position = f.tell()
        if ignore_startswith is not None:
            line = f.readline()
            while line.startswith(ignore_startswith):
                current_position = f.tell()
                line = f.readline()
            f.seek(current_position)

        start = f.tell()

        size -= start  # Adjust file size to remove header.

        # Estimate the line length using first 100 lines
        line_length = np.empty((100))
        prev = start
        for i in range(100):
            _ = f.readline()
            line_length[i] = f.tell() - prev
            prev = f.tell()

        line_length = np.mean(line_length)
        approx_num_lines = size / line_length

        # Go back to start
        f.seek(start)

        # Add the first line to the index.
        l1 = f.readline().rstrip().split(delimiter)
        chrom1, pos1 = (l1[chrom_col], l1[pos_col])
        chrom1 = chrom1.lstrip("chr")
        pos1 = int(pos1)
        _db_insert(cur, chrom1, pos1, start)
        f.seek(start)

        # Compute the seek jump size.
        target_num_lines = index_rate * approx_num_lines
        seek_jump = size / target_num_lines

        # Now we will jump of seek_jump, get to a fresh line (because we will
        # probably land in the middle of a line).
        current_position = f.tell() + seek_jump

        # We need to remember the previous position to make sure the file is
        # sorted.
        # Tuple of chromosome, position.
        prev = (-1, -1)

        while current_position + seek_jump < size:
            # Jump in the file.
            f.seek(current_position)
            # Throw away partial line.
            f.readline()

            tell = f.tell()
            line = f.readline().rstrip().split(delimiter)
            chrom, pos = (line[chrom_col], line[pos_col])
            if chrom.startswith("chr"):
                chrom = chrom[3:]

            pos = int(pos)

            if prev[0] == chrom and prev[1] == pos:
                # We need to keep searching.
                current_position += seek_jump
                continue

            try:
                _db_insert(cur, chrom, pos, tell)
            except sqlite3.IntegrityError:
                con.close()
                raise Exception("File is not sorted.")

            prev = (chrom, pos)

    # Commit the index.
    con.commit()
    con.close()

    return idx_fn


def get_index(fn):
    """Restores the index for a given file or builds it if the index was not
       previously created.

    :param fn: The filname of the file to index.
    :type fn: str

    :returns: The index dictionary corresponding to the input file.
    :rtype: dict

    """

    global _registered_connexions

    con = sqlite3.connect(_get_index_fn(fn))
    cur = con.cursor()

    _registered_connexions.append(con)

    return cur


def goto(f, cursor, chrom, pos):
    """Given a file and the cursor to an index, go to the genomic coordinates.

    :param f: An open file.
    :type f: file

    :param cursor: A sqlite3 cursor to the index database.
    :param idx: :py:class:`sqlite3.Cursor`

    :param chrom: The queried chromosome.
    :param pos: The queried position on the chromosome.

    :returns: True if the position was found and the cursor moved, False if
              the queried chromosome, position wasn't found.
    :rtype: bool

    """

    # Encode the target.
    #   - Find the chromosome code.
    cursor.execute("SELECT ckey FROM key WHERE chrom=?", (chrom, ))
    try:
        ckey = cursor.fetchone()[0]
    except TypeError:
        # Couldn't find chromosome key. This means it is not indexed.
        return Exception("Chromosome '{}' not found in the index.".format(
            chrom, 
        ))

    #   - Encode the position
    target = ckey + pos

    # Find the two nearest index anchors.
    sql = ("SELECT code, seek, abs(code - {target}) AS abs_distance,"
           " code - {target} AS distance "
           "FROM idx ORDER BY abs_distance ASC LIMIT 2")
    sql = sql.format(target=target)

    cursor.execute(sql)
    hits = cursor.fetchall()

    assert len(hits) == 2, "Internal error in the index."

    # If we have a sandwich hit, we will have one positive and one negative
    # distance.
    if hits[0][3] != hits[1][3]:
        top, bottom = sorted(hits, key=lambda x: x[1])
        return goto_fine(f, chrom, pos, top[1], bottom[1])

    # We didn't get a sandwich hit. Check if they are both after the hit. If 
    # this is the case, we will search from the top of the file.
    if hits[0][3] + hits[0][3] >= 0:
        return goto_fine(f, chrom, pos, 0, min(hits[0][1], hits[1][1]))

    # If this wasn't before everything, maybe it is after everything.
    cur.execute("SELECT max(seek) FROM idx")
    max_seek = cur.fetchone()[0]
    return goto_fine(f, chrom, pos, max_seek, float("+infinity"))

    # If we're still here then I don't know what happened.
    raise Exception("Couldn't use index, no anchors found.")


def goto_fine(f, chrom, pos, top, bottom):
    # TODO implement binary search.

    # Get the meta information.
    idx_filename = _get_index_fn(f.name)
    con = sqlite3.connect(idx_filename)
    cursor = con.cursor()
    try:
        delim, chrom_col, pos_col = cursor.execute(
            "SELECT * FROM meta"
        ).fetchone()
    finally:
        con.close()

    # Go to the start of the plausible region.
    f.seek(top)

    while f.tell() < bottom:
        # Iterate through the file.
        here = f.tell()
        line = f.readline().rstrip().split(delim)

        if len(line) == 1:  # This is 1 because of the split of ""
            # End of file.
            return False

        this_chrom = line[chrom_col]
        this_chrom = this_chrom.lstrip("chr")

        this_pos = int(line[pos_col])

        # We got it!
        if this_chrom == chrom and this_pos == pos:
            f.seek(here)
            return True

    return False


def _get_index_fn(fn):
    """Generates the index filename from the path to the indexed file.

    :param fn: The name of the file to index.
    :type fn: str

    """

    if not os.path.isfile(fn):
        raise Exception("File '{}' does not exist.".format(
            fn
        ))

    return os.path.abspath("{}.gtidx".format(fn))


def close_connexions():
    for con in _registered_connexions:
        con.close()


atexit.register(close_connexions)
