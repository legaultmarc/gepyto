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

import logging
import re
import os

import numpy as np

MAGIC_NUMBER = 10 ** 9  # This should be larger than any indexed position.


class EndOfFile(Exception):
    pass


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

    size = os.path.getsize(fn)  # Total filesize
    start = 0  # Byte position of the meat of the file (data).
    with open(fn, "r") as f:
        if skip_lines > 0:
            for i in range(skip_lines):
                # Skip header lines if needed.
                f.readline()

        current_position = f.tell()
        # Skip the lines that start with the user provided string.
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
        # We take 100 sample positions in the file.
        for i, jump in enumerate(np.linspace(0, 0.9 * size, 100)):
            f.seek(start + int(jump))
            f.readline()  # Throw away chunk.
            line_start = f.tell()
            f.readline()
            line_length[i] = f.tell() - line_start

        line_length = np.mean(line_length)
        approx_num_lines = size / line_length

        # Go back to start
        f.seek(start)

        # Prepare the variables to keep track of the position encoding.
        encoding = {
            "current_key": 0,
            "chromosomes": {}
        }
        index = []  # List of (code, seek) tuples.

        def encode_locus(chrom, pos):
            chrom = str(chrom)
            if chrom not in encoding["chromosomes"]:
                encoding["current_key"] += 1
                encoding["chromosomes"][chrom] = encoding["current_key"]
            
            return encoding["chromosomes"][chrom] * MAGIC_NUMBER + pos


        def get_locus(line):
            if not line:
                raise EndOfFile()

            line = line.rstrip().split(delimiter)
 
            chrom = line[chrom_col]
            if chrom.startswith("chr"):
                chrom = chrom[3:]
            pos = int(line[pos_col])

            return chrom, pos

        # Add the first line to the index.
        chrom1, pos1 = get_locus(f.readline())
        index.append((encode_locus(chrom1, pos1), start))

        # Compute the seek jump size.
        target_num_lines = index_rate * approx_num_lines
        seek_jump = size / target_num_lines

        # We start indexing here.
        current_position = f.tell()

        while current_position + seek_jump < size + start:
            # Jump in the file.
            f.seek(current_position + seek_jump)
            # Throw away partial line.
            f.readline()

            # We need to make sure this is a unique line.
            try:
                cur_chrom, cur_pos = get_locus(f.readline())

                tell = f.tell()
                next_chrom, next_pos = get_locus(f.readline())
                while next_chrom == cur_chrom and next_pos == cur_pos:
                    tell = f.tell()
                    next_chrom, next_pos = get_locus(f.readline())

                # We found "new" content.
                # First we make sure that the file looks sorted. Then, we
                # index it.
                code = encode_locus(next_chrom, next_pos)
                if index[-1][0] > code:
                    raise Exception("This file is not sorted.")

                index.append((code, tell))
                current_position = tell

            except EndOfFile:
                break  # Reached the end of the file.

    # Write the index pickle.
    index = np.array(index)
    np.save(idx_fn, index)

    return idx_fn


def get_index(fn):
    """Restores the index for a given file or builds it if the index was not
       previously created.

    :param fn: The filname of the file to index.
    :type fn: str

    :returns: The index dictionary corresponding to the input file.
    :rtype: dict

    """

    return np.load(fn + ".npy")


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
    sql = ("SELECT code, seek, min(abs(code - {target})), "
           "code - {target} as dist"
           " FROM idx GROUP BY dist<=0")
    sql = sql.format(target=target)

    cursor.execute(sql)
    hits = cursor.fetchall()

    # Handle the case where the position is directly in the index.
    direct_hit = False
    for hit in hits:
        if hit[2] == 0:
            f.seek(hit[1])
            return True

    # If we have a sandwich hit, we will have one positive and one negative
    # distance.
    if len(hits) == 2:
        top, bottom = hits if hits[0][3] <= 0 else hits[::-1]
        return goto_fine(f, chrom, pos, top[1], bottom[1])

    assert len(hits) == 1, "Invalid indexing."
    hits = hits[0]
    # We didn't get a sandwich hit. Check if we have a hit before or after the
    # nearest anchor.
    if hits[3] < 0:
        # We hit after the last index.
        return goto_fine(f, chrom, pos, hits[1], float("+infinity"))
    else:
        # We hit before the first index.
        # Because the first line is always indexed, we just return false.
        return False

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
