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

try:
    import cPickle as pickle
except ImportError:
    import pickle

import bisect
import logging
import re
import os
import functools

import numpy as np


MAGIC_NUMBER = 10 ** 9  # This should be larger than any indexed position.


class EndOfFile(Exception):
    pass


class ChromosomeNotIndexed(Exception):
    pass


def _get_locus(line, chrom_col, pos_col, delimiter):
    if not line:
        raise EndOfFile()

    line = line.rstrip().split(delimiter)

    chrom = line[chrom_col]
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    pos = int(line[pos_col])

    return chrom, pos


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

    assert chrom_col != pos_col

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

        # Bind some parameters so that function calls look nicer.
        get_locus = functools.partial(
            _get_locus,
            chrom_col=chrom_col,
            pos_col=pos_col,
            delimiter=delimiter
        )

        # Add the first line to the index.
        chrom1, pos1 = get_locus(f.readline())
        index.append((encode_locus(chrom1, pos1), start))

        # Compute the seek jump size.
        target_num_lines = index_rate * approx_num_lines
        seek_jump = size / target_num_lines

        # We start indexing here.
        current_position = f.tell()

        if index_rate == 1:
            logging.debug("Full indexing mode.")
            tell = f.tell()
            line = f.readline()
            while line:
                chrom, pos = get_locus(line)
                code = encode_locus(chrom, pos)

                if index[-1][0] > code:
                    raise Exception("This file is not sorted.")
                elif index[-1][0] == code:
                    pass
                else:
                    index.append((code, tell))

                tell = f.tell()
                line = f.readline()

        else:
            logging.debug("Sparse indexing mode.")
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

    index = np.array(index)

    # Create a dict containing the relevant information to be able to find the
    # chromosome and position columns.
    info = {"chrom_col": chrom_col, "pos_col": pos_col, "delimiter": delimiter,
            "chrom_codes": encoding["chromosomes"]}

    pickle_string = pickle.dumps(info)

    # Write the index pickle (and numpy matrix).
    with open(idx_fn, "wb") as f:
        f.write(pickle_string)
        np.save(f, index)

    return idx_fn


def get_index(fn):
    """Restores the index for a given file or builds it if the index was not
       previously created.

    :param fn: The filname of the file to index.
    :type fn: str

    :returns: The numpy array representing the actual index.
    :rtype: :py:class:`numpy.ndarray`

    """

    # Read the information pickle part.
    # We use the numpy format definition to know when to stop:
    # https://github.com/numpy/numpy/blob/master/doc/neps/npy-format.rst
    indexed_filename = fn
    fn = _get_index_fn(indexed_filename)
    with open(fn, "rb") as f:
        i = 0
        chunk = None
        while chunk != b"\x93NUMPY":
            f.seek(i)
            chunk = f.read(6)
            if not chunk:
                raise Exception("Invalid format for the index.")
            i += 1

        pickle_length = f.tell() - 6
        f.seek(0)
        info = pickle.loads(f.read(pickle_length))

        index = np.load(f)

    # This is a class that will behave like the (info, index) tuple.
    # It is only used to make it printing the index object prettier.
    class Index(object):
        def __init__(self, info, index, filename):
            self.filename = filename
            self.info = info
            self.index = index

        def __getitem__(self, key):
            if key == 0:
                return self.info
            elif key == 1:
                return self.index
            else:
                raise IndexError()

        def __iter__(self):
            return (i for i in (self.info, self.index))

        def __repr__(self):
            return "<{} object for file '{}'>".format(self.__class__.__name__,
                                                      indexed_filename)

    return Index(info, index, fn)


def goto(f, index, chrom, pos):
    """Given a file, a locus and the index, go to the genomic coordinates.

    :param f: An open file.
    :type f: file

    :param index: This is actually a tuple. The first element is an information
                  dict containing the delimiter, chromosome column and position
                  column. The second element is a numpy matrix containing the
                  encoded loci and the "tell" positions.
    :param index: tuple

    :param chrom: The queried chromosome.
    :param pos: The queried position on the chromosome.

    :returns: True if the position was found and the cursor moved, False if
              the queried chromosome, position wasn't found.
    :rtype: bool

    """

    # Type checks for the parameters.
    chrom = str(chrom)
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    pos = int(pos)

    info, index = index

    # Encode the locus.
    chrom_code = info["chrom_codes"].get(str(chrom))
    if chrom_code is None:
        raise ChromosomeNotIndexed(
            "Chromosome '{}' is not in the index.".format(chrom)
        )

    code = chrom_code * MAGIC_NUMBER + pos
    logging.debug("Looking for code: {}".format(code))

    # Find the boundaries for the locus.
    boundary = bisect.bisect_left(index[:, 0], code)

    # The hit is at the end of the index.
    if boundary == index.shape[0]:
        logging.debug("Looking from the end of the index.")
        return goto_fine(f, chrom, pos, index[-1][1], float("+infinity"), info)

    if index[boundary, 0] == code:
        # We got a direct hit.
        logging.debug("Direct hit on locus.")
        f.seek(index[boundary, 1])
        return True

    # We always index the first position, so nothing should come before unless
    # it's a direct hit.
    if boundary == 0:
        logging.debug("Locus before first index.")

    left = index[boundary - 1, 1]
    left_code = index[boundary, 0]
    try:
        right = index[boundary + 1, 1]
        right_code = index[boundary + 1, 0]
    except IndexError:
        right = float("+infinity")  # Right boundary not in index, use the EOF.
        right_code = float("+infinity")

    logging.debug("Found boundaries: {}:{} and {}:{}. Using linear search to "
                  "find the exact position.".format(left_code, left,
                                                    right_code, right))
    return goto_fine(f, chrom, pos, left, right, info)


def goto_fine(f, chrom, pos, left, right, info):

    # Go to the start of the plausible region.
    f.seek(left)

    # Bind some parameters to the line parser.
    get_locus = functools.partial(
        _get_locus,
        chrom_col=info["chrom_col"],
        pos_col=info["pos_col"],
        delimiter=info["delimiter"],
    )

    # Iterate through the file.
    while True:
        here = f.tell()

        try:
            this_chrom, this_pos = get_locus(f.readline())

            # We got it!
            if this_chrom == chrom and this_pos == pos:
                f.seek(here)
                return True

        except EndOfFile:
            logging.debug("Reached end of file without finding the locus.")
            return False

        if here > right:
            logging.debug("Didn't find the locus before the right boundary.")
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
