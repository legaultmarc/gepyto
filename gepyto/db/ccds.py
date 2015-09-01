"""
Work with the CCDS database containing canonical isoforms of proteins along
with their coding sequence.
Because of the large size of this database, using this module for the first
time will require an installation step.
"""

# This file is part of gepyto.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.

from __future__ import division

__author__ = "Marc-Andre Legault"
__copyright__ = ("Copyright 2014 Marc-Andre Legault and Louis-Philippe "
                 "Lemieux Perreault. All rights reserved.")
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


import os
import zlib
import operator
import sqlite3
import functools
import contextlib
import multiprocessing
import collections
import logging
logger = logging.getLogger(__name__)

import six
from six.moves.urllib.request import urlretrieve, urlopen
from six.moves import input, configparser

from Bio import bgzf  # Dependency for BioPython.

from .. import settings
from ..reference import GenericReference
from ..utils.progress import Progress
from ..structures.region import Region


ccds_fields = ("chromosome", "nc_accession", "gene", "gene_id", "ccds_id",
               "ccds_status", "cds_strand", "cds_from", "cds_to",
               "cds_locations", "match_type")

indexed_fields = set(("nc_accesssion", "gene", "gene_id", "ccds_id"))

CCDSRecord = collections.namedtuple("CCDSRecord", ccds_fields)


class CCDS(object):
    def __init__(self, install=False):
        # Try to detect the installation parameters.
        state = _load_installation()
        if not state:
            logger.info("Could not find a CCDS installation. Use install=True "
                        "if you want to install a local copy.")
            if install:
                local_install(confirm=False)
            else:
                raise Exception("CCDS not installed.")

        self.build_info = state["build_info"]
        self.con = state["index"][0]
        self.cur = state["index"][1]
        self.fasta = state["fastas"]
        self.ccds = open(state["ccds_filename"], "r")

    @staticmethod
    def _protein_to_ccds(start, end):
        """Convert the amino acids coordinates to CCDS.

        :param start: the protein starting position (protein coordinates)
        :param end: the protein ending position (protein coordinates)

        :type start: int
        :type end: int

        :returns: the same region, but in CCDS coordinates.
        :rtype: tuple

        """
        return start * 3, end * 3 + 2

    def _ccds_to_genomic(self, record, start, end):
        """converts the CCDS coordinates to genomic.

        :param record: the record for this CCDS.
        :param start: the CCDS starting position (CCDS coordinates)
        :param end: the CCDS ending position (CCDS coordinates)

        :type record: CCDSRecord
        :type start: int
        :type end: int

        :returns: The same region, but in genomic coordinates. This region
                  might be a non contiguous region.
        :rtype: gepyto.structures.region.Region

        """
        # The chromosome
        chrom = record.chromosome

        # Creating a region for the asked CDS segment (CDS coordinate)
        cds_region = Region(chrom=chrom, start=start, end=end)

        # The final region
        final_regions = []

        # Checking for the strand (if '-', we need to cycle from the end)
        cds_locations = sorted(record.cds_locations, key=lambda x: x[0])
        if record.cds_strand == "-":
            cds_locations = cds_locations[::-1]

        # The first segment starting position in CCDS coordinates
        cds_start = 1

        # Cycling through the CDS segment (genomic coordinates)
        for genomic_start, genomic_end in cds_locations:
            # Getting the current segment ending position (CDS coordinates)
            cds_end = cds_start + genomic_end - genomic_start

            # Creating a region for this specif region (CDS coordinates)
            curr_cds_region = Region(chrom=chrom, start=cds_start, end=cds_end)

            if curr_cds_region.start > cds_region.end:
                # We're pass the segment
                break

            if not cds_region.overlaps_with(curr_cds_region):
                # The regions are not overlapping
                continue

            # Getting the overlapping region of the two CDS region
            overlap = cds_region.intersect(curr_cds_region)

            # Computing the CDS region (genomic coordinates)
            cds_genomic_start = genomic_start + overlap.start - cds_start
            cds_genomic_end = genomic_start + overlap.end - cds_start

            if record.cds_strand == "-":
                # The correction is different if this is the first segment
                if len(final_regions) == 0:
                    correction = -(cds_genomic_start - genomic_start)
                else:
                    correction = genomic_end - cds_genomic_end

                cds_genomic_start += correction
                cds_genomic_end += correction

            # Adding this genomic region
            final_regions.append(Region(
                chrom=chrom,
                start=cds_genomic_start,
                end=cds_genomic_end,
            ))

            # Getting the next segment staring position (CDS coordinates)
            cds_start = cds_end + 1

        return Region.from_regions(final_regions)

    def protein_to_genomic(self, record, start, end):
        """Convert the amino acids coordinates to genomic.

        """
        # TODO IMPLEMENT ME
        raise NotImplementedError()

    def get_record(self, **kwargs):
        """Get a record using the col=value syntax."""
        if len(kwargs) != 1:
            raise TypeError("You need to provide a single col=value to query "
                            "the database.")

        field, value = list(kwargs.items())[0]
        if field not in indexed_fields:
            msg = "Can't query on '{}'. Valid columns are: '{}'".format(
                field, indexed_fields
            )
            raise TypeError(msg)
        sql = "SELECT * FROM ccds_index WHERE {}=?;".format(field)
        self.cur.execute(sql, (value, ))
        results = []
        for tu in self.cur:
            tell = tu[-1]
            self.ccds.seek(tell)
            line = self.ccds.readline().rstrip().split("\t")

            # Parse the positions into integers.
            positions = []
            for pos in line[-2].lstrip("[").rstrip("]").split(","):
                start, end = [int(i) for i in pos.split("-")]
                positions.append((start, end))

            line[-2] = positions

            results.append(
                CCDSRecord(*line)
            )

        return results

    def close(self):
        self.con.close()
        self.ccds.close()


def local_install(directory=None, build=settings.BUILD, confirm=True):
    if directory is None:
        # Check if user specified a file in the configuration file.
        config_filename = settings.get_config()
        config = configparser.RawConfigParser()
        config.read(config_filename)

        try:
            directory = config.get("databases", "CCDS_PATH")
        except Exception:
            directory = ""

        if not directory:
            directory = ""
    else:
        logger.warning("You specified the path to an existing CCDS. "
                       "Note that gepyto will not remember this path and that "
                       "you should specify it explicitly every time you use "
                       "this module. If you want gepyto to remember it, edit "
                       "the configuration file available at: '{}'".format(
                           settings.get_config()
                       ))

    if directory == "":
        # If the user didn't specify anything, we use the default, which
        # is a GEPYTO_ROOT/CCDS.
        directory = os.path.join(settings.GEPYTO_ROOT, "CCDS")
        if not config.has_section("databases"):
            config.add_section("databases")

        config.set("databases", "CCDS_PATH", directory)
        with open(config_filename, "w") as f:
            config.write(f)

        os.makedirs(directory)

    build = build.lower()
    if build not in ("grch37", "grch38"):
        raise ValueError("Build {} not recognized. Compatible builds are: "
                         "GRCh37 and GRCh38.".format(build))

    # Get a CCDS version compatible with the build.
    if build == "grch37":
        base = "ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/15"
        urls = {
            "build_info": os.path.join(base, "BuildInfo.20131129.txt"),
            "ccds": os.path.join(base, "CCDS.20131129.txt"),
            "nucleotides": os.path.join(
                base, "CCDS_nucleotide.20131129.fna.gz"
            ),
            "proteins": os.path.join(base, "CCDS_protein.20131129.faa.gz")
        }
        version = "15"

    elif build == "grch38":
        base = "ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/18"
        urls = {
            "build_info": os.path.join(base, "BuildInfo.20150512.txt"),
            "ccds": os.path.join(base, "CCDS.20150512.txt"),
            "nucleotides": os.path.join(
                base, "CCDS_nucleotide.20150512.fna.gz"
            ),
            "proteins": os.path.join(base, "CCDS_protein.20150512.faa.gz")
        }
        version = "18"

    logger.info("Installing CCDS version '{}' for build '{}' in '{}'.".format(
        version, build, directory
    ))

    _download_all(urls, version, directory, confirm=confirm)


def _reporthook(num_blocks, chunk_size, total_size, bar):
    if bar["bar"] is None:
        bar["bar"] = Progress(total_size)

    if num_blocks * chunk_size < total_size:
        bar["bar"].update(num_blocks * chunk_size)
    else:
        bar["bar"].update(total_size)


def _download_all(urls, version, directory, confirm=True):
    """Download all specified files to the directory."""
    if confirm:
        msg = ("A copy of the CCDS database will be downloaded for a total "
               "of ~35MB at {}.\nIs this ok? [y/N] ").format(directory)
        confirm = input(msg)
        if confirm != "y":
            logger.warning("CCDS download cancelled by user.")
            return

    path = os.path.join(directory, version)
    os.makedirs(path)
    procs = []
    for key in ("nucleotides", "proteins"):
        url = urls[key]
        logger.info(
            "Starting a download and BGZF re-encoding for {}".format(url)
        )
        p = multiprocessing.Process(
            target=_download_recompress,
            args=(url, os.path.join(path, key) + ".fa.gz")
        )
        procs.append(p)
        p.start()

    for key in ("ccds", "build_info"):
        url = urls[key]
        # Regular (unzipped file).
        logger.info("Downloading: {}".format(url))
        bar = {"bar": None}
        hook = functools.partial(_reporthook, bar=bar)
        urlretrieve(
            url,
            os.path.join(path, key) + ".txt",
            reporthook=hook
        )
        bar["bar"].finish()

    p.join()


def _download_recompress(url, destination):
    """Download a gzipped file and recompress using bgzf."""
    # Start a download in a subprocess and convert from gzip to bgzf before
    # writing to disk.

    # The contextlib is needed for PY2.
    with contextlib.closing(urlopen(url)) as f:
        data = f.read()
        uncompressed = zlib.decompress(data, zlib.MAX_WBITS | 32)
        with bgzf.BgzfWriter(destination) as out:
                out.write(uncompressed)


def _load_installation():
    config_filename = settings.get_config()
    config = configparser.RawConfigParser()
    config.read(config_filename)

    try:
        path = config.get("databases", "CCDS_PATH")
    except Exception:
        return

    version = os.listdir(path)
    if len(version) == 0:
        # No local CCDS installation.
        return
    elif len(version) > 1:
        logger.warning("Multiple CCDS versions are installed. Restore the one "
                       "you wish to use manually.")
        return

    # There is only one installed version so we can restore it safely.
    version = version[0]
    path = os.path.join(path, version)

    # Parse the build info.
    try:
        with open(os.path.join(path, "build_info.txt"), "r") as f:
            build_info = [i.split() for i in f.readlines()]
            build_info[0][0] = build_info[0][0].lstrip("#")
            build_info = dict(zip(*build_info))
    except IOError:
        logger.debug("Could not parse the build_info.txt file.")
        return

    # Build an index in memory.
    ccds_filename = os.path.join(path, "ccds.txt")
    try:
        with open(ccds_filename, "r") as f:
            con, cur = _get_ccds_index(f)
    except IOError:
        logger.debug("Could not parse the ccds.txt file.")
        return

    fastas = {
        "proteins": os.path.join(path, "proteins.fa.gz"),
        "nucleotides": os.path.join(path, "nucleotides.fa.gz"),
    }
    for k, filename in fastas.items():
        if not os.path.isfile(filename):
            logger.debug("Could not find the '{}' file.".format(filename))
            return

    return {
        "build_info": build_info,
        "index": (con, cur),
        "fastas": fastas,
        "ccds_filename": ccds_filename
    }


def _get_ccds_index(file_object):
    con = sqlite3.connect(":memory:")
    cur = con.cursor()

    # Create the table.
    cur.execute("""
        CREATE TABLE ccds_index (
            nc_accession TEXT,
            gene TEXT,
            gene_id INTEGER,
            ccds_id TEXT,
            tell INTEGER
        );
    """)

    header = file_object.readline().lstrip("#").split()

    for i, field in enumerate(header):
        assert field == ccds_fields[i]

    header = dict(zip(header, range(len(header))))

    _data = []
    tell = file_object.tell()
    line = file_object.readline()
    while line:
        line = line.rstrip().split()
        _data.append((
            line[header["nc_accession"]],
            line[header["gene"]],
            int(line[header["gene_id"]]),
            line[header["ccds_id"]],
            tell
        ))

        tell = file_object.tell()
        line = file_object.readline()

    cur.executemany("INSERT INTO ccds_index VALUES (?, ?, ?, ?, ?)", _data)
    return con, cur
