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
import functools
import contextlib
import multiprocessing
import logging
logger = logging.getLogger(__name__)

import six
from six.moves.urllib.request import urlretrieve, urlopen
from six.moves import input, configparser

from Bio import bgzf  # Dependency for BioPython.

from .. import settings
from ..utils.progress import Progress


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
        uncompressed = zlib.decompress(data, zlib.MAX_WBITS|32)
        with bgzf.BgzfWriter(destination) as out:
            for line in uncompressed.splitlines():
                if line.startswith(">"):
                    line = line.split("|")[0]

                out.write(line + "\n")


def _restore_state():
    pass


_restore_state()
