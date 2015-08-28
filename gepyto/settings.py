
# Default settings.

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


import os
import re

try:
    import ConfigParser as configparser
except ImportError:
    import configparser


BUILD = ""
REFERENCE_PATH = ""
GEPYTO_ROOT = ""
DEBUG = False

CHROM_REGEX = re.compile(r"([0-9]{1,2}|MT|X|Y)")


def get_config():
    """Return a path to the gepyto configuration file."""
    # Check if the configuration file exists.
    config_file = os.path.join(GEPYTO_ROOT, "gepytorc.ini")
    if not os.path.isfile(config_file):
        _create_default_config(config_file)

    return config_file


def _create_default_config(fn):
    """Create the default configuration file."""

    config = configparser.RawConfigParser()

    config.add_section("gepytoConfiguration")
    config.set("gepytoConfiguration", "DEBUG", False)
    config.set("gepytoConfiguration", "BUILD", "GRCh37")
    config.set("gepytoConfiguration", "REFERENCE_PATH", "")

    config.add_section("databases")
    config.set("databases", "CCDS_PATH", "")

    with open(fn, "w") as f:
        config.write(f)


def _init_reference(config):
    global BUILD
    global REFERENCE_PATH

    BUILD = config.get("gepytoConfiguration", "BUILD")
    REFERENCE_PATH = config.get("gepytoConfiguration", "REFERENCE_PATH")

    # The REFERENCE_PATH can also be set as an environment variable.
    if REFERENCE_PATH == "" and os.environ.get("REFERENCE_PATH"):
        REFERENCE_PATH = os.environ.get("REFERENCE_PATH")


def _init_settings():
    # Create the directory where the configuration file will be.
    global DEBUG
    global GEPYTO_ROOT
    try:
        GEPYTO_ROOT = os.path.expanduser(os.environ["GEPYTO_ROOT"])
    except KeyError:
        GEPYTO_ROOT = None

    if GEPYTO_ROOT is None:
        GEPYTO_ROOT = os.path.abspath(os.path.join(
            os.path.expanduser("~"),
            ".gepyto"
        ))


    if not os.path.isdir(GEPYTO_ROOT):
        os.mkdir(GEPYTO_ROOT)

    # Read the configuration file.
    config = configparser.RawConfigParser()
    config.read(get_config())

    # Load settings related to the reference genome.
    _init_reference(config)

    # Check if debug mode.
    DEBUG = config.get("gepytoConfiguration", "DEBUG").lower() == "true"


_init_settings()
