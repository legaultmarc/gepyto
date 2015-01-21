
# Default settings.

# This file is part of genometools.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


__author__ = "Marc-Andre Legault"
__copyright__ = ("Copyright 2014 Marc-Andre Legault and Louis-Philippe Lemieux "
                 "Perreault. All rights reserved.")
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


import os

try:
    import ConfigParser as configparser
except ImportError:
    import configparser


BUILD = ""
REFERENCE_PATH = ""

def _create_default_config(fn):
    # Create the default configuration file.

    config = configparser.RawConfigParser()

    config.add_section("GenometoolsConfiguration")
    config.set("GenometoolsConfiguration", "BUILD", "GRCh37")
    config.set("GenometoolsConfiguration", "REFERENCE_PATH", "")

    with open(fn, "w") as f:
        config.write(f)


def _init_reference(config):
    global BUILD
    global REFERENCE_PATH

    BUILD = config.get("GenometoolsConfiguration", "BUILD")
    REFERENCE_PATH = config.get("GenometoolsConfiguration", "REFERENCE_PATH")

    # The REFERENCE_PATH can also be set as an environment variable.
    if REFERENCE_PATH == "" and os.environ.get("REFERENCE_PATH"):
        REFERENCE_PATH = os.environ.get("REFERENCE_PATH")


def _init_settings():
    # Create the directory where the configuration file will be.
    config_dir = os.path.abspath(os.path.join(
        os.path.expanduser("~"),
        ".gtconfig"
    ))

    if not os.path.isdir(config_dir):
        os.mkdir(config_dir)

    # Check if the configuration file exists.
    config_file = os.path.join(config_dir, "gtrc.ini")
    if not os.path.isfile(config_file):
        _create_default_config(config_file)

    # Read the configuration file.
    config = configparser.RawConfigParser()
    config.read(config_file)

    # Load settings related to the reference genome.
    _init_reference(config)


_init_settings()

