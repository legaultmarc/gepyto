#!/usr/bin/env python

# How to build source distribution
# python setup.py sdist --format bztar
# python setup.py sdist --format gztar
# python setup.py sdist --format zip


import os
from setuptools import setup


setup(
    name="gepyto",
    version="0.1",
    description=("Utilities and tools to interface with genomic databases "
                 "and for genome bioinformatics."),
    author="Marc-Andre Legault",
    author_email="legaultmarc@gmail.com",
    url="https://github.com/legaultmarc",
    license="CC BY-NC 4.0",
    packages=["gepyto", "gepyto.annotation", "gepyto.db", "gepyto.structures",
              "gepyto.tests", "gepyto.utils", "gepyto.formats"],
    package_data={"gepyto.db": ["data/*", ], },
    classifiers=["Operating System :: Linux",
                 "Programming Language :: Python",
                 "Programming Language :: Python :: 2.7",
                 "Programming Language :: Python :: 3.4"],
    test_suite="gepyto.tests.test_suite",
    install_requires=["numpy >= 1.8.1", "requests >= 2.4.3", "pandas >= 0.15"],
)

