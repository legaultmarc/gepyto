[![Build Status](https://travis-ci.org/legaultmarc/genometools.svg?branch=master)](https://travis-ci.org/legaultmarc/genometools)

[![Documentation Status](https://readthedocs.org/projects/genometools/badge/?version=latest)](https://readthedocs.org/projects/genometools/?badge=latest)

# Introduction

This is a small module that contains reusable code. I aim to use it as a 
submodule for other projects.

# Features

- Python objects for Variants including Indels and SNPs.
- Python objects for Genes and Transcripts.
- Abstraction of the APPRIS transcript annotation database.
- Abstraction of the Ensembl API for variant queries in a range and from an rs
  number.
- Throttling for Ensembl database queries.

# Installation
## For users

First, download genometools using either the "Download ZIP" button or by using 
the ``git clone https://github.com/legaultmarc/genometools.git`` command.
After downloading (and extracting if needed), navigate to the directory and
run the command:

```shell
python setup.py install
```

This will install the package into your current Python (so if you want to use
it in a virtual environment, make sure you activate it first).

You can also test the package using:

```shell
python setup.py test
```

Note that the test coverage is fairly low for now. Don't hesitate to contact us
to report problems with the installation.

## For developers

You could use the _"For users"_ instructions, but if you will be changing
things it will be easier to simply add the package to your python path.

1. ``git clone git@github.com:legaultmarc/genometools.git``
2. Find the absolute path to the _genometools_ root directory and add it to your
   path by using ``export PYTHONPATH="${PYTHONPATH}:your_absolute_path"``.

You can either use the _export_ command in every terminal you use to work on
_genometools_ or you can add the line from step 2. to your `~/.bash_profile` or
equivalent.

It is also still a good idea to run ``python setup.py test`` to make sure the
tests pass.

A script automating the previously described steps would look like this:

```shell
git clone git@github.com:legaultmarc/genometools.git
gtpath=$(find $PWD -maxdepth 1 -name genometools)
echo 'export PYTHONPATH="${PYTHONPATH}:'${gtpath}'"' >> ~/.bash_profile
```

Note that this will work on most Linux and Mac OS versions.

# WIP

- Checking of reference alleles for Variants using a provided genome reference.
- Annotation module.
- Get all genes from region in Ensembl DB (utils).

# Demos

This project will use IPython notebooks for feature demonstrations. For now, 
there is one that can be viewed on [NBViewer](http://nbviewer.ipython.org/github/legaultmarc/genometools/blob/master/demos/Variant%20Annotation.ipynb) or from the `demos` 
directory. It demonstrates how this module can be used to fetch information 
about genes and variants and how to annotate the latter.

