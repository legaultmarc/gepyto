[![Build Status](https://travis-ci.org/legaultmarc/gepyto.svg?branch=master)](https://travis-ci.org/legaultmarc/gepyto)
[![Documentation Status](https://readthedocs.org/projects/gepyto/badge/?version=latest)](https://readthedocs.org/projects/gepyto/?badge=latest)


# Introduction

The name comes from **ge**nome **py**thon **to**ols. Originally, this package
was called `genometools`, but was renamed to avoid confusion as another popular
project uses that name.

The fundamental goal of ``gepyto`` is to provide bioinformaticians with a set
of tools to make developing scripts faster. This means that most features are
included in ``gepyto`` because they provide a clear and eloquent way of
expressing common genomics tasks through Python code. As an example, fetching
gene or variant information from public databases is abstracted away through
extensible Python objects letting users focus on data manipulation rather than
losing time writing code to handle HTTP requests and database queries. 

# Demos

This project will use IPython notebooks for feature demonstrations. They are
available both in the `demos` directory or through nbviewer:

- [Sequence analysis](http://nbviewer.ipython.org/github/legaultmarc/gepyto/blob/master/demos/Sequence%20analysis.ipynb)

We didn't write demos yet, but we have the nice ``gepyto.formats`` module which
adds functionality to parse both ``Impute2`` and ``SeqXML`` files.

# Installation
## For users

First, download gepyto using either the "Download ZIP" button or by using the
``git clone https://github.com/legaultmarc/gepyto.git`` command. After
downloading (and extracting if needed), navigate to the directory and run the
command:

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

<a name="devs_install">
## For developers

See also: [Instructions for contributors](CONTRIBUTING.markdown).

You could use the _"For users"_ instructions, but if you will be changing
things it will be easier to simply add the package to your python path. To
do this, fork the repository, clone it and add it to your `PYTHONPATH`.

1. ``git clone git@github.com:your_user/gepyto.git``
2. Find the absolute path to the _gepyto_ root directory and add it to your
   path by using ``export PYTHONPATH="${PYTHONPATH}:your_absolute_path"``.

You can either use the _export_ command in every terminal you use to work on
_gepyto_ or you can add the line from step 2. to your `~/.bash_profile` or
equivalent.

It is also still a good idea to run ``python setup.py test`` to make sure the
tests pass.

A script automating the previously described steps would look like this:

```shell
git clone git@github.com:your_user/gepyto.git
gtpath=$(find $PWD -maxdepth 1 -name gepyto)
echo 'export PYTHONPATH="${PYTHONPATH}:'${gtpath}'"' >> ~/.bash_profile
```

Don't forget to change `your_user` with your Github username.

Note that this will work on most Linux and Mac OS versions.

# About

This project was intitated by a bioinformatician at the 
[Beaulieu-Saucier Pharmacogenomics Center](http://www.pharmacogenomics.ca/) and 
a student of the [StatGen](http://statgen.org/) lab.

Both are located at the Montreal Heart Institute in Canada.

![StatGen and PGX Logo](https://raw.github.com/legaultmarc/gepyto/master/docs/_static/logo_statgen_pgx.png)

