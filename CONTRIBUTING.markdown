# Instructions for developers

Thank you for your interest in ``genometools``. To get started, you can follow
the instructions from the [README](README.markdown#devs_install). After you
will be ready to start working on genometools!

## Fork, clone, commit and all that

To get more information on how to collaborate using git/Github, you can visit 
Github's very good 
[documentation](https://help.github.com/categories/collaborating/) on the 
subject.

## Quality

To make sure to have a high quality package, we will only merge code that is
thoroughly tested and documented. Documentation is achieved by using the
``sphinx`` docstrings (See Coding guidelines section). After that, if you
created a new module, add a _RestructuredText_ file in the ``docs`` folder by
following the pattern of what is already in there. Basically, generating the
documentation is as easy as:

```
.. automodule:: genometools.a_package.a_module
    :members:
```

Writing good tests is also a prerequisite for merging. All the tests from the
``genometools.tests`` directory will be automatically included in the build
process. They can be invoked using the ``setup.py test`` command.

We use the dafault Python ``unittest`` module. Reading the official
documentation is probably the best way to understand how to write tests. You
can also look at previously implemented _TestCase_s.

Also note that the tests should pass on both Python 2.7 and Python 3.4.

## Coding guidelines

We mostly follow the [Google Python Style Guide](https://google-styleguide.googlecode.com/svn/trunk/pyguide.html).

The most important things to note are:

* 80 characters maximum line length
* Use 4 spaces for indentation
* Two blank lines between top-level definitions and one blank line between 
  method definitions.
* Use [format](https://docs.python.org/2/library/functions.html#format) to
  manipulate strings.
* Naming: CamelCase for classes, CAPS\_FOR\_GLOBALS, underscore\_notation for 
  everything else.

Also for dependencies, we try to use the standard library if convenient. If a
third-party library is needed, we usually just import it. The user will need
to install it manually if he uses the corresponding `genometools` functionality.

We also document all public functions for use with Sphinx autodoc. Reading
a little bit of our current code should allow you to familiarize with the system
easily, but the general structure looks like this:

```python
def a_cool_function(a):
    """A short description for the function.

    :param a: Description of the first argument.
    :type a: int

    :returns: What this function returns.
    :rtype: float

    A longer description, notes or details.

    """

    return float(a)

```

You will need to add your module to the documentation in the `docs` directory.
Refer to the [Sphinx documentation](http://sphinx-doc.org/) and the existing
`genometools` documentation to get more information.

## Design

The idea behind this package was to provide bioinformaticists working in Python
a general purpose framework to work with "omics" data. It originated from our
realisation that we kept writing the same functions accross projects and that
doing it well, once was probably better. Hence, we try to make genometools
easy to use, small and modular so we can reuse it often accross projects.

If you add nice things to the package, it would be a good idea to write 
thorough tests, documentation and [IPython Notebook](http://ipython.org/notebook.html) 
demonstrations in the ``demos`` directory.

