
Database module (``db``)
==========================

Ensembl
--------

`Ensembl <http://www.ensembl.org/>`_ provides a very useful 
`REST API <http://rest.ensembl.org/>`_. We added an interface that does some
throttling and manages the JSON response from the API. 

.. automodule:: genometools.db.ensembl
    :members:

Appris
-------

The `Appris <http://appris.bioinfo.cnio.es/appris.html>`_ database aims to
annotate alternative splice isoforms. We provide an interface to this database
that can annotate :py:class:`genometools.structures.genes.Transcript` objects.

.. automodule:: genometools.db.appris
    :members:

Index
------

Index is a lower level implementation of an indexing data structure. It was
designed to be able to handle any text delimited file with a chromosome column
and a position column.

Indexing is fast because not all the file is read: jumps of a fixed size are 
used. The jump size is estimated by the ``index_rate`` parameter which is
an estimated coverage of the indexed file. Lower rates will take less time to
index and create smaller index files, but will result in slower lookups.

The structure of the index is a pickled python dictionary using chromosomes
as keys and lists of ``(position, file seek)`` as values.

.. automodule:: genometools.db.index
    :members:

