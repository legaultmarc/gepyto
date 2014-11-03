
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

