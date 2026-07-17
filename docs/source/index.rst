CoRE MOF Tools
==============

CoRE MOF Tools is the Python interface for querying and processing the CoRE MOF
database and for running the associated curation, validation, descriptor, and
property-prediction workflows.

Start with :doc:`installation`, then run ``coremof doctor`` to see which features
are available in your environment. Most workflows accept a CIF path and either
return a dictionary/dataframe or write documented files to an output directory.

.. warning::

   CSD downloads and MOSAEC require a licensed CSD installation. Zeo++ and MOFid
   are external programs. Heat-capacity prediction requires the ensemble files
   from the full GitHub repository. A successful base installation does not imply
   that these optional components are present.

.. toctree::
   :maxdepth: 2
   :caption: User guide

   installation
   quickstart
   features
   troubleshooting
   references

.. toctree::
   :maxdepth: 2
   :caption: API reference

   CoREMOF

Project links
-------------

* `GitHub repository <https://github.com/Chung-Research-Group/CoRE-MOF-Tools>`_
* `CoRE MOF web application <https://mof-db.pusan.ac.kr/>`_
* `Issue tracker <https://github.com/Chung-Research-Group/CoRE-MOF-Tools/issues>`_

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
