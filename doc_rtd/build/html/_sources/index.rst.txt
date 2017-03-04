.. amplicon_sequencing_pipeline_doc documentation master file, created by
   sphinx-quickstart on Fri Mar  3 22:29:34 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Introduction
============

This document describes the process of going from raw 16S or ITS data to
processed data (OTU tables, oligotypes, etc.). This is orchestrated by
the script ``raw2otu.py``. The platform is designed in such a manner that
the user interfaces with a single script, ``Master.py``, which takes as
input the path to the folder containing the dataset. Each dataset folder
must contain a machine-readable text file, called a summary file, which
gives instructions to ``Master.py``. The format of these summary files, and
all files required for 16S or ITS processing, are described in :ref:`preparing_dataset`.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   quickstart
   preparing-dataset
   running
   output
   otu-tables
   source-code
   troubleshooting
   


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
