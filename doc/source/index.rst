.. amplicon_sequencing_pipeline_doc documentation master file, created by
   sphinx-quickstart on Fri Mar  3 22:29:34 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Amplicon sequencing pipeline
============================

This document describes the process of going from raw 16S or ITS data to
processed data (OTU tables, oligotypes, etc.) using the scripts in the
Alm lab's `processing pipeline`_. 
Most of the 16S and ITS processing steps are orchestrated by
the script ``raw2otu.py``. However, the larger pipeline platform is designed in such 
a manner that the user interfaces with a single script, ``Master.py``, which takes 
as input the path to a folder containing the dataset and a machine-readable file
called a summary file. The summary file tells ``Master.py`` what type of
processing to do (e.g. 16S or ITS) and what steps to do for each type of processing.
Currently, the pipeline can only process 16S and ITS data. 
The format of these summary files, and all files required for 16S or ITS processing, 
are described in :ref:`preparing_dataset`.

.. _processing pipeline: https://github.com/thomasgurry/amplicon_sequencing_pipeline

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
   
