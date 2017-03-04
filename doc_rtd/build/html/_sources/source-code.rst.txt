
Source code
===========

The following code can all be in ``/home/ubuntu/scripts``. Some of these
scripts are used by the pipeline, others are not currently incorporated
but may still prove useful to your 16S-related work.

Python modules
--------------

-  Formatting.py - Module to house miscellaneous formatting methods,
   e.g. conversion from classic dense format to BIOM format, OTU table
   transposition, etc.

-  QualityControl.py - Methods for quality control diagnostics on a
   dataset.

-  preprocessing_16S.py - Methods and wrappers for raw 16S sequence
   data processing.

-  Taxonomy.py - Methods for taxonomy-related feature extraction and
   analytics. Includes functions for things like: -adding latin names to
   a GreenGenes-referenced OTU table -collapsing abundances at different
   taxonomic levels

-  Phylogeny.py - Methods for phylogenetic feature extraction, e.g.
   left/right (LR) abundance ratios at each node of a phylogenetic tree.

-  Analytics.py - Generic statistical analysis tools, e.g. Wilcoxon
   tests across all available taxa.

-  Regressions.py - Performs different types of regressions *en masse*.

-  PipelineFilesInterface.py - Methods for reading and moving around raw
   data files and for reading specific groups of attributes from the
   summary file.

Scripts and routines
--------------------

-  Master.py - Master script that calls relevant processing pipelines,
   e.g. raw2otu.py.

-  raw2otu.py - Pipeline for converting raw 16S FASTQ sequence files to
   OTU tables. Handles parallelization requirements in these processing
   steps automatically. Takes as input a directory that contains a
   summary file and the raw data.

-   dbotu.py - Module for doing distribution-based OTU calling.

-  rdp_classify.py - Wrapper for calling the RDP classifier jar file,
   found at ``/home/ubuntu/tools/RDPTools/classifier.jar``.
