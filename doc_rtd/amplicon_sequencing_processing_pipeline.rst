
Running the pipeline
====================

Processing your data
--------------------

Once you have placed all relevant files and folders, including the
summary file, into a single folder, you can call the script Master.py
with the path to this folder as input. Master.py will parse the summary
file and launch the appropriate scripts to process your request. Suppose
the path to the folder containing the dataset is
/home/ubuntu/dataset\_folder, and contains a summary file where you
specified the dataset ID as ’myDataset’, you would run the following
command from anywhere:

::

    python ~/scripts/Master.py -i /home/ubuntu/dataset_folder

If you prefer to run the program in the background, you can add an & to
the end of the command. Alternatively, you can use your favorite
terminal multi-plexer (e.g. screen or tmux) to run the code while
continuing your work. Using nohup ensures that any hangups do not
interrupt the processing:

::

    nohup python ~/scripts/Master.py -i /home/ubuntu/dataset_folder &

Processing happens in a folder within the proc folder, with path
/home/ubuntu/proc/
myDataset\_proc\_16S or /home/ubuntu/proc/myDataset\_proc\_ITS,
depending on the amplicon being analyzed. Final results are put in
/home/ubuntu/processing\_results/. Log files are put in the
/home/ubuntu/logs folder, with two files created for each dataset:
stderr\_datasetID.log for error and warning messages and
stdout\_datasetID.log for various progress messages. Certain std out and
std err messages also go to the ``stderr_master.log`` and
``stdout_master.log`` in the directory from which you called the
``python Master.py`` command. Make sure to check both of these places
when troubleshooting!

**Please remove your files from the ~/proc folder once you have checked
your processing results!** The ~/proc folder can easily fill up - it
contains the original raw files are and additional copies corresponding
to each processing step!

Underlying function calls
-------------------------

The following table describes the underlying scripts and functions used
in each step of 16S and ITS data processing. Python scripts are
passed-down Alm lab scripts (for the most part).

All code is in ``/home/ubuntu/scripts`` and we encourage users to take a
look!

+---------------------------------------+------------------------------------------------------------------+
| **Processing**                        | **Function/script call**                                         |
+=======================================+==================================================================+
| Merging                               | ``usearch8 -fastq_mergepairs``                                   |
+---------------------------------------+------------------------------------------------------------------+
| De-multiplexing                       | ``2.split_by_barcodes.py``                                       |
+---------------------------------------+------------------------------------------------------------------+
| Primer trimming                       | ``1.remove_primers.py``                                          |
+---------------------------------------+------------------------------------------------------------------+
| Quality filtering (truncation)        | ``usearch8 -fastq_filter -fastq_truncqual``                      |
+---------------------------------------+------------------------------------------------------------------+
| Quality filtering (expected errors)   | ``usearch8 -fastq_filter -fastq_maxee``                          |
+---------------------------------------+------------------------------------------------------------------+
| Length trimming (FASTQ)               | ``usearch8 -fastq_filter -fastq_trunclen``                       |
+---------------------------------------+------------------------------------------------------------------+
| Length trimming (FASTA)               | ``usearch8 -fastx_truncate -trunclen``                           |
+---------------------------------------+------------------------------------------------------------------+
| Dereplication                         | ``3.dereplicate.py``                                             |
+---------------------------------------+------------------------------------------------------------------+
| De novo clustering                    | ``usearch8 -cluster_otus -otu_radius_pct``                       |
+---------------------------------------+------------------------------------------------------------------+
| Closed-reference mapping              | ``usearch8 -usearch_global -db GG_database_file -strand both``   |
+---------------------------------------+------------------------------------------------------------------+
| RDP taxonomy assignment               | ``rdp_classify.py``                                              |
+---------------------------------------+------------------------------------------------------------------+
| Distribution-based clustering         | ``dbotu.py`` ``call_otus()``                                     |
+---------------------------------------+------------------------------------------------------------------+

Pipeline output files and directories
=====================================

The pipeline outputs different OTU tables and corresponding
representative sequences. All final outputs can be found in the
processing\_results folder, under a sub-directory labeled
myDataset\_results. Files in this directory are labeled systematically,
usually with the format
myDataset.file\_description.otu\_similarity.file\_type.

Top directory
-------------

Files in the top results directory are as follows:

-  myDataset.otu\_seqs.N.fasta: FASTA file with the representative
   sequences for the denovo OTUs, clustered at N%.

-  myDataset.otu\_seqs.dbOTU.fasta: FASTA file with the representative
   sequences for the distribution-based OTUs.

-  myDataset.otu\_table.N.denovo: OTU table with N% denovo OTUs labeled
   denovo1, denovo2, ... in the rows and samples in the columns.

-  myDataset.otu\_table.N.dbOTU: OTU table with distribution-based OTUs
   labeled dbotu1, dbotu2, ... in the rows and samples in the columns.

-  myDataset.otu\_table.N.denovo\_oligotypes: OTU table with N% denovo
   OTUs separated into unique oligotypes. Each OTU is labeled denovo1.1,
   denovo1.2, denovo2.1, denovo2.2, denovo2.3, etc. The first number
   corresponds to the parent denovo OTU number; the second is the
   oligotype number. Oligotypes are calculated as each unique sequence
   within an OTU cluster.

-  myDataset.raw\_dereplicated.fasta: FASTA file with all unique
   sequences in the dataset. Only sequences which appear more times than
   the MIN\_COUNT specified in the summary file are included (default
   value for MIN\_COUNT is 10).

-  summary\_file.txt: Updated summary file containing original
   processing request and resulting file names.

Also within each dataset’s processing\_results directory, 3
sub-directories are created:

Quality control
---------------

Within the results folder, there is a subfolder called quality\_control,
which contains various plots diagnostic of dataset quality. Currently,
the pipeline outputs:

-  Histogram showing distribution of read lengths, taken from the first
   100,000 reads in the raw FASTQ file.

-  Bar chart showing number of reads per sample.

-  File showing percentage of reads thrown out at each processing step
   (processing\_summary.txt).

Note that the information in these files is not always accurate - you
should probably do more thorough quality control yourself.

RDP
---

This folder contains the RDP-assigned OTU tables. It has two files:

-  myDataset.otu\_table.N.denovo.rdp\_assigned: OTU table with denovo
   OTUs assigned Latin names with RDP. OTUs are in rows and samples are
   in columns. OTU names are of the format:

   ::

        k__kingdom;p__phylum;c__class;o__order;f__family;g__genus;s__species;
        d__denovoID

   where the denovoID corresponds to the respective sequence in
   ../myDataset.otu\_seqs.N.fasta.

-  myDataset.otu\_table.dbOTU.rdp\_assigned: OTU table with
   distribution-based OTUs assigned Latin names with RDP. OTUs are in
   rows and samples are in columns. OTU names are of the format:

   ::

        k__kingdom;p__phylum;c__class;o__order;f__family;g__genus;s__species;
        d__dbotuID

   where the dbotuID corresponds to the respective sequence in
   ../myDataset.otu\_seqs.dbOTU.fasta.

GG
--

This folder contains the Green Genes-assigned OTU table. It has multiple
files, which are described further in the following section.

-  myDataset.otu\_table.N.gg.consensusM: closed-reference OTU table with
   dereplicated sequences mapped to Green Genes using usearch. OTUs are
   in rows and samples are in columns. OTU names are of the format:

   ::

        k__kingdom;p__phylum;c__class;o__order;f__family;g__genus;s__species;
        d__derepID

   where the derepID corresponds to the sequence in
   ../myDataset.raw\_dereplicated.fasta with the same ID number.

Description of pipeline OTU tables
==================================

OTU calling
-----------

The pipeline produces various OTU tables. All tables are constructed
from the set of raw dereplicated reads. These are a set of reads in
FASTA format that appear in the file datasetID.raw\_dereplicated.fasta
and correspond to the set of unique sequences present in the raw data.

.. figure:: figs/OTU_calling.png
   :alt: Schematic of OTU calling pipeline.
   :width: 80.0%

   Schematic of OTU calling pipeline.

*de novo* OTU tables
~~~~~~~~~~~~~~~~~~~~

The dereplicated reads are clustered to within the specified sequence
similarity percentage cut-off (OTU\_SIMILARITY, default: 97) using
usearch. Individual reads are assigned either as the OTU centroid or as
a match. Centroids are labeled as ’OTU\_ID.0’ and matches count from
’OTU\_ID.1’ onwards. These are separate oligotypes within the OTU
’OTU\_ID’. All oligotype counts are then collapsed to their respective
OTU, resulting in a fully *de novo* OTU table which can be found in the
filename datasetID.otu\_table.\*.denovo where ’’ gets replaced by the
OTU similarity cut-off.

Distribution-based OTU tables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A counts table is made from the dereplicated reads and dereplication
map. The counts table and dereplicated reads are used as inputs to
``dbotu.py``\ ’s ``call_otus()`` function. Briefly, sequences are
compared against existing OTU representative sequences in descending
abundance order. If the candidate sequence’s abundance is within the
``ABUNDANCE_CRITERIA`` of an existing OTU and its count distribution
across samples is not significantly different from the OTU’s, then the
candidate sequence is grouped with that OTU. You can learn more about
the distribution-based clustering algorithm at
``http://almlab.mit.edu/dbotu3.html``.

Intermediate files such as the membership file (which lists which
sequences belong to which OTU) and the table of sequence counts (used as
an input to the ``dbotu.py`` function call) are in the respective
``/home/ubuntu/proc/`` folder for your dataset.

Closed-reference OTU table
~~~~~~~~~~~~~~~~~~~~~~~~~~

There are two types of database-referencing that are outputed from the
pipeline by default. The first is assignments from the Ribosomal
Database Project (RDP), which returns probabilistic assignments. The
default probability cut-off below which an assignment is labeled as
unidentified at a given taxonomic level is 0.5, but this can be set
using a summary file parameter RDP\_CUTOFF. This OTU table can be found
in the results sub-folder called ’RDP’. The full classifications can be
found in ``/home/ubuntu/proc/``, in ``RDP_classifications.denovo.txt``
or ``RDP_classifications.dbOTU.txt``.

The dereplicated reads are also aligned to a standard database
(GreenGenes in the case of 16S sequences and UNITE in the case of ITS
sequences). In the case of GreenGenes, the database is determined based
on the specified OTU similarity cut-off: e.g. 97\_otus.fasta for
OTU\_SIMILARITY set to 97). The alignment is performed using usearch,
and considers the top 10 hits. Consensus assignments are then produced
for the top 1, top 3, top 5 and top 10 hits (where a taxonomic level is
only assigned a latin name if the top N hits from GreenGenes agree), and
the corresponding OTU tables are output. Thus, the OTU table called
datasetID.otu\_table.\*.gg.consensus10 contains latin names which are
formed from a minimum consensus of the top 10 hits for each taxonomic
level, where ’’ gets replaced by the OTU similarity cut-off. Levels are
left unidentified (e.g. ’s\_\_’) if the consensus requirement is not
met.

Note that the database referencing process is one of the slower steps in
the pipeline, so if you only care about RDP, you can skip the GreenGenes
assignments steps by setting the parameter GG\_ALIGN to False in the
summary file. Similarly, for ITS sequences, you can skip the UNITE
assignments steps by setting the parameter UNITE\_ALIGN to False.

Open-reference OTU table
~~~~~~~~~~~~~~~~~~~~~~~~

DEPRECATED: The pipeline no longer produces open-reference OTUs.

If any reads in the dereplicated sequences do not align to
GreenGenes/UNITE to within the specified similarity cut-off, they are
clustered with the desired similarity cut-off using usearch, and
appended to the GreenGenes/UNITE closed-reference table to produce a
full, open-reference OTU table (which combines GreenGenes-referenced and
*de novo* OTUs) at the desired similarity cut-off, with filename
datasetID.otu\_table.gg.\*.open\_ref.

Oligotypes
~~~~~~~~~~

The OTU calling process automatically assigns separate IDs to unique
sequences within OTUs (termed ’oligotypes’ by some). These offer a
greater level of granularity. For example, if sequences A, B and C are
all the same OTU, called OTU1 (and the centroid is sequence A), these
would be relabeled ’OTU1.0’, ’OTU1.1’ and ’OTU1.2’. Note that these
unique ‘oligotypes’ are generated by any fancy information-based
algorithms - they are simply the unique sequences within each OTU.

| There is currently one such table output: the full *de novo* oligotype
  table (where sequences are labeled according to a fully *de novo* OTU
  clustering process) -
| myDataset.otu\_table.N.denovo\_oligotypes

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

-  preprocessing\_16S.py - Methods and wrappers for raw 16S sequence
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

-  | rdp\_classify.py - Wrapper for calling the RDP classifier jar file,
   | found at /home/ubuntu/tools/RDPTools/classifier.jar.

Troubleshooting
===============

Standard error (stderr) and standard out (stdout) are directed to the
directory /home/ubuntu/logs. Errors will appear in the stderr file for
the corresponding datasets. There may also be useful errors in the
``stderr_master.log`` and ``stdout_master.log`` files that are created
in the directory from which you ran the ``python Master.py`` call. In
addition, file outputs from each step of the pipeline can be found in
/home/ubuntu/proc/datasetID\_proc\_16S, which can help to diagnose
errors. Common sources of errors or difficulties include:

-  Bad formatting of the summary, barcodes, or primers file. There
   should be no white spaces, and columns should be tab-delimited. If
   you created the file in Excel, it may have carriage return or other
   non-linux formatting characters that introduce difficulties. Safest
   is to go from a previous summary file template, or to create it
   manually.

-  Typos in a 16S attribute key-value pair in the summary file.

-  Incorrect ASCII encoding specified. Check whether 33 or 64 using
   usearch8 -fastq\_chars yourFASTQfile.fastq. Default is 33, so if left
   unspecified and the file is ASCII base 64, quality trimming will
   fail.

-  For ITS sequences, the RDP classifier often runs out of RAM when
   loading the UNITE database. This is usually what happened if you find
   RDP\_classification.txt to be empty in the relevant folder in
   /home/ubuntu/proc. If you rerun the pipeline several times, have
   checked everything else and keep encountering this problem, you may
   need a node with more RAM to do your processing.

-  If you specify a file containing a list of FASTQ/FASTA files for the
   raw reads, make sure the file paths are relative to the directory
   hosting the summary file.

Appendix
========

List of 16S and ITS attributes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+-----------------------+-------------------------------------------------------------------------------+
| **Attribute**         | **Description**                                                               |
+=======================+===============================================================================+
| RAW\_FASTQ\_FILE      | Raw FASTQ file name/path within the dataset directory                         |
+-----------------------+-------------------------------------------------------------------------------+
| RAW\_FASTA\_FILE      | Raw FASTA file name/path if raw data is in FASTA format                       |
+-----------------------+-------------------------------------------------------------------------------+
| RAW\_FASTQ\_FILES     | For demultiplexed datasets where samples are separated                        |
+-----------------------+-------------------------------------------------------------------------------+
|                       | into separate FASTQ files. Filename of two column file                        |
+-----------------------+-------------------------------------------------------------------------------+
|                       | containing FASTQ filenames in first column and sample IDs                     |
+-----------------------+-------------------------------------------------------------------------------+
|                       | in the second column.                                                         |
+-----------------------+-------------------------------------------------------------------------------+
| RAW\_FASTA\_FILES     | For demultiplexed datasets where samples are separated                        |
+-----------------------+-------------------------------------------------------------------------------+
|                       | into separate FASTA files. Filename of two column file                        |
+-----------------------+-------------------------------------------------------------------------------+
|                       | containing FASTA filenames in first column and sample IDs                     |
+-----------------------+-------------------------------------------------------------------------------+
|                       | in the second column.                                                         |
+-----------------------+-------------------------------------------------------------------------------+
| ASCII\_ENCODING       | ASCII quality encoding in FASTQ. Supports either                              |
+-----------------------+-------------------------------------------------------------------------------+
|                       | ’ASCII\_BASE\_33’ or ’ASCII\_BASE\_64’. Set to 33 if unspecified.             |
+-----------------------+-------------------------------------------------------------------------------+
| PRIMERS\_FILE         | Filename/path to primers file.                                                |
+-----------------------+-------------------------------------------------------------------------------+
|                       | required: If primers have already been removed, specify ’None’.               |
+-----------------------+-------------------------------------------------------------------------------+
| BARCODES\_MAP         | Filename/path to barcodes map file.                                           |
+-----------------------+-------------------------------------------------------------------------------+
|                       | Tab-delimited file contains sampleIDs in first column                         |
+-----------------------+-------------------------------------------------------------------------------+
|                       | and barcode sequences in second column.                                       |
+-----------------------+-------------------------------------------------------------------------------+
|                       | required: If barcodes have already been removed, specify ’None’.              |
+-----------------------+-------------------------------------------------------------------------------+
| BARCODES\_MODE        | ’1’ = barcodes in sequence ID,                                                |
+-----------------------+-------------------------------------------------------------------------------+
|                       | ’2’ = barcodes in sequences themselves.                                       |
+-----------------------+-------------------------------------------------------------------------------+
|                       | Required if BARCODES\_MAP is not None.                                        |
+-----------------------+-------------------------------------------------------------------------------+
| BARCODES\_SEPARATOR   | Separator character. See description in Case 2 below.                         |
+-----------------------+-------------------------------------------------------------------------------+
| METADATA\_FILE        | Filename/path to metadata file.                                               |
+-----------------------+-------------------------------------------------------------------------------+
| MERGE\_PAIRS          | If need to merge paired-end reads, set to “True”.                             |
+-----------------------+-------------------------------------------------------------------------------+
| FWD\_SUFFIX           | Filename suffix of files with forward reads.                                  |
+-----------------------+-------------------------------------------------------------------------------+
|                       | Should include filename extension. If not specified, defaults to \_1.fastq    |
+-----------------------+-------------------------------------------------------------------------------+
| REV\_SUFFIX           | Filename suffix of files with reverse reads.                                  |
+-----------------------+-------------------------------------------------------------------------------+
|                       | Should include filename extension. If not specified, defaults to \_2.fastq    |
+-----------------------+-------------------------------------------------------------------------------+
| PROCESSED             | True/False flag for whether data have already been processed.                 |
+-----------------------+-------------------------------------------------------------------------------+
|                       | required: Set to ’False’ for processing to proceed.                           |
+-----------------------+-------------------------------------------------------------------------------+
| TRIM\_LENGTH          | Length to which all sequences should be trimmed.                              |
+-----------------------+-------------------------------------------------------------------------------+
|                       | Defaults to 101 if unspecified.                                               |
+-----------------------+-------------------------------------------------------------------------------+
| QUALITY\_TRIM         | Minimum quality score allowed.                                                |
+-----------------------+-------------------------------------------------------------------------------+
|                       | Sequences are truncated at the first base having quality                      |
+-----------------------+-------------------------------------------------------------------------------+
|                       | score less than value. Defaults to 25 if unspecified.                         |
+-----------------------+-------------------------------------------------------------------------------+
|                       | If set to None, no quality filtering or trimming will be performed.           |
+-----------------------+-------------------------------------------------------------------------------+
|                       | If both QUALITY\_TRIM and MAX\_ERRORS are included in summary file,           |
+-----------------------+-------------------------------------------------------------------------------+
|                       | MAX\_ERRORS will be ignored (even if QUALITY\_TRIM = None).                   |
+-----------------------+-------------------------------------------------------------------------------+
| MAX\_ERRORS           | Maximum expected errors allowed.                                              |
+-----------------------+-------------------------------------------------------------------------------+
|                       | After length trimming, sequences with more than MAX\_ERRORS                   |
+-----------------------+-------------------------------------------------------------------------------+
|                       | expected errors are discarded. If not specified or if a TRIM\_QUALITY value   |
+-----------------------+-------------------------------------------------------------------------------+
|                       | is specified, defaults to quality trimming behavior, above.                   |
+-----------------------+-------------------------------------------------------------------------------+
| MIN\_COUNT            | Minimum sequence count in dereplication across                                |
+-----------------------+-------------------------------------------------------------------------------+
|                       | all samples. Defaults to 10 if unspecified (i.e. sequences with fewer than    |
+-----------------------+-------------------------------------------------------------------------------+
|                       | 10 occurrences in the entire dataset will not be considered downstream).      |
+-----------------------+-------------------------------------------------------------------------------+
| OTU\_SIMILARITY       | Integer specifying the percent similarity desired                             |
+-----------------------+-------------------------------------------------------------------------------+
|                       | in OTU clustering. Defaults to 97 if unspecified.                             |
+-----------------------+-------------------------------------------------------------------------------+
| RDP\_CUTOFF           | Desired probability cut-off for Ribosomal Database Project                    |
+-----------------------+-------------------------------------------------------------------------------+
|                       | assignments. Assignments at each taxonomic level will be evaluated and        |
+-----------------------+-------------------------------------------------------------------------------+
|                       | those with a lower probability than this cutoff will be                       |
+-----------------------+-------------------------------------------------------------------------------+
|                       | labeled as unidentified. Defaults to 0.5 if unspecified.                      |
+-----------------------+-------------------------------------------------------------------------------+
| GG\_ALIGN             | Specific to 16S sequences. True/False flag for whether                        |
+-----------------------+-------------------------------------------------------------------------------+
|                       | GreenGenes alignments are desired. Defaults to ’True’ if unspecified.         |
+-----------------------+-------------------------------------------------------------------------------+
| UNITE\_ALIGN          | Specific to ITS sequences. True/False flag for whether                        |
+-----------------------+-------------------------------------------------------------------------------+
|                       | UNITE alignments are desired. Defaults to ’True’ if unspecified.              |
+-----------------------+-------------------------------------------------------------------------------+
| DBOTU                 | Whether to perform distribution-based OTU calling.                            |
+-----------------------+-------------------------------------------------------------------------------+
|                       | Defaults to False if unspecified.                                             |
+-----------------------+-------------------------------------------------------------------------------+
| ABUNDANCE\_CRITERIA   | Abundance criteria for distribution-based OTU calling.                        |
+-----------------------+-------------------------------------------------------------------------------+
|                       | Defaults to 10 if unspecified.                                                |
+-----------------------+-------------------------------------------------------------------------------+
| DISTANCE\_CRITERIA    | Distance criteria for distribution-based OTU calling.                         |
+-----------------------+-------------------------------------------------------------------------------+
|                       | Defaults to 0.1 if unspecified.                                               |
+-----------------------+-------------------------------------------------------------------------------+
| DBOTU\_PVAL           | P value cutoff for distribution-based OTU calling.                            |
+-----------------------+-------------------------------------------------------------------------------+
|                       | Defaults to 0.0005 if unspecified.                                            |
+-----------------------+-------------------------------------------------------------------------------+
| OUTDIR                | Full path to processing directory. Defaults                                   |
+-----------------------+-------------------------------------------------------------------------------+
|                       | to /home/ubuntu/proc/ if not specified.                                       |
+-----------------------+-------------------------------------------------------------------------------+

.. [1]
   http://almlab.mit.edu/dbotu3.html
