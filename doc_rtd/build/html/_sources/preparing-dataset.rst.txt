.. _preparing_dataset:

Preparing your dataset for processing
=====================================

Raw data files
--------------

The first step to process your data is to understand what raw data you have. 

* Do you have FASTQ or FASTA data? 
* Are your sequences already de-multiplexed with one file per sample, or will you need to split sequences by barcodes? 
* Do you have unmerged forward and reverse reads that you’ll need to merge? 
* Are the primers still in the sequences?
* Are all reads already all trimmed to a certain length? 

Every sequencing center provides a different kind of “raw data”, so make sure you know 
what you’re starting with! The pipeline takes many different kinds of inputs
and can perform many different processing steps, so it’s important to
know what you’ll be needing.

Summary file
------------

This file, named summary_file.txt, is a machine-readable,
**tab-delimited** file that must accompany any dataset directory when
uploaded to the cloud. It orchestrates all processing that will happen,
and is where you specify each of your processing requests. It should be
found in the highest directory level for the dataset directory in the S3
bucket. It is a text file with descriptors for the data and paths to all
relevant datafiles within the directory. It can include a True or False
flag for whether any associated raw 16S/ITS data has already been
processed. It is case-sensitive.

The order in which items are listed between the lines ``#16S_start`` and
``#16S_end`` (for 16S) and between lines ``#ITS_start`` and ``#ITS_end`` (for
ITS) does not matter.

Note that any white space in the summary file should correspond to a
single tab character.

Summary file format and options
-------------------------------

The first line in your summary file should be:

::

    DATASET_ID  myDataset

16S and ITS processing parameters are specified by attributes placed on
separate lines between ``#16S_start`` and ``#16S_end`` and/or
``#ITS_start`` and ``#ITS_end``. All white spaces in the summary file
should be tabs.

Required attributes are:
 
* ``PRIMERS_FILE``,
* ``BARCODES_MAP``,
* ``PROCESSED``,
*  and one of:
   
   * ``RAW_FASTQ_FILE``,
   * ``RAW_FASTQ_FILES``,
   * ``RAW_FASTA_FILE``, or
   * ``RAW_FASTA_FILES``.

For processing to occur, ``PROCESSED`` should be ``False`` (case-sensitive).

The following section will go through all available summary file
attributes, and are summarized in `List of 16S and ITS attributes`_. These
options are presented in roughly the same order as they are processed.

.. _input-files:

Input file(s)
~~~~~~~~~~~~~

The pipeline takes as input many kinds of raw data:

**one FASTQ file.** This file should contain all sequences for all of
your samples. Sequences may still contain barcodes, or they can have
had barcodes removed and FASTQ headers replaced with sample
identifiers. Specify with ``RAW_FASTQ_FILE``.

**multiple FASTQ files.** Each of these files should contain sequences
for one sample only. Specify with ``RAW_FASTQ_FILES``.

**one FASTA file.** This file should contain all sequences for all of
your samples, and should have a sample identifier in the FASTA header
line. Specify with ``RAW_FASTA_FILE``.

**multiple FASTA files.** Each of these files should contain sequences
for one sample only. Specify with ``RAW_FASTA_FILES``.

If you have one sequence file, ``RAW_FASTQ_FILE`` or ``RAW_FASTA_FILE``
should refer to the FASTQ/A file name. 

If you have multiple sequence files, ``RAW_FASTQ_FILES`` or ``RAW_FASTA_FILES`` 
should refer to a text file that relates each FASTQ/A file to its sample identifier. 
This text file is tab-delimited with the FASTQ/A file name in the first column and
the corresponding sample ID in the second column. This text file should
not contain a header.

An example ``fastq2sid.txt`` file mapping each FASTQ file to its sample is:

::

   SRR1324.fastq    sample1
   SRR1325.fastq    sample2
   SRR1326.fastq    sample3
   ...

Note that all file paths provided should be *relative to the summary
file directory*. If your FASTQ to sample ID file map is in a
subdirectory called ``filemaps/`` and your sequence files are in a
subdirectory called ``datafiles/``, then ``RAW_FASTQ_FILES`` would be
``filemaps/fastq_filemap.txt``, and the entries in ``fastq_filemap.txt`` 
would include the ``datafiles/`` prefix (e.g. ``datafiles/SRR1324.fastq``, etc.)

Merging
~~~~~~~

If you have unmerged paired-end reads, you can merge them by setting
``MERGE_PAIRS`` to ``True`` (case sensitive). Merging can be performed
in both the case where your reads are not de-multiplexed (i.e. you have
one file with *all* of your forward reads for all of your samples, and
one file with all of the reverse reads) and when they are (i.e. you have
two files per sample: one with the forward reads and one with the
reverse reads).

If you need to merge reads, you also need to specify the file name
suffixes for both the forward and reverse read files. These are
specified in ``FWD_SUFFIX`` and ``REV_SUFFIX`` and default to
``_R1.fastq`` and ``_R2.fastq``, respectively. If you have more
complicated file names, either rename them to have consistent suffixes,
or talk to Claire to see if you can incorporate more complicated regex
matching into the code.

Note that the files specified either in ``RAW_FASTQ_FILE`` or in the
``fastq_filemap.txt`` (specified in ``RAW_FASTQ_FILES``) should refer to
the full name of the FASTQ file(s) containing the forward reads.

See `Case 4: multiple demultiplexed raw paired-end FASTQ files of 16S sequences which need merging`_ for an example.

.. _demultiplexing:

De-multiplexing
~~~~~~~~~~~~~~~

If your FASTQ file still contains the barcodes in the sequences, you
will need to include ``BARCODES_MAP`` and ``BARCODES_MODE`` in your
summary file. **\*\*Note that BARCODES_MAP is a required
attribute.** If you do not need to de-multiplex your sequences,
``BARCODES_MAP`` should be ``None`` (case-sensitive).

``BARCODES_MAP`` refers to a tab-delimited file which has the sample
identifiers in the first column and the corresponding barcode sequence
in the second column. This file does not have a header.

An example ``BARCODES_MAP`` file could be::

   C01     AGAGACAT
   C03     AGAGATGT
   C05     AGATGTAG
   C07     AGCGATCT
   C09     AGCTCTAG
   C11     AGTACGAG

You may also specify a ``BARCODES_MODE``, which specifies where the
barcodes are to found in the FASTQ file. If the barcodes are still in
the sequences, ``BARCODES_MODE`` should be 2. If the barcodes are in the
FASTQ sequence header, ``BARCODES_MODE`` should be 1. ``BARCODES_MODE``
defaults to 2.

See `Case 1: raw FASTQ file of 16S sequences, still includes primers and barcodes`_ for an example.

Sometimes, the ’raw’ data has already had primers and barcodes removed
but still has all samples in the same FASTQ file. In this case,
``BARCODES_MAP`` should be ``None`` and the sample IDs must be listed in
the sequence header lines of the FASTQ file. If there is text other than
the sample ID in the header, you need to specify the first non-sample ID
character in ``BARCODES_SEPARATOR``. For example, sequences in these
kinds of files are often labeled like:

::

    @sample1_seq1
    ...<rest of fastq record>
    @sample1_seq2
    ...<rest of fastq record>
    @sample2@_seq1
    ...<rest of fastq record>
    @sample3_seq1
    ...<rest of fastq record>
    @sample2_seq2
    ...<rest of fastq record>

In this case, the barcodes separator would be an underscore (``_``),
which is the default.

Primer trimming
~~~~~~~~~~~~~~~

If you need to remove primers from your sequences, you can specify
``PRIMERS_FILE``, a text file with your primer sequences. **\*\*Note
that PRIMERS_FILE is a required attribute.** If you do not need to
remove primers from your sequences, ``PRIMERS_FILE`` should be ``None``
(case sensitive).

Your primers file should have each primer on its own line and no header::

     CCTACGGGAGGCAGCAG
     ATTACCGCGGCTGCT

The pipeline does not currently remove reverse primers. If your
sequences still contain reverse primers, you can remove them yourself or
trim your sequences to a length shorter than the start of your reverse
primer.

Quality filtering
~~~~~~~~~~~~~~~~~

There are two ways to quality filter your sequences. One is based on the
number of expected errors in your sequence, and the other truncates
reads after a certain quality is encountered. You can learn more about
these approaches by reading the USEARCH documentation:
http://www.drive5.com/usearch/manual/readqualfiltering.html

**To truncate reads after a base with a certain quality is encountered**,
use the ``QUALITY_TRIM`` option. A default value that is often used is
25. This step is performed *before* length trimming.

**To discard reads based on their number of expected errors**, use the
``MAX_ERRORS`` option. A default value that is often used is 2 (i.e.
reads with more than 2 expected errors are discarded). This step is
performed *after* length trimming

If nothing is specified, the pipeline defaults to ``QUALITY_TRIM`` of
25. If both ``MAX_ERRORS`` and ``QUALITY_TRIM`` are specified, quality
filtering by truncation is performed (i.e. ``MAX_ERRORS`` is ignored).

You may also need to specify the encoding of the quality scores.
``ASCII_ENCODING`` can be either ``ASCII_BASE_33`` (default) or
``ASCII_BASE_64``. You can check the encoding of your file using
usearch: ``usearch -fastq_chars yourFASTQfile.fastq``

Length trimming
~~~~~~~~~~~~~~~

By default, the pipeline trims all reads to 101 base pairs before
dereplication and clustering. You can specify a different length by
using ``TRIM_LENGTH``. Any reads which are shorter than the specified
length are discarded.

Dereplication
~~~~~~~~~~~~~

In the dereplication step, unique sequences are identified and the
samples from which they came are tracked (sometimes referred to as
“provenancing”). By default, unique sequences which are present fewer
than 10 times in the entire dataset are discarded. If you want to change
this number, specify it with ``MIN_COUNT``. (e.g. if ``MIN_COUNT`` is 2,
only singleton sequences are discarded).

OTU calling
~~~~~~~~~~~

You can specify the similarity used to define OTUs in the
``OTU_SIMILARITY`` attribute. The default value is 97, corresponding to
97% OTUs.

By default, the pipeline clusters OTUs using both *de novo* and
closed-reference approaches. If you specify an OTU similarity that does
not have a corresponding Green Genes reference file, closed-reference
clustering will not be performed. OTU similarities supported by Green
Genes closed-reference mapping are: 61, 64, 67, 70, 73, 76, 79, 82, 85,
88, 91, 94, 97, and 99%. The database files used for this mapping can be
found in ``/home/ubuntu/databases/gg_13_5_otus/rep_set_latin/``.

The pipeline assigns taxonomies to *de novo* OTUs using the naive-Bayes
RDP classifier. By default, the confidence cutoff is 0.5. You can
specify a different value with the ``RDP_CUTOFF`` attribute.

Distribution-based OTU calling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The pipeline also performs distribution-based OTU calling [1]_. You
can set the abundance, distance and p value criteria in the summary
file attributes ``DISTANCE_CRITERIA``, ``ABUNDANCE_CRITERIA``, and ``DBOTU_PVAL``.

Distribution-based clustering is **not** performed by default. You can
turn it on by setting the summary file attribute ``DBOTU`` to ``True``.

Sample summary files
--------------------

Case 1: raw FASTQ file of 16S sequences, still includes primers and barcodes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The simplest case is if you have the following files: a raw FASTQ file;
a file specifying the map between barcode sequences and IDs; and a file
specifying the primers used. Your summary file would look something like
this:

::

    DATASET_ID  myDataset

    #16S_start
    RAW_FASTQ_FILE      myData.fastq
    ASCII_ENCODING      ASCII_BASE_33
    PRIMERS_FILE        primers.txt
    BARCODES_MAP        barcodes_map.txt
    BARCODES_MODE       2
    METADATA_FILE       metadata.txt
    PROCESSED           False
    #16S_end

Note that you **must** also specify the place where barcodes are to be
found, i.e. either in the "``>``” sequence ID lines (mode 1) or in
the sequences themselves (mode 2). The ``PROCESSED`` flag tells the
processing instance that the dataset needs to be processed into OTU
tables.

Your ``barcodes_map.txt`` file would look something like this:

::

    S1      ATCGCTAGTA
    S2      TCGCTATATA
    S3      TCTACAGCGT
    S4      CGTACTCAGT

And your ``primers.txt`` file could be::

    CCTACGGGAGGCAGCAG
    ATTACCGCGGCTGCT

Case 2: raw FASTQ file of ITS sequences, primers and barcodes have been removed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the case where the ’raw’ data has already had primers and barcodes
removed (but is not yet de-multiplexed, i.e. all samples are still in
the same FASTQ file), the sample IDs must be listed in the sequence ID
lines of the FASTQ file. When the pipeline removes barcodes itself and
replaces them with sample IDs, individual sequence reads for a given
sampleID will be annotated as sampleID;1, sampleID;2, etc., where we
note here that the BARCODES\_SEPARATOR is ’;’. However, in a dataset
where the barcodes have previously been removed, you will have to look
into the FASTQ file to check the ’separator’ character. Your summary
file would look something like this:

.. code:: bash

    DATASET_ID  myDataset

    #ITS_start
    RAW_FASTQ_FILE     myData.fastq
    ASCII_ENCODING     ASCII_BASE_33
    PRIMERS_FILE       None
    BARCODES_MAP       None
    BARCODES_SEPARATOR ;
    METADATA_FILE      metadata.txt
    PROCESSED          False
    #ITS_end

Case 3: multiple demultiplexed raw FASTQ or FASTA files of 16S sequences, each file corresponding to a single sample
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes sequencing data are available in a demultiplexed form, where
the reads for each sample are split into separate files. Many datasets
in the SRA, for example, are available in this form. In this case, you
can create a **two-column, tab-delimited file** where the first column
lists the filename and the second column lists the corresponding sample
ID. Note that paths should be relative paths *within the current
directory*, e.g. ``datafiles/file1.txt`` for files in a folder called ``datafiles``
within the current directory. In the summary file, the ``RAW_FASTQ_FILE``
line becomes ``RAW_FASTQ_FILES`` (plural), and instead refers to this
filename. If your files are FASTA rather than FASTQ, simply use
``RAW_FASTA_FILES`` (also plural). For a filename fastq_filemap.txt, your
summary file would look something like this:

::

    DATASET_ID  myDataset

    #16S_start
    RAW_FASTQ_FILES     fastq_filemap.txt
    ASCII_ENCODING      ASCII_BASE_33
    PRIMERS_FILE        primers.txt
    METADATA_FILE       metadata.txt
    PROCESSED           False
    PRIMERS_FILE        None
    BARCODES_MAP        None
    #16S_end

And your fastq_filemap.txt file would look something like this (note
that white spaces in the following example correspond to a single tab
character):

::

    SRR10001.fastq  S1
    SRR10002.fastq  S2
    SRR10003.fastq  S3
    SRR10004.fastq  S4


Case 4: multiple demultiplexed raw paired-end FASTQ files of 16S sequences which need merging
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If your 16S FASTQ files are split into forward and reverse paired-end
reads, the pipeline can merge them for you. Specify ``MERGE_PAIRS`` in the
summary file, and also include the filename suffixes corresponding to
forward and reverse reads. If your forward read fastq files were named
``sampleID_L001_R1.fastq`` and your reverse read fastq files were named
``sampleID_L001_R2.fastq``, your summary file would look something like
this:

::

    DATASET_ID  myDataset

    #16S_start
    RAW_FASTQ_FILES fastq_filemap.txt
    PRIMERS_FILE    None
    BARCODES_MAP    None
    MERGE_PAIRS     True
    FWD_SUFFIX      _L001_R1.fastq
    REV_SUFFIX      _L001_R2.fastq
    PROCESSED       False
    #16S_end

And your fastq_filemap would look like

::

    S1_L001_R1.fastq  S1
    S2_L001_R1.fastq  S2
    S3_L001_R1.fastq  S3
    S4_L001_R1.fastq  S4


If you have de-multiplexed files (as in this example), the file names in
your fastq_filemap.txt file should be the forward read fastq files.

If instead you have non-demultiplexed sequences (i.e. two fastq files,
one containing your forward reads and one containing your reverse
reads), ``RAW_FASTQ_FILE`` should point to the file containing the forward
reads.

.. _attribute-list:

List of 16S and ITS attributes
------------------------------

+-----------------------+-------------------------------------------------------------------------------+
| **Attribute**         | **Description**                                                               |
+=======================+===============================================================================+
| RAW\_FASTQ\_FILE      | Raw FASTQ file name/path within the dataset directory. If ``MERGE_PAIRS`` is  |
|                       | ``True``, this should be the full name of the forward read file.              |
+-----------------------+-------------------------------------------------------------------------------+
| RAW\_FASTA\_FILE      | Raw FASTA file name/path if raw data is in FASTA format                       |
+-----------------------+-------------------------------------------------------------------------------+
| RAW\_FASTQ\_FILES     | For demultiplexed datasets where samples are separated                        |
|                       | into separate FASTQ files. Filename of two column file                        |
|                       | containing FASTQ filenames in first column and sample IDs                     |
|                       | in the second column. If ``MERGE_PAIRS`` is ``True``, these file              |
|                       | names should be the full forward read file names.                             |
+-----------------------+-------------------------------------------------------------------------------+
| RAW\_FASTA\_FILES     | For demultiplexed datasets where samples are separated                        |
|                       | into separate FASTA files. Filename of two column file                        |
|                       | containing FASTA filenames in first column and sample IDs                     |
|                       | in the second column.                                                         |
+-----------------------+-------------------------------------------------------------------------------+
| ASCII\_ENCODING       | ASCII quality encoding in FASTQ. Supports either                              |
|                       | ``ASCII_BASE_33`` or ``ASCII_BASE_64``. Set to 33 if unspecified.             |
+-----------------------+-------------------------------------------------------------------------------+
| PRIMERS\_FILE         | | Filename/path to primers file.                                              |
|                       | | **required**: If primers have already been removed, specify ``None``.       |
+-----------------------+-------------------------------------------------------------------------------+
| BARCODES\_MAP         | Filename/path to barcodes map file.                                           |
|                       | Tab-delimited file contains sampleIDs in first column                         |
|                       | and barcode sequences in second column.                                       |
|                       | | **required**: If barcodes have already been removed, specify ``None``.      |
+-----------------------+-------------------------------------------------------------------------------+
| BARCODES\_MODE        | | ``1`` = barcodes in sequence ID,                                            |
|                       | | ``2`` = barcodes in sequences themselves.                                   |
|                       | | ``3`` = barcodes in separate index file (beta)                              |
|                       | | Required if ``BARCODES_MAP`` is not None.                                   |
+-----------------------+-------------------------------------------------------------------------------+
| BARCODES\_SEPARATOR   | Separator character. See description in                                       |
|                       | `De-multiplexing`_                                                            |
+-----------------------+-------------------------------------------------------------------------------+
| METADATA\_FILE        | Filename/path to metadata file.                                               |
+-----------------------+-------------------------------------------------------------------------------+
| MERGE\_PAIRS          | If need to merge paired-end reads, set to ``True``.                           |
+-----------------------+-------------------------------------------------------------------------------+
| FWD\_SUFFIX           | Filename suffix of files with forward reads.                                  |
|                       | Should include filename extension. If not specified, defaults to ``_1.fastq`` |
+-----------------------+-------------------------------------------------------------------------------+
| REV\_SUFFIX           | Filename suffix of files with reverse reads.                                  |
|                       | Should include filename extension. If not specified, defaults to ``_2.fastq`` |
+-----------------------+-------------------------------------------------------------------------------+
| PROCESSED             | True/False flag for whether data have already been processed.      `          |
|                       | required: Set to ``False`` for processing to proceed.                         |
+-----------------------+-------------------------------------------------------------------------------+
| TRIM\_LENGTH          | Length to which all sequences should be trimmed.                              |
|                       | Defaults to 101 if unspecified.                                               |
+-----------------------+-------------------------------------------------------------------------------+
| QUALITY\_TRIM         | Minimum quality score allowed.                                                |
|                       | Sequences are truncated at the first base having quality                      |
|                       | score less than value. Defaults to 25 if unspecified.                         |
|                       | If set to None, no quality filtering or trimming will be performed.           |
|                       | If both ``QUALITY_TRIM`` and ``MAX_ERRORS`` are included in summary file,     |
|                       | ``MAX_ERRORS`` will be ignored (even if ``QUALITY_TRIM`` = ``None``).         |
+-----------------------+-------------------------------------------------------------------------------+
| MAX\_ERRORS           | Maximum expected errors allowed.                                              |
|                       | After length trimming, sequences with more than ``MAX_ERRORS``                |
|                       | expected errors are discarded. If not specified or if a ``QUALITY_TRIM`` value|
|                       | is specified, defaults to quality trimming behavior, above.                   |
+-----------------------+-------------------------------------------------------------------------------+
| MIN\_COUNT            | Minimum sequence count in dereplication across                                |
|                       | all samples. Defaults to 10 if unspecified (i.e. sequences with fewer than    |
|                       | 10 occurrences in the entire dataset will not be considered downstream).      |
+-----------------------+-------------------------------------------------------------------------------+
| OTU\_SIMILARITY       | Integer specifying the percent similarity desired                             |
|                       | in OTU clustering. Defaults to 97 if unspecified.                             |
+-----------------------+-------------------------------------------------------------------------------+
| RDP\_CUTOFF           | Desired probability cut-off for Ribosomal Database Project                    |
|                       | assignments. Assignments at each taxonomic level will be evaluated and        |
|                       | those with a lower probability than this cutoff will be                       |
|                       | labeled as unidentified. Defaults to 0.5 if unspecified.                      |
+-----------------------+-------------------------------------------------------------------------------+
| GG\_ALIGN             | Specific to 16S sequences. True/False flag for whether                        |
|                       | GreenGenes alignments are desired. Defaults to ``True`` if unspecified.       |
+-----------------------+-------------------------------------------------------------------------------+
| UNITE\_ALIGN          | Specific to ITS sequences. True/False flag for whether                        |
|                       | UNITE alignments are desired. Defaults to ``True`` if unspecified.            |
+-----------------------+-------------------------------------------------------------------------------+
| DBOTU                 | Whether to perform distribution-based OTU calling.                            |
|                       | Defaults to ``False`` if unspecified.                                         |
+-----------------------+-------------------------------------------------------------------------------+
| ABUNDANCE\_CRITERIA   | Abundance criteria for distribution-based OTU calling.                        |
|                       | Defaults to 10 if unspecified.                                                |
+-----------------------+-------------------------------------------------------------------------------+
| DISTANCE\_CRITERIA    | Distance criteria for distribution-based OTU calling.                         |
|                       | Defaults to 0.1 if unspecified.                                               |
+-----------------------+-------------------------------------------------------------------------------+
| DBOTU\_PVAL           | P value cutoff for distribution-based OTU calling.                            |
|                       | Defaults to 0.0005 if unspecified.                                            |
+-----------------------+-------------------------------------------------------------------------------+
| OUTDIR                | Full path to processing directory. Defaults                                   |
|                       | to /home/ubuntu/proc/ if not specified.                                       |
+-----------------------+-------------------------------------------------------------------------------+



.. [1]
   http://almlab.mit.edu/dbotu3.html
