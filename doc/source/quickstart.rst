==========
Quickstart
==========

Prepare your data
-------------------

To process your data, put all of your files in the same directory
and write your summary file. The summary file contains all of the
instructions for processing that the pipeline will use.

Summary file
~~~~~~~~~~~~

The summary file should be called ``summary_file.txt`` and be 
tab-delimited.

Summary file attributes tell the processing script what to do.
Attributes should be in all-caps in the first column, and
associated values in the second column.
Attributes should be between ``#16S_start`` and ``#16S_end``
or ``#ITS_starts`` and ``#ITS_end`` lines.
Everything in the summary file is case-sensitive. 
Your summary file must have the following attributes:

* ``PRIMERS_FILE``,
* ``BARCODES_MAP``,
* ``PROCESSED``,
*  and one of:
   
   * ``RAW_FASTQ_FILE``,
   * ``RAW_FASTQ_FILES``,
   * ``RAW_FASTA_FILE``, or
   * ``RAW_FASTA_FILES``.

Your summary file could look like this:

::

    DATASET_ID  myDataset

    #16S_start
    RAW_FASTQ_FILES     fastq_filemap.txt
    PROCESSED           False
    PRIMERS_FILE        None
    BARCODES_MAP        None
    #16S_end

The script considers everything relative to the summary file,
so any relative file paths must be valid with respect to where
your summary file is.

All of the available summary file attributes and their associated
default values are in :ref:`attribute-list`.

Data
~~~~

The pipeline can handle FASTQ and FASTA files. It can also handle
un-demultiplexed data (i.e. reads for all samples are in one file) or
de-multiplexed data (i.e. reads for each sample are in separate files).
It can also handle un-merged paired-end data.

If you have de-multiplexed data, you need to specify a file map
with two columns. The first column has the data files and the 
second column has the corresponding sample IDs (see
:ref:`input-files` for more).

If you need to de-multiplex your data (e.g. remove the barcodes),
you should provide a barcodes map with sample IDs in the first
column and barcodes in the second column (see :ref:`demultiplexing`
for more).

In general, none of your auxiliary files should have headers (e.g.
primers, barcodes, and data-to-sample maps).

Run the pipeline
--------------------

From anywhere, you can run:

.. code:: bash

	python ~/scripts/Master.py -i ~/path/to/summary/file/directory

If your data will take a while to process, we recommend using
a screen or tmux session.

Troubleshooting and logs
------------------------

Check both ``stderr_master.log`` and ``stdout_master.log`` in your current directory and the 
respective file in ``~/logs`` (e.g. ``stderr_datasetID_proc_16s.log``).
If you used nohup, also look at nohup.out in your current directory.

The most common problem is messed up white space in one of your files,
for example:

- extra tabs in summary files
- ^M character in your barcodes or fastq file map (has to do with Mac OS's newline character)

The next most common problem is typos:

- check your file paths
- check the summary file attributes and values

Sometimes the processing can mess up because your data is doing
something unexpected. An easy way to troubleshoot these problems 
is to look at file sizes in ``~/proc/yourdataset_proc_16S``.
File suffixes indicate what's been done to each file:

* ``.pt`` - primer trimmed
* ``.qt`` - quality trimmed
* ``.sb`` - split by barcodes 
* ``.lt`` - length trimmed

Getting help
------------

If you're still getting an error, you can email Claire and Thomas
for help with your problem. You can also email the alm-comp [at] mit.edu
email list to see if others have encountered something similar.