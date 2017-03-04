
Running the pipeline
====================

Processing your data
--------------------

Once you have placed all relevant files and folders, including the
summary file, into a single folder, you can call the script Master.py
with the path to this folder as input. Master.py will parse the summary
file and launch the appropriate scripts to process your request. Suppose
the path to the folder containing the dataset is
``/home/ubuntu/dataset_folder``, and contains a summary file where you
specified the dataset ID as ’myDataset’, you would run the following
command from anywhere:

.. code:: bash

    python ~/scripts/Master.py -i /home/ubuntu/dataset_folder

If you prefer to run the program in the background, you can add an ``&`` to
the end of the command. Alternatively, you can use your favorite
terminal multi-plexer (e.g. ``screen`` or ``tmux``) to run the code while
continuing your work. Using ``nohup`` ensures that any hangups do not
interrupt the processing:

.. code:: bash

    nohup python ~/scripts/Master.py -i /home/ubuntu/dataset_folder 

Processing happens in a folder within the proc folder, with path
``/home/ubuntu/proc/myDataset_proc_16S`` or ``/home/ubuntu/proc/myDataset_proc_ITS``,
depending on the amplicon being analyzed. Final results are put in
``/home/ubuntu/processing_results/``. Log files are put in the
``/home/ubuntu/logs folder``, with two files created for each dataset:
``stderr_datasetID.log`` for error and warning messages and
``stdout_datasetID.log`` for various progress messages. Certain std out and
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
