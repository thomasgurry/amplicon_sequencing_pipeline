
Troubleshooting
===============

Standard error (stderr) and standard out (stdout) are directed to the
directory ``/home/ubuntu/logs``. Errors will appear in the ``stderr`` file for
the corresponding datasets. There may also be useful errors in the
``stderr_master.log`` and ``stdout_master.log`` files that are created
in the directory from which you ran the ``python Master.py`` call. In
addition, file outputs from each step of the pipeline can be found in
``/home/ubuntu/proc/datasetID_proc_16S``, which can help to diagnose
errors. Common sources of errors or difficulties include:

-  Bad formatting of the summary, barcodes, or primers file. There
   should be no white spaces, and columns should be tab-delimited. If
   you created the file in Excel, it may have carriage return or other
   non-linux formatting characters that introduce difficulties. Safest
   is to go from a previous summary file template, or to create it
   manually.

   Sometimes emacs or other editors will insert multiple tab characters
   where you only pressed tab once. Make sure that all white spaces
   correspond to **one** tab character only.

   Moving text files between OS's can also lead to errors. For example,
   the Mac newline character is a weird ``^M`` carriage return. Make
   sure that any of your files replace these with normal newline (``\n``)
   characters. One way to check this is to ``less`` your file
   (e.g. :code:`$ less <your_file>`).

-  Typos in a 16S attribute key-value pair in the summary file.

-  Incorrect ASCII encoding specified. Check whether 33 or 64 using
   ``usearch8 -fastq_chars yourFASTQfile.fastq``. Default is 33, so if left
   unspecified and the file is ASCII base 64, quality trimming will
   fail.

-  For ITS sequences, the RDP classifier often runs out of RAM when
   loading the UNITE database. This is usually what happened if you find
   ``RDP_classification.txt`` to be empty in the relevant folder in
   ``/home/ubuntu/proc``. If you rerun the pipeline several times, have
   checked everything else and keep encountering this problem, you may
   need a node with more RAM to do your processing.

-  If you specify a file containing a list of FASTQ/FASTA files for the
   raw reads, make sure the file paths are relative to the directory
   hosting the summary file.
