Pipeline output files and directories
=====================================

The pipeline outputs different OTU tables and corresponding
representative sequences. All final outputs can be found in the
``processing_results`` folder, under a sub-directory labeled
``myDataset_results``. Files in this directory are labeled systematically,
usually with the format
``myDataset.file_description.otu_similarity.file_type``.

Top directory
-------------

Files in the top results directory are as follows:

-  **myDataset.otu_seqs.N.fasta**: FASTA file with the representative
   sequences for the denovo OTUs, clustered at N%.

-  **myDataset.otu_seqs.dbOTU.fasta**: FASTA file with the representative
   sequences for the distribution-based OTUs.

-  **myDataset.otu_table.N.denovo**: OTU table with N% denovo OTUs labeled
   denovo1, denovo2, ... in the rows and samples in the columns.

-  **myDataset.otu_table.N.dbOTU**: OTU table with distribution-based OTUs
   labeled dbotu1, dbotu2, ... in the rows and samples in the columns.

-  **myDataset.otu_table.N.denovo_oligotypes**: OTU table with N% denovo
   OTUs separated into unique oligotypes. Each OTU is labeled denovo1.1,
   denovo1.2, denovo2.1, denovo2.2, denovo2.3, etc. The first number
   corresponds to the parent denovo OTU number; the second is the
   oligotype number. Oligotypes are calculated as each unique sequence
   within an OTU cluster.

-  **myDataset.raw_dereplicated.fasta**: FASTA file with all unique
   sequences in the dataset. Only sequences which appear more times than
   the ``MIN_COUNT`` specified in the summary file are included (default
   value for ``MIN_COUNT`` is 10).

-  **summary_file.txt**: Updated summary file containing original
   processing request and resulting file names.

Also within each dataset’s ``processing_results`` directory, 3
sub-directories are created:

Quality control
---------------

Within the results folder, there is a subfolder called ``quality_control``,
which contains various plots diagnostic of dataset quality. Currently,
the pipeline outputs:

-  Histogram showing distribution of read lengths, taken from the first
   100,000 reads in the raw FASTQ file.

-  Bar chart showing number of reads per sample.

-  File showing percentage of reads thrown out at each processing step
   (``processing_summary.txt``).

Note that the information in these files is not always accurate - you
should probably do more thorough quality control yourself.

RDP
---

This folder contains the RDP-assigned OTU tables. It has two files:

-  **myDataset.otu_table.N.denovo.rdp_assigned**: OTU table with denovo
   OTUs assigned Latin names with RDP. OTUs are in rows and samples are
   in columns. OTU names are of the format:

   ::

        k__kingdom;p__phylum;c__class;o__order;f__family;g__genus;s__species;
        d__denovoID

   where the denovoID corresponds to the respective sequence in
   ``../myDataset.otu_seqs.N.fasta``.

-  **myDataset.otu_table.dbOTU.rdp_assigned**: OTU table with
   distribution-based OTUs assigned Latin names with RDP. OTUs are in
   rows and samples are in columns. OTU names are of the format:

   ::

        k__kingdom;p__phylum;c__class;o__order;f__family;g__genus;s__species;
        d__dbotuID

   where the dbotuID corresponds to the respective sequence in
   ``../myDataset.otu_seqs.dbOTU.fasta``.

GG
--

This folder contains the Green Genes-assigned OTU table. It has multiple
files, which are described further in the following section.

-  **myDataset.otu_table.N.gg.consensusM**: closed-reference OTU table with
   dereplicated sequences mapped to Green Genes using usearch. OTUs are
   in rows and samples are in columns. OTU names are of the format:

   ::

        k__kingdom;p__phylum;c__class;o__order;f__family;g__genus;s__species;
        d__derepID--GGggid

   where the derepID corresponds to the sequence in
   ``../myDataset.raw_dereplicated.fasta`` with the same ID number.
