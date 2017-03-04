
Description of pipeline OTU tables
==================================

OTU calling
-----------

The pipeline produces various OTU tables. All tables are constructed
from the set of raw dereplicated reads. These are a set of reads in
FASTA format that appear in the file ``datasetID.raw_dereplicated.fasta``
and correspond to the set of unique sequences present in the raw data.

.. image:: OTU_calling.png
   :alt: Schematic of OTU calling pipeline.
   :width: 80.0%

   Schematic of OTU calling pipeline.

*de novo* OTU tables
~~~~~~~~~~~~~~~~~~~~

The dereplicated reads are clustered to within the specified sequence
similarity percentage cut-off (``OTU_SIMILARITY``, default: 97) using
usearch. Individual reads are assigned either as the OTU centroid or as
a match by usearch's output. Centroids are labeled as ``OTU_ID.0`` and matches count from
``OTU_ID.1`` onwards. These are separate "oligotypes" within the OTU
``OTU_ID``. These are written out to an oligotype table in 
``datasetID.otu_table.\*.denovo_oligotypes``.
Note that this definition of 'oligotype' is simply the unique sequences within
each OTU cluster.

All oligotype counts are then collapsed to their respective parent
OTU, resulting in a fully *de novo* OTU table which can be found in the
filename ``datasetID.otu_table.\*.denovo where`` ``\*`` gets replaced by the
OTU similarity cut-off.

Distribution-based OTU tables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A counts table is made from the dereplicated reads and dereplication
map. The counts table and dereplicated reads are used as inputs to
``dbotu.py``'s ``call_otus()`` function. Briefly, sequences are
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
using a summary file parameter ``RDP_CUTOFF``. This OTU table can be found
in the results sub-folder called ``RDP``. The full classifications can be
found in ``/home/ubuntu/proc/``, in ``RDP_classifications.denovo.txt``
or ``RDP_classifications.dbOTU.txt``.

The dereplicated reads are also aligned to a standard database
(GreenGenes in the case of 16S sequences and UNITE in the case of ITS
sequences). In the case of GreenGenes, the database is determined based
on the specified OTU similarity cut-off: e.g. ``97_otus.fasta`` for
``OTU_SIMILARITY`` set to 97). The alignment is performed using usearch,
and considers the top 10 hits. Consensus assignments are then produced
for the top 1, top 3, top 5 and top 10 hits (where a taxonomic level is
only assigned a latin name if the top N hits from GreenGenes agree), and
the corresponding OTU tables are output. Thus, the OTU table called
``datasetID.otu_table.\*.gg.consensus10`` contains latin names which are
formed from a minimum consensus of the top 10 hits for each taxonomic
level, where ’\*’ gets replaced by the OTU similarity cut-off. Levels are
left unidentified (e.g. ``s__``) if the consensus requirement is not
met.

Note that the database referencing process is one of the slower steps in
the pipeline, so if you only care about RDP, you can skip the GreenGenes
assignments steps by setting the parameter ``GG_ALIGN`` to False in the
summary file. Similarly, for ITS sequences, you can skip the UNITE
assignments steps by setting the parameter ``UNITE_ALIGN`` to False.

Open-reference OTU table
~~~~~~~~~~~~~~~~~~~~~~~~

**DEPRECATED: The pipeline no longer produces open-reference OTUs.**

If any reads in the dereplicated sequences do not align to
GreenGenes/UNITE to within the specified similarity cut-off, they are
clustered with the desired similarity cut-off using usearch, and
appended to the GreenGenes/UNITE closed-reference table to produce a
full, open-reference OTU table (which combines GreenGenes-referenced and
*de novo* OTUs) at the desired similarity cut-off, with filename
datasetID.otu\_table.gg.\*.open\_ref.

Intermediate Files
------------------

A description of the potentially useful intermediate files that you
can find in ``~/proc/``. TODO.