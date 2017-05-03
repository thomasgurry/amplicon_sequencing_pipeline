# 16S Amplicon Sequencing Processing Pipeline

This repo contains scripts to process raw 16S or ITS sequencing data into processed data like dereplicated sequences, OTU tables, etc.
The platform is designed in such a manner that the user interfaces with a single script, `Master.py`, which takes as input the path to 
the folder containing the dataset raw files. Each dataset folder must contain a machine-readable text file, called a summary file, 
which gives instructions to `Master.py`.

You can read more about how to use the pipeline and what it does in the [documentation](http://amplicon-sequencing-pipeline.readthedocs.io/en/latest/).

