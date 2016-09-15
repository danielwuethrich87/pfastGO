pfastGO
=======================

pfastGO is an easy to use bioinformatics pipeline that determines the GO term of bacterial proteins. The pipeline supports multithreading and is able to annotated a bacterial genome with 3000 CDS, using 16 cpus within 10 minutes.

It uses blastp to find homologs of already annotated proteins in the manually curated database swissprot. Homologs are selected using the machine learning algorithm DBSCAN, from which the GO terms are transferred. 
To expand the analysis to proteins, which have no closely related annotated homologs, it also performs pfamscan to find conserved domains in the proteins sequences. The GO terms of the found conserved domains are also added to the analysis.
Finally, the GO terms are mapped to EC numbers, KEGG ids and Metacyc reactions to set the annotated proteins into the context of pathways.

As input a multi fasta file (.faa file) with the amino acid sequences  is needed.

#Requirements:

Linux 64 bit system

-python (version 2.7)
-python packages:
NetworkX (available https://networkx.github.io/)
scikit-learn (available http://scikit-learn.org/stable/)

-GNU parallel*
-blast+*
-hmmsearch (version 3.1b2)*
-pfam_scan*

*software will be installed with the script: install.sh

#Installation: 
sh install.sh

#Usage: 
sh phastGo.sh <protein_file> <sample_id> <cores>
 
  <protein_file>	multi fasta file that contains amino acid sequences of proteins (path)
  <sample_id>		short identifier for a sample (string)
  <cores>		Number cpus to run (int)

