pfastGO
=======================

pfastGO is an easy to use bioinformatics pipeline that determines the GO term of bacterial proteins. The pipeline supports multithreading and is able to annotated a bacterial genome with 3000 CDS, using 16 cpus within 10 minutes.<br />

It uses blastp to find homologs of already annotated proteins in the manually curated database swissprot. Homologs are selected using the machine learning algorithm DBSCAN, from which the GO terms are transferred.<br /> 
To expand the analysis to proteins, which have no closely related annotated homologs, it also performs pfamscan to find conserved domains in the proteins sequences. The GO terms of the found conserved domains are also added to the analysis.<br />
Finally, the GO terms are mapped to EC numbers, KEGG ids and Metacyc reactions to set the annotated proteins into the context of pathways.<br />

As input a multi fasta file (.faa file) with the amino acid sequences  is needed.<br />

#Requirements:

-Linux 64 bit system<br />

-python (version 2.7)<br />
-python package:NetworkX (available https://networkx.github.io/)<br />
-python package:scikit-learn (available http://scikit-learn.org/stable/)<br />

Following software will be installed automatically with the script: install.sh<br />

-GNU parallel<br />
-blast+<br />
-hmmsearch (version 3.1b2)<br />
-pfam_scan<br />

#Installation: 
sh install.sh<br />

#Usage: 
sh phastGo.sh protein_file sample_id cores<br />
 
protein_file: multi fasta file that contains amino acid sequences of proteins (path)<br />
sample_id: short identifier for a sample (string)<br />
cores: number cpus to run (int)<br />

