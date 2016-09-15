pfastGO
=======================

pfastGO is an easy to use bioinformatics pipeline that annotates bacterial proteins with GO terms. The pipeline supports multithreading and is able to annotated a bacterial genome with 3000 CDS, using 16 threads in less then 10 minutes.<br />

In the pfastGO pipeline the curated database swissprot is searched using blastp for homolgs of the input protein sequences. The best homologs are selected using the machine learning algorithm DBSCAN. The GO terms of the selected homologs are transferred to the proteins.<br /> 
Additionally, pfamscan is used to find conserved domains in the protein sequences, allowing the annotation of proteins of which no annotated homolog exists. The GO terms of the found conserved domains are also added to the analysis.<br />
Finally, the GO terms are mapped to EC numbers, KEGG ids and Metacyc reactions allowing to set the annotated proteins into the context of pathways.<br />

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

wget https://github.com/danielwuethrich87/pfastGO/archive/master.zip 
unzip master.zip
cd pfastGO-master
sh install.sh<br />

#Usage:

sh phastGo.sh protein_file sample_id cores<br />
 
protein_file: multi fasta file that contains amino acid sequences of proteins (path)<br />
sample_id: short identifier for a sample (string)<br />
cores: number of parallel threads to run (int)<br />

