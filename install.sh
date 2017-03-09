mkdir -p Software/databases/blast

wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.31/ncbi-blast-2.2.31+-x64-linux.tar.gz -P Software
tar xvzf Software/ncbi-blast-2.2.31+-x64-linux.tar.gz -C Software
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_bacteria.dat.gz -P Software/databases/blast
echo "Reformating blast database ..."
python script/convert_swissprot_to_fasta.py Software/databases/blast/uniprot_sprot_bacteria.dat.gz  > Software/databases/blast/uniprot_sprot_bacteria.fasta
Software/ncbi-blast-2.2.31+/bin/makeblastdb -parse_seqids -dbtype prot -in Software/databases/blast/uniprot_sprot_bacteria.fasta

mkdir -p Software/databases/pfam

wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz -P Software
tar xvzf Software/hmmer-3.1b2-linux-intel-x86_64.tar.gz -C Software
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/PfamScan.tar.gz -P Software
tar xvzf Software/PfamScan.tar.gz -C Software
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz -P Software/databases/pfam
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz -P Software/databases/pfam
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/active_site.dat.gz -P Software/databases/pfam

gunzip Software/databases/pfam/* 
Software/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmpress Software/databases/pfam/Pfam-A.hmm

mkdir -p Software/databases/mapping

wget http://geneontology.org/external2go/pfam2go -P Software/databases/mapping
wget http://www.geneontology.org/ontology/go.obo -P Software/databases/mapping

wget http://ftp.gnu.org/gnu/parallel/parallel-20160822.tar.bz2 -P Software
tar -jxvf Software/parallel-20160822.tar.bz2 -C Software
