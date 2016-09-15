#!/bin/bash

export working_dir=$PWD
export cores=$3
export sample_id=$2
export protein_file=$1
export temp="$working_dir"/temp/"$sample_id"
export output_path="$working_dir"/results/"$sample_id"
export pfastGO_location=$(dirname $0)

export PATH=$PATH:"$pfastGO_location"/Software/ncbi-blast-2.2.31+/bin/
export PATH=$PATH:"$pfastGO_location"/Software/parallel-20160822/src/
export PATH=$PATH:"$pfastGO_location"/Software/PfamScan/
export PATH=$PATH:"$pfastGO_location"/Software/hmmer-3.1b2-linux-intel-x86_64/src/

is_command_installed () {
if which $1 &>/dev/null; then
    echo "$1 is installed in" $(which $1)
else
    echo "ERROR: $1 not found"
fi
}

echo
echo "Checking software ..."
echo

is_command_installed parallel
is_command_installed blastp
is_command_installed pfam_scan.pl
is_command_installed python
is_command_installed hmmsearch

if [ which parallel &>/dev/null ]&&[ which blastp &>/dev/null ]&&[ which pfam_scan.pl &>/dev/null ]&&[ which python &>/dev/null ]&&[ which hmmsearch &>/dev/null ];

then

if [ -r "$protein_file" ]&&[ -n "$sample_id" ]&&[ "$cores" -eq "$cores" ];

then

mkdir -p $temp

echo $(date)" ::: Splitting input into "$cores" parallel jobs"
$(
python << END

counter=0
input_file = [n for n in open("$protein_file",'r').read().replace("\r","").split("\n") if len(n)>0]
for line in input_file:
	if line[0:1]==">":
		counter+=1
		if counter>int("$cores"):
			counter=1
		output_file= open("$temp"+"/"+str(counter)+'_proteins.faa', 'a')
    		output_file.write(line+"\n")

	else:
		output_file.write(line+"\n")

END
2>&1)

echo $(date)" ::: Running Blast search"
parallel -j "$cores" 'blastp -db $pfastGO_location/Software/databases/blast/uniprot_sprot_bacteria.fasta -num_threads 1 -query {1} -out {1}.blast.tab -outfmt "6 qseqid sseqid stitle qlen slen length pident nident mismatch gaps evalue bitscore" -max_target_seqs 250' ::: "$temp"/*.faa

echo $(date)" ::: Running pfam domain search"
parallel -j "$cores" 'pfam_scan.pl -dir $pfastGO_location/Software/databases/pfam/ -cpu 1 -as -fasta {1} > {1}.pfam.tsv' ::: "$temp"/*.faa

echo $(date)" ::: Annotating sequences"
parallel -j "$cores" 'python $pfastGO_location/script/parse_blast_sepperate.py {1} {1}.blast.tab {1}.pfam.tsv > {1}.pfastgo_result.tab' ::: "$temp"/*.faa

mkdir -p "$output_path"

echo $(date)" ::: Concatenating results"
grep -v -h "^#" "$temp"/*.faa.pfam.tsv | grep -v "^$" > "$output_path"/"$sample_id".pfam.tsv
cat "$temp"/*.faa.blast.tab > "$output_path"/"$sample_id".blast.tab
cat "$temp"/*.faa.pfastgo_result.tab > "$output_path"/"$sample_id".pfastgo_result.tab

echo $(date)" ::: pfastGO has finished the annotation"

else

echo " "
echo "ERROR: Incorrect input!"
echo "pfastGO version 0.1 by Daniel Wüthrich (danielwue@hotmail.com)"
echo " "
echo "Usage: "
echo "  sh phastGo.sh <protein_file> <sample_id> <cores>"
echo " "
echo "  <protein_file>  multi fasta file that contains amino acid sequences of proteins (path)"
echo "  <sample_id>     short identifier for a sample (string)"
echo "  <cores>         Number cpus to run (int)"
echo " "

fi

else

echo
echo "ERROR: Not all required softwre was found"
echo

fi

