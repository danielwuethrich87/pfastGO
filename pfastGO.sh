#$ -q all.q
#$ -e $JOB_ID.pfam.err
#$ -o $JOB_ID.pfam.out
#$ -cwd #executes from the current directory and safes the ouputfiles there
#$ -pe smp 16


date
export working_dir=$PWD
export temp=/scratch/"$JOB_ID"."$SGE_TASK_ID".all.q
export PERL5LIB=/home/dwuethrich/Application/pfamscan/PfamScan
export PATH=$PATH:/home/dwuethrich/Application/pfamscan/PfamScan
module add SequenceAnalysis/HMM-Profile/hmmer/3.1b2;
module add Blast/ncbi-blast/2.2.31+;




for i in FAM1476c1p FAM1476c2p FAM10921p FAM13875p FAM14217p FAM18098p FAM18149p FAM18814p FAM19036p FAM19038p FAM19132p FAM20860p FAM21731p FAM22155p FAM22284p FAM22367p FAM8105p FAM8627p

do

mkdir -p "$working_dir"/results/"$i"
mkdir -p "$working_dir"/temp/"$i"

export protein_file=$working_dir/../../annotation/"$i"/"$i"_07142016.faa
export output_path=$working_dir/temp/"$i"/

$(
python << END

counter=0
input_file = [n for n in open("$protein_file",'r').read().replace("\r","").split("\n") if len(n)>0]
for line in input_file:
	if line[0:1]==">":
		counter+=1
		if counter>int("$NSLOTS"):
			counter=1
		output_file= open("$output_path"+str(counter)+'_proteins.faa', 'a')
    		output_file.write(line+"\n")

	else:
		output_file.write(line+"\n")

END
2>&1) 



parallel -j "$NSLOTS" '

blastp -db /data5/users/dwuethrich/bacteria_swiss_prot_19_08_2016/format/blastdbs/uniprot_sprot_bacteria.fasta -num_threads 1 -query {1} -out {1}.blast.tab -outfmt "6 qseqid sseqid qlen slen length pident nident mismatch gaps evalue bitscore" -max_target_seqs 250

pfam_scan.pl -dir ~/Application/pfamscan/PfamScan/pfam_hmm/ -cpu 1 -as -fasta {1} > {1}.pfam.tsv

python script/parse_blast_sepperate.py  {1} {1}.blast.tab {1}.pfam.tsv > {1}.pfastgo_result.tab

' ::: "$working_dir"/temp/"$i"/*.faa

grep -v -h "^#" "$working_dir"/temp/"$i"/*.faa.pfam.tsv | grep -v "^$" > "$working_dir"/results/"$i"/"$i".pfam.tsv
cat "$working_dir"/temp/"$i"/*.faa.blast.tab > "$working_dir"/results/"$i"/"$i".blast.tab
cat "$working_dir"/temp/"$i"/*.faa.pfastgo_result.tab > "$working_dir"/results/"$i"/"$i".pfastgo_result.tab

rm -rf "$working_dir"/temp/

done

date

#pfam_scan.pl -dir ~/Application/pfamscan/PfamScan/pfam_hmm/ -cpu "$NSLOTS" -as -fasta "$working_dir"/../../annotation/"$i"/"$i"_*.fa
#awk -v path="$working_dir/temp/$i/" '/>/ {counter+=1; OUT=path"protein_"counter".fa"}; OUT{print >OUT}' "$working_dir"/../../annotation/"$i"/"$i"_*.faaa > out.tsv
