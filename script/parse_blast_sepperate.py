import operator
import csv
import sys
import math
from Bio.Blast import NCBIXML
import networkx as nx
import numpy as np
from sklearn.cluster import DBSCAN
import re
inputOptions = sys.argv[1:]

# python parse_blast_sepperate.py ../test/longest_proteins.fasta ../test/blast/results/longest_proteins/longest_proteins_swissprot_blast.csv ../test/pfam_only/out.tsv
p_value=float("1e-6")

def main():

	sys.stderr.write("pfastGO:Reading FASTA\n")
	sequences={}
	sequences=read_fasta(inputOptions[0])
	sys.stderr.write("pfastGO:Reading BLAST output\n")
	sequences=read_blast(sequences,inputOptions[1]) # add first blast (swiss prot)
	sys.stderr.write("pfastGO:Reading PFAM output\n")
	sequences=add_interpro_information(sequences,inputOptions[2])
	sys.stderr.write("pfastGO:Building GO graph\n")
	sequences=add_full_path_to_protein(sequences)

	#print "Locus id\tProkka\tUniprot homolog\tGOs\tmetacyc\tEC\tKegg"
	for seq_id in sorted(sequences.keys()):
		print seq_id+"\t"+sequences[seq_id].prokka_name+"\t"+sequences[seq_id].name+"\t"+sequences[seq_id].most_downstream_GO_string()+"\t"+sequences[seq_id].GO_string()+"\t"+sequences[seq_id].metacyc_string()+"\t"+sequences[seq_id].ec_string()+"\t"+sequences[seq_id].kegg_string()
	sys.stderr.write("pfastGO:Finished\n")


def read_fasta(input_file):

	sequences={}
	input_file = [n for n in open(input_file,'r').read().replace("\r","").split("\n") if len(n)>0]
	for line in input_file:
		if line[0:1]==">":
			seq_id = line[1:].split(" ")[0]
			sequences[seq_id]=(Sequence(seq_id))
			sequences[seq_id].prokka_name=line[1:].replace(seq_id+" ","")
		else:
			sequences[seq_id].sequence+=line


	return sequences


def read_blast(sequences,input_file):


	seq_id2blast_hits={}
	input_file = [n for n in open(input_file,'r').read().replace("\r","").split("\n") if len(n)>0]
	for line in input_file:

		
		seq_id=line.split("\t")[0]

		if (seq_id in seq_id2blast_hits.keys()) == bool(0):
			seq_id2blast_hits[seq_id]=[]
			
		hsp=Hsp()	
		hsp.hit_def=line.split("\t")[1]
		hsp.gos=str(hsp.hit_def).split("|GO=")[1].split("|")[0].split(";")
		hsp.name = hsp.hit_def.split("|name=")[1].split(" |")[0].split("|")[0]
		hsp.align_length=int(line.split("\t")[4])
		hsp.gaps=int(line.split("\t")[8])
		hsp.score=float(line.split("\t")[10])
		hsp.identities=int(line.split("\t")[6])
		hsp.e_value=float(line.split("\t")[9])
		#if cal_HSSP_curve(hsp) >=5: # 5 are 90% of enzymes correct (http://www.sciencedirect.com/science/article/pii/S0022283602000165 (fig. 3))
			#print str(cal_HSSP_curve(hsp))+"\t"+line
		if hsp.e_value<=p_value:
			seq_id2blast_hits[seq_id].append(hsp)

	for seq_id in seq_id2blast_hits.keys():

		blast_hits = sorted(seq_id2blast_hits[seq_id], key=lambda x: x.score, reverse=True)

		sequences=assign_go_terms(blast_hits,seq_id,sequences)


	return sequences


def assign_name(blast_hits,sequence,first_clusters):
	# Assigning the name of the best blast hit of cluster 0 and with a specific name
	for hsp,cluster in zip(blast_hits,first_clusters):
		if sequence.name.lower()=="hypothetical protein" and cluster==first_clusters[0]:# changes from "cluster==0"
				
			if len(hsp.gos)>=1:		
				sequence.name=hsp.name

	return sequence


def assign_go_terms(blast_hits,seq_id,sequences):

	if len(blast_hits)>0:
		sequence_length=len(sequences[seq_id].sequence)
		X=list()

		for hsp in blast_hits:
			X.append([float(hsp.align_length)/float(sequence_length), float(hsp.identities)/float(hsp.align_length)])

		X.append([0,0])	
		X.append([1,1])
		X=np.array(X)
		eps_value=0.2
		db = DBSCAN(eps=eps_value, min_samples=1).fit(X)
		clusters = db.labels_

		sequences[seq_id]=add_go_term_by_cluster(blast_hits,sequences[seq_id],clusters)



	return sequences


def add_go_term_by_cluster(blast_hits,sequence,clusters):
	# The GO terms of the clusters build using DB scan are summed up. The to go term of the best 20 or all blast hits of the best cluster are collected.

	cluster_with_best_hit=clusters[0] # first hit is best blast hit
	GOs_of_best_cluster=[]

	counter=1
	for hsp,cluster in zip(blast_hits,clusters[:len(clusters)-2]):
		if (cluster == cluster_with_best_hit and counter <= 20):
			counter+=1
			GOs_of_best_cluster.extend(hsp.gos)


	sequence.GO_terms=list(set(GOs_of_best_cluster))
	sequence=assign_name(blast_hits,sequence,clusters[:len(clusters)-2])


	return sequence


def add_interpro_information(sequences,input_file):
	seq_id2interpro = read_interpro(input_file,sequences)
	for seq_id in sequences.keys():

		for interpro_result in seq_id2interpro[seq_id]:
			sequences[seq_id].GO_terms=list(set(sequences[seq_id].GO_terms).union(set(interpro_result.GO_terms)))


		sequences=alter_names_of_unknown(seq_id,sequences,seq_id2interpro)

		#if sequences[seq_id].name.lower().find("transposase")!=-1: # add tranposase GO
		#	sequences[seq_id].GO_terms.append("GO:0004803")		
		sequences[seq_id].GO_terms=list(set(sequences[seq_id].GO_terms)) # remove duplicate GO

	return sequences


def alter_names_of_unknown(seq_id,sequences,seq_id2interpro):

	sorted_interpro_results = sorted(seq_id2interpro[seq_id], key=lambda x: x.p_value, reverse=False)
	for interpro_result in sorted_interpro_results:
		if sequences[seq_id].name.lower()=="hypothetical protein":
			if interpro_result.p_value<1 and interpro_result.name.lower().find("protein of unknown function")==-1 and interpro_result.name.lower().find("uncharacterised protein")==-1:
				sequences[seq_id].name= interpro_result.name
	return sequences


def read_interpro(input_file,sequences):

	pfam2go=read_pfam2go()
	seq_id2interpro={}
	for seq_id in sequences.keys():
		seq_id2interpro[seq_id]=list()

	input_file = [n for n in open(input_file,'r').read().replace("\r","").split("\n") if len(n)>0]
	for line in input_file:
		line=re.sub(' +','\t',line)
		if line[0:1]!='#':
			seq_id=line.split("\t")[0]

			interpro_result=Interpro_result(line.split("\t")[0])
			interpro_result.p_value=float(line.split("\t")[12])
			pfam_id=line.split("\t")[5].split(".")[0]
			if (pfam_id in pfam2go.keys()) == bool(1):
				interpro_result.GO_terms=pfam2go[pfam_id]
			interpro_result.name=line.split("\t")[6]+" "+line.split("\t")[7]+" protein"
			
			if interpro_result.p_value<=p_value:
				seq_id2interpro[seq_id].append(interpro_result)
	return seq_id2interpro


def add_full_path_to_protein(sequences):

	go_graph=GO_graph()
	for seq_id in sequences.keys():


		sequences[seq_id].most_downstream_GO_terms=select_most_downstream_gos(go_graph,sequences[seq_id].GO_terms)		
		sequences[seq_id].GO_terms=add_full_go_path(go_graph,sequences[seq_id].GO_terms)

	
	return sequences


def read_pfam2go():
	pfam2go={}
	input_file = [n for n in open('script/pfam2go','r').read().replace("\r","").split("\n") if len(n)>0]
	for line in input_file:
		if line[0:1]!='!':
			pfam_id=line.split('Pfam:')[1].split(' ')[0]
			if (pfam_id in pfam2go.keys())==bool(0):
				pfam2go[pfam_id]=[]
			pfam2go[pfam_id].append('GO:'+line.split('; GO:')[1])	

	return pfam2go	
	
		
def add_full_go_path(go_graph,GO_terms):
	new_go_terms=[]
	for go in GO_terms:		
		if go!='':			
			for basic_go in ("GO:0008150","GO:0003674","GO:0005575"):
				if go_graph.obsolete[go]!="true": # go with no newer version are excluded
					gos_of_all_paths=list()
					for path in nx.all_simple_paths(go_graph.graph,go_graph.obsolete[go],basic_go):
						gos_of_all_paths=list(set(gos_of_all_paths).union(set(path)))
						
					new_go_terms=list(set(new_go_terms).union(set(gos_of_all_paths)))

	return new_go_terms
	
def select_most_downstream_gos(go_graph,GO_terms):
	most_downstream_go_terms=[]
	for go in GO_terms:
		path_to_other_go_exists=bool(0)
		for go_compare in GO_terms:
			if go!='' and go_compare!='' and go!=go_compare:
				if nx.has_path(go_graph.graph,go_graph.obsolete[go],go_graph.obsolete[go_compare]):	
					path_to_other_go_exists=bool(1)
				
		if path_to_other_go_exists==bool(0):
			most_downstream_go_terms.append(go)
		
	return most_downstream_go_terms


def cal_HSSP_curve(hsp):
	hssp_distance=0
	align_length=float(hsp.align_length-hsp.gaps)
	idendity=float(hsp.identities)/float(align_length)*float(100)
	if align_length<=11:
		hssp_distance=idendity-100
	elif align_length>450:
		hssp_distance=idendity-19.5
	else:
		hssp_distance=idendity-480*pow(align_length,(-0.32*(1+(math.exp(-align_length/1000)))))

	return hssp_distance


class Hsp:
	def __init__(self):
		self.hit_def=""
		self.gos=""
		self.name = ""
		self.align_length=""
		self.gaps=""
		self.identities=""
		self.score=""
		self.e_value=1
		
class Sequence:

	def __init__(self, seq_id):
		self.seq_id = seq_id
		self.sequence = ""
		self.GO_terms = list()
		self.most_downstream_GO_terms = list()
		self.ec_numbers = list()
		self.metacyc_reactions = list()
		self.kegg_reactions = list()
		self.name = "hypothetical protein"
		self.prokka_name = ""		





	def most_downstream_GO_string(self):
		to_print=""
		for GO_term in set(self.most_downstream_GO_terms):
			to_print+=", "+GO_term

		if len(to_print)>0:
			to_print=to_print[2:]
		return to_print

	def GO_string(self):
		to_print=""
		for GO_term in set(self.GO_terms):
			to_print+=", "+GO_term

		if len(to_print)>0:
			to_print=to_print[2:]
		return to_print

	def ec_string(self):
		go_graph=GO_graph()
		to_print=""

		for GO_term in set(self.GO_terms):
			
			if (GO_term in go_graph.go_to_ec.keys()) == bool(1):
				for ec in go_graph.go_to_ec[GO_term]:
					to_print+=", "+"EC:"+ec
		if len(to_print)>0:
			to_print=to_print[2:]
		return to_print

	def metacyc_string(self):

		go_graph=GO_graph()
		to_print=""

		for GO_term in set(self.GO_terms):
			
			if (GO_term in go_graph.go_to_metacyc.keys()) == bool(1):
				for metacyc in go_graph.go_to_metacyc[GO_term]:
					to_print+=", "+metacyc
		if len(to_print)>0:
			to_print=to_print[2:]
		return to_print

	def kegg_string(self):

		go_graph=GO_graph()
		to_print=""

		for GO_term in set(self.GO_terms):
			
			if (GO_term in go_graph.go_to_kegg.keys()) == bool(1):
				for kegg in go_graph.go_to_kegg[GO_term]:
					to_print+=", "+kegg
		if len(to_print)>0:
			to_print=to_print[2:]
		return to_print


class Interpro_result:
	def __init__(self, seq_id):
		self.seq_id = seq_id
		self.database = ""
		self.p_value = 1
		self.GO_terms = list()
		self.name = ""	


class GO_graph:	

	ec_to_go={}
	go_to_kegg={}
	go_to_metacyc={}
	go_to_ec={}
	graph={}
	obsolete={}

	def __init__(self):
		if len(GO_graph.ec_to_go.keys())==0:

			obo_file = [n for n in open("script/go.obo",'r').read().replace("\r","").split("\n") if len(n)>0]

			edges=list()
			go_names={}
			name_space={}
			obsolete={}
	
			go_to_kegg={}
			go_to_ec={}
			ec_to_go={}
			go_to_metacyc={}
			for line in obo_file:

				#------ec and metacyc ids
				if line[0:9]=="xref: EC:":
					go_to_ec[go_term].append(line[9:])
					if (line[9:] in ec_to_go.keys()) == bool(0):
						ec_to_go[line[9:]]=list()
					ec_to_go[line[9:]].append(go_term)

				if line[0:14]=="xref: MetaCyc:":	
					go_to_metacyc[go_term].append(line[14:])

				if line[0:11]=="xref: KEGG:":	
					go_to_kegg[go_term].append(line[11:])	



				#------ec and metacyc ids

				if line[0:4]=="id: ":
					alt_ids=list()
					go_term = line.split(" ")[1]
					obsolete[go_term]=go_term

					go_to_ec[go_term]=list()
					go_to_metacyc[go_term]=list()
					go_to_kegg[go_term]=list()
				if line[0:11]=="namespace: ": 
					name_space_name=line[11:]
					name_space[go_term]=name_space_name
				if line[0:6]=="name: ":
					name=line[6:]
					go_names[go_term]=name

				if line[0:8]=="alt_id: ":
					alt_go=line[8:]
					obsolete[alt_go]=go_term

				if line[0:6]=="is_a: ":
					edge=(go_term,line.split(" ")[1])
					edges.append(edge)
					for z in alt_ids:
						edge=(z,line.split(" ")[1])
						edges.append(edge)
			
				if line[0:17]=="is_obsolete: true": # with this obsolete go terms with no replacements can be replaced
					obsolete[go_term]="true"
				if line[0:13]=="replaced_by: ":
					obsolete[go_term]=line[13:]
				if line[0:10]=="consider: ":
					obsolete[go_term]=line[10:]
		
			graph=nx.DiGraph()
			graph.add_edges_from(edges)

			GO_graph.ec_to_go=ec_to_go
			GO_graph.go_to_kegg=go_to_kegg
			GO_graph.go_to_metacyc=go_to_metacyc
			GO_graph.go_to_ec=go_to_ec
			GO_graph.graph=graph
			GO_graph.obsolete=obsolete
main()		
