#!/usr/bin/env python

import numpy as np
import subprocess
import sys
import os
import sys
import gzip

inputOptions = sys.argv[1:]

#usage: file1


def main():

	info={}
	info["sequence"]=""
	info["prot_ID"]=""
	info["RefSeq"]=""
	info["NCBI_TaxID"]=""
	info["organism"]=""
	info["sequence_length"]=0
	info["protein_name"]=""
	info["GO"]=""
	info["KEGG"]=""
	info["EC"]=""

 
	with gzip.open(inputOptions[0],'r') as f:
    		for line in f:
			construct_fasta(line.replace("\n",""),info)
			

def construct_fasta(line,info):
	
	if line[0:2]=="DE":
		if line.find("Full=")!=-1 and info["protein_name"]=="":
			info["protein_name"]+=line.split("Full=")[1].split("{")[0].replace(";","")
		if line.find("EC=")!=-1:
			info["EC"]+=line.split("EC=")[1].split(" ")[0].split(";")[0]+";"

	if line[0:23]=="GN   OrderedLocusNames=":
		to_remove_part=line.split("=")[1].split(";")[0].split(" {")[0]
		to_remove_part=" "+to_remove_part
		info["protein_name"]=info["protein_name"].replace(to_remove_part,"")


	if line[0:10]=="DR   KEGG;":
		if info["KEGG"]=="":
			info["KEGG"]+=line.split("DR   KEGG; ")[1].split(";")[0]+";"
	if line[0:2]=="ID":
		info["prot_ID"]+=line.split("   ")[1]
	if line[0:8]=="DR   GO;":
		info["GO"]+=line.split("GO; ")[1].split(";")[0]+";"
	if line[0:13]=="DR   RefSeq; ":
		if info["RefSeq"]=="":
			info["RefSeq"]+=(line.split("DR   RefSeq; ")[1]+";").replace(".;",";").replace("; ",";")
	if line[0:2]=="OX":
		info["NCBI_TaxID"]+=line.split("NCBI_TaxID=")[1].replace(";","").split(" ")[0]
	if line[0:2]=="OS":
		info["organism"]+=line.split("   ")[1].split(" (")[0]
	if line[0:2]=="SQ":
		info["sequence_length"]+=int(line.split("   ")[2].split(" ")[0])
	if line[0:5]=="     ":
		info["sequence"]+=line.replace(" ","")
	if line[0:2]=="//":
		assert(info["sequence_length"]==len(info["sequence"]))
		

		info_string="|name="+info["protein_name"]+"|species="+info["organism"]+"|NCBI_TaxID="+info["NCBI_TaxID"]+ "|GO="+info["GO"]+ "|RefSeq="+info["RefSeq"]+ "|KEGG="+info["KEGG"]+ "|EC="+info["EC"]+"|seq_length="+str(info["sequence_length"])+"|"

		print ">"+info["prot_ID"]+ " "+info_string.replace(";|","|").replace(".|","|").replace(" |","|")
		i=0
		print_string=""
		while (i<=len(info["sequence"])):
			print_string+= info["sequence"][i:i+60]+"\n"
			i=i+60
		assert(info["sequence_length"]==len(print_string.replace("\n","")))

		print print_string

		#cleaning up
		for key in info.keys():
			info[key]=""
		info["sequence_length"]=0


main()
