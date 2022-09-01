#!/bin/bash

#This is a script to run the first couple of steps of Change-O in a folder containing Cell Ranger results (in separate folders). Should be generalizable for now.
root=$(pwd)
cd from_cellranger
for d in *; do
	#Preparing AIRR-formatted data:
	cd $root/from_cellranger/$d/vdj_b;
	AssignGenes.py igblast -s filtered_contig.fasta -b ~/yard/share/igblast/ --organism human --loci ig --format airr;
	MakeDb.py igblast -i filtered_contig_igblast.fmt7 -s filtered_contig.fasta -r ~/yard/share/germlines/imgt/human/vdj/imgt_human_*.fasta --10x filtered_contig_annotations.csv --extended;

	#Separating heavy and light chains:
	ParseDb.py select -d filtered_contig_igblast_db-pass.tsv -f locus -u "IGH" --logic all --regex --outname heavy
	ParseDb.py select -d filtered_contig_igblast_db-pass.tsv -f locus -u "IG[LK]" --logic all --regex --outname light 

done
cd $root
exit 0;

