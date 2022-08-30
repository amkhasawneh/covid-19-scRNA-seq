#!/bin/bash

#This is a script to run the first couple of steps of Change-O in a folder containing Cell Ranger results (in separate folders). Should be generalizable for now.
root=$(pwd)
for d in *; do
	cd $root/$d/outs/per_sample_outs/*/vdj_b;
	AssignGenes.py igblast -s filtered_contig.fasta -b ~/yard/share/igblast/ --organism human --loci ig --format  blast;
	MakeDb.py igblast -i filtered_contig_igblast.fmt7 -s filtered_contig.fasta -r ~/yard/share/germlines/imgt/human/vdj/imgt_human_*fasta --10x filtered_contig_annotations.csv --extended;
done

exit 0;

