#!/bin/bash

set -eux;

if [ ! -d from_cellranger ]; then mkdir from_cellranger; fi

for i in [^fs]*[^.sh$];	do
	mkdir from_cellranger/"$i"; 
	mkdir from_cellranger/"$i"/{count,vdj_{b,t}}; 
	cp -vr "$i"/outs/per_sample_outs/"$i"/count/sample_filtered_feature_bc_matrix from_cellranger/"$i"/count/; 
	cp -v "$i"/outs/per_sample_outs/"$i"/count/sample_molecule_info.h5 from_cellranger/"$i"/count/; 
	cp -v "$i"/outs/multi/count/raw_molecule_info.h5 from_cellranger/"$i"/count/; 
	cp -v "$i"/outs/multi/count/raw_feature_bc_matrix.h5 from_cellranger/"$i"/count/; 
	cp -v "$i"/outs/per_sample_outs/"$i"/web_summary.html from_cellranger/"$i"/;
	for j in t b; do
		cp -v "$i"/outs/per_sample_outs/"$i"/vdj_$j/airr_rearrangement.tsv from_cellranger/"$i"/vdj_$j/; 
		cp -v "$i"/outs/per_sample_outs/"$i"/vdj_$j/filtered_contig_annotations.csv from_cellranger/"$i"/vdj_$j/; 
		cp -v "$i"/outs/per_sample_outs/"$i"/vdj_$j/clonotypes.csv from_cellranger/"$i"/vdj_$j/;
	done


done

exit 0;
