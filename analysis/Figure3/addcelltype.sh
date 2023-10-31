#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source trapdebug
## AC
scrnah5adaddmetadatamerge.sh -d "$outdir" -b snRNA_AC -m AC_celltype_table.txt -k cluster2 -- snRNA_AC_scVI.h5ad
scrnah5adumapby.sh -W 5 -H 5 -l group1 -l group2 -l group3 -l group4 -l group5 -l group6 -l cluster2 -t -d "$outdir" -b "$bname" -- snRNA_AC.h5ad
## RGC
scrnah5adaddmetadatamerge.sh -d "$outdir" -b snRNA_RGC -m RGC_celltype_table.txt -k group1 -c RGC_clusternum -- snRNA_RGC_scVI.h5ad
scrnah5adumapby.sh -W 3.5 -H 3.5 -l group1 -l cluster2 -l RGC_clusternum -t -d "$outdir" -b "$bname" -- snRNA_RGC.h5ad
