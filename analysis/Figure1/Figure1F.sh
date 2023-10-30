#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source trapdebug
scrnah5adrawcounts2bulkby -d "$outdir" -b HRCA_snRNA_scRNA_bulk -k sampleid -k dataset -k majorclass -H bulkid -- HRCA_snRNA_scRNA.h5ad
scrnah5ad2bulkdegcondtypebydeseq2 -e deseq2 -d "$outdir" -b HRCA_snRNA_scRNA -k dataset -1 snRNA -2 scRNA -l majorclass -r 10 -c 10 -H 5 -W 8 -s 100 -- HRCA_snRNA_scRNA_bulk.h5ad
deseq2enhancedvolcano.sh -H 4 -W 6 -q 0.05 -l 1 -d "$outdir" -b "$bname" --xlim 10 --ylim 300 -s 'RHO' -p right -- HRCA_snRNA_scRNA_full_bylabel_dataset_snRNAVsscRNA.txt.gz
