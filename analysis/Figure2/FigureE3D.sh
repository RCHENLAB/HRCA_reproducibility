#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source trapdebug
scrnah5adrawcounts2bulkby -d "$outdir" -b HRCA_snRNA_BC -k donor -k celltype -H bulkid -f sampleid -c 2000 --HRCA_snRNA_BC_rawcounts.h5ad 
scrnah5ad2degcontrastpairbydeseq2 -e deseq2 -d "$outdir" -b HRCA_snRNA_BC -k celltype -1 GB -2 BB -p donor -r 10 -H 5 -W 8 -s 100 -- HRCA_snRNA_BC.h5ad
