#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source trapdebug
scrnah5adrawcounts2bulkby -d "$outdir" -b HRCA_snRNA_scRNA_bulk -k sampleid -k dataset -k majorclass -H bulkid -- HRCA_snRNA_scRNA.h5ad
scrnah5adhvg2metaneighborus -d "$outdir" -b HRCA_snRNA_scRNA_metaneighbor -s sampleid -u dataset -c majorclass -W 6 -H 6 -m 5 -- HRCA_snRNA_scRNA_bulk.h5ad
