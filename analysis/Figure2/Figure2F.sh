#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source trapdebug
scrnah5adrawcounts2bulkby -d "$outdir" -b BC5A_D -k sampleid -k species -k celltype -H bulkid -- BC5A_D.h5ad
scrnah5adhvg2metaneighborus -d "$outdir" -b BC5A_D -s sampleid -u species -c celltype -W 6 -H 6 -m 5 -- BC5A_D.h5ad
