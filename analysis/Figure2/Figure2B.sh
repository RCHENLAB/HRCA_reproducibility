#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source trapdebug
scrnah5adumapby.sh -W 4 -H 3 -l group1 -l group2 -l group3 -t -d "$outdir" -b "$bname" -- snRNA_BC_scVI.h5ad
