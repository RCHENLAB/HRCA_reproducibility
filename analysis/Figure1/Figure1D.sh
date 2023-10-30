#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source trapdebug
scrnah5adfiles2scviwkfl -d "$outdir" -e u_scvi -t 2 -c Figure1D_config.yaml -- snRNA.h5ad
