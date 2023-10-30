#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source trapdebug
saturncrossspecieswkfl -d "$outdir" -e u_saturn -t 4 -c Figure2C_config.yaml -- Figure2C_saturn.yaml

