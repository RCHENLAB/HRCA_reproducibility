#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source trapdebug
zcat Figure2E_cast.txt.gz | cut -f 2- | rplotvenndiagramovlp.sh -o Figure2E.pdf --cast -k 1 -W 3.5 -H 3.5
