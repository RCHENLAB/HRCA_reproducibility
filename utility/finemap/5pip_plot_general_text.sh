#!/bin/sh
command=human_meta/scripts/finemap/5pip_plot_general_text
file=$1
software/R-4.0.0/bin/Rscript --vanilla  ${command}.R  $file
