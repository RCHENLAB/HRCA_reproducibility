#!/bin/sh
dir="human_meta/scripts/finemap_eQTL"
software/R-4.0.0/bin/R --no-save < ${dir}/5pip_plot_merged.R > ${dir}/5pip_plot_merged.out 2> ${dir}/5pip_plot_merged.err
