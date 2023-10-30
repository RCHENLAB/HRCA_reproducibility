#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

args=$(getopt -o ho:W:H:m:f -l help,outfile:,width:,height:,main:,fixorder --name "$0" -- "$@") || exit "$?"
eval set -- "$args"

absdir=$(dirname $(readlink -f "$0"))
scriptname=$(basename "$0" .sh)

width=5
height=4
main=
fixorder=F
while true
do
	case "$1" in
		-h|--help)
			exec cat "$absdir/$scriptname.txt"
			;;
		-o|--outfile)
			outfile=$2
			shift 2
			;;
		-W|--width)
			width=$2
			shift 2
			;;
		-H|--height)
			height=$2
			shift 2
			;;
		-m|--main)
			main=$2
			shift 2
			;;
		-f|--fixorder)
			fixorder=T
			shift
			;;
		--)
			shift
			break
			;;
		*)
			echo "$0: not implemented option: $1" >&2
			exit 1
			;;
	esac
done

[[ $outfile ]] || { echo "$scriptname.sh: -o|--outfile must be specified."; exit 1; }

outdir=$(dirname "$outfile")
bname=$(basename.sh -s .pdf -s .png "$outfile")
mkdir -p "$outdir" && Rvanilla.sh \
	-e "outfile='$outfile'" \
	-e "width=$width" \
	-e "height=$height" \
	-e "main='$main'" \
	-e "fixorder=$fixorder" \
	-e "source('$absdir/R/$scriptname.R')"
