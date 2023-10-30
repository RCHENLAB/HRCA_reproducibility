#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

args=$(getopt -o ho:W:H:m:k: -l help,outfile:,width:,height:,main:,melt,cast,keyword: --name "$0" -- "$@") || exit "$?"
eval set -- "$args"

absdir=$(dirname $(readlink -f "$0"))
scriptname=$(basename "$0" .sh)

width=5
height=5
main=
mode=melt
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
		--melt)
			mode=melt
			shift
			;;
		--cast)
			mode=cast
			shift
			;;
		-k|--keyword)
			keyword=$2
			shift 2
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

[[ $outfile ]] || { echo "$scriptname.sh: -o|--outfile must be specified!"; exit 1; }

if [[ $mode = cast ]]
then
	[[ $keyword ]] || { echo "$scriptname.sh: -k|--keyword is required for the --cast mode!"; exit 1; }
fi

mkdir -p "$(dirname "$outfile")" && exec R --slave --no-save --no-restore --no-init-file \
	-e "width=$width" \
	-e "height=$height" \
	-e "outfile='$outfile'" \
	-e "main='$main'" \
	-e "keyword='$keyword'" \
	-e "source('$absdir/R/$mode/$scriptname.R')"
