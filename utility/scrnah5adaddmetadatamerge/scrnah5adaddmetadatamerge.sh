#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

args=$(getopt -o hd:b:m:k:c:e: -l help,outdir:,bname:,metadata:,key:,character:,condaenv: --name "$0" -- "$@") || exit "$?"
eval set -- "$args"

absdir=$(dirname $(readlink -f "$0"))
scriptname=$(basename "$0" .sh)

outdir=.
while true
do
	case "$1" in
		-h|--help)
			exec cat "$absdir/$scriptname.txt"
			;;
		-d|--outdir)
			outdir=$2
			shift 2
			;;
		-b|--bname)
			bname=$2
			shift 2
			;;
		-m|--metadata)
			metadata=$2
			shift 2
			;;
		-k|--key)
			key=$2
			shift 2
			;;
		-c|--character)
			character+=("$2")
			shift 2
			;;
		-e|--condaenv)
			condaenv=$2
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

[[ $bname ]] || { echo "$scriptname.sh: -b|--bname must be specified!"; exit 1; }
[[ $metadata ]] || { echo "$scriptname.sh: -m|--metadata must be specified!"; exit 1; }
[[ $key ]] || { echo "$scriptname.sh: -k|--key must be specified!"; exit 1; }

function cmd {
if [[ $condaenv ]]
then
	source "$(conda info --base)/etc/profile.d/conda.sh"
	conda activate "$condaenv"
fi

local f=$1
set -x
mkdir -p "$outdir" && pycmd.sh \
	-e "f='$f'" \
	-e "outdir='$outdir'" \
	-e "bname='$bname'" \
	-e "metadata='$metadata'" \
	-e "key='$key'" \
	-e "character=$(basharr2pylist.sh -c -- "${character[@]}")" \
	-s "$absdir/python/$scriptname.py"
}

if (($#))
then
	cmd "$@"
fi
