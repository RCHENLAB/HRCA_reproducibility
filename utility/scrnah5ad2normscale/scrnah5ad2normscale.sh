#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

args=$(getopt -o hd:b:rse:t: -l help,outdir:,bname:,raw,scale,useraw,condaenv:,targetsum: --name "$0" -- "$@") || exit "$?"
eval set -- "$args"

absdir=$(dirname $(readlink -f "$0"))
scriptname=$(basename "$0" .sh)

outdir=.
raw=False
useraw=False
scale=False
targetsum=None
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
		-r|--raw)
			raw=True
			shift
			;;
		--useraw)
			useraw=True
			shift
			;;
		-s|--scale)
			scale=True
			shift
			;;
		-e|--condaenv)
			condaenv=$2
			shift 2
			;;
		-t|--targetsum)
			targetsum=$2
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

function cmd {
if [[ $condaenv ]]
then
	if [[ $MAMBA_ROOT_PREFIX ]]
	then
		source "$MAMBA_ROOT_PREFIX/etc/profile.d/micromamba.sh"
		micromamba activate "$condaenv"
	else
		source "$(conda info --base)/etc/profile.d/conda.sh"
		conda activate "$condaenv"
	fi
fi

set -xe
local f=$1
hdf5ls.sh "$f"
mkdir -p "$outdir" && exec pycmd.sh \
	-e "f='$f'" \
	-e "outdir='$outdir'" \
	-e "bname='$bname'" \
	-e "raw=$raw" \
	-e "useraw=$useraw" \
	-e "scale=$scale" \
	-e "targetsum=$targetsum" \
	-s "$absdir/python/$scriptname.py"
}

if (($#))
then
	cmd "$@"
fi
