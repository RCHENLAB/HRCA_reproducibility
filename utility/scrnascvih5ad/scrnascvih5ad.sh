#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

args=$(getopt -o hd:b:k:e:l:t:n:f:s:p:g: -l help,outdir:,bname:,batchkey:,condaenv:,nlayer:,nlatent:,ntop:,flavor:,seed:,epoch:,normcounts,gpu: --name "$0" -- "$@") || exit "$?"
eval set -- "$args"

absdir=$(dirname $(readlink -f "$0"))
scriptname=$(basename "$0" .sh)

outdir=.
nlayer=2
nlatent=30
ntop=2000
flavor=seurat
seed=12345
epoch=None
normcounts=False
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
		-k|--batchkey)
			batchkey=$2
			shift 2
			;;
		-e|--condaenv)
			condaenv=$2
			shift 2
			;;
		-l|--nlayer)
			nlayer=$2
			shift 2
			;;
		-t|--nlatent)
			nlatent=$2
			shift 2
			;;
		-n|--ntop)
			ntop=$2
			shift 2
			;;
		-f|--flavor)
			flavor=$2
			shift 2
			;;
		-s|--seed)
			seed=$2
			shift 2
			;;
		-p|--epoch)
			epoch=$2
			shift 2
			;;
		--normcounts)
			normcounts=True
			shift
			;;
		-g|--gpu)
			gpu=$2
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
[[ $batchkey ]] || { echo "$scriptname.sh: -k|--batchkey must be specified!"; exit 1; }

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

if ! [[ $gpu ]]
then
	host=$(hostname)
	if [[ $host == mhgcp-g00.grid.bcm.edu ]]
	then
		gpu=$((RANDOM%2))
	elif [[ $host == mhgcp-g01.grid.bcm.edu ]]
	then
		gpu=$((RANDOM%4))
	else
		ngpu=$(ngpu.sh)
		gpu=$((RANDOM%${ngpu}))
	fi
fi
export CUDA_VISIBLE_DEVICES=$gpu

set -xe
local f=$(abspath.sh "$1")
hdf5ls.sh "$f"
mkdir -p "$outdir" && cd "$outdir" && pycmd.sh \
	-e "f='$f'" \
	-e "bname='$bname'" \
	-e "batchkey='$batchkey'" \
	-e "nlayer=$nlayer" \
	-e "nlatent=$nlatent" \
	-e "ntop=$ntop" \
	-e "flavor='$flavor'" \
	-e "seed=$seed" \
	-e "epoch=$epoch" \
	-e "normcounts=$normcounts" \
	-s "$absdir/python/$scriptname.py"
}

if (($#))
then
	cmd "$@"
fi
