#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

args=$(getopt -o hd:b:W:H:q:l:s:p: -l help,outdir:,bname:,width:,height:,qthr:,log2fcthr:,labelfc:,labelqvalue:,labelsymbol:,xlim:,ylim:,selectlab:,legendposition:,pointsize:,labsize: --name "$0" -- "$@") || exit "$?"
eval set -- "$args"

absdir=$(dirname $(readlink -f "$0"))
scriptname=$(basename "$0" .sh)

outdir=.
width=8
height=8
qthr=0.05
log2fcthr=0.1
labelfc=log2FoldChange
labelqvalue=padj
labelsymbol=symbol
ylim=-1
xlim=-10
legendposition=top
pointsize=1
labsize=3
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
		-W|--width)
			width=$2
			shift 2
			;;
		-H|--height)
			height=$2
			shift 2
			;;
		-q|--qthr)
			qthr=$2
			shift 2
			;;
		-l|--log2fcthr)
			log2fcthr=$2
			shift 2
			;;
		--labelfc)
			labelfc=$2
			shift 2
			;;
		--labelqvalue)
			labelqvalue=$2
			shift 2
			;;
		--labelsymbol)
			labelsymbol=$2
			shift 2
			;;
		--xlim)
			xlim=$2
			shift 2
			;;
		--ylim)
			ylim=$2
			shift 2
			;;
		-s|--selectlab)
			selectlab+=("$2")
			shift 2
			;;
		-p|--legendposition)
			legendposition=$2
			shift 2
			;;
		--pointsize)
			pointsize=$2
			shift 2
			;;
		--labsize)
			labsize=$2
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

[ "$bname" ] || { echo "$scriptname.sh: -b|--bname must be specified!"; exit 1; }

function cmd {
source traptimes
SECONDS=0
local f=$1
set -ex

mkdir -p "$outdir" && Rvanilla.sh \
	-e "f='$f'" \
	-e "outdir='$outdir'" \
	-e "bname='$bname'" \
	-e "width=$width" \
	-e "height=$height" \
	-e "qthr=$qthr" \
	-e "log2fcthr=$log2fcthr" \
	-e "labelfc='$labelfc'" \
	-e "labelqvalue='$labelqvalue'" \
	-e "labelsymbol='$labelsymbol'" \
	-e "xlim=$xlim" \
	-e "ylim=$ylim" \
	-e "selectlab=$(basharr2Rvec.sh -c -- "${selectlab[@]}")" \
	-e "legendposition='$legendposition'" \
	-e "pointsize=$pointsize" \
	-e "labsize=$labsize" \
	-e "source('$absdir/R/$scriptname.R')"
set +x
echo "Time elapsed: $(date -u -d @$SECONDS +'%H:%M:%S')"
}

if (($#))
then
	cmd "$@"
fi
