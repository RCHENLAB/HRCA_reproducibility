#!/bin/sh

main="/storage/chen/home/jw29"

export PYTHONPATH="/storage/chen/home/jw29/software/anaconda3/lib/python3.8/site-packages:/storage/chen/home/jw29/software/anaconda3/lib/python3.8:/storage/chen/home/jw29/software/anaconda3/bin"
#date="4-12-23"
#spe="mm"
date=$1
spe=$2
${main}/software/anaconda3/bin/python enhancer_validation/scripts/9barcode_mean.py $date $spe

#spe="hm"
#${main}/software/anaconda3/bin/python enhancer_validation/scripts/9barcode_mean.py $date $spe
