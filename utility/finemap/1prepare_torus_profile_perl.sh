#!/bin/sh

export PERL5LIB=/storage/chen/home/jw29/software/perl/lib/perl5/Math-CDF-0.1/lib/perl5/x86_64-linux/:$PERL5LIB
dir="human_meta/scripts/finemap"

/storage/chen/Software/perl-5.10.1/perl ${dir}/1prepare_torus_profile.pl $1  $2 $3 $4
