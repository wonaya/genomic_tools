#!/bin/bash

wavelets="/work/02114/wonaya/software/hotspot-distr-v3/hotspot-deploy/bin/wavelets"

usage="\n
waveTime <signal_bedgraph> <level> <chromosome_file>\n
<signal_bedgraph> is the only required argument.\n
<signal_bedgraph> is a 4-column bedGraph file (signal in column 4)\n
<chromosome_file> contains chromosomes to be processed, in first column. By default, use all chromosomes in the signal file.\n
<level> level of smoothing; by default 3. If the resolution of the input file is x, then the results are \n
        smoothed out to a scale of (2^lvl)*x.\n"
        
if [ $# -gt 4 ] || [ $# -lt 2 ]; then
    echo -e $usage
    exit
fi

# Our temp directory
tmpd=./tmp$$

# Density file (signal_bedgraph)
den_orig=$1
den=`basename $den_orig .bed`.sorted.wig
sort -k1,1 -k2,2g -o $den $den_orig

## Default values.
#filter=LA8
filter=Haar
halfsite=75
lvl=$2
echo $3
chrs=$(cut -f1 $3)
pid=$$
mkdir -p $tmpd
outpl=""
for chr in $chrs
do
    echo $chr 1>&2
    egrep $chr $den > $tmpd/${chr}.den
    tmptxt=$tmpd/tmp.`basename $den .bed`.$chr.J$lvl.$pid.txt
    tmpsmth=`echo $tmptxt | sed s/txt$/smooth.txt/`
	## Get scores to smooth from 4th column of bedGraph file
    cut -f4 $tmpd/${chr}.den > $tmptxt
    $wavelets --level $lvl --to-stdout --boundary reflected --filter $filter $tmptxt > $tmpsmth
done

cat $tmpd/*.smooth.txt > `basename $den`.J$lvl.smoothed.txt
rm -Rf $tmpd
paste $den `basename $den`.J$lvl.smoothed.txt | cut -f 1,2,3,5 > `basename $den`.J$lvl.$filter.smoothed.wig
