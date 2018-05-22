#!/bin/bash
while getopts :r: TEST; do
	case $TEST in
	r) OPT_R=$OPTARG
	;;	
	esac
done
mkdir -p GenomeStats
GenomeStats.py -I $OPT_R > GenomeStats/$(basename $OPT_R)".stats.txt"
mdust $OPT_R > GenomeStats/$(basename $OPT_R)".dust.fasta"
repeat-match -n 50 GenomeStats/$(basename $OPT_R)".dust.fasta" > GenomeStats/$(basename $OPT_R)".repeat.fasta"
cd GenomeStats
GenomeStats.py -I $(basename $OPT_R)".dust.fasta" > $(basename $OPT_R)"genomestats_afterdust.txt"
cd ../
