#!/bin/bash

while getopts :a:1:2: TEST; do
  case $TEST in

    #Aligner
    a) OPT_A=$OPTARG
    ;;
    #Paired Reads one
    1) OPT_1=$OPTARG
    ;;
    #Paired Reads two
    2) OPT_2=$OPTARG
    ;;
    esac
done
#Get pwd for starting directory
pwd=`pwd`


#Check that both files are given to command line
if [ -z $OPT_1 ] || [ -z $OPT_2 ] 
then 
	printf "${bold} One or more input files not given to command line with [-1] or [-2]"
	exit 1
fi

#Check that both provided files exist in directory.
if [ -e $OPT_1 ]  
then
	:
else
	echo "${bold} Input file one [-1] not detected in directory"
	exit 1
fi

if [ -e $OPT_2 ] 
then
	:
else 
	echo "${bold} Input file two [-2] not detected in current working directory"
	exit 1
fi

if [ "$OPT_A" = "spades" ] || [ "$OPT_A" = "iva" ] || [ "$OPT_A" = "trinity" ] || [ "$OPT_A" = "allpaths" ] || [ "$OPT_A" = "celera" ] || [ "$OPT_A" = "abyss" ] || [ "$OPT_A" = "vicuna" ] || [ "$OPT_A" = "idba" ] || [ "$OPT_A" = "test" ] || [ "$OPT_A" = "velvet" ]
	then
		:	
	else
		echo "The spelling of your aligner looks wrong"
		echo "Try [spades] / [iva] / [trinity] / [allpath] / [celera] / [abyss] / [vicuna] / [idba] / [velvet]"
		exit 1
	fi

pythContigs(){
	Contigs2OneLine.py $refContig
	echo ${refContig%.fa*}
	ln -s ${refContig%.fa*}"_oneline.fa" ./
	RemoveContigs.py *.fa
}

alignReads(){
	mkdir -p BowtieIndexes
	mkdir -p BowtieOutput
	bowtie2-build *contig*_oneline.fa BowtieIndexes/${PWD##*/}
	bowtie2 -p 15 -x BowtieIndexes/${PWD##*/} -1 ../../$OPT_1 -2 ../../$OPT_2 -S BowtieOutput/${PWD##*/}"_aligned.sam" 2>&1 | tee ${PWD##*/}"_bowtie.stats"
	refLen=$(awk '{print $2}' ../${pwd##*/}_weesam.stats | grep -v RefLength)
	echo "!!!!! $refLen !!!!"
}

samtobam(){
	echo "Sorting sam file: BowtieOutput/'${PWD##*/}'_aligned.sam"
	samtools view -bS BowtieOutput/${PWD##*/}"_aligned.sam" | samtools sort -o BowtieOutput/${PWD##*/}"_aligned.bam"
	rm -f BowtieOutput/${PWD##*/}"_aligned.sam"

	samtools index BowtieOutput/${PWD##*/}"_aligned.bam"

	weeSAMv1.3 -b BowtieOutput/${PWD##*/}"_aligned.bam" -out BowtieOutput/${PWD##*/}"_weesam.stats"
}

weeSamStats(){
	#Always call after bowtie method, need a variablei
	echo "Running weesam with a ref length of $refLen"
	WeeSamContigs.py BowtieOutput/${PWD##*/}"_weesam.stats" $refLen
	mv WeeSamContigs_Stats.txt BowtieOutput/
}

blast(){
	mkdir BlastOutput
	blastn -query *contig*_oneline.fa -evalue 1e-5 -num_alignments 1 -num_threads 12 -db nt -out BlastOutput/${PWD##*/}_blast.txt
}

##	Spades	##
if [ $OPT_A = "spades" ]
then
	mkdir -p Alignment/SpadesContigs
	cd Alignment/SpadesContigs
	refContig=$pwd/SpadesOutput/contigs.fasta
	pythContigs
	alignReads
	samtobam
	weeSamStats
	blast
fi
##	 IDBA	  ##
if [ $OPT_A = "idba" ]
then
	mkdir -p Alignment/IDBAContigs
	cd Alignment/IDBAContigs
	refContig=$pwd/IDBAOutput/contig.fa
	pythContigs
	alignReads
	samtobam
	weeSamStats
	blast
fi
#	Abyss
if [ $OPT_A = "abyss" ]
then
	mkdir -p Alignment/ABySSContigs
	cd Alignment/ABySSContigs
	refContig=$pwd/ABySSOutput/${pwd##*/}-contigs.fa
	pythContigs
	alignReads
	samtobam
	weeSamStats
	blast
fi

#VICUNA
if [ $OPT_A = "vicuna" ]
then
	mkdir -p Alignment/VicunaContigs
	cd Alignment/VicunaContigs
	refContig=$pwd/VicunaOutput/contig.fasta
	pythContigs
	alignReads
	samtobam
	weeSamStats
	blast
fi

#IVA
if [ $OPT_A = "iva" ] 
then
	mkdir -p Alignment/IVAContigs
	cd Alignment/IVAContigs
	refContig=$pwd/IVAOutput/contigs/contigs.fasta
	pythContigs
	alignReads
	samtobam
	weeSamStats
	blast
fi 

# Velvet
if [ $OPT_A = "velvet" ] 
then
	mkdir -p Alignment/VelvetContigs
	cd Alignment/VelvetContigs
	refContig=$pwd/${pwd##*/}_VelvetAssembly/contigs.fa
	pythContigs
	alignReads
	samtobam
	weeSamStats
	blast
fi

