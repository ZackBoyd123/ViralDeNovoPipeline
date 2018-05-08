#!/bin/bash

while getopts :a:1:2:r: TEST; do
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
    r) OPT_R=$OPTARG
    ;;
    esac
done
#Get pwd for starting directory
pwd=`pwd`
substr=,
if [ "${OPT_A/$substr}" = "$OPT_A" ] 
then
	:
else
	echo "Multiple Aligners Selected"
	IFS=',' read -a array <<< $OPT_A
fi
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

if [ -n $OPT_A ] && [ -z $array ] 
then
	if [ "$OPT_A" = "spades" ] || [ "$OPT_A" = "iva" ] || [ "$OPT_A" = "trinity" ] || [ "$OPT_A" = "allpaths" ] || [ "$OPT_A" = "celera" ] || [ "$OPT_A" = "abyss" ] || [ "$OPT_A" = "vicuna" ] || [ "$OPT_A" = "idba" ] || [ "$OPT_A" = "test" ] || [ "$OPT_A" = "velvet" ]
		then
			:	
		else
			echo "The spelling of your aligner looks wrong"
			echo "Try [spades] / [iva] / [trinity] / [allpath] / [celera] / [abyss] / [vicuna] / [idba] / [velvet]"
			exit 1
		fi
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
	refLen=$(cat ../${pwd##*/}_weesam.stats | grep -v RefLength | awk '{sum += $2} END {print sum}')
	#refLen=$(awk '{print $2}' ../${pwd##*/}_weesam.stats | grep -v RefLength)
	echo "!!!!! $refLen !!!!"
}

samtobam(){
	echo "Sorting sam file: BowtieOutput/'${PWD##*/}'_aligned.sam"
	samtools view -bS BowtieOutput/${PWD##*/}"_aligned.sam" | samtools sort -o BowtieOutput/${PWD##*/}"_aligned.bam"
	rm -f BowtieOutput/${PWD##*/}"_aligned.sam"

	samtools index BowtieOutput/${PWD##*/}"_aligned.bam"

	weeSAMv1.3 -b BowtieOutput/${PWD##*/}"_aligned.bam" -out BowtieOutput/${PWD##*/}"_weesam.stats"
	echo "I finished sam2bam stuff at $(date)"
}

weeSamStats(){
	#Always call after bowtie method, need a variablei
        refLen=$(cat ../${pwd##*/}_weesam.stats | grep -v RefLength | awk '{sum += $2} END {print sum}')
	WeeSamContigs.py BowtieOutput/${PWD##*/}"_weesam.stats" $refLen
	mv WeeSamContigs_Stats.txt BowtieOutput/
	echo "I finished WeeSamstuff at $(date)"
}

blast(){
	mkdir BlastOutput
	blastdb=BlastDB
	newblast=$(basename $OPT_R)
	if [ -d ../"$blastdb" ]
	then
		echo "im in here"
		:
	else
		mkdir -p ../$blastdb
		echo "no blast ref db, making one now"
		makeblastdb -in "$OPT_R" -out ../BlastDB/"$newblast" -dbtype nucl
	fi
	blastn -query *contig*_oneline.fa -evalue 1e-5 -num_alignments 1 -num_threads 12 -db ../BlastDB/"$newblast" -out BlastOutput/${PWD##*/}_blast2ref.txt -outfmt 6
	sort -u -k1,1 BlastOutput/${PWD##*/}_blast2ref.txt > BlastOutput/${PWD##*/}_blast2ref_unique.txt
	#refLen=$(cat ../${pwd##*/}_weesam.stats | grep -v RefLength | awk '{sum += $2} END {print sum}')
	blastn -query *contig*_oneline.fa -evalue 1e-5 -num_alignments 1 -num_threads 12 -db nt -out BlastOutput/${PWD##*/}_blastHits.txt -outfmt 6
	sort -u -k1,1 BlastOutput/${PWD##*/}_blastHits.txt > BlastOutput/${PWD##*/}_blastHits_unique.txt
	mkdir -p ../TSV_Files
	ConvertToTSV -1 BowtieOutput/WeeSamConitgs_Stats.txt -2 *contig*_oneline_stats.txt -3 BlastOutput/${PWD##*/}_blast2ref_unique.txt -4 BlastOutput/${PWD##*/}_blastHits_unique.txt
	mv *.tsv ../TSV_Files
	echo "I finished blast stuff at $(date)"
}

##	Spades	##
if [ $OPT_A = "spades" ] || [[ " ${array[@]} " =~ " spades "  ]]
then
	echo "I started in spades $(date)"
	mkdir -p Alignment/SpadesContigs
	cd Alignment/SpadesContigs
	refContig=$pwd/SpadesOutput/contigs.fasta
	#pythContigs
	#alignReads
	#samtobam
	#weeSamStats
	blast
	echo "I finished in spades at $(date)"
	cd ../../
fi
##	 IDBA	  ##
if [ $OPT_A = "idba" ] || [[ " ${array[@]} " =~ " idba " ]]
then
	echo "I started in IDBA $(date)"
	mkdir -p Alignment/IDBAContigs
	cd Alignment/IDBAContigs
	refContig=$pwd/IDBAOutput/contig.fa
	pythContigs
	alignReads
	samtobam
	weeSamStats
	blast
	echo "I finished in IDBA $(date)"
	cd ../../
fi
#	Abyss
if [ $OPT_A = "abyss" ] || [[ " ${array[@]} " =~ " abyss " ]]
then
	echo "I started in abyss $(date)"
	mkdir -p Alignment/ABySSContigs
	cd Alignment/ABySSContigs
	refContig=$pwd/ABySSOutput/${pwd##*/}-contigs.fa
	pythContigs
	alignReads
	samtobam
	weeSamStats
	blast
	echo "I finished in abyss $(date)"
	cd ../../
fi

#VICUNA
if [ $OPT_A = "vicuna" ] || [[ " ${array[@]} " =~ " vicuna " ]]
then
	echo "I started in vicuna $(date)"
	mkdir -p Alignment/VicunaContigs
	cd Alignment/VicunaContigs
	refContig=$pwd/VicunaOutput/contig.fasta
	pythContigs
	alignReads
	samtobam
	weeSamStats
	blast
	echo "I finished in vicuna $(date)"
	cd ../../
fi

#IVA
if [ $OPT_A = "iva" ] || [[ " ${array[@]} " =~ " iva "  ]]  
then
	echo "I started in iva $(date)"
	mkdir -p Alignment/IVAContigs
	cd Alignment/IVAContigs
	refContig=$pwd/IVAOutput/contigs/contigs.fasta
	pythContigs
	alignReads
	samtobam
	weeSamStats
	blast
	echo "I finished in iva $(date)"
	cd ../../
fi 

# Velvet
if [ $OPT_A = "velvet" ] || [[ " ${array[@]} " =~ " velvet "  ]]  
then
	echo "I started in velvet $(date)"
	mkdir -p Alignment/VelvetContigs
	cd Alignment/VelvetContigs
	refContig=$pwd/${pwd##*/}_VelvetAssembly/contigs.fa
	pythContigs
	alignReads
	samtobam
	weeSamStats
	blast
	echo "I finished in velvet $(date)"
	cd ../../
fi

