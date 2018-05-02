#!/bin/bash

while getopts :a:1:2:t:U:r: TEST; do
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
    #Trimming program
    t) OPT_T=$OPTARG
    ;;
    #Reference
    r) OPT_R=$OPTARG
    ;;
    esac
done 

foldername=${PWD##*/}
pwd=`pwd`
#Make a bold font for warnings
bold=$(tput bold)
normal=$(tput sgr0)

##	Check variables are in place before running script	##
#Trimming variables
if [ -z $OPT_T ] 
then
	echo "no trimming programme selected with [-t]"
	exit 1
fi

if [ -n $OPT_T ] 
then
	if [ "$OPT_T" = "trimgalore" ] || [ "$OPT_T" = "trimmomatic" ] 
	then
		:	
	else
		echo "The spelling of your selected trimming programme is wrong"
		echo "Try [trimgalore]"
		exit 1
	fi
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

#Check for reference file given to command line
if [ -z $OPT_R ] 
then
	echo "No reference file given with [-r]"	
	exit 1
else
	#Get first char of reference $OPT_R to check for full path
	firstchar="${OPT_R:0:1}"
	
	if [ "$firstchar" == "/" ] || [ "$firstchar" == "~" ]
	then
		
		#If file exists and bowtie indexes dont build indexes, else exit script
        	if [ -e $OPT_R ]
		then
			if [ -d BowtieIndexes ] 
			then
				echo "Bowtie2 Indexes already present within directory"
			else
				printf "Reference file detected...\n ${bold} No Bowtie2 indexes present, building now ${normal}\n\n"
				mkdir -p BowtieIndexes
				bowtie2-build "$OPT_R" BowtieIndexes/${PWD##*/}	
			fi
		else
			echo "Couldn't find specifed reference file in given path"
			exit 1
		fi

	else
        	echo "Please specify reference [-r] as a full "/" or relative "~" path"
	        exit 1
	fi

	
	
fi

if [ $OPT_T = "trimgalore" ] 
then
	if [ ! -e "${OPT_1%.fastq}_val_1.fq" ] && [ ! -e "${OPT_2%.fastq}_val_2.fq" ]
	then
		echo "Running selected program $OPT_T"
		echo "Paired reads detected, trimming"
		trim_galore --paired "$OPT_1" "$OPT_2" --length 80

		mkdir -p ReportsAndUntrimmed
		for i in $(ls | grep trimming)
		do
			mv $i ReportsAndUntrimmed      
		done

		OPT_1="${OPT_1%.fastq}_val_1.fq"
		OPT_2="${OPT_2%.fastq}_val_2.fq"
  
	else
		echo "Trimmed file already exists"
		OPT_1="${OPT_1%.fastq}_val_1.fq"
		OPT_2="${OPT_2%.fastq}_val_2.fq"

	 fi

fi


mkdir -p Alignment
bowtie2 -p 15 -x BowtieIndexes/${PWD##*/} -1 $OPT_1 -2 $OPT_2 -S Alignment/${PWD##*/}".sam"
echo "Sorting sam file: Alignment/'${PWD##*/}'.sam"
samtools view -bS Alignment/${PWD##*/}".sam" | samtools sort -o Alignment/${PWD##*/}".bam"
rm -f Alignment/${PWD##*/}.sam

samtools index Alignment/${PWD##*/}".bam" 

weeSAMv1.3 -b Alignment/${PWD##*/}".bam" -out Alignment/${PWD##*/}"_weesam.stats"
