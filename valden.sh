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
    #Unpaired Reads
    U) OPT_U=$OPTARG
    ;;
    #Reference
    r) OPT_R=$OPTARG
    ;;
    esac
done 

# help options
if [ $# -eq 0 ];
then
	echo "+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Author: Zack Boyd																			|
| ----------------------------------------------------															|
|																					|
| This script will take an input, raw fastq file, trim it using trim galore, or trimmomatic.										|
|																					|
| Trimmed fastq files are then assembled using [spades] / [iva] /  [abyss] / [vicuna] / [idba] / [velvet] / [mira] / [masurca]						|
| If multiple assemblers are to be used, they must be presented in a comma seperate list with no whitespace, e.g iva,idba,velvet,vicuna					|
| 																					|
| The contig files generated from assembling the reads are automatically analysed using quast.										|
|	Note: If no reference file is specified the script will still run, but quast analysis will not occur.								|
| 																					| 
| 																					|
| 																					|
| Flags required for the script are as follows: 															|				
| [-1] : First raw fastq file of a paired reads sample.															|
| [-2] : Second raw fastq file of a paired reads sample.														|
| 																					|
| [-a] : Which aligner to use. Specified as [spades] / [iva] / [spadesnok] / [mira] / [masurca] / [abyss] / [vicuna] / [idba] / [velvet]				|					
|		Multiple aligners specified as one string seperated by commas and no whitespace.									|
|																					|
| [-r] : A reference genome specified from root \"/\" or relative to home \"~/\"											|
+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------+"
	exit 1
	

fi



#Get foldername
foldername=${PWD##*/}
pwd=`pwd`
#Make a bold font for warnings
bold=$(tput bold)
normal=$(tput sgr0)
#Make a substring variable to look for in OPT_A
substr=,
#Check if substr exists in OPT_A
if [ "${OPT_A/$substr}" = "$OPT_A" ]
then
	:
else
	printf "\t${bold} Multiple Aligners Selected...\n${normal}"
	IFS=',' read -a array <<< $OPT_A
fi

pathArray=()
quastName=()

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
		echo "Try [trimgalore] / [trimmomatic]"
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

#Check for vicuna config file.
if [ $OPT_A = "vicuna" ] && [ ! -e "vicuna_config.txt" ] 
then
	echo "No vicuna config file detected in current working directory"
	exit 1
else
	:	
fi

if [[ " ${array[@]} " =~ " vicuna " ]] && [ ! -e "vicuna_config.txt" ]
then
	echo "No vicuna config file in directory"
	exit 1
else
	:
fi
#Check for reference file given to command line
if [ -z $OPT_R ] 
then
	echo "No reference file given with [-r]"	
	exit 1
	printf "	${bold}!	Warning No reference specifed, script will still continue but QUAST analysis won't occur	!${normal}\n\n"

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
				:
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

#Aligner variables
if [ -n $OPT_A ] && [ -z $array ]
then
	if [ "$OPT_A" = "spades" ] || [ "$OPT_A" = "iva" ] || [ "$OPT_A" = "trinity" ] || [ "$OPT_A" = "allpaths" ] || [ "$OPT_A" = "celera" ] || [ "$OPT_A" = "abyss" ] || [ "$OPT_A" = "vicuna" ] || [ "$OPT_A" = "idba" ] || [ "$OPT_A" = "test" ] || [ "$OPT_A" = "velvet" ] || [ "$OPT_A" = "mira" ] || [ "$OPT_A" = "spadesnok" ] || [ "$OPT_A" = "masurca" ] 
	then
		:	
	else
		echo "The spelling of your aligner looks wrong"
		echo "Try [spades] / [iva] / [spadesnok] / [mira] / [] / [abyss] / [vicuna] / [idba] / [velvet] / [masurca]"
		exit 1
	fi
fi

##	Trimming Work	 ##
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

		RAW_1=$OPT_1
		RAW_2=$OPT_2
		OPT_1="${OPT_1%.fastq}_val_1.fq"
		OPT_2="${OPT_2%.fastq}_val_2.fq"
  
	else
		echo "Trimmed file already exists"
		RAW_1=$OPT_1
		RAW_2=$OPT_2
		OPT_1="${OPT_1%.fastq}_val_1.fq"
		OPT_2="${OPT_2%.fastq}_val_2.fq"

	 fi


else
	echo "unrecognised trim program selected: try [trimgalore] / [trimomatic]"

fi

# Function to print time to a file
scriptrun(){
	# Check how many lines are present in RunTime file
	var=$(cat ScriptRunTime.txt | wc -l)
	echo $(date) >> "ScriptRunTime.txt"

	# If more than two lines in file, rewrite file
	if [ "$var" -ge "2" ]
	then
		echo "Rewriting to run time logs, previous runs detected"
		> ScriptRunTime.txt
		echo $(date) >> "ScriptRunTime.txt"
	else
		:	
	fi
}
#skipAlignment="True"
if [[ $skipAlignment == "True" ]]
then
	:
else
	out=$(basename $OPT_R)
	out=${out%.fa*}
	# Get Genome stats here.
	mkdir -p GenomeStats
	GenomeStats.py -I $OPT_R > GenomeStats/$out".stats.txt"
	mdust $OPT_R > GenomeStats/$out".dust.fasta"
	repeat-match -n 50 GenomeStats/$out".dust.fasta" > GenomeStats/$out".repeat.fasta"
	exact-tandems_edited $OPT_R 50 > GenomeStats/$out".tandems.fasta"
	cd GenomeStats
	GenomeStats.py -I $out".dust.fasta" > $out".genomestats_afterdust.txt"
	cd ../
	########################

	# Get alignment stats here
	if [ ! -f Alignment/${PWD##*/}.bam ]
	then
		mkdir -p Alignment
		bowtie2 -p 15 -x BowtieIndexes/${PWD##*/} -1 $OPT_1 -2 $OPT_2 -S Alignment/${PWD##*/}".sam"
		echo "Sorting sam file: Alignment/'${PWD##*/}'.sam"
		samtools view -bS Alignment/${PWD##*/}".sam" | samtools sort -o Alignment/${PWD##*/}".bam"
		rm -f Alignment/${PWD##*/}.sam
		samtools index Alignment/${PWD##*/}".bam" 
		samtools depth -d 1000000 Alignment/${PWD##*/}".bam" > Alignment/${PWD##*/}"_coverage.bam"
		weeSAMv1.3 -b Alignment/${PWD##*/}".bam" -out Alignment/${PWD##*/}"_weesam.stats"
	else
		:
	fi
fi
########################


##	Aligner work	##
echo "Aligner: $OPT_A selected, running...."
echo "Script run time being logged, check 'ScriptRunTime.txt' in output directories"

# Add in all the methods which help with validating contigs.

pythContigs(){
	echo "Sorting Contigs!"
	Contigs2OneLine.py $refContig
	echo ${refContig%.fa*}
	file_to_link=$(echo $(basename $refContig) | cut -d. -f1)
	echo $file_to_link
	ln -s $(dirname $refContig)"/"$file_to_link"_oneline.fa" ./
	RemoveContigs.py *.fa*
}

alignReads(){
	echo "Aligning reads"
	mkdir -p BowtieIndexes
	mkdir -p BowtieOutput
	bowtie2-build *_oneline.fa BowtieIndexes/${PWD##*/}
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
	blastn -query *_oneline.fa -evalue 1e-5 -num_alignments 1 -num_threads 12 -db ../BlastDB/"$newblast" -out BlastOutput/${PWD##*/}_blast2ref.txt -outfmt "6 qseqid sseqid length qlen qstart qend sstart send mismatch gapopen evalue bitscore"
	#sort -k1,1 BlastOutput/${PWD##*/}_blastTemp.txt > BlastOutput/${PWD##*/}_blast2ref.txt
	sort -u -k1,1 BlastOutput/${PWD##*/}_blast2ref.txt > BlastOutput/${PWD##*/}_blast2ref_unique.txt
	#refLen=$(cat ../${pwd##*/}_weesam.stats | grep -v RefLength | awk '{sum += $2} END {print sum}')
	blastn -query *_oneline.fa -evalue 1e-5 -num_alignments 1 -num_threads 12 -db nt -out BlastOutput/${PWD##*/}_blastHits.txt -outfmt 6
	sort -u -k1,1 BlastOutput/${PWD##*/}_blastHits.txt > BlastOutput/${PWD##*/}_blastHits_unique.txt
	echo "I finished blast stuff at $(date)"
}

makecsv(){
	mkdir -p ../CSV_Files
	if [ ! -e ../CSV_Files/All_data.csv ]
	then
		echo ",Total#Contigs,Longest Contig,Contigs <300bp,Total#Reads,Total Mapped,Mean # Mapped,Median # Mapped,<5 Mapped,0 Mapped,Assembly Length,N50,NG50,Contigs blast 2 ref,Blast 1e-5 Hits" > ../CSV_Files/All_data.csv
	else
		:
	fi

	ConvertToTSV.py -1 BowtieOutput/WeeSamContigs_Stats.txt -2 *_oneline_stats.txt -3 BlastOutput/${PWD##*/}_blast2ref_unique.txt -4 BlastOutput/${PWD##*/}_blastHits_unique.txt
        tail -1 *.csv >> ../CSV_Files/All_data.csv
	GenomeStats.py -I $OPT_R --repeats ../../GenomeStats/$(basename $OPT_R)".repeat.fasta" --coverage ../${pwd##*/}_coverage.bam --contig BlastOutput/${PWD##*/}_blast2ref.txt --tandems ../../GenomeStats/$(basename $OPT_R)".tandem.fasta"
        mv *.csv ../CSV_Files

}

reaper(){
	echo "marking duplicates"
	newFold=${PWD##*/}
	java -jar ~orto01r/programs/picard.jar MarkDuplicates I=BowtieOutput/${PWD##*/}_aligned.bam O=BowtieOutput/${PWD##*/}_dupes.bam M=BowtieOutput/dupe_stats.txt
	echo "doing reapr"
	mkdir -p REAPROutput
	cd REAPROutput
	echo "checking fasta..."
	reapr facheck ../*_oneline.fa new_assembly
	reapr smaltmap -y 0.9 -s 2 -k 13 new_assembly.fa ../../../$OPT_1 ../../../$OPT_2 ../BowtieOutput/$newFold"_dupes.bam"
	reapr pipeline new_assembly.fa ../BowtieOutput/$newFold"_dupes.bam" outdir 
	cd ../
}

contents(){
	pythContigs
        alignReads
        samtobam
        weeSamStats
        blast
        reaper
#	Removing this for now as it isn't working 12/7/18
#	makecsv
}
##	Spades	##
if [ $OPT_A = "spades" ] || [[ " ${array[@]} " =~ " spades "  ]]
then
	mkdir -p SpadesOutput
	cd SpadesOutput	
	scriptrun	
	cd ../
	mkdir -p SpadesOutput
	(time spades.py --careful --cov-cutoff auto -1 $OPT_1 -2 $OPT_2 -k 77,99,127 -o SpadesOutput/) 2>&1 | tee spades.log.txt
	cd SpadesOutput
	scriptrun
	cd ../
	refContig=$pwd/SpadesOutput/contigs.fasta
	mkdir -p Alignment/SpadesContigs
	cd Alignment/SpadesContigs
	contents
	echo "I finished in spades at $(date)"
	cd ../../

	quastpath=SpadesOutput/contigs.fasta
	pathArray+=("$quastpath")
	quastName+=("Spades")
fi
## Spades No K
if [ $OPT_A = "spadesnok" ] || [[ " ${array[@]} " =~ " spadesnok " ]]
then
	mkdir -p SpadesOutputNoKMER
        (time spades.py -1 $OPT_1 -2 $OPT_2 -o SpadesOutputNoKMER/) 2>&1 | tee spadesnok.log.txt
	refContig=$pwd/SpadesOutputNoKMER/contigs.fasta
	mkdir -p Alignment/SpadesOutputNoKMER
	cd Alignment/SpadesOutputNoKMER
	contents
	echo "Finished in spades at $(date)"
	quastpath=$refContig
	pathArray+=("$quastpath")
	quastName+=("SpadesNoK")
fi

##	 IDBA	  ##
if [ $OPT_A = "idba" ] || [[ " ${array[@]} " =~ " idba " ]]
then
	echo "${bold} Begining IDBA analysis....${normal}"
	mkdir -p IDBAOutput
	cd IDBAOutput
	scriptrun
	cd ../
	fq2fa --merge "$OPT_1" "$OPT_2" "idba.merged.fa"
	(time idba -r "idba.merged.fa" -o IDBAOutput/) 2>&1 | tee idba.log.txt
	cd IDBAOutput
	scriptrun 
	cd ../
	
	refContig=$pwd/IDBAOutput/contig.fa
	mkdir -p Alignment/IDBAContigs
        cd Alignment/IDBAContigs
	contents
        echo "I finished in IDBA at $(date)"

        cd ../../


	quastpath=IDBAOutput/contig.fa
	pathArray+=("$quastpath")
	quastName+=("IDBA")
fi
#	Abyss
if [ $OPT_A = "abyss" ] || [[ " ${array[@]} " =~ " abyss " ]]
then
	mkdir -p ABySSOutput
	cd ABySSOutput
	scriptrun
	(time abyss-pe name="$foldername" k=96 in="../$OPT_1 ../$OPT_2") 2>&1 | tee abyss.log.txt
	scriptrun
	cd ../

	refContig=$pwd/ABySSOutput/${pwd##*/}-contigs.fa
	mkdir -p Alignment/ABySSContigs
        cd Alignment/ABySSContigs
	contents
        echo "I finished in ABySS at $(date)"

        cd ../../
	

	quastpath=ABySSOutput/"${PWD##*/}-contigs.fa"
	pathArray+=("$quastpath")
	quastName+=("ABySS")
fi

# VICUNA
if [ $OPT_A = "vicuna" ] || [[ " ${array[@]} " =~ " vicuna " ]]
then
	mkdir -p VicunaOutput
	currentpath=$pwd"/vicuna_config.txt"
	otherpath=$pwd

	#Export vars for python script
	export currentpath
	export otherpath
	cd VicunaOutput 

	#Sym link to trimmed fastq files in directory above
	ln -s ../$OPT_1 ./
	ln -s ../$OPT_2 ./
	
	#Edit the config file, output: vicuna_edit_config.txt
	EditVicunaConfig.py

	#Run assembler
	scriptrun
	(time OMP_NUM_THREADS=12 ~hugh01j/bin/VICUNA_v1.3/executable/vicuna-omp.static.linux64 vicuna_config.txt) 2>&1 | tee vicuna.log.txt
	scriptrun
	cd ../

	refContig=$pwd/VicunaOutput/contig.fasta
	mkdir -p Alignment/VicunaContigs
        cd Alignment/VicunaContigs
	contents
        echo "I finished in Vicuna at $(date)"
        cd ../../


	quastpath=VicunaOutput/contig.fasta
	pathArray+=("$quastpath")
	quastName+=("Vicuna")
fi

#IVA
if [ $OPT_A = "iva" ] || [[ " ${array[@]} " =~ " iva "  ]] 
then
	mkdir -p IVAOutput
	cd IVAOutput
	scriptrun
	(time iva -f ../$OPT_1 -r ../$OPT_2 contigs) 2>&1 | tee iva.log.txt
	scriptrun
	cd ../

	refContig=$pwd/IVAOutput/contigs/contigs.fasta
	mkdir -p Alignment/IVAContigs
        cd Alignment/IVAContigs
	contents
        cd ../../


	echo "Finished in iva $(date)"
	quastpath=IVAOutput/contigs/contigs.fasta
	pathArray+=("$quastpath")
	quastName+=("IVA")

fi

# Velvet
if [ $OPT_A = "velvet" ] || [[ " ${array[@]} " =~ " velvet "  ]]
then
	echo "Running Velvet..."
	(time -p sh -c 'velveth ${PWD##*/}"_VelvetAssembly" 99 -shortPaired -fastq -separate '$OPT_1' '$OPT_2' ; velvetg ${PWD##*/}"_VelvetAssembly"') 2>&1 | tee velvet.log.txt

	refContig=$pwd/${pwd##*/}_VelvetAssembly/contigs.fa
	mkdir -p Alignment/VelvetContigs
        cd Alignment/VelvetContigs
	contents
        echo "I finished in spades at $(date)"
        cd ../../

	quastpath=${PWD##*/}_VelvetAssembly/contigs.fa
	pathArray+=("$quastpath")
	quastName+=("Velvet")

fi

if [ $OPT_A = "mira" ] || [[ " ${array[@]} " =~ " mira " ]]
then
	ln -s $OPT_1 ${OPT_1%.fq}.fastq
	ln -s $OPT_2 ${OPT_2%.fq}.fastq

	printf "project = MIRA\njob = est,denovo,accurate\nparameters = -GE:not=15 -NW:cmrnl=warn SOLEXA_SETTINGS\n\nreadgroup = testgrp\ndata = ${OPT_1%.fq}.fastq ${OPT_2%.fq}.fastq\ntechnology = solexa\n" > mira.manifest 
	(time mira mira.manifest) 2>&1 | tee mira.log.txt
	cd MIRA_assembly/MIRA_d_results
	mv MIRA_out.unpadded.fasta contig.fa
	refContig=$pwd/MIRA_assembly/MIRA_d_results/contig.fa
	cd ../../
	mkdir -p Alignment/MiraContigs
	cd Alignment/MiraContigs
	contents 

	echo "Finised mira at: $(date)"
	cd ../../

	quastpath=MIRA_assembly/MIRA_d_results/contig.fa
	pathArray+=("$quastpath")
	quastName+=("Mira")

fi

if [ $OPT_A = "masurca" ] || [[ "${array[@]} " =~ " masurca " ]]
then
	mkdir -p MasurcaOutput
	cd MasurcaOutput
	printf "DATA\nPE= pe 150 20 $pwd/$RAW_1 $pwd/$RAW_2 \nEND" > masurca.config
	(time /home1/boyd01z/SoftwareInstall/MaSuRCA-3.2.6/bin/masurca masurca.config) 
	(time ./assemble.sh) 2>&1 | tee masurca.log.txt 
	cd ../
	mkdir -p Alignment/MasurcaContigs/
	cd Alignment/MasurcaContigs/
	refContig=$pwd/MasurcaOutput/CA/final.genome.scf.fasta
	contents

	echo "Finished in Masurca at $(date)"
	cd ../../
	quastpath=$refContig
	pathArray+=("$quastpath")
	quastName+=("Masurca")

fi

if [ -z $OPT_R ]
then	
	echo "Programme terminating"
else
	echo "Im here, doing quast analysis"
	
	if [ -z $array ]
	then
		echo "${bold}Generating quast stats for one aligner specified${normal}"
		mkdir -p "QuastOutput_$OPT_A"
		quast.py -o "QuastOutput_$OPT_A" -R $OPT_R $quastpath

	else
		var=""
		printf -v var "%s," "${quastName[@]}"
		var=${var%, }
		var=$(echo $var)
		var=${var%,}
		list=$(echo ${pathArray[@]})

		echo "Quast contigs ${pathArray[@]}"
		echo "${bold}Generating quast stats for $var ....${normal}"
		echo "quast.py -l $var $list -R $OPT_R --reads1 $OPT_1 --reads2 $OPT_2"
		quast.py -l $var $list -R $OPT_R --reads1 $OPT_1 --reads2 $OPT_2
	fi
fi
