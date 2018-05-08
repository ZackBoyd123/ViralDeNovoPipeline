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
| Trimmed fastq files are then assembled using [spades] / [iva] / [trinity] / [allpath] / [celera] / [abyss] / [vicuna] / [idba] / [velvet]				|
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
| [-a] : Which aligner to use. Specified as [spades] / [iva] / [trinity] / [allpath] / [celera] / [abyss] / [vicuna] / [idba] / [velvet]				|					
|		Multiple aligners specified as one string seperated by commas and no whitespace.									|
|																					|
| Optional flags for the script are as follows:																|
| [-r] : A reference genome specified from root \"/\" or relative to home \"~/\"												|
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
	if [ "$OPT_A" = "spades" ] || [ "$OPT_A" = "iva" ] || [ "$OPT_A" = "trinity" ] || [ "$OPT_A" = "allpaths" ] || [ "$OPT_A" = "celera" ] || [ "$OPT_A" = "abyss" ] || [ "$OPT_A" = "vicuna" ] || [ "$OPT_A" = "idba" ] || [ "$OPT_A" = "test" ] || [ "$OPT_A" = "velvet" ]
	then
		:	
	else
		echo "The spelling of your aligner looks wrong"
		echo "Try [spades] / [iva] / [trinity] / [allpath] / [celera] / [abyss] / [vicuna] / [idba] / [velvet]"
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

		OPT_1="${OPT_1%.fastq}_val_1.fq"
		OPT_2="${OPT_2%.fastq}_val_2.fq"
  
	else
		echo "Trimmed file already exists"
		OPT_1="${OPT_1%.fastq}_val_1.fq"
		OPT_2="${OPT_2%.fastq}_val_2.fq"

	 fi


#elif [ $OPT_T = "trimmomatic" ]
#then
#	if [ ! -e "$foldername""_trimmomatic_1P.fq" ] && [ ! -e "$foldername""_trimmomatic_2P.fq" ]
#	then
#		# Got default settings from:http://www.usadellab.org/cms/?page=trimmomatic "
#		echo "Running selected program $OPT_T"
#		java -jar /home1/boyd01z/SoftwareInstall/Trinityrnaseq-v2.6.6/trinity-plugins/Trimmomatic-0.36/trimmomatic.jar PE -threads 12 $OPT_1 $OPT_2 -baseout "$foldername""_trimmomatic.fq" ILLUMINACLIP:/home1/boyd01z/SoftwareInstall/Trinityrnaseq-v2.6.6/trinity-plugins/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
#		OPT_1=$foldername"_trimmomatic_1P.fq"
#		OPT_2=$foldername"_trimmomatic_2P.fq"
#
#	else
#		echo "Trimmed files already generated by valden detected in directory, running....."
#		OPT_1=$foldername"_trimmomatic_1P.fq"
#                OPT_2=$foldername"_trimmomatic_2P.fq"		
#
#	fi
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

# Get Genome stats here.
mkdir -p GenomeStats
GenomeStats.py $OPT_R > GenomeStats/genomestats.txt
mdust $OPT_R > GenomeStats/${PWD##*/}".dust.fasta"
#python3 /home1/boyd01z/MScProjectWork/RunScripts/FixFasta.py ${PWD##*/}".dust.fasta"

printf "!!!Stats After Dust!!!\n"
cd GenomeStats
GenomeStats.py ${pwd##*/}".dust.fasta" > genomestats_afterdust.txt
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
	weeSAMv1.3 -b Alignment/${PWD##*/}".bam" -out Alignment/${PWD##*/}"_weesam.stats"
else
	:
fi
########################


##	Aligner work	##
echo "Aligner: $OPT_A selected, running...."
echo "Script run time being logged, check 'ScriptRunTime.txt' in output directories"
#echo $OPT_A
#printf '%s\n' "${array[@]}"
#Testing stuff
if [ $OPT_A = "test" ] || [[ " ${array[@]} " =~ " test " ]]
then
	testvar="IAMTESTVAR"
	#var=""
	echo "test in array"
	if [[ " ${array[@]} " =~ " cat " ]]
	then
		echo "test and idba"
		pathArray+=("IDBAOutput/out/contig.fa")
		quastName+=("IDBA")
		echo "${pathArray[0]}"
		echo "\"$OPT_A\" in array"
	fi
	if [[ " ${array[@]} " =~ " dog " ]]
	then
		echo "test and idba and iva"
		pathArray+=("IVAOutput/contigs.fa")
		quastName+=("IVA")
	#	printf -v var "%s," "${quastName[@]}"
	fi
	
	#var=${var%,}
	#echo "New quast thing $var"
	#echo "Quast contigs ${pathArray[@]}"
	#exit 1
else
	echo "no test in array"
fi

# Add in all the methods which help with validating contigs.

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
	ConvertToTSV.py -1 BowtieOutput/WeeSamContigs_Stats.txt -2 *contig*_oneline_stats.txt -3 BlastOutput/${PWD##*/}_blast2ref_unique.txt -4 BlastOutput/${PWD##*/}_blastHits_unique.txt
	tail -1 *.csv >> ../CSV_Files/All_data.csv
	mv *.csv ../CSV_Files
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
}
##	Spades	##
if [ $OPT_A = "spades" ] || [[ " ${array[@]} " =~ " spades "  ]]
then
	mkdir -p SpadesOutput
	cd SpadesOutput	
	scriptrun	
	cd ../
	(time python2 ~fawc01h/Documents/SPAdes-3.6.1-Linux/bin/spades.py --careful --cov-cutoff auto -1 $OPT_1 -2 $OPT_2 -k 77,99,127 -o SpadesOutput/) 2>&1 | tee spades.log.txt
	cd SpadesOutput
	scriptrun
	cd ../
	refContig=$pwd/SpadesOutput/contigs.fasta
	mkdir -p Alignment/SpadesContigs
	cd Alignment/SpadesContigs
	pythContigs
	alignReads
	samtobam
	weeSamStats
	makecsv
	blast
	echo "I finished in spades at $(date)"
	cd ../../

	quastpath=SpadesOutput/contigs.fa
	pathArray+=("$quastpath")
	quastName+=("Spades")
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
        pythContigs
        alignReads
        samtobam
        weeSamStats
	makecsv
        blast
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
        pythContigs
        alignReads
        samtobam
        weeSamStats
	makecsv
        blast
        echo "I finished in ABySS at $(date)"
        cd ../../


	quastpath=ABySSOutput/"${PWD##*/}-contigs.fa"
	pathArray+=("$quastpath")
	quastName+=("ABySS")
fi

# Celera
#if [ $OPT_A = "celera" ] || [[ " ${array[@]} " =~ " celera " ]] 
#then
#	mkdir -p CeleraOutput
#	frgname=${PWD##*/}
#	cd CeleraOutput
#	scriptrun
#	
#	# Generate FRG files needed for celera assembly
#	fastqToCA -libraryname celera -insertsize 500 50 -technology illumina -mates "../$OPT_1","../$OPT_2" > "$frgname"".frg"
#	echo "FRG File generated"	
#
#	#Begin alignment
#	echo "Running celera aligner...."
#	time runCA -p "$frgnname" -d "CeleraRun" "$frgname"."frg"
#
#	scriptrun
#	cd ../
#fi

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
	(time ~hugh01j/bin/VICUNA_v1.3/executable/vicuna-omp.static.linux64 vicuna_config.txt) 2>&1 | tee vicuna.log.txt
	scriptrun
	cd ../

	refContig=$pwd/VicunaOutput/contig.fasta
	mkdir -p Alignment/VicunaContigs
        cd Alignment/VicunaContigs
        pythContigs
        alignReads
        samtobam
        weeSamStats	
	makecsv
        blast
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
        pythContigs
        alignReads
        samtobam
        weeSamStats
	makecsv
        blast
        cd ../../


	echo "Finished in iva $(date)"
	quastpath=IVAOutput/contigs/contigs.fa
	pathArray+=("$quastpath")
	quastName+=("IVA")

fi
#Allpaths-lg
#if [ $OPT_A = "allpaths" ] || [[ " ${array[@]} " =~ " allpaths "  ]] 
#then	
#	#Make directories which allpaths needs
#	mkdir -p AllPathsOutput/${PWD##*/}/Data
#
#	#Generate in_groups.csv and in_libs.csv
#	printf "group_name, library_name, file_name \n1, Illumina_1, $pwd"/"*.fastq \n2, Illumina_2, $pwd"/"*.fastq" > "in_groups.csv"
#	printf "library_name, project_name, organism_name, type,paired, frag_size, frag_stddev, insert_size, insert_stddev, read_orientation, genomic_start, genomic_end \nIllumina_1, AllPaths, ${PWD##*/}, fragment, 1, 377, 137, , , inward, , \nIllumina_2, AllPaths, ${PWD##*/}, jumping, 1, , , 143, 15, outward, , " > "in_libs.csv"
#
#	#Move them to AllPathsOutput
#	mv in_groups.csv AllPathsOutput
#	mv in_libs.csv AllPathsOutput
#	cd AllPathsOutput
#
#	#Get name of directory one above, which is the name of the folder created in AllPathsOutput
#	allpathsname=`awk -F "/" '{print $(NF-1)}'<<< $PWD`
#	
#	#Prepare files for allpaths run
#        PrepareAllPathsInputs.pl DATA_DIR=$PWD""/"$allpathsname" PICARD_TOOLS_DIR=/usr/bin/picard-tools PLOIDY=1
#	cd $allpathsname
#	mv frag_reads_orig.* jump_reads_orig.* ploidy Data/
#	cd ../
#
#	#Run Allpaths
#	RunAllPathsLG PRE=./ REFERENCE_NAME="$allpathsname" DATA_SUBDIR=Data RUN=run SUBDIR=small OVERWRITE=True VAPI_WARN_ONLY=True
#	cd ../
#	
#	quastpath=AllPathsOutput/${PWD##*/}/Data/run/ASSEMBLIES/small/final.contigs.fasta
#	pathArray+=$quastpath

#fi

# Velvet
if [ $OPT_A = "velvet" ] || [[ " ${array[@]} " =~ " velvet "  ]]
then
	echo "Running Velvet..."
#	echo $OPT_1
#	echo $OPT_2

	(time -p sh -c 'velveth ${PWD##*/}"_VelvetAssembly" 99 -shortPaired -fastq -separate '$OPT_1' '$OPT_2' ; velvetg ${PWD##*/}"_VelvetAssembly"') 2>&1 | tee velvet.log.txt

	refContig=$pwd/${pwd##*/}_VelvetAssembly/contigs.fa
	mkdir -p Alignment/VelvetContigs
        cd Alignment/VelvetContigs
        pythContigs
        alignReads
        samtobam
        weeSamStats
	makecsv
        blast
        echo "I finished in spades at $(date)"
        cd ../../



	quastpath=${PWD##*/}_VelvetAssembly/contigs.fa
	pathArray=("$quastpath")
	quastName=("Velvet")

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
		printf -v var "%s, " "${quastName[@]}"
		var=${var%, }
		var=$(echo \"$var\")
		list=$(echo ${pathArray[@]})

		echo "Quast contigs ${pathArray[@]}"
		#Finish this bit vvv
		echo "${bold}Generating quast stats for $var ....${normal}"
		echo "quast.py -l $var $list -R $OPT_R --reads1 $OPT_1 --reads2 $OPT_2"
		quast.py -l $var $list -R $OPT_R --reads1 $OPT_1 --reads2 $OPT_2
	fi
fi
