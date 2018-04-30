GenomeStats.py $OPT_R
mdust $OPT_R > ${PWD##*/}".dust.fasta"
#python3 /home1/boyd01z/MScProjectWork/RunScripts/FixFasta.py ${PWD##*/}".dust.fasta"

printf "!!!Stats After Dust!!!\n"

GenomeStats.py ${PWD##*/}".dust.fasta"
#rm -f ${PWD##*/}.dust.fasta
