### _valden_  
A bash / python pipeline used to benchmark multiple NGS assembler programmes.  
  
#### Run the script as follows:  
``valden.sh -a assembler -t trimgalore -1 sample_R1.fq -2 sample_R2.fq -r reference.fa``  
  
The above example demonstrates one assembler being benchmarked, however, valden has the capability to benchmark as many as you would like, as long as the code has been added to the pipeline. The pipeline currently supports the following assemblers:  
`` spades, vicuna, mira, masurca, iva, idba, abyss and velvet``  
#### Example benchmarking all the assemblers.  
`` valden.sh -a spades,vicuna,mira,masurca,iva,idba,abyss,velvet -1 sample_R1.fq -2 sample_R2.fq -r reference.fa -t trimgalore``    
  
### Scripts  
What each script does and how it is run.  
#### _GenomeStats.py_  
``user@host: GenomeStats.py -I input.fasta --coverage depth.txt --repeats repeats.txt --tandems tandems.txt --multi-file    --contig blast.txt ``  
GenomeStats.py generates two output files containing information required by ``Final_contig_plot.r`` to generate its output pdf file. These files will be called:  
``input_stats.csv & input_multi_stats.csv``  
``input_stats.csv`` contains the plotting information for the genome plot and the coverage plot and ``input_multi_stats`` contains the information on where the contigs blast to the reference.  
The command line options are as follows:  
``-I:`` An input fasta file.   
``--coverage:`` A file containing information on the read coverage, usually generated by samtools depth.  
``--repeats:`` A file telling the script where repetitive regions in the genome occur, generated by ``repeat-match``.    
``--tandems:`` Where tandem repeats occur in the genome, generated by exact-tandems.    
``--multi-file:`` If this is specified the script will treat all .txt files in the directory as files containing coordinates where contigs blast to a reference. Only specify this if you are sure each .txt file was produced by blastn.    
``--contig:`` Blast file containing coordinates of where contigs blast to a reference sequence.    
The only required flag is ``-I`` if this is specified the script will run and produces information on the specified fasta file.  
#### _Final_contig_plot.r_  
``user@host: ContigPlot.r file1 file2 out.pdf``  
Generates a genome plot, coverage plot and a plot of where contigs blast to a reference. file1 and file2 are generated by GenomeStats.py and out.pdf is the name of the output R plot.  
#### _Contigs2OneLine.py_  
``user@host: Contigs2OneLine.py input output``  
This script will convert all multiline sequences in a file to one line.  Requires two arguments, an input file and an output file name.  
#### _RemoveContigs.py_  
``user@host: RemoveContigs.py input output``  
Takes an input fasta file containing multiple sequences and prints stats of that file. The script will also write two new files, one with all contigs less than 300bp and one with all contigs greater than 300bp.  
The statistics generated are:  
```
Longest contig in file.  
Number of contigs in file.  
Number of contigs < 300bp in file.  
```
#### _WeeSamContigs.py_
``user@host: input genome_length[int]``  
Analyses the output from weeSAM and converts it to a file of relevant information.  
The information in the ouput file is:  
```
Total # Reads mapped,   # Reads mapped to longest contig,   Mean # reads mapped,    Median # reads mapped  
Total assembly length,  N50,                                Genome Length,          NG50,   LG100  
# Contigs < 5 reads mapped, Total # contigs in file   
```  
The input file is a file generated from weeSAM and the genome length is the length of genome in base pairs.  
#### _EditVicunaConfig.py_  
Edits specific lines in the vicuna config file to tell vicuna where fastq files are. The script requires specific exported bash variables, so it's best not to run this outside of the bash script.  
 




