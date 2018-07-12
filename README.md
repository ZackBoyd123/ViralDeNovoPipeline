### _valden_  
A bash / python pipeline used to benchmark multiple NGS assembler programmes.  
  
#### Run the script as follows:  
``valden.sh -a assembler -t trimgalore -1 sample_R1.fq -2 sample_R2.fq -r reference.fa``  
  
The above example demonstrates one assembler being benchmarked, however, valden has the capability to benchmark as many as you would like, as long as the code has been added to the pipeline. The pipeline currently supports the following assemblers:  
`` spades, vicuna, mira, masurca, iva, idba, abyss and velvet``  
#### Example benchmarking all the assemblers.  
`` valden.sh -a spades,vicuna,mira,masurca,iva,idba,abyss,velvet -1 sample_R1.fq -2 sample_R2.fq -r reference.fa -t trimgalore``  

