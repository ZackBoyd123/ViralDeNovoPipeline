#!/usr/bin/python3
import itertools
import sys


oldFile = sys.argv[1] 
discardFile = open(oldFile.split(".",1)[0]+"_discarded.fa", "w")



original = sys.stdout
sys.stdout = open(oldFile.split(".",1)[0]+"_edit.fa", "w")


with open(oldFile) as file:

    lst= []
    totalContigs = 0
    contigsOver300 = 0
    longestContig = ""
    longestLength = ""


    for l1,l2 in itertools.zip_longest(*[file]*2):

        
        if str(l1).startswith(">"):
            totalContigs += 1
        if not len(str(l2)) < 300:
            contigsOver300 += 1
            print(l1+l2, end='')
        else:
            discardFile.write(l1+l2)
        if len(l2) > len(str(longestLength)):
            longestContig = l1
            longestLength = str(l2)

    longestContig = longestContig.strip("\n")
sys.stdout = original
print("Longest contig in file:"+"\t"+"'"+longestContig+"'"+"\t"+"At length:"+"\t"+str(len(longestLength)))
print("Total Contigs in File:",oldFile+"\t",totalContigs)
print("Contigs removed from file as they weren't longer than specified threshold:","\t",str(int(totalContigs)-int(contigsOver300)))
print("New edited file name:"+"\t"+oldFile.split(".",1)[0]+"_edit.fasta")
sys.stdout = open(oldFile.split(".",1)[0]+"_stats.txt", "w")
print("Total Contigs in File:",oldFile+"\t",totalContigs)
print("Contigs removed from file as they weren't longer than specified threshold:","\t",str(int(totalContigs)-int(contigsOver300)))
print("Longest contig in file:","\t","'",longestContig,"'","\t","At length:","\t",str(len(longestLength)))
