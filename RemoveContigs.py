#!/usr/bin/python3
import itertools
import sys

# This script is only going to work if contigs are on one line each.  I always use this script in the pipeline after ive converted
# the contig fasta file. 
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
	# If line no. one starts with a > increment total contigs by one.         
        if str(l1).startswith(">"):
            totalContigs += 1
	# If contig is above 300 write that sequence plus its header to a new file.
        if not len(str(l2)) < 300:
            contigsOver300 += 1
            print(l1+l2, end='')
	# If the contig is < 300 write it to a different file.
        else:
            discardFile.write(l1+l2)
        if len(l2) > len(str(longestLength)):
            longestContig = l1
            longestLength = str(l2)

    longestContig = longestContig.replace("\n","")
# Print some stats to an output file. 
sys.stdout = original
print("Longest contig in file:"+"\t"+"'"+longestContig+"'"+"\t"+"At length:"+"\t"+str(len(longestLength)))
print("Total Contigs in File:",oldFile+"\t",totalContigs)
print("Contigs removed from file as they weren't longer than specified threshold:","\t",str(int(totalContigs)-int(contigsOver300)))
print("New edited file name:"+"\t"+oldFile.split(".",1)[0]+"_edit.fasta")
sys.stdout = open(oldFile.split(".",1)[0]+"_stats.txt", "w")
print("Total Contigs in File:",oldFile+"\t",totalContigs)
print("Contigs removed from file as they weren't longer than specified threshold:","\t",str(int(totalContigs)-int(contigsOver300)))
print("Longest contig in file:","\t","'",longestContig,"'","\t","At length:","\t",str(len(longestLength)))
