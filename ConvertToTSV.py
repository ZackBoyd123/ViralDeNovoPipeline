#!/usr/bin/python3
import sys
import csv
import argparse
import re

parser = argparse.ArgumentParser(description="")
parser.add_argument("-1","--WeeSamFile",help="first input file",required=True)
parser.add_argument("-2","--ContigFile",help="first input file",required=True)
parser.add_argument("-3","--BlastFile",help="first input file",required=True)
parser.add_argument("-4","--BlastHits",help="first input file",required=True)
args = parser.parse_args()

blastFile = args.BlastFile
blastHits = args.BlastHits
weeSam = args.WeeSamFile
contigFile = args.ContigFile
suffix = blastFile.rsplit("/")[-1]
suffix=suffix.split("_")[0][:-7]
#print(blastFile,blastHits,contigFile,weeSam)

with open(weeSam) as f1:
    readsList = []
    mappedList = []
    assembly = []
    for i, line in enumerate(f1):
        if i == 1:
            readsList.append(line)
        if i == 4:
            assembly.append(line)
        if i == 7:
            mappedList.append(line)


        #print(line)
    readsList = "".join(readsList)
    readsList = readsList.strip("\n")
    readsList = readsList.split(",")


    mappedList = "".join(mappedList)
    mappedList = mappedList.strip("\n")
    mappedList = mappedList.split(",")
    zeromap = mappedList[2]
    zeromap = int(zeromap)
    mappedList = mappedList[0]
    #print(mappedList)

    #print(assembly)
    assembly = "".join(assembly)
    assembly = assembly.strip("\n")
    assembly = assembly.split(",")
    N50 = assembly[1]
    NG50 = assembly[3]
    assembly = assembly[0]




    f1.close()

with open(contigFile) as f2:
    contigList = []
    longList = []
    removed = []
    for i, line in enumerate(f2):
        if i == 0:
            contigList.append(line)
        if i == 1:
            removed.append(line)

        if i == 2:
            longList.append(line)

    contigList = "".join(contigList)
    contigList = contigList.strip("\n")
    contigList = contigList.split("\t")
    contigList = contigList[1]
    contigList = int(contigList)
    #print(contigList)


    ###
    longList = "".join(longList)
    longList = longList.strip("\n")
    #if "'" in longList:
    #    #print("Im in here")
    #    result = re.search("'(.*)'",longList)
    #    longContig = result.group(0)
    #else:
    #    longList = longList.split("\t")
    #    longContig = longList[3]
    #    longContig = longContig.strip(" ")


    longList = longList.split("\t")
    #print(longList)
    longContig = longList[-1]
    longContig = longContig.strip(" ")
    print(longContig)

    removed = "".join(removed)
    removed = removed.strip("\n")
    removed = removed.split("\t")
    removed = removed[1]
    removed = removed.strip(" ")
    #print(removed)
    f2.close()

with open(blastFile) as f3:
    unique = 0
    for line in f3:
        unique += 1
    f3.close()
with open(blastHits) as f4:
    total = 0
    for line in f4:
        total += 1


sys.stdout=open(suffix+"_data.csv","w")
#print("\t"+"Total # Contigs"+"\t"+"Longest Contig"+"\t"+"Contig < 300bp"+"\t"+"Total # Reads"+"\t"+"Total Mapped to longest contig"+"\t"+"Mean # mapped"+"\t"+"Median # mapped"+"\t"+"<5 mapped"+"\t"+"0 mapped"+"\t"+"Assembly Length"+"\t"+"N50"+"\t"+"NG50"+"\t"+"Contigs blast 2 ref"+"\t"+"Blast hits")
#print(suffix,end="\t")
#print(str(contigList)+"\t"+str(longContig)+"\t"+str(removed),end="\t")
#print(*readsList,sep="\t",end="\t")
#print(mappedList+"\t"+str(contigList-zeromap)+"\t"+str(assembly)+"\t"+str(N50)+"\t"+str(NG50)+"\t"+str(unique)+"\t"+str(total))

print(","+"Total # Contigs"+","+"Longest Contig"+","+"Contig < 300bp"+","+"Total # Reads"+","+"Total Mapped to longest contig"+","+"Mean # mapped"+","+"Median # mapped"+","+"<5 mapped"+","+"0 mapped"+","+"Assembly Length"+","+"N50"+","+"NG50"+","+"Contigs blast 2 ref"+","+"Blast hits")
print(suffix,end=",")
print(str(contigList)+","+str(longContig)+","+str(removed),end=",")
print(*readsList,sep=",",end=",")
print(mappedList+","+str(contigList-zeromap)+","+str(assembly)+","+str(N50)+","+str(NG50)+","+str(unique)+","+str(total))
