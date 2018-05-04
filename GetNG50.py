#!/bin/python3

import sys
import csv


file = sys.argv[1]
genomeLength = sys.argv[2]
outfile = "UniqueContigsNG50.txt"


with open(file) as file:
    data = csv.reader(file,delimiter="\t")
    totalLength = []
    for line in data:
        totalLength.append(line[3])
        #print(line[0:4])



    file.close()
#total number of contigs mapped to genome
print(len(totalLength))

totalLength = [int(i)for i in totalLength]
totalLength = sorted(totalLength)
totalLength = reversed(totalLength)
totalLength = list(totalLength)
#print(totalLength)


total = 0
ng50Total = []


for i in totalLength:
    total += i

    if total > genomeLength / 2:
        ng50Total.append(i)
        #print(i)

sys.stdout = open(outfile,"w")
print("NG50 total:"+"\t"+str(ng50Total[0]))

