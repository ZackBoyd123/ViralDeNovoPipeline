#!/usr/bin/python3
import sys
import csv
import statistics

file = sys.argv[1]
genomeLength = int(sys.argv[2])


with open (file) as file:
    data = csv.reader(file, delimiter="\t")
    next(data,None)

    totalContigs = 0
    longestContig = 0
    readsmapped = []
    refLength = []

    for line in data:
        totalContigs += 1
        toAppend = line[2]
        readsmapped.append(toAppend)
        refAppend = line[1]
        refLength.append(refAppend)







file.close()



refLength = [int(i)for i in refLength]
#Where in the list the longest contig occurs
#print(refLength.index(max(refLength)))
#Max contig length
#print(max(refLength))
readsmapped = [int(i)for i in readsmapped]
x = readsmapped[refLength.index(max(refLength))]

#print("Reads mapped to longest contig:","\t",x)

sys.stdout = open("WeeSamContig_stats.txt","w")
readsmapped = sorted(readsmapped)


print("\t","Total Number of reads Mapped","\t","Reads Mapped to Longest Contig","\t","Mean number of reads mapped","\t","Median Number of Reads Mapped","\t","Total Number of Contigs in File")

print("\t",sum(readsmapped),"\t",x,"\t",statistics.mean(readsmapped),"\t",statistics.median(readsmapped),"\t",totalContigs)

print("\n","\t","Total Assembly length","\t","N50","\t","Genome Length","\t","NG50")

assemblyLength = sum(refLength)


y=(sorted(refLength))
halfAssembly = (assemblyLength / 2)

new = list(reversed(y))

total = 0
ng50total = 0
removeList=[]

n50list = []
ng50list = []

for i in new:
    total += i

    if total > halfAssembly:
        n50list.append(i)

for i in new:
    ng50total += i

    if ng50total > genomeLength / 2:
        ng50list.append(i)




print("\t",assemblyLength,"\t",n50list[0],"\t",genomeLength,"\t",ng50list[0])



