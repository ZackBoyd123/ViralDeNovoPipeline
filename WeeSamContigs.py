#!/usr/bin/python3
import sys
import csv
import statistics

file = sys.argv[1] 
genomeLength = sys.argv[2]


with open (file) as file:
    data = csv.reader(file, delimiter="\t")
    next(data,None)

    totalContigs = 0
    longestContig = 0
    readsmapped = []
    refLength = []
    contigbelow5 = 0
    contigzero = 0

    for line in data:
        totalContigs += 1
        toAppend = line[2]
        readsmapped.append(toAppend)
        refAppend = line[1]
        refLength.append(refAppend)
        if int(toAppend) <= 5:
            contigbelow5 += 1
        elif int(toAppend) == 0:
            contigzero += 1








file.close()



refLength = [int(i)for i in refLength]
#Where in the list the longest contig occurs
#print(refLength.index(max(refLength)))
#Max contig length
#print(max(refLength))
readsmapped = [int(i)for i in readsmapped]
x = readsmapped[refLength.index(max(refLength))]

#print("Reads mapped to longest contig:","\t",x)


readsmapped = sorted(readsmapped)

sys.stdout = open("WeeSamContigs_Stats.txt","w")
print("Total Number of reads Mapped"+","+"Reads Mapped to Longest Contig"+","+"Mean number of reads mapped"+","+"Median Number of Reads Mapped")

print(str(sum(readsmapped))+","+str(x)+","+str(statistics.mean(readsmapped))+","+str(statistics.median(readsmapped)))

print("\n"+"Total Assembly length"+","+"N50"+","+"Genome Length"+","+"NG50")


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

    if ng50total > int(genomeLength) / 2:
        ng50list.append(i)




print(str(assemblyLength)+","+str(n50list[0])+","+str(genomeLength)+","+str(ng50list[0]))

print("\n"+"Number of contigs with 5 or less reads mapped"+","+"Number of contigs with no reads mapped")
print(str(contigbelow5)+","+str(contigzero))



