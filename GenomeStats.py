#!/usr/bin/python3
import sys
import statistics
import csv
import argparse

parser = argparse.ArgumentParser(description="Give me -I (--repeats,--coverage,--contig)")
parser.add_argument("-I","--input", help= "Input fasta file for genome analysis. If this file was the result of mdust \
                                           some extra analyis will occur. This input file is required for all stages \
                                          of analysis", required=True)
parser.add_argument("--repeats", help="A file containing repeat coordinates on a genome, produced by mummer. If you \
                                      specifiy this file, then you also need to specifiy a coverage file and a file \
                                      containing coordianted where contigs blast to a reference.")
parser.add_argument("--coverage", help="Coverage file containing a genome position and it's associated coverage")
parser.add_argument("--contig",help="File containing coordinates where a contig blasts to your ref file")
args = parser.parse_args()

inputFasta = str(args.input)
goodtogo = False
inputRepeats = ""
inputCoverage = ""
inputContig = ""

if args.repeats == None:
    inputRepeats_exists = False
else:
    inputRepeats_exists = True
    inputRepeats = str(args.repeats)

if inputRepeats_exists:
    if args.coverage == None:
        print("You specified --repeats, which means you must also specify --coverage and --conitg")
        exit(1)
    if args.contig == None:
        print("You specified --repeats, which means you must also specify --coverage and --conitg")
        exit(1)
    else:
        inputCoverage = str(args.coverage)
        inputContig = str(args.contig)
        goodtogo = True

with open(inputFasta) as file:
    seq=[]
    for line in file:
        if line.startswith(">"):
            pass
        else:
            seq.append(line)

    seq="".join(seq)
    seq=seq.replace("\n","")
    seq=seq.rstrip("\n")
    file.close()

print('Genome Length'+"\t"+ str(len(seq))+"\n"+"_______________________________________")
print('G(%):'+"\t"+ str((int(seq.count("G"))/int(len(seq))*100)))
print('C(%):'+"\t"+ str((int(seq.count("C"))/int(len(seq))*100)))
print("T(%):"+"\t"+ str((int(seq.count("T"))/int(len(seq))*100)))
print("A(%):"+"\t"+ str((int(seq.count("A"))/int(len(seq))*100))+"\n"+"_______________________________________")
print('Total CpG Sites:'+"\t"+str(seq.count("CG")))
gccontent = ((int(seq.count("G"))+ int(seq.count("C"))) / int(len(seq))) * 100
print("GC Content:"+"\t"+str(gccontent)+"\n")

seq = list(seq)
testlist = []
if seq.count("N") > 1:
    nlist = []
    for i, j in enumerate(seq):
        if j == "N":
            nlist.append(i)

    totNumOfN = 1
    lengthlist =[nlist[0]]

    for i in nlist:
        try:
            if int(i) - int(nlist[nlist.index(i)+1]) == -1:
                pass
            else:
                lengthlist.append(i)
                lengthlist.append(nlist[nlist.index(i)+1])

        except IndexError:
            pass
    lengthlist.append(nlist[-1])

    testlist = [lengthlist[i:i+2]for i in range(0, len(lengthlist), 2)]

    repeatRegions=[]
    for i in testlist:

        repeatRegions.append((list(range(int(i[0]),int(i[1])))))

    newNRegion = [item for sublist in repeatRegions for item in sublist]

    lst=[]
    for i, j in enumerate(seq):
        lst.append(i)

    finalList =[]

    for i in testlist:

        i = i[1] - i[0]
        finalList.append(i)

    finalList = sorted(finalList, reverse=True)


    print("_______________________________________")
    print("Number of N windows:", "\t", str(len(testlist)))
    print("Number of Ns in file:"+"\t"+str(seq.count("N")))
    print("Median length of N windows:","\t",str(statistics.median(finalList)))
    print("Mean length of N windows:","\t",str(int(seq.count("N"))/int(len(testlist))))


if goodtogo:
    repeatList=[]
    totlines = 0
    with open(inputRepeats) as f2:
        data = csv.reader(f2, delimiter="\t")
        for i,line in enumerate(data):
            totlines += 1
            if i != 1:
                repeatList.append(line)

        f2.close()

    del repeatList[0]
    repeatList = [item for sublist in repeatList for item in sublist]
    repeatList = ",".join(repeatList)
    repeatList = " ".join(repeatList.split())
    repeatList = repeatList.split(" ")

    sys.stdout=open(inputFasta.split(".")[0]+"_stats.csv","w")
    print("Figure"+","+"X-Min"+","+"X-Max"+","+"Y-Min"+","+"Y-Max"+","+"Plot")
    print("Genome"+","+str(0)+","+str(len(seq))+","+str(0)+","+str(3)+","+"GenomePlot")
    for i in testlist:
        print("Low Entropy"+","+str(i[0])+","+str(i[1])+","+str(0)+","+str(3)+","+"GenomePlot")


    def chunks(l, n):
        for i in range(0, len(l), n):
            yield l[i:i + n]

    x = (chunks(repeatList, 3))
    if totlines >= 3:
        for i in x:
            i[2] = i[2].replace(",", "")
            if i[1].endswith("r"):
                i[1] = i[1].replace("r","")
                print("RepeatForward"+","+str(i[0])+","+str(int(i[0])+int(i[2]))+","+str(0)+","+str(3)+","+"GenomePlot")
                print("RepeatReverse"+","+str(i[1])+","+str(int(i[1])-int(i[2]))+","+str(0)+","+str(3)+","+"GenomePlot")
            else:
                print("RepeatForward"+","+str(i[0])+","+str(int(i[0])+int(i[2]))+","+str(0)+","+str(3)+","+"GenomePlot")
                print("RepeatForward"+","+str(i[1])+","+str(int(i[1])+int(i[2]))+","+str(0)+","+str(3)+","+"GenomePlot")


    print("A%" + "," + (
                    "%.2f" % (int(seq.count("A")) / int(len(seq)) * 100) + "," + str(0) + "," + str(0) + "," + str(
                0) + "," + "PiePlot"))
    print("T%" + "," + (
                    "%.2f" % (int(seq.count("T")) / int(len(seq)) * 100) + "," + str(0) + "," + str(0) + "," + str(
                0) + "," + "PiePlot"))
    print("C%"+","+("%.2f" % (int(seq.count("C"))/int(len(seq))*100)+","+str(0)+","+str(0)+","+str(0)+","+"PiePlot"))
    print("G%" + "," + (
                "%.2f" % (int(seq.count("G")) / int(len(seq)) * 100) + "," + str(0) + "," + str(0) + "," + str(
            0) + "," + "PiePlot"))





    coveragelist=[]
    linesinfile = 0
    with open(inputCoverage) as f3:
        data = csv.reader(f3, delimiter="\t")
        for line in data:
            linesinfile += 1
            coveragelist.append(line[1:3])

        f3.close()

    for i,j in enumerate(coveragelist):
        if i % 200 == 0:
            print("Coverage"+","+str(j[0])+","+str(j[1])+","+str(0)+","+str(0)+","+"CoveragePlot")

    #sys.stdout= sys.__stdout__

    contigList=[]
    qlenList=[]
    with open(inputContig) as f4:
        data = csv.reader(f4, delimiter="\t")
        for line in data:
            contigList.append(line[2:8])



        f4.close()

    for i in contigList:

        if int(i[0]) == int(i[1]):
            if int(i[4]) >= int(i[5]):
                print("PerfectAlignReverse" + "," + str(i[4]) + "," + str(i[5]) + "," + str(4.1) + "," + str(5) + "," + "ContigPlot")
                contigList.remove(i)

            else:
                print("PerfectAlignForward" + "," + str(i[4]) + "," + str(i[5]) + "," + str(4.1) + "," + str(5) + "," + "ContigPlot")
                contigList.remove(i)

    #print(len(contigList))
    for i in contigList:

        if int(i[4]) >= int(i[5]):
            print("AlignReverse" + "," + str(i[4]) + "," + str(i[5]) + "," + str(1) + "," + str(2) + "," + "ContigPlot")
        else:
            print("AlignForward"+","+str(i[4])+","+str(i[5])+","+str(2.1)+","+str(3)+","+"ContigPlot")


