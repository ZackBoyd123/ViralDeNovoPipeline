#!/usr/bin/python3
import sys
import statistics


with open(sys.argv[1]) as file:

    seq=[]

    for line in file:
        if line.startswith(">"):
            pass

        else:
            seq.append(line)

    seq="".join(seq)
    gcount = 0
    ccount = 0
    acount = 0
    tcount = 0
    ncount = 0

    total = 0
    for i in seq:
        total += 1

        if i.startswith("G"):
            gcount += 1
        if i.startswith("C"):
            ccount += 1
        if i.startswith("T"):
            tcount += 1
        if i.startswith("A"):
            acount += 1
        if i.startswith("N"):
            ncount += 1


    dinucleotide = [seq[i:i+2] for i in range(0,len(seq),2)]

    gcdinuc = 0
    gc=0
    for i in dinucleotide:
        gc +=1
        if i == "GC":
            gcdinuc += 1



    print('Genome Length'+"\t"+ str(total)+"\n"+"_______________________________________")

    print('G(%):'+"\t"+ str((int(gcount)/int(total)*100)))
    print('C(%):'+"\t"+ str((int(ccount)/int(total)*100)))
    print("T(%):"+"\t"+ str((int(tcount)/int(total)*100)))
    print("A(%):"+"\t"+ str((int(acount)/int(total)*100))+"\n"+"_______________________________________")
    print('Total GC dinucleotides:'+"\t"+str(gcdinuc))
    gccontent = ((int(gcount)+ int(ccount)) / int(total)) * 100
    print("GC Content:"+"\t"+str(gccontent)+"\n")





    seq = list(seq)
    if ncount > 1:
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
                    totNumOfN += 1
            except IndexError:
                pass
    
        lengthlist.append(nlist[-1])
        #print(nlist)
        #print(lengthlist)
        
    
        testlist = [lengthlist[i:i+2]for i in range(0, len(lengthlist), 2)]
        #print(testlist)
        finalList =[]
        for i in testlist:
            i = i[1] - i[0]
            finalList.append(i)
    
        finalList = sorted(finalList, reverse=True)

    
        print("_______________________________________")
        print("Number of N windows:", "\t", totNumOfN)
        print("Number of Ns in file:"+"\t"+str(ncount))
        print("Median length of N windows:","\t",str(statistics.median(finalList)))
        print("Mean length of N windows:","\t",str(statistics.mean(finalList)))
