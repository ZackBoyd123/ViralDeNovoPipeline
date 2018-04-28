#!/usr/bin/python3

import re
import sys
import os
fileIn = os.environ['currentpath']
fileOut = "vicuna_config.txt"
path = os.environ['otherpath']

sys.stdout = open(fileOut, "w")

for line in open(fileIn, 'r').readlines():
    line = re.sub(r'outputDIR.+', r'outputDIR\t{0}/VicunaOutput/'.format(path), line)
    line = re.sub(r'pFqDir.+', r'pFqDir\t{0}/VicunaOutput'.format(path), line)
    print(line)
