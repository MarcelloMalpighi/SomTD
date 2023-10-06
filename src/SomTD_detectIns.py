#!/usr/bin/env python3
### 2023.09.28
###### a script for detecting insertions from collected clipped read pairs

from SomTD_ReadGTP import *
import pysam
import os
from optparse import OptionParser
parser=OptionParser()
parser.add_option('-i', '--input', type = 'str', dest = 'inputFileLoc', help = "input")
parser.add_option('-o', '--output', type = 'str', dest = 'output', help = "output")
(options, args) = parser.parse_args()
inputFileLoc,output = options.inputFileLoc,options.output


inputFile=pysam.AlignmentFile(inputFileLoc,"r")
testList=[]
interList=[]
NR = 1
for i in inputFile:
    interList.append(i)
    if NR % 12 == 0:
        testList.append(interList)
        interList = []
    NR = NR + 1

inputFile.close()

### test output
ref2NameList = [] # generate (1) accessory 2ref read name list
resultDict = {} # generate (2) read information list, key: read name, value: read information lists including 1 or 2 arboutputlist
resultFile = open(output + ".bed","w")
discNameFile = open(output + ".clip2disc.fragment.names","w")
NR = 1
for i in testList:
    testGTP=ReadGTP(i,150,TElenDict)
    testGTP.arbMode()
    if testGTP.arbDisc:
        _=discNameFile.write(testGTP.alignList[0].query_name + "\n")
    for j in testGTP.arbOutputList:
        if j[3] in resultDict.keys():
            resultDict[j[3]].append(j)
        else:
            resultDict[j[3]] = [j]
        if j[4].endswith("2"):
            ref2NameList.append(j[3])
    NR = NR + 1

ref2NameList = list(set(ref2NameList))
discNameFile.close()
for i in sorted(resultDict.keys()):
    for j in resultDict[i]:
        _=resultFile.write("\t".join(list(map(str,j))) + "\n")

resultFile.close()

# generate (3) insertion information dict (4) read to insertion dict
os.system("sort -k 1,1 -k 2,2n " + output + ".bed > " + output + ".sorted.bed")
os.system("bedtools merge -d 100 -c 4,4,5,6,7,8,9,10,11,12 -o count,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse -i " + output + ".sorted.bed > " + output + ".sorted.merged.bed")
read2InsDict = {} # key: read name, value: insertion name
insInfoDict = {} # key insertion name, value: insertion information
mergedFile = open( output + ".sorted.merged.bed","r")
for i in mergedFile:
    interList = i.strip("\n").split("\t")
    insInfoDict[str(interList[0]) + "_" + str(interList[1]) + "_" + str(interList[2])] = interList
    for readName in interList[4].split(","):
        if readName in read2InsDict.keys():
            read2InsDict[readName].append(str(interList[0]) + "_" + str(interList[1]) + "_" + str(interList[2]))
        else:
            read2InsDict[readName] = [str(interList[0]) + "_" + str(interList[1]) + "_" + str(interList[2])]

mergedFile.close()

### choose better location with more conclusive reads for accessory 2ref reads
def getConcluNum(insName):
    """
    function for getting the number of conclusive read of a insertion
    """
    concluNum = 0
    for i in insInfoDict[insName][6].split(","):
        if i == "CONCLUSIVE":
            concluNum += 1
    return concluNum

for i in ref2NameList:
    ins1,ins2=read2InsDict[i][0],read2InsDict[i][1]
    if getConcluNum(ins1) >= getConcluNum(ins2):
        for j in range(len(insInfoDict[ins2][4].split(","))):
            if insInfoDict[ins2][4].split(",")[j] == i:
                delID = insInfoDict[ins2][5].split(",")[j]
        for j in range(len(resultDict[i])):
            if resultDict[i][j][4] == delID:
                del resultDict[i][j]
                break
    else:
        for j in range(len(insInfoDict[ins1][4].split(","))):
            if insInfoDict[ins1][4].split(",")[j] == i:
                delID = insInfoDict[ins1][5].split(",")[j]
        for j in range(len(resultDict[i])):
            if resultDict[i][j][4] == delID:
                del resultDict[i][j]
                break

newOutputFile=open( output+ ".new.bed","w")
for i in sorted(resultDict.keys()):
    for j in resultDict[i]:
        _=newOutputFile.write("\t".join(list(map(str,j))) + "\n")
newOutputFile.close()

os.system("sort -k 1,1 -k 2,2n " + output + ".new.bed > " + output + ".new.sorted.bed")
os.system("bedtools merge -d 100 -c 4,4,5,6,7,8,9,10,11,12 -o count,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse -i " + output + ".new.sorted.bed > " + output + ".new.sorted.merged.bed")

newMergedInsFile = open(output + ".new.sorted.merged.bed")
outputMergedInsFile = open(output + ".new.sorted.merged.conclu.morethan1.bed","w")
for i in newMergedInsFile:
    interList = i.strip("\n").split("\t")
    concluMark = False
    for j in interList[6].split(","):
        if j.startswith("CONCLUSIVE"): # at lease one conclusive
            concluMark =True
    headNum,tailNum,headTruncNum,tailTruncNum,wrongNum,allNum,truncMark=0,0,0,0,0,len(interList[6].split(",")),True # not all single end truncation
    for k in interList[10].split(","):
        if k == "head":
            headNum += 1
        elif k == "tail":
            tailNum += 1
        elif k == "headTruncation":
            headTruncNum += 1
        elif k== "tailTruncation":
            tailTruncNum += 1
        else:
            wrongNum += 1
    if wrongNum == allNum or allNum == headTruncNum or allNum == tailTruncNum:
        truncMark = False
    if concluMark and truncMark and allNum >= 2:
        _=outputMergedInsFile.write(i)

newMergedInsFile.close()
outputMergedInsFile.close()
