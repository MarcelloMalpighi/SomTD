#!/usr/bin/env python3
### 2023.07.25
### a script for combining alignments of Genome primary alignments, TE primary alignments, 
### Genome supplementary alignments(all 2048 alignments including not clipped pairs), TE supplementary alignments
### (only clipped read pairs), constructed polyA primary pseudoalignments and constructed polyA supplementary pseudoalignments
### format:
### @header
### G primary alignment
### G supplementary alignment
### T primary alignment
### T supplementary alignment
### P primary pseudoalignment
### P supplementary pseudoalignment
### 2023.07.26
### try to combine polyA constructed alignment
### 2023.08.15
### fix the bug that pseudoSAM alignments without polyA tail are not tagged with XP
### 2023.09.12
### XG for Gsup and curated Gsup, XT for Ppri, Psup and curated Psup are added
### use SomTD_ReadGTP instead of SomTD_ReadGT and replace checkSoftClipped with checkClipped
### 2023.09.13
### add standard for choosing Gsup by clipped part number and relative location with mate read and Tsup
### by soft-clipped part number

import pysam
import os
from optparse import OptionParser
from SomTD_ReadGTP import complseq
from SomTD_ReadGTP import checkClipped
from SomTD_ReadGTP import checkClippedNum

def getMateName(readName):
    """
    function for getting mate read name
    """
    mateName = readName[:-1] + "2" if readName.endswith("1") else readName[:-1] + "1"
    return mateName



parser=OptionParser()
parser.add_option('-g', dest='GpriFileLoc')
parser.add_option('-t', dest='TpriFileLoc')
parser.add_option('-G', dest='GsupFileLoc')
parser.add_option('-T', dest='TsupFileLoc')
parser.add_option('-c', dest='cutadaptFileLoc')
parser.add_option('-o', dest='output')
(options, args) = parser.parse_args()
GpriFileLoc=options.GpriFileLoc
TpriFileLoc=options.TpriFileLoc
GsupFileLoc=options.GsupFileLoc
TsupFileLoc=options.TsupFileLoc
cutadaptFileLoc=options.cutadaptFileLoc
outputClipFileLoc=options.output + ".combine.clip.pair.sam"

GpriFile=pysam.AlignmentFile(GpriFileLoc,"r")
TpriFile=pysam.AlignmentFile(TpriFileLoc,"r")
GsupFile=pysam.AlignmentFile(GsupFileLoc,"r")
TsupFile=pysam.AlignmentFile(TsupFileLoc,"r")
cutadaptFile=open(cutadaptFileLoc,"r")
os.system("grep -h \"^@\" " + GpriFileLoc + " " + TpriFileLoc + " " + GsupFileLoc + " " + TsupFileLoc + " " + cutadaptFileLoc + " | sort | uniq > " + outputClipFileLoc)
outputFile=open(outputClipFileLoc,"a")

### read in alignments and combine to a dict
GDict={}
TDict={}
GSupDict={}
TSupDict={}
seqPhredDict={}
outputDict={}

### Genome alignments
for i in GpriFile: # Genome primary alignments
    GDict[i.query_name + "_" + i.get_tag("XO")] = [i]
    seqPhredDict[i.query_name + "_" + i.get_tag("XO")] = [i.flag,i.query_sequence,i.tostring().split("\t")[10]]

# Genome supplementary alignments
for i in GsupFile: # generate Genome supplementary alignments dict
    if i.query_name + "_" + i.get_tag("XO") in GSupDict.keys():
        GSupDict[i.query_name + "_" + i.get_tag("XO")].append(i)
    else:
        GSupDict[i.query_name + "_" + i.get_tag("XO")] = [i]

for i in GSupDict.keys(): # add best Genome supplementary alignments
    if i in GDict.keys():
        GSupNum = len(GSupDict[i])
        selectedID = 0
        if GSupNum != 1: # more than one Genome supplementary alignment
            for j in range(GSupNum):
                if checkClippedNum(GSupDict[i][j].cigartuples,10) == 1 and abs(GSupDict[i][j].reference_start - GDict[getMateName(i)][0].reference_start) <= 700 and GSupDict[i][j].reference_name == GDict[getMateName(i)][0].reference_name: # good: one clipped part and close to mate read
                    selectedID = j
            if selectedID == 0:
                for j in range(GSupNum):
                    if checkClippedNum(GSupDict[i][j].cigartuples,10) == 1: # good: one clipped part, maybe mate read aligned to reference insertion
                        selectedID = j
            if selectedID == 0:
                for j in range(GSupNum): # good: better accordance on location
                    if abs(GSupDict[i][j].reference_start - GDict[getMateName(i)][0].reference_start) <= 700 and GSupDict[i][j].reference_name == GDict[getMateName(i)][0].reference_name:
                        selectedID = j
        GsupAdd = GSupDict[i][selectedID]
        if checkClipped(GsupAdd.cigartuples,10):
            GsupAdd.set_tag("XG","CLIP","Z")
        else:
            GsupAdd.set_tag("XG","UNCLIP","Z")
        GDict[i].append(GsupAdd)

for i in sorted(GDict.keys()): # deal with primary alignments without supplementary alignments
    if len(GDict[i]) == 1:
        newList=GDict[i][0].tostring().split("\t")[0:12]
        if int(newList[1]) & 16 == 16: # keep original direction of the read
            newList[9] = complseq(newList[9])
            newList[10] = "".join(list(reversed(newList[10])))
        newList.append("XO:Z:" + GDict[i][0].get_tag("XO") + "\tXG:Z:UNCLIP")
        newList[1:9] = '2052','*','0','0','*','*','0','0' # clumsy single mode simulation
        GDict[i].append("\t".join(newList))
    outputDict[i] = []
    for j in GDict[i]:
        if type(j) is str:
            outputDict[i].append(j)
        else:
            outputDict[i].append(j.tostring())

GpriFile.close()
GsupFile.close()

### TE alignments
for i in TpriFile: # TE primary alignments
    if checkClipped(i.cigartuples,10):
        i.set_tag("XT","CLIP","Z")
    else:
        i.set_tag("XT","UNCLIP","Z")
    TDict[i.query_name + "_" + i.get_tag("XO")] = [i]

for i in TsupFile: # generate Genome supplementary alignments dict
    if i.query_name + "_" + i.get_tag("XO") in TSupDict.keys():
        TSupDict[i.query_name + "_" + i.get_tag("XO")].append(i)
    else:
        TSupDict[i.query_name + "_" + i.get_tag("XO")] = [i]

for i in TSupDict.keys(): # add best TE supplementary alignments
    if i in TDict.keys():
        TSupNum = len(TSupDict[i])
        selectedID = 0
        if TSupNum != 1: # more than one TE supplementary alignment
            for j in range(TSupNum):
                if checkClippedNum(TSupDict[i][j].cigartuples,10) == 1 and abs(TSupDict[i][j].reference_start - TDict[getMateName(i)][0].reference_start) <= 700 and TSupDict[i][j].reference_name == TDict[getMateName(i)][0].reference_name: # good: one clipped part and close to mate read
                    selectedID = j
            if selectedID == 0:
                for j in range(TSupNum):
                    if checkClippedNum(TSupDict[i][j].cigartuples,10) == 1: # good: one clipped part, maybe clipped read and mate read come from reference insertion and new insertion respectively
                        selectedID = j
            if selectedID == 0:
                for j in range(TSupNum): # good: better accordance on location
                    if abs(TSupDict[i][j].reference_start - TDict[getMateName(i)][0].reference_start) <= 700 and TSupDict[i][j].reference_name == TDict[getMateName(i)][0].reference_name:
                        selectedID = j
        TsupAdd = TSupDict[i][selectedID]
    if checkClipped(TsupAdd.cigartuples,10):
        TsupAdd.set_tag("XT","CLIP","Z")
    else:
        TsupAdd.set_tag("XT","UNCLIP","Z")
    TDict[i].append(TsupAdd)

for i in sorted(TDict.keys()): # deal with primary alignments without supplementary alignments
    if len(TDict[i]) == 1: 
        newList=TDict[i][0].tostring().split("\t")[0:12]
        if int(newList[1]) & 16 == 16: # keep original direction of the read
            newList[9] = complseq(newList[9])
            newList[10] = "".join(list(reversed(newList[10])))
        newList.append("XO:Z:" + TDict[i][0].get_tag("XO") + "\tXT:Z:UNCLIP")
        newList[1:9] = '2052','*','0','0','*','*','0','0' # clumsy single mode simulation
        TDict[i].append("\t".join(newList))
    for j in TDict[i]:
        if type(j) is str:
            outputDict[i].append(j)
        else:
            outputDict[i].append(j.tostring())

TpriFile.close()
TsupFile.close()

### deal with cutadapt results, sequence and phred values should be added
NR = 1
for i in cutadaptFile:
    if i[0] != "@":
        lineList = i.strip("\n").split("\t")
        seqPhredList=seqPhredDict[lineList[0]]
        nameList = lineList[0].split("_")
        lineList[0]=nameList[0]
        lineList.append("XO:Z:" + str(nameList[1]))
        if int(lineList[1]) & 16 == int(seqPhredList[0]) & 16:
            outputDict["_".join(nameList)].append("\t".join(lineList[0:9]) + "\t" + seqPhredList[1] + "\t" + seqPhredList[2] + "\t" + lineList[10] + "\t" + lineList[9])
        elif int(lineList[1]) & 4 == 4:
            if int(seqPhredList[0]) & 16 == 16:
                outputDict["_".join(nameList)].append("\t".join(lineList[0:9]) + "\t" + complseq(seqPhredList[1]) + "\t" + "".join(list(reversed(seqPhredList[2]))) + "\t" + lineList[10] + "\t" + lineList[9])
            else:
                outputDict["_".join(nameList)].append("\t".join(lineList[0:9]) + "\t" + seqPhredList[1] + "\t" + seqPhredList[2] + "\t" + lineList[10] + "\t" + lineList[9])
        else:
            outputDict["_".join(nameList)].append("\t".join(lineList[0:9]) + "\t" + complseq(seqPhredList[1]) + "\t" + "".join(list(reversed(seqPhredList[2]))) + "\t" + lineList[10] + "\t" + lineList[9])
        NR = NR + 1

cutadaptFile.close()

for i in sorted(outputDict.keys()): # add null information for reads without cutadapt results
    if len(outputDict[i]) == 4:
        genomePriList = outputDict[i][0].split("\t")
        if int(genomePriList[1]) & 16 == 16:
            genomePriList[9] = complseq(genomePriList[9])
            genomePriList[10] = "".join(list(reversed(genomePriList[10])))
        outputDict[i].append("\t".join([genomePriList[0],"4","*","0","0","*","*","0","0",genomePriList[9],genomePriList[10],"XO:Z:" + str(i.split("_")[1]),"XP:Z:NULL"]))
        outputDict[i].append("\t".join([genomePriList[0],"2052","*","0","0","*","*","0","0",genomePriList[9],genomePriList[10],"XO:Z:" + str(i.split("_")[1]),"XP:Z:NULL"]))


### output results, only the first supplementary alignment is output
for i in sorted(outputDict.keys()):
    for j in outputDict[i]:
        _=outputFile.write(j + "\n")

outputFile.close()