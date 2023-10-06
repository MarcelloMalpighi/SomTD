#!/usr/bin/env python3
### 2023.07.28
### a script for converting cutadapt results to pseudoSAM format
### 2023.07.31
### deal with 2 kinds of overlapping
### 2023.09.29
### fix the bug of ignoring deleting 0bp genome part output

from optparse import OptionParser

def convertDirec(direc):
    """
    function for converting FOR and REV to flag 0 and 16
    """
    if direc == "FOR":
        return 0
    elif direc == "REV":
        return 16
    else:
        return None


parser=OptionParser()
parser.add_option('-t', dest='TEtaileCutadaptFileLoc')
parser.add_option('-a', dest='polyAcutadaptFileLoc')
parser.add_option('-r', type = 'int', dest='readLen')
parser.add_option('-o', dest='output')
(options, args) = parser.parse_args()
TEtailCutadaptFileLoc=options.TEtaileCutadaptFileLoc
polyAcutadaptFileLoc=options.polyAcutadaptFileLoc
readLen=options.readLen
output=options.output

TEtailCutadaptFile=open(TEtailCutadaptFileLoc,"r")
polyAcutadaptFile=open(polyAcutadaptFileLoc,"r")
outputFile=open(output+".cutadapt.txt","w")

outputFile.write("@SQ\tSN:polyA\tLN:100\n")
outputDict={}

NR = 1
for i in polyAcutadaptFile:
    if NR % 4 == 1:
        lineList=i.strip("\n").split("_")
        readName = lineList[0][1:] + "_" + lineList[1][-1]
        polyAflag = convertDirec(lineList[2])
        polyAlen=len(lineList[3].split("=")[1])
        if not outputDict.__contains__(readName): # determine whether the other side polyA exist
            outputDict[readName] = [[str(polyAflag),"polyA","1","60","","*","0","0","XP:Z:polyA," + str(polyAflag) + ",1," + str(polyAlen)],["2052","*","0","0","*","*","0","0","XP:Z:NULL"]] # remove @ in the read name, sequence and phred could be retrieved from genome alignment
            secondCombine = 0
        else:
            outputDict[readName][1] = [str(2048 + polyAflag),"polyA","1","60","","*","0","0","XP:Z:polyA," + str(polyAflag) + ",1," + str(polyAlen)]
            secondCombine = 1
    elif NR % 4 == 2:
        genomeLen = len(i.strip("\n"))
        if secondCombine == 0:
            if genomeLen != 0:
                outputDict[readName][0][4] = str(polyAlen) + "M" + str(genomeLen) + "S"
            else:
                outputDict[readName][0][4] = str(polyAlen) + "M"
        else:
            if genomeLen != 0:
                outputDict[readName][1][4] = str(polyAlen) + "M" + str(genomeLen) + "S"
            else:
                outputDict[readName][1][4] = str(polyAlen) + "M"
    NR = NR + 1

NR = 1
for i in TEtailCutadaptFile:
    if NR % 4 == 1:
        lineList=i.strip("\n").split("_")
        readName = lineList[0][1:] + "_" + lineList[1][-1]
        TEflag = convertDirec(lineList[2])
        TEtailLen = len(lineList[3].split("=")[1])
        TEtype = str(lineList[3].split("=")[0])
        polyAlen=len(lineList[4].split("=")[1])
        if outputDict.__contains__(readName): # Deal with overlapping records, replace polyA results with TEtail results for the latter is more reliable
            for j in range(len(outputDict[readName])): # replace
                if int(outputDict[readName][j][0]) & 16 == TEflag & 16 or int(outputDict[readName][j][0]) & 4 == 4:
                    outputDict[readName][j] = [str(TEflag),"polyA","1","60","","*","0","0","XP:Z:" + TEtype + "," + str(TEflag) + ",1," + str(TEtailLen) + ";polyA," + str(TEflag) + "," + str(TEtailLen+1) + "," + str(TEtailLen + polyAlen)]
                    replaceMark = j
                    break
        else: # add
            outputDict[readName]=[[str(TEflag),"polyA","1","60","","*","0","0","XP:Z:" + TEtype + "," + str(TEflag) + ",1," + str(TEtailLen) + ";polyA," + str(TEflag) + "," + str(TEtailLen+1) + "," + str(TEtailLen + polyAlen)],["2052","*","0","0","*","*","0","0","XP:Z:NULL"]]
            replaceMark = "add"
    elif NR % 4 == 2:
        genomeLen = len(i.strip("\n")) # sequence which is beyond TE tail and polyA
        if replaceMark == "add": # add
            if genomeLen != 0:
                outputDict[readName][0][4] = str(TEtailLen) + "S" + str(polyAlen) + "M" + str(genomeLen) + "S" # no need to change CIGAR according to the strandness, i.e., polyA or polyT because SAM shows forward strand result
            else:
                outputDict[readName][0][4] = str(TEtailLen) + "S" + str(polyAlen) + "M"
        else: # replace
            if genomeLen != 0:
                outputDict[readName][j][4] = str(TEtailLen) + "S" + str(polyAlen) + "M" + str(genomeLen) + "S"
            else:
                outputDict[readName][j][4] = str(TEtailLen) + "S" + str(polyAlen) + "M"
    NR = NR + 1

### output results
for i in sorted(outputDict.keys()):
    for j in outputDict[i]:
        outputFile.write(i + "\t" + "\t".join(j) + "\n")
