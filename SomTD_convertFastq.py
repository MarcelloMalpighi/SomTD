#!/usr/bin/env python3
### 2022.10.25
### 2023.04.11
### change lowFI to SomTD
### 2023.06.30
### improve names of OptionParser parameters
### 2023.07.25
### change name of reads because .1 and .2 has been added in previous steps
### 2023.09.24
### update for using in both clip and disc situations



# import
import pysam
from optparse import OptionParser
from SomTD_ReadGT import complseq

# read in parameters
parser=OptionParser()
parser.add_option('-i', dest='input')
parser.add_option('-o', dest='output')
parser.add_option('-t', dest='dataType', default='clip')
(options, args) = parser.parse_args()
inputFileLoc=options.input
if options.dataType == "clip":
    outputFileLoc1=options.output + ".genome.clip.pair.fastq.R1"
    outputFileLoc2=options.output + ".genome.clip.pair.fastq.R2"
elif options.dataType == "disc":
    outputFileLoc1=options.output + ".genome.disc.pair.fastq.R1"
    outputFileLoc2=options.output + ".genome.disc.pair.fastq.R2"

# output fastq of overlapping reads
inputFile=pysam.AlignmentFile(inputFileLoc,"r") # input soft-clipped read pairs genome sam file
outputFile1=open(outputFileLoc1,"w") # output soft-clipepd read pairs genome fastq file
outputFile2=open(outputFileLoc2,"w")
for read in inputFile:
    if read.get_tag("XO") == "1":
        outputFile1.write("@"+read.query_name + " XO:Z:1\n")
        if read.flag & 16 == 16:
            outputFile1.write(complseq(read.query_sequence)+"\n"+"+\n"+read.qual[::-1]+"\n")
        else:
            outputFile1.write(read.query_sequence+"\n"+"+\n"+read.qual+"\n")
    else:
        outputFile2.write("@"+read.query_name+" XO:Z:2\n")
        if read.flag & 16 == 16:
            outputFile2.write(complseq(read.query_sequence)+"\n"+"+\n"+read.qual[::-1]+"\n")
        else:
            outputFile2.write(read.query_sequence+"\n"+"+\n"+read.qual+"\n")

# close files
inputFile.close()
outputFile1.close()
outputFile2.close()