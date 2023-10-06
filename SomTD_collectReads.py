#!/usr/bin/env python3
### 2023.04.17
### a script for collecting soft-clipped and discordant read pairs in genome alignments
### 2023.06.13
### changed for receiving output from samtools
### 2023.06.19
### add the standness condition for discordant fragments
### 2023.06.30
### combine soft-clipped and discordant scripts
### change names of output files
### improve names of OptionParser parameters
### 2023.07.26
### alignments XO tag which indicate R1 R2 has been added in the previous step, so no need
### to add \1 or \2 to the read name
### 2023.08.02
### fix a bug about outpuing UNCLIP instead of DISC for the first read when the second read is soft-clipped
### 2023.09.12
### use SomTD_ReadGTP instead of SomTD_ReadGT and replace checkSoftClipped with checkClipped
### 2023.09.13
### mate read of soft-clipped read should be mapped now for convenience of following operations
### new result is the same as original testSup.genome.clip.pair.sam, so use the original one



# import
import pysam
from optparse import OptionParser
from SomTD_ReadGTP import checkClipped

# read in parameters
parser=OptionParser()
parser.add_option('-o', dest='output')
parser.add_option('-c', dest='cutoff')
parser.add_option('-m', dest='meanLen')
parser.add_option('-s', dest='stdLen')
(options, args) = parser.parse_args()
outputClipFileLoc=options.output + ".genome.clip.pair.sam"
outputDiscFileLoc=options.output + ".genome.disc.pair.sam"
cutoff=int(options.cutoff)
meanLen=int(options.meanLen)
stdLen=int(options.stdLen)

# read in files
inputFile=pysam.AlignmentFile("-","r") # input genome sam file
outputClipFile=pysam.AlignmentFile(outputClipFileLoc,"w",template=inputFile) # output soft-clipped read pairs genome sam file
outputDiscFile=pysam.AlignmentFile(outputDiscFileLoc,"w",template=inputFile) # output discordant read pairs genome sam file

# collect discordant read pairs
NR = 1
for read in inputFile:
    if NR % 2 == 1:
        if checkClipped(read.cigartuples,cutoff=cutoff) and not read.mate_is_unmapped:# soft-clipped check
            read.set_tag("XG","CLIP","Z")
            checkTag = 1 # soft-clipped read pair
        elif read.is_paired and not read.is_unmapped and not read.mate_is_unmapped and (read.reference_id != read.next_reference_id or abs(read.reference_start - read.next_reference_start) >= meanLen + 3 * stdLen or (read.is_forward == read.mate_is_forward and read.is_reverse == read.mate_is_reverse)): # discordant check: chromosome comparison in 3stdLen and strandness is not necessary because python will not process the following conditions if the previous one is False
        # 2 reads are both mapped
        # 2 reads without proper soft-clipped parts
        # 2 reads are aligned to 2 different chromosomes or
        # 2 reads are aligned to the same chromosomes but the disctance between 2 reads is more than mean fragment length + 3 std or
        # 2 reads are aligned to the same chromosomes but the standness is incompatible
            read.set_tag("XG","DISC","Z")
            checkTag = 2 # discordant read pair
        else:
            read.set_tag("XG","UNCLIP","Z")
            checkTag = 0 # other read pair
        interline = read
    else:
        if checkClipped(read.cigartuples,cutoff=cutoff) and not read.mate_is_unmapped:
            read.set_tag("XG","CLIP","Z")
            if interline.get_tag("XG") != "CLIP":
                interline.set_tag("XG","UNCLIP","Z")
            outputClipFile.write(interline)
            outputClipFile.write(read)
        elif checkTag == 1:
            read.set_tag("XG","UNCLIP","Z")
            outputClipFile.write(interline)
            outputClipFile.write(read)
        elif checkTag == 2:
            read.set_tag("XG","DISC","Z")
            outputDiscFile.write(interline)
            outputDiscFile.write(read)
    NR = NR + 1

# close files
inputFile.close()
outputClipFile.close()
outputDiscFile.close()