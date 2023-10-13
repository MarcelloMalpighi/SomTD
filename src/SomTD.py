#!/usr/bin/env python3
### 2023.02.05
### 2023.03.28
### modified for a temporary version which can normally run with both single and double mode
### 2023.04.11
### particulars modified: SomTD name changed; help content added; selfFilter no longer receive a parameter;
### 2023.06.01
### modify bwa mem parameters
### 2023.06.30
### change to a python version
### 2023.07.21
### main script of SomTD in a python style
### change to a totally new pipeline
### 2023.09.28
### construct an intact version



import sys,os,datetime
from optparse import OptionParser



##################################
####### receive parameters #######
##################################

parser=OptionParser()
parser.add_option('--input1', type = 'str', dest = 'input1', help = "input file, bam/sam/fastq, please use .bam .sam .fq/.fastq as a suffix, mandatory") # need to be changed for adapting fastq
parser.add_option('--input2', type = 'str', dest = 'input2', help = "input file, if fastq, please use .fq/.fastq as a suffix, optional") # need to be changed for adapting fastq
parser.add_option('-o', '--outputPath', type = 'str', dest = 'outputPath', help = "output path, directory name will be output name, directory will be generated autonomously if not exist, mandatory")
parser.add_option('-c', '--cutoff', type = 'str', dest = 'cutoff', help = 'minimum soft-clipped length, limit itself is included, optional, default: 10')

parser.add_option('-f', '--frag', type = 'str', dest = 'fragLen', help = 'expected fragment length, mandatory')
parser.add_option('--std', type = 'str', dest = 'fragStd', help = 'standard deviation of fragment length, mandatory')
parser.add_option('-r', '--readLen', type = 'str', dest = 'readLen', help = 'read length, mandatory')

parser.add_option('-g', '--Greference', type = 'str', dest = 'Greference', help = 'genome reference sequence, fa, mandatory')
parser.add_option('-t', '--Treference', type = 'str', dest = 'Treference', help = 'transposon reference sequence, fa, mandatory')
parser.add_option('--TreferenceRecom', type = 'str', dest = 'TreferenceRecom', help = 'reversed complemented transposon reference sequence, fa, mandatory')
parser.add_option('-G', '--Gindex', type = 'str', dest = 'Gindex', help = 'genome reference sequence bwa index, mandatory') # need to be changed to optional
parser.add_option('-T', '--Tindex', type = 'str', dest = 'Tindex', help = 'transposon reference sequence bwa index, mandatory')
parser.add_option('-p', '--parallel', type = 'str', dest = 'parallel', help = 'number of threads, optional, default: 2') # need to be changed for a global control
parser.add_option('-m', '--memory', type = 'str', dest = 'memory', help = 'memory per thread used, optional, defalut: 2G') # need to be changed for a global memory control
# parser.add_option('-d', '--distRef', type = 'str', dest = 'distRef', help = 'minimum distance from reference insertions, nonreference insertions with a shorter one will be ignored, optional, default: 0')
# parser.add_option('-s', '--selfFilter', dest = 'selfFilter', action = 'store_true', help = 'filter insertions shared by all input samples, optional, default: disabled ')

(options, args) = parser.parse_args()
input1,input2,outputPath,Greference,Treference,Gindex,Tindex,fragLen,fragStd,readLen,TreferenceRecom = options.input1,options.input2,str(options.outputPath),str(options.Greference),str(options.Treference),str(options.Gindex),str(options.Tindex),str(options.fragLen),str(options.fragStd),str(options.readLen),str(options.TreferenceRecom)
cutoff = str(options.cutoff) if options.cutoff else "10"
parallel = str(options.parallel) if options.parallel else "2"
memory = str(options.memory) if options.memory else "2G"
# distRef = str(options.distRef) if options.distRef else "0"
# selfFilter = options.selfFilter if options.selfFilter else False

if input1 and outputPath and Greference and Treference and Gindex and Tindex:
    print("SomTD begins, time: "+ str(datetime.datetime.now()), file=sys.stderr)
else:
    print("Error: lack of mandatory parameters", file=sys.stderr)
    parser.print_help(file=sys.stderr)
    exit()

# check for paths
os.system("SomTDPath=`which SomTD | awk '{split($0,a,\"/\");for(i=2;i<=length(a)-2;i++){printf \"/%s\",a[i]}}'`") # get SomTD location
os.system("if [ ! -d " + outputPath + " ]; then mkdir " + outputPath + ";fi") # generate output directory if not exist
os.system("cd " + outputPath)
output=outputPath.split("/")[-1] # get output name



##################################
########## alignGenome ###########
##################################

if (input1.endswith("fq") or input1.endswith("fastq")) and (input2.endswith("fq") or input2.endswith("fastq")):
    fq1FileName,fq2FileName = output + ".R1.fq",output + ".R2.fq"
    fq1File,fq2File,input1File,input2File = open(fq1FileName,"w"),open(fq2FileName,"w"),open(input1,"r"),open(input2,"r")
    NR = 1
    for i in input1File:
        if NR%4==1:
            fq1File.write(i.strip("\n") + " XO:Z:1\n")
        else:
            fq1File.write(i)
        NR += 1
    NR = 1
    for i in input2File:
        if NR%4==1:
            fq2File.write(i.strip("\n") + " XO:Z:2\n")
        else:
            fq2File.write(i)
        NR += 1
    markdupSortedFileName = output+".markdup.sorted.bam"

    print(" ".join(["bwa","mem","-v","2","-t",parallel,"-C","-R","\"@RG\\tID:" + output + "\\tLB:" + output + "\\tPL:" + output + "\\tPU:" + output + "\\tSM:" + output + "\"",Gindex,fq1FileName,fq2FileName,"|","samblaster","-r","2>","markdup.sorted.log","|","samtools","view","-b","-u","-@",parallel,"|","samtools","sort","-@",parallel,"-m",memory,"-o",markdupSortedFileName]))
    os.system(" ".join(["bwa","mem","-v","2","-t",parallel,"-C","-R","\"@RG\\tID:" + output + "\\tLB:" + output + "\\tPL:" + output + "\\tPU:" + output + "\\tSM:" + output + "\"",Gindex,fq1FileName,fq2FileName,"|","samblaster","-r","2>","markdup.sorted.log","|","samtools","view","-b","-u","-@",parallel,"|","samtools","sort","-@",parallel,"-m",memory,"-o",markdupSortedFileName]))

    print(" ".join(["samtools","view","-u","-F","3328",markdupSortedFileName,"|","samtools","sort","-n","-@",parallel,"-m",memory,"-O","sam","|","SomTD_collectReads.py","-o",output,"-c",cutoff,"-m",fragLen,"-s",fragStd]))
    os.system(" ".join(["samtools","view","-u","-F","3328",markdupSortedFileName,"|","samtools","sort","-n","-@",parallel,"-m",memory,"-O","sam","|","SomTD_collectReads.py","-o",output,"-c",cutoff,"-m",fragLen,"-s",fragStd]))
else:
    print("SAM/BAM or other data format are not supported now",file=sys.stderr)
    exit()


##################################
########## collectReads ##########
##################################
# collect soft-clipped and discordant read pairs
os.system("SomTD_convertFastq.py -i " + output +".genome.clip.pair.sam -o " + output)
clipFq1FileName,clipFq2FileName = output + ".genome.clip.pair.fastq.R1",output + ".genome.clip.pair.fastq.R2"
print("SomTD collectReads done!", file=sys.stderr)
os.system(" ".join(["samtools view -u -F 1280 -f 2048",markdupSortedFileName,"| samtools sort -n -@",parallel,"-m",memory,"-O sam >",output+".genome.clip.sup.sam"]))



##################################
############ alignTE #############
##################################
# align to TE consensus sequence with bwa mem
os.system(" ".join(["bwa","mem","-v","2","-t",parallel,"-C","-R","\"@RG\\tID:" + output + "\\tLB:" + output + "\\tPL:" + output + "\\tPU:" + output + "\\tSM:" + output + "\"",Tindex,clipFq1FileName,clipFq2FileName,"|","samtools","view","-h","-F","2304",">",output + ".te.clip.ori.sam"]))
######### here to write one -k seems to influence greatly
### get TE supplementary alignments
os.system(" ".join(["samtools","view","-h","-F","2304",output + ".te.clip.ori.sam",">",output + ".te.clip.pri.sam"]))
os.system(" ".join(["samtools","view","-h","-F","256","-f","2048",output+".te.clip.ori.sam",">",output+".te.clip.sup.sam"]))
print("SomTD alignTE done!", file=sys.stderr)

##################################
############ cutadapt ############
##################################

# forward insertion tails
print(" ".join(['cutadapt','-j',parallel,'-e 0.1 -O 5 -g',"file:" + Treference,"--action=trim --untrimmed-output /dev/null --rename='{id}_XO:Z:1_FOR_{adapter_name}={match_sequence}'",clipFq1FileName,'-o','cutadapt.forwardTEtail.fastq.R1','>','cutadapt.TE.log']))
os.system(" ".join(['cutadapt','-j',parallel,'-e 0.1 -O 5 -g',"file:" + Treference,"--action=trim --untrimmed-output /dev/null --rename='{id}_XO:Z:1_FOR_{adapter_name}={match_sequence}'",clipFq1FileName,'-o','cutadapt.forwardTEtail.fastq.R1','>','cutadapt.TE.log']))
os.system(" ".join(['cutadapt','-j',parallel,'-e 0.1 -O 5 -g',"file:"+ Treference ,"--action=trim --untrimmed-output /dev/null --rename='{id}_XO:Z:2_FOR_{adapter_name}={match_sequence}'",clipFq2FileName,'-o','cutadapt.forwardTEtail.fastq.R2','>>', 'cutadapt.TE.log']))
os.system("cat cutadapt.forwardTEtail.fastq.R1 cutadapt.forwardTEtail.fastq.R2 > cutadapt.forwardTEtail.fastq")
os.system(" ".join(["cutadapt", '-j', parallel, "-e 0.1 -O 5 -g XA{100} --action=trim --rename='{id}_polyA={match_sequence}' --untrimmed-output /dev/null -o cutadapt.forwardTEtail.polyA.fastq cutadapt.forwardTEtail.fastq >> cutadapt.TE.log"]))

# reverse insertion tails
os.system(" ".join(['cutadapt','-j',parallel,'-e 0.1 -O 5 -a',"file:" + TreferenceRecom,"--action=trim --untrimmed-output /dev/null --rename='{id}_XO:Z:1_REV_{adapter_name}={match_sequence}'",clipFq1FileName,'-o','cutadapt.reverseTEtail.fastq.R1 >> cutadapt.TE.log']))
os.system(" ".join(['cutadapt','-j',parallel,'-e 0.1 -O 5 -a','file:'+TreferenceRecom,"--action=trim --untrimmed-output /dev/null --rename='{id}_XO:Z:2_REV_{adapter_name}={match_sequence}'",clipFq2FileName,'-o','cutadapt.reverseTEtail.fastq.R2 >> cutadapt.TE.log']))
os.system("cat cutadapt.reverseTEtail.fastq.R1 cutadapt.reverseTEtail.fastq.R2 > cutadapt.reverseTEtail.fastq")
os.system(" ".join(['cutadapt -j',parallel,"-e 0.1 -O 5 -a T{100}X --action=trim --rename='{id}_polyA={match_sequence}' --untrimmed-output /dev/null -o cutadapt.reverseTEtail.polyA.fastq cutadapt.reverseTEtail.fastq >> cutadapt.TE.log"]))

os.system("cat cutadapt.forwardTEtail.polyA.fastq cutadapt.reverseTEtail.polyA.fastq > cutadapt.TEtail.polyA.fastq")
os.system("rm -rf cutadapt.forwardTEtail.fastq.R1 cutadapt.forwardTEtail.fastq.R2 cutadapt.forwardTEtail.fastq cutadapt.forwardTEtail.polyA.fastq cutadapt.reverseTEtail.fastq.R1 cutadapt.reverseTEtail.fastq.R2 cutadapt.reverseTEtail.fastq cutadapt.reverseTEtail.polyA.fastq")

# forward insertion polyA
os.system(' '.join(["cutadapt",'-j',parallel,"-e 0.1 -O 10 -g XA{100} --action=trim --rename='{id}_XO:Z:1_FOR_polyA={match_sequence}' --untrimmed-output /dev/null -o cutadapt.forwardPolyA.fastq.R1",clipFq1FileName,"> cutadapt.polyA.log"]))
os.system(' '.join(['cutadapt','-j',parallel,"-e 0.1 -O 10 -g XA{100} --action=trim --rename='{id}_XO:Z:2_FOR_polyA={match_sequence}' --untrimmed-output /dev/null -o cutadapt.forwardPolyA.fastq.R2",clipFq2FileName, "> cutadapt.polyA.log"]))

# reverse insertion polyA
os.system(' '.join(["cutadapt",'-j',parallel,"-e 0.1 -O 10 -a T{100}X --action=trim --rename='{id}_XO:Z:1_REV_polyA={match_sequence}' --untrimmed-output /dev/null -o cutadapt.reversePolyA.fastq.R1",clipFq1FileName,">> cutadapt.polyA.log"]))
os.system(' '.join(['cutadapt','-j',parallel,"-e 0.1 -O 10 -a T{100}X --action=trim --rename='{id}_XO:Z:2_REV_polyA={match_sequence}' --untrimmed-output /dev/null -o cutadapt.reversePolyA.fastq.R2",clipFq2FileName,">> cutadapt.polyA.log"]))

os.system("cat cutadapt.forwardPolyA.fastq.R1 cutadapt.forwardPolyA.fastq.R2 cutadapt.reversePolyA.fastq.R1 cutadapt.reversePolyA.fastq.R2 > cutadapt.polyA.fastq")
os.system("rm -rf cutadapt.forwardPolyA.fastq.R1 cutadapt.forwardPolyA.fastq.R2 cutadapt.reversePolyA.fastq.R1 cutadapt.reversePolyA.fastq.R2")
os.system(' '.join(["SomTD_convertCutadapt.py","-t","cutadapt.TEtail.polyA.fastq","-a","cutadapt.polyA.fastq","-r",readLen, "-o",output]))


#########################################
############ detectInsertion ############
#########################################

os.system(' '.join(["SomTD_combineAlign.py","-g",output + ".genome.clip.pair.sam","-t",output+".te.clip.pri.sam","-G",output + ".genome.clip.sup.sam", "-T",output + ".te.clip.sup.sam","-c",output +".cutadapt.txt","-o",output]))
### deal with testAuto
os.system(' '.join(["SomTD_detectIns.py","-i",output + ".combine.clip.pair.sam","-o",output]))