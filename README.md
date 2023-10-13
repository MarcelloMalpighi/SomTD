# SomTD  
Retrotransposons contribute approximately 40 % of the human genome and subfamilies of ALU, LINE1 and SVA elements remain actively mobile. Nonetheless, the detection of de novo transposable element (TE) insertions poses a significant challenge due to chimera artifacts. Additionally, evaluating the insertion rate with bulk or single-cell sequencing data presents certain challenges.  
SomTD is a tool designed for detecting de novo TE insertions and evaluating the insertion rate utilizing both traditional, rule-based algorithms and a convolutional neural network (CNN) model. It is applicable for both bulk and single-cell sequencing data. Notable features of SomTD are as follows:  
- SomTD prioritizes split read pairs, subsequently extracting discordant read pairs as required. The supplementary alignments of the split read pairs are utilized for pinpointing the location because TE part of the clipped read may be identified as the primary alignment, potentially resulting in a misleading insertion location around a reference insertion (fig1).  
- A lightweight CNN, is applied to extract every suitable read pair, contrasting previous machine/deep learning applications focusing on the insertion level for insertion detection or genotyping (fig2). This approach allows for the detection of weak signals of rare insertions by minimizing data loss and distinguishing chimera artifacts via features difficult to discern.  
- SomTD estimates the insertion rate based on the cumulative sum of variant allele fraction including rare insertions in bulk sequencing data based on its elaborated sensitivity and accuracy, which remains comparable across bulk and single-cell sequencing data.  
<img src="https://github.com/MarcelloMalpighi/SomTD/blob/main/SomTD_fig1.png" height="197px" width="320px"/>  
<img src="https://github.com/MarcelloMalpighi/SomTD/blob/main/SomTD_fig2.png" height="292px" width="587px"/>  

## Dependencies  
1. bedtools  
2. bwa  
3. cutadapt  
4. pysam  
5. pytorch  
6. samblaster  
7. samtools  
## Run SomTD  
```
SomTD
Usage: SomTD.py [options]

Options:
  -h, --help            show this help message and exit
  --input1=INPUT1       input file, bam/sam/fastq, please use .bam .sam
                        .fq/.fastq as a suffix, mandatory
  --input2=INPUT2       input file, if fastq, please use .fq/.fastq as a
                        suffix, optional
  -o OUTPUTPATH, --outputPath=OUTPUTPATH
                        output path, directory name will be output name,
                        directory will be generated autonomously if not exist,
                        mandatory
  -c CUTOFF, --cutoff=CUTOFF
                        minimum soft-clipped length, limit itself is included,
                        optional, default: 10
  -f FRAGLEN, --frag=FRAGLEN
                        expected fragment length, mandatory
  --std=FRAGSTD         standard deviation of fragment length, mandatory
  -r READLEN, --readLen=READLEN
                        read length, mandatory
  -g GREFERENCE, --Greference=GREFERENCE
                        genome reference sequence, fa, mandatory
  -t TREFERENCE, --Treference=TREFERENCE
                        transposon reference sequence, fa, mandatory
  --TreferenceRecom=TREFERENCERECOM
                        reversed complemented transposon reference sequence,
                        fa, mandatory
  -G GINDEX, --Gindex=GINDEX
                        genome reference sequence bwa index, mandatory
  -T TINDEX, --Tindex=TINDEX
                        transposon reference sequence bwa index, mandatory
  -p PARALLEL, --parallel=PARALLEL
                        number of threads, optional, default: 2
  -m MEMORY, --memory=MEMORY
                        memory per thread used, optional, defalut: 2G
```
