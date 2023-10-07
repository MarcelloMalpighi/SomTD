# SomTD  
Retrotransposons contribute approximately 40 % of the human genome and subfamilies of ALU, LINE1 and SVA elements remain actively mobile. Nonetheless, the detection of de novo transposable element (TE) insertions poses a significant challenge due to the interference of chimera artifacts. Additionally, evaluating the insertion rate with bulk or single-cell sequencing data presents certain challenges.  
SomTD is a tool designed for detecting de novo TE insertions and evaluating the insertion rate utilizing both traditional, rule-based algorithms and a convolutional neural network (CNN) model. It is applicable for both bulk and single-cell sequencing data. Notable features of SomTD are as follows:  
- SomTD prioritizes split read pairs, subsequently extracting discordant read pairs as required. The supplementary alignments of the split read pairs are utilized for pinpointing the location because TE part of the clipped read may be identified as the primary alignment, potentially resulting in a misleading insertion location around a reference insertion.  
- A lightweight CNN, is applied at the read pair level to extract every suitable read pair, contrasting previous machine or deep learning applications that focused on the insertion level for detection or genotyping of insertions. It is capable of distinguishing chimera artifacts via features that are difficult to discern through conventionally designed rules. Other tools must discard some read pairs relating to a potential insertion to generate a fixed-shape graph at the insertion level for CNN, which results in data loss. By contrast, SomTD applies deep learning at the read pair level, ensuring input shape consistency, and eliminating data loss. The attributes of CNN are also well-suited for managing the multi-channel information furnished by a read pair.  
- SomTD introduces a new metric to assess transposition burden, which remains comparable across bulk and single-cell sequencing data. The enhanced sensitivity of CNN at the read pair level allows for precise calculation of this metric in bulk sequencing data.  
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
  -d DISTREF, --distRef=DISTREF
                        minimum distance from reference insertions,
                        nonreference insertions with a shorter one will be
                        ignored, optional, default: 0
  -s, --selfFilter      filter insertions shared by all input samples,
                        optional, default: disabled
```
