# Tutorial2 Understand the Raw Data

In this tutorial you're gonna learn:
1. Experienmentally:
    1) The raw data for RNAseq data.
    2) Obtaining the raw data.
2. Bioinformatics:
    3) Raw data quality check.
    4) Raw data cleaning and preprocessing. 
    

## 1). RNAseq data
The RNA sequencing technology was designed to answer one of the fundamental questions in biology - 
which genomic loci are expressed at a specific time and what's the expression level of it? 

## 2). Obtaining the raw data 
NGS sequencing technology introduction, Illustration.
Where to download,
Database,
GEO.

## 3). Raw data quality check
In this project we're gonna download four fastq file from EBI RNA-seq trainning bootcamp. The information about the data is [here](ftp://ftp.ebi.ac.uk/pub/training/Train_online/RNA-seq_exercise/zebrafish-rna-seq.pdf). Basically it's a pairend RNAseq data of zebrafish at two time points. Download the data into any folder you'd like. 
```Shell
mkdir /mnt/d/UCSD/RNAseq   # I will put the data in this path
wget ftp://ftp.ebi.ac.uk/pub/training/Train_online/RNA-seq_exercise/2cells_1.fastq  
wget ftp://ftp.ebi.ac.uk/pub/training/Train_online/RNA-seq_exercise/2cells_2.fastq
wget ftp://ftp.ebi.ac.uk/pub/training/Train_online/RNA-seq_exercise/6h_1.fastq
wget ftp://ftp.ebi.ac.uk/pub/training/Train_online/RNA-seq_exercise/6h_2.fastq

ls  # check all files are shown in your folder 
```
The RNAseq raw reads file are all stored in [FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format) files, which is the most commonly used format to store reads information (another format is FASTA, you could learn [here](https://bioinformatics.stackexchange.com/questions/14/what-is-the-difference-between-fasta-fastq-and-sam-file-formats)). It includes:
1. @ followed by the read ID and possibly sequencing information
2. sequenced bases ATCG...
3. + (perhaps followed by the read ID again, or some other description)
4. quality scores for each base of the sequence, should be of the same length as 2.

<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/fastq.png">
For pairend sequencing data, there are always two fastq files for the same sample, one is forward and the other is the backward. They should have the same number of reads, and usually the same file size. (With this property you could do a quick check at your data)

### Raw Reads quality check with FastQC
First make sure you've downloaded fastqc, you could check using type `fastqc -h` in the console, and if successfully installed you should see a manual page of the tool. 

To run fastqc is really simple:
```Shell
cd /mnt/d/UCSD/RNAseq
fastqc -h                                                               # get the help page of fastqc

# run fastqc and the output file will be generated in the same folder 
fastqc 2cells_1.fastq  2cells_2.fastq  6h_1.fastq  6h_2.fastq           
```
:question:What if you want to put the result in other folders? Use `fastqc -h` to find the anwser!

### Understand Fastqc results:
After you successfully run the code above, you should get one `.html` file for each fastq file with the same sample name. Open these `html` file use a web browser on your computer (If you don't know where to find the files, please review the file system mapping we've talked about in tutorial 1.). In general, we only care about few aspects of the whole report, I will use 2cells_2_fastqc.html as an example.

#### 1)Per Base sequence quality
Show the plot of quality of each base for all reads (Red: poor; Orange: reasonable; Green: good), this is a good per base quality graph because almost all base quality is above orange area. 


<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/fastqc_per_base_seq_quality.png">

#### 2)Per sequence quality score
Show the number of mean quality for all reads. What to look for: The distribution of average read quality should be fairly tight in the	
upper range of the plot.

<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/baseQ.png">

#### 3)Per base sequence content
Show the percentage of A/C/T/G for each position of all reads. What to look for: For RNA library first 10~15 bases, there will be non-uniform distribution of ATCG because of the library kit used. However, for DNA sequencing, usually all nucleotides should have similar percentage across all positions.

<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/basecontent.png">

#### 4) Per sequence GC content
Show the distribution of GC content for all reads. What to look for: whether the distribution is a normal distribution with the peak pf the curve at the mean GC content for the organism be sequenced. 
<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/GC.png">

#### 5) Sequence duplication levels
This plot means percentage of reads of asgiven sequence in the file which are present a given number of times in the file. There are two sources of duplication:
- Biased PCR enrichment (some sequence might be favored by PCR)
- The overrepresented sequence in your library <br>

Your RNA sequencing library will contain highly abundant transcripts. It is expected that duplicated reads will be observed for high abundance transcripts. The RNA-Seq data below was flagged as Failed by FastQC even though the duplication is expected in this case.

<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/dup.png">

#### 6) Overrepresented sequences:
fastqc will look at the first 50 bp of all sequences, see if this 50-mer is overrepresented compare to an expected ratio of all reads, eg 0.1%. The table will give you useful information, even includes possible sources. 
<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/over.png">

:star2:Your tasks: Understand what other plots of the report indicates. If there is a problem please refer to our [reference manual](http://chagall.med.cornell.edu/RNASEQcourse/FastQC_TutorialAndFAQ_080717.pdf). Understand why there are some difference in fastqc results with respect to DNA library and RNA library. 


## 4). Preprocessing of the raw data 
Our next step is to clean the raw data. Why we need to clean our data, and what are those unwanted contaminats in RNAseq data? 
In essence, we want our RNAseq library to reflect a faithful representation of the real situation in a cell. Which is to say, we're trying to aviod these signals for the downstream analysis:
- artificial DNA, adapters, primers.
- contamination from other species, like bateria.
- rRNA usually counted to a very large percent of all transcriptome, which might suppress the library amplication efficacy. Remove them use kit while doing experiement or remove bioinformatically after you get the raw reads.

Here we will learn how to remove the adapters of the raw reads. Retaining adapters might larged influence the reads alignment step (too many bases would be identified as mismatch as the adapters are foreign sequence to the organism genome). 

```Shell
# We will use fastp to automatically detect and remove adpaters. 

# check if successfully installed fastp, if not please use conda install fastp
fastp -h  

# run fastp
mkdir fastp  # to store fastp results

fastp -i 2cells_1.fastq -I 2cells_2.fastq -o fastp/2cells_R1_clean.fastq -O fastp/2cells_R2_clean.fastq -h 2cell_fastp.html
fastp -i 6h_1.fastq -I 6h_2.fastq -o fastp/6h_R1_clean.fastq -O fastp/6h_R2_clean.fastq -h 6h_fastp.html
```

It will finish with few seconds and the results would print out to the screen:
```Shell
Read1 before filtering:
total reads: 786742
total bases: 59792392
Q20 bases: 57669779(96.45%)
Q30 bases: 43262473(72.3545%)

Read2 before filtering:
total reads: 786742
total bases: 59792392
Q20 bases: 57447802(96.0788%)
Q30 bases: 40784290(68.2098%)

Read1 after filtering:
total reads: 770487
total bases: 58522540
Q20 bases: 56667730(96.8306%)
Q30 bases: 42631358(72.846%)

Read2 aftering filtering:
total reads: 770487
total bases: 58522540
Q20 bases: 56891294(97.2126%)
Q30 bases: 40538007(69.269%)

Filtering result:
reads passed filter: 1540974
reads failed due to low quality: 31554
reads failed due to too many N: 956
reads failed due to too short: 0
reads with adapter trimmed: 2790
bases trimmed due to adapters: 69036

Duplication rate: 11.508%

Insert size peak (evaluated by paired-end reads): 118

JSON report: fastp.json
HTML report: fastp.html

fastp -i 2cells_1.fastq -I 2cells_2.fastq -o fastp/2cells_R1_clean.fastq -O fastp/2cells_R2_clean.fastq
fastp v0.20.0, time used: 6 seconds
```

You could check the result in 2cell_fastp.html and 6h_fastp.html. This is very detailed and pretty much like what you've learned in fastqc report. 

Now you also got two clean pairend fastq files! We're ready for the next step. 
```Shell
fastp/2cells_R1_clean.fastq  
fastp/2cells_R2_clean.fastq  
fastp/6h_R1_clean.fastq
fastp/6h_R2_clean.fastq
```

