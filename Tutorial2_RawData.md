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
After you successfully run the code above, you should get one `.html` file for each fastq file with the same sample name. Open these `html` file use a web browser on your computer (If you don't know where to find the files, please review the file system mapping we've talked about in tutorial 1.). In general, we only care about I will use 2cells_2_fastqc.html as an example.

#### 1)Per Base sequence quality
Show the plot of quality of each base for all reads (Red: poor; Orange: reasonable; Green: good)
<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/fastqc_per_base_seq_quality.png">

#### 2)


## 4). Preprocessing of the raw data f


