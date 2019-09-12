# Tutorial 1 Preparation  
In this tutorial, you're gonna learn how to 
  1) get the working environment ready.
  2) Learn Linux basics.
  3) Learn R basics.
  4) Download and install softwares that will be used in the later sessions.
  
## Working environment

Most bioinformatics analysis are deal with super large datasets (over few gigabytes even hundreds gigabytes), which is impracticle for manual checking and modification. 
In this case, most bioinformatics person put linux/Unix in their daily toolbox. With a lot commands that especially designed for large text file,
you're able to sort / extract / subsititute / get the head / get the tail / search super fast without loading the whole file into your memory 
( as you will load the big file into your memory if it is opened in a windows / mac text reader ).

Beyond that, Linux OS is also 
- convinent and compatible with mutiple versions of the same software.
- easy to build up pipelines for batch files processing. 
- super cool text processing softwares beyond your imagine. 
- Great habitat for a lot of open source softwares ( especially the ones for bioinformatics ). 

# Prepare a linux working environment
## My system is windows:
For windows users, your options would be: 
1) Download one linux system distribution from windows store ( Ubuntu recommended ).
2) Log in a known server ( if your lab or other sources could provide you an account to an linux system ). 
3) Download Cygwin.
4) Download a virtual machine and install an linux system. 

## Linux basics:
3.2.1.1 Lab 1a
check the your present directory
pwd
check history
history
pipe history to grep to search for the cd command
history | grep cd
put history into a history.txt file
history > history.txt
make a directory called data
mkdir data
change into data directory
cd data
move history.txt file into data directory
mv ../history.txt ./
check manual page of wget command
man wget
redirect wget maunual page output into a file called wget.txt
man wget > wget.txt
return the lines that contain output in the wget.txt file
cat wget.txt | grep output
grep -i output wget.txt
Compress wget.txt file
gzip wget.txt
View Compressed file
cat wget.txt.qz
zcat wget.txt.qz
zcat wget.txt.qz | less

# R basics

# Software installation
http://www.bio-info-trainee.com/4030.html
https://www.jianshu.com/p/6e493a1e4240?utm_campaign=maleskine&utm_content=note&utm_medium=seo_notes&utm_source=recommendation
conda install -c bioconda samtools=1.5
conda install -c bioconda htseq=0.7.2
conda install -c bioconda hisat2=2.0.5
conda install -c bioconda fastqc=0.11.5
conda install -c jfear sratoolkit=2.8.1
conda install star
