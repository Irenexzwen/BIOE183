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
- super cool text processing softwares and text editor. 
- Great habitat for a lot of open source softwares ( especially the ones for bioinformatics ). 

# Prepare a linux working environment
## My system is windows:
For windows users, your options would be: 
1) Download one linux system distribution from windows store ( Ubuntu recommended ).
2) Log in a known server ( if your lab or other sources could provide you an account to an linux system ). 
3) Download a virtual machine and install an linux system. 

The first option would be recommended as the windows has it's own Windows Subsystem Language embedded. 
For the first 5min of the [video](https://www.youtube.com/watch?v=xzgwDbe7foQ), you will learn how to set up a linux distribution "Ubuntu" of your windows 10 system. 

Basically the steps include:
1) search "features" in windows10, open "Turn Windows features on or off".
2) On "Windows Features," tick the Windows Subsystem for Linux.
3) Choose "restart your computer right now" after you click ok. 
4) Download and install Ubuntu from the [windows store](https://www.microsoft.com/en-us/p/ubuntu/9nblggh4msv6?ocid=AID2000142_aff_7593_159229&tduid=%28ir__ka1wk3nrygkfr3k2kk0sohzixu2xg0e9pm2nab1600%29%287593%29%28159229%29%28%29%28UUwpUdUnU49661YYwYg%29&rtc=1&irgwc=1&irclickid=_ka1wk3nrygkfr3k2kk0sohzixu2xg0e9pm2nab1600&activetab=pivot:overviewtab).
5) Run ubuntu and setup your user name and password.

Once you finished that, you should see something like this:
<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/ubuntu.png">
Now you're officially logged into a linux system on your own computer! Next, we're going to learn more linux skills that will enable you to walk around your computer file system as easy as your familiar system. 


## My system is MacOS:


## Linux basics:
The file system structure of linux is hierarchical, everything is a file (document, .exe, script, .mp3 etc) at a specific layer of the system. Unlike windows, which store files on different drives (C: D: G: etc). The "/" is the root of the whole system, and any other files are stored in different levels of the tree. 
<img src="https://www.google.com/url?sa=i&source=images&cd=&cad=rja&uact=8&ved=2ahUKEwjrlcjDnNHkAhXCHjQIHY22B_0Qjhx6BAgBEAI&url=https%3A%2F%2Fwww.ques10.com%2Fp%2F17415%2Fexplain-file-system-hierarchyfhs-of-linux-1%2F&psig=AOvVaw1qhs83vWwVP7OA8-KPkSKv&ust=1568580928757395"> 
*Credit to: https://nepalisupport.wordpress.com/2016/06/29/linux-file-system-hierarchy/.*
 

#### 1. file system
The file system 
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
