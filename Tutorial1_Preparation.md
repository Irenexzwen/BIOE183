# Tutorial 1 Preparation  
In this tutorial, you will learn to 
  1) set up Linux environment on your laptop (for Windows users),
  2) learn Linux basics,
  3) download and install softwares useful for genomic data analysis and required for this course.
  
## Working environment

Most bioinformatics analysis deals with very large datasets (from hundreds of GBs to several TB), which is impractical for manual checking and modification. 
In this case, most bioinformatics person put linux/Unix in their daily toolbox. With a lot commands that especially designed for large text file,
you're able to sort / extract / subsititute / get the head / get the tail / search super fast without loading the whole file into your memory 
( as you will load the big file into your memory if it is opened in a windows / mac text reader ).

Beyond that, Linux OS is also 
- convinent and compatible with mutiple versions of the same software.
- easy to build up pipelines for batch files processing. 
- super cool text processing softwares and text editor. 
- Great habitat for a lot of open source softwares ( especially the ones for bioinformatics ). 

# Prepare a linux working environment
## My system is MacOS:
You can directly type Linux commands on your computer. Go to Linux basics.

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

# Linux basics:
The file system structure of linux is hierarchical, everything is a file (document, .exe, script, .mp3 etc) at a specific layer of the system. Unlike windows, which store files on different drives (C: D: G: etc). The "/" is the root of the whole system, and any other files are stored in different levels of the tree. Important folders have been listed below, about their naming tradition one can refer to [here](http://www.linuxstories.net/linux-directory-structure-file-system-structure/).
<img src="https://github.com/WGLab/dragonstar2019/blob/master/day1_linux/img/directory.png"> 

### 1) directory navigation
Home directory:
Everytime you open "Ubuntu" you will see the shell prompt like:

`user_name @ servername: current_directory $ COMMAND_TYPE_HERE`.

Every new user will have their own home directory: `/home/user_name`, your home directory are often short for `~`,which you could check it here:
<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/linux.png">

Relative path and absolute path:
```
## The following characters have special meanings in linux system.
.     # current work directory
..    # parent directory
-     # previous work directory
~     # home directory
```

Navigate between directories:
```Shell
cd \mnt\c       # change directory to \mnt\c
cd ..           # change to parent directory
cd .            # stay where you are 
pwd             # Check the absolute path of current directory
mkdir A         # Make a new directory at the current path called A.
rmdir -r A      # remove a A at the current path. (here -r means, remove A in an recursively fashion.)
```

### 2) Check files in current folder:
```Shell
ls            # list files in current folder, not include hidden files
ls -a         # list files in current folder, include hidden files
ls -t         # list files in current folder, sorted by time
man ls        # check all the features of tool `ls`
```

### 3) Copy, rename and delete a file or folder:
```Shell
cp A B     # copy file A to file B in current folder, you could also change B to a absolute path
rm A       # remove file A 
rm *.txt   # remove all txt files (* is a wildcard character)
rm -rf A   # remove directory A without inquiry ( could be dangerous some time )
mv A B     # change the name of file(folder) A to B 
mv A C     # move A into folder C ( if C is a folder) 
```

### 4) Check text file and modification:
```Shell
less m.txt              # check file m.txt, scroll down a file using ↑ and ↓ (most commonly used)
cat m.txt               # open file m.txt from the beginning to end
tac m.txt               # open file m.txt from the end to beginning
head -n 10 m.txt        # only open the first 10 lines of m.txt
tail -n 10 m.txt        # only open the last 10 lines of m.txt
```
To learn more about linux basic operations, you could check these resources:
- bash command notebook from UCSD people! Thanks Jessica Zhou and Bill Greenwald (be sure to use jaccob engineering email to login): [here](https://colab.research.google.com/drive/1RRxtCSvTgXnLZWnfXvHpVu8G1-yR3-bl)
- Linux tutorial from MIT: [here](http://math.mit.edu/services/help/new/unix.php).
- Unix/Linux tutorial from Berkeley: [here](https://people.ischool.berkeley.edu/~kevin/unix-tutorial/toc.html).

*Corresponding between windows file path and Ubuntu file path system.*

```Shell
\mnt\c   <-----> C:\
\mnt\d   <-----> D:\
```


## Software installation with miniconda
Next we're going to download and install some softwares for bioinformatics analysis, especially for this project - RNAseq analysis.
Most people use [Anaconda](https://en.wikipedia.org/wiki/Anaconda_(Python_distribution)) as a package management tool for python. However, it's power is much beyond that, here we're gonna use miniconda (lightweight version of conda without preinstalled some popular packages).

### Download installation script on your favorite place on your computer:
For example, I'm gonna download and install it on D:\Miniconda on my own computer:
```Shell
$ cd /mnt/d
$ mkdir Miniconda
$ cd Miniconda
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh  # for Windows users
$ bash Miniconda3-latest-Linux-x86_64.sh # for Windows users 
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh  # for Mac users
$ bash Miniconda3-latest-MacOSX-x86_64.sh # for Mac users 

```
`wget` download the script into the folder, `bash x.sh` is excute x.sh file. During the installation you are required to follow the instructions prompted from Miniconda.  

If you install conda with default setting, the excutable conda usually locates in somewhere else, for example:
<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/conda.png">

Now we need to add that path into our environment default path so that you could use `conda` anywhere.
```Shell
export PATH="/path/to/your/miniconda3/bin/:$PATH"
conda   # running conda, now you should see the conda is successfully installed.
```

### Use bioconda to install softwares we want:
[Bioconda](https://daler.github.io/bioconda-docs/) is a channel for the conda package manager specializing in bioinformatics softwares. A full list could be found [here](https://anaconda.org/bioconda/repo). It also provides convient way for version checking on different softwares as well as virtual environment for each individual project or different analysis pipelines. 

Now, we're going to install softwares that will be useful for RNAseq analysis. 

#### 1)Step1, set up a virtual environment specially for this RNAseq analysis project:
```Shell
conda create -n rnaseq python=3         # create an environment called "rnaseq", and install python3 .
conda init bash 
```
Close your current bash and restart ubuntu. Now you could see `(base)` has been had at the beginning of the orginal line. Then you switch to RNAseq environment and install packages only in this environment.

```Shell
conda activate rnaseq                   # activate a new environment called rnaseq
conda deactivate                 		# deactivate 
```

Add bioconda source:
```Shell
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
```

More options for conda:
```Shell
conda env list         # check all virtual environment on your machine
conda list             # list out all the packages installed in this conda environment
conda update conda     # update conda
conda update A         # update package A
conda search A         # check whether package A is in the source channel.

conda env remove -n ENV_NAME    # remove a virtual environment.
```

Install packages: 
```Shell
conda install fastp=0.20.0 
conda install fastqc=0.11.8
conda install subread=1.6.4
conda install samtools=1.9
conda install htseq=0.11.2
conda install hisat2=2.1.0
conda install star=2.7.2b
···
