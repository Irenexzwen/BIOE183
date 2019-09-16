# Tutorial3 Mapping and Quantification
In this tutorial you're gonna learn how to map the RNAseq data to a zebrafish reference genome and quantify each gene's expression level. 

## 1.Introduction of Reads mapping
http://chagall.med.cornell.edu/RNASEQcourse/Slides_July2019_Day2.pdf

## 2.Different mapping strategy 
### 2.1 Alignment based
### 2.2 Alignment free

## 3.Get reference sequences ready
You could download reference genome / transcriptome / gtf files of your familiar species from [ENSEMBLE](https://uswest.ensembl.org/info/data/ftp/index.html).
If you are analyzing Human or Mouse, you could also try Genecode.

```Shell
# Download zebrafish genome reference sequence 
wget ftp://ftp.ensembl.org/pub/release-97/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna_rm.primary_assembly.fa.gz
gzip -d Danio_rerio.GRCz11.dna_rm.primary_assembly.fa.gz  # decompress .gz file 

# Download zebrafish transctiptome (cDNA) sequence
wget ftp://ftp.ensembl.org/pub/release-97/fasta/danio_rerio/cdna/Danio_rerio.GRCz11.cdna.all.fa.gz
gzip -d Danio_rerio.GRCz11.dna_rm.primary_assembly.fa.gz

# Download zebrafish annotation gtf file
ftp://ftp.ensembl.org/pub/release-97/gtf/danio_rerio/Danio_rerio.GRCz11.97.chr.gtf.gz
gzip -d Danio_rerio.GRCz11.97.chr.gtf.gz
```
## 4.Align to the genome and quatification
We will use STAR to align reads to the whole genome.
```Shell
conda install star

# check if STAR has been successfully installed.
STAR -h 

# build index
STAR --runThreadN 3 --runMode genomeGenerate \
--genomeDir /mnt/d/UCSD/RNAseq/STAR_index \                    
--genomeFastaFiles /mnt/d/UCSD/RNAseq/reference/genome/Danio_rerio.GRCz11.dna_rm.primary_assembly.fa \  
--sjdbGTFfile /mnt/d/UCSD/RNAseq/reference/gtf/Danio_rerio.GRCz11.97.chr.gtf \
--sjdbOverhang 74                                                
```

`--runThreadN` Threads you use to run on your computer.
`--genomeDir` The place you want to put your reference.
`--genomeFastaFiles` The genome fasta file we're just downloaded and uncompressed.
`--sjdbGTFfile` Gtf file.
`--sjdbOverhang` Usually equals read length minus 1.

It takes X min to finish on my computer. The total Memory usage peak at 10g during this process. If this memory requirement is beyond your computer, you could download the pre-computed index from [here]()


## 5.Align to the trasctriptome and quantification
We will use kallisto to do pseudoalignment. 
```Shell

# Install and check if kallisto is ready
conda install kallisto
kallisto -h

# build index
mkdir kallisto_index && cd kallisto_index
kallisto index -i zebrafish_transcriptome.idx Danio_rerio.GRCz11.cdna.all.fa.gz

# quantification
kallisto quant -i /mnt/d/UCSD/RNAseq/reference/transcriptome/zebrafish_transcriptome.idx -o /mnt/d/UCSD/RNAseq/quant/2cell /mnt/d/UCSD/RNAseq/fastp/2cells_R1_clean.fastq /mnt/d/UCSD/RNAseq/fastp/2cells_R2_clean.fastq 

kallisto quant -i /mnt/d/UCSD/RNAseq/reference/transcriptome/zebrafish_transcriptome.idx -o /mnt/d/UCSD/RNAseq/quant/6h /mnt/d/UCSD/RNAseq/fastp/6h_R1_clean.fastq /mnt/d/UCSD/RNAseq/fastp/6h_R2_clean.fastq

# 
```



