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


It might take some time to finish the alignment, and the total Memory usage peak at 10g during this process. If this memory requirement is beyond your computer, you could download the pre-computed index from our resouce page. 


After we build the index, we're gonna map our reads towards the genome.
```Shell
for i in F_head1 F_head2 F_midgut1 F_midgut2
do
printf "STAR 
        --runThreadN 20 
        --genomeDir /DS/reference/genome_STARidx 
        --readFilesIn /${i}_clean_R1.fq /${i}*/${i}_clean_R2.fq 
        --outSAMtype BAM SortedByCoordinate 
        --outFileNamePrefix ./${i}_STAR_genome" > ${i}_STAR.sh
        bash ${i}_STAR.sh 
done
```
This will give you two important results:
1) A sorted [bam file](https://support.illumina.com/help/BS_App_RNASeq_Alignment_OLH_1000000006112/Content/Source/Informatics/BAM-Format.html) based on the coordinates.
2) A final \*Log.final.out file that includes mapping results information: total mapping ratio / unique mapping ratio / number of mapped reads etc. 
 
#### Quantify gene expression level use FeatureCount
Once you get the bam file (which records each reads align to which specific locations of the genome), you may want to summarize reads abundance for each gene.   
```Shell
for i in F_head1 F_head2 F_midgut1 F_midgut2
do
        printf "
        /software/featureCounts/subread-1.6.4-Linux-x86_64/bin/featureCounts -p -a                                /dataOS/wenxingzhao/class/BIOE183/DS/reference/Drosophila_melanogaster.BDGP6.22.97.chr.gtf -T 10 -o ${i}_count.txt /dataOS/wenxingzhao/class/BIOE183/DS/quant/STAR/${i}*.bam" > ${i}.sh
        bash ${i}.sh 
done

```


## 5.Align to the trasctriptome and quantification using alignment-free method
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
```

Now you've finished the kallisto alignment step, the results of the abundance of different genes are summarized in the file 
`abundance.tsv`. Be careful that the two samples's abundance file are of the same name. 

Next we will use the gene expression table (from kallisto) of the two sample to do differential expression analysis. 




