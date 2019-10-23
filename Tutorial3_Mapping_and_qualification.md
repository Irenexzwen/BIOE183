# Tutorial3 Mapping and Quantification
In this tutorial you're gonna learn how to map the RNAseq data to a drosophila reference genome and quantify each gene's expression level. 

## 1.Introduction of Reads mapping
http://chagall.med.cornell.edu/RNASEQcourse/Slides_July2019_Day2.pdf

## 2.Different mapping strategy 
### 2.1 Alignment based
### 2.2 Alignment free
More info here: [Alignment-free sequence comparison: benefits, applications, and tools](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1319-7)

## 3.Get raw data and reference sequences ready
### 3.1 Raw data cleaning and preprocessing
We will use RNAseq data from [FlyAtlas2 database](http://flyatlas.gla.ac.uk/FlyAtlas2/index.html), which collects hundreds of RNAseq data of drosophila melanogaster. You can search by gene, category or tissue.

The raw RNAseq data can be found on the [resource page](https://github.com/Irenexzwen/BIOE183/blob/master/README.md) for the homework tutorials. If you'd like, you can analyze and clean the raw data using the techniques you learned from Tutorial 2.

Otherwise, we have conveniently provided the cleaned data for you here:
- [[female_head1_R1_clean]](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_raw/F_head1_clean_R1.fq),[[female_head1_R2_clean]](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_raw/F_head1_clean_R2.fq)
                   
- [[female_head2_R1_clean]](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_raw/F_head2_clean_R1.fq),[[female_head2_R2_clean]](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_raw/F_head2_clean_R2.fq)

- [[female_midgut1_R1_clean]](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_raw/F_midgut1_clean_R1.fq),[[female_midgut1_R2_clean]](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_raw/F_midgut1_clean_R2.fq)
                     
- [[female_midgut2_R1_clean]](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_raw/F_midgut2_clean_R1.fq),[[female_midgut2_R2_clean]](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_raw/F_midgut2_clean_R2.fq)

### 3.2 Reference sequences
You could download reference genome / transcriptome / gtf files of your familiar species from [ENSEMBLE](https://uswest.ensembl.org/info/data/ftp/index.html).
If you are analyzing Human or Mouse, you could also try Genecode.

```Shell
# Download drosophila genome reference sequence 
wget ftp://ftp.ensembl.org/pub/release-97/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa.gz
gzip -d Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa.gz  # decompress .gz file 

# Download drosophila transctiptome (cDNA) sequence
wget ftp://ftp.ensembl.org/pub/release-97/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.22.cdna.all.fa.gz
gzip -d Drosophila_melanogaster.BDGP6.22.cdna.all.fa.gz

# Download drosophila annotation gtf file
wget ftp://ftp.ensembl.org/pub/release-97/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.22.97.chr.gtf.gz
gzip -d Drosophila_melanogaster.BDGP6.22.97.chr.gtf.gz
```
## 4.Align to the genome and quantification
We will use STAR to align reads to the whole genome.
```Shell
conda install star

# check if STAR has been successfully installed.
STAR -h 

# build index
STAR --runThreadN 10 --runMode genomeGenerate \
--genomeDir /DS/reference/genome_STARidx \
--genomeFastaFiles /DS/reference/Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa \
--sjdbGTFfile /DS/reference/Drosophila_melanogaster.BDGP6.22.97.chr.gtf \
--sjdbOverhang 100 #reads length minus 1                                               
```
(Change the file paths to reflect where the respective files are located on your computer.)

`--runThreadN` Number of threads you use to run on your computer.  

`--genomeDir` The place you want to put your reference. (Make sure the location already exists; the command does not make a new folder.)

`--genomeFastaFiles` The genome fasta file we've just downloaded and uncompressed.  

`--sjdbGTFfile` Gtf file.  

`--sjdbOverhang` Usually equals read length minus 1.  

It might take some time to finish the alignment, and the total Memory usage peak at 10g during this process. If this memory requirement is beyond your computer, you could download the pre-computed index from our resource page. 

After we build the index, we're gonna map our reads towards the genome.
```Shell
for i in F_head1 F_head2 F_midgut1 F_midgut2
do
printf "STAR 
        --runThreadN 20 
        --genomeDir /DS/reference/genome_STARidx 
        --readFilesIn ${i}_clean_R1.fq ${i}_clean_R2.fq 
        --outSAMtype BAM SortedByCoordinate 
        --outFileNamePrefix ./${i}_STAR_genome" > ${i}_STAR.sh
        bash ${i}_STAR.sh 
done
```
This will give you two important results:
1) A sorted [bam file](https://support.illumina.com/help/BS_App_RNASeq_Alignment_OLH_1000000006112/Content/Source/Informatics/BAM-Format.htm) based on the coordinates.
2) A final \*Log.final.out file that includes mapping results information: total mapping ratio / unique mapping ratio / number of mapped reads etc. 
 
#### Quantify gene expression level using FeatureCount
Once you get the bam file (which records each reads align to which specific locations of the genome), you may want to summarize reads abundance for each gene.   
```Shell
for i in F_head1 F_head2 F_midgut1 F_midgut2
do
        printf "
        featureCounts -p -a                                 
        /DS/reference/Drosophila_melanogaster.BDGP6.22.97.chr.gtf -T 10 -o ${i}_count.txt 
        /DS/quant/STAR/${i}*.bam" > ${i}.sh
        bash ${i}.sh 
done

```
`-p` Check validity of paired-end distance.  

`-a` Name of an annotation file. GTF/GFF format by default.

To compare different featureCounts results into a whole, you could utilize some Unix techniques you've learned previously.
```Shell
paste <(less F_head1_count.txt|cut -f1,7-|sed 1d) <(less F_head2_count.txt|cut -f1,7-|sed 1d) <(less F_midgut1_count.txt|cut -f1,7-|sed 1d) <(less F_midgut2_count.txt|cut -f1,7-|sed 1d)|cut -f 1,2,4,6,8 > all_sample_count.txt
```

#### Normalization using bamCoverage
bamCoverage takes an alignment of reads or fragments as input (BAM file) and generates a coverage track (bigWig or bedGraph) as output. The coverage is calculated as the number of reads per bin, where bins are short consecutive counting windows of a defined size.

bamCoverage offers normalization by:
- RPKM = Reads Per Kilobase per Million mapped reads
- CPM = Counts Per Million mapped reads, same as CPM in RNA-seq
- BPM = Bins Per Million mapped reads, same as TPM in RNA-seq
- RPGC = reads per genomic content (1x normalization)

You may need to install the deeptools package which has the bamCoverage tool:
```
conda install -c bioconda deeptools
```

To run bamCoverage with BPM (TPM) normalization:
```
bamCoverage -b filename.bam -o filename.bw --normalizeUsing BPM
```

## 5.Align to the transcriptome and quantification using alignment-free method
We will use kallisto to do pseudoalignment. 
```Shell

# Install and check if kallisto is ready
conda install kallisto
kallisto -h

# build index
mkdir kallisto_index && cd kallisto_index
kallisto index -i drosophila_transcriptome.idx PATHTO/Drosophila_melanogaster.BDGP6.22.cdna.all.fa

# quantification
kallisto quant -i drosophila_transcriptome.idx -o /quant_kallisto/F_head1 PATH_TO_clean_R1.fq PATH_TO_clean_R2.fq
```

Now you've finished the kallisto alignment step, the results of the abundance of different genes are summarized in the file 
`abundance.tsv`. Be careful that the four samples's abundance file are of the same name. 

