# BENG183

## Homework notes:
1) [Working environment setup](https://github.com/Irenexzwen/BIOE183/blob/master/Tutorial1_Preparation.md)
2) [Raw Data QC and Cleaning](https://github.com/Irenexzwen/BIOE183/blob/master/Tutorial2_RawData.md)
3) [Mapping and quantification](https://github.com/Irenexzwen/BIOE183/blob/master/Tutorial3_Mapping_and_qualification.md)
4) [Differential analysis](https://github.com/Irenexzwen/BIOE183/blob/master/Tutorial4_DE.md)

## Raw Sample files:
### 1)Sample for Raw Data cleaning illustration
We will use a dataset from a previous [EBI training courses](https://www.ebi.ac.uk/training/online/course/ebi-next-generation-sequencing-practical-course/rna-sequencing/rna-seq-analysis-transcriptome). This data is derived from sequencing of mRNA from zebrafish embryos in two diﬀerent developmental stages. Sequencing was performed on the Illumina platform
and generated 76bp paired-end sequence data using poly-(A)+ selected RNA. 

Download data from ftp://ftp.ebi.ac.uk/pub/training/Train_online/RNA-seq_exercise/
- [2cells_1.fastq](ftp://ftp.ebi.ac.uk/pub/training/Train_online/RNA-seq_exercise/2cells_1.fastq) 
- [2cells_2.fastq](ftp://ftp.ebi.ac.uk/pub/training/Train_online/RNA-seq_exercise/2cells_2.fastq)
- [6h_post_fertilisation_R1.fastq](ftp://ftp.ebi.ac.uk/pub/training/Train_online/RNA-seq_exercise/6h_1.fastq)
- [6h_post_fertilisation_R2.fastq](ftp://ftp.ebi.ac.uk/pub/training/Train_online/RNA-seq_exercise/6h_2.fastq)

### 2)Sample for the homework
We will use RNAseq data from [FlyAtlas2 database](http://flyatlas.gla.ac.uk/FlyAtlas2/index.html), which collects hundreds of RNAseq data of drosophila melanogaster. You can search by gene, category or tissue. Here we downloaded 4 samples (female_head x 2, female_midgut x 2).
- [[female_head1_R1_raw]](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_raw/2089N0008_GCCAAT_L001_R1_001.fastq.gz),[[female_head1_R2_raw]](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_raw/2089N0008_GCCAAT_L001_R2_001.fastq.gz)

- [[female_head2_R1_raw]](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_raw/2089N0007_ACAGTG_L001_R1_001.fastq.gz), [[female_head2_R2_raw]](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_raw/2089N0007_ACAGTG_L001_R2_001.fastq.gz )

- [[female_midgut1_R1_raw]](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_raw/160825_D00261_0358_BC9JAVANXX_1_TP-D7-009_TP-D5-001_1.fastq.gz),[[female_midgut1_R2_raw]](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_raw/160825_D00261_0358_BC9JAVANXX_1_TP-D7-009_TP-D5-001_2.fastq.gz )

- [[female_midgut2_R1_raw]](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_raw/160825_D00261_0358_BC9JAVANXX_1_TP-D7-010_TP-D5-008_1.fastq.gz),[[female_midgut2_R2_raw]](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_raw/160825_D00261_0358_BC9JAVANXX_1_TP-D7-010_TP-D5-008_2.fastq.gz)

## Clean Data
- [[female_head1_R1_clean]](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_raw/F_head1_clean_R1.fq),[[female_head1_R2_clean]](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_raw/F_head1_clean_R2.fq)
                   
- [[female_head2_R1_clean]](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_raw/F_head2_clean_R1.fq),[[female_head2_R2_clean]](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_raw/F_head2_clean_R2.fq)

- [[female_midgut1_R1_clean]](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_raw/F_midgut1_clean_R1.fq),[[female_midgut1_R2_clean]](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_raw/F_midgut1_clean_R2.fq)
                     
- [[female_midgut2_R1_clean]](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_raw/F_midgut2_clean_R1.fq),[[female_midgut2_R2_clean]](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_raw/F_midgut2_clean_R2.fq)

## Reference genome.fa / transcriptome.fa / gtf
We usually download the reference data from [ensemble](https://uswest.ensembl.org/info/data/ftp/index.html). You search "drosophila" and choose DNA / cDNA / gtf, then you use a `wget` to download.
- drosophila genome: 
ftp://ftp.ensembl.org/pub/release-97/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa.gz

- drosophila transcriptome: 
ftp://ftp.ensembl.org/pub/release-97/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.22.cdna.all.fa.gz

- drosophila gtf: 
ftp://ftp.ensembl.org/pub/release-97/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.22.97.chr.gtf.gz

## Index files
Pre-computed index files: [download here](https://drive.google.com/open?id=1CT-iZ2PzPwWa44KxhjYsv4pD5GnnJStg)

## Mapping .bam file
- [female_head1_bam](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_quant/F_head1_STAR_genomeAligned.sortedByCoord.out.bam)
- [female_head2_bam](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_quant/F_head2_STAR_genomeAligned.sortedByCoord.out.bam)
- [female_midgut1_bam](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_quant/F_midgut1_STAR_genomeAligned.sortedByCoord.out.bam)
- [female_midgut2_bam](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_quant/F_midgut2_STAR_genomeAligned.sortedByCoord.out.bam)

## FeatureCounts Count table
- [female_head1_count](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_quant/F_head1_count.txt)
- [female_head2_count](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_quant/F_head2_count.txt)
- [female_midgut1_count](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_quant/F_midgut1_count.txt)
- [female_midgut2_count](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_quant/F_midgut2_count.txt)

- [all sample](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_quant/all_sample_count.txt)

## Recourses:
#### 1.Tutorials:
[1] Weill Cornell Medical Colledge: http://chagall.med.cornell.edu/RNASEQcourse/

#### 2.Software manuals:
- [Bioconda](https://www.youtube.com/watch?v=lGa9PCSH5IU) starting from 3:30. 
- [fastQC manual](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf)
- [fastp manual](https://github.com/OpenGene/fastp)

