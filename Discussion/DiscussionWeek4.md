# Discussion session Week 4
Mapping reads with STAR and visualization of mapped data with IGV

### 1. Download required data/software
1) IGV browser (https://software.broadinstitute.org/software/igv/download)

2) reference genome of Drosophila
```bash
wget http://homer.ucsd.edu/zeyang/BENG183/Drosophila_melanogaster.BDGP6.22.dna.toplevel.2R.fa
```

3) annotations of genes (information of where exons and introns are so that reads can be mapped even if there is splicing)
```bash
wget ftp://ftp.ensembl.org/pub/release-97/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.22.97.chr.gtf.gz
gzip -d Drosophila_melanogaster.BDGP6.22.97.chr.gtf.gz
```

4) fastq files
```bash
wget http://homer.ucsd.edu/zeyang/BENG183/lib_002_mapped.1.fastq
wget http://homer.ucsd.edu/zeyang/BENG183/lib_002_mapped.2.fastq
```

### 2. Build STAR index for fast mapping
```bash
mkdir ./Drosophila_2R_star/
STAR --runMode genomeGenerate --genomeDir ./Drosophila_2R_star/ --genomeFastaFiles Drosophila_melanogaster.BDGP6.22.dna.toplevel.2R.fa --sjdbGTFfile Drosophila_melanogaster.BDGP6.22.97.chr.gtf --sjdbOverhang 99 --genomeSAindexNbases 11
```

### 3. Mapping reads with STAR
```bash
STAR --runThreadN 1 --genomeDir ./Drosophila_2R_star/ --readFilesIn lib_002_mapped.1.fastq lib_002_mapped.2.fastq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix Drosophila_mapped_
```

### 4. Visualize mapped data in IGV
Generate index file for BAM file before loading into IGV
```bash
samtools index Drosophila_mapped_Aligned.sortedByCoord.out.bam
```

### 5. Quantify gene expression
```bash
featureCounts -p -a Drosophila_melanogaster.BDGP6.22.97.chr.gtf -o counts.txt Drosophila_mapped_Aligned.sortedByCoord.out.bam
```

