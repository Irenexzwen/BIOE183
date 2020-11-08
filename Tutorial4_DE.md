# Tutorial4 Analysis of differentially expressed genes

## Preparations
### R basics
If you would like to brush up R basics, find the following resources:
[R basics](https://www.quora.com/What-are-some-good-resources-for-learning-R-1)
[- R bioconductor](https://www.coursera.org/learn/bioconductor)
[- previous discussion session](https://github.com/Irenexzwen/BIOE183/blob/master/Discussion/DiscussionTutorial_ClusterAnalysis.md#2-r-basics).

R is a language and environment for statistical computing and graphics. It is wildly used among bioinformaticians not only because it's pretty easy to work with table data, but also because bioconductor (a pool of packages target for high-throughput biological assays) was released on R platform.

RStudio is recommended for this tutorial even though it is completely doable to run all commands in R. To install RStudio, please refer to [the previous discussion notes](https://github.com/Irenexzwen/BIOE183/blob/master/Discussion/DiscussionTutorial_ClusterAnalysis.md#1-prepare-the-r-working-environment). 

### Install DESeq2 by pasting the following commands in R
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
```

### Download [raw counts file](http://homer.ucsd.edu/zeyang/BENG183/combine_sample_raw_counts.txt) used for this tutorial and Homework 4
This file has been processed in a way so that it only contains gene lengths and raw counts from the featureCounts output. The raw counts for all the four samples that you mapped in the previous tutorial have been combined into this file. Other columns in featureCounts output were left out since they will not be used in the following analysis. You can use "wget" in Linux to download the file:
```bash
wget http://homer.ucsd.edu/zeyang/BENG183/combine_sample_raw_counts.txt
```

## Replicates & experiments quality check
After obtaining reads counts for every gene, we can first check which of these four samples have similar/different gene expression profiles. To compare gene expression across samples, the first step is always normalization. The following codes compute TPM values from raw counts:
```R
# First, add the folder where you store the downloaded raw counts to the paths that R searches for
setwd("/path/to/your/data/")
table <- read.table("combine_sample_raw_counts.txt", header = T,stringsAsFactors = F,row.names = 1) # read in raw count matrix
counts <- table[, c(2:5)] # store raw counts
gene.length <- table[, 1] # store gene lengths
library.size <- colSums(counts) # compute total counts of reads for each sample
counts.normalized.by.length <- counts[,c(1:4)]/gene.length # compute raw counts normalized for gene length

# compute TPM values, refer to the formula to better understand why it is done this way
tpm <- 1e6*t(t(counts.normalized.by.length[c(1:nrow(counts.normalized.by.length)),])/colSums(counts.normalized.by.length))
head(tpm) # take a look at the TPM values
colSums(tpm) # make sure that TPM values add up to the same number for each sample
# convert TPM values to log2 scale
log2.tpm <- log2(tpm+1)
```
After obtaining TPM values, you are now able to compare the gene expression across samples. To check the correlation between replicates, you can use function "plot" in R to visualize the expression of all the genes as a scatter plot. Here is a toy example showing how you can use function "plot" to do that. Replace "X" and "Y" by the columns in table "log2.tpm" to plot for log2 TPM values of replicates.
```R
X <- c(1,2,3,4,5)
Y <- c(2,4,6,8,10)
plot(X, Y, xlab="X sample replicate 1", ylab="X sample replicate 2")
```
<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/Rplot_toy.png">

To compare multiple replicates from different experiments, let's see how these samples cluster with each other. Refer to [previous discussion notes](https://github.com/Irenexzwen/BIOE183/blob/master/Discussion/DiscussionTutorial_ClusterAnalysis.md#3-cluster-analysis-in-r) for more information about how to do clustering in R. The following commands conduct hierarchical clustering on the log2 TPM values of the four samples. 
```R
d <- dist(t(log2.tpm), method = "euclidean") # compute distance
hc1 <- hclust(d, method = "complete") # conduct hierarchical clustering using complete linkage
plot(hc1) # plot dendrogram
```

## Differential analysis using [DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

### Read in raw counts output from featureCounts
```R
setwd("/path/to/your/data/")
table <- read.table("combine_sample_raw_counts.txt", header = T,stringsAsFactors = F,row.names = 1)
counts <- table[, c(2:5)] # load the raw counts and leave out gene lengths
head(counts) # look at what the table looks like
```

### Run DESeq2
```R
# Load DESeq2 package
library(DESeq2)

# Create a data frame that includes the information of replicates for the samples
labels <- factor(c("head","head","midgut","midgut")) # specify that the first two samples are replicates of "head" and the last two samples are replicates of "midgut"
labels.df <- DataFrame(condition = labels, row.names = colnames(counts))

# Load the raw counts, replicate information into DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = labels.df,
                              design = ~ condition)

# Run DESeq2
dds <- DESeq(dds)

# Show results that include fold changes and p-values
res <- results(dds)
res.df <- as.data.frame(res)
head(res.df)
```

### Visualize differentially expressed genes
DESeq2 provides a convenient visualization function "plotMA" to quickly check where differentially expressed genes (red dots) are distributed with respect to their expression levels and fold changes. 
```R
plotMA(res)
```
<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/Rplot_DESeq2Plot.png">

You can also visualize the actual expression of differential genes in the samples. To do this, we first sort out differentially expressed genes from the outputs of DESeq2. To identify differential genes, the best practice is to set cutoffs for both fold change and significance value. Here, we use column “log2FoldChange” that includes log2 transformed fold changes, and column “padj” that contains adjusted p-values corrected for multiple testing. To learn more about the approach DESeq2 used for multiple testing correction, please refer to the [Benjamini–Hochberg procedure](https://en.wikipedia.org/wiki/False_discovery_rate#Benjamini%E2%80%93Hochberg_procedure). 
```R
res.df[is.na(res.df[,2]),2] <- 0 # replace null in log2 fold change by 0
res.df[is.na(res.df[,6]),6] <- 1 # replace null in adjusted p-values by 1
DE.bools <- abs(res.df[,2]) > 1 & res.df[,6] < 0.05 # filter for fold change greater than 2 and adjusted p-value less than 0.05
```
Then we can use R function "heatmap()" to plot a heatmap of the log2 TPM values across samples.
```R
heatmap(log2.tpm[DE.bools,])
```

### Find functional enrichment of differentially expressed genes
We can use [Metascape](http://metascape.org/gp/index.html#/main/step1) to find the biological pathways that differential genes are enriched for. Input to Metascape is simply a list of genes (either gene names or gene ids). Let's first output the differentially expressed genes into a text file.
```R
DE.list <- row.names(res.df)[DE.bools] # extract gene IDs for differential genes
write.table(as.data.frame(DE.list),file="~/DE.txt",quote=F,sep=",",row.names=F) # save genes in a text file in your home directory
```
Now let's go the [Metascape](http://metascape.org/gp/index.html#/main/step1) website and conduct functional enrichment analysis for these genes. 

**Step 1: upload the text file that include the differential genes**

<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/Metascape_step1.png">

**Step 2: select Drosophila (D. melanogaster) as input sepcies and analysis species**

<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/Metascape_step2.png">

**Step 3: hit "Express Analysis" and wait for the results (usually take less than 2 minutes)**

**Step 4: see the analysis reports and read the top pathways that the differential genes are enriched for under "Bar Graph Summary" section**

<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/Metascape_reports.png">

