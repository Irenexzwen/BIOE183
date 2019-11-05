# Tutorial4 Analysis of differentially expressed genes

## Preparations
### R basics
If you would like to brush up R basics, find the following resources:
**resources**: 
[- R basics](https://www.quora.com/What-are-some-good-resources-for-learning-R-1)
[- R bioconductor](https://www.coursera.org/learn/bioconductor)
[- R basics from previous discussion session](https://github.com/Irenexzwen/BIOE183/blob/master/Discussion/DiscussionTutorial_ClusterAnalysis.md#2-r-basics)

R is a language and environment for statistical computing and graphics. It is wildly used among bioinformaticians not only because it's pretty easy to work with table data, but also because bioconductor (a pool of packages target for high-throughput biological assays) was released on R platform.

RStudio is recommended for this tutorial even though it is completely doable to run all commands in R. To install RStudio, please refer to [the previous discussion notes](https://github.com/Irenexzwen/BIOE183/blob/master/Discussion/DiscussionTutorial_ClusterAnalysis.md#1-prepare-the-r-working-environment). 

### Install required R packages
```R
# install package "dplyr"
install.packages("dplyr")

# install package "DESeq2"
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")

# install package "ggplot2" for plotting results
install.packages("ggplot2")

# install package "pheatmap" to plot heatmaps
install.packages("pheatmap")

# install package "ggfortify" to conduct PCA analysis
install.packages("ggfortify")
```

### Download raw counts file used for this tutorial
Download the raw counts for the four samples you mapped in the previous tutorial [here](http://sysbio.ucsd.edu/public/wenxingzhao/CourseFall2019/DS_quant/all_sample_count.txt).

## Experiment/sample quality check
After obtaining the mapped reads and reads counts for every gene, we can first check whether these four samples have very different gene expression profiles. To compare tens of thousands genes automatically, we can use the following visualization methods:
1) Heatmap:
In R you could simply plot a heatmap use `pheatmap` package:
```R
library(pheatmap)
setwd("/path/to/your/data/")
table <- read.table("combine_sample_raw_counts.txt", header = T,stringsAsFactors = F,row.names = 1)
pheatmap::pheatmap(log(table+1),show_rownames = F)
```

2) PCA plot
```R
ggplot2::autoplot(ggfortify::prcomp(t(table)),label=TRUE)
```
Principal component analysis (PCA) is the most simple but direct way to visualize a high dimentional data in a low dimensional space.
If you're not familiar with PCA, please make sure understand it use any resource you could find! 

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
labels <- c("head","head","midgut","midgut") # showing the first two samples are replicates and the last two samples are replicates
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
Visualize top DE genes (sorted by adjusted pvalue)
```R

deheat <- table[rownames(head(degenes,20)),] %>% as.matrix
pheatmap::pheatmap(log(deheat+1))
```

