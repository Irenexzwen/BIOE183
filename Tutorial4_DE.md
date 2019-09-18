# Tutorial4 Differential expressed gene analysis

## What is differentially expressed genes?
- basic principle
- biological / technical replicates

## DE analysis without biological replicates

- Null hypothesis
- Fold Change estimates
- which distribution

## DE analysis with biological replicates
- Principles
- Introduction to how DESeq2 or edgeR works

[short reference](http://chagall.med.cornell.edu/RNASEQcourse/Slides_July2019_Day4.pdf)


## DE analysis in R

### R basics 
R is a language and environment for statistical computing and graphics. 
It is wildly used among bioinformaticians not only because it's pretty easy to work with table data, 
but also because bioconductor(a pool of packages target for high-throughput biological assays) was released on R platform. 

It's always good for you to learn some R basics. Use more google if any task sounds hard to you! 

**resources**: 
[- R basics](https://www.quora.com/What-are-some-good-resources-for-learning-R-1)
[- R bioconductor](https://www.coursera.org/learn/bioconductor)


### DE analysis with DESeq2
```R
# google how to install these packages if you haven't installed yet
requrie(dplyr)
require(DESeq2)
require(ggplot2)
require(pheatmap)
require(ggfortify)

# read in the summarize raw count table

table <- read.table("D:/UCSD/class/all_sample_count.txt",header = T,stringsAsFactors = F,row.names = 1)

# have a glance of what the four sample looks like, using a PCA plot
ggplot2::autoplot(prcomp(t(table)),label=TRUE)


# To know what's going on please see the function defination below.
degenes <- Desq2(table[,c(1,2)],table[,c(3,4)],c("head","midgut")) %>% as.data.frame() 


# See the differential expressed genes in four samples using heatmap

deheat <- table[rownames(head(degenes,20)),] %>% as.matrix
pheatmap::pheatmap(log(deheat+1))

Desq2 <- function (count1, count2, group_pattern, thred = 20) 
{
    table <- data.frame(cbind(count1, count2))
    pat1 <- group_pattern[1]
    pat2 <- group_pattern[2]
    groups <- rep(pat1, length(colnames(table)))
    groups[grep(pat2, colnames(table))] <- pat2
    counts <- table[rowSums(table) > 0, ]
    nGenes <- length(counts[, 1])
    coverage <- colSums(counts)/nGenes
    counts <- counts[, coverage > 10]
    groups <- groups[coverage > 10]
    groups <- factor(groups, levels = c(pat1, pat2))
    coverage <- coverage[coverage > 1]
    library(DESeq2)
    library("BiocParallel")
    dds <- DESeqDataSetFromMatrix(counts, DataFrame(groups), 
        ~groups)
    register(MulticoreParam(thred))
    dds <- DESeq(dds, parallel = TRUE)
    res <- results(dds)
    find.significant.genes <- function(de.result, alpha = 0.05) {
        filtered <- de.result[(de.result$padj < alpha) & !is.infinite(de.result$log2FoldChange) & 
            !is.nan(de.result$log2FoldChange) & !is.na(de.result$padj), 
            ]
        sorted <- filtered[order(filtered$padj), c(1, 2, 6)]
    }
    de2.genes <- find.significant.genes(res)
    return(de2.genes)
}

```

After we load the data, we first want to see whether these four samples have very different different profile. 
Since there are tens of thousands genes we can't check them manually. We usually use some visualization method to do that:
1) Heatmap:
In R you could simply plot a heatmap use `pheatmap` package:
```R
require(pheatmap)
table <- read.table("D:/UCSD/class/all_sample_count.txt",header = T,stringsAsFactors = F,row.names = 1)
pheatmap::pheatmap(log(table+1),show_rownames = F)
```
<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/all_gene_heatmap.png">

2) PCA plot
```R
ggplot2::autoplot(ggfortify::prcomp(t(table)),label=TRUE)
```
Principal component analysis (PCA) is the most simple but direct way to visualize a high dimentional data in a low dimensional space.
If you're not familiar with PCA, please make sure understand it use any resource you could find! 

<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/PCA_DS.png">

3) Visualize top DE genes (sorted by adjusted pvalue)
```R

deheat <- table[rownames(head(degenes,20)),] %>% as.matrix
pheatmap::pheatmap(log(deheat+1))
```
<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/topDE.png">
