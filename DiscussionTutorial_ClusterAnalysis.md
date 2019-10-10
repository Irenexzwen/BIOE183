# Discussion Tutorial: Cluster Analysis Tools in R
In this tutorial, you will:
  1) [Set up R environment on your laptop](#1-prepare-the-r-working-environment),
  2) [Learn the basics of the R language](#2-r-basics),
  3) [Perform cluster analysis in R](#3-cluster-analysis-in-r).
  
## What is R?
R is a language and environment for statistical computing and graphics. It provides a wide variety of statistical (linear and nonlinear modelling, classical statistical tests, time-series analysis, classification, clustering, …) and graphical techniques, and is highly extensible. It includes:
- an effective data handling and storage facility,
- a suite of operators for calculations on arrays, in particular matrices,
- a large, coherent, integrated collection of intermediate tools for data analysis,
- graphical facilities for data analysis and display either on-screen or on hardcopy, and
- a well-developed, simple and effective programming language which includes conditionals, loops, user-defined recursive functions and input and output facilities.

More information and documentation about R can be found on their [website](https://www.r-project.org/).

# 1) Prepare the R working environment
First, download and install R on your computer from the website [here](https://cran.rstudio.com/). R is supported on Windows, Mac, and Linux platforms.

Next, download and install [RStudio Desktop](https://rstudio.com/products/rstudio/download/). RStudio provides a convenient GUI to work with R and includes a console, syntax-highlighting editor that supports direct code execution, and a variety of robust tools for plotting, viewing history, debugging and managing your workspace.

**Note:** For Windows users, it is recommended that you also install [RTools](https://cran.rstudio.com/bin/windows/Rtools/), as it is necessary for building R packages in Windows.

Once everything is installed, simply find and launch RStudio from your computer. 

# 2) R basics
R is a command line driven program. Within the console, the user enters commands at the prompt (> by default) and each command is executed one at a time.

The workspace is your current R working environment and includes any user-defined objects (vectors, matrices, data frames, lists, functions). At the end of an R session, the user can save an image of the current workspace that is automatically reloaded the next time R is started.

You can also write an R script, which is simply a text file containing a set of commands and comments. The script can be saved and used later to re-execute the saved commands.

## Operators in R
### Arithmetic Operators:
|Operator|Description   |
|--------|--------------|
|+       |Addition      |
|-       |Subtraction   |
|*       |Multiplication|
|/       |Division      |
|^ or ** |Exponentiation|

### Logical Operators:
|Operator|Description             |
|--------|------------------------|
|>       |Greater than            |
|>=      |Greater than or equal to|
|==      |Exactly equal to        |
|!=      |Not equal to            |

## Data Types
R has a wide variety of data types including scalars, vectors (numerical, character, logical), matrices, data frames, and lists.

### Creating New Variables
Use the assignment operator <- to create new variables.
```R
# An example of computing the mean with variables

mydata$sum <- mydata$x1 + mydata$x2
mydata$mean <- (mydata$x1 + mydata$x2)/2
```

### Importing Data
Importing data into R is fairly simple. R offers options to import many file types, from CSVs to databases. For example, this is how to import a CSV into R.
```R
# first row contains variable names, comma is separator
# assign the variable id to row names
# note the / instead of \ on mswindows systems

mydata <- read.table("c:/mydata.csv", header=TRUE,
   sep=",", row.names="id")
```

## Packages
Packages are collections of R functions, data, and compiled code in a well-defined format. The directory where packages are stored is called the library. R comes with a standard set of packages. Others are available for download and installation. Once installed, they have to be loaded into the session to be used.
```R
install.packages("packagename") # install a package
.libPaths()                     # get library location
library()                       # see all packages installed
library("packagename")          # load a package
search()                        # see packages currently loaded
```

## Useful Keyboard Shortcuts for RStudio
|Description|Windows/Linux|Mac|
|-----------|-------------|---|
|Clear console|Ctrl+L|Ctrl+L|
|Interrupt currently executing command|Esc|Esc|
|Run current line/selection (in a script)|Ctrl+Enter|Command+Enter|
|Comment/uncomment current line/selection|Ctrl+Shift+C|Command+Shift+C|

More shortcuts are listed [here](https://support.rstudio.com/hc/en-us/articles/200711853-Keyboard-Shortcuts).

# 3) Cluster Analysis in R
Now, we will use R to perform two types of clustering methods on a dataset:
  1) [K-Means Clustering](#k-means-clustering)
  2) [Hierarchical Clustering](#hierarchical-clustering)

For convenience, all of the code used in this tutorial has been compiled into a single R script [here](https://github.com/Irenexzwen/BIOE183/blob/master/ClusteringTutorial.R).

## Package Requirements
For this tutorial, you will need to load the following packages:
```R
library(tidyverse)        # data manipulation
library(cluster)          # clustering algorithms
library(factoextra)       # clustering algorithms & visualization
library(cluster.datasets) # sample datasets for cluster analysis
library(gridExtra)        # for plotting multiple plots in grid
library(dendextend)       # for comparing two dendrograms
```
Note: You may need to install the packages first by running `install.packages("packagename")` before you can load them.

## Data Preparation
To perform a cluster analysis in R, generally, the data should be prepared as follows:
  1) Rows are observations (individuals) and columns are variables
  2) Any missing value in the data must be removed or estimated.
  3) The data must be standardized (i.e., scaled) to make variables comparable. Recall that, standardization consists of transforming the variables such that they have mean zero and standard deviation one.

Here, we'll use a dataset from the "cluster.datasets" package called `acidosis.patients`, which contains measures of various compounds in cerebrospinal fluid and blood for 40 acidosis patients.
```R
data(acidosis.patients)  # load the dataset
ds <- acidosis.patients  # import data into variable
```
To remove any missing values that might be present in the data, type this:
```R
ds <- na.omit(ds)
```
Next, we scale/standardize the data using the R function `scale`:
```R
ds <- scale(ds)
```

## Clustering Distance Measures
To classify observations into groups, one needs to compute the distance or (dis)similarity between each pair of observations. The result of this computation is called a dissimilarity or distance matrix. There are many methods to calculate distance measures, including *Euclidean*, *Manhattan*, and *correlation-based* distances.

Within R it is simple to compute and visualize the distance matrix using the functions `get_dist` and `fviz_dist` from the `factoextra` R package. This starts to illustrate which patients have large dissimilarities (red) versus those that appear to be fairly similar (teal).
- `get_dist`: for computing a distance matrix between the rows of a data matrix. The default distance computed is the Euclidean; however, `get_dist` also supports other methods such as Manhattan, Pearson, Spearman, etc.
- `fviz_dist`: for visualizing a distance matrix
```R
distance <- get_dist(ds, method = "euclidean")
fviz_dist(distance, order=FALSE, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
```
<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/ClusteringTutorial_DistPlot.png">

## K-Means Clustering
K-means clustering is the most commonly used unsupervised machine learning algorithm for partitioning a given data set into a set of *k* groups (i.e. *k* clusters), where *k* represents the number of groups pre-specified by the analyst. The basic idea behind k-means clustering consists of defining clusters so that the total intra-cluster variation (known as total within-cluster variation) is minimized.

The K-means algorithm can be summarized as follows:
  1) Specify the number of clusters (*K*) to be created (by the analyst)
  2) Randomly select *k* objects from the data set as the initial cluster centers or means
  3) Assign each observation to their closest centroid, based on the Euclidean distance between the object and the centroid
  4) For each of the *k* clusters, update the cluster centroid by calculating the new mean values of all the data points in the cluster. The centroid of a *K*th cluster is a vector of length *p* containing the means of all variables for the observations in the *k*th cluster; *p* is the number of variables.
  5) Iteratively minimize the total within sum of squares. That is, iterate steps 3 and 4 until the cluster assignments stop changing or the maximum number of iterations is reached. By default, the R software uses 10 as the default value for the maximum number of iterations.

We can compute k-means in R with the `kmeans` function. Here will group the data into two clusters (`centers = 2`). The `kmeans` function also has an `nstart` option that attempts multiple initial configurations and reports on the best one. For example, adding `nstart = 25` will generate 25 initial configurations. This approach is often recommended.
```R
k2 <- kmeans(ds, centers = 2, nstart = 25)
str(k2)  # print the output
##List of 9
## $ cluster     : Named int [1:40] 2 2 2 2 2 2 2 2 2 1 ...
##  ..- attr(*, "names")= chr [1:40] "1" "2" "3" "4" ...
## $ centers     : num [1:2, 1:6] -0.851 0.15 1.66 -0.293 -1.816 ...
##  ..- attr(*, "dimnames")=List of 2
##  .. ..$ : chr [1:2] "1" "2"
##  .. ..$ : chr [1:6] "ph.cerebrospinal.fluid" "ph.blood" "hco3.cerebrospinal.fluid" "hco3.blood" ...
## $ totss       : num 234
## $ withinss    : num [1:2] 23.5 102.3
## $ tot.withinss: num 126
## $ betweenss   : num 108
## $ size        : int [1:2] 6 34
## $ iter        : int 1
## $ ifault      : int 0
## - attr(*, "class")= chr "kmeans"
```
The output of `kmeans` is a list with several bits of information. The most important being:
- `cluster`: A vector of integers (from 1:k) indicating the cluster to which each point is allocated.
- `centers`: A matrix of cluster centers.
- `totss`: The total sum of squares.
- `withinss`: Vector of within-cluster sum of squares, one component per cluster.
- `tot.withinss`: Total within-cluster sum of squares, i.e. sum(withinss).
- `betweenss`: The between-cluster sum of squares, i.e. totss - tot.withinss.
- `size`: The number of points in each cluster.

We can view our results by using `fviz_cluster` to provide a nice illustration of the clusters. If there are more than two dimensions (variables) `fviz_cluster` will perform principal component analysis (PCA) and plot the data points according to the first two principal components that explain the majority of the variance.
```R
fviz_cluster(k2, data = ds)
```
<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/ClusteringTutorial_plot1.png">

Because the number of clusters (*k*) must be set before we start the algorithm, it is often advantageous to use several different values of *k* and examine the differences in the results. We can execute the same process for 3, 4, and 5 clusters, and the results are shown in the figure:
```R
k3 <- kmeans(ds, centers = 3, nstart = 25)
k4 <- kmeans(ds, centers = 4, nstart = 25)
k5 <- kmeans(ds, centers = 5, nstart = 25)

# plots to compare
p1 <- fviz_cluster(k2, geom = "point", data = ds) + ggtitle("k = 2")
p2 <- fviz_cluster(k3, geom = "point",  data = ds) + ggtitle("k = 3")
p3 <- fviz_cluster(k4, geom = "point",  data = ds) + ggtitle("k = 4")
p4 <- fviz_cluster(k5, geom = "point",  data = ds) + ggtitle("k = 5")

grid.arrange(p1, p2, p3, p4, nrow = 2) # arrange plots on a grid
```
<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/ClusteringTutorial_plot2.png">

### Determining Optimal Clusters
Since the number of clusters is set manually, how does the analyst know what is the optimal number to use? There are several methods, two of which we address here:

#### Elbow Method:
Recall that, the basic idea behind cluster partitioning methods, such as k-means clustering, is to define clusters such that the total intra-cluster variation (known as total within-cluster variation or total within-cluster sum of square) is minimized. Thus, we can use the following algorithm to define the optimal clusters:
  1) Compute clustering algorithm (e.g., k-means clustering) for different values of *k*. For instance, by varying *k* from 1 to 10 clusters.
  2) For each *k*, calculate the total within-cluster sum of square ("wss").
  3) Plot the curve of wss according to the number of clusters *k*.
  4) The location of a bend (elbow) in the plot is generally considered as an indicator of the appropriate number of clusters.

We can implement this in R with the function `fviz_nbclust` by specifying `method = "wss"`. The results suggest that 3 is the optimal number of clusters as it appears to be the bend in the "elbow".
```R
set.seed(123) # set seed for random number generator to ensure reproducibility
fviz_nbclust(ds, kmeans, method = "wss")
```

<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/ClusteringTutorial_elbowplot.png">

#### Average Silhouette Method:
In short, the average silhouette approach measures the quality of a clustering. That is, it determines how well each object lies within its cluster. A high average silhouette width indicates a good clustering. The average silhouette method computes the average silhouette of observations for different values of *k*. The optimal number of clusters *k* is the one that maximizes the average silhouette over a range of possible values for *k*.

Again, we can implement this in R with the function `fviz_nbclust` by specifying `method = "silhouette"`.
```R
fviz_nbclust(ds, kmeans, method = "silhouette")
```
<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/ClusteringTutorial_silhouetteplot.png">

### Extracting Results
With these approaches suggesting 3 as the optimal number of clusters, we can perform the final analysis and extract the results using 3 clusters.
```R
final <- kmeans(ds, 3, nstart = 25)
print(final)
## K-means clustering with 3 clusters of sizes 5, 21, 14
##
## Cluster means:
##   ph.cerebrospinal.fluid    ph.blood hco3.cerebrospinal.fluid hco3.blood co2.cerebrospinal.fluid  co2.blood
## 1             -0.9504782  1.96494447               -1.9117531 -1.9340627              -1.8599395 -1.6311646
## 2             -0.3050189 -0.41230383               -0.1343910 -0.1713011              -0.2012997 -0.2759981
## 3              0.7969848 -0.08331014                0.8843555  0.9476883               0.9662137  0.9965559
##
## Clustering vector:
##  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 
##  2  2  2  2  2  2  2  2  2  1  2  2  1  1  1  1  2  3  2  2  3  3  3  3  3  3  3  3  3  3  3  3  2  2  2  2  3  2  2  2 
##
## Within cluster sum of squares by cluster:
## [1] 18.72524 29.97298 25.44634
##  (between_SS / total_SS =  68.3 %)
##
## Available components:
##
## [1] "cluster"      "centers"      "totss"        "withinss"     "tot.withinss" "betweenss"    "size"         "iter"        
## [9] "ifault" 
```
We can visualize the results using `fviz_cluster`:
```R
fviz_cluster(final, data = ds)
```
<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/ClusteringTutorial_finalkmeansplot.png">

And we can extract the clusters and add to our initial data to do some descriptive statistics at the cluster level:
```R
acidosis.patients %>%
  mutate(Cluster = final$cluster) %>%
  group_by(Cluster) %>%
  summarise_all("mean")
## # A tibble: 3 x 7
##   Cluster ph.cerebrospinal.fluid ph.blood hco3.cerebrospinal.fluid hco3.blood co2.cerebrospinal.fluid co2.blood
##     <int>                  <dbl>    <dbl>                    <dbl>      <dbl>                   <dbl>     <dbl>
## 1       1                   44.0     58.5                     12.8       8.38                    24.2      19.0
## 2       2                   46.4     38.0                     21.3      21.4                     43.3      33.5
## 3       3                   50.5     40.8                     26.2      29.7                     56.7      47.1
```

## Hierarchical Clustering
Hierarchical clustering is an alternative approach for identifying groups in a dataset. Unlike the k-means approach, it does not require us to pre-specify the number of clusters to be generated. Furthermore, hierarchical clustering has an added advantage over k-means clustering in that it results in an attractive tree-based representation of the observations, called a dendrogram.

### Types of Hierarchical Clustering
Hierarchical clustering can be divided into two main types: *agglomerative* and *divisive*.
  1) **Agglomerative clustering**: Also known as AGNES (Agglomerative Nesting), it works in a *bottom-up* manner. That is, each object is initially considered as a single-element cluster (leaf). At each step of the algorithm, the two clusters that are the most similar are combined into a new bigger cluster (nodes). This procedure is iterated until all points are member of just one single big cluster (root). The result is a tree which can be plotted as a dendrogram.
  2) **Divisive hierarchical clustering**: Also known as DIANA (Divise Analysis), it works in a *top-down* manner. The algorithm is an inverse order of AGNES. It begins with the root, in which all objects are included in a single cluster. At each step of iteration, the most heterogeneous cluster is divided into two. The process is iterated until all objects are in their own cluster.

Agglomerative clustering is good at identifying small clusters, while divisive hierarchical clustering is good at identifying large clusters.

### Cluster Linkage Methods
Previously, in the k-means approach, we measure the (dis)similarity of observations using distance measures (i.e., Euclidean distance, Manhattan distance, etc.). Now, for hierarchical clustering, we need to measure the (dis)similarity between two *clusters* of observations. This can be done using a number of different cluster agglomeration methods (linkage methods). More information about these methods can be found [here](https://support.minitab.com/en-us/minitab/19/help-and-how-to/modeling-statistics/multivariate/how-to/cluster-observations/methods-and-formulas/linkage-methods/).

### Implementing in R
There are different functions available in R for computing hierarchical clustering (HC). The commonly used functions are:
- `hclust` (in stats package) and `agnes` (in cluster package) for agglomerative HC
- `diana` (in cluster package) for divisive HC

#### Agglomerative HC
We can perform agglomerative HC with `hclust`. First we compute the dissimilarity values with `dist` and then feed these values into `hclust` and specify the agglomeration method to be used (i.e. “complete”, “average”, “single”, “ward.D”). We can then plot the dendrogram.
```R
d <- dist(ds, method = "euclidean")   # Dissimilarity matrix
hc1 <- hclust(d, method = "complete") # Hierarchical clustering using Complete Linkage
plot(hc1, cex = 0.6, hang = -1)       # Plot the resulting dendrogram
```
<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/ClusteringTutorial_dendro1.png">

Alternatively, we can use the `agnes` function. These functions behave very similarly; however, with the `agnes` function you can also get the *agglomerative coefficient*, which measures the amount of clustering structure found (values closer to 1 suggest strong clustering structure).
```R
hc2 <- agnes(ds, method = "complete") # Compute with agnes
hc2$ac                                # Agglomerative coefficient
## [1] 0.9064081
```

This number allows us to find certain hierarchical clustering methods that can identify stronger clustering structures. Here we see that Ward’s method identifies the strongest clustering structure of the four methods assessed.
```R
# Assess different linkage methods
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")
# function to compute coefficient
ac <- function(x) {
  agnes(ds, method = x)$ac
}
map_dbl(m, ac)
##   average    single  complete      ward 
## 0.8220451 0.7245571 0.9064081 0.9319864
```
We can visualize the dendrogram similar to before:
```R
hc3 <- agnes(ds, method = "ward")
pltree(hc3, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 
```
<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/ClusteringTutorial_dendro2.png">

#### Divisive HC
The R function `diana` provided by the cluster package allows us to perform divisive hierarchical clustering. `diana` works similar to `agnes`, and provides a similar value that measures the amount of clustering structure found called the *divisive coefficient*.
```R
# compute divisive hierarchical clustering
hc4 <- diana(ds)
hc4$dc # Divisive coefficient
## [1] 0.8955153

# plot the dendrogram
pltree(hc4, cex = 0.6, hang = -1, main = "Dendrogram of diana")
```
<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/ClusteringTutorial_dendro3.png">

### Interpreting Dendrograms
In the dendrograms above, each leaf corresponds to one observation. As we move up the tree, observations that are similar to each other are combined into branches, which are themselves fused at a higher height.

The (dis)similarity between two observations is indicated by the height of the fusion, provided on the vertical axis. The higher the height of the fusion, the less similar the observations are.

Using the function `cutree`, we can "cut" a dendrogram at a certain height to identify clusters. Different numbers of clusters are obtained by changing the height of the cut. It plays the same role as the *k* value in k-means clustering.
```R
hc5 <- hclust(d, method = "ward.D2") # Ward's method
sub_grp <- cutree(hc5, k = 3)        # Cut tree into 3 groups
table(sub_grp)                       # Number of members in each cluster
## sub_grp
##  1  2  3 
## 29  5  6
```
It’s possible to draw the dendrogram with a border around the clusters. The argument `border` is used to specify the border colors for the rectangles:
```R
plot(hc5, cex = 0.6)
rect.hclust(hc5, k = 3, border = 2:4)
```
<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/ClusteringTutorial_dendro3.png">

To use `cutree` with `agnes` and `diana` you can perform the following:
 ```R
 # Cut agnes() tree into 3 groups
hc_a <- agnes(ds, method = "ward")
cutree(as.hclust(hc_a), k = 3)

# Cut diana() tree into 3 groups
hc_d <- diana(ds)
cutree(as.hclust(hc_d), k = 3)
```

# References
The information in this tutorial is adapted from the following material:
1) https://www.statmethods.net/r-tutorial/index.html
2) http://uc-r.github.io/kmeans_clustering
3) https://uc-r.github.io/hc_clustering
