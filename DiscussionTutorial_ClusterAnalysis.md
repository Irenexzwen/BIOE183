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
```
# An example of computing the mean with variables

mydata$sum <- mydata$x1 + mydata$x2
mydata$mean <- (mydata$x1 + mydata$x2)/2
```

### Importing Data
Importing data into R is fairly simple. R offers options to import many file types, from CSVs to databases. For example, this is how to import a CSV into R.
```
# first row contains variable names, comma is separator
# assign the variable id to row names
# note the / instead of \ on mswindows systems

mydata <- read.table("c:/mydata.csv", header=TRUE,
   sep=",", row.names="id")
```

## Packages
Packages are collections of R functions, data, and compiled code in a well-defined format. The directory where packages are stored is called the library. R comes with a standard set of packages. Others are available for download and installation. Once installed, they have to be loaded into the session to be used.
```
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
## Package Requirements
For this tutorial, you will need to load the following packages:
```
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
```
data(acidosis.patients)  # load the dataset
ds <- acidosis.patients  # import data into variable
```
To remove any missing values that might be present in the data, type this:
```
ds <- na.omit(ds)
```
Next, we scale/standardize the data using the R function `scale`:
```
ds <- scale(ds)
```

## Clustering Distance Measures
To classify observations into groups, one needs to compute the distance or (dis)similarity between each pair of observations. The result of this computation is called a dissimilarity or distance matrix. There are many methods to calculate distance measures, including *Euclidean*, *Manhattan*, and *correlation-based* distances.

Within R it is simple to compute and visualize the distance matrix using the functions `get_dist` and `fviz_dist` from the `factoextra` R package. This starts to illustrate which patients have large dissimilarities (red) versus those that appear to be fairly similar (teal).
- `get_dist`: for computing a distance matrix between the rows of a data matrix. The default distance computed is the Euclidean; however, `get_dist` also supports other methods such as Manhattan, Pearson, Spearman, etc.
- `fviz_dist`: for visualizing a distance matrix
```
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
```
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
```
fviz_cluster(k2, data = ds)
```
<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/ClusteringTutorial_plot1.png">

Because the number of clusters (*k*) must be set before we start the algorithm, it is often advantageous to use several different values of *k* and examine the differences in the results. We can execute the same process for 3, 4, and 5 clusters, and the results are shown in the figure:
```
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

## Determining Optimal Clusters
Since the number of clusters is set by the analyst, how does the analyst know what is the optimal number of clusters? There are several methods, two of which we address here:

### Elbow Method:
Recall that, the basic idea behind cluster partitioning methods, such as k-means clustering, is to define clusters such that the total intra-cluster variation (known as total within-cluster variation or total within-cluster sum of square) is minimized. Thus, we can use the following algorithm to define the optimal clusters:
  1) Compute clustering algorithm (e.g., k-means clustering) for different values of *k*. For instance, by varying *k* from 1 to 10 clusters.
  2) For each *k*, calculate the total within-cluster sum of square ("wss").
  3) Plot the curve of wss according to the number of clusters *k*.
  4) The location of a bend (elbow) in the plot is generally considered as an indicator of the appropriate number of clusters.

We can implement this in R with the function `fviz_nbclust` by specifying `method = "wss"`. The results suggest that 3 is the optimal number of clusters as it appears to be the bend in the "elbow".
```
set.seed(123) # set seed for random number generator to ensure reproducibility
fviz_nbclust(ds, kmeans, method = "wss")
```

<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/ClusteringTutorial_elbowplot.png">

### Average Silhouette Method:
In short, the average silhouette approach measures the quality of a clustering. That is, it determines how well each object lies within its cluster. A high average silhouette width indicates a good clustering. The average silhouette method computes the average silhouette of observations for different values of *k*. The optimal number of clusters *k* is the one that maximizes the average silhouette over a range of possible values for *k*.

Again, we can implement this in R with the function `fviz_nbclust` by specifying `method = "silhouette"`.
```
fviz_nbclust(ds, kmeans, method = "silhouette")
```
<img src="https://github.com/Irenexzwen/BIOE183/blob/master/images/ClusteringTutorial_silhouetteplot.png">


**More to be added soon...**

# References
This tutorial is adapted from the following material:
1) https://www.statmethods.net/r-tutorial/index.html
2) http://uc-r.github.io/kmeans_clustering
3) https://uc-r.github.io/hc_clustering
