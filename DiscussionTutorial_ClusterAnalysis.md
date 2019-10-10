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
ds <- acidosis.patients
```
To remove any missing values that might be present in the data, type this:
```
ds <- na.omit(ds)
```
Next, we scale/standardize the data using the R function `scale`:
```
ds <- scale(ds)
```

**More to be added soon...**

# References
This tutorial is adapted from the following material:
1) https://www.statmethods.net/r-tutorial/index.html
2) http://uc-r.github.io/kmeans_clustering
3) https://uc-r.github.io/hc_clustering
