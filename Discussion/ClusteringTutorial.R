### Load Packages ###
library(tidyverse)        # data manipulation
library(cluster)          # clustering algorithms
library(factoextra)       # clustering algorithms & visualization
library(cluster.datasets) # sample datasets for cluster analysis
library(gridExtra)        # for plotting multiple plots in grid
library(dendextend)       # for comparing two dendrograms


### Data Preparation ###
data(acidosis.patients) # load the dataset
ds <- acidosis.patients # import data into variable
ds <- na.omit(ds)       # remove missing values from data
ds <- scale(ds)         # scale/standardize the data


### Compute and Visualize Distance Matrix ###
distance <- get_dist(ds, method = "euclidean")
fviz_dist(distance, order=FALSE, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))


### K-means Clustering ###
# k = 2 clusters
k2 <- kmeans(ds, centers = 2, nstart = 25)
str(k2) # print the output

fviz_cluster(k2, data = ds) # visualize on cluster plot

# try different numbers of clusters
k3 <- kmeans(ds, centers = 3, nstart = 25)
k4 <- kmeans(ds, centers = 4, nstart = 25)
k5 <- kmeans(ds, centers = 5, nstart = 25)

# plots to compare
p1 <- fviz_cluster(k2, geom = "point", data = ds) + ggtitle("k = 2")
p2 <- fviz_cluster(k3, geom = "point",  data = ds) + ggtitle("k = 3")
p3 <- fviz_cluster(k4, geom = "point",  data = ds) + ggtitle("k = 4")
p4 <- fviz_cluster(k5, geom = "point",  data = ds) + ggtitle("k = 5")

grid.arrange(p1, p2, p3, p4, nrow = 2) # arrange plots on a grid

# determine optimal number of clusters
set.seed(123) # set seed for random number generator to ensure reproducibility
fviz_nbclust(ds, kmeans, method = "wss")        # using the "elbow method"
fviz_nbclust(ds, kmeans, method = "silhouette") # using the "silhouette method"

# Compute k-means clustering with k = 3
final <- kmeans(ds, 3, nstart = 25)
print(final)

fviz_cluster(final, data = ds)

# extract clusters and add to initial data to do some descriptive statistics at the cluster level
acidosis.patients %>%
  mutate(Cluster = final$cluster) %>%
  group_by(Cluster) %>%
  summarise_all("mean")


### Hierarchical Clustering ###
d <- dist(ds, method = "euclidean") # Dissimilarity matrix

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete")
plot(hc1, cex = 0.6, hang = -1)

# Compute with agnes using Complete Linkage
hc2 <- agnes(ds, method = "complete")
hc2$ac # Agglomerative coefficient

# Assess different linkage methods
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")
# function to compute coefficient
ac <- function(x) {
  agnes(ds, method = x)$ac
}
map_dbl(m, ac)

# Visualize the dendrogram from agnes
hc3 <- agnes(ds, method = "ward")
pltree(hc3, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 

# compute divisive hierarchical clustering
hc4 <- diana(ds)
hc4$dc # Divisive coefficient
pltree(hc4, cex = 0.6, hang = -1, main = "Dendrogram of diana") #plot dendrogram

# Cutting a dendrogram
hc5 <- hclust(d, method = "ward.D2") # Ward's method
sub_grp <- cutree(hc5, k = 3)        # Cut tree into 3 groups
table(sub_grp)                       # Number of members in each cluster

# Visualize dendrogram with borders around clusters
plot(hc5, cex = 0.6)
rect.hclust(hc5, k = 3, border = 2:4)

# Using cutree with agnes and diana
# Cut agnes() tree into 3 groups
hc_a <- agnes(ds, method = "ward")
cutree(as.hclust(hc_a), k = 3)

# Cut diana() tree into 3 groups
hc_d <- diana(ds)
cutree(as.hclust(hc_d), k = 3)
