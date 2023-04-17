# Calculation of enterotypes using different methods

library(tidyverse)
# Example of how to calculate enterotypes from https://enterotype.embl.de/enterotypes.html

# Get table of relative abundances from https://enterotype.embl.de/MetaHIT_SangerSamples.genus.txt 
# all samples add up to 1
dat <- read.table("MetaHIT_SangerSamples.genus.txt", header=T, row.names=1, dec=".", sep="\t")
dat <- dat[-1, ]

# Function to determine the Jensen-Shannon Divergence, a distance metric

dist.JSD <- function(inMatrix, pseudocount = 0.000001, ...) {
  KLD <- function(x, y) sum(x * log(x/y))
  JSD <- function(x, y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix, 1:2, function(x) ifelse (x == 0, pseudocount, x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i, j] <- JSD(as.vector(inMatrix[ ,i]),
                             as.vector(inMatrix[ ,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix) -> resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}

# Apply function to data
data.dist <- dist.JSD(dat)

# Use partitioning around medoids (PAM) clustering to cluster the abundance profiles
# Similar to k-means but supports any distance measure
# x is a distance matrix and k the number of clusters

library(cluster)

pam.clustering <- function(x, k) {
  clust <- as.vector(pam(as.dist(x), k, diss = TRUE)$clustering)
  return(clust)
}

# Get clusters
data.cluster <- pam.clustering(data.dist, k = 3)

# Assess optimal number of clusters
require(clusterSim)
nclusters <- index.G1(t(dat), data.cluster, d = data.dist, centrotypes = "medoids")

# Try 1-20 clusters and plot optimal CH index
nclusters <- NULL

for (k in 1:20) { 
  if (k == 1) {
    nclusters[k] = NA 
  } else {
    data.cluster_temp <- pam.clustering(data.dist, k)
    nclusters[k] <- index.G1(t(dat),data.cluster_temp,  d = data.dist, centrotypes = "medoids")
  }
}

plot(nclusters, type="h", xlab="k clusters", ylab="CH index", main="Optimal number of clusters")

# Optimal cluster number is 3

# Cluster validation

obs.silhouette <- mean(silhouette(data.cluster, data.dist)[ , 3])
cat(obs.silhouette) #0.1899451

dat <- noise.removal(dat, percent=0.01)

## Multi-dimentional scaling plots to visualise enterotypes
## Between class analysis (plot 1)
library(ade4)
obs.pca <- dudi.pca(data.frame(t(dat)), scannf = F, nf = 10)
obs.bet <- bca(obs.pca, fac = as.factor(data.cluster), scannf = F, nf = k-1) 
#dev.new()
s.class(obs.bet$ls, fac = as.factor(data.cluster), grid = F, sub = "Between-class analysis")

# Principal coordinate analysis (plot 2)
obs.pcoa <- dudi.pco(data.dist, scannf = F, nf = 3)
#dev.new()
s.class(obs.pcoa$li, fac = as.factor(data.cluster), grid = F, sub = "Principal coordinate analysis")


#### Nutriperso data

load("nutriperso_phyloseq.rda")
nutri.rel  <- transform_sample_counts(nutri, function(x) x / sum(x) )
otus <- otu_table(nutri.rel)

data.dist <- dist.JSD(otus)

# Get clusters (need to specify how many)
data.cluster <- pam.clustering(data.dist, k = 3)

# Assess optimal number of clusters
library(clusterSim)
nclusters <- index.G1(t(otus), data.cluster, d = data.dist, centrotypes = "medoids")

# Try 1-20 clusters and plot optimal CH index
nclusters <- NULL

for (k in 1:20) { 
  if (k == 1) {
    nclusters[k] = NA 
  } else {
    data.cluster_temp <- pam.clustering(data.dist, k)
    nclusters[k] <- index.G1(t(otus), data.cluster_temp, d = data.dist, centrotypes = "medoids")
  }
}

plot(nclusters, type="h", xlab="k clusters", ylab="CH index", main="Optimal number of clusters")
# Optimal number of clusters is 2

data.cluster <- pam.clustering(data.dist, k = 2)

# Cluster validation

obs.silhouette <- mean(silhouette(data.cluster, data.dist)[ , 3])
cat(obs.silhouette) #0.1899451

otus <- noise.removal(otus, percent=0.01)

## Multi-dimentional scaling plots to visualise enterotypes
## Between class analysis (plot 1)
library(ade4)
obs.pca <- dudi.pca(data.frame(t(otus)), scannf = F, nf = 10)
obs.bet <- bca(obs.pca, fac = as.factor(data.cluster), scannf = F, nf = k-1) 
#dev.new()
s.class(obs.bet$ls, fac = as.factor(data.cluster), grid = F, sub = "Between-class analysis")

# Principal coordinate analysis (plot 2)
obs.pcoa <- dudi.pco(data.dist, scannf = F, nf = 3)
#dev.new()
s.class(obs.pcoa$li, fac = as.factor(data.cluster), grid = F, sub = "Principal coordinate analysis")

# Using the mbOmic bioconductor package and data table
library(BiocManager)
BiocManager::install("mbOmic")
library(mbOmic)
library(data.table)

# From https://bioconductor.org/packages/release/bioc/vignettes/mbOmic/inst/doc/enterotyping.html
# Data at http://enterotypes.org/ref_samples_abundance_MetaHIT.txt

dat <- read.delim("ref_samples_abundance_MetaHIT.txt")
dat <- impute::impute.knn(as.matrix(dat), k = 100)
dat <- as.data.frame(dat$data + 0.001) 
setDT(dat, keep.rownames = TRUE)
dat

# Convert to bSet object and do enterotyping
dat <- bSet(b = dat)
res <- estimate_k(dat)
res
ret <- enterotyping(dat, res$verOptCluster) 
ret


### Get enterotypes for our Nutriperso data. Load phyloseq object
library(phyloseq)

# First transform OTU table to relative abundance
nutri.rel  <- transform_sample_counts(nutri, function(x) x / sum(x) )
otus <- otu_table(nutri.rel)
data.dist <- dist.JSD(otus)

# Get clusters
data.cluster <- pam.clustering(data.dist, k=3)

# Assess optimal number of clusters
require(clusterSim)
nclusters = index.G1(t(otus), data.cluster, d = data.dist, centrotypes = "medoids")

# Try 1-20 clusters and plot optimal CH index
nclusters=NULL

for (k in 1:20) { 
  if (k == 1) {
    nclusters[k] = NA 
  } else {
    data.cluster_temp = pam.clustering(data.dist, k)
    nclusters[k] = index.G1(t(otus), data.cluster_temp, d = data.dist,
                          centrotypes = "medoids")
  }
}

plot(nclusters, type="h", xlab="k clusters", ylab="CH index",main="Optimal number of clusters")
data.cluster

# Package method. Transform to relative abudance
nutri.rel  <- transform_sample_counts(nutri, function(x) x / sum(x) )
otus1 <- data.frame(otus)
library(impute)
otus1 <- impute.knn(as.matrix(otus1), k = 100)
otus1 <- as.data.frame(otus1$data + 0.001) 
setDT(otus1, keep.rownames = TRUE)

otus1

# Convert to bSet object and do enterotyping
otus1 <- bSet(b = otus1)
res <- estimate_k(otus1)
res
ret <- enterotyping(otus1, res$verOptCluster) 
ret

### With phyloseq? Has distance methods
# https://joey711.github.io/phyloseq/distance.html
# Get the JSD on the phyloseq object
nutri.rel  <- transform_sample_counts(nutri, function(x) x / sum(x) )
jsdist <- phyloseq::distance(nutri.rel, method="jsd") # or JSD(physeq)
# Get clusters
clust <- as.vector(pam(jsdist, 3, diss = TRUE)$clustering)
clust <- pam(jsdist, 3, diss = TRUE)$clustering


jsMDS <- ordinate(nutri, "MDS", distance = jsdist)
p <- plot_ordination(nutri, jsMDS, color="casnutpkt", shape = clust)




# Get clusters
data.cluster <- pam.clustering(jsdist, k=3)

# Assess optimal number of clusters
require(clusterSim)
nclusters <- index.G1(t(otu_table(nutri.rel)), data.cluster, d = jsdist, centrotypes = "medoids")
