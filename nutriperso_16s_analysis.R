# Diversity measures and descriptive analysis
load("nutriperso_phyloseq.rda")
library(phyloseq)
library(tidyverse)

# Remind size and participant variables
colnames(sample_data(nutri))
ntaxa(nutri)
nsamples(nutri)


# Abundance of taxa
plot_bar(nutri, fill = "Phylum") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Subset the phyloseq object otherwise plots too hard to read
nutri0 <- subset_samples(nutri, casnutpkt == "TEMOIN" & RUN == "RUN1")
nsamples(nutri0) #69 samples

# Plots. Note: full data too large to visualise or see case/control differences
# see: http://joey711.github.io/phyloseq/plot_bar-examples#enterotypes_dataset_examples
# Default plot is abundance with stacked OTUs
plot_bar(nutri0, fill = "Phylum") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Remind participant variables
colnames(sample_data(nutri))

# Group by variables (hard to see because many OTUs)
plot_bar(nutri0, x="tabac", fill = "Phylum")

# Facetted (too many samples to read)
plot_bar(nutri0, fill = "Phylum", facet_grid = ~ casnutpkt)

# Heatmap with clustering based on ordination methods
# https://joey711.github.io/phyloseq/plot_heatmap-examples.html
# Too many low OTU values; no obvious differences between samples
library(scales)
plot_heatmap(nutri0, na.value = "black")

# Plot tree coloured by phylum (too complex to read)
plot_tree(nutri, color="Phylum", size="abundance")
plot_tree(nutri, ladderize="left", color="Class")
plot_tree(nutri, color="Phylum") + coord_polar(theta="y")


# Alpha and beta diversity (within and between sample variability)
### Measures of alpha diversity by status. By default gives all measures
plot_richness(nutri, x="casnutpkt", color = "casnutpkt") + theme_bw()
plot_richness(nutri, x="tabac", color = "casnutpkt")
plot_richness(nutri, x="casnutpkt", measures = c("Shannon", "Simpson"))

# Jittered
plot_richness(nutri, x="casnutpkt", color = "casnutpkt", measures = c("Shannon", "Simpson")) +
  theme_bw() + theme(legend.title = element_blank(), axis.title.x = element_blank()) +
  geom_point(colour = "white") + theme(legend.position="none") +
  geom_jitter()

# Estimate alpha diversity for epidemiological models (based on diversity() from vegan)
div <- estimate_richness(nutri, measures = c("Shannon", "Simpson"))

# Beta-diversity analysis and enterotypes. Ordination plots (for samples)
# http://joey711.github.io/phyloseq/plot_ordination-examples.html

# Plot of samples
ord <- ordinate(nutri, "NMDS", "bray")
plot_ordination(nutri, ord, type="samples", color="casnutpkt") + 
  geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() + theme(legend.title = element_blank())

plot_ordination(nutri, ord, type="samples", color="casnutpkt", shape="diabete_groupe")
# More info on enterotypes at https://enterotype.embl.de/enterotypes.html

# Can also plot ordination of OTUs
plot_ordination(nutri, ord, type="taxa", color="Phylum", title="taxa") + theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") 



# Get enterotypes (for full detail, see enterotypes.R)
load("nutriperso_phyloseq.rda")
nutri.rel  <- transform_sample_counts(nutri, function(x) x / sum(x) )
otus <- otu_table(nutri.rel)

data.dist <- dist.JSD(otus)
data.dist1 <- phyloseq::distance(nutri.rel, method="jsd")

# Do a test run on our dataset with number of clusters set to 3
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

# Optimal number of clusters is 2; reset to 2
data.cluster <- pam.clustering(data.dist, k = 2)

# Cluster validation (index 3 is the sil_width for each obs)
obs.silhouette <- mean(silhouette(data.cluster, data.dist)[ , 3])
cat(obs.silhouette) #0.1209738 (acceptable)

# Remove low abundance genera
noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe -> Matrix
  bigones <- rowSums(Matrix) * 100 / (sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)
}

otus <- noise.removal(otus, percent=0.01)


## Multi-dimentional scaling plots to visualise enterotypes
## Get data for BCA and PCA
library(ade4)
obs.pca <- dudi.pca(data.frame(t(otus)), scannf = F, nf = 10)
obs.bet <- bca(obs.pca, fac = as.factor(data.cluster), scannf = F, nf = k-1) 

## Between class analysis (plot 1)
dev.new()
s.class(obs.bet$ls, fac = as.factor(data.cluster), grid = F, sub = "Between-class analysis")
# Doesn't work with only 2 clusters

# Principal coordinate analysis (plot 2)
obs.pcoa <- dudi.pco(data.dist, scannf = F, nf = 3)
#dev.new()
s.class(obs.pcoa$li, fac = as.factor(data.cluster), grid = F, sub = "Principal coordinate analysis")

  

