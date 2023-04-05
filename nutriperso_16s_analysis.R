# Diversity measures and descriptive analysis
load("nutriperso_phyloseq.rda")
library(phyloseq)
library(tidyverse)

# Remind size and participant variables
colnames(sample_data(nutri))
ntaxa(nutri)
nsamples(nutri)

### Measures of alpha diversity by status. By default gives all measures
plot_richness(nutri, x="casnutpkt", color = "casnutpkt")
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
  theme_bw() + theme(legend.title = element_blank())
plot_ordination(nutri, ord, type="samples", color="casnutpkt", shape="diabete_groupe")
# More info on enterotypes at https://enterotype.embl.de/enterotypes.html

# Can also plot ordination of OTUs
plot_ordination(nutri, ord, type="taxa", color="Phylum", title="taxa")



### Other plots of descriptive data
### Subset the phyloseq object otherwise plots too hard to read
nutri0 <- subset_samples(nutri, casnutpkt == "TEMOIN" & RUN == "RUN1")
nsamples(nutri0) #69 samples

# Plots. Note: full data too large to visualise or see case/control differences
# see: http://joey711.github.io/phyloseq/plot_bar-examples#enterotypes_dataset_examples
# Default plot is abundance with stacked OTUs
plot_bar(nutri0, fill = "Phylum")

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

