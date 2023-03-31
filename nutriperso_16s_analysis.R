# Analysis of Nutriperso Phyloseq object
# Descriptive analyses and diversity (alpha and beta)
load("nutriperso_phyloseq.rda")
library(phyloseq)

<<<<<<< HEAD
# Plots. see:
# http://joey711.github.io/phyloseq/plot_bar-examples#enterotypes_dataset_examples
# Default plot is abundance
plot_bar(physeq1)
plot_bar(physeq1, fill = "Phylum") #+ theme(axis.labels.x = element_blank())
=======
# Plots. see: http://joey711.github.io/phyloseq/plot_bar-examples#enterotypes_dataset_examples
# Default plot is abundance with stacked OTUs
plot_bar(nutri)
plot_bar(nutri, fill = "Phylum") #+ theme(axis.labels.x = element_blank())
>>>>>>> 540025a9495e676ff66aa33c420b02147885bdfa

# Remind partipant variables
colnames(sample_data(nutri))

# Group by variables (hard to see because many OTUs)
plot_bar(nutri, x="casnutpkt", fill = "Phylum")
plot_bar(nutri, x="diabete_groupe", fill = "Phylum")

# Facetted (too many samples to read)
plot_bar(nutri, fill = "Phylum", facet_grid = ~ casnutpkt)

# Heatmap with clustering based on ordination methods
# https://joey711.github.io/phyloseq/plot_heatmap-examples.html
# Too many low OTU values; no obvious differences between samples
library(scales)
plot_heatmap(nutri, na.value = "black")

# Plot tree with labelled diabetes group (too complex to read)
plot_tree(nutri, color="Phylum", shape="diabete_groupe", size="abundance")
plot_tree(nutri)

# Plot a tree, mapping colour to case-control status or class
plot_tree(nutri, ladderize="left", color="casnutpkt")
plot_tree(nutri, ladderize="left", color="Class")

plot_tree(physeq1, color="casnutpkt", ladderize="left") + coord_polar(theta="y")

# Subsets
ntaxa(nutri)
nsamples(nutri)

# Measures of alpha diversity. By default gives all measures
plot_richness(nutri, x="diabete_groupe")
plot_richness(nutri, x="casnutpkt", color = "casnutpkt")
plot_richness(nutri, x="diabete_groupe", measures = c("Shannon", "Simpson"))

# Estimate alpha diversity for epidemiological models (based on diversity() from vegan)
estimate_richness(nutri, measures = c("Shannon", "Simpson"))

# Beta-diversity analysis and enterotypes
# Ordination plots
# Plot of OTUs
# http://joey711.github.io/phyloseq/plot_ordination-examples.html
ord <- ordinate(nutri, "NMDS", "bray")
plot_ordination(nutri, ord, type="taxa", color="Phylum", title="taxa")

# Plot of samples
<<<<<<< HEAD
plot_ordination(physeq1, GP.ord, type="samples", color="casnutpkt") + 
  theme_bw() + theme(legend.title = element_blank())
plot_ordination(physeq1, GP.ord, type="samples", color="casnutpkt", shape="diabete_groupe")

# Measures of alpha diversity
plot_richness(physeq1, x="diabete_groupe")
plot_richness(physeq1, x="casnutpkt", color = "casnutpkt")
# Just selected measures
plot_richness(physeq1, x="casnutpkt", color = "casnutpkt", measures = c("Shannon", "Simpson")) +
  theme_bw() + theme(legend.title = element_blank(), axis.title.x = element_blank()) +
  geom_point(colour = "white") + theme(legend.position="none") +
  geom_jitter()
=======
plot_ordination(nutri, ord, type="samples", color="casnutpkt", shape="diabete_groupe")

>>>>>>> 540025a9495e676ff66aa33c420b02147885bdfa

# More info on enterotypes at https://enterotype.embl.de/enterotypes.html


# Subset controls (need to change names)
NPfr.cont <- subset_samples(NPfr, casnutpkt == "TEMOIN" & RUN == "RUN1")
nsamples(NPfr.cont) #69 samples
plot_bar(NPfr.cont, fill = "Phylum")
