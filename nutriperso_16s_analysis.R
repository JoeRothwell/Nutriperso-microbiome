library(phyloseq)

# Plots. see:
# http://joey711.github.io/phyloseq/plot_bar-examples#enterotypes_dataset_examples
# Default plot is abundance
plot_bar(physeq1)
plot_bar(physeq1, fill = "Phylum") #+ theme(axis.labels.x = element_blank())

# Group by a variable
plot_bar(physeq1, x="casnutpkt", fill = "Phylum")
plot_bar(physeq1, x="diabete_groupe", fill = "Phylum")

# Facetted (doesn't work well)
plot_bar(physeq1, fill = "Phylum", facet_grid = ~ diabete_groupe)

# Heatmap (too complex)
plot_heatmap(physeq1)

# Plot tree with labelled diabetes group (too complex to read)
plot_tree(physeq1, color="Phylum", shape="diabete_groupe", size="abundance")
plot_tree(physeq1)

# Subsets
ntaxa(physeq1)
nsamples(physeq1)

# Ordination plots
# Plot of OTUs
# http://joey711.github.io/phyloseq/plot_ordination-examples.html
GP.ord <- ordinate(physeq1, "NMDS", "bray")
plot_ordination(physeq1, GP.ord, type="taxa", color="Phylum", title="taxa")

# Plot of samples
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

# Plot a tree, mapping colour to case-control status or class
plot_tree(physeq1, ladderize="left", color="casnutpkt")
plot_tree(physeq1, ladderize="left", color="Class")

plot_tree(physeq1, color="casnutpkt", ladderize="left") + coord_polar(theta="y")


# Subset controls (need to change names)
NPfr.cont <- subset_samples(NPfr, casnutpkt == "TEMOIN" & RUN == "RUN1")
nsamples(NPfr.cont) #69 samples
plot_bar(NPfr.cont, fill = "Phylum")
