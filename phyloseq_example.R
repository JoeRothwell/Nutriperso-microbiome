# Phyloseq examples
# From https://joey711.github.io/phyloseq/preprocess.html
library(phyloseq)
packageVersion("phyloseq")

# Load example data
data(GlobalPatterns)

# Get basic info on tables
ntaxa(GlobalPatterns)
nsamples(GlobalPatterns)
sample_names(GlobalPatterns)[1:5]
rank_names(GlobalPatterns)
sample_variables(GlobalPatterns)
otu_table(GlobalPatterns)[1:5, 1:5]
tax_table(GlobalPatterns)[1:5, 1:4]
phy_tree(GlobalPatterns)
taxa_names(GlobalPatterns)[1:10]

# Get sums of abundances for taxa
Sums <- taxa_sums(GlobalPatterns)

# Get most abundant taxa
myTaxa <- names(sort(Sums, decreasing = TRUE)[1:10])

# Prune taxa and draw tree
ex1 <- prune_taxa(myTaxa, GlobalPatterns)
plot(phy_tree(ex1), show.node.label = TRUE)


plot_tree(ex1, color = "SampleType", label.tips = "Phylum", ladderize = "left", justify = "left" , size = "Abundance")


# Pre-processing. Transform to relative abundance to allow the filtering of 
# OTUs below a given threshold
GPr  = transform_sample_counts(GlobalPatterns, function(x) x / sum(x) )
GPfr = filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE)

GP.chl = subset_taxa(GlobalPatterns, Phylum=="Chlamydiae")
GP.chl = prune_samples(sample_sums(GP.chl)>=20, GP.chl)


