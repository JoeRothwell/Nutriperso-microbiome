# Script to get phylogenetic trees from seqeuence data for phyloseq 
# The UniFrac metric uses it for beta diversity calculations, as does the phylogenetic alpha diversity metric, Faith’s PD

# Info: https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html#construct_phylogenetic_tree
# Also: https://rachaellappan.github.io/16S-analysis/index.html

# Read in fna file and
# Use msa package to align 16s sequences (uses command line tools clustal and muscle)
# Only muscle worked (takes a while)

library(msa)
dat1 <- readDNAStringSet("NUTRIPERSO_assembled_350_OTU.fna")
retained.otu <- dat1[rownames(otumat1)]
otu.align <- msa(retained.otu, type = "dna", method = "Muscle")

# Convert to a seqinr object
library(seqinr)
otu.align1 <- msaConvert(otu.align, type = "seqinr::alignment")
otu.dist <- dist.alignment(otu.align1, "identity")

# Neighbour-joining tree estimation in ape
library(ape)

# Here is a random tree in ape:
#random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
#plot(random_tree)

# Tree of our data
otuTree <- nj(otu.dist)
plot(otuTree)

# Can also use phyloseq (and phangorn)
library(phyloseq)
plot_tree(otuTree, size="Abundance", color="Sample", base.spacing=0.03)

# Other things investigated but not used
# Read in the OTU sequences from the fna file with seqinr (not used)
seqs <- read.fasta("NUTRIPERSO_assembled_350_OTU.fna")

# Phangorn produces more sophisticated trees
library(phangorn)
# Need to fit a model to the data and then extract $tree
# Also see: https://fuzzyatelin.github.io/bioanth-stats/module-24/module-24.html
# Also: https://bioinformatics.stackexchange.com/questions/7019/how-to-create-phylogenetic-trees-from-fasta-files-in-python-or-r/7021
# Not sure how to do it this way

# Reading in data with ape
dat <- read.dna("NUTRIPERSO_assembled_350_OTU.fna", format = "fasta")