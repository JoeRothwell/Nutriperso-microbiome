# Microbiome metadata
library(readr)
library(readxl)
library(haven)
library(janitor)
library(tidyverse)

# Case-control data (initially sent to INRA)
dat <- read_xlsx("NP_RawSeq_ID.xlsx")
# E3N variables
#micro <- read_xls("microbiome_20210106.xls")
micro <- read_sas("nutriperso_20210304.sas7bdat")

# Analysis of 16s tables. Two approaches are used, OTU and ASV. 
# There are two tables per method: the distributions (.txt) and the taxonomic 
# assignments for each OTU (see emails from Patricia, 18th May and 11th Oct)
# Explanation of OTUs and ASVs at the beginning of this guide:
# https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html#abstract

# Start with OTU approach (distribution and taxonomies in separate tables)
otus <- read.delim("NUTRIPERSO_assembled_350_TAB.txt")
tax  <- read_tsv("NUTRIPERSO_assembled_350_OTU.tax", col_names = F)
# Parsing failure rows 190 and 507 (concatenated cols), to remove

# Tidy taxonomy table and filter by threshold >0.5 or >0.7
thr <- 0.7
taxmat <- tax %>% remove_constant() %>% select(-X2) %>% slice(-c(190, 507)) %>%
  filter(X5 > thr & X8 > thr & X11 > thr & X14 > thr & X17 > thr & X20 > thr) %>%
  select(-X5, -X8, -X11, -X14, -X17, -X20) %>% as.matrix

nrow(taxmat) #632 OTUs remaining (or 584)
colnames(taxmat) <- c("otu", "Domain", "Phylum", "Class", "Order", "Family", "Genus")

# Get overlap between kept OTUs above threshold and all in otu matrix
vec <- intersect(taxmat[, 1], otus$OTU)

# Prepare data for phyloseq. See: http://joey711.github.io/phyloseq/import-data.html
# Also see https://www.yanh.org/2021/01/01/microbiome-r/
# Also see Semi's script microbiome_nfbc1

# Remove OTU names and add back as rownames. Subset OTUs by their names in taxmat
otumat <- otus[, -1]
rownames(otumat) <- otus[, 1]
otumat1 <- otumat[vec, ]

# Add colnames for taxonomies
taxmat <- taxmat[, -1]
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
rownames(taxmat) <- rownames(otumat1)

# Make phyloseq object using constructors like otu_table
library(phyloseq)
OTU = otu_table(otumat1, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)

# Combine objects
physeq <- phyloseq(OTU, TAX)
physeq
plot_bar(physeq, fill = "Phylum")

# Create sample data
rownames(dat) <- str_remove(dat$ID, pattern = "-")
sampdata1 <- dat[colnames(otumat1), ]
# Make phyloseq object
sampledata <- sample_data(data.frame(sampdata1, row.names = sample_names(physeq),
                                     stringsAsFactors = F))

# Add sample data and phylogenetic tree
physeq1 <- merge_phyloseq(physeq, sampledata, otuTree)
physeq1


# Now just need the phylogenetic tree for downstream analysis and to add this to phyloseq
# Some info here: https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html#construct_phylogenetic_tree
# Also: https://rachaellappan.github.io/16S-analysis/index.html . From this site:
# the UniFrac metric uses it for beta diversity calculations, as does the phylogenetic  
# alpha diversity metric Faith’s PD

# Use msa package for to align 16s sequences (uses command line tools clustal and muscle)
# Only muscle worked
library(msa)
dat1 <- readDNAStringSet("NUTRIPERSO_assembled_350_OTU.fna")
retained.otu <- dat1[rownames(otumat1)]
otu.align <- msa(retained.otu, type = "dna", method = "Muscle")

# Convert to a seqinr object
otu.align1 <- msaConvert(otu.align, type = "seqinr::alignment")
otu.dist <- dist.alignment(otu.align1, "identity")
# Neighbour-joining tree estimation in ape
otuTree <- nj(otu.dist)
plot(otuTree)


# Other things investigated but not used
# Here is a random tree in ape:
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
plot(random_tree)

# Read in the OTU sequences from the fna file with seqinr (not used)
library(seqinr)
seqs <- read.fasta("NUTRIPERSO_assembled_350_OTU.fna")

library(phangorn)
# Need to fit a model to the data and then extract $tree
# Also see: https://fuzzyatelin.github.io/bioanth-stats/module-24/module-24.html
# Also: https://bioinformatics.stackexchange.com/questions/7019/how-to-create-phylogenetic-trees-from-fasta-files-in-python-or-r/7021
# Not sure how to do it this way

# Reading in data with ape
library(ape)
dat <- read.dna("NUTRIPERSO_assembled_350_OTU.fna", format = "fasta")

# Use msa package for alignment (uses command line tools clustal and muscle)


# Plots. see:
# http://joey711.github.io/phyloseq/plot_bar-examples#enterotypes_dataset_examples
# Default plot is abundance
plot_bar(physeq1)
plot_bar(physeq1, fill = "Phylum")

# Group by a variable
plot_bar(physeq1, x="casnutpkt", fill = "Phylum")
plot_bar(physeq1, x="diabete_groupe", fill = "Phylum")

# Facetted (doesn't work well)
plot_bar(physeq1, fill = "Phylum", facet_grid = ~ diabete_groupe)

# Subsets
ntaxa(physeq1)
nsamples(physeq1)

# Pre-processing. Transform to relative abundances (scale 0-1), filter
NPr  = transform_sample_counts(physeq1, function(x) x / sum(x) )
NPfr = filter_taxa(NPr, function(x) mean(x) > 1e-5, TRUE)

# Subset controls
NPfr.cont <- subset_samples(NPfr, casnutpkt == "TEMOIN" & RUN == "RUN1")
nsamples(NPfr.cont) #69 samples
plot_bar(NPfr.cont, fill = "Phylum")



# Import with mothur
library(mothur)
dat <- import_mothur(mothur_shared_file = "NUTRIPERSO_assembled_350_TAB.txt",
                     mothur_constaxonomy_file = "NUTRIPERSO_assembled_350_ASV.tax")

import("mothur", "NUTRIPERSO_assembled_350_TAB.txt")

### ASV approach. Data and taxonomies in the same table, with taxonomies at the end
asvs <- read.delim("NUTRIPERSO_assembled_350_TAB_ASV.txt")
tax <- asvs %>% select(OTU0:X1.0.5)

#tax1 <- read_tsv("NUTRIPERSO_assembled_350_ASV.tax", col_names = F) 


# Aligned sequences
seqs <- read.delim("NUTRIPERSO_assembled_350_OTU.fna")
taxmat <- tax1 %>% remove_constant() %>% 
  select(-X1, -X5, -X8, -X11, -X14, -X16, -X17, -X19, -X20) %>% as.matrix

# Initial 16S data (abandonned)
#dat16s <- read_xlsx("RelAb_OTUs_Tab_NP_Lepage.xlsx")
