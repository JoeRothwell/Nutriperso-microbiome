# Preparation of 16s microbiome data and import into a phyloseq object
library(readr)
library(janitor)
library(tidyverse)

# Preparation of participant data ----
# Data for case-control initially sent to INRA
library(readxl)
dat <- read_xlsx("NP_RawSeq_ID.xlsx")

# E3N variables
micro <- read_xls("microbiome_20210106.xls")
library(haven)
# micro <- read_sas("nutriperso_20210304.sas7bdat")

# Join 2 sets of participant data
dat$ident <- as.character(dat$ident)
dat <- left_join(dat, micro, by = "ident")

# 16s data
# Two sets of data are provided, OTU and ASV. Explanation at:
# https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html#abstract
# There are two tables per method: the distributions (.txt) and the taxonomic assignments for each OTU 
# See emails from Patricia, 18th May and 11th Oct



### Preparation of OTU data ----
# Start with OTU approach (distribution and taxonomies in separate tables)
# For ASV see below
otus <- read.delim("NUTRIPERSO_assembled_350_TAB.txt")
tax  <- read_tsv("NUTRIPERSO_assembled_350_OTU.tax", col_names = F)
# Parsing failure rows 190 and 507 (concatenated cols), to remove

# Tidy taxonomy table and filter by threshold >0.5 or >0.7
thr <- 0.7
taxmat <- tax %>% 
  remove_constant() %>% 
  select(-X2) %>% dplyr::slice(-c(190, 507)) %>%
  filter(X5 > thr & X8 > thr & X11 > thr & X14 > thr & X17 > thr & X20 > thr) %>%
  select(-X5, -X8, -X11, -X14, -X17, -X20) %>% as.matrix

nrow(taxmat) #632 OTUs remaining (or 584 for 0.7 threshold)
colnames(taxmat) <- c("otu", "Domain", "Phylum", "Class", "Order", "Family", "Genus")

# Get overlap between kept OTUs above threshold and all in otu matrix
vec <- intersect(taxmat[, 1], otus$OTU)

# Prepare otu and tax data for phyloseq. See: http://joey711.github.io/phyloseq/import-data.html
# Also see https://www.yanh.org/2021/01/01/microbiome-r/
# Also see Semi's script microbiome_nfbc1.R

# Remove OTU names and add back as rownames. Subset OTUs by their names in taxmat
otumat <- otus[, -1]
rownames(otumat) <- otus[, 1]
otumat1 <- otumat[vec, ]

# Add colnames for taxonomies
taxmat <- taxmat[, -1]
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
rownames(taxmat) <- rownames(otumat1)



### Get phylogenetic tree ----
# Some info here: https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html#construct_phylogenetic_tree
# Also: https://rachaellappan.github.io/16S-analysis/index.html . From this site:
# the UniFrac metric uses it for beta diversity calculations, as does the phylogenetic  
# alpha diversity metric, Faith’s PD

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
# seqs <- read.fasta("NUTRIPERSO_assembled_350_OTU.fna")

# Phangorn produces more sophisticated trees
library(phangorn)
# Need to fit a model to the data and then extract $tree
# Also see: https://fuzzyatelin.github.io/bioanth-stats/module-24/module-24.html
# Also: https://bioinformatics.stackexchange.com/questions/7019/how-to-create-phylogenetic-trees-from-fasta-files-in-python-or-r/7021
# Not sure how to do it this way

# Reading in data with ape
# fas <- read.dna("NUTRIPERSO_assembled_350_OTU.fna", format = "fasta")



### Make phyloseq object ----
# (using constructor functions)
OTU = otu_table(otumat1, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)

# Combine objects
physeq <- phyloseq(OTU, TAX)
physeq
plot_bar(physeq, fill = "Phylum")

# Get particpant data (sample data)
rownames(dat) <- str_remove(dat$ID, pattern = "-")
sampdata1 <- dat[colnames(otumat1), ]
sampledata <- sample_data(data.frame(sampdata1, row.names = sample_names(physeq),
                                     stringsAsFactors = F))

# Add sample data and phylogenetic tree
nutri <- merge_phyloseq(physeq, sampledata, otuTree)
nutri

### Overview of data
ntaxa(nutri)
nsamples(nutri)
sample_names(nutri)[1:5]

### Pre-process data. Transform to relative abundance and filter taxa
# nutrir  <- transform_sample_counts(nutri, function(x) x / sum(x) )
# nutrifr <- filter_taxa(nutri, function(x) mean(x) > 1e-5, TRUE)
# ntaxa(nutrifr) # Only 258 taxa remain

### Need to use prune_samples as well?

save(nutri, file = "nutriperso_phyloseq.rda")


### Preparation of ASV data ----
# Read in ASVs, removing taxonomies tagged onto end of table
asvs <- read.delim("NUTRIPERSO_assembled_350_TAB_ASV.txt") %>% select(-(OTU0:X1.0.5))

# Read in taxonomies
tax <- read_tsv("NUTRIPERSO_assembled_350_ASV.tax", col_names = F) 

# Tidy taxonomy table and filter by threshold >0.5 or >0.7
thr <- 0.7
taxmat <- tax %>% 
  remove_constant() %>% 
  select(-X2) %>% #dplyr::slice(-c(190, 507)) %>%
  filter(X5 > thr & X8 > thr & X11 > thr & X14 > thr & X17 > thr & X20 > thr) %>%
  select(-X5, -X8, -X11, -X14, -X17, -X20) %>% as.matrix

nrow(taxmat) #2219 ASVs remaining for 0.7 threshold
colnames(taxmat) <- c("otu", "Domain", "Phylum", "Class", "Order", "Family", "Genus")

# Get overlap between kept OTUs above threshold and all in otu matrix
vec <- intersect(taxmat[, 1], asvs$OTU)

# In the asv data, remove OTU names and add back as rownames. Subset OTUs by their names in taxmat
asvmat <- asvs[, -1]
rownames(asvmat) <- asvs[, 1]
asvmat1 <- asvmat[vec, ]

# Add colnames for taxonomies
taxmat <- taxmat[, -1]
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
rownames(taxmat) <- rownames(asvmat1)



# Get phylogenetic tree. Used by Faith's PD (alpha) and UniFrac (beta-diversity) metrics
# Use msa package to align 16s sequences (uses command line tools clustal and muscle)
# Only muscle worked (takes a while)
dat1 <- readDNAStringSet("NUTRIPERSO_assembled_350_OTU.fna")
retained.otu <- dat1[rownames(asvmat1)]
library(msa)
otu.align <- msa(retained.otu, type = "dna", method = "Muscle")

# Convert to a seqinr object
library(seqinr)
otu.align1 <- msaConvert(otu.align, type = "seqinr::alignment")
otu.dist <- dist.alignment(otu.align1, "identity")

# Neighbour-joining tree estimation in ape
library(ape)
