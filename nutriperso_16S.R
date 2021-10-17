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

# Start with OTU approach (distribution and taxonomies in separate tables)
otus <- read.delim("NUTRIPERSO_assembled_350_TAB.txt")
tax  <- read_tsv("NUTRIPERSO_assembled_350_OTU.tax", col_names = F)
# Parsing failure rows 190 and 507 (concatenated cols), to remove

# Tidy taxonomy table and filter by threshold >0.5 or >0.7
thr <- 0.7
taxmat <- tax %>% remove_constant() %>% select(-X2) %>% slice(-c(190, 507)) %>%
  filter(X5 > thr & X8 > thr & X11 > thr & X14 > thr & X17 > thr & X20 > thr) %>%
  select(-X5, -X8, -X11, -X14, -X16, -X17, -X19, -X20) %>% as.matrix

nrow(taxmat) #632 OTUs remaining
colnames(taxmat) <- c("otu", "Domain", "Phylum", "Class", "Order", "Family", "Genus")

# Get overlap between kept OTUs above threshold and all in otu matrix
vec <- intersect(taxmat[, 1], otus$OTU)

# Prepare data for phyloseq. See: http://joey711.github.io/phyloseq/import-data.html
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

# Add sample data
physeq1 <- merge_phyloseq(physeq, sampledata)
physeq1
# Get phylogenetic tree?

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




# What is the .fna file for?
seqs <- read.delim("NUTRIPERSO_assembled_350_OTU.fna")

# Import with mothur
dat <- import_mothur(mothur_shared_file = "NUTRIPERSO_assembled_350_TAB.txt",
                     mothur_constaxonomy_file = "NUTRIPERSO_assembled_350_ASV.tax")

import("mothur", "NUTRIPERSO_assembled_350_TAB.txt")

### ASV approach. Data and taxonomies in the same table, with taxonomies at the end
asvs <- read.delim("NUTRIPERSO_assembled_350_TAB_ASV.txt")
tax <- asvs %>% select(OTU0:X1.0.5)

#tax1 <- read_tsv("NUTRIPERSO_assembled_350_ASV.tax", col_names = F) 


# Aligned sequences
#seqs <- read.delim("NUTRIPERSO_assembled_350_OTU.fna")
taxmat <- tax1 %>% remove_constant() %>% 
  select(-X1, -X5, -X8, -X11, -X14, -X16, -X17, -X19, -X20) %>% as.matrix

# Initial 16S data (abandonned)
#dat16s <- read_xlsx("RelAb_OTUs_Tab_NP_Lepage.xlsx")
