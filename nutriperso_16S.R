# Microbiome metadata
library(readr)
library(readxl)
library(haven)
library(janitor)
library(tidyverse)

# E3N data
micro <- read_xls("microbiome_20210106.xls")
micro1 <- read_sas("nutriperso_20210304.sas7bdat")
# Initial 16S data
dat16s <- read_xlsx("RelAb_OTUs_Tab_NP_Lepage.xlsx")

# Case-control data (initially sent to INRA)
meta <- read_xlsx("NP_RawSeq_ID.xlsx")

# Analysis of OTU tables
# Two approaches are used, OTU and ASV. 
# There are two tables per method: the distributions (.txt) and the taxonomic 
# assignments for each OTU (.tax, .fna, see email from Patricia, 18th May)

### ASV approach (distribution and taxonomies in the same table)
asvs <- read.delim("NUTRIPERSO_assembled_350_TAB_ASV.txt")
tax <- asvs %>% select(OTU0:X1.0.5)

### OTU approach (distribution and taxonomies in separate tables)
otus <- read.delim("NUTRIPERSO_assembled_350_TAB.txt")
tax1 <- read_tsv("NUTRIPERSO_assembled_350_ASV.tax", col_names = F) 
#seqs <- read.delim("NUTRIPERSO_assembled_350_OTU.fna")
taxmat <- tax1 %>% remove_constant() %>% 
  select(-X1, -X5, -X8, -X11, -X14, -X16, -X17, -X19, -X20) %>% as.matrix
ncol(taxmat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
rownames(taxmat) <- tax1$X1

# Prepare data for phyloseq
# http://joey711.github.io/phyloseq/import-data.html
# Using OTU table
otumat <- otus[, -1]
rownames(otumat) <- otus[, 1]

# Make phyloseq object
library(phyloseq)
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)

# Create sample data
sampledata = meta[, -1]
rownames(sampledata) <- str_replace(meta$ID, "-", "")

physeq = phyloseq(OTU, TAX)
physeq
plot_bar(physeq, fill = "Family")


dat <- import_mothur(mothur_shared_file = "NUTRIPERSO_assembled_350_TAB.txt",
                     mothur_constaxonomy_file = "NUTRIPERSO_assembled_350_ASV.tax")

import("mothur", "NUTRIPERSO_assembled_350_TAB.txt")
