# Adapted from Semi's script

library(ape)
library(ggplot2)
library(gplots)
library(lme4)
library(phangorn)
library(plotly)
library(tidyverse)
library(vegan)
library(VennDiagram)
library(foreign)
library(haven)
library(Hmisc)
library(expss)
library(biomformat)
library(wesanderson)
library(xtable)
library(ggpubr)
library(car)
library(RVAideMemoire)
library(yingtools2)
library(gridExtra)
library(ppcor)
library(gamlss)
library(phyloseq)
library(microbiome)
library(knitr)
library(DESeq2)
library(MiRKAT)
library(gamlss)
library(reshape2)

setwd("C:/Users/zouiouichs/Documents/nfbc_2018/final_databases")
set.seed(5) 

# ------------------------------------------- BUILDING MICROBIOME DATASET -------------------------------------------

# ------------------------------------------- DATA IMPORT ------------------------------------------- 

# import biom file (from QIIME)
dat <- read_biom("C:/Users/zouiouichs/Documents/nfbc_2018/sample_updated/QIIME1/finland_sample_gg.biom")

# import taxonomy file
tax <- observation_metadata(dat)

# import OTU table
# in our case, there was a problem with the ID of some participants, this is why we have to import another OTU table
# all core-metrics for non-QC samples determined based on sample name with NFB (585)
OTU <- read.csv("C:/Users/zouiouichs/Documents/nfbc_2018/sample_updated/OTU.csv")
OTU$ID <- OTU$sample_name
# import file with the correct participant IDs
sample_id<-read.table("C:/Users/zouiouichs/Documents/nfbc_2018/sample/sampleid2.txt",
                      sep='\t', header=T, comment="")
sample_id$ID <- sample_id$X.SampleID
sample_id$X.SampleID <- NULL
# merge the OTU table with the correct IDs
OTU<- right_join(x = OTU, y = sample_id, by = "ID")
OTU$X.SampleID <- NULL
# extract the correct ID from the full string
OTU$ID <- substr(OTU$X.realid, 7, 15)
OTU <- OTU %>%
  dplyr::select(ID, everything())
# delete duplicates
which(duplicated(OTU$ID))
OTU <- OTU[!duplicated(OTU$ID), ]
OTU$X.realid <- NULL
OTU$sample_name <- NULL

# import phylogenetic tree
NJ.tree<-read.tree("C:/Users/zouiouichs/Documents/nfbc_2018/sample/tree.nwk")

# import metadata
load(file = "metadata.R")

# ------------------------------------------- DATA CLEANING/DATA PREPARATION ------------------------------------------- 

# OTU table
row.names(OTU) = OTU$ID
OTU.clean = OTU[,-which(names(OTU) %in% "ID")]
row.names(metadata) = metadata$ID
OTU.clean = OTU.clean[order(row.names(OTU.clean)),]
metadata = metadata[order(row.names(metadata)),]
# remove all the participants that don't occur in the metadata
OTU.clean = OTU.clean[row.names(OTU.clean) %in% row.names(metadata),]
# remove missing values from OTU table
OTU.clean.filtered <- OTU.clean
OTU.clean.filtered[is.na(OTU.clean.filtered)] <- 0
OTU.clean.filtered <- OTU.clean.filtered[-which(rowSums(OTU.clean.filtered)==0),]
OTU.clean.filtered <- OTU.clean.filtered[,-which(colSums(OTU.clean.filtered)==0)]

# taxonomy table
# remove all the OTUs that don't occur in our OTU.clean data set
tax.clean = tax[row.names(tax) %in% colnames(OTU.clean),]
# rename a column in R
colnames(tax.clean)[colnames(tax.clean)=="taxonomy1"] <- "Kingdom"
colnames(tax.clean)[colnames(tax.clean)=="taxonomy2"] <- "Phylum"
colnames(tax.clean)[colnames(tax.clean)=="taxonomy3"] <- "Class"
colnames(tax.clean)[colnames(tax.clean)=="taxonomy4"] <- "Order"
colnames(tax.clean)[colnames(tax.clean)=="taxonomy5"] <- "Family"
colnames(tax.clean)[colnames(tax.clean)=="taxonomy6"] <- "Genus"
colnames(tax.clean)[colnames(tax.clean)=="taxonomy7"] <- "Species"
# remove confidence column
tax.clean = tax.clean[,-which(names(tax.clean) %in% "Confidence")]
# remove species column
tax.clean = tax.clean[,-which(names(tax.clean) %in% "Species")]

# metadata
# create categorical variables based on continuous variables
# for insulin resistance
# homa
metadata$homa_quartiles <- NA
metadata$homa_quartiles[which(metadata$HOMAIR_46<quantile(metadata$HOMAIR_46, prob=0.25, na.rm=T) )] <- 0
metadata$homa_quartiles[which(metadata$HOMAIR_46>=quantile(metadata$HOMAIR_46, prob=0.25, na.rm=T)&metadata$HOMAIR_46<quantile(metadata$HOMAIR_46, prob=0.5, na.rm=T) )] <- 1
metadata$homa_quartiles[which(metadata$HOMAIR_46>=quantile(metadata$HOMAIR_46, prob=0.5, na.rm=T)&metadata$HOMAIR_46<quantile(metadata$HOMAIR_46, prob=0.75, na.rm=T) )] <- 2
metadata$homa_quartiles[which(metadata$HOMAIR_46>=quantile(metadata$HOMAIR_46, prob=0.75, na.rm=T) )] <- 3
metadata <- apply_labels(metadata, 
                         homa_quartiles	=	"Insulin resistance variable based on quartiles of HOMAIR_46",
                         homa_quartiles=c("Q1"=0, "Q2"=1,"Q3"=2,"Q4"=3))
describe(metadata$homa_quartiles)
table(metadata$homa_quartiles)
# visceral fat
metadata$visceral_quartiles <- NA
metadata$visceral_quartiles[which(metadata$C6646C_BI_005<quantile(metadata$C6646C_BI_005, prob=0.25, na.rm=T) )] <- 0
metadata$visceral_quartiles[which(metadata$C6646C_BI_005>=quantile(metadata$C6646C_BI_005, prob=0.25, na.rm=T)&metadata$C6646C_BI_005<quantile(metadata$C6646C_BI_005, prob=0.5, na.rm=T) )] <- 1
metadata$visceral_quartiles[which(metadata$C6646C_BI_005>=quantile(metadata$C6646C_BI_005, prob=0.5, na.rm=T)&metadata$C6646C_BI_005<quantile(metadata$C6646C_BI_005, prob=0.75, na.rm=T) )] <- 2
metadata$visceral_quartiles[which(metadata$C6646C_BI_005>=quantile(metadata$C6646C_BI_005, prob=0.75, na.rm=T) )] <- 3
metadata <- apply_labels(metadata, 
                         visceral_quartiles	=	"Visceral fat area variable based on quartiles of C6646C_BI_005",
                         visceral_quartiles=c("Q1"=0, "Q2"=1,"Q3"=2,"Q4"=3))

# for inflammation
# hsCRP - high sensitivity c-reactive protein
metadata$crp_quartiles <- NA
metadata$crp_quartiles[which(as.numeric(metadata$C6646B_Serum1_hsCRP.x)<quantile(as.numeric(metadata$C6646B_Serum1_hsCRP.x), prob=0.25, na.rm=T) )] <- 0
metadata$crp_quartiles[which(as.numeric(metadata$C6646B_Serum1_hsCRP.x)>=quantile(as.numeric(metadata$C6646B_Serum1_hsCRP.x), prob=0.25, na.rm=T)&as.numeric(metadata$C6646B_Serum1_hsCRP.x)<quantile(as.numeric(metadata$C6646B_Serum1_hsCRP.x), prob=0.5, na.rm=T) )] <- 1
metadata$crp_quartiles[which(as.numeric(metadata$C6646B_Serum1_hsCRP.x)>=quantile(as.numeric(metadata$C6646B_Serum1_hsCRP.x), prob=0.5, na.rm=T)&as.numeric(metadata$C6646B_Serum1_hsCRP.x)<quantile(as.numeric(metadata$C6646B_Serum1_hsCRP.x), prob=0.75, na.rm=T) )] <- 2
metadata$crp_quartiles[which(as.numeric(metadata$C6646B_Serum1_hsCRP.x)>=quantile(as.numeric(metadata$C6646B_Serum1_hsCRP.x), prob=0.75, na.rm=T) )] <- 3
metadata <- apply_labels(metadata, 
                         crp_quartiles	=	"HsCRP variable based on quartiles of C6646B_Serum1_hsCRP.x",
                         crp_quartiles=c("Q1"=0, "Q2"=1,"Q3"=2,"Q4"=3))
# S-BIL - bilirubin
metadata$bilirubin_quartiles <- NA
metadata$bilirubin_quartiles[which(metadata$C6646B_serum_gel11_12<quantile(metadata$C6646B_serum_gel11_12, prob=0.25, na.rm=T) )] <- 0
metadata$bilirubin_quartiles[which(metadata$C6646B_serum_gel11_12>=quantile(metadata$C6646B_serum_gel11_12, prob=0.25, na.rm=T)&metadata$C6646B_serum_gel11_12<quantile(metadata$C6646B_serum_gel11_12, prob=0.5, na.rm=T) )] <- 1
metadata$bilirubin_quartiles[which(metadata$C6646B_serum_gel11_12>=quantile(metadata$C6646B_serum_gel11_12, prob=0.5, na.rm=T)&metadata$C6646B_serum_gel11_12<quantile(metadata$C6646B_serum_gel11_12, prob=0.75, na.rm=T) )] <- 2
metadata$bilirubin_quartiles[which(metadata$C6646B_serum_gel11_12>=quantile(metadata$C6646B_serum_gel11_12, prob=0.75, na.rm=T) )] <- 3
metadata <- apply_labels(metadata, 
                         bilirubin_quartiles	=	"Bilirubin variable based on quartiles of C6646B_serum_gel11_12",
                         bilirubin_quartiles=c("Q1"=0, "Q2"=1,"Q3"=2,"Q4"=3))
# S-TESTO - testosterone
metadata$testo_quartiles <- NA
metadata$testo_quartiles[which(metadata$C6646B_serum1_01.x<quantile(metadata$C6646B_serum1_01.x, prob=0.25, na.rm=T) )] <- 0
metadata$testo_quartiles[which(metadata$C6646B_serum1_01.x>=quantile(metadata$C6646B_serum1_01.x, prob=0.25, na.rm=T)&metadata$C6646B_serum1_01.x<quantile(metadata$C6646B_serum1_01.x, prob=0.5, na.rm=T) )] <- 1
metadata$testo_quartiles[which(metadata$C6646B_serum1_01.x>=quantile(metadata$C6646B_serum1_01.x, prob=0.5, na.rm=T)&metadata$C6646B_serum1_01.x<quantile(metadata$C6646B_serum1_01.x, prob=0.75, na.rm=T) )] <- 2
metadata$testo_quartiles[which(metadata$C6646B_serum1_01.x>=quantile(metadata$C6646B_serum1_01.x, prob=0.75, na.rm=T) )] <- 3
metadata <- apply_labels(metadata, 
                         testo_quartiles	=	"Testosterone variable based on quartiles of C6646B_serum1_01.x",
                         testo_quartiles=c("Q1"=0, "Q2"=1,"Q3"=2,"Q4"=3))
# S-SHBG - sex hormone-binding globulin
metadata$shbg_quartiles <- NA
metadata$shbg_quartiles[which(metadata$C6646B_serum1_02.x<quantile(metadata$C6646B_serum1_02.x, prob=0.25, na.rm=T) )] <- 0
metadata$shbg_quartiles[which(metadata$C6646B_serum1_02.x>=quantile(metadata$C6646B_serum1_02.x, prob=0.25, na.rm=T)&metadata$C6646B_serum1_02.x<quantile(metadata$C6646B_serum1_02.x, prob=0.5, na.rm=T) )] <- 1
metadata$shbg_quartiles[which(metadata$C6646B_serum1_02.x>=quantile(metadata$C6646B_serum1_02.x, prob=0.5, na.rm=T)&metadata$C6646B_serum1_02.x<quantile(metadata$C6646B_serum1_02.x, prob=0.75, na.rm=T) )] <- 2
metadata$shbg_quartiles[which(metadata$C6646B_serum1_02.x>=quantile(metadata$C6646B_serum1_02.x, prob=0.75, na.rm=T) )] <- 3
metadata <- apply_labels(metadata, 
                         shbg_quartiles	=	"Sex hormone-binding globulin variable based on quartiles of C6646B_serum1_02.x",
                         shbg_quartiles=c("Q1"=0, "Q2"=1,"Q3"=2,"Q4"=3))
# S-SHBG - sex hormone-binding globulin
metadata$hba1c_quartiles <- NA
metadata$hba1c_quartiles[which(metadata$C6646B_blood_edta2_1<quantile(metadata$C6646B_blood_edta2_1, prob=0.25, na.rm=T) )] <- 0
metadata$hba1c_quartiles[which(metadata$C6646B_blood_edta2_1>=quantile(metadata$C6646B_blood_edta2_1, prob=0.25, na.rm=T)&metadata$C6646B_blood_edta2_1<quantile(metadata$C6646B_blood_edta2_1, prob=0.5, na.rm=T) )] <- 1
metadata$hba1c_quartiles[which(metadata$C6646B_blood_edta2_1>=quantile(metadata$C6646B_blood_edta2_1, prob=0.5, na.rm=T)&metadata$C6646B_blood_edta2_1<quantile(metadata$C6646B_blood_edta2_1, prob=0.75, na.rm=T) )] <- 2
metadata$hba1c_quartiles[which(metadata$C6646B_blood_edta2_1>=quantile(metadata$C6646B_blood_edta2_1, prob=0.75, na.rm=T) )] <- 3
metadata <- apply_labels(metadata, 
                         hba1c_quartiles	=	"Glycated haemoglobin variable based on quartiles of C6646B_blood_edta2_1",
                         hba1c_quartiles=c("Q1"=0, "Q2"=1,"Q3"=2,"Q4"=3))

# for the confounders
# gender
metadata$gender <- as.factor(metadata$gender)
# smoking status
metadata$smoking <- as.factor(metadata$smoke)

# remove all the participants that don't occur in the filtered OTU table
metadata.filtered <- metadata
metadata.filtered = metadata.filtered[row.names(metadata.filtered) %in% row.names(OTU.clean.filtered),]

# preparing metadata for analyses
# some analyses (permanova for example) cannot handle or account for NAs or blank. 
# subset to only samples with complete metadata before running the analyses if necessary.
meta.permanova <- metadata.filtered[!is.na(metadata.filtered$BMI_46_M_cat),]
meta.permanova <- meta.permanova[!is.na(meta.permanova$smoking),]
meta.permanova <- meta.permanova[!is.na(meta.permanova$gender),]
meta.permanova <- meta.permanova[!is.na(meta.permanova$homa_quartiles),]
meta.permanova <- meta.permanova[!is.na(meta.permanova$visceral_quartiles ),]
meta.permanova <- meta.permanova[!is.na(meta.permanova$bilirubin_quartiles  ),]
meta.permanova <- meta.permanova[!is.na(meta.permanova$testo_quartiles   ),]
meta.permanova <- meta.permanova[!is.na(meta.permanova$shbg_quartiles    ),]
meta.permanova <- meta.permanova[!is.na(meta.permanova$hba1c_quartiles    ),]
# remove all the participants that don't occur in the filtered metadata
OTU.permanova= OTU.clean.filtered[row.names(OTU.clean.filtered) %in% row.names(meta.permanova),]
# create a subset of the metadata by selecting columns of interest in the filtered metadata
meta<- dplyr::select(meta.permanova,homa_quartiles,visceral_quartiles,crp_quartiles,bilirubin_quartiles,testo_quartiles,shbg_quartiles,BMI_46_M_cat, hba1c_quartiles, gender
                     ,smoking, HOMAIR_46 , C6646C_BI_005 , C6646B_Serum1_hsCRP.x , C6646B_serum_gel11_12 , C6646B_serum1_01.x , C6646B_serum1_02.x, C6646B_blood_edta2_1 )
meta$homa <- as.factor(meta$homa_quartiles)
meta$visceral <- as.factor(meta$visceral_quartiles)
meta$crp <- as.factor(meta$crp_quartiles)
meta$bilirubin <- as.factor(meta$bilirubin_quartiles)
meta$testo <- as.factor(meta$testo_quartiles)
meta$shbg <- as.factor(meta$shbg_quartiles)
meta$hba1c <- as.factor(meta$hba1c_quartiles)

# crp has a lot more missing values so treating it apart
meta.permanova.crp <- meta.permanova[!is.na(meta.permanova$crp_quartiles),]
OTU.permanova.crp= OTU.clean.filtered[row.names(OTU.clean.filtered) %in% row.names(meta.permanova.crp),]
meta.crp<- dplyr::select(meta.permanova.crp,homa_quartiles,visceral_quartiles,crp_quartiles,bilirubin_quartiles,testo_quartiles,shbg_quartiles,BMI_46_M_cat,gender
                         ,smoking, HOMAIR_46 , C6646C_BI_005 , C6646B_Serum1_hsCRP.x , C6646B_serum_gel11_12 , C6646B_serum1_01.x , C6646B_serum1_02.x, hba1c_quartiles, C6646B_blood_edta2_1 )
meta.crp$homa <- as.factor(meta.crp$homa_quartiles)
meta.crp$visceral <- as.factor(meta.crp$visceral_quartiles)
meta.crp$crp <- as.factor(meta.crp$crp_quartiles)
meta.crp$bilirubin <- as.factor(meta.crp$bilirubin_quartiles)
meta.crp$testo <- as.factor(meta.crp$testo_quartiles)
meta.crp$shbg <- as.factor(meta.crp$shbg_quartiles)
meta.crp$hba1c <- as.factor(meta.crp$hba1c_quartiles)

# ------------------------------------------- BUILDING PHYLOSEQ OBJECTS ------------------------------------------- 

# phyloseq object 1 (all except crp)
OTU.UF.perma = otu_table(as.matrix(OTU.permanova), taxa_are_rows=FALSE)
tax.permanova = tax.clean[row.names(tax.clean) %in% colnames(OTU.permanova),]
tax.UF.perma = tax_table(as.matrix(tax.permanova))
meta.UF.perma = sample_data(meta)
physeq.permanova = phyloseq(OTU.UF.perma, tax.UF.perma, meta.UF.perma)
phys = merge_phyloseq(physeq.permanova, NJ.tree)
phys
# summarize the contents of a phyloseq object
summarize_phyloseq(phys)

# phyloseq object 2 (crp)
OTU.UF.perma.crp = otu_table(as.matrix(OTU.permanova.crp), taxa_are_rows=FALSE)
tax.permanova.crp = tax.clean[row.names(tax.clean) %in% colnames(OTU.permanova.crp),]
tax.UF.perma.crp = tax_table(as.matrix(tax.permanova.crp))
meta.UF.perma.crp = sample_data(meta.crp)
physeq.permanova.crp = phyloseq(OTU.UF.perma.crp, tax.UF.perma.crp, meta.UF.perma.crp)
phys.crp = merge_phyloseq(physeq.permanova.crp, NJ.tree)
rank_names(phys.crp)
table(tax_table(phys.crp)[, "Phylum"], exclude = NULL)
phys.crp <- subset_taxa(phys.crp, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
prevdf.crp = apply(X = otu_table(phys.crp),
                   MARGIN = ifelse(taxa_are_rows(phys.crp), yes = 1, no = 2),
                   FUN = function(x){sum(x > 0)})
prevdf.crp = data.frame(Prevalence = prevdf.crp,
                        TotalAbundance = taxa_sums(phys.crp),
                        tax_table(phys.crp))
plyr::ddply(prevdf.crp, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
filterPhyla = c("p__Acidobacteria", "p__Caldithrix","p__Euryarchaeota","p__Synergistetes","p__TM7","p__Elusimicrobia")
ps.crp = subset_taxa(phys.crp, !Phylum %in% filterPhyla)
prevdf1.crp = subset(prevdf.crp, Phylum %in% get_taxa_unique(ps.crp, "Phylum"))
prevalenceThreshold = 0.1 * nsamples(ps.crp)
keepTaxa.crp = rownames(prevdf1.crp)[(prevdf1.crp$Prevalence >= prevalenceThreshold)]
ps2.crp = prune_taxa(keepTaxa.crp, ps.crp)
length(get_taxa_unique(ps.crp, taxonomic.rank = "Genus"))
ps1.crp = tax_glom(ps.crp, "Genus", NArm = TRUE)
psra.crp = transform_sample_counts(ps.crp, function(x){x / sum(x)})
ps.crp <- prune_samples(sample_sums(ps.crp) > 1000, ps.crp)
ps.rarified.crp <- rarefy_even_depth(ps.crp, sample.size = 29000)
OTU1.crp=as(otu_table(ps.crp), "matrix")
OTUdf.crp = as.data.frame(OTU1.crp)
bi.dist.crp=vegdist(OTU1.crp, distance="binomial")
bc.dist.crp=vegdist(OTU1.crp, distance="bray")
jc.dist.crp=vegdist(OTU1.crp, distance="jaccard", binary = TRUE)
wUF.dist.crp = UniFrac(ps.crp, weighted=TRUE, normalized=TRUE)
uwUF.dist.crp = UniFrac(ps.crp, weighted=FALSE, normalized=TRUE)

# ------------------------------------------- FILTERING ------------------------------------------- 

# taxonomic filtering
# show available ranks in the dataset
rank_names(phys)
# create table, number of features for each phyla
table(tax_table(phys)[, "Phylum"], exclude = NULL)
# the following ensures that features with ambiguous phylum annotation are also removed 
phys <- subset_taxa(phys, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
# compute prevalence (define here as the number of samples in which a taxon appears at least once) of each feature, store as data.frame
prevdf = apply(X = otu_table(phys),
               MARGIN = ifelse(taxa_are_rows(phys), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(phys),
                    tax_table(phys))
# compute the total and average prevalences of the features in each phylum
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
# define phyla to filter
filterPhyla = c("p__Acidobacteria", "p__Caldithrix","p__Euryarchaeota","p__Synergistetes","p__TM7","p__Elusimicrobia")
# filter entries with unidentified phylum
ps = subset_taxa(phys, !Phylum %in% filterPhyla)

# prevalence filtering
# subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
# include a guess for parameter
geom_hline(yintercept = 0.01, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
facet_wrap(~Phylum) + theme(legend.position="none")
# define prevalence threshold as 10% of total samples
prevalenceThreshold = 0.1 * nsamples(ps)
prevalenceThreshold
# execute prevalence filter, using `prune_taxa()` function ; e.g. present in ???10% of the total number of samples
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps)

# abundance value transformation
# begin by defining a custom plot function, plot_abundance(), that uses phyloseq's function to define a relative abundance graphic
plot_abundance = function(physeq,title = "",
                          Facet = "Order", Color = "Phylum"){
  # arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(physeq, Phylum %in% c("p__Firmicutes"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x= "gender",y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}
# the transformation in this case converts the counts from each sample into their frequencies, often referred to as proportions or relative abundances
# transform to relative abundance
psra = transform_sample_counts(ps, function(x){x / sum(x)})
plotBefore = plot_abundance(ps,"")
plotAfter = plot_abundance(psra,"")
# combine each plot into one graphic
# the graphic shows the comparison of original abundances (top panel) and relative abundances (lower)
grid.arrange(nrow = 2,  plotBefore, plotAfter)

# rarefaction
# remove samples with fewer than 1000 reads
ps <- prune_samples(sample_sums(ps) > 1000, ps)
# subsampling of the read counts of each sample to a common read depth
ps.rarified <- rarefy_even_depth(ps, sample.size = 29000)
OTU1=as(otu_table(ps.rarified), "matrix")
OTUdf = as.data.frame(OTU1)

# ------------------------------------------- GENERATING DISTANCES AND DIVERSITY METRICS ------------------------------------------- 

# calculate distance and save as a matrix
bi.dist=vegdist(OTU1, distance="binomial")
bc.dist=vegdist(OTU1, distance="bray")
jc.dist=vegdist(OTU1, distance="jaccard", binary = TRUE)
wUF.dist = UniFrac(ps.rarified, weighted=TRUE, normalized=TRUE)
uwUF.dist = UniFrac(ps.rarified, weighted=FALSE, normalized=TRUE)

# create alpha diversity metrics
tab <- alpha(ps.rarified, index = "all")
htmlTable(head(tab))
# merge alpha diversity metrics with metadata
meta=data.frame(sample_data(ps.rarified))
tab$ID=row.names(tab)
meta$ID=row.names(meta)
meta<- right_join(x = meta, y = tab, by = "ID")
OTUdf1 = OTUdf
OTUdf1$ID=row.names(OTUdf1)
meta_asv<- right_join(x = meta, y = OTUdf1, by = "ID")

# ------------------------------------------- ANALYSES ------------------------------------------- 

# study sample characteristics
table(meta$homa_quartiles, meta$BMI_46_M_cat)

describe(meta$BMI_46_M_cat)

describe(meta$gender[meta$BMI_46_M_cat==1])
describe(meta$gender[meta$BMI_46_M_cat==2])
describe(meta$gender[meta$BMI_46_M_cat==3])

describe(meta$smoking[meta$BMI_46_M_cat==1])
describe(meta$smoking[meta$BMI_46_M_cat==2])
describe(meta$smoking[meta$BMI_46_M_cat==3])

describe(meta$HOMAIR_46[meta$BMI_46_M_cat==1])
sd(meta$HOMAIR_46[meta$BMI_46_M_cat==1])
describe(meta$HOMAIR_46[meta$BMI_46_M_cat==2])
sd(meta$HOMAIR_46[meta$BMI_46_M_cat==2])
describe(meta$HOMAIR_46[meta$BMI_46_M_cat==3])
sd(meta$HOMAIR_46[meta$BMI_46_M_cat==3])

describe(as.numeric(meta$C6646B_Serum1_hsCRP.x[meta$BMI_46_M_cat==1]))
sd(as.numeric(meta$C6646B_Serum1_hsCRP.x[meta$BMI_46_M_cat==1]), na.rm = T)
describe(as.numeric(meta$C6646B_Serum1_hsCRP.x[meta$BMI_46_M_cat==2]))
sd(as.numeric(meta$C6646B_Serum1_hsCRP.x[meta$BMI_46_M_cat==2]), na.rm = T)
describe(as.numeric(meta$C6646B_Serum1_hsCRP.x[meta$BMI_46_M_cat==3]))
sd(as.numeric(meta$C6646B_Serum1_hsCRP.x[meta$BMI_46_M_cat==3]), na.rm = T)

describe(meta$C6646B_serum_gel11_12[meta$BMI_46_M_cat==1])
sd(meta$C6646B_serum_gel11_12[meta$BMI_46_M_cat==1])
describe(meta$C6646B_serum_gel11_12[meta$BMI_46_M_cat==2])
sd(meta$C6646B_serum_gel11_12[meta$BMI_46_M_cat==2])
describe(meta$C6646B_serum_gel11_12[meta$BMI_46_M_cat==3])
sd(meta$C6646B_serum_gel11_12[meta$BMI_46_M_cat==3])

describe(meta$C6646B_serum1_01.x[meta$BMI_46_M_cat==1])
sd(meta$C6646B_serum1_01.x[meta$BMI_46_M_cat==1])
describe(meta$C6646B_serum1_01.x[meta$BMI_46_M_cat==2])
sd(meta$C6646B_serum1_01.x[meta$BMI_46_M_cat==2])
describe(meta$C6646B_serum1_01.x[meta$BMI_46_M_cat==3])
sd(meta$C6646B_serum1_01.x[meta$BMI_46_M_cat==3])

describe(meta$C6646B_serum1_02.x[meta$BMI_46_M_cat==1])
sd(meta$C6646B_serum1_02.x[meta$BMI_46_M_cat==1])
describe(meta$C6646B_serum1_02.x[meta$BMI_46_M_cat==2])
sd(meta$C6646B_serum1_02.x[meta$BMI_46_M_cat==2])
describe(meta$C6646B_serum1_02.x[meta$BMI_46_M_cat==3])
sd(meta$C6646B_serum1_02.x[meta$BMI_46_M_cat==3])

describe(meta$C6646B_blood_edta2_1[meta$BMI_46_M_cat==1])
sd(meta$C6646B_blood_edta2_1[meta$BMI_46_M_cat==1])
describe(meta$C6646B_blood_edta2_1[meta$BMI_46_M_cat==2])
sd(meta$C6646B_blood_edta2_1[meta$BMI_46_M_cat==2])
describe(meta$C6646B_blood_edta2_1[meta$BMI_46_M_cat==3])
sd(meta$C6646B_blood_edta2_1[meta$BMI_46_M_cat==3])

# correlation
# prepare the data
meta.cor<- dplyr::select(meta.permanova, HOMAIR_46 , C6646C_BI_005 , C6646B_Serum1_hsCRP.x , C6646B_serum_gel11_12 , C6646B_serum1_01.x ,
                         C6646B_serum1_02.x, C6646B_blood_edta2_1, BMI_46_M )
str(meta.cor)
meta.cor$BMI <- meta.cor$BMI_46_M
meta.cor$HOMA <- meta.cor$HOMAIR_46
meta.cor$Testosterone  <- meta.cor$C6646B_serum1_01.x
meta.cor$SHBG <- meta.cor$C6646B_serum1_02.x
meta.cor$HbA1c <- meta.cor$C6646B_blood_edta2_1
meta.cor$CRP <- as.numeric(meta.cor$C6646B_Serum1_hsCRP.x)
meta.cor$Bilirubin <- as.numeric(meta.cor$C6646B_serum_gel11_12)
meta.cor<- dplyr::select(meta.cor, BMI, HOMA, Testosterone, SHBG, HbA1c, CRP, Bilirubin)
head(meta.cor)
# calculate correlation matrix
cormat <- cor(meta.cor, use="pairwise.complete.obs")
head(cormat)
melted_cormat <- melt(cormat)
head(melted_cormat)
# upper triangle
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
upper_tri <- get_upper_tri(cormat)
# reorganized heatmap
reorder_cormat <- function(cormat){
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
print(ggheatmap)

ggheatmap + 
  geom_text(aes(Var2, Var1, label = round(value, 2)), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

# alpha diversity
# check distribution to see if roughly normal, if yes, ANOVA (categorical) or t-test (continuous), if no Kruskal-Wallis (categorical) or Wilcoxon (continuous)
shapiro.test(meta$diversity_shannon)
shapiro.test(meta$observed)

# description of the potential counfounders
# BMI (shannon significant, observed significant)
kruskal.test(list(meta$diversity_shannon[as.factor(meta$BMI_46_M_cat==1)],meta$diversity_shannon[as.factor(meta$BMI_46_M_cat==2)],
                  meta$diversity_shannon[as.factor(meta$BMI_46_M_cat==3)]))
kruskal.test(list(meta$observed[as.factor(meta$BMI_46_M_cat==1)],meta$observed[as.factor(meta$BMI_46_M_cat==2)],
                  meta$observed[as.factor(meta$BMI_46_M_cat==3)]))
# smoking (shannon non significant, observed significant)
kruskal.test(list(meta$diversity_shannon[meta$smoking=="never smoker"], meta$diversity_shannon[meta$smoking=="former smoker"],
                  meta$diversity_shannon[meta$smoking=="current smoker"]))
kruskal.test(list(meta$observed[meta$smoking=="never smoker"], meta$observed[meta$smoking=="former smoker"],
                  meta$observed[meta$smoking=="current smoker"]))
# gender (shannon non significant, observed non significant)
kruskal.test(list(meta$diversity_shannon[meta$gender=="Male"], meta$diversity_shannon[meta$gender=="Female"]))
kruskal.test(list(meta$observed[meta$gender=="Male"], meta$observed[meta$gender=="Female"]))

# anova with categorical variables
# shannon
lm.shannon.homa=lm(diversity_shannon ~ homa + as.factor(BMI_46_M_cat) + smoking + gender, data=meta)
Anova(lm.shannon.homa)
summary(lm.shannon.homa)
lm.shannon.crp=lm(diversity_shannon ~ crp + as.factor(BMI_46_M_cat) + smoking + gender, data=meta)
Anova(lm.shannon.crp)
summary(lm.shannon.crp)
lm.shannon.bilirubin=lm(diversity_shannon ~ bilirubin + as.factor(BMI_46_M_cat) + smoking + gender, data=meta)
Anova(lm.shannon.bilirubin)
summary(lm.shannon.bilirubin)
lm.shannon.testo=lm(diversity_shannon ~ testo + as.factor(BMI_46_M_cat) + smoking + gender, data=meta)
Anova(lm.shannon.testo)
summary(lm.shannon.testo)
lm.shannon.shbg=lm(diversity_shannon ~ shbg + as.factor(BMI_46_M_cat) + smoking + gender, data=meta)
Anova(lm.shannon.shbg)
summary(lm.shannon.shbg)
lm.shannon.hba1c=lm(diversity_shannon ~ hba1c + as.factor(BMI_46_M_cat) + smoking + gender, data=meta)
Anova(lm.shannon.hba1c)
summary(lm.shannon.hba1c)

# observed ASVs
lm.observed.homa=lm(observed ~ homa + as.factor(BMI_46_M_cat) + smoking + gender, data=meta)
Anova(lm.observed.homa)
summary(lm.observed.homa)
lm.observed.crp=lm(observed ~ crp + as.factor(BMI_46_M_cat) + smoking + gender, data=meta)
Anova(lm.observed.crp)
summary(lm.observed.crp)
lm.observed.bilirubin=lm(observed ~ bilirubin + as.factor(BMI_46_M_cat) + smoking + gender, data=meta)
Anova(lm.observed.bilirubin)
summary(lm.observed.bilirubin)
lm.observed.testo=lm(observed ~ testo + as.factor(BMI_46_M_cat) + smoking + gender, data=meta)
Anova(lm.observed.testo)
summary(lm.observed.testo)
lm.observed.shbg=lm(observed ~ shbg + as.factor(BMI_46_M_cat) + smoking + gender, data=meta)
Anova(lm.observed.shbg)
summary(lm.observed.shbg)
lm.observed.hba1c=lm(observed ~ hba1c + as.factor(BMI_46_M_cat) + smoking + gender, data=meta)
Anova(lm.observed.hba1c)
summary(lm.observed.hba1c)

# general linear model using the non-normal quasipoisson method with continuous variables
# shannon
qp.shannon.homa = glm(diversity_shannon ~ HOMAIR_46 + as.factor(BMI_46_M_cat) + gender + smoking, data=meta, family="quasipoisson")
summary(qp.shannon.homa)
qp.shannon.crp = glm(diversity_shannon ~ as.numeric(C6646B_Serum1_hsCRP.x)   + as.factor(BMI_46_M_cat) + gender + smoking, data=meta, family="quasipoisson")
summary(qp.shannon.crp)
qp.shannon.bili = glm(diversity_shannon ~ C6646B_serum_gel11_12    + as.factor(BMI_46_M_cat) + gender + smoking, data=meta, family="quasipoisson")
summary(qp.shannon.bili)
qp.shannon.testo = glm(diversity_shannon ~ C6646B_serum1_01.x     + as.factor(BMI_46_M_cat) + gender + smoking, data=meta, family="quasipoisson")
summary(qp.shannon.testo)
qp.shannon.shbg = glm(diversity_shannon ~ C6646B_serum1_02.x      + as.factor(BMI_46_M_cat) + gender + smoking, data=meta, family="quasipoisson")
summary(qp.shannon.shbg)
qp.shannon.hba1c = glm(diversity_shannon ~ C6646B_blood_edta2_1    + as.factor(BMI_46_M_cat) + gender + smoking, data=meta, family="quasipoisson")
summary(qp.shannon.hba1c)
qp.shannon.homa.crp = glm(diversity_shannon ~ HOMAIR_46 + as.numeric(C6646B_Serum1_hsCRP.x) + as.factor(BMI_46_M_cat) + gender + smoking, data=meta, family="quasipoisson")
summary(qp.shannon.homa.crp)

# observed ASVs
qp.observed.homa = glm(observed ~ HOMAIR_46 + as.factor(BMI_46_M_cat) + gender + smoking, data=meta, family="quasipoisson")
summary(qp.observed.homa)
qp.observed.crp = glm(observed ~ as.numeric(C6646B_Serum1_hsCRP.x)   + as.factor(BMI_46_M_cat) + gender + smoking, data=meta, family="quasipoisson")
summary(qp.observed.crp)
qp.observed.bili = glm(observed ~ C6646B_serum_gel11_12    + as.factor(BMI_46_M_cat) + gender + smoking, data=meta, family="quasipoisson")
summary(qp.observed.bili)
qp.observed.testo = glm(observed ~ C6646B_serum1_01.x     + as.factor(BMI_46_M_cat) + gender + smoking, data=meta, family="quasipoisson")
summary(qp.observed.testo)
qp.observed.shbg = glm(observed ~ C6646B_serum1_02.x      + as.factor(BMI_46_M_cat) + gender + smoking, data=meta, family="quasipoisson")
summary(qp.observed.shbg)
qp.observed.hba1c = glm(observed ~ C6646B_blood_edta2_1    + as.factor(BMI_46_M_cat) + gender + smoking, data=meta, family="quasipoisson")
summary(qp.observed.hba1c)
qp.observed.homa.crp = glm(observed ~ HOMAIR_46 + as.numeric(C6646B_Serum1_hsCRP.x) + as.factor(BMI_46_M_cat) + gender + smoking, data=meta, family="quasipoisson")
summary(qp.observed.homa.crp)

# beta diversity

# PERMANOVA

homa.bi = adonis(bi.dist ~ homa+ as.factor(BMI_46_M_cat)+ gender + smoking, data = meta, permutations = 1000, by = NULL)
homa.bi

# MiRKAT
# construct kernel matrix from distance matrix
k.bi = D2K(as.matrix(bi.dist))
k.jc = D2K(as.matrix(jc.dist))
k.wUF = D2K(as.matrix(wUF.dist))
k.uwUF = D2K(as.matrix(uwUF.dist))
# treat crp apart because of missing values
k.bi.crp = D2K(as.matrix(bi.dist.crp))
k.jc.crp = D2K(as.matrix(jc.dist.crp))
k.wUF.crp = D2K(as.matrix(wUF.dist.crp))
k.uwUF.crp = D2K(as.matrix(uwUF.dist.crp))

# with continuous variables
# homa
MiRKAT(y = meta$HOMAIR_46, Ks = k.bi, X = cbind(as.factor(meta$BMI_46_M_cat), meta$gender, meta$smoking)
       , out_type = "C", method = "davies", nperm=1000)
MiRKAT(y = meta$HOMAIR_46, Ks = k.jc, X = cbind(as.factor(meta$BMI_46_M_cat), meta$gender, meta$smoking)
       , out_type = "C", method = "davies", nperm=1000)
MiRKAT(y = meta$HOMAIR_46, Ks = k.wUF, X = cbind(as.factor(meta$BMI_46_M_cat), meta$gender, meta$smoking)
       , out_type = "C", method = "davies", nperm=1000)
MiRKAT(y = meta$HOMAIR_46, Ks = k.uwUF, X = cbind(as.factor(meta$BMI_46_M_cat), meta$gender, meta$smoking)
       , out_type = "C", method = "davies", nperm=1000)
# bilirubin
MiRKAT(y = meta$C6646B_serum_gel11_12  , Ks = k.bi, X = cbind(as.factor(meta$BMI_46_M_cat), meta$gender, meta$smoking)
       , out_type = "C", method = "davies", nperm=1000)
MiRKAT(y = meta$C6646B_serum_gel11_12  , Ks = k.jc, X = cbind(as.factor(meta$BMI_46_M_cat), meta$gender, meta$smoking)
       , out_type = "C", method = "davies", nperm=1000)
MiRKAT(y = meta$C6646B_serum_gel11_12  , Ks = k.wUF, X = cbind(as.factor(meta$BMI_46_M_cat), meta$gender, meta$smoking)
       , out_type = "C", method = "davies", nperm=1000)
MiRKAT(y = meta$C6646B_serum_gel11_12  , Ks = k.uwUF, X = cbind(as.factor(meta$BMI_46_M_cat), meta$gender, meta$smoking)
       , out_type = "C", method = "davies", nperm=1000)
# testo
MiRKAT(y = meta$C6646B_serum1_01.x  , Ks = k.bi, X = cbind(as.factor(meta$BMI_46_M_cat), meta$gender, meta$smoking)
       , out_type = "C", method = "davies", nperm=1000)
MiRKAT(y = meta$C6646B_serum1_01.x  , Ks = k.jc, X = cbind(as.factor(meta$BMI_46_M_cat), meta$gender, meta$smoking)
       , out_type = "C", method = "davies", nperm=1000)
MiRKAT(y = meta$C6646B_serum1_01.x  , Ks = k.wUF, X = cbind(as.factor(meta$BMI_46_M_cat), meta$gender, meta$smoking)
       , out_type = "C", method = "davies", nperm=1000)
MiRKAT(y = meta$C6646B_serum1_01.x  , Ks = k.uwUF, X = cbind(as.factor(meta$BMI_46_M_cat), meta$gender, meta$smoking)
       , out_type = "C", method = "davies", nperm=1000)
# shbg
MiRKAT(y = meta$C6646B_serum1_02.x  , Ks = k.bi, X = cbind(as.factor(meta$BMI_46_M_cat), meta$gender, meta$smoking)
       , out_type = "C", method = "davies", nperm=1000)
MiRKAT(y = meta$C6646B_serum1_02.x  , Ks = k.jc, X = cbind(as.factor(meta$BMI_46_M_cat), meta$gender, meta$smoking)
       , out_type = "C", method = "davies", nperm=1000)
MiRKAT(y = meta$C6646B_serum1_02.x  , Ks = k.wUF, X = cbind(as.factor(meta$BMI_46_M_cat), meta$gender, meta$smoking)
       , out_type = "C", method = "davies", nperm=1000)
MiRKAT(y = meta$C6646B_serum1_02.x  , Ks = k.uwUF, X = cbind(as.factor(meta$BMI_46_M_cat), meta$gender, meta$smoking)
       , out_type = "C", method = "davies", nperm=1000)
# crp
MiRKAT(y = meta.crp$C6646B_Serum1_hsCRP.x   , Ks = k.bi.crp, X = cbind(as.factor(meta.crp$BMI_46_M_cat), meta.crp$gender, meta.crp$smoking)
       , out_type = "C", method = "davies", nperm=1000)
MiRKAT(y = meta.crp$C6646B_Serum1_hsCRP.x  , Ks = k.jc.crp, X = cbind(as.factor(meta.crp$BMI_46_M_cat), meta.crp$gender, meta.crp$smoking)
       , out_type = "C", method = "davies", nperm=1000)
MiRKAT(y = meta.crp$C6646B_Serum1_hsCRP.x  , Ks = k.wUF.crp, X = cbind(as.factor(meta.crp$BMI_46_M_cat), meta.crp$gender, meta.crp$smoking)
       , out_type = "C", method = "davies", nperm=1000)
MiRKAT(y = meta.crp$C6646B_Serum1_hsCRP.x , Ks = k.uwUF.crp, X = cbind(as.factor(meta.crp$BMI_46_M_cat), meta.crp$gender, meta.crp$smoking)
       , out_type = "C", method = "davies", nperm=1000)
# hba1c
MiRKAT(y = meta$C6646B_blood_edta2_1  , Ks = k.bi, X = cbind(as.factor(meta$BMI_46_M_cat), meta$gender, meta$smoking)
       , out_type = "C", method = "davies", nperm=1000)
MiRKAT(y = meta$C6646B_blood_edta2_1  , Ks = k.jc, X = cbind(as.factor(meta$BMI_46_M_cat), meta$gender, meta$smoking)
       , out_type = "C", method = "davies", nperm=1000)
MiRKAT(y = meta$C6646B_blood_edta2_1  , Ks = k.wUF, X = cbind(as.factor(meta$BMI_46_M_cat), meta$gender, meta$smoking)
       , out_type = "C", method = "davies", nperm=1000)
MiRKAT(y = meta$C6646B_blood_edta2_1  , Ks = k.uwUF, X = cbind(as.factor(meta$BMI_46_M_cat), meta$gender, meta$smoking)
       , out_type = "C", method = "davies", nperm=1000)
# crp+homa
MiRKAT(y = meta.crp$C6646B_Serum1_hsCRP.x   , Ks = k.bi.crp, X = cbind(meta.crp$HOMAIR_46, as.factor(meta.crp$BMI_46_M_cat), meta.crp$gender, meta.crp$smoking)
       , out_type = "C", method = "davies", nperm=1000)
MiRKAT(y = meta.crp$C6646B_Serum1_hsCRP.x  , Ks = k.jc.crp, X = cbind(meta.crp$HOMAIR_46, as.factor(meta.crp$BMI_46_M_cat), meta.crp$gender, meta.crp$smoking)
       , out_type = "C", method = "davies", nperm=1000)
MiRKAT(y = meta.crp$C6646B_Serum1_hsCRP.x  , Ks = k.wUF.crp, X = cbind(meta.crp$HOMAIR_46, as.factor(meta.crp$BMI_46_M_cat), meta.crp$gender, meta.crp$smoking)
       , out_type = "C", method = "davies", nperm=1000)
MiRKAT(y = meta.crp$C6646B_Serum1_hsCRP.x , Ks = k.uwUF.crp, X = cbind(meta.crp$HOMAIR_46, as.factor(meta.crp$BMI_46_M_cat), meta.crp$gender, meta.crp$smoking)
       , out_type = "C", method = "davies", nperm=1000)
#homa+crp
MiRKAT(y = meta.crp$HOMAIR_46   , Ks = k.bi.crp, X =cbind(meta.crp$C6646B_Serum1_hsCRP.x, as.factor(meta.crp$BMI_46_M_cat), meta.crp$gender, meta.crp$smoking)
       , out_type = "C", method = "davies", nperm=1000)
MiRKAT(y = meta.crp$HOMAIR_46  , Ks = k.jc.crp, X = cbind(meta.crp$C6646B_Serum1_hsCRP.x, as.factor(meta.crp$BMI_46_M_cat), meta.crp$gender, meta.crp$smoking)
       , out_type = "C", method = "davies", nperm=1000)
MiRKAT(y = meta.crp$HOMAIR_46  , Ks = k.wUF.crp, X = cbind(meta.crp$C6646B_Serum1_hsCRP.x, as.factor(meta.crp$BMI_46_M_cat), meta.crp$gender, meta.crp$smoking)
       , out_type = "C", method = "davies", nperm=1000)
MiRKAT(y = meta.crp$HOMAIR_46 , Ks = k.uwUF.crp, X = cbind(meta.crp$C6646B_Serum1_hsCRP.x, as.factor(meta.crp$BMI_46_M_cat), meta.crp$gender, meta.crp$smoking)
       , out_type = "C", method = "davies", nperm=1000)

# ASVs
# there are a number of way to decrease the number of OTUs you're looking at:
# 1. don't use OTUs, add together genus or family groups and test if all or some of these taxa differ across variables of interest
# 2. apply an abundance cutoff such as only looking at OTUs/taxa that are at least 1% abundance in at least one sample
# 3. apply a frequency cutoff such as only looking at OTUs/taxa that occur in at least 10% of samples
# 4. combine 2 and 3

# gamlss package using zero-inflated beta distribution 

# agglomerate taxa, for generating relative abundance tables only /!\
# taxonomic agglomeration groups all the "leaves" in the hierarchy that descend from the user-prescribed agglomerating rank, this is sometimes called 'glomming'
# how many genera would be present after filtering?
ps.ra = tax_glom(psra, "Genus", NArm = TRUE)
ps.ra
# replace actual sequences with OTU ids
taxa_names(ps.ra)
n_seqs <- seq(ntaxa(ps.ra))
len_n_seqs <- nchar(max(n_seqs))
taxa_names(ps.ra) <- paste("Seq", formatC(n_seqs, 
                                          width = len_n_seqs, 
                                          flag = "0"), sep = "_")
taxa_names(ps.ra)
# generate a vector containing the full taxonomy path for all OTUs
wholetax <- do.call(paste, c(as.data.frame(tax_table(ps.ra))
                             [c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")], 
                             sep = "."))  # to distinguish from "_" within tax ranks
# turn the otu_table into a data.frame
otu_export <- as.data.frame(otu_table(ps.ra))
tmp <- names(otu_export)
tmp
# paste wholetax and OTU_ids together
for(i in 1:length(tmp)){
  names(tmp)[i] = paste(wholetax[i], tmp[i], sep = "__")
}
# overwrite old names
names(otu_export) <- names(tmp)
head(otu_export)[5]

metara=data.frame(sample_data(ps.ra))
metara$ID=row.names(metara)
otu_export$ID=row.names(otu_export)
metara<- right_join(x = metara, y = otu_export, by = "ID")
row.names(metara)=metara$ID
metara$bmi = as.factor(metara$BMI_46_M_cat)
otu_export$ID<-NULL
names(otu_export)

# remove ASVs with only 0
which(colSums(otu_export)==0)
OTUdfra <- otu_export[,-(which(colSums(otu_export) == 0))]

# solve deviance problem
dev <- function(x){
  model1<-gamlss(OTUdfra[,x] ~ metara$bmi+ metara$gender+ metara$smoking, family=BEZI, nu.formula = ~ metara$bmi+ metara$gender+ metara$smoking, sigma.formula = ~ metara$bmi+ metara$gender+ metara$smoking)
  model2<-gamlss(OTUdfra[,x] ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking, family=BEZI, nu.formula = ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking, sigma.formula = ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking)
  D0 <- deviance(model1)
  D1 <- deviance(model2)
  chi <- D0-D1
}
nb.dev <- as.matrix(c(1:ncol(OTUdfra)))
dev.res.all <-apply(nb.dev, 1, dev)
which(dev.res.all < 0)
OTUdfra <- OTUdfra[,-37]

# checking dimensions of the dataframe (rows columns)
dim(metara)
dim(OTUdfra)

lr <- function(x){
  model1<-gamlss(OTUdfra[,x] ~ metara$bmi+ metara$gender+ metara$smoking, family=BEZI, nu.formula = ~ metara$bmi+ metara$gender+ metara$smoking, sigma.formula = ~ metara$bmi+ metara$gender+ metara$smoking)
  model2<-gamlss(OTUdfra[,x] ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking, family=BEZI, nu.formula = ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking, sigma.formula = ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking)
  lr.res<-LR.test(model1,model2, print=F)$p.val
}

nb <- as.matrix(c(1:ncol(OTUdfra)))

lr.res.all <-apply(nb, 1, lr)

p.fdr <- p.adjust(as.numeric(lr.res.all), "fdr")
which(p.fdr < 0.05)
homa.fdr <-OTUdfra[,c(30 , 54 , 79 , 86 , 94 , 96 ,104 ,105, 113, 129 ,140)]
names(homa.fdr)

# homa
gamlss.homa1<-gamlss(homa.fdr[,1] ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking, 
                     family=BEZI, 
                     nu.formula = ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking, 
                     sigma.formula = ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking)
gamlss.homa1.sum<-as.data.frame(summary(gamlss.homa1))
gamlss.homa1.sum$fdr <- p.adjust(gamlss.homa1.sum[,4], "fdr")
gamlss.homa1.sum

gamlss.homa2<-gamlss(homa.fdr[,2] ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking, 
                     family=BEZI, 
                     nu.formula = ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking, 
                     sigma.formula = ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking)
gamlss.homa2.sum<-as.data.frame(summary(gamlss.homa2))
gamlss.homa2.sum$fdr <- p.adjust(gamlss.homa2.sum[,4], "fdr")
gamlss.homa2.sum

gamlss.homa3<-gamlss(homa.fdr[,3] ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking, 
                     family=BEZI, 
                     nu.formula = ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking, 
                     sigma.formula = ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking)
gamlss.homa3.sum<-as.data.frame(summary(gamlss.homa3))
gamlss.homa3.sum$fdr <- p.adjust(gamlss.homa3.sum[,4], "fdr")
gamlss.homa3.sum

gamlss.homa4<-gamlss(homa.fdr[,4] ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking, 
                     family=BEZI, 
                     nu.formula = ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking, 
                     sigma.formula = ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking)
gamlss.homa4.sum<-as.data.frame(summary(gamlss.homa4))
gamlss.homa4.sum$fdr <- p.adjust(gamlss.homa4.sum[,4], "fdr")
gamlss.homa4.sum

gamlss.homa5<-gamlss(homa.fdr[,5] ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking, 
                     family=BEZI, 
                     nu.formula = ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking, 
                     sigma.formula = ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking)
gamlss.homa5.sum<-as.data.frame(summary(gamlss.homa5))
gamlss.homa5.sum$fdr <- p.adjust(gamlss.homa5.sum[,4], "fdr")
gamlss.homa5.sum

gamlss.homa6<-gamlss(homa.fdr[,6] ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking, 
                     family=BEZI, 
                     nu.formula = ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking, 
                     sigma.formula = ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking)
gamlss.homa6.sum<-as.data.frame(summary(gamlss.homa6))
gamlss.homa6.sum$fdr <- p.adjust(gamlss.homa6.sum[,4], "fdr")
gamlss.homa6.sum

gamlss.homa7<-gamlss(homa.fdr[,7] ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking, 
                     family=BEZI, 
                     nu.formula = ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking, 
                     sigma.formula = ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking)
gamlss.homa7.sum<-as.data.frame(summary(gamlss.homa7))
gamlss.homa7.sum$fdr <- p.adjust(gamlss.homa7.sum[,4], "fdr")
gamlss.homa7.sum

gamlss.homa8<-gamlss(homa.fdr[,8] ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking, 
                     family=BEZI, 
                     nu.formula = ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking, 
                     sigma.formula = ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking)
gamlss.homa8.sum<-as.data.frame(summary(gamlss.homa8))
gamlss.homa8.sum$fdr <- p.adjust(gamlss.homa8.sum[,4], "fdr")
gamlss.homa8.sum

gamlss.homa9<-gamlss(homa.fdr[,9] ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking, 
                     family=BEZI, 
                     nu.formula = ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking, 
                     sigma.formula = ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking)
gamlss.homa9.sum<-as.data.frame(summary(gamlss.homa9))
gamlss.homa9.sum$fdr <- p.adjust(gamlss.homa9.sum[,4], "fdr")
gamlss.homa9.sum

gamlss.homa10<-gamlss(homa.fdr[,10] ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking, 
                     family=BEZI, 
                     nu.formula = ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking, 
                     sigma.formula = ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking)
gamlss.homa10.sum<-as.data.frame(summary(gamlss.homa10))
gamlss.homa10.sum$fdr <- p.adjust(gamlss.homa10.sum[,4], "fdr")
gamlss.homa10.sum

gamlss.homa11<-gamlss(homa.fdr[,11] ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking, 
                     family=BEZI, 
                     nu.formula = ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking, 
                     sigma.formula = ~ metara$HOMAIR_46+ metara$bmi+ metara$gender+ metara$smoking)
gamlss.homa11.sum<-as.data.frame(summary(gamlss.homa11))
gamlss.homa11.sum$fdr <- p.adjust(gamlss.homa11.sum[,4], "fdr")
gamlss.homa11.sum
