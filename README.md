# Data analysis for the Nutriperso 16s microbiome project

Nutriperso is a study of 297 cases of diabetes and 297 controls from the E3N cohort who gave saliva samples at baseline. These saliva samples were analysed by 16s microbiome profiling. This repository is the codebase for data analysis.

This analysis mainly uses the *Phyloseq* Bioconductor package. The script nutriperso_data_prep.R creates the phyloseq objects, from 16s data, assigned taxonomies and participant data. 16s data are available in both .otu and .asv format. For the OTU approach, a phylogenetic tree is created from aligned sequence data, and this is added to the phyloseq object.

Other scripts perform more downstream analysis: descriptive analysis of 16s data and calculation of alpha and beta diversity, calculation of enterotypes and epidemiological models. 
