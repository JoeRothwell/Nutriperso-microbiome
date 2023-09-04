# Data analysis for the Nutriperso 16s microbiome project

Nutriperso is a study of 147 cases of diabetes and 147 controls from the E3N cohort. Saliva samples, taken from participants at baseline, were analysed by 16s rRNA microbiome profiling. This repository is the codebase for data analysis.

This analysis of 16s rRNA data mainly uses the *Phyloseq* Bioconductor package, although other packages are used that do more specific things. The script nutriperso_data_prep.R creates the phyloseq objects from 16s data, assigned taxonomies and participant data. 16s data are available in both .otu and .asv format. For the OTU approach, a phylogenetic tree is created from aligned sequence data, and this is added to the phyloseq object.

Other scripts perform more downstream analysis: descriptive analysis of 16s data and calculation of alpha and beta diversity, calculation of enterotypes and epidemiological models. 
