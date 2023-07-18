# Epidemiological models for diversity. First load phyloseq object
load("nutriperso_phyloseq.rda")

# Join alpha diversity data from nutriperso_16s_analysis.R. Make ID variables consistent
library(tidyverse)
library(phyloseq)
dat1 <- data.frame(sample_data(nutri)) %>% mutate(ID = str_remove(ID, pattern = "-"))

# Estimate alpha diversity for epidemiological models (based on diversity() from vegan) and join
div <- estimate_richness(nutri, measures = c("Shannon", "Simpson"))
div1 <- div %>% rownames_to_column(var = "ID")
divdat <- left_join(dat1, div1)


# Get enterotypes and add to diversity table
nutri.rel  <- transform_sample_counts(nutri, function(x) x / sum(x) )
jsdist <- phyloseq::distance(nutri.rel, method="jsd") # or JSD(physeq)
# Get clusters

# Optimal number of clusters is 2 (see enterotypes script)

#clust <- as.vector(pam(jsdist, 2, diss = TRUE)$clustering)
library(cluster)
clust <- pam(jsdist, 2, diss = TRUE)$clustering
#clust1 <- data.frame(clust) %>% rownames_to_column("ID")

divdat$enterotype <- clust
divdat %>% select(ID, Shannon, Simpson, enterotype) %>% write_csv("Nutriperso_diversity.csv")

# Hack to code case and control as 1 and 0
divdat$ct <- as.numeric(fct_rev(divdat$casnutpkt)) - 1

# Logistic model to get OR for diversity
# Unadjusted models
library(survival)
mod1 <- clogit(ct ~ scale(Shannon) + strata(matchnutpkt), divdat)
mod2 <- clogit(ct ~ scale(Simpson) + strata(matchnutpkt), divdat)
summary(mod1)

library(broom)
tidy(mod1, exponentiate = TRUE, conf.int = TRUE)
tidy(mod2, exponentiate = TRUE, conf.int = TRUE)

# Adjusted models (removed tabac as gave unstable coefficients)
mod3 <- clogit(ct ~ scale(Shannon) + RUN + imcq10 + alcoolq8 + METsTotalQ8 + SUCRES + strata(matchnutpkt), divdat)
mod4 <- clogit(ct ~ scale(Simpson) + RUN + imcq10 + alcoolq8 + METsTotalQ8 + SUCRES + strata(matchnutpkt), divdat)

tidy(mod3, exponentiate = TRUE, conf.int = TRUE)
tidy(mod4, exponentiate = TRUE, conf.int = TRUE)


# Enterotypes. First read in data calculated in enterotypes.R
nutr <- read.csv("Nutriperso_diversity.csv")


