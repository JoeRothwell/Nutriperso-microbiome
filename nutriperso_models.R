# Epidemiological models for diversity

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

varlist <- c("Tr_int1", "Tr_int2", "Tr_int3", "Tr_int4", "tabac", "Pb_dents1", "Pb_dents_dt1", 
             "dent_app_o", "dent_app_n", "dent_app_dt", "dent_app_com", "dent_app_bas", 
             "dent_app_haut","dent_app_tj_o", "dent_app_tj_n", "dent_app_sup_o", "dent_app_sup_n", 
             "RUN", "diabete_groupe", "imc_oms_salive_k", "diabete_groupe1")

dat <- dat %>% mutate(across((varlist), as.factor))

# Get participant data (sample data)
rownames(dat) <- str_remove(dat$ID, pattern = "-")
sampdata1 <- dat[colnames(otumat1), ]
sampledata <- sample_data(data.frame(sampdata1, row.names = sample_names(physeq),
                                     stringsAsFactors = F))


# Join alpha diversity data from nutriperso_16s_analysis.R. Make ID variables consistent
dat1 <- dat %>% mutate(ID = str_remove(ID, pattern = "-"))
div1 <- div %>% rownames_to_column() %>% rename(rowname = "ID")
divdat <- left_join(dat1, div1)

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
library(survival)
mod3 <- clogit(ct ~ scale(Shannon) + RUN + imcq10 + alcoolq8 + METsTotalQ8 + SUCRES + strata(matchnutpkt), divdat)
mod4 <- clogit(ct ~ scale(Simpson) + RUN + imcq10 + alcoolq8 + METsTotalQ8 + SUCRES + strata(matchnutpkt), divdat)

library(broom)
tidy(mod3, exponentiate = TRUE, conf.int = TRUE)
tidy(mod4, exponentiate = TRUE, conf.int = TRUE)


