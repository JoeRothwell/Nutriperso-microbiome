# Epidemiological models

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

# Join 2 sets of participant data
div1 <- div %>% rownames_to_column() %>% rename(rowname = "ID")
divdat <- left_join(dat, div)

# Get particpant data (sample data)
rownames(dat) <- str_remove(dat$ID, pattern = "-")
sampdata1 <- dat[colnames(otumat1), ]
sampledata <- sample_data(data.frame(sampdata1, row.names = sample_names(physeq),
                                     stringsAsFactors = F))
