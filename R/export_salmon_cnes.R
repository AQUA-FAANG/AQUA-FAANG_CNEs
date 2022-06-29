library(CNEr)
library(tidyverse)

cneFiles <- list.files("/mnt/orca/projects/AQUA-FAANG/species/CNEs/",
                       pattern = "RDS", full.names = T)
names(cneFiles) <- str_match(cneFiles, "Salmon_(.+)\\.RDS")[,2]

cneClass <- lapply(cneFiles, readRDS)
