library(tidyverse)
library(stringr)

setwd("~/plantcad-zygophyllum-orthologs/")
cwd_path <- file.path("05-additional-genes")

s9 <- read_csv("data/T_mongolica/liu2023_tableS9_TPS.csv")

writeLines(
  str_replace(s9$gene_id, "\\.\\d+", ""),
  file.path(cwd_path, "tps_ids.txt")
)
