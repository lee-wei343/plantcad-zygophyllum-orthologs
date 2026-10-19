library(tidyverse)

setwd("~/plantcad-zygophyllum-orthologs/05-additional-genes/")

busco_Zsta <- data.table::fread("out/stapffii_proteom/run_embryophyta_odb10/full_table.tsv")

tps_orthos <- read_tsv("tps_orthos.tsv") |>
  mutate(Zsta_tx = paste0(Zsta, ".t1"))

tps_orthos$in_BUSCO <- tps_orthos$Zsta_tx %in% busco_Zsta$Sequence

write_tsv(tps_orthos, file = "tps_orthos_in_busco.tsv")
