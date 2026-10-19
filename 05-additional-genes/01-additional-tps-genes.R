

# setup -------------------------------------------------------------------

library(stringr)
library(txdbmaker)
library(tidyverse)
setwd("~/plantcad-zygophyllum-orthologs/")
source("R/extract-helpers.R")
source("04-top15/R/plot-helpers.R")


# set ---------------------------------------------------------------------

s9 <- read_csv("data/T_mongolica/liu2023_tableS9_TPS.csv")

writeLines(str_replace(s9$gene_id, "\\.\\d+", ""), "tps_ids.txt")

s10 <- read_csv("data/T_mongolica/liu2023_tableS10_TAG.csv")

zs_v_tm <- read_tsv("02-orthofinder/Results_Sep29/Orthologues/Orthologues_Z_stapffii/Z_stapffii__v__T_mongolica.tsv")

tm_gff_path <- file.path("data", "T_mongolica", "T_mongolica_genes.gff3")

zs_v_tm_uniq <- zs_v_tm |>
  mutate(
    n_zs = if_else(Z_stapffii == "", 0, str_count(Z_stapffii, ",") + 1),
    n_tm = if_else(T_mongolica == "", 0, str_count(T_mongolica, ",") + 1)
  ) |>
  mutate(
    relation = case_when(
      n_zs == 0 | n_tm == 0 ~ "absent",
      n_zs == 1 & n_tm == 1 ~ "1:1",
      n_zs == 1 | n_tm == 1 ~ "1:many",
      .default = "many:many"
    )
  ) |>
  mutate(
    one_to_one = relation == "1:1"
  )


zs_v_tm_long <- zs_v_tm_uniq |>
  separate_longer_delim("Z_stapffii", delim = ", ") |>
  separate_longer_delim("T_mongolica", delim = ", ")

tps_orthos <- s9 |>
  mutate(gene_id = str_replace(gene_id, "\\.\\d+", "")) |>
  left_join(
    zs_v_tm_long,
    by = c("gene_id" = "T_mongolica"))

# read in BUSCO gene names from BUSCO proteom run to filter on non-BUSCO genes








# functions ---------------------------------------------------------------

add_metadata_col <- function(exons_by_tx, col, col_name, match_by) {
  idx <- match(match_by, names(exons_by_tx))
  if (anyNA(idx)) stop("some IDs in match_by are not in names(exons_by_tx)")
  
  mcols(exons_by_tx)[[col_name]] <- rep(NA, length(exons_by_tx))
  mcols(exons_by_tx)[[col_name]][idx] <- col
  exons_by_tx
}


# main --------------------------------------------------------------------

exons_by_tx <- get_exons_by_tx(gff_path = tm_gff_path)
exons_by_tx <- exons_by_tx[names(exons_by_tx) %in% s9$gene_id]

exons_by_tx <- add_metadata_col(exons_by_tx,
                                col      = s9$tps_subfamily,
                                col_name = "tps_subfamily",
                                match_by = s9$gene_id)


map(exons_by_tx, plot_transcript_simple)


# s9_sub <- s9 |>
#   dplyr::select(gene_id, name, tps_subfamily, chrom) |>
#   mutate(gene_id = str_replace(gene_id, "\\.[0-9]+", ""))
# 
# ids <- s9_sub$gene_id
# matches <- match(ids, tm_gff$ID)
# tm_gff$tps_subfamily <- NA
# tm_gff$tps_subfamily[matches] <- s9_sub$tps_subfamily
# 
# 
# 
# genes <- tm_gff[tm_gff$type == "gene"]
# 
# s9_sub <- s9 |>
#   select(gene_id, name, tps_subfamily, chrom) |>
#   mutate(gene_id = str_replace(gene_id, "\\.[0-9]+", ""))
