library(tidyverse)

setwd("~/plantcad-zygophyllum-orthologs/")

msic <- read_csv("06-model-comparison/msic_long.csv", show_col_types = FALSE)
# -------------------------------------------------------------------------

get_blastp_by_tx <- function(dir_path) {
  blastp <- list.files(dir_path, pattern = "_with_.*\\.tsv$") |>
    set_names(basename) |>
    map(
      \(f) read_tsv(file.path(dir_path, f), show_col_types = FALSE)
    ) |>
    list_rbind(names_to = "file") |>
    mutate(
      at_tx = str_extract(string = file, pattern = "^[^_]+"),
      zs_tx = str_extract(string = file,
                          pattern = "_with_(.+)\\.tsv", group = 1)
    )
  
  best <- blastp |>
    slice_max(bitscore, n = 1, with_ties = FALSE, by = c(at_tx, zs_tx))
  
  pident_by_tx <- bind_rows(
    best |> select(tx = at_tx, pident),
    best |> select(tx = zs_tx, pident)
  )
}

top15 <- get_blastp_by_tx(file.path("04-top15", "out", "blastp"))
tps <- get_blastp_by_tx(file.path("05-additional-genes", "blastp"))

pident_by_tx <- row_bind(top15, tps)

msic_blastp <- msic |>
  left_join(pident_by_tx, by = "tx", relationship = "many-to-one")

write_csv(msic_blastp, file = "06-model-comparison/msic_blastp.csv")
