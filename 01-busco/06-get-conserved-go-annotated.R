library(org.At.tair.db)
library(tidyverse)

keep_only_complete_genes <- function(tsv) {
  n_cols <- length(str_split_1(read_lines(tsv, skip = 3, n_max = 1), "\t"))
  cols <- if (n_cols == 10) {
    c("busco_id", "status", "sequence", "gene_start", "gene_end", "strand",
      "score", "length", "orthodb_url", "description")
  } else {
    c("busco_id", "status", "sequence", "score", "length", "orthodb_url", "description")
  }
  read_tsv(
    tsv,
    comment = "#",
    col_names = cols,
  ) |>
    filter(status == "Complete") |>
    mutate(tair = str_remove(sequence, "\\.\\d+$"))
}

keep_only_go_and_complete_genes <- function(tsv) {
   
  busco <- keep_only_complete_genes(tsv)
  
  # take only the BUSCO genes with experimental data, which are at least somewhat known
  experimental_evidence <- c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP")
  evident_go_annotation <- AnnotationDbi::select(
    org.At.tair.db,
    keys = busco$tair,
    columns = c("SYMBOL", "GENENAME", "GO", "EVIDENCE"),
    keytype = "TAIR"
  ) |>
    filter(EVIDENCE %in% experimental_evidence) |>
    distinct(TAIR, SYMBOL, GENENAME) |>
    filter(!is.na(GENENAME))
  
  evident_go_busco <- busco |>
    inner_join(evident_go_annotation, by = c("tair" = "TAIR"))
  
 
  return(evident_go_busco)
}


A_thaliana <- keep_only_go_and_complete_genes(
  tsv = "../out/busco/busco_A_thaliana/run_embryophyta_odb10/full_table.tsv"
  )

O_sativa <- keep_only_complete_genes(
  tsv = "../out/busco/busco_O_sativa/run_embryophyta_odb10/full_table.tsv"
)

P_patens <- keep_only_complete_genes(
  tsv = "../out/busco/busco_P_patens/run_embryophyta_odb10/full_table.tsv"
)

Z_stapfii <- keep_only_complete_genes(
  tsv = "../out/busco/busco_Z_stapfii/run_embryophyta_odb10/full_table.tsv"
)

T_mongolica <- keep_only_complete_genes(
  tsv = "../out/busco/busco_T_mongolica/run_embryophyta_odb10/full_table.tsv"
)

conserved <- A_thaliana |>
  select(busco_id, description, tair, SYMBOL, GENENAME) |>
  semi_join(O_sativa, join_by(busco_id)) |>
  semi_join(P_patens, join_by(busco_id)) |>
  semi_join(T_mongolica, join_by(busco_id)) |>
  inner_join(
    Z_stapfii |>
      select(busco_id, contig = sequence, gene_start, gene_end, strand),
    join_by(busco_id)
  )
