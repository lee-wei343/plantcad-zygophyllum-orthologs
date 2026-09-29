library(org.At.tair.db)
library(tidyverse)

setwd("~/plantcad-zygophyllum-orthologs/")

keep_only_complete_genes <- function(tsv) {
  n_cols <- length(str_split_1(read_lines(tsv, skip = 3, n_max = 1), "\t"))
  # if protein - else genome
  cols <- if (n_cols == 10) {
    c(
      "busco_id", "status", "sequence", "gene_start", "gene_end", "strand",
      "score", "length", "orthodb_url", "description"
    )
  } else {
    c("busco_id", "status", "sequence", "score", "length", "orthodb_url", "description")
  }
  read_tsv(
    tsv,
    comment = "#",
    col_names = cols,
    show_col_types = FALSE
  ) |>
    filter(status %in% c("Complete", "Duplicated")) |>
    # if protein - else genome
    # agic = Arabidopsis Genome Initiative Code
    mutate(agic = if (n_cols == 10) paste(sequence, gene_start) else str_remove(sequence, "\\.\\d+$")) |>
    filter(n_distinct(agic) == 1, .by = busco_id) |>
    slice_max(score, n = 1, by = busco_id, with_ties = F)
}

keep_only_go_and_complete_genes <- function(tsv) {
  busco <- keep_only_complete_genes(tsv)

  # take only the BUSCO genes with experimental data, which are at least somewhat known
  experimental_evidence <- c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP")
  evident_go_annotation <- AnnotationDbi::select(
    org.At.tair.db,
    keys = busco$agic,
    columns = c("SYMBOL", "GENENAME", "GO", "EVIDENCE"),
    keytype = "TAIR"
  ) |>
    filter(EVIDENCE %in% experimental_evidence) |>
    filter(!is.na(GENENAME)) |>
    summarise(
      SYMBOL = paste(unique(SYMBOL), collapse = "/"),
      GENENAME = first(GENENAME),
      .by = TAIR
    )

  evident_go_busco <- busco |>
    inner_join(evident_go_annotation, by = c("agic" = "TAIR"))


  return(evident_go_busco)
}


A_thaliana <- keep_only_go_and_complete_genes(
  tsv = "01-busco/out/busco/busco_A_thaliana/run_embryophyta_odb10/full_table.tsv"
)

O_sativa <- keep_only_complete_genes(
  tsv = "01-busco/out/busco/busco_O_sativa/run_embryophyta_odb10/full_table.tsv"
)

P_patens <- keep_only_complete_genes(
  tsv = "01-busco/out/busco/busco_P_patens/run_embryophyta_odb10/full_table.tsv"
)

Z_stapfii <- keep_only_complete_genes(
  tsv = "01-busco/out/busco/busco_Z_stapfii/run_embryophyta_odb10/full_table.tsv"
)

T_mongolica <- keep_only_complete_genes(
  tsv = "01-busco/out/busco/busco_T_mongolica/run_embryophyta_odb10/full_table.tsv"
)

conserved <- A_thaliana |>
  select(busco_id, description, agic, SYMBOL, GENENAME) |>
  semi_join(O_sativa, join_by(busco_id)) |>
  semi_join(P_patens, join_by(busco_id)) |>
  semi_join(T_mongolica, join_by(busco_id)) |>
  inner_join(
    Z_stapfii |>
      select(busco_id, contig = sequence, gene_start, gene_end, strand),
    join_by(busco_id)
  )

stopifnot(!anyDuplicated(conserved$busco_id))


# pots <- c("AT1G71220", "AT5G54440", "AT5G15540", "AT4G14180", "AT3G24320", "AT3G48110", "AT4G02070", "AT3G20440", "AT1G44900", "AT4G17380", "AT1G11290", "AT4G01690", "AT1G67840", "AT5G01920", "AT5G05680", "AT5G54080", "AT5G41480", "AT3G59380", "AT1G64350", "AT3G47860", "AT5G22370", "AT4G17050", "AT1G21710", "AT4G25080", "AT3G05210", "AT2G48120", "AT4G29170")
# pets <- c("AT1G67840", "AT5G01920", "AT5G41480", "AT3G47860", "AT4G17050", "AT1G21710", "AT2G48120", "AT4G29170")


# conserved_well_annotated <- conserved |>
#   filter(agic %in% pots)
# 
# top_candidates <- conserved |>
#   filter(agic %in% pets)

write_tsv(conserved, "01-busco/out/conserved_go_annotated.tsv")
