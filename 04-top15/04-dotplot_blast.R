library(zoo)
library(tidyverse)
library(stringr)
library(patchwork)
library(Rsamtools)
library(txdbmaker)
library(rtracklayer)
library(GenomicRanges)
library(Biostrings)

setwd("~/plantcad-zygophyllum-orthologs/")
source("R/extract-helpers.R")
source("04-top15/R/plot-helpers.R")


# settings ----------------------------------------------------------------

top15_dir <- "04-top15/out"
blastn_dir <- file.path(top15_dir, "blastn")
blastp_dir <- file.path(top15_dir, "blastp")
plot_dir <- file.path(top15_dir, "dot_plot")
dir.create(blastn_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(blastp_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)


# color blind friendly palette --------------------------------------------

palette <- palette.colors(palette = "Okabe-Ito", names = TRUE)
class_colors <- c(
  "low" = palette[["blue"]],
  "neutral" = palette[["gray"]],
  "high" = palette[["vermillion"]]
)


# load --------------------------------------------------------------------

tx_pair_keys <- read_tsv(
  "03-candidate-transcripts/out/conserved_go_annotated_braker_pairchecked.tsv"
  ) |>
  mutate(
    agic_tx = paste0(agic, ".1"),
    bgc_tx = paste0(braker_gene_id, ".t1")
  ) |>
  dplyr::select(agic_tx, bgc_tx)

at_tx_list <- readRDS(file = "04-top15/out/at_tx_list.rds")
zs_tx_list <- readRDS(file = "04-top15/out/zs_tx_list.rds")

busco_top15 <- read_tsv("04-top15/out/busco_top15.tsv")

# filter 42 genes to the 15 actually used
tx_pair_keys <- tx_pair_keys |>
  filter(
    agic_tx %in% paste0(busco_top15$agic, ".1"),
    bgc_tx %in% paste0(busco_top15$braker_gene_id, ".t1")
  )

at_genome_fa <- load_genome(
  genome_path = "data/A_thaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
)$fasta

zs_genome <- load_genome(
  genome_path = "data/Z_stapfii/Z_stapfii_genome.fa",
  gff_path = "data/Z_stapfii/braker.gff3"
)
zs_genome_fa <- zs_genome$fasta


at_tx_db <- makeTxDbFromGFF(
  "data/A_thaliana/Arabidopsis_thaliana.TAIR10.63.gff3",
  format = "gff3"
  )

# TxDb of the BRAKER gff: transcript names come from Name, not transcript_id
zs_gff <- zs_genome$gff
zs_gff$transcript_id <- NULL
zs_gff$Name <- ifelse(zs_gff$type == "mRNA", zs_gff$ID, NA_character_)
zs_tx_db <- makeTxDbFromGRanges(zs_gff)

at_cds_by_tx <- cdsBy(at_tx_db, by = "tx", use.names = TRUE)
zs_cds_by_tx <- cdsBy(zs_tx_db, by = "tx", use.names = TRUE)

# functions ---------------------------------------------------------------

process_pair <- function(agic_tx, bgc_tx) {
  pair_name <- paste0(agic_tx, "_with_", bgc_tx)
  message("processing ", pair_name)
  
  blast_dotplot(
    at_genome_fa = at_genome_fa,
    zs_genome_fa = zs_genome_fa,
    at_tx_list = at_tx_list,
    zs_tx_list = zs_tx_list,
    agic_tx = agic_tx,
    bgc_tx = bgc_tx,
    type = "blastn",
    outdir = blastn_dir,
    plot_dir = plot_dir
  )
  
  
  blast_dotplot(
    at_genome_fa = at_genome_fa,
    zs_genome_fa = zs_genome_fa,
    at_tx_list = at_tx_list,
    zs_tx_list = zs_tx_list,
    agic_tx = agic_tx,
    bgc_tx = bgc_tx,
    type = "blastp",
    outdir = blastp_dir,
    plot_dir = plot_dir,
    at_cds_by_tx = at_cds_by_tx,
    zs_cds_by_tx = zs_cds_by_tx
  )
  
}


# run ---------------------------------------------------------------------

walk2(tx_pair_keys$agic_tx, tx_pair_keys$bgc_tx, process_pair)
