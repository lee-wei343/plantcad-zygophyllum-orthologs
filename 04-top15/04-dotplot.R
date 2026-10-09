library(zoo)
library(tidyverse)
library(stringr)
library(patchwork)
library(Rsamtools)

setwd("~/plantcad-zygophyllum-orthologs/")
source("04-top15/R/plot-helpers.R")


# settings ----------------------------------------------------------------

top15_dir <- "04-top15/out"
blast_dir <- file.path(top15_dir, "blastn")
plot_dir <- file.path(top15_dir, "dot_plot")
dir.create(blast_dir, showWarnings = FALSE, recursive = TRUE)
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
  select(agic_tx, bgc_tx)

at_tx_list <- readRDS(file = "04-top15/out/at_tx_list.rds")
zs_tx_list <- readRDS(file = "04-top15/out/zs_tx_list.rds")

# filter 42 genes to the 15 actually used
tx_pair_keys <- tx_pair_keys |>
  filter(
    agic_tx %in% names(at_file_ids_pcad1_index),
    bgc_tx %in% names(zs_file_ids_pcad1_index)
  )

at_genome_path <- "data/A_thaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
if (!file.exists(paste0(at_genome_path, ".fai"))) {
  indexFa(at_genome_path)
}
at_genome_fa <- FaFile(at_genome_path)

zs_genome_path <- "data/Z_stapfii/Z_stapfii_genome.fa"
if (!file.exists(paste0(zs_genome_path, ".fai"))) {
  indexFa(zs_genome_path)
}
zs_genome_fa <- FaFile(zs_genome_path)


# functions ---------------------------------------------------------------

process_pair <- function(agic_tx, bgc_tx) {
  pair_name <- paste0(agic_tx, "_with_", bgc_tx)
  message("processing ", pair_name)

  # sequences + blastn
  at_fa <- write_fasta(
    getSeq(
      at_genome_fa,
      at_tx_list[[agic_tx]]$current_tx
    ),
    blast_dir,
    agic_tx
  )
  zs_fa <- write_fasta(
    getSeq(
      zs_genome_fa,
      zs_tx_list[[bgc_tx]]$current_tx
    ),
    blast_dir,
    bgc_tx
  )
  
  blastn_tsv <- file.path(blast_dir, paste0(pair_name, ".tsv"))
  hsps <- run_blastn(
    at_fa,
    zs_fa,
    blastn_tsv
  )

  # plot
  # TODO: write plot_dotplot_marginal() in R/plot-helpers.R
  # p <- plot_dotplot_marginal(hsps, msic_at, msic_zs)
  # ggsave(file.path(plot_dir, paste0(pair_name, ".png")), p,
  #        dpi = 300, width = 10, height = 10)
  
  message("boop")
  
  a <- read_tsv(file = blastn_tsv)
  
}

a <- read_tsv(file = "04-top15/out/blastn/AT1G71220.1_with_g20167.t1.tsv")


# run ---------------------------------------------------------------------

walk2(tx_pair_keys$agic_tx[1], tx_pair_keys$bgc_tx[1], process_pair)
