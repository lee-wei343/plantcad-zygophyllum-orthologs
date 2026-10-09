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

at_pcad1_dir <- file.path(top15_dir, "a_thaliana_top15_PCAD1_results")
at_pcad2s_dir <- file.path(top15_dir, "a_thaliana_top15_PCAD2S_results")
zs_pcad1_dir  <- file.path(top15_dir, "z_stapffii_top15_PCAD1_results")
zs_pcad2s_dir <- file.path(top15_dir, "z_stapffii_top15_PCAD2S_results")



# color blind friendly palette --------------------------------------------

palette <- palette.colors(palette = "Okabe-Ito", names = TRUE)
class_colors <- c("low" = palette[["blue"]],
                  "neutral" = palette[["gray"]],
                  "high" = palette[["vermillion"]])


# load --------------------------------------------------------------------

tx_pair_keys <- read_tsv("03-candidate-transcripts/out/conserved_go_annotated_braker_pairchecked.tsv") |>
  mutate(agic_tx = paste0(agic, ".1"),
         bgc_tx = paste0(braker_gene_id, ".t1")) |>
  select(agic_tx, bgc_tx)

at_tx_list <- readRDS(file = "04-top15/out/at_tx_list.rds")

at_file_ids_pcad1 <- readLines("04-top15/out/a_thaliana_top15_PCAD1_results/fileIDs.txt")
at_file_ids_pcad1_index <- make_tx_index(at_file_ids_pcad1)

at_file_ids_pcad2s <- readLines("04-top15/out/a_thaliana_top15_PCAD2S_results/fileIDs.txt")
at_file_ids_pcad2s_index <- make_tx_index(at_file_ids_pcad2s)


zs_tx_list <- readRDS(file = "04-top15/out/zs_tx_list.rds")

zs_file_ids_pcad1 <- readLines("04-top15/out/z_stapffii_top15_PCAD1_results//fileIDs.txt")
zs_file_ids_pcad1_index <- make_tx_index(zs_file_ids_pcad1)

zs_file_ids_pcad2s <- readLines("04-top15/out/z_stapffii_top15_PCAD2S_results/fileIDs.txt")
zs_file_ids_pcad2s_index <- make_tx_index(zs_file_ids_pcad2s)

# filter 42 genes to the 15 actually used
tx_pair_keys <- tx_pair_keys |>
  filter(agic_tx %in% names(at_file_ids_pcad1_index),
         bgc_tx %in% names(zs_file_ids_pcad1_index))

# -------------------------------------------------------------------------

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

  # MSIC
  msic_at <- bind_rows(
    AT_PCAD1  = read_msic(at_pcad1_dir,  "a_thaliana_top15_PCAD1_scores_table",  at_file_ids_pcad1_index,  agic_tx),
    AT_PCAD2S = read_msic(at_pcad2s_dir, "a_thaliana_top15_PCAD2S_scores_table", at_file_ids_pcad2s_index, agic_tx),
    .id = "run"
  )
  msic_zs <- bind_rows(
    ZS_PCAD1  = read_msic(zs_pcad1_dir,  "z_stapffii_top15_PCAD1_scores_table",     zs_file_ids_pcad1_index,  bgc_tx),
    ZS_PCAD2S = read_msic(zs_pcad2s_dir, "z_stapffii_top15_PLANTCAD2S_scores_table", zs_file_ids_pcad2s_index, bgc_tx),
    .id = "run"
  )

  # sequences + blastn
  at_fa <- write_fasta(getSeq(at_genome_fa, at_tx_list[[agic_tx]]$current_tx), blast_dir, agic_tx)
  zs_fa <- write_fasta(getSeq(zs_genome_fa, zs_tx_list[[bgc_tx]]$current_tx), blast_dir, bgc_tx)
  hsps  <- run_blastn(at_fa, zs_fa, file.path(blast_dir, paste0(pair_name, ".tsv")))

  # plot
  # TODO: write plot_dotplot_marginal() in R/plot-helpers.R
  # p <- plot_dotplot_marginal(hsps, msic_at, msic_zs)
  # ggsave(file.path(plot_dir, paste0(pair_name, ".png")), p,
  #        dpi = 300, width = 10, height = 10)

}


# run ---------------------------------------------------------------------

walk2(tx_pair_keys$agic_tx, tx_pair_keys$bgc_tx, process_pair)
