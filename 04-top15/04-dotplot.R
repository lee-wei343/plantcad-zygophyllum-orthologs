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
  
  blastn_tsv_path <- file.path(blast_dir, paste0(pair_name, ".tsv"))
  blastn_tsv <- run_blastn(
    at_fa,
    zs_fa,
    blastn_tsv_path
  )
  
  current_q_tx_list <- at_tx_list[[agic_tx]]
  rel_qstart <- current_q_tx_list$start
  rel_qend <- current_q_tx_list$end
  q_exons <- flip_rel_exons(strand = current_q_tx_list$strand,
                 start = rel_qstart,
                 end = rel_qend,
                 exon_starts = current_q_tx_list$relative_exon_start,
                 exon_ends = current_q_tx_list$relative_exon_end)
  rel_q_exon_start <- q_exons$exon_starts
  rel_q_exon_end <- q_exons$exon_ends
  
  
  current_s_tx_list <- zs_tx_list[[bgc_tx]]
  rel_sstart <- current_s_tx_list$start
  rel_send <- current_s_tx_list$end
  s_exons <- flip_rel_exons(strand = current_s_tx_list$strand,
                 start = rel_sstart,
                 end = rel_send,
                 exon_starts = current_s_tx_list$relative_exon_start,
                 exon_ends = current_s_tx_list$relative_exon_end)
  rel_s_exon_start <- s_exons$exon_starts
  rel_s_exon_end <- s_exons$exon_ends

  dotplot <- dotplot_bn(
    blastn_tsv = blastn_tsv,
    rel_qstart = rel_qstart,
    rel_qend = rel_qend,
    rel_sstart = rel_sstart,
    rel_send = rel_send,
    rel_q_exon_start = rel_q_exon_start,
    rel_q_exon_end = rel_q_exon_end,
    rel_s_exon_start = rel_s_exon_start,
    rel_s_exon_end = rel_s_exon_end
  )
  ggsave(file.path(plot_dir, paste0(pair_name, ".png")), dotplot,
         dpi = 300, width = 10, height = 10)
  
  
}


# run ---------------------------------------------------------------------

walk2(tx_pair_keys$agic_tx, tx_pair_keys$bgc_tx, process_pair)
