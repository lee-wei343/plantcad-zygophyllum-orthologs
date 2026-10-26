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

out_dir <- "05-additional-genes"
blastn_dir <- file.path(out_dir, "blastn")
blastp_dir <- file.path(out_dir, "blastp")
plot_dir <- file.path(out_dir, "dot_plot")
dir.create(blastn_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(blastp_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

at_genome_path <- "data/A_thaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
at_gff_path <- "data/A_thaliana/Arabidopsis_thaliana.TAIR10.63.gff3"

zs_genome_path <- "data/Z_stapfii/Z_stapfii_genome.fa"
zs_gff_path <- "data/Z_stapfii/braker.gff3"

flank_bp <- 1000L


# color blind friendly palette --------------------------------------------

palette <- palette.colors(palette = "Okabe-Ito", names = TRUE)
class_colors <- c(
  "low" = palette[["blue"]],
  "neutral" = palette[["gray"]],
  "high" = palette[["vermillion"]]
)


# load --------------------------------------------------------------------

tx_pair_keys <- read_tsv(
  "05-additional-genes/tps_orthos_in_busco.tsv"
  ) |>
  mutate(
    agic_tx = paste0(Atha, ".1"),
    bgc_tx = paste0(Zsta, ".t1")
  ) |>
  dplyr::select(agic_tx, bgc_tx)


at_genome <- load_genome(
  genome_path = at_genome_path,
  gff_path = at_gff_path
)

zs_genome <- load_genome(
  genome_path = zs_genome_path,
  gff_path = zs_gff_path
)


# make x_tx_by_cds for blastp ---------------------------------------------

at_tx_db <- makeTxDbFromGFF(
  at_gff_path,
  format = "gff3",
  )

# TxDb of the BRAKER gff: transcript names come from Name, not transcript_id
zs_gff <- zs_genome$gff
zs_gff$transcript_id <- NULL
zs_gff$Name <- ifelse(zs_gff$type == "mRNA", zs_gff$ID, NA_character_)
zs_tx_db <- makeTxDbFromGRanges(zs_gff)

at_cds_by_tx <- cdsBy(at_tx_db, by = "tx", use.names = TRUE)
zs_cds_by_tx <- cdsBy(zs_tx_db, by = "tx", use.names = TRUE)


# tx lists ----------------------------------------------------------------

at_exons_by_tx <- exonsBy(at_tx_db, by = "tx", use.names = TRUE)
at_exons_by_tx <- at_exons_by_tx[names(at_exons_by_tx) %in% tx_pair_keys$agic_tx]
at_exon_spans <- add_flanks(at_exons_by_tx, flank_bp, at_genome$fasta)

at_tx_list <- make_tx_list(
  exon_spans = at_exon_spans,
  exons_by_tx = at_exons_by_tx,
  headers = make_header(at_exon_spans, at_exons_by_tx, flank_bp)
)

zs_exons_by_tx <- split_exons_by_tx(zs_genome$gff)
zs_exons_by_tx <- zs_exons_by_tx[names(zs_exons_by_tx) %in% tx_pair_keys$bgc_tx]
zs_exon_spans <- add_flanks(zs_exons_by_tx, flank_bp, zs_genome$fasta)

zs_tx_list <- make_tx_list(
  exon_spans = zs_exon_spans,
  exons_by_tx = zs_exons_by_tx,
  headers = make_header(zs_exon_spans, zs_exons_by_tx, flank_bp, seq_label = "ctg")
)

stopifnot(
  all(tx_pair_keys$agic_tx %in% names(at_tx_list)),
  all(tx_pair_keys$bgc_tx %in% names(zs_tx_list))
)


# load --------------------------------------------------------------------

tx_pair_keys <- read_tsv("03-candidate-transcripts/out/conserved_go_annotated_braker_pairchecked.tsv") |>
  mutate(agic_tx = paste0(agic, ".1"),
         bgc_tx = paste0(braker_gene_id, ".t1")) |>
  select(agic_tx, bgc_tx)

at_file_ids_pcad1 <- readLines("04-top15/out/a_thaliana_top15_PCAD1_results/fileIDs.txt")
at_file_ids_pcad1_index <- make_tx_index(at_file_ids_pcad1)

at_file_ids_pcad2s <- readLines("04-top15/out/a_thaliana_top15_PCAD2S_results/fileIDs.txt")
at_file_ids_pcad2s_index <- make_tx_index(at_file_ids_pcad2s)


zs_file_ids_pcad1 <- readLines("04-top15/out/z_stapffii_top15_PCAD1_results//fileIDs.txt")
zs_file_ids_pcad1_index <- make_tx_index(zs_file_ids_pcad1)

zs_file_ids_pcad2s <- readLines("04-top15/out/z_stapffii_top15_PCAD2S_results/fileIDs.txt")
zs_file_ids_pcad2s_index <- make_tx_index(zs_file_ids_pcad2s)

# filter 42 genes to the 15 actually used
tx_pair_keys <- tx_pair_keys |>
  filter(agic_tx %in% names(at_file_ids_pcad1_index),
         bgc_tx %in% names(zs_file_ids_pcad1_index))

# -------------------------------------------------------------------------

at_genome_fa <- load_genome(
  genome_path = "data/A_thaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
)$fasta

zs_genome_fa <- load_genome(
  genome_path = "data/Z_stapfii/Z_stapfii_genome.fa"
)$fasta

for (i in seq_len(nrow(tx_pair_keys))) {
  
  # thaliana
  current_agic_tx <- tx_pair_keys$agic_tx[i]
  
  at_agic_tx_pcad1 <- read_csv(
    paste0(
      "04-top15/out/a_thaliana_top15_PCAD1_results/a_thaliana_top15_PCAD1_scores_table",
      at_file_ids_pcad1_index[current_agic_tx],
      ".csv"
    )
  )
  
  at_agic_tx_pcad2s <- read_csv(
    paste0(
      "04-top15/out/a_thaliana_top15_PCAD2S_results/a_thaliana_top15_PCAD2S_scores_table",
      at_file_ids_pcad2s_index[current_agic_tx],
      ".csv"
    )
  )
  
  at_agic_tx_pcad2l <- read_csv(
    "04-top15/out/a_thaliana_top15_PCAD2L_results/msic.csv"
  ) |>
    mutate(agic = str_extract(header, "^[^|]+")) |>
    filter(agic == current_agic_tx) |>
    dplyr::rename(MSIC = msic)
  
  at_agic_tx_bot1s <- read_csv(
    "04-top15/out/a_thaliana_top15_Botanic1S_results/msic.csv"
    ) |>
    mutate(agic = str_extract(header, "^[^|]+")) |>
    filter(agic == current_agic_tx) |>
    dplyr::rename(MSIC = msic)
  
  
  # Z_stapffii
  current_bgc_tx <- tx_pair_keys$bgc_tx[i]
  
  zs_bgc_tx_pcad1 <- read_csv(
    paste0(
      "04-top15/out/z_stapffii_top15_PCAD1_results/z_stapffii_top15_PCAD1_scores_table",
      zs_file_ids_pcad1_index[current_bgc_tx],
      ".csv"
    )
  )
  
  zs_bgc_tx_pcad2s <- read_csv(
    paste0(
      "04-top15/out/z_stapffii_top15_PCAD2S_results/z_stapffii_top15_PLANTCAD2S_scores_table",
      zs_file_ids_pcad2s_index[current_bgc_tx],
      ".csv"
    )
  )
  
  zs_bgc_tx_pcad2l <- read_csv(
    "04-top15/out/z_stapffii_top15_PCAD2L_result/msic.csv"
  ) |>
    mutate(bgc = str_extract(header, "^[^|]+")) |>
    filter(bgc == current_bgc_tx) |>
    dplyr::rename(MSIC = msic)
  
  zs_bgc_tx_bot1s <- read_csv(
    "04-top15/out/z_stapffii_top15_Botanic1S_result/msic.csv"
  ) |>
    mutate(bgc = str_extract(header, "^[^|]+")) |>
    filter(bgc == current_bgc_tx) |>
    dplyr::rename(MSIC = msic)
  
  
  at_agic_tx_pcad1 <- add_rolling_mean(add_class(add_pos(at_agic_tx_pcad1)))
  at_agic_tx_pcad2s <- add_rolling_mean(add_class(add_pos(at_agic_tx_pcad2s)))
  at_agic_tx_pcad2l <- add_rolling_mean(add_class(add_pos(at_agic_tx_pcad2l)))
  at_agic_tx_bot1s <- add_rolling_mean(add_class(add_pos(at_agic_tx_bot1s)))
  
  
  zs_bgc_tx_pcad1 <- add_rolling_mean(add_class(add_pos(zs_bgc_tx_pcad1)))
  zs_bgc_tx_pcad2s <- add_rolling_mean(add_class(add_pos(zs_bgc_tx_pcad2s)))
  zs_bgc_tx_pcad2l <- add_rolling_mean(add_class(add_pos(zs_bgc_tx_pcad2l)))
  zs_bgc_tx_bot1s <- add_rolling_mean(add_class(add_pos(zs_bgc_tx_bot1s)))

  msic_table_at <- list(
    PCAD1 = at_agic_tx_pcad1,
    PCAD2S = at_agic_tx_pcad2s,
    PCAD2L = at_agic_tx_pcad2l,
    BOT1S = at_agic_tx_bot1s
  )

  msic_table_zs <- list(
    PCAD1 = zs_bgc_tx_pcad1,
    PCAD2S = zs_bgc_tx_pcad2s,
    PCAD2L = zs_bgc_tx_pcad2l,
    BOT1S = zs_bgc_tx_bot1s
    )
  

  msic_long_at <- bind_rows(msic_table_at, .id = "run") |>
    mutate(run = factor(run, levels = names(msic_table_at)))

  msic_long_zs <- bind_rows(msic_table_zs, .id = "run") |>
    mutate(run = factor(run, levels = names(msic_table_zs)))

  x_range  <- c(1, max(msic_long_at$pos, msic_long_zs$pos))
  x_scale  <- scale_x_continuous(expand = expansion(mult = 0.01))

  msic_plot_at <- ggplot(msic_long_at, aes(x = pos, y = MSIC)) +
    geom_col(aes(fill = class), width = 1L) +
    geom_hline(yintercept = 0) +
    geom_line(aes(y = rolling_mean), linewidth = 0.3) +
    scale_fill_manual(values = class_colors) +
    facet_grid(run ~ .) +
    x_scale +
    coord_cartesian(xlim = x_range, ylim = c(-1, 1)) +
    labs(y = "MSIC", x = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_blank())


  msic_plot_zs <- ggplot(msic_long_zs, aes(x = pos, y = MSIC)) +
    geom_col(aes(fill = class), width = 1L) +
    geom_hline(yintercept = 0) +
    geom_line(aes(y = rolling_mean), linewidth = 0.3) +
    scale_fill_manual(values = class_colors) +
    facet_grid(run ~ .) +
    x_scale +
    coord_cartesian(xlim = x_range, ylim = c(-1, 1)) +
    labs(y = "MSIC", x = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_blank())

  current_at_tx_list <- at_tx_list[[current_agic_tx]]
  tx_plot_at <- plot_transcript_simple(
    start = current_at_tx_list$start,
    end = current_at_tx_list$end,
    strand = current_at_tx_list$strand,
    exon_starts = current_at_tx_list$relative_exon_start,
    exon_ends = current_at_tx_list$relative_exon_end,
    title = current_at_tx_list$header
  ) +
    x_scale +
    coord_cartesian(xlim = x_range)

  current_zs_tx_list <- zs_tx_list[[current_bgc_tx]]
  tx_plot_zs <- plot_transcript_simple(
    start = current_zs_tx_list$start,
    end = current_zs_tx_list$end,
    strand = current_zs_tx_list$strand,
    exon_starts = current_zs_tx_list$relative_exon_start,
    exon_ends = current_zs_tx_list$relative_exon_end,
    title = current_zs_tx_list$header
  ) +
    x_scale +
    coord_cartesian(xlim = x_range)

  combined_plot <- tx_plot_at / msic_plot_at / tx_plot_zs / msic_plot_zs +
    plot_layout(heights = c(1, 8, 1, 8), guides = "collect")

  combined_plot

  multi_path <- "04-top15/out/multi_plot2"
  dir.create(multi_path, showWarnings = FALSE)
  ggsave(paste0(multi_path, "/", current_agic_tx, "_with_", current_bgc_tx, ".png"),
         plot = combined_plot,
         dpi = 300,
         height = 7,
         width = 12)
}