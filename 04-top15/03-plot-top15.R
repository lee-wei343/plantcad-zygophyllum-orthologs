library(zoo)
library(tidyverse)
library(stringr)
library(patchwork)
library(pwalign)
library(Rsamtools)

setwd("~/plantcad-zygophyllum-orthologs/")

source("R/extract-helpers.R")
source("04-top15/R/plot-helpers.R")

# color-blind palette -----------------------------------------------------

# colorblind friendly
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
    plot_layout(heights = c(1, 4, 1, 4), guides = "collect")

  combined_plot

  multi_path <- "04-top15/out/multi_plot2"
  dir.create(multi_path, showWarnings = FALSE)
  ggsave(paste0(multi_path, "/", current_agic_tx, "_with_", current_bgc_tx, ".png"),
         plot = combined_plot,
         dpi = 300,
         height = 7,
         width = 12)
  
  # current_at_tx <- at_tx_list[[current_agic_tx]]$current_tx
  # current_zs_tx <- zs_tx_list[[current_bgc_tx]]$current_tx
  # nuc_sub_matr <- nucleotideSubstitutionMatrix(match = 2,
  #                                              mismatch = -1,
  #                                              baseOnly = TRUE)
  # 
  # pw_align <- pairwiseAlignment(
  #   pattern = getSeq(at_genome_fa, current_at_tx),
  #   subject = getSeq(zs_genome_fa, current_zs_tx),
  #   type = "local",
  #   substitutionMatrix = nuc_sub_matr,
  #   gapOpening = 3,
  #   gapExtension = 1
  # )
  # print(pw_align)
}

  
