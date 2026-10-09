library(zoo)
library(tidyverse)
library(stringr)
library(patchwork)
library(pwalign)
library(Rsamtools)

setwd("~/plantcad-zygophyllum-orthologs/")


# color-blind palette -----------------------------------------------------

# colorblind friendly
palette <- palette.colors(palette = "Okabe-Ito", names = TRUE)
class_colors <- c("low" = palette[["blue"]],
                  "neutral" = palette[["gray"]],
                  "high" = palette[["vermillion"]])


# functions ---------------------------------------------------------------

add_pos <- function(df) {
  df$pos <- seq_len(nrow(df))
  return(df)
}
add_class <- function(df) {
  df <- df |>
    mutate(class = case_when(
      is.na(MSIC) ~ NA_character_,
      MSIC > 0.5 ~ "high",
      MSIC < -0.5 ~ "low",
      .default = "neutral"
    )) |>
    mutate(class = factor(class, levels = c("high", "neutral", "low")))
  return(df)
}
add_rolling_mean <- function(df) {
  rolling <- rollapply(data = df$MSIC, # change MSIC to IC
                       width = list(-30:29),
                       FUN = mean,
                       partial = TRUE)
  df$rolling_mean <- (rolling) 
  return(df)
}

make_tx_index <- function(id_file) {
  indexed <- list()
  for (i in seq_along(id_file)) {
    number <- str_split_1(id_file[i], ": ")[1]
    header <- str_split_1(id_file[i], ": ")[2]
    id <- str_split_1(header, "\\|")[1]
    
    indexed[id] <- number
  }
  return(indexed)
}

plot_transcript <- function(start, end, strand, exon_starts, exon_ends, title) {
  # flip minus strand, so both line up 5' to '3 
  if (strand == "-") {
    flipped_starts <- start + end - exon_ends
    flipped_ends <- start + end - exon_starts
    exon_starts <- flipped_starts
    exon_ends <- flipped_ends
  }

  ggplot() +
    geom_segment(
      aes(x = start, xend = end, y = 1, yend = 1)
    ) +
    geom_rect(
      aes(xmin = exon_starts, xmax = exon_ends, ymin = 0.8, ymax = 1.2),
      fill = palette[["gray"]]
    ) +
    scale_y_continuous(breaks = NULL) +
    labs(y = NULL,
         x = NULL,
         title = title) +
    theme_minimal()
}

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

for (i in seq_len(nrow(tx_pair_keys))) {
  
  # BREAK HERE
  if (i != 1) break
  
  
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
  
  at_agic_tx_pcad1 <- add_rolling_mean(add_class(add_pos(at_agic_tx_pcad1)))
  at_agic_tx_pcad2s <- add_rolling_mean(add_class(add_pos(at_agic_tx_pcad2s)))

  zs_bgc_tx_pcad1 <- add_rolling_mean(add_class(add_pos(zs_bgc_tx_pcad1)))
  zs_bgc_tx_pcad2s <- add_rolling_mean(add_class(add_pos(zs_bgc_tx_pcad2s)))

  msic_table_at <- list(
    AT_PCAD1 = at_agic_tx_pcad1,
    AT_PCAD2S = at_agic_tx_pcad2s
  )

  msic_table_zs <- list(
    ZS_PCAD1 = zs_bgc_tx_pcad1,
    ZS_PCAD2S = zs_bgc_tx_pcad2s
    )

  msic_long_at <- bind_rows(msic_table_at, .id = "run") |>
    mutate(run = factor(run, levels = names(msic_table_at)))

  msic_long_zs <- bind_rows(msic_table_zs, .id = "run") |>
    mutate(run = factor(run, levels = names(msic_table_zs)))

  x_range  <- c(1, max(msic_long_at$pos))
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
  tx_plot_at <- plot_transcript(
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
  tx_plot_zs <- plot_transcript(
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

  multi_path <- "04-top15/out/multi_plot"
  dir.create(multi_path, showWarnings = FALSE)
  ggsave(paste0(multi_path, "/", current_agic_tx, "_with_", current_bgc_tx, ".png"),
         plot = combined_plot,
         dpi = 300,
         height = 7,
         width = 12)
  
  current_at_tx <- at_tx_list[[current_agic_tx]]$current_tx
  current_zs_tx <- zs_tx_list[[current_bgc_tx]]$current_tx
  nuc_sub_matr <- nucleotideSubstitutionMatrix(match = 2,
                                               mismatch = -1,
                                               baseOnly = TRUE)

  pw_align <- pairwiseAlignment(
    pattern = getSeq(at_genome_fa, current_at_tx),
    subject = getSeq(zs_genome_fa, current_zs_tx),
    type = "local",
    substitutionMatrix = nuc_sub_matr,
    gapOpening = 3,
    gapExtension = 1
  )
  print(pw_align)
  
  break
  
# dotplot AI generated! ---------------------------------------------------

  
  current_at_tx <- at_tx_list[[current_agic_tx]]$current_tx
  current_zs_tx <- zs_tx_list[[current_bgc_tx]]$current_tx
  seq_at <- as.character(getSeq(at_genome_fa, current_at_tx))
  seq_zs <- as.character(getSeq(zs_genome_fa, current_zs_tx))
  
  window    <- 30   # window length along the diagonal
  threshold <- 21   # minimum matching bases per window (70%)
  
  # Convert bases to numbers: A=1, C=2, G=3, T=4
  # Anything else (e.g. N) gets 5 in at and 6 in zs, so N never matches N
  codes_at <- match(strsplit(seq_at, "")[[1]], c("A", "C", "G", "T"), nomatch = 5)
  codes_zs <- match(strsplit(seq_zs, "")[[1]], c("A", "C", "G", "T"), nomatch = 6)
  
  n_at <- length(codes_at)
  n_zs <- length(codes_zs)
  
  # Every diagonal, identified by its offset zs_pos - at_pos
  diagonals <- seq(-(n_at - 1), n_zs - 1)
  segment_list <- vector("list", length(diagonals))
  
  for (j in seq_along(diagonals)) {
    offset <- diagonals[j]
    
    # at positions that lie on this diagonal inside the grid
    first_at <- max(1, 1 - offset)
    last_at  <- min(n_at, n_zs - offset)
    diagonal_length <- last_at - first_at + 1
    if (diagonal_length < window) next
    
    at_positions <- first_at:last_at
    is_match <- codes_at[at_positions] == codes_zs[at_positions + offset]
    
    # Matches per window, via a running total
    running_total  <- cumsum(c(0, is_match))
    window_matches <- running_total[(window + 1):(diagonal_length + 1)] -
      running_total[1:(diagonal_length - window + 1)]
    
    passes <- window_matches >= threshold
    if (!any(passes)) next
    
    # Merge neighbouring passing windows into one segment
    runs       <- rle(passes)
    run_ends   <- cumsum(runs$lengths)
    run_starts <- run_ends - runs$lengths + 1
    keep       <- runs$values
    
    segment_start <- at_positions[run_starts[keep]]
    segment_end   <- at_positions[run_ends[keep]] + window - 1
    
    segment_list[[j]] <- data.frame(
      at_start = segment_start,
      at_end   = segment_end,
      zs_start = segment_start + offset,
      zs_end   = segment_end + offset
    )
  }
  
  segments <- bind_rows(segment_list)
  nrow(segments)
  
  
  min_length <- 20   # drop segments shorter than this
  
  long_segments <- segments %>%
    mutate(segment_length = at_end - at_start + 1) %>%
    filter(segment_length >= min_length)
  
  nrow(long_segments)
  
  dot_plot <- ggplot(long_segments) +
    geom_segment(aes(x = at_start, y = zs_start, xend = at_end, yend = zs_end),
                 colour = "darkcyan", linewidth = 0.6) +
    coord_fixed(xlim = c(1, n_at), ylim = c(1, n_zs)) +
    labs(x = "A. thaliana position (bp)",
         y = "Z. stapffii position (bp)",
         title = ) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "grey40", fill = NA))
  
  dot_plot_path <- "04-top15/out/dot_plot"
  dir.create(dot_plot_path, showWarnings = FALSE)
  ggsave(paste0(dot_plot_path, "/", current_agic_tx, "_with_", current_bgc_tx, ".png"),
         plot = dot_plot,
         dpi = 300,
         height = 7,
         width = 8)
}

  
