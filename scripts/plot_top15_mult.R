library(zoo)
library(tidyverse)
library(stringr)
library(patchwork)

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

tx_pair_keys <- read_tsv("out/conserved_go_annotated_braker_pairchecked.tsv") |>
  mutate(agic_tx = paste0(agic, ".1"),
         bgc_tx = paste0(braker_gene_id, ".t1")) |>
  select(agic_tx, bgc_tx)

at_tx_list <- readRDS(file = "out/top15/at_tx_list.rds")

at_file_ids_pcad1 <- readLines("out/top15/a_thaliana_top15_PCAD1_results/fileIDs.txt")
at_file_ids_pcad1_index <- make_tx_index(at_file_ids_pcad1)

at_file_ids_pcad2s <- readLines("out/top15/a_thaliana_top15_PCAD2S_results/fileIDs.txt")
at_file_ids_pcad2s_index <- make_tx_index(at_file_ids_pcad2s)


zs_tx_list <- readRDS(file = "out/top15/zs_tx_list.rds")

zs_file_ids_pcad1 <- readLines("out/top15/z_stapffii_top15_PCAD1_results//fileIDs.txt")
zs_file_ids_pcad1_index <- make_tx_index(zs_file_ids_pcad1)

zs_file_ids_pcad2s <- readLines("out/top15/z_stapffii_top15_PCAD2S_results/fileIDs.txt")
zs_file_ids_pcad2s_index <- make_tx_index(zs_file_ids_pcad2s)

# filter 42 genes to the 15 actually used
tx_pair_keys <- tx_pair_keys |>
  filter(agic_tx %in% names(at_file_ids_pcad1_index),
         bgc_tx %in% names(zs_file_ids_pcad1_index))

# -------------------------------------------------------------------------

for (i in seq_len(nrow(tx_pair_keys))) {
  current_agic_tx <- tx_pair_keys$agic_tx[i]
  at_agic_tx_pcad1 <- read_csv(
    paste0(
      "out/top15/a_thaliana_top15_PCAD1_results/a_thaliana_top15_PCAD1_scores_table",
      at_file_ids_pcad1_index[current_agic_tx],
      ".csv"
    )
  )
  at_agic_tx_pcad2s <- read_csv(
    paste0(
      "out/top15/a_thaliana_top15_PCAD2S_results/a_thaliana_top15_PCAD2S_scores_table",
      at_file_ids_pcad2s_index[current_agic_tx],
      ".csv"
    )
  )
  
  current_bgc_tx <- tx_pair_keys$bgc_tx[i]
  zs_bgc_tx_pcad1 <- read_csv(
    paste0(
      "out/top15/z_stapffii_top15_PCAD1_results/z_stapffii_top15_PCAD1_scores_table",
      zs_file_ids_pcad1_index[current_bgc_tx],
      ".csv"
    )
  )
  zs_bgc_tx_pcad2s <- read_csv(
    paste0(
      "out/top15/z_stapffii_top15_PCAD2S_results/z_stapffii_top15_PLANTCAD2S_scores_table",
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
  
  multi_path <- "out/top15/multi_plot"
  dir.create(multi_path)
  ggsave(paste0(multi_path, "/", current_agic_tx, "_with_", current_bgc_tx, ".png"),
         plot = combined_plot,
         dpi = 300,
         height = 7,
         width = 12)
}

# -------------------------------------------------------------------------

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

# flip minus strand, so both line up 5' to '3 
if (strand == "-") {
  flipped_starts <- start + end - exon_ends
  flipped_ends <- start + end - exon_starts
  exon_starts <- flipped_starts
  exon_ends <- flipped_ends
}

ggplot() +
  geom_segment(
    aes(x = start, xend = end, y = 1.5, yend = 1)
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


ggplot() +
  geom_segment(
    aes(
      x = current_zs_tx_list$start,
      xend = current_zs_tx_list$end,
      y = 1.5,
      yend = 1.5
    )
  ) +
  geom_segment(
    aes(
      x = current_at_tx_list$start,
      xend = current_at_tx_list$end,
      y = -1.5,
      yend = -1.5
    )
  ) +
  geom_rect(
    aes(xmin = current_at_tx_list$relative_exon_start,
        xmax = current_at_tx_list$relative_exon_end,
        ymin = -1,
        ymax = -2),
    fill = palette[["gray"]]
  ) +
  geom_rect(
    aes(
      xmin = current_zs_tx_list$relative_exon_start,
      xmax = current_zs_tx_list$relative_exon_end,
      ymin = 1,
      ymax = 2
    ),
    fill = palette[["gray"]]
  ) +
  scale_y_continuous(breaks = NULL) +
  labs(y = NULL,
       x = NULL,
       title = NULL) +
  theme_minimal()
