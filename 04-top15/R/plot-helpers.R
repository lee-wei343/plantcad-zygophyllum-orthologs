# helper script
#
# Expects the calling script to have loaded: tidyverse, zoo, Biostrings
# (via Rsamtools), and to define `palette` (Okabe-Ito) before plotting.

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
  rolling <- rollapply(data = df$MSIC, # TODO: change MSIC to IC
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

plot_transcript_simple <- function(start, end, strand, exon_starts, exon_ends,
                                   title
                                   ) {
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

plot_transcript_arrows <- function(start, end, strand, exon_starts, exon_ends,
                                   title
                                   ) {
  if (strand == "+") {
    direction <- 1
  } else {
    direction <- -1
  }
  
  arrow_spacing <- 1000
  arrow_length <- 20
  arrow_starts <- seq(
    from = start,
    to = end,
    by = arrow_spacing
  )
  arrows <- data.frame(
    x = arrow_starts,
    xend = arrow_starts + (direction * arrow_length)
  )
  
  ggplot() +
    geom_segment(
      aes(x = start, xend = end, y = 1, yend = 1)
    ) +
    geom_segment(
      data = arrows,
      aes(x = x, xend = xend, y = 1, yend = 1),
      arrow = arrow(length = unit(2, "mm"))
    ) +
    geom_rect(
      aes(xmin = exon_starts, xmax = exon_ends, ymin = 0.8, ymax = 1.2),
      fill = "blue"
    ) +
    scale_y_continuous(breaks = NULL) +
    labs(x = "position (bp)",
         y = NULL,
         title = title) +
    theme_minimal()
}


# MSIC --------------------------------------------------------------------

read_msic <- function(results_dir, file_prefix, index, tx_id) {
  path <- file.path(results_dir, paste0(file_prefix, index[[tx_id]], ".csv"))
  read_csv(path, show_col_types = FALSE) |>
    add_pos() |>
    add_class() |>
    add_rolling_mean()
}


# blastn ------------------------------------------------------------------

blast_cols <- c("qstart", "qend", "sstart", "send",
                "pident", "length", "evalue", "bitscore")

write_fasta <- function(seq, dir, name) {
  path <- file.path(dir, paste0(name, ".fa"))
  writeXStringSet(seq, path)
  path
}

run_blastn <- function(query_fa, subject_fa, out_tsv,
                       env = "plantcad-zygophyllum-orthologs") {
  outfmt <- paste("6", paste(blast_cols, collapse = " "))
  
  system2(
    "conda",
    c(
      "run", "-n", env, "blastn",
      "-query", query_fa, "-subject", subject_fa, "-task", "blastn",
      "-outfmt", shQuote(outfmt),
      "-out", out_tsv
      )
    )
  hsps <- read_tsv(out_tsv, col_names = blast_cols, show_col_types = FALSE)
  write_tsv(hsps, out_tsv)
  hsps
}


# improved dotplot --------------------------------------------------------

dotplot_bn <- function(
    blastn_tsv, rel_qstart, rel_qend, rel_sstart, rel_send, rel_q_exon_start,
    rel_q_exon_end, rel_s_exon_start, rel_s_exon_end
) {
  palette <- palette.colors(palette = "Okabe-Ito", names = TRUE)
  
  main <- ggplot(blastn_tsv) +
    # TODO: correct the exon starts and ends for the minus strands
    annotate("rect",
             xmin = -Inf, xmax = Inf,
             ymin = rel_q_exon_start, ymax = rel_q_exon_end,
             fill = palette[["vermillion"]], alpha = 0.3
    ) +
    annotate("rect",
             xmin = rel_s_exon_start, xmax = rel_s_exon_end,
             ymin = -Inf, ymax = Inf,
             fill = palette[["blue"]], alpha = 0.3
    ) +
    geom_segment(aes(y = qstart, yend = qend, x = sstart, xend = send)) +
    coord_cartesian(xlim = c(rel_sstart, rel_send), ylim = c(rel_qstart, rel_qend)) +
    theme(legend.position = "none") +
    labs(
      y = "A. thaliana (bp)",
      x = "Z. stapffii (bp)"
    )
  
  top <- ggplot(blastn_tsv) +
    geom_segment(aes(x = sstart, xend = send, y = pident, yend = pident)) +
    xlim(rel_sstart, rel_send) +
    guides(x = "none") +
    labs(x = NULL, y = "% identity") +
    theme(
      legend.position = "none",
    ) +
    theme_minimal()
  
  right <- ggplot(blastn_tsv) +
    geom_segment(aes(y = qstart, yend = qend, x = pident, xend = pident)) +
    guides(y = "none") +
    theme(
      legend.position = "none",
    ) +
    labs(y = NULL, x = "% identity") +
    theme_minimal()
  
  dotplot <- top + plot_spacer() + main + right +
    plot_layout(
      ncol = 2,
      widths = c(6, 1),
      heights = c(1, 6),
      guides = "collect"
    ) &
    theme(plot.margin = margin(0, 0, 0, 0))
} 


