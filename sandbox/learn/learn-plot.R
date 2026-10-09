library(ggplot2)

toy_exons <- data.frame(
  exon_number = 1:4,
  start = c(100, 400, 800, 1200),
  end   = c(250, 520, 1000, 1350)
)
toy_strand <- "+"

plot_transcript <- function(start, end, strand, exon_starts, exon_ends, title) {
  if (strand == "+") {
    direction <- 1
  } else {
    direction <- -1
  }
  
  arrow_spacing <- 200
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
      arrow = arrow(length = unit(5, "mm"))
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

plot_transcript(title = "bob", start = 0, end = 2000, strand = toy_strand, exon_starts = toy_exons$start, exon_ends = toy_exons$end)
