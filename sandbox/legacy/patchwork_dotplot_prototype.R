library(ggplot2)
library(patchwork)
library(GenomicRanges)

# -------------------------------------------------------------------------
cblastn <- read_tsv(
  "~/plantcad-zygophyllum-orthologs/04-top15/out/blastn/AT1G71220.1_with_g20167.t1.tsv"
) |>
  filter(length > 40)

current_q_tx_list <- at_tx_list[["AT1G71220.1"]]
rel_qstart <- current_q_tx_list$start
rel_qend <- current_q_tx_list$end
rel_q_exon_start <- current_q_tx_list$relative_exon_start
rel_q_exon_end <- current_q_tx_list$relative_exon_end

current_s_tx_list <- zs_tx_list[["g20167.t1"]]
rel_sstart <- current_s_tx_list$start
rel_send <- current_s_tx_list$end
rel_s_exon_start <- current_s_tx_list$relative_exon_start
rel_s_exon_end <- current_s_tx_list$relative_exon_end

dotplot_bn <- function(
    blastn_tsv, rel_qstart, rel_qend, rel_sstart, rel_send, rel_q_exon_start,
    rel_q_exon_end, rel_s_exon_start, rel_s_exon_end
    ) {
  palette <- palette.colors(palette = "Okabe-Ito", names = TRUE)
  
  main <- ggplot(blastn_tsv) +
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

plot <- dotplot_bn(
  blastn_tsv = cblastn,
  rel_qstart = rel_qstart,
  rel_qend = rel_qend,
  rel_sstart = rel_sstart,
  rel_send = rel_send,
  rel_q_exon_start = rel_q_exon_start,
  rel_q_exon_end = rel_q_exon_end,
  rel_s_exon_start = rel_s_exon_start,
  rel_s_exon_end = rel_s_exon_end
  )
plot
