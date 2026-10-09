library(ggplot2)
library(patchwork)
library(GenomicRanges)

# -------------------------------------------------------------------------
cblastn <- read_tsv(
  "~/plantcad-zygophyllum-orthologs/04-top15/out/blastn/AT1G71220.1_with_g20167.t1.tsv"
) |>
  filter(length > 50)

current_q_tx_list <- at_tx_list[["AT1G71220.1"]]
qstart <- current_q_tx_list$start
qend <- current_q_tx_list$end

current_s_tx_list <- zs_tx_list[["g20167.t1"]]
sstart <- current_s_tx_list$start
send <- current_s_tx_list$end


main <- ggplot(mtcars, aes(wt, mpg, colour = factor(cyl))) +
  annotate("rect",
    xmin = 1, xmax = 2,
    ymin = -Inf, ymax = Inf,
    fill = "blue", alpha = 0.2
  ) +
  annotate("rect",
    xmin = -Inf, xmax = Inf,
    ymin = 20, ymax = 25,
    fill = "red", alpha = 0.2
  ) +
  geom_point() +
  coord_cartesian(xlim = c(sstart, send), ylim = c(qstart, qend)) +
  theme(legend.position = "none") +
  labs(
    y = "A. thaliana (bp)",
    x = "Z. stapffii (bp)"
  )

main

top <- ggplot(mtcars, aes(wt, fill = factor(cyl))) +
  geom_density() +
  xlim(sstart, send) +
  theme_void() +
  theme(legend.position = "none")

right <- ggplot(mtcars, aes(y = mpg, fill = factor(cyl))) +
  geom_histogram(bins = 10) +
  ylim(qstart, qend) +
  theme_void() +
  theme(legend.position = "none")

top + plot_spacer() + main + right +
  plot_layout(
    ncol = 2,
    widths = c(4, 1),
    heights = c(1, 4),
    guides = "collect"
  ) &
  theme(plot.margin = margin(0, 0, 0, 0))
