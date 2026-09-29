library(Rsamtools)
library(rtracklayer)
library(txdbmaker)

flank_bp <- 1000L
write_fasta <- TRUE
make_plots <- TRUE

plot_transcript <- function(start, end, strand, exon_starts, exon_ends, title) {
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

# Arabidopsis thaliana ----------------------------------------------------

# at = Arabidopsis thaliana
# Arabidopsis Genome Initiative Code
busco <- readr::read_tsv("01-busco/out/conserved_go_annotated.tsv")

at_genome_path <- "data/A_thaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
if (!file.exists(paste0(at_genome_path, ".fai"))) {
  indexFa(at_genome_path)
}
at_genome <- FaFile(at_genome_path)

at_gff_path <- "data/A_thaliana/Arabidopsis_thaliana.TAIR10.63.gff3"
at_gff <- import(at_gff_path)

at_txdb <- makeTxDbFromGFF(file = at_gff_path)
at_exons_by_tx <- exonsBy(at_txdb, by = "tx", use.name = TRUE)

at_exons_by_busco_tx <- at_exons_by_tx[
  names(at_exons_by_tx) %in% paste0(busco$agic, ".1")
  ]

at_exon_spans_pre <- unlist(range(at_exons_by_busco_tx))

at_exon_spans <- at_exon_spans_pre + flank_bp
# stoping if contig or chromosome ends before the flank range
stopifnot(all(at_exon_spans_pre + flank_bp * 2 == width(at_exon_spans)))

# If there is nothing to trim, they remain same legnth
at_exon_spans_trim <- trim(at_exon_spans)
stopifnot(identical(at_exon_spans, at_exon_spans_trim))

at_tx_seq <- getSeq(at_genome, at_exon_spans_trim)
stopifnot(
  length(at_exon_spans_trim) == length(at_tx_seq),
  identical(width(at_exon_spans_trim), width(at_tx_seq))
)

at_tx_id <- names(at_exon_spans_trim)
at_chrom <- as.character(seqnames(at_exon_spans_trim))
at_strand <- as.character(strand(at_exon_spans_trim))
at_range <- paste0(
  start(at_exon_spans_trim), "-",
  end(at_exon_spans_trim)
)
at_n_exons_by_tx <- lengths(at_exons_by_busco_tx)

headers <- paste0(
  at_tx_id, "|",
  "chrom:", at_chrom, "|",
  at_strand, "|",
  at_range, "|",
  "n_exons:", at_n_exons_by_tx, "|",
  "flank:", flank_bp
)

names(at_tx_seq) <- headers

if (write_fasta) {
  writeXStringSet(at_tx_seq, filepath = "out/tx_z_stapffii.fa")
}

