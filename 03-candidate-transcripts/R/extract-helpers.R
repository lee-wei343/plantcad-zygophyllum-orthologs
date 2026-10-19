# functions for script 01-extract-transcripts.R

load_genome <- function(genome_path, gff_path) {
  if (!file.exists(paste0(genome_path, ".fai"))) {
    indexFa(genome_path)
  }
  genome <- FaFile(genome_path)
  
  gff <- import(gff_path)
  
  list("fasta" = genome, "gff" = gff)
}

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

add_flanks <- function(exons_by_tx, flank_bp, genome_fa) {
  exon_spans_pre <- unlist(range(exons_by_tx))
  
  contig_lengths <- seqlengths(seqinfo(genome_fa))
  seqlengths(exon_spans_pre) <- contig_lengths[seqlevels(exon_spans_pre)]
  stopifnot(!any(is.na(seqlengths(exon_spans_pre))))
  
  exon_spans <- exon_spans_pre + flank_bp
  
  # If there is nothing to trim, they remain same legnth
  exon_spans_trim <- trim(exon_spans)
  stopifnot(identical(exon_spans, exon_spans_trim))
  exon_spans_trim
}

make_header <- function(exon_spans, exons_listed, flank_bp) {
  tx_ids <- names(exon_spans)
  chrom <- as.character(seqnames(exon_spans))
  strand <- as.character(strand(exon_spans))
  range <- paste0(
    start(exon_spans), "-",
    end(exon_spans)
  )
  n_exons_by_tx <- lengths(exons_listed)
  
  headers <- paste0(
    tx_ids, "|",
    "chrom:", chrom, "|",
    strand, "|",
    range, "|",
    "n_exons:", n_exons_by_tx, "|",
    "flank:", flank_bp
  )
}