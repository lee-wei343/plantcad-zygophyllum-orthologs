# shared helpers for loading genomes and extracting flanked transcripts
#
# Used by 03-candidate-transcripts/, 04-top15/ and 05-additional-genes/.
# Expects the calling script to have loaded Rsamtools, rtracklayer, txdbmaker
# and ggplot2.


# load --------------------------------------------------------------------

# GFF is null if no gff given
load_genome <- function(genome_path, gff_path = NULL) {
  if (!file.exists(paste0(genome_path, ".fai"))) {
    indexFa(genome_path)
  }
  genome <- FaFile(genome_path)

  gff <- NULL
  if (!is.null(gff_path)) {
    gff <- import(gff_path)
  }

  list("fasta" = genome, "gff" = gff)
}

# Works for GFFs with transcript names (A. thaliana, T. mongolica)
get_exons_by_tx <- function(gff_path) {
  tx_db <- makeTxDbFromGFF(gff_path, format = "gff3")
  exonsBy(tx_db, by = "tx", use.names = TRUE)
}

# TxDb can't read BRAKER transcript names
split_exons_by_tx <- function(gff) {
  exon_rows <- gff[gff$type == "exon"]
  split(exon_rows, exon_rows$transcript_id)
}


# extract -----------------------------------------------------------------

add_flanks <- function(exons_by_tx, flank_bp, genome_fa) {
  exon_spans_pre <- unlist(range(exons_by_tx))

  contig_lengths <- seqlengths(seqinfo(genome_fa))
  seqlengths(exon_spans_pre) <- contig_lengths[seqlevels(exon_spans_pre)]
  stopifnot(!any(is.na(seqlengths(exon_spans_pre))))

  exon_spans <- exon_spans_pre + flank_bp

  exon_spans_trim <- trim(exon_spans)
  stopifnot(identical(exon_spans, exon_spans_trim))
  exon_spans_trim
}

make_header <- function(exon_spans, exons_by_tx, flank_bp,
                        seq_label = "chrom", extra = NULL) {
  tx_ids <- names(exon_spans)
  chrom <- as.character(seqnames(exon_spans))
  strand <- as.character(strand(exon_spans))
  range <- paste0(
    start(exon_spans), "-",
    end(exon_spans)
  )
  n_exons_by_tx <- lengths(exons_by_tx[tx_ids])

  id <- tx_ids
  if (!is.null(extra)) {
    stopifnot(length(extra) == length(tx_ids))
    id <- paste0(tx_ids, "|", extra)
  }

  paste0(
    id, "|",
    seq_label, ":", chrom, "|",
    strand, "|",
    range, "|",
    "n_exons:", n_exons_by_tx, "|",
    "flank:", flank_bp
  )
}


# plot --------------------------------------------------------------------

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

plot_transcripts <- function(exon_spans, exons_by_tx, headers, out_dir) {
  stopifnot(length(headers) == length(exon_spans))
  tx_list <- list()

  for (i in seq_along(exon_spans)) {
    tx_id <- names(exon_spans)[i]
    current_tx <- exon_spans[tx_id]
    region_start <- start(current_tx)

    current_tx_exons <- exons_by_tx[[tx_id]]
    relative_exon_start <- start(current_tx_exons) - region_start + 1
    relative_exon_end   <- end(current_tx_exons) - region_start + 1

    tx_plot <- plot_transcript_arrows(
      start = 1,
      end = width(current_tx),
      strand = as.character(strand(current_tx)),
      exon_starts = relative_exon_start,
      exon_ends = relative_exon_end,
      title = headers[i]
    )

    tx_list[[tx_id]] <- list(current_tx = current_tx,
                             start = 1,
                             end = width(current_tx),
                             strand = as.character(strand(current_tx)),
                             relative_exon_start = relative_exon_start,
                             relative_exon_end = relative_exon_end,
                             header = headers[i])

    ggsave(
      filename = paste0(tx_id, ".png"),
      plot = tx_plot,
      path = out_dir,
      width = 10,
      height = 2.5,
      dpi = 300,
      create.dir = TRUE
    )
  }

  invisible(tx_list)
}
