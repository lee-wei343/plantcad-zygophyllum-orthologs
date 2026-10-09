library(Rsamtools)
library(rtracklayer)

setwd("~/plantcad-zygophyllum-orthologs/")


# loading -----------------------------------------------------------------

genome_path <- "data/Z_stapfii/Z_stapfii_genome.fa"
if (!file.exists(paste0(genome_path, ".fai"))) {
  indexFa(genome_path)
}
genome_fa <- FaFile(genome_path)

gff <- import("data/Z_stapfii/braker.gff3")
genes <- gff[gff$type == "gene"]

conserved <- data.table::fread("out/busco/conserved.tsv")

# quick look --------------------------------------------------------------

contig_lengths <- seqlengths(seqinfo(genome_fa))

# n contigs
length(contig_lengths)

which.max(contig_lengths)
max(contig_lengths)

table(gff$type)


# BUSCO cords to GRanges --------------------------------------------------

# GRanges always start as lower, BUSCO gave higher as start when minus strand
hit_start <- pmin(conserved$gene_start, conserved$gene_end)
hit_end <- pmax(conserved$gene_start, conserved$gene_end)

busco_hits <- GRanges(
  seqnames = conserved$contig,
  ranges = IRanges(start = hit_start, end = hit_end),
  strand = conserved$strand
)

# BUSCO cord to Braker gene -----------------------------------------------

# start and end of BUSCO and Braker might differ, but they will overlap
hits <- findOverlaps(busco_hits, genes)

# Check if sussy
# Each BUSCO matches one Braker
# It does not:
# matches no Braker;
# overlapps two braker genes;
# two BUSCO match to one Braker gene.
stopifnot(
  length(busco_hits) == length(hits),
  !anyDuplicated(queryHits(hits)),
  !anyDuplicated(subjectHits(hits)),
  identical(queryHits(hits), seq_along(busco_hits))
)

conserved$braker_gene_id <- genes$ID[subjectHits(hits)]

#conserved[, c("agic", "braker_gene_id")]

# get exons for each transcript -------------------------------------------

exon_rows <- gff[gff$type == "exon"]
head(mcols(exon_rows))

exons_by_tx <- split(exon_rows, exon_rows$transcript_id)

# I take only the first transcripts (.t1 not .t2 etc.) of a gene
tx_ids <- paste0(conserved$braker_gene_id, ".t1")
stopifnot(all(tx_ids %in% names(exons_by_tx)))

cons_exons_by_tx <- exons_by_tx[tx_ids]
tx_spans <- unlist(range(cons_exons_by_tx))
stopifnot(identical(names(tx_spans), tx_ids))

# add 1000bp flanks; keep inside contig -----------------------------------

contig_lengths <- seqlengths(seqinfo(genome_fa))

seqlengths(tx_spans) <- contig_lengths[seqlevels(tx_spans)]
stopifnot(!any(is.na(seqlengths(tx_spans))))

flank_bp <- 1000L
tx_spans_flanked <- tx_spans + flank_bp
# stoping if contig or chromosome ends before the flank range
stopifnot(all(width(tx_spans) + (2 * flank_bp) == width(tx_spans_flanked)))

# I could also allow trimming, and note it in FASTA header
# If there is nothing to trim, they remain same legnth
tx_spans_flanked_trimmed <- trim(tx_spans_flanked)
stopifnot(identical(tx_spans_flanked, tx_spans_flanked_trimmed))

# getSeq ------------------------------------------------------------------

flanked_seqs <- getSeq(genome_fa, tx_spans_flanked_trimmed)
stopifnot(
  length(flanked_seqs) == length(tx_spans_flanked_trimmed),
  identical(width(flanked_seqs), width(tx_spans_flanked_trimmed))
)

transcript_id <- names(tx_spans_flanked_trimmed)
A_tha_locus <- conserved$agic
contig <- as.character(seqnames(tx_spans_flanked_trimmed))
strand <- as.character(strand(tx_spans_flanked_trimmed))
range <- paste0(
  start(tx_spans_flanked_trimmed),
  "-",
  end(tx_spans_flanked_trimmed)
)

headers <- paste0(
  transcript_id, "|",
  A_tha_locus, "|",
  contig, "|",
  strand, "|",
  range, "|",
  "flank:", flank_bp
)

names(flanked_seqs) <- headers

writeXStringSet(flanked_seqs, "out/gene_region.Z_stapfii.fa")


# plot transcript exon fingerprint ----------------------------------------

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

for (tx_id in tx_ids) {

  current_tx <- tx_spans_flanked_trimmed[tx_id]
  region_start <- start(current_tx)
  current_tx_exons <- cons_exons_by_tx[[tx_id]]
  relative_exon_start <- start(current_tx_exons) - region_start + 1
  relative_exon_end <- end(current_tx_exons) - region_start + 1
  header <- headers[startsWith(headers, tx_id)]
  tx_plot <- plot_transcript(
    start = 1,
    end = width(current_tx),
    strand = as.character(strand(current_tx)),
    exon_starts = relative_exon_start,
    exon_ends = relative_exon_end,
    title = header
  )
  
  ggsave(filename = paste0(tx_id, ".png"),
         plot = tx_plot,
         path = "out/Z_stapfii_tx_plots",
         width = 10,
         height = 2.5,
         dpi = 300,
         create.dir = TRUE
  )
}
