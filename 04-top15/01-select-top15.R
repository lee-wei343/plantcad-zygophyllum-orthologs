library(Rsamtools)
library(rtracklayer)
library(txdbmaker)
library(ggplot2)
library(pwalign)

setwd("~/plantcad-zygophyllum-orthologs/")

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


# Z stapffii --------------------------------------------------------------

# zs = Zygophyllum stapffii
zs_genome_path <- "data/Z_stapfii/Z_stapfii_genome.fa"
if (!file.exists(paste0(zs_genome_path, ".fai"))) {
  indexFa(zs_genome_path)
}
zs_genome_fa <- FaFile(zs_genome_path)

zs_gff <- import("data/Z_stapfii/braker.gff3")
zs_genes <- zs_gff[zs_gff$type == "gene"]


# quick look --------------------------------------------------------------

# contig_lengths <- seqlengths(seqinfo(genome_fa))
# 
# # n contigs
# length(contig_lengths)
# 
# which.max(contig_lengths)
# max(contig_lengths)
# 
# table(gff$type)


# BUSCO cords to GRanges --------------------------------------------------

# GRanges always start as lower, BUSCO gave higher as start when minus strand
hit_start <- pmin(busco$gene_start, busco$gene_end)
hit_end <- pmax(busco$gene_start, busco$gene_end)

busco_hits <- GRanges(
  seqnames = busco$contig,
  ranges = IRanges(start = hit_start, end = hit_end),
  strand = busco$strand
)

# BUSCO cord to Braker gene -----------------------------------------------

# start and end of BUSCO and Braker might differ, but they will overlap
hits <- findOverlaps(busco_hits, zs_genes)

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

busco$braker_gene_id <- zs_genes$ID[subjectHits(hits)]

# get exons for each transcript -------------------------------------------

zs_exon_rows <- zs_gff[zs_gff$type == "exon"]

zs_exons_by_tx <- split(zs_exon_rows, zs_exon_rows$transcript_id)

# I take only the first transcripts (.t1 not .t2 etc.) of a gene
zs_tx_ids <- paste0(busco$braker_gene_id, ".t1")
stopifnot(all(zs_tx_ids %in% names(zs_exons_by_tx)))


zs_exons_by_busco_tx <- zs_exons_by_tx[zs_tx_ids]
zs_tx_spans <- unlist(range(zs_exons_by_busco_tx))
stopifnot(identical(names(zs_tx_spans), zs_tx_ids))

# add 1000bp flanks; keep inside contig -----------------------------------

zs_contig_lengths <- seqlengths(seqinfo(zs_genome_fa))

seqlengths(zs_tx_spans) <- zs_contig_lengths[seqlevels(zs_tx_spans)]
stopifnot(!any(is.na(seqlengths(zs_tx_spans))))

zs_exon_spans <- zs_tx_spans + flank_bp
# stoping if contig or chromosome ends before the flank range
stopifnot(all(width(zs_tx_spans) + (2 * flank_bp) == width(zs_exon_spans)))

# I could also allow trimming, and note it in FASTA header
# If there is nothing to trim, they remain same legnth
zs_exon_spans_trim <- trim(zs_exon_spans)
stopifnot(identical(zs_exon_spans, zs_exon_spans_trim))

# Pairwise align ----------------------------------------------------------

# Pairwise global alignment
dna_matrix <- nucleotideSubstitutionMatrix(
  match = 2,
  mismatch = -1,
  baseOnly = TRUE
  )

alignment <- list()
pw_align_scores <- numeric(nrow(busco))
pid_align_scores <- numeric(nrow(busco))
for (i in seq_len(nrow(busco))) {
  at_tx_id <- paste0(busco$agic[i], ".1")
  at_tx <- at_exon_spans_trim[at_tx_id]
  at_tx_seq <- getSeq(at_genome, at_tx)
  
  zs_tx_id <- paste0(busco$braker_gene_id[i], ".t1")
  zs_tx <- zs_exon_spans_trim[zs_tx_id]
  zs_tx_seq <- getSeq(zs_genome_fa, zs_tx)
  
  global_align <- pairwiseAlignment(
    at_tx_seq,
    zs_tx_seq,
    type = "global",
    substitutionMatrix = dna_matrix,
    gapOpening = 5, gapExtension = 1
    )
  alignment[[paste0(at_tx_id, "_vs_", zs_tx_id)]] <- global_align
    
  pw_align_scores[i] <- score(global_align)
  pid_align_scores[i] <- pid(global_align)
}

busco$pw_align_scores <- pw_align_scores
busco$pid_align_scores <- pid_align_scores

busco_score_15 <- busco |>
  dplyr::slice_max(pw_align_scores, n = 15)

dir.create(path = "04-top15/out")




# -------------------------------------------------------------------------

at_exon_spans_trim <- at_exon_spans_trim[names(at_exon_spans_trim) %in% paste0(busco_score_15$agic, ".1")]

at_tx_seqs <- getSeq(at_genome, at_exon_spans_trim)
stopifnot(
  length(at_exon_spans_trim) == length(at_tx_seqs),
  identical(width(at_exon_spans_trim), width(at_tx_seqs))
)

at_tx_ids <- names(at_exon_spans_trim)
at_chrom  <- as.character(seqnames(at_exon_spans_trim))
at_strand <- as.character(strand(at_exon_spans_trim))
at_range  <- paste0(
  start(at_exon_spans_trim), "-",
  end(at_exon_spans_trim)
)

at_n_exons_by_tx <- lengths(at_exons_by_busco_tx[at_tx_ids])

headers <- paste0(
  at_tx_ids, "|",
  "chrom:", at_chrom, "|",
  at_strand, "|",
  at_range, "|",
  "n_exons:", at_n_exons_by_tx, "|",
  "flank:", flank_bp
)
stopifnot(length(headers) == length(at_tx_seqs))

names(at_tx_seqs) <- headers

if (write_fasta) {
  writeXStringSet(at_tx_seqs, filepath = "04-top15/out/tx_a_thaliana.fa")
}

# -------------------------------------------------------------------------

at_tx_list <- list()
current_tx <- at_exon_spans_trim[at_tx_id]
region_start <- start(current_tx)

current_tx_exons <- at_exons_by_busco_tx[[at_tx_id]]
relative_exon_start <- start(current_tx_exons) - region_start + 1
relative_exon_end   <- end(current_tx_exons) - region_start + 1

header <- headers[at_tx_ids == at_tx_id]
at_tx_list["A"] <- list(current_tx = current_tx,
                             start = 1,
                             end = width(current_tx),
                             relative_exon_start,
                             relative_exon_end,
                             header = header)


# -------------------------------------------------------------------------


at_tx_list <- list()
if (make_plots) {
  for (at_tx_id in at_tx_ids) {
    current_tx <- at_exon_spans_trim[at_tx_id]
    region_start <- start(current_tx)
    
    current_tx_exons <- at_exons_by_busco_tx[[at_tx_id]]
    relative_exon_start <- start(current_tx_exons) - region_start + 1
    relative_exon_end   <- end(current_tx_exons) - region_start + 1
    
    header <- headers[at_tx_ids == at_tx_id]
    
    tx_plot <- plot_transcript(
      start = 1,
      end = width(current_tx),
      strand = as.character(strand(current_tx)),
      exon_starts = relative_exon_start,
      exon_ends = relative_exon_end,
      title = header
    )
    
    at_tx_list[[at_tx_id]] <- list(current_tx = current_tx,
                             start = 1,
                             end = width(current_tx),
                             strand = as.character(strand(current_tx)),
                             relative_exon_start = relative_exon_start,
                             relative_exon_end = relative_exon_end,
                             header = header)
    
    
    ggsave(
      filename = paste0(at_tx_id, ".png"),
      plot = tx_plot,
      path = "04-top15/out/a_thaliana_tx_plots",
      width = 10,
      height = 2.5,
      dpi = 300,
      create.dir = TRUE
    )
  }
}


# getSeq ------------------------------------------------------------------

zs_exon_spans_trim <- zs_exon_spans_trim[names(zs_exon_spans_trim) %in% paste0(busco_score_15$braker_gene_id, ".t1")]

zs_tx_seqs <- getSeq(zs_genome_fa, zs_exon_spans_trim)
stopifnot(
  length(zs_tx_seqs) == length(zs_exon_spans_trim),
  identical(width(zs_tx_seqs), width(zs_exon_spans_trim))
)

zs_tx_ids <- names(zs_exon_spans_trim)

zs_gene_ids <- sub("\\.t1$", "", zs_tx_ids)
busco_row_for_tx <- match(zs_gene_ids, busco_score_15$braker_gene_id)
A_tha_locus <- busco_score_15$agic[busco_row_for_tx]
stopifnot(!anyNA(A_tha_locus))

zs_contig <- as.character(seqnames(zs_exon_spans_trim))
zs_strand <- as.character(strand(zs_exon_spans_trim))
zs_range  <- paste0(
  start(zs_exon_spans_trim), "-",
  end(zs_exon_spans_trim)
)
zs_n_exons_by_tx <- lengths(zs_exons_by_busco_tx[zs_tx_ids])

zs_headers <- paste0(
  zs_tx_ids, "|",
  A_tha_locus, "|",
  "ctg:", zs_contig, "|",
  zs_strand, "|",
  zs_range, "|",
  "n_exons:", zs_n_exons_by_tx, "|",
  "flank:", flank_bp
)

stopifnot(length(zs_headers) == length(zs_tx_seqs))

names(zs_tx_seqs) <- zs_headers

if (write_fasta) {
  writeXStringSet(zs_tx_seqs, filepath = "04-top15/out/tx_z_stapffii.fa")
}

zs_tx_list <- list()
if (make_plots) {
  for (zs_tx_id in zs_tx_ids) {
    current_tx <- zs_exon_spans_trim[zs_tx_id]
    region_start <- start(current_tx)
    
    current_tx_exons <- zs_exons_by_busco_tx[[zs_tx_id]]
    relative_exon_start <- start(current_tx_exons) - region_start + 1
    relative_exon_end   <- end(current_tx_exons) - region_start + 1
    
    header <- zs_headers[zs_tx_ids == zs_tx_id]
    
    tx_plot <- plot_transcript(
      start = 1,
      end = width(current_tx),
      strand = as.character(strand(current_tx)),
      exon_starts = relative_exon_start,
      exon_ends = relative_exon_end,
      title = header
    )
    
    zs_tx_list[[zs_tx_id]] <- list(current_tx = current_tx,
                             start = 1,
                             end = width(current_tx),
                             strand = as.character(strand(current_tx)),
                             relative_exon_start = relative_exon_start,
                             relative_exon_end = relative_exon_end,
                             header = header)
    
    ggsave(
      filename = paste0(zs_tx_id, ".png"),
      plot = tx_plot,
      path = "04-top15/out/z_stapffii_tx_plots",
      width = 10,
      height = 2.5,
      dpi = 300,
      create.dir = TRUE
    )
  }
}


# save the tx_lists -------------------------------------------------------

saveRDS(object = at_tx_list, "04-top15/out/at_tx_list.rds")
saveRDS(object = zs_tx_list, "04-top15/out/zs_tx_list.rds")

# -------------------------------------------------------------------------

top15_align_ids <- paste0(busco_score_15$agic, ".1", "_vs_", busco_score_15$braker_gene_id, ".t1")

top15_alignment <- alignment[top15_align_ids]

dir.create(path = "04-top15/out/palignments")


for (i in seq_len(length(top15_alignment))) {
  filename <- names(top15_alignment)[i]
  filepath <- paste0("04-top15/out/palignments/", filename)
  writePairwiseAlignments(top15_alignment[[i]],
                        file = filepath)
}
