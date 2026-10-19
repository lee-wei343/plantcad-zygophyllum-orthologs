library(Rsamtools)
library(rtracklayer)
library(txdbmaker)
library(ggplot2)
library(pwalign)

setwd("~/plantcad-zygophyllum-orthologs/")
source("03-candidate-transcripts/R/extract-helpers.R")


# setup -------------------------------------------------------------------

flank_bp <- 1000L
write_fasta <- TRUE
make_plots <- TRUE


# generally load ----------------------------------------------------------

busco <- readr::read_tsv("01-busco/out/conserved_go_annotated.tsv")

# Arabidopsis thaliana ----------------------------------------------------

# at = Arabidopsis thaliana
# Arabidopsis Genome Initiative Code


# load AT -----------------------------------------------------------------

at_genome_path <- "data/A_thaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
at_gff_path <- "data/A_thaliana/Arabidopsis_thaliana.TAIR10.63.gff3"


at_genome <- load_genome(genome_path = at_genome_path, gff_path = at_gff_path)

at_txdb <- makeTxDbFromGFF(file = at_gff_path)
at_exons_by_tx <- exonsBy(at_txdb, by = "tx", use.name = TRUE)


# Filter AT ---------------------------------------------------------------

at_exons_by_busco_tx <- at_exons_by_tx[
  names(at_exons_by_tx) %in% paste0(busco$agic, ".1")
]


# get Seq with header -----------------------------------------------------

# GRanges
at_ex_by_tx_flank <- add_flanks(
  exons_by_tx = at_exons_by_busco_tx,
  flank_bp = flank_bp,
  genome_fa = at_genome$fasta
  )

# DNA-string set
at_tx_seqs <- getSeq(at_genome$fasta, at_ex_by_tx_flank)
stopifnot(
  length(at_ex_by_tx_flank) == length(at_tx_seqs),
  identical(width(at_ex_by_tx_flank), width(at_tx_seqs))
)

headers <- make_header(
  exon_spans = at_ex_by_tx_flank,
  exons_listed = at_exons_by_busco_tx,
  flank_bp = flank_bp
)
names(at_tx_seqs) <- headers


# write FASTA -------------------------------------------------------------

if (write_fasta) {
  writeXStringSet(
    at_tx_seqs,
    filepath = "03-candidate-transcripts/out/tx_a_thaliana.fa"
    )
}


# Plot Seqs ---------------------------------------------------------------

if (make_plots) {
  for (at_tx_id in names(at_ex_by_tx_flank)) {
    current_tx <- at_ex_by_tx_flank[at_tx_id]
    region_start <- start(current_tx)
    current_tx_exons <- at_exons_by_busco_tx[[at_tx_id]]
    relative_exon_start <- start(current_tx_exons) - region_start + 1
    relative_exon_end <- end(current_tx_exons) - region_start + 1
    header <- headers[startsWith(headers, at_tx_id)]
    
    tx_plot <- plot_transcript(
      start = 1,
      end = width(current_tx),
      strand = as.character(strand(current_tx)),
      exon_starts = relative_exon_start,
      exon_ends = relative_exon_end,
      title = header
    )

    ggsave(
      filename = paste0(at_tx_id, ".png"),
      plot = tx_plot,
      path = "03-candidate-transcripts/out/a_thaliana_tx_plots",
      width = 10,
      height = 2.5,
      dpi = 300,
      create.dir = TRUE
    )
  }
}


# Z stapffii --------------------------------------------------------------

# zs = Z stapffii
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

# Check if suspicious
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

# getSeq ------------------------------------------------------------------

zs_tx_seqs <- getSeq(zs_genome_fa, zs_exon_spans_trim)
stopifnot(
  length(zs_tx_seqs) == length(zs_exon_spans_trim),
  identical(width(zs_tx_seqs), width(zs_exon_spans_trim))
)


zs_tx_ids <- names(zs_exon_spans_trim)
A_tha_locus <- busco$agic
zs_contig <- as.character(seqnames(zs_exon_spans_trim))
zs_strand <- as.character(strand(zs_exon_spans_trim))
zs_range <- paste0(
  start(zs_exon_spans_trim),
  "-",
  end(zs_exon_spans_trim)
)
zs_n_exons_by_tx <- lengths(zs_exons_by_busco_tx)

zs_headers <- paste0(
  zs_tx_ids, "|",
  A_tha_locus, "|",
  "ctg:", zs_contig, "|",
  zs_strand, "|",
  zs_range, "|",
  "n_exons:", zs_n_exons_by_tx, "|",
  "flank:", flank_bp
)

names(zs_tx_seqs) <- zs_headers

if (write_fasta) {
  writeXStringSet(zs_tx_seqs, "03-candidate-transcripts/out/tx_z_stapffii.fa")
}

if (make_plots) {
  for (zs_tx_id in zs_tx_ids) {
    current_tx <- zs_exon_spans_trim[zs_tx_id]
    region_start <- start(current_tx)
    current_tx_exons <- zs_exons_by_busco_tx[[zs_tx_id]]
    relative_exon_start <- start(current_tx_exons) - region_start + 1
    relative_exon_end <- end(current_tx_exons) - region_start + 1
    header <- zs_headers[startsWith(zs_headers, zs_tx_id)]
    tx_plot <- plot_transcript(
      start = 1,
      end = width(current_tx),
      strand = as.character(strand(current_tx)),
      exon_starts = relative_exon_start,
      exon_ends = relative_exon_end,
      title = header
    )

    ggsave(
      filename = paste0(zs_tx_id, ".png"),
      plot = tx_plot,
      path = "03-candidate-transcripts/out/z_stapfii_tx_plots",
      width = 10,
      height = 2.5,
      dpi = 300,
      create.dir = TRUE
    )
  }
}


# Pairwise align ----------------------------------------------------------

# Pairwise global alignment
dna_matrix <- nucleotideSubstitutionMatrix(
  match = 2,
  mismatch = -1,
  baseOnly = TRUE
)

pw_align_scores <- numeric(nrow(busco))
pid_align_scores <- numeric(nrow(busco))
for (i in seq_len(nrow(busco))) {
  at_tx <- at_ex_by_tx_flank[paste0(busco$agic[i], ".1")]
  at_tx_seq <- getSeq(at_genome$fasta, at_tx)

  zs_tx <- zs_exon_spans_trim[paste0(busco$braker_gene_id[i], ".t1")]
  zs_tx_seq <- getSeq(zs_genome_fa, zs_tx)

  global_align <- pairwiseAlignment(
    at_tx_seq,
    zs_tx_seq,
    type = "global",
    substitutionMatrix = dna_matrix,
    gapOpening = 5, gapExtension = 1
  )
  pw_align_scores[i] <- score(global_align)
  pid_align_scores[i] <- pid(global_align)
}

busco$pw_align_scores <- pw_align_scores
busco$pid_align_scores <- pid_align_scores


readr::write_tsv(busco, file = "03-candidate-transcripts/out/conserved_go_annotated_braker.tsv")
