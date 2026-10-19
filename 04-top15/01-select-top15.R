library(Rsamtools)
library(rtracklayer)
library(txdbmaker)
library(ggplot2)
library(pwalign)

setwd("~/plantcad-zygophyllum-orthologs/")

flank_bp <- 1000L
save_fasta <- TRUE
make_plots <- TRUE

source("R/extract-helpers.R")

# Arabidopsis thaliana ----------------------------------------------------

# at = Arabidopsis thaliana
# Arabidopsis Genome Initiative Code
busco <- readr::read_tsv("01-busco/out/conserved_go_annotated.tsv")

at_genome <- load_genome(
  genome_path = "data/A_thaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
)

at_exons_by_tx <- get_exons_by_tx(
  gff_path = "data/A_thaliana/Arabidopsis_thaliana.TAIR10.63.gff3"
)

at_exons_by_busco_tx <- at_exons_by_tx[
  names(at_exons_by_tx) %in% paste0(busco$agic, ".1")
  ]

at_exon_spans_trim <- add_flanks(
  exons_by_tx = at_exons_by_busco_tx,
  flank_bp = flank_bp,
  genome_fa = at_genome$fasta
)


# Z stapffii --------------------------------------------------------------

# zs = Zygophyllum stapffii
zs_genome <- load_genome(
  genome_path = "data/Z_stapfii/Z_stapfii_genome.fa",
  gff_path = "data/Z_stapfii/braker.gff3"
)

zs_genes <- zs_genome$gff[zs_genome$gff$type == "gene"]


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

zs_exons_by_tx <- split_exons_by_tx(zs_genome$gff)

# I take only the first transcripts (.t1 not .t2 etc.) of a gene
zs_tx_ids <- paste0(busco$braker_gene_id, ".t1")
stopifnot(all(zs_tx_ids %in% names(zs_exons_by_tx)))


zs_exons_by_busco_tx <- zs_exons_by_tx[zs_tx_ids]

# add 1000bp flanks; keep inside contig -----------------------------------

# I could also allow trimming, and note it in FASTA header
zs_exon_spans_trim <- add_flanks(
  exons_by_tx = zs_exons_by_busco_tx,
  flank_bp = flank_bp,
  genome_fa = zs_genome$fasta
)
stopifnot(identical(names(zs_exon_spans_trim), zs_tx_ids))

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
  at_tx_seq <- getSeq(at_genome$fasta, at_tx)

  zs_tx_id <- paste0(busco$braker_gene_id[i], ".t1")
  zs_tx <- zs_exon_spans_trim[zs_tx_id]
  zs_tx_seq <- getSeq(zs_genome$fasta, zs_tx)

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

dir.create(path = "04-top15/out", showWarnings = FALSE, recursive = TRUE)
readr::write_tsv(busco_score_15, "04-top15/out/busco_top15.tsv")


# top 15 A thaliana -------------------------------------------------------

at_exon_spans_trim <- at_exon_spans_trim[names(at_exon_spans_trim) %in% paste0(busco_score_15$agic, ".1")]

at_tx_seqs <- getSeq(at_genome$fasta, at_exon_spans_trim)
stopifnot(
  length(at_exon_spans_trim) == length(at_tx_seqs),
  identical(width(at_exon_spans_trim), width(at_tx_seqs))
)

headers <- make_header(
  exon_spans = at_exon_spans_trim,
  exons_by_tx = at_exons_by_busco_tx,
  flank_bp = flank_bp
)
names(at_tx_seqs) <- headers

if (save_fasta) {
  writeXStringSet(at_tx_seqs, filepath = "04-top15/out/tx_a_thaliana.fa")
}

at_tx_list <- list()
if (make_plots) {
  at_tx_list <- plot_transcripts(
    exon_spans = at_exon_spans_trim,
    exons_by_tx = at_exons_by_busco_tx,
    headers = headers,
    out_dir = "04-top15/out/a_thaliana_tx_plots"
  )
}


# top 15 Z stapffii -------------------------------------------------------

zs_exon_spans_trim <- zs_exon_spans_trim[names(zs_exon_spans_trim) %in% paste0(busco_score_15$braker_gene_id, ".t1")]

zs_tx_seqs <- getSeq(zs_genome$fasta, zs_exon_spans_trim)
stopifnot(
  length(zs_tx_seqs) == length(zs_exon_spans_trim),
  identical(width(zs_tx_seqs), width(zs_exon_spans_trim))
)

zs_tx_ids <- names(zs_exon_spans_trim)

zs_gene_ids <- sub("\\.t1$", "", zs_tx_ids)
busco_row_for_tx <- match(zs_gene_ids, busco_score_15$braker_gene_id)
A_tha_locus <- busco_score_15$agic[busco_row_for_tx]
stopifnot(!anyNA(A_tha_locus))

zs_headers <- make_header(
  exon_spans = zs_exon_spans_trim,
  exons_by_tx = zs_exons_by_busco_tx,
  flank_bp = flank_bp,
  seq_label = "ctg",
  extra = A_tha_locus
)
names(zs_tx_seqs) <- zs_headers

if (save_fasta) {
  writeXStringSet(zs_tx_seqs, filepath = "04-top15/out/tx_z_stapffii.fa")
}

zs_tx_list <- list()
if (make_plots) {
  zs_tx_list <- plot_transcripts(
    exon_spans = zs_exon_spans_trim,
    exons_by_tx = zs_exons_by_busco_tx,
    headers = zs_headers,
    out_dir = "04-top15/out/z_stapffii_tx_plots"
  )
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
