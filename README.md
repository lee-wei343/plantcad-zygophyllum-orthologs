# plantcad-zygophyllum-orthologs

How confident are plant DNA language models on a genome they have never seen?
This project compares the masked-nucleotide prediction confidence of **PlantCAD** (PlantCAD1, PlantCAD2-small/-large) and **Botanic1** (small/large) on ortholog genes of the Namib desert shrub _Zygophyllum stapffii_, which is outside every model's training data, and of _Arabidopsis thaliana_, which is inside all of them.

Confidence is measured per base as the maximum scaled information content (MSIC), where 1 means low entropy and the model prediction matches the reference base.

The models were run outside this project, indicated by markdown files.

## Key findings

- Botanic1 performance doesn't drop, but it is less confident than PlantCaduceus.

## Data

Genomes, annotations and results are not included. The scripts expect them under `data/<Species>/`:

- _A. thaliana_ (TAIR10/Araport11), _O. sativa_ (IRGSP-1.0) and
  _P. patens_ (Phypa V3): Ensembl Plants release 63
- _T. mongolica_: genome and annotation from Liu et al. (2023)
- _Z. stapffii_: unpublished PacBio assembly and BRAKER annotation
  (P. Karunarathne); not publicly available

## Software

BUSCO 6.1.0, OrthoFinder 3.1.5, BLAST 2.17.0, R 4.5.2 with Bioconductor 3.22 (Biostrings, GenomicRanges, rtracklayer, txdbmaker, pwalign, ape) and the tidyverse.
