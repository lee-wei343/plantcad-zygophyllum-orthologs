library(tidyverse)

# Check BUSCOxBraker ------------------------------------------------------


busco <- read_tsv("out/conserved_go_annotated_braker.tsv")
ortho_zs_at <- read_tsv("02-orthofinder/Results_Sep29/Orthologues/Orthologues_Z_stapffii/Z_stapffii__v__A_thaliana.tsv")

ortho_pairs <- ortho_zs_at |>
  separate_longer_delim(Z_stapffii, ", ") |>
  separate_longer_delim(A_thaliana, ", ")

ortho_pairs <- ortho_pairs |>
  mutate(Z_stapffii = sub("\\.t[0-9]+$", "", Z_stapffii))

ortho_keys <- paste(ortho_pairs$A_thaliana, ortho_pairs$Z_stapffii)
busco_keys <- paste(busco$agic, busco$braker_gene_id)

res <- busco_keys %in% ortho_keys

busco$buscopair_in_orthofinderpair <- res

write_tsv(busco, "out/conserved_go_annotated_braker_pairchecked.tsv")
