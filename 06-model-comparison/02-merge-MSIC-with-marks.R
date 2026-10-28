library(tidyverse)
library(stringr)

# functions ---------------------------------------------------------------

get_msic_by_tx_llmpipe <- function(tx, csv_path, id_file_path, suffix = "csv") {
  at_file_ids_pcad1 <- readLines(id_file_path)
  at_file_ids_pcad1_index <- make_tx_index(at_file_ids_pcad1)
  
  at_agic_tx_pcad1 <- read_csv(
    paste0(
      csv_path,
      at_file_ids_pcad1_index[tx],
      ".", suffix
    ),
    show_col_types = FALSE,
    progress = FALSE
  ) |>
    rename(msic = MSIC)
}

get_msic_by_tx_silvan <- function(dir, tx) {
  filepath <- file.path(dir, "msic.csv")
  msic_by_tx <- read_csv(
    filepath,
    show_col_types = FALSE,
    progress = FALSE
    ) |>
    mutate(tx = str_extract(header, "^[^|]+")) |>
    filter(tx == .env$tx)
}

path_constructor <- function(model, species, group) {
  if (group == "top15") {
    folder <- file.path("04-top15", "out")
  } else if (group == "tps") {
    folder <- file.path("05-additional-genes", "out")
  }
  
  
  id_file_path <- file.path(folder,
                            paste0(species, "_", group, "_", model, "_results"),
                            "fileIDs.txt")
  if (!file.exists(id_file_path)) {
    id_file_path <- NULL
  }
  
  out_dir <- file.path(folder,
                    paste0(species, "_", group, "_", model, "_results"))
  
  csv_path <- file.path(
    out_dir,
    paste0(species, "_", group, "_", model, "_scores_table")
  )
  
  list(id_file_path = id_file_path, out_dir = out_dir, csv_path = csv_path)
  
}


get_msic <- function(species, groups, models, keys) {
  compact(
    map(
      set_names(species),
      \(specie) {
        compact(
          map(
            set_names(groups),
            \(group) {
              compact(
                map(
                  set_names(models),
                  \(model) {
                    
                    path <- path_constructor(model, specie, group)
                    compact(
                      map(
                        set_names(keys[[group]][[specie]]),
                        \(tx) {
                          if (!is.null(path$id_file_path)) {
                            get_msic_by_tx_llmpipe(
                              tx,
                              path$csv_path,
                              path$id_file_path
                            )
                          } else if (file.exists(
                            file.path(path$out_dir, "msic.csv"))
                          ) {
                            get_msic_by_tx_silvan(
                              dir = path$out_dir,
                              tx = tx
                            )
                          } else {
                            NULL
                          }
                        }
                      )
                    )
                  }
                )
              )
            }
          )
        )
      }
    )
  )
}

rbind_models <- function(grp) {
  list_rbind(map(grp, \(m) list_rbind(m, names_to = "tx")), names_to = "model")
}

# setup -------------------------------------------------------------------

setwd("~/plantcad-zygophyllum-orthologs/")

source("04-top15/R/plot-helpers.R")

tps_marks <- readRDS('06-model-comparison/tps_marks.rds')
top15_marks <- readRDS('06-model-comparison/top15_marks.rds')

tx_pair_keys_tps <- read_tsv(
  "05-additional-genes/tps_orthos_in_busco.tsv"
) |>
  mutate(
    agic_tx = paste0(Atha, ".1"),
    bgc_tx = paste0(Zsta, ".t1")
  ) |>
  dplyr::select(agic_tx, bgc_tx)

tx_pair_keys_top15 <- read_tsv("04-top15/out/busco_top15.tsv") |>
  mutate(agic_tx = paste0(agic, ".1"),
         bgc_tx = paste0(braker_gene_id, ".t1")) |>
  dplyr::select(agic_tx, bgc_tx)

keys <- list(
  top15 = list(
    a_thaliana = tx_pair_keys_top15$agic_tx,
    z_stapffii = tx_pair_keys_top15$bgc_tx
  ),
  tps = list(
    a_thaliana = tx_pair_keys_tps$agic_tx,
    z_stapffii = tx_pair_keys_tps$bgc_tx
  )
)



# -------------------------------------------------------------------------

species <- c("a_thaliana", "z_stapffii")
groups <- c("top15", "tps")
models <- c("PCAD1", "PCAD2S", "PCAD2L", "Botanic1S", "Botanic1L")

msic <- get_msic(
  species,
  groups,
  models,
  keys)



long <- msic |>
  map(\(sp) list_rbind(map(sp, rbind_models), names_to = "group")) |>
  list_rbind(names_to = "species") |>
  mutate(position = row_number(), .by = c(species, group, model, tx))

long <- mutate(long, position = row_number(), .by = c(model, tx))

marks <- bind_rows(
  top15_marks$at_pos,
  top15_marks$zs_pos,
  tps_marks$at_pos,
  tps_marks$zs_pos
  )


joined <- long |>
  left_join(
    marks,
    by = c("tx" = "tx_id", "position" = "pos")
    ) |>
  mutate(species = case_when(
    str_detect(tx, "^AT") ~ "a_thaliana",
    str_detect(tx, "^g") ~ "z_stapffii",
    TRUE ~ NA
  ),
  is_busco = if_else(group == "top15", TRUE, FALSE)) |>
  select(
    tx, species, position, model, msic, in_exon, in_cds, is_busco, group,
    ref, p_a, p_t, p_c, p_g, p_ref,
    p_max, ic, everything()
    )


joined |>
  # filter(model == "PCAD1", in_exon == TRUE, species == "a_thaliana") |>
  filter(group == "tps") |>
  View()

setwd("06-model-comparison/")

write_csv(joined, "msic_long.csv")
