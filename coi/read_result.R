# Read and modify model results, because it is becoming too much for the same script 

rm(list = ls())
setwd("coi")

# LOAD FUNCTIONS ---------------------------------------------------------------

source("load_FinBOL_GBOL.R")
source("functions.R")

# SET PARAMETERS ---------------------------------------------------------------

short <- F
shorttxt <- ""
undshort <- ""
if(short){
  shorttxt <- "short"
  undshort <- "_short"
}
keep_na <- F
natxt <- ""
if(keep_na) natxt <- "_keepNA"

# LOAD DATA --------------------------------------------------------------------

# Load test and train data 
data <- load_FinBOL_GBOL()

# Replace labels of taxa unique to test
# Sort test data in alphabetical order
data_true <- get_data_true(data$test, data$train) |>
  dplyr::arrange(ID)

class_order_match <- data$train |> 
  dplyr::select(Class, Order) |> 
  dplyr::distinct()

train_tax <- 
  read.table("data/train_tax.tsv", 
             header = T, 
             fill = T) |> 
  dplyr::distinct(kingdom, phylum, class, order, family, genus, species) |> 
  dplyr::rename(Class = kingdom, 
                Order = phylum, 
                Family = class, 
                Subfamily = order, 
                Tribe = family, 
                Genus = genus, 
                Species = species)


# CORRECT COLUMN NAMES AND ORDER -----------------------------------------------

ranks <- colnames(data_true)[-1]
correct_cols <- c("ID", 
                  ranks, 
                  sapply(ranks, function(x) paste0("Prob_", x), USE.NAMES = F))

# READ RESULTS -----------------------------------------------------------------

# BayesANT
result_BayesANT <-
  read.table(paste0("results/bayesant/bayesant_test", shorttxt, "_nt_aln_label.txt"),
             header = T) |>
  tibble::rownames_to_column("ID") |>
  dplyr::arrange(ID)
result_BayesANT <- arrange_columns(result_BayesANT, correct_cols)

# RDP
result_RDP <-
  read.table(paste0("results/rdp/rdp_test", shorttxt, "_nt_aln_label.txt"),
             header = T) |>
  tibble::rownames_to_column("ID") |>
  dplyr::arrange(ID)
result_RDP <- arrange_columns(result_RDP, correct_cols)

# PROTAX
# result_PROTAX <- 
#   read.table(paste0("models/results_protax_finbol_gbol/protax_finbol_gbol_",
#                     level, undshort, ".txt"),
#              header = T)
# result_PROTAX[is.na(result_PROTAX)] <- 0
# result_PROTAX <- rename_PROTAX_output(result_PROTAX)
# # Quick check that order of columns is identical to BayesANT
# # all(stringr::str_to_lower(colnames(result_PROTAX)) ==
# #       stringr::str_to_lower(colnames(result_BayesANT)))
# # Rename Protax columns
# colnames(result_PROTAX) <- colnames(result_BayesANT)

# EPA-ng Taxonomy tree
result_epatax <-
  read.table(paste0("results/epang_taxtree/epang_taxtree_test", shorttxt,
                    "_nt_all.tsv"),
             header = T) |>
  dplyr::arrange(ID)
result_epatax <- rename_unk_output(result_epatax)
result_epatax <- arrange_columns(result_epatax, correct_cols)

# Sintax
result_SINTAX <-
  read.table(paste0("results/sintax/sintax_test", shorttxt, "_nt_16.tsv"),
             header = TRUE) |> 
  dplyr::arrange(ID)
result_SINTAX <- arrange_columns(result_SINTAX, correct_cols)

# EPA-ng phylogenetic tree 
result_epaphyl <- 
  read.table(paste0("results/epang_phyltree/epang_phyltree_test",
                    shorttxt, "_finbol-gbol.txt"),
             header = T, sep = "\t")  |>
  dplyr::arrange(ID)
result_epaphyl$class[which(result_epaphyl$class == "DISTANT")] <- "unk"
result_epaphyl$species <- sapply(result_epaphyl$species, 
                                 function(x) stringr::str_replace(x, " ", "_"),
                                 USE.NAMES = F)
result_epaphyl <- rename_unk_output(result_epaphyl)
result_epaphyl <- arrange_columns(result_epaphyl, correct_cols)

# MycoAI-CNN 
result_aicnn <-
  read.table(paste0("results/mycoai_cnn/mycoai_cnn_test", shorttxt,
                    "_nt_gpu.tsv"),
             header = T) |>
  dplyr::arrange(ID)
result_aicnn <- 
  result_aicnn |> 
  dplyr::left_join(class_order_match, 
                   by = c("order" = "Order")) |> 
  dplyr::mutate(Prob_Class = Prob_order)
result_aicnn <- arrange_columns(result_aicnn, correct_cols)

# MycoAI-BERT 
result_aibert <- 
  read.table(paste0("results/mycoai_bert/mycoai_bert_test", shorttxt, 
                    "_nt_gpu.tsv"), 
             header = T) |> 
  dplyr::arrange(ID)
result_aibert <- 
  result_aibert |> 
  dplyr::left_join(class_order_match, 
                   by = c("order" = "Order")) |> 
  dplyr::mutate(Prob_Class = Prob_order)
result_aibert <- arrange_columns(result_aibert, correct_cols)

# BLAST top hit 
result_blast_top <-
  read.table(paste0("results/blast/blast_top_hit_test", shorttxt, "_nt_16.tsv"),
             header = T) |> 
  dplyr::arrange(ID) |> 
  dplyr::mutate(species = purrr::map_chr(species, 
                                         ~stringr::str_extract(.x, "^[^;]+")))
result_blast_top <- arrange_columns(result_blast_top, correct_cols)

## HERE COMES THE ALGORITHMS WHICH HAVE TWO OPTIONS: (1) KEEP NA AS NA or (2) USE NA AS NOVEL 

# BLAST threshold
result_blast_thresh <-
  read.table(paste0("results/blast/blast_thresh_test", shorttxt, "_nt_16.tsv"),
             header = T) |> 
  dplyr::arrange(ID)
result_blast_thresh <- arrange_columns(result_blast_thresh, correct_cols)
if(!keep_na){
  result_blast_thresh[is.na(result_blast_thresh)] <- "unk"
  result_blast_thresh[,correct_cols[which(stringr::str_starts(correct_cols, "Prob_"))]] <- 1
  result_blast_thresh <- rename_unk_output(result_blast_thresh)
}

# DNA-barcoder
result_dnabarcoder <-
  read.table(paste0("results/dnabarcoder/test", shorttxt, "_nt_4.tsv"),
             header = T) |> 
  dplyr::arrange(ID)
if(keep_na){
  # Unidentified = NA 
  result_dnabarcoder[result_dnabarcoder == "unidentified"] <- NA
  result_dnabarcoder <- 
    result_dnabarcoder |>
    tidyr::pivot_longer(cols = all_of(stringr::str_to_lower(ranks)), 
                        names_to = "rank", values_to = "taxon") |> 
    dplyr::mutate(prob = dplyr::if_else(is.na(taxon), NA, 1)) |> 
    tidyr::pivot_wider(names_from = rank, 
                       values_from = c(taxon, prob)) |> 
    dplyr::rename_with(~gsub("taxon_", "", .), dplyr::starts_with("taxon_"))
}else{
  result_dnabarcoder <-
    result_dnabarcoder |>
    rename_unk_output(unktxt = "unidentified") |>
    cbind(purrr::map_dfc(ranks, ~ dplyr::tibble(!!paste0("Prob_", .x) := 1)))
}
result_dnabarcoder <- arrange_columns(result_dnabarcoder, correct_cols)

# IDTAXA 
result_idtaxa <-
 read.table(paste0("results/idtaxa/idtaxa_test", shorttxt, "_nt_4.tsv"), 
            fill = T, header = T) |> 
  dplyr::arrange(ID)
result_idtaxa <- arrange_columns(result_idtaxa, correct_cols)
if(keep_na){
  # Unclassified = NA 
  result_idtaxa <- 
    result_idtaxa |> 
    dplyr::mutate(across(all_of(ranks), 
                         ~stringr::str_replace(., "^unclassified_.*", ""))) |> 
    tidyr::pivot_longer(cols = all_of(ranks), 
                        names_to = "rank", 
                        values_to = "taxon") |> 
    tidyr::pivot_longer(cols = paste0("Prob_", ranks), 
                        names_to = "Prob_rank", 
                        names_pattern = "Prob_(.*)", 
                        values_to = "Prob") |> 
    dplyr::filter(rank == Prob_rank) |> 
    dplyr::mutate(Prob = dplyr::if_else(taxon == "", NA, Prob)) |> 
    dplyr::select(-Prob_rank) |> 
    tidyr::pivot_wider(names_from = rank, 
                       values_from = c(taxon, Prob)) |> 
    dplyr::rename_with(~gsub("taxon_", "", .), dplyr::starts_with("taxon_")) |> 
    dplyr::mutate(across(all_of(ranks), 
                         ~dplyr::na_if(., "")))
}else{
  result_idtaxa <-
    result_idtaxa |>
    dplyr::mutate(across(all_of(ranks),
                         ~dplyr::na_if(., ""))) |>
    dplyr::mutate(across(all_of(ranks),
                         ~dplyr::coalesce(stringr::str_replace(., "^unclassified_.*", "unk"), "unk"))) |> 
    tidyr::pivot_longer(cols = all_of(correct_cols[which(stringr::str_starts(correct_cols, "Prob_"))]), 
                        names_to = "rank", values_to = "prob") |> 
    dplyr::group_by(ID) |> 
    dplyr::mutate(prob = dplyr::if_else(is.na(prob), dplyr::lag(prob, default = FALSE), prob)) |>
    tidyr::fill(prob, .direction = "down") |> 
    dplyr::ungroup() |> 
    tidyr::pivot_wider(names_from = rank, values_from = prob)
  result_idtaxa <- rename_unk_output(result_idtaxa)
}

# Crest4
result_crest4 <- 
  read.table(paste0("results/crest4/crest4_test", shorttxt, "_nt_4.tsv"), 
             header = T)
result_crest4[result_crest4 == "unclassified"] <- NA
result_crest4[paste0("Prob_", ranks)] <- 1
result_crest4 <- arrange_columns(result_crest4, correct_cols)


# Mystery result ------
# # Brendan's mystery model
# # Special keep_na case: 
# # Unknown = novel (has probability), NA = NA (no probability)
# # if the rank above is unknown, then the whole row should be unknown with the same probability.
# # When not keep_na, the same rule applies for all unknown and NAs. 
# result_mystery <- read.table("models/mystery/mystery_model_finbol_gbol.txt",
#                              header = T) |> 
#   dplyr::arrange(ID)
# result_mystery <- arrange_columns(result_mystery, correct_cols)
# if(keep_na){
#   result_mystery <- 
#     result_mystery |> 
#     dplyr::mutate(across(all_of(ranks), 
#                          ~stringr::str_replace(., "^unknown_.*", "unk"))) |> 
#     tidyr::pivot_longer(cols = all_of(ranks),
#                         names_to = "rank",
#                         values_to = "taxon") |>
#     tidyr::pivot_longer(cols = paste0("Prob_", ranks),
#                         names_to = "Prob_rank",
#                         names_pattern = "Prob_(.*)",
#                         values_to = "Prob") |> 
#     dplyr::filter(rank == Prob_rank) |> 
#     dplyr::group_by(ID) |> 
#     dplyr::mutate(taxon = dplyr::if_else(is.na(taxon) & any(taxon == "unk"), "unk", taxon)) |> 
#     dplyr::mutate(Prob_filled = dplyr::if_else(taxon == "unk" & is.na(Prob), dplyr::lag(Prob, default = NA), Prob)) |> 
#     tidyr::fill(Prob_filled, .direction = "down") |> 
#     dplyr::mutate(Prob = dplyr::if_else(taxon == "unk" & is.na(Prob), Prob_filled, Prob)) |> 
#     dplyr::select(-Prob_filled) |> 
#     dplyr::ungroup() |> 
#     dplyr::select(-Prob_rank) |> 
#     tidyr::pivot_wider(names_from = rank,
#                        values_from = c(taxon, Prob)) |> 
#     dplyr::rename_with(~gsub("taxon_", "", .), dplyr::starts_with("taxon_")) 
# }else{
#   result_mystery <- 
#     result_mystery |> 
#     dplyr::mutate(across(all_of(ranks), 
#                          ~stringr::str_replace(., "^unknown_.*", "unk"))) |> 
#     tidyr::pivot_longer(cols = all_of(ranks),
#                         names_to = "rank",
#                         values_to = "taxon") |>
#     tidyr::pivot_longer(cols = paste0("Prob_", ranks),
#                         names_to = "Prob_rank",
#                         names_pattern = "Prob_(.*)",
#                         values_to = "Prob") |> 
#     dplyr::filter(rank == Prob_rank) |> 
#     dplyr::group_by(ID) |> 
#     dplyr::mutate(taxon = dplyr::if_else(is.na(taxon), "unk", taxon)) |>  
#     dplyr::mutate(Prob = dplyr::if_else(is.na(Prob), dplyr::lag(Prob, default = NA), Prob)) |> 
#     tidyr::fill(Prob, .direction = "down") |> 
#     dplyr::ungroup() |> 
#     dplyr::select(-Prob_rank) |> 
#     tidyr::pivot_wider(names_from = rank,
#                        values_from = c(taxon, Prob)) |> 
#     dplyr::rename_with(~gsub("taxon_", "", .), dplyr::starts_with("taxon_")) 
# }
# result_mystery <- rename_unk_output(result_mystery)

# UGLY FIX FOR BLAST MISSING DATA ----------------------------------------------

add_blast_thresh <- data.frame(ID = data_true$ID[which(!(data_true$ID %in% result_blast_thresh$ID))])
add_blast_thresh[,correct_cols[-1]] <- NA

result_blast_thresh <- rbind(result_blast_thresh, add_blast_thresh) |> dplyr::arrange(ID)

add_blast_top <- data.frame(ID = data_true$ID[which(!(data_true$ID %in% result_blast_top$ID))])
add_blast_top[,correct_cols[-1]] <- NA

result_blast_top <- rbind(result_blast_top, add_blast_top) |> dplyr::arrange(ID)

results <- list(
  "BayesANT" = result_BayesANT,
  #  "PROTAX" = result_PROTAX,
  "EPA-ng taxtree" = result_epatax,
  "RDP" = result_RDP,
  "BLAST top hit" = result_blast_top, 
  "BLAST threshold" = result_blast_thresh, 
  "SINTAX" = result_SINTAX,
  # "Mystery" = result_mystery,
  "EPA-ng phyltree" = result_epaphyl,
  "DNABarcoder" = result_dnabarcoder,
  "IDTAXA" = result_idtaxa, 
  "MycoAI-CNN" = result_aicnn, 
  "MycoAI-BERT" = result_aibert, 
  "Crest4" = result_crest4
)

results <- lapply(results, function(x) as.data.frame(x))

saveRDS(results, paste0("result_list", natxt, ".rds"))

# CHECK TAXONOMY ---------------------------------------------------------------

mismatch_taxonomy <- check_taxonomy(taxdf = train_tax, results = results, ranks = ranks)

library(tidyverse)

mismatch_taxonomy |> 
  as_tibble() |> 
  filter(str_starts(model, "MycoAI") | model == "SINTAX") |> 
  mutate(intid = row_number()) |> 
  left_join(train_tax, by = "Species", suffix = c("_pred", "_train")) |> 
  relocate(intid, .before = Class_pred) |> 
  relocate(model, .after = intid) |> 
  relocate(Species, .after = model) |>  
  pivot_longer(4:15, 
               names_to = c("rank", "set"), 
               values_to = "taxa", 
               names_sep = "_") |> 
  pivot_wider(names_from = set, values_from = taxa) |> 
  filter(pred != train)

 mismatch_taxonomy |> 
   as_tibble() |> 
   filter(str_starts(model, "MycoAI") | model == "SINTAX")






