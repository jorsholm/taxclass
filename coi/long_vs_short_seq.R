# Compare model performance using long and short sequences 

rm(list = ls())
setwd("coi")

# LOAD FUNCTIONS ---- 

source("load_FinBOL_GBOL.R")
source("functions.R")

# LOAD DATA ---- 

# Load test and train data 
data_long <- load_FinBOL_GBOL(level = "species", short = F)
data_short <- load_FinBOL_GBOL(level = "species", short = T)

# Replace labels of taxa unique to test
# Sort test data in alphabetical order
data_true <- get_data_true(data_long$test, data_long$train) |>
  dplyr::arrange(ID)

# READ RESULTS ---- 

# BayesANT 

BayesANT_long <- 
  read.table("models/results_bayesant_finbol_gbol/bayesant_test_finbol-gbol_species.txt", 
             header = T) |> 
  dplyr::arrange(ID)

BayesANT_short <- 
  read.table("models/results_bayesant_finbol_gbol/bayesant_testshort_finbol-gbol_species.txt", 
             header = T) |> 
  dplyr::arrange(ID)

# RDP 

RDP_long <- 
  read.table("models/results_rdp_finbol_gbol/rdp_test_finbol-gbol_species.txt", 
             header = T) |> 
  dplyr::arrange(ID)

RDP_short <- 
  read.table("models/results_rdp_finbol_gbol/rdp_testshort_finbol-gbol_species.txt", 
             header = T) |> 
  dplyr::arrange(ID)

# PROTAX 

protax_long <- 
  read.table("models/results_protax_finbol_gbol/protax_finbol_gbol_species.txt", 
             header = T) |> 
  dplyr::arrange(ID)
protax_long[is.na(protax_long)] <- 0
protax_long <- rename_PROTAX_output(protax_long)
colnames(protax_long) <- colnames(BayesANT_long) 

protax_short <- 
  read.table("models/results_protax_finbol_gbol/protax_finbol_gbol_species_short.txt", 
             header = T) |> 
  dplyr::arrange(ID)
protax_short[is.na(protax_short)] <- 0 
protax_short <- rename_PROTAX_output(protax_short)
colnames(protax_short) <- colnames(BayesANT_long)

# EPA-ng taxtree 

epatax_long <- 
  read.table("models/results_epa_finbol_gbol/epa_test_finbol-gbol_species.txt", 
             header = T) |> 
  dplyr::arrange(ID)
epatax_long <- rename_PROTAX_output(epatax_long)
colnames(epatax_long) <- colnames(BayesANT_long)

epatax_short <- 
  read.table("models/results_epa_finbol_gbol/epa_testshort_finbol-gbol_species.txt", 
             header = T) |> 
  dplyr::arrange(ID)
epatax_short <- rename_PROTAX_output(epatax_short)
colnames(epatax_short) <- colnames(BayesANT_long)

# sintax 

sintax_long <- 
  read.table("models/results_sintax_finbol_gbol/sintax_finbol_gbol_species.txt", 
             header = T) |> 
  dplyr::arrange(ID)
sintax_long <- rename_PROTAX_output(sintax_long)
colnames(sintax_long) <- colnames(BayesANT_long)

sintax_short <- 
  read.table("models/results_sintax_finbol_gbol/sintax_finbol_gbol_species_short.txt", 
             header = T) |> 
  dplyr::arrange(ID)
sintax_short <- rename_PROTAX_output(sintax_short)
colnames(sintax_short) <- colnames(BayesANT_long)

# PREPARE RESULTS ----

results <- list(
  "BayesANT_long" = BayesANT_long, 
  "BayesANT_short" = BayesANT_short, 
  "PROTAX_long" = protax_long, 
  "PROTAX_short" = protax_short, 
  "EPA-ngtaxtree_long" = epatax_long, 
  "EPA-ngtaxtree_short" = epatax_short, 
  "RDP_long" = RDP_long, 
  "RDP_short" = RDP_short, 
  "SINTAX_long" = sintax_long, 
  "SINTAX_short" = sintax_short
)

## Identify novel taxa 

observed_everywhere <- !grepl("_new$", data_true[, ncol(data_true)])

id_novel <- list()
for(rank in colnames(data_true)[-1]){
  id_novel[[rank]] <- which(stringr::str_ends(data_true[,rank], paste0(rank, "_new")))
}

id_observed <- list()
for(rank in colnames(data_true)[-1]){
  id_observed[[rank]] <- which(!stringr::str_ends(data_true[,rank], "_new"))
}

id_all <- lapply(data_true, function(x) 1:length(x))

# ANALYSIS ----

calibrations <- get_calibration(results, data_true, observed_everywhere)

accuracy <- get_accuracy_from_cal(calibrations)

library(tidyverse)

acc_comp <- 
  accuracy |> 
  separate(model, into = c("model", "sequence_length"), 
           sep = "_") |> 
  pivot_wider(names_from = sequence_length, 
              values_from = accuracy) |> 
  mutate(accuracy_drop = short - long) 

acc_comp$rank <- factor(acc_comp$rank, levels = names(data_true)[-1])

comp_plot <- 
  acc_comp |> 
  ggplot() + 
  geom_line(aes(x = rank, 
                y = accuracy_drop, 
                color = model, 
                group = model)) + 
  geom_point(aes(x = rank, 
                y = accuracy_drop, 
                color = model, 
                group = model)) +
  facet_wrap(~set) + 
  theme_bw() + 
  geom_hline(yintercept = 0, 
             linetype = "dashed")

ggsave(plot = comp_plot, 
       filename = "../plots/accuracy_change.pdf", 
       width = 300, 
       height = 130, 
       units = "mm")
