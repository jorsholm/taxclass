
rm(list = ls())

# LOAD FUNCTIONS ---------------------------------------------------------------

source("load_FinBOL_GBOL.R")
source("functions.R")

# LOAD DATA GENUS --------------------------------------------------------------

# Load test and train data 
data <- load_FinBOL_GBOL(level = "genus")

# Replace labels of taxa unique to test 
# Sort test data in alphabetical order 
data_true <- get_data_true(data$test, data$train) |>
  dplyr::arrange(ID)

# READ RESULTS GENUS -----------------------------------------------------------

# BayesANT 
result_BayesANT <- 
  read.table("models/results_bayesant_finbol_gbol/bayesant_test_finbol-gbol_genus.txt",
             header = T) |> 
  dplyr::arrange(ID)

# RDP 
result_RDP <- 
  read.table("models/results_rdp_finbol_gbol/rdp_test_finbol-gbol_genus.txt", 
             header = T) |> 
  dplyr::arrange(ID)

# PROTAX 
result_PROTAX <- read.table("models/results_protax_finbol_gbol/protax_finbol_gbol_genus.txt", 
                     header = T)
result_PROTAX[is.na(result_PROTAX)] <- 0
result_PROTAX <- rename_PROTAX_output(result_PROTAX)
# Quick check that order of columns is identical to BayesANT 
# all(stringr::str_to_lower(colnames(result_PROTAX)) ==
#       stringr::str_to_lower(colnames(result_BayesANT)))
# Rename Protax columns 
colnames(result_PROTAX) <- colnames(result_BayesANT)

# EPA-ng Taxonomy tree
result_epatax <- read.table("models/results_epa_finbol_gbol/epa_test_finbol-gbol_genus.txt",
                            header = T) |> 
  dplyr::arrange(ID)
result_epatax <- rename_PROTAX_output(result_epatax)
# all(stringr::str_to_lower(colnames(result_epatax)) ==
#       stringr::str_to_lower(colnames(result_BayesANT)))
colnames(result_epatax) <- colnames(result_BayesANT)

# PLOT RESULTS GENUS -----------------------------------------------------------

results <- list(
  "BayesANT" = result_BayesANT, 
  "PROTAX" = result_PROTAX, 
  "EPA-ng taxtree" = result_epatax, 
  "RDP" = result_RDP
)

# # Are all colnames in results identical? 
# length(unique(lapply(results, function(x) colnames(x)))) == 1
# # Are taxonomic columns identical to data_true? 
# all(colnames(data_true) == colnames(results[[1]])[1:ncol(data_true)])

plot_comparison(results, data_true, level = "Genus")

# LOAD DATA SPECIES ------------------------------------------------------------

# Load test and train data 
data <- load_FinBOL_GBOL(level = "species")

# Replace labels of taxa unique to test 
# Sort test data in alphabetical order 
data_true <- get_data_true(data$test, data$train) |>
  dplyr::arrange(ID)

# READ RESULTS SPECIES ---------------------------------------------------------

# BayesANT 
result_BayesANT <- 
  read.table("models/results_bayesant_finbol_gbol/bayesant_test_finbol-gbol_species.txt",
             header = T) |> 
  dplyr::arrange(ID)

# RDP 
result_RDP <- 
  read.table("models/results_rdp_finbol_gbol/rdp_test_finbol-gbol_species.txt", 
             header = T) |> 
  dplyr::arrange(ID)

# PROTAX 
result_PROTAX <- read.table("models/results_protax_finbol_gbol/protax_finbol_gbol_species.txt", 
                            header = T)
result_PROTAX[is.na(result_PROTAX)] <- 0
result_PROTAX <- rename_PROTAX_output(result_PROTAX)
# Quick check that order of columns is identical to BayesANT 
# all(stringr::str_to_lower(colnames(result_PROTAX)) ==
#       stringr::str_to_lower(colnames(result_BayesANT)))
# Rename Protax columns 
colnames(result_PROTAX) <- colnames(result_BayesANT)

# EPA-ng Taxonomy tree
result_epatax <- read.table("models/results_epa_finbol_gbol/epa_test_finbol-gbol_species.txt",
                            header = T) |> 
  dplyr::arrange(ID)
result_epatax <- rename_PROTAX_output(result_epatax)
# all(stringr::str_to_lower(colnames(result_epatax)) ==
#       stringr::str_to_lower(colnames(result_BayesANT)))
colnames(result_epatax) <- colnames(result_BayesANT)

# PLOT RESULT SPECIES ----------------------------------------------------------

results <- list(
  "BayesANT" = result_BayesANT, 
  "PROTAX" = result_PROTAX, 
  "EPA-ng taxtree" = result_epatax, 
  "RDP" = result_RDP
)

# # Are all colnames in results identical?
# length(unique(lapply(results, function(x) colnames(x)))) == 1
# # Are taxonomic columns identical to data_true?
# all(colnames(data_true) == colnames(results[[1]])[1:ncol(data_true)])

plot_comparison(results, data_true, level = "Species")

# How well-calibrated are different models across ranks? 
plotlist_rank_calibrations <- list()

for(i in 1:length(results)){
  plotlist_rank_calibrations[[i]] <- plot_ranks(results[[i]], data_true, 
                                                names(results)[i])
}

ggpubr::ggarrange(plotlist = plotlist_rank_calibrations, 
                  ncol = length(results)/2, nrow = 2)

id_observed <- !grepl("_new$", data_true[, ncol(data_true)])

plot_accuracies(results, data_true, id_observed)

bayesant_cal <- plot_cal_unknown(results$BayesANT, data_true, model_name = "BayesANT", id_observed, legend = "none")
protax_cal <- plot_cal_unknown(results$PROTAX, data_true, model_name = "PROTAX", id_observed, legend = "none") 
rdp_cal <- plot_cal_unknown(results$RDP, data_true, model_name = "RDP", id_observed, legend = "none") 
epatax_cal <- plot_cal_unknown(results$`EPA-ng taxtree`, data_true, model_name = "EPA-ng", id_observed, legend = "none") 

cal_plotlist <- list(bayesant_cal, protax_cal, rdp_cal, epatax_cal)

ggpubr::ggarrange(plotlist = cal_plotlist, ncol = 1)

