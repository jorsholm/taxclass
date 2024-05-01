
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

# Sintax
result_SINTAX <- read.table("models/results_sintax_finbol_gbol/sintax_finbol_gbol_genus.txt",
                            header = TRUE)
result_SINTAX <- rename_PROTAX_output(result_SINTAX)
result_SINTAX <- dplyr::arrange(result_SINTAX, ID)
# all(stringr::str_to_lower(colnames(result_SINTAX)) ==
#       stringr::str_to_lower(colnames(result_BayesANT)))
colnames(result_SINTAX) <- colnames(result_BayesANT)

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

# SINTAX
result_SINTAX <- read.table("models/results_sintax_finbol_gbol/sintax_finbol_gbol_species.txt",
                            header = T)
result_SINTAX[is.na(result_SINTAX)] <- 0
result_SINTAX <- rename_PROTAX_output(result_SINTAX)
result_SINTAX <- dplyr::arrange(result_SINTAX, ID)
# Quick check that order of columns is identical to BayesANT
# all(stringr::str_to_lower(colnames(result_SINTAX)) ==
#       stringr::str_to_lower(colnames(result_BayesANT)))
# Rename SINTAX columns
colnames(result_SINTAX) <- colnames(result_BayesANT)

# PREPARE AND CHECK RESULT -----------------------------------------------------

results <- list(
  "BayesANT" = result_BayesANT,
  "PROTAX" = result_PROTAX,
  "EPA-ng taxtree" = result_epatax,
  "RDP" = result_RDP,
  "SINTAX" = result_SINTAX
)

id_observed <- !grepl("_new$", data_true[, ncol(data_true)])

# Some quick checks 
if(!length(unique(lapply(results, function(x) colnames(x)))) == 1) print("One of the results data frames have wrong column names.")
if(!all(sapply(results, function(x) all(x[,1] == data_true[,1])))) print("One of the results data frames is not sorted by ID.")
if(!all(colnames(data_true) == colnames(results[[1]])[1:ncol(data_true)])) print("Column names containing taxonomic information not identical to data_true")

# PLOT CALIBRATION CURVES ------------------------------------------------------

calibrations <- get_calibration(results, data_true, id_observed)

calibrations$rank <- factor(calibrations$rank, 
                            levels = colnames(data_true)[-1])
calibrations$set <- factor(calibrations$set, 
                           levels = c("All", "Observed", "Novel"))

#### Plot all calibration curves (all ranks, taxa sets and models)
p_cal <- plot_calibration(calibrations) 

ggplot2::ggsave(plot = p_cal, 
                filename = "../plots/calibration_COI_full.png", 
                width = 210, 
                height = 297, 
                units = "mm")

# PLOT ACCURACIES --------------------------------------------------------------

plot_accuracies(results, data_true, id_observed)


