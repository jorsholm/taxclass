
rm(list = ls())

# LOAD FUNCTIONS AND PACKAGES --------------------------------------------------

library(tidyverse)
source("functions.R")

case <- "coi"

ranks <- c("Class", 
           "Order",
           "Family",
           "Subfamily",
           "Tribe",
           "Genus",
           "Species")

# LOAD DATA --------------------------------------------------------------------

# Load test and train data
data <- load_train_test(train_file = paste0(case, "/data/train_tax.tsv"), 
                        test_file = paste0(case, "/data/test_tax.tsv"),
                        ranks = ranks)

# Replace labels of taxa unique to test
# Sort test data in alphabetical order
data_true <- get_data_true(data$test, data$train, 
                           ranks = ranks) |>
  dplyr::arrange(ID)

# CORRECT COLUMN NAMES AND ORDER -----------------------------------------------

correct_cols <- c("ID",
                  ranks,
                  sapply(ranks, function(x) paste0("Prob_", x), USE.NAMES = F))

# READ RESULTS -----------------------------------------------------------------

nt_long <-  readRDS("coi/result_list_nt.rds")
aa_long <- readRDS("coi/result_list_aa.rds")

# IDENTIFY TAXA SETS -----------------------------------------------------------

# Observed on all ranks ("Observed species")
observed_everywhere <- !grepl("_new$", data_true[, ncol(data_true)])

# Novel at each rank ("Novel taxa")
id_novel <- list()
for(rank in colnames(data_true)[-1]){
  id_novel[[rank]] <- which(stringr::str_ends(data_true[,rank], paste0(rank, "_new")))
}

# Observed at each rank ("Observed taxa")
id_observed <- list()
for(rank in colnames(data_true)[-1]){
  id_observed[[rank]] <- which(!stringr::str_ends(data_true[,rank], "_new"))
}

id_all <- lapply(data_true, function(x) 1:length(x))

# TODO: Some of our results have more sequences than data_true
# Here, I remove them as a temporary solution 
nt_long <- lapply(nt_long, function(x) x[which(x[,1] %in% data_true[,1]),])
aa_long <- lapply(aa_long, function(x) x[which(x[,1] %in% data_true[,1]),])

# Some quick checks of results structure 
if(!length(unique(lapply(aa_long, function(x) colnames(x)))) == 1) print("One of the results data frames have wrong column names.")
if(!all(sapply(aa_long, function(x) all(x[,1] == data_true[,1])))) print("One of the results data frames is not sorted by ID.")
if(!all(colnames(data_true) == colnames(aa_long[[1]])[1:ncol(data_true)])) print("Column names containing taxonomic information not identical to data_true")


# CALC ACCURACIES --------------------------------------------------------------

# TODO: calculate accuracy without calibrations 
nt_calibrations <- 
  get_calibration(nt_long, data_true, observed_everywhere) |> 
  mutate(rank = factor(rank, levels = ranks), 
         set = factor(set, levels = c("All", "Observed", "Novel")), 
         model = factor(model, levels = algs_sorted))

aa_calibrations <- 
  get_calibration(aa_long, data_true, observed_everywhere) |> 
  mutate(rank = factor(rank, levels = ranks), 
         set = factor(set, levels = c("All", "Observed", "Novel")), 
         model = factor(model, levels = algs_sorted))

### If a sequence belongs to a new taxon (on any rank), how well is it predicted on higher ranks? 
nt_accuracies <- get_accuracy_from_cal(nt_calibrations) |> 
  mutate(rank = factor(rank, levels = ranks), 
         set = factor(set, levels = c("All", "Observed", "Novel")))
aa_accuracies <- get_accuracy_from_cal(aa_calibrations) |> 
  mutate(rank = factor(rank, levels = ranks), 
         set = factor(set, levels = c("All", "Observed", "Novel")))

tab <- 
  bind_rows(nt_accuracies |> mutate(datatype = "nt"), 
          aa_accuracies |> mutate(datatype = "aa")) |> 
  mutate(cname = map2_chr(set, datatype, ~paste0(.y, "_", .x))) |> 
  select(-set, -datatype) |> 
  pivot_wider(names_from = cname, values_from = accuracy) |> 
  filter(!if_all(c(aa_All, aa_Observed, aa_Novel), is.na)) |> 
  select(model, rank, nt_All, aa_All, 
         nt_Observed, aa_Observed, 
         nt_Novel, aa_Novel)

print(xtable::print.xtable(xtable::xtable(tab), type = "latex", include.rownames = F))

