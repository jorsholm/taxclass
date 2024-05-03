
rm(list = ls())

# LOAD FUNCTIONS ---------------------------------------------------------------

source("load_FinBOL_GBOL.R")
source("functions.R")

# LOAD DATA --------------------------------------------------------------------

level <- "species"
short <- TRUE
shorttxt <- ""
undshort <- ""
if(short){
  shorttxt <- "short"
  undshort <- "_short"
} 

# Load test and train data 
data <- load_FinBOL_GBOL(level = level)

# Replace labels of taxa unique to test
# Sort test data in alphabetical order
data_true <- get_data_true(data$test, data$train) |>
  dplyr::arrange(ID)

# READ RESULTS -----------------------------------------------------------------

# BayesANT
result_BayesANT <-
  read.table(paste0("models/results_bayesant_finbol_gbol/bayesant_test", 
                    shorttxt, "_finbol-gbol_", level, ".txt"),
             header = T) |>
  dplyr::arrange(ID)

# RDP
result_RDP <-
  read.table(paste0("models/results_rdp_finbol_gbol/rdp_test", 
                    shorttxt, "_finbol-gbol_", level, ".txt"),
             header = T) |>
  dplyr::arrange(ID)

# PROTAX
result_PROTAX <- 
  read.table(paste0("models/results_protax_finbol_gbol/protax_finbol_gbol_",
                    level, undshort, ".txt"),
             header = T)
result_PROTAX[is.na(result_PROTAX)] <- 0
result_PROTAX <- rename_PROTAX_output(result_PROTAX)
# Quick check that order of columns is identical to BayesANT
# all(stringr::str_to_lower(colnames(result_PROTAX)) ==
#       stringr::str_to_lower(colnames(result_BayesANT)))
# Rename Protax columns
colnames(result_PROTAX) <- colnames(result_BayesANT)

# EPA-ng Taxonomy tree
result_epatax <- 
  read.table(paste0("models/results_epa_finbol_gbol/epa_test", 
                    shorttxt, "_finbol-gbol_", level,".txt"),
             header = T) |>
  dplyr::arrange(ID)
result_epatax <- rename_PROTAX_output(result_epatax)
# all(stringr::str_to_lower(colnames(result_epatax)) ==
#       stringr::str_to_lower(colnames(result_BayesANT)))
colnames(result_epatax) <- colnames(result_BayesANT)

# Sintax
result_SINTAX <- 
  read.table(paste0("models/results_sintax_finbol_gbol/sintax_finbol_gbol_",
                    level, undshort, ".txt"),
                            header = TRUE)
result_SINTAX <- rename_PROTAX_output(result_SINTAX)
result_SINTAX <- dplyr::arrange(result_SINTAX, ID)
# all(stringr::str_to_lower(colnames(result_SINTAX)) ==
#       stringr::str_to_lower(colnames(result_BayesANT)))
colnames(result_SINTAX) <- colnames(result_BayesANT)

# PREPARE AND CHECK RESULT -----------------------------------------------------

results <- list(
  "BayesANT" = result_BayesANT,
  "PROTAX" = result_PROTAX,
  "EPA-ng taxtree" = result_epatax,
  "RDP" = result_RDP,
  "SINTAX" = result_SINTAX
)

#### Identify novel taxa across all ranks 

id_novel <- list()
for(rank in colnames(data_true)[-1]){
  id_novel[[rank]] <- which(str_ends(data_true[,rank], paste0(rank, "_new")))
}

id_observed <- list()
for(rank in colnames(data_true)[-1]){
  id_observed[[rank]] <- which(!str_ends(data_true[,rank], "_new"))
}

id_all <- lapply(data_true, function(x) 1:length(x))

# Some quick checks 
if(!length(unique(lapply(results, function(x) colnames(x)))) == 1) print("One of the results data frames have wrong column names.")
if(!all(sapply(results, function(x) all(x[,1] == data_true[,1])))) print("One of the results data frames is not sorted by ID.")
if(!all(colnames(data_true) == colnames(results[[1]])[1:ncol(data_true)])) print("Column names containing taxonomic information not identical to data_true")

# PLOT CALIBRATION CURVES ------------------------------------------------------

calibrations <- get_calibration(results, data_true, id_observed$Species)

calibrations$rank <- factor(calibrations$rank, 
                            levels = colnames(data_true)[-1])
calibrations$set <- factor(calibrations$set, 
                           levels = c("All", "Observed", "Novel"))

#### Plot all calibration curves (all ranks, taxa sets and models)
p_cal <- plot_calibration(calibrations) 

ggplot2::ggsave(plot = p_cal, 
                filename = paste0("../plots/calibration_COI", undshort, ".png"), 
                width = 210, 
                height = 297, 
                units = "mm")

# PLOT ACCURACIES --------------------------------------------------------------

### If a sequence belongs to a new taxon (on any rank), how well is it predicted on higher ranks? 
accuracies <- get_accuracy_from_cal(calibrations)

accuracies$rank <- factor(accuracies$rank,
                          levels = colnames(data_true)[-1])
accuracies$set <- factor(accuracies$set, 
                           levels = c("All", "Observed", "Novel"))

#### Plot accuracy across all models and taxa sets 
p_acc <- plot_accuracies(accuracies)

ggplot2::ggsave(plot = p_acc, 
                filename = paste0("../plots/accuracy_COI", undshort, ".png"), 
                width = 210, 
                height = 150, 
                units = "mm")



#### Marginal accuracy 
marg_accuracies_novel <- marginal_accuracy(results, data_true, id_novel)

marg_accuracies_obs <- marginal_accuracy(results, data_true, id_observed)

marg_accuracies_all <- marginal_accuracy(results, data_true, id_all)

### Conditional accuracy 
cond_accuracies_novel <- conditional_accuracy(results, data_true, id_novel)

cond_accuracies_obs <- conditional_accuracy(results, data_true, id_observed)

cond_accuracies_all <- conditional_accuracy(results, data_true, id_all)

accuracy_df <- rbind(cbind(marg_accuracies_all, "set" = "all"),
                     cbind(marg_accuracies_obs, "set" = "observed"),
                     cbind(marg_accuracies_novel, "set" = "novel")
                     ) |>
  dplyr::left_join(rbind(cbind(cond_accuracies_all, "set" = "all"), 
                         cbind(cond_accuracies_obs, "set" = "observed"), 
                         cbind(cond_accuracies_novel, "set" = "novel")
                         )) |>
  dplyr::rename(marginal = marg_accuracy, 
                conditional = cond_accuracy) |> 
  tidyr::pivot_longer(c(3,5), 
                      names_to = "measure", 
                      values_to = "accuracy")

accuracy_df$rank <- factor(accuracy_df$rank, 
                           levels = colnames(data_true)[-1])
accuracy_df$set <- factor(accuracy_df$set, 
                          levels = c("all", "observed", "novel"))
accuracy_df$measure <- factor(accuracy_df$measure, 
                              levels = c("marginal", "conditional"))

## TODO : Put numbers in cond accuracy denominator
cond_marg_plot <- 
  ggplot2::ggplot() + 
  ggplot2::theme_bw() + 
  ggplot2::theme(aspect.ratio = 1) + 
  ggplot2::geom_line(data = accuracy_df, 
                     mapping = ggplot2::aes(x = rank, 
                                            y = accuracy, 
                                            color = model, 
                                            group = model)) + 
  ggplot2::geom_point(data = accuracy_df, 
                     mapping = ggplot2::aes(x = rank, 
                                            y = accuracy, 
                                            color = model)) + 
  ggplot2::facet_grid(measure~set)

ggplot2::ggsave(plot = cond_marg_plot, 
                filename = paste0("../plots/marg_cond_accuracy", undshort, ".pdf"), 
                width = 310, 
                height = 200, 
                units = "mm")

# How many of the novel taxa are predicted as novel *on any taxonomic rank* (in practice the correct one or above)
a <- get_novelty_accuracy(results, id_novel)

novel_accuracy <- a$acc
novel_count <- a$count
partitioned_novelty <- a$part_novelty

novel_accuracy <- 
  novel_accuracy |> 
  dplyr::mutate(rank = purrr::map_chr(rank, ~paste0(.x, " (", novel_count[[.x]], ")")))

novel_accuracy$rank <- factor(novel_accuracy$rank, 
                              levels = unique(novel_accuracy$rank))

p <- ggplot2::ggplot() +
  ggplot2::theme_bw() +
  ggplot2::theme(aspect.ratio = 1) + 
  ggplot2::ylim(0, 100) +
  ggplot2::theme(axis.title.x = element_blank()) + 
  ggplot2::ylab("Accuracy %")

# Proportion novel taxa correctly predicted novel (on any taxonomic level, i.e. correct or above)
# truly_novel <- p + 
#   ggplot2::geom_line(data = novel_accuracy, 
#                      mapping = ggplot2::aes(x = rank, 
#                                             y = novel_accuracy, 
#                                             color = model, 
#                                             group = model)) + 
#   ggplot2::geom_point(data = novel_accuracy, 
#                       mapping = ggplot2::aes(x = rank,
#                                              y = novel_accuracy,
#                                              color = model)) + 
#   ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   ggplot2::labs(color = "Model") + 
#   ggplot2::facet_wrap(~paste("Correctly predicted novel (on any rank) /\n the number of novel taxa on that rank"))



# How many were predicted new (on any level) that were actually observed?
# TODO: scale this to the number of true novel sequences 
# wrong_any <- p + 
#   ggplot2::geom_line(data = novel_accuracy, 
#                      mapping = ggplot2::aes(x = rank, 
#                                             y = false_novelty, 
#                                             color = model, 
#                                             group = model)) + 
#   ggplot2::geom_point(data = novel_accuracy, 
#                       mapping = ggplot2::aes(x = rank,
#                                              y = false_novelty,
#                                              color = model)) + 
#   ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   ggplot2::labs(color = "Model") + 
#   ggplot2::facet_wrap(~paste("Wrongly predicted as novel (on any rank) /\n all seq with known id on that rank"))

# Partition the correctly predicted novel taxa 
partitioned_novelty <- 
  partitioned_novelty |> 
  dplyr::filter(novelty_rank != "None") |> 
  dplyr::mutate(novelty_rank = if_else(novelty_rank == rank, 
                                       "Correct rank", novelty_rank))

partitioned_novelty$novelty_rank <- factor(partitioned_novelty$novelty_rank, 
                                           levels = c(colnames(data_true)[-1], "Correct rank"))
partitioned_novelty$rank <- factor(partitioned_novelty$rank, 
                                   levels = colnames(data_true)[-1])

color_alternatives <- c("#eff3ff", "#bdd7e7", "#6baed6", "#3182bd", "#08519c")

part_novelty <- p + 
  ggplot2::geom_bar(data = partitioned_novelty, 
                    mapping = ggplot2::aes(x = rank, 
                                           y = perc, 
                                           fill = novelty_rank), 
                    position = "stack", 
                    stat = "identity") + 
  ggplot2::facet_wrap(~model) + 
  ggplot2::scale_fill_manual(limits = c(colnames(data_true)[-1], "Correct rank"), 
                             values = c(color_alternatives[1:length(colnames(data_true)[-1])], "grey"))

ggpubr::ggarrange(plotlist = list(accuracy_plot, part_novelty), 
                  nrow = 2)

