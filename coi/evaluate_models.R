
rm(list = ls())
setwd("coi")

# LOAD FUNCTIONS ---------------------------------------------------------------

source("load_FinBOL_GBOL.R")
source("functions.R")

# LOAD DATA --------------------------------------------------------------------

level <- "species"
short <- F
shorttxt <- ""
undshort <- ""
if(short){
  shorttxt <- "short"
  undshort <- "_short"
} 

# Load test and train data 
data <- load_FinBOL_GBOL(level = level, short = short)

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

# Some quick checks 
if(!length(unique(lapply(results, function(x) colnames(x)))) == 1) print("One of the results data frames have wrong column names.")
if(!all(sapply(results, function(x) all(x[,1] == data_true[,1])))) print("One of the results data frames is not sorted by ID.")
if(!all(colnames(data_true) == colnames(results[[1]])[1:ncol(data_true)])) print("Column names containing taxonomic information not identical to data_true")

# PLOT CALIBRATION -------------------------------------------------------------

calibrations <- get_calibration(results, data_true, observed_everywhere)
calibrations$rank <- factor(calibrations$rank, 
                            levels = colnames(data_true)[-1])
calibrations$set <- factor(calibrations$set, 
                           levels = c("All", "Observed", "Novel"))

#### Plot all calibration curves (all ranks, taxa sets and models) #####
p_cal <- plot_calibration(calibrations) 

ggplot2::ggsave(plot = p_cal, 
                filename = paste0("../plots/calibration_COI", undshort,
                                  "_", level, ".pdf"), 
                width = 210, 
                height = 297, 
                units = "mm")

#### Calibration plot for publication #### 
obs_perc <- round(sum(observed_everywhere)/nrow(data_true) * 100, 1)
novel_perc <- round(sum(!observed_everywhere)/nrow(data_true) * 100, 1)

calibrations_for_pub <- 
  calibrations |> 
  dplyr::filter(set != "All")

calibrations_for_pub$set <- as.character(calibrations_for_pub$set)
calibrations_for_pub$set[which(calibrations_for_pub$set == "Observed")] <- paste0("Observed (", obs_perc, "%)")
calibrations_for_pub$set[which(calibrations_for_pub$set == "Novel")] <- paste0("Novel (", novel_perc, "%)")
calibrations_for_pub$set <- factor(calibrations_for_pub$set, levels = c(paste0("Observed (", obs_perc, "%)"), 
                                                                        paste0("Novel (", novel_perc, "%)")))

calibration_plot_publ <- 
  plot_calibration(calibrations_for_pub) +
  ggplot2::facet_grid(set~model) + 
  theme(panel.grid.minor = element_blank(), 
        legend.position = "bottom")

ggplot2::ggsave(plot = calibration_plot_publ, 
                filename = paste0("../plots/publication/calibration", undshort, ".pdf"), 
                width = 235, 
                height = 130, 
                units = "mm")
ggplot2::ggsave(plot = calibration_plot_publ, 
                filename = paste0("../plots/publication/calibration", undshort, ".png"), 
                width = 235, 
                height = 130, 
                units = "mm")

#### Binned calibration #####
binned_calibrations <- get_binned_calibration(results, data_true, observed_everywhere)

binned_calibrations$rank <- factor(binned_calibrations$rank, 
                                   levels = colnames(data_true)[-1])
binned_calibrations$set <- factor(binned_calibrations$set, 
                                   levels = c("All", "Observed", "Novel"))

binned_calibrations <- 
  binned_calibrations |> 
  dplyr::mutate(bin = purrr::map_dbl(bin, function(x) 
                                     stringr::str_remove(x, "\\(") |> 
                                       stringr::str_split(",") |> 
                                       unlist() |> 
                                       head(1) |> 
                                       as.numeric())) |> 
  dplyr::mutate(bin = bin + 0.05)

# Genus-level plot 
binned_genus <- 
  binned_calibrations |> 
  dplyr::filter(set != "All") |> 
  dplyr::filter(rank == "Genus") |> 
  ggplot() + 
  ggplot2::ylim(0, 1) +
  ggplot2::xlim(0, 1.05) +
  ggplot2::theme_bw() +
  ggplot2::theme(aspect.ratio = 1) +
  ggplot2::geom_abline(intercept = 0, slope = 1, color = "grey") +
  ggplot2::geom_bar(ggplot2::aes(x = bin, 
                                 y = correct, 
                                 fill = log(count)), 
                    stat = "identity", 
                    alpha = 0.5) + 
  ggplot2::facet_grid(set~model) + 
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) + 
  ggplot2::ylab("Prop. correct") + 
  ggplot2::xlab("Pred. probability") 

ggsave(plot = binned_genus, 
       filename = "../plots/binned_genus.pdf", 
       width = 200, 
       height = 100, 
       units = "mm")

# Plot all ranks for all models 
plotlist <- list()

for(m in names(results)){
  
  plotlist[[m]] <- 
    ggplot2::ggplot() + 
    ggplot2::ylim(0, 1) +
    ggplot2::xlim(0, 1.05) +
    ggplot2::theme_bw() +
    ggplot2::theme(aspect.ratio = 1) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "grey") +
    ggplot2::geom_bar(data = binned_calibrations |> dplyr::filter(model == m), 
                      ggplot2::aes(x = bin, 
                                   y = correct, 
                                   fill = log(count)), 
                      stat = "identity", 
                      alpha = 0.5) + 
    ggplot2::facet_grid(set~rank) + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) + 
    ggplot2::ylab("Prop. correct") + 
    ggplot2::xlab("Pred. probability") + 
    ggplot2::ggtitle(m)
}

pdf(paste0("../plots/binned_calibration", undshort, "_", level, ".pdf"), 
    width = 10)
lapply(plotlist, print)
dev.off()

# PLOT ACCURACIES --------------------------------------------------------------

### If a sequence belongs to a new taxon (on any rank), how well is it predicted on higher ranks? 
accuracies <- get_accuracy_from_cal(calibrations)

accuracies$rank <- factor(accuracies$rank,
                          levels = colnames(data_true)[-1])
accuracies$set <- factor(accuracies$set, 
                           levels = c("All", "Observed", "Novel"))

#### Plot across all models and species sets ####
p_acc <- plot_accuracies(accuracies)

ggplot2::ggsave(plot = p_acc, 
                filename = paste0("../plots/accuracy_COI", undshort, 
                                  "_", level, ".pdf"), 
                width = 210, 
                height = 150, 
                units = "mm")

#### Marginal and conditional accuracy #####
marg_accuracies_novel <- marginal_accuracy(results, data_true, id_novel)

marg_accuracies_obs <- marginal_accuracy(results, data_true, id_observed)

marg_accuracies_all <- marginal_accuracy(results, data_true, id_all)

cond_accuracies_novel <- conditional_accuracy(results, data_true, id_novel)

cond_accuracies_obs <- conditional_accuracy(results, data_true, id_observed)

cond_accuracies_all <- conditional_accuracy(results, data_true, id_all)

### Plot marginal and conditional accuracy 
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
                      values_to = "accuracy") |> 
  dplyr::mutate(denom_cond = dplyr::if_else(measure == "marginal", 0, denom_cond)) |> 
  dplyr::mutate(denom_cond = dplyr::na_if(denom_cond, 0))

accuracy_df$rank <- factor(accuracy_df$rank, 
                           levels = colnames(data_true)[-1])
accuracy_df$set <- factor(accuracy_df$set, 
                          levels = c("all", "observed", "novel"))
accuracy_df$measure <- factor(accuracy_df$measure, 
                              levels = c("marginal", "conditional"))

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
  ggplot2::facet_grid(measure~set) + 
  ggplot2::geom_text(data = accuracy_df |> dplyr::filter(measure == "conditional", 
                                                         set == "novel", 
                                                         !(model %in% c("RDP", "SINTAX"))), 
                     mapping = ggplot2::aes(x = rank, y = accuracy, color = model, 
                                            label = denom_cond), 
                     nudge_y = 0.05, size = 3)

ggplot2::ggsave(plot = cond_marg_plot, 
                filename = paste0("../plots/marg_cond_accuracy", undshort, 
                                  "_", level, ".pdf"), 
                width = 310, 
                height = 200, 
                units = "mm")

##### For publ 
accuracy_for_publ <- 
  accuracy_df |>
  dplyr::filter(measure == "marginal") |> 
  dplyr::mutate(set = stringr::str_to_title(set))

accuracy_for_publ$set <- factor(accuracy_for_publ$set, 
                                levels = c("All", "Observed", "Novel"))

marginal_accuracy_publ <- 
  ggplot2::ggplot() + 
  ggplot2::theme_bw() + 
  ggplot2::theme(aspect.ratio = 1) + 
  ggplot2::geom_line(data = accuracy_for_publ, 
                     mapping = ggplot2::aes(x = rank, 
                                            y = accuracy * 100, 
                                            color = model, 
                                            group = model)) + 
  ggplot2::geom_point(data = accuracy_for_publ, 
                      mapping = ggplot2::aes(x = rank, 
                                             y = accuracy * 100, 
                                             color = model)) + 
  ggplot2::facet_wrap(~set) + 
  ggplot2::ylab("Accuracy (%)") + 
  ggplot2::labs(color = "Model") + 
  theme(panel.grid.minor = element_blank(), 
        axis.title.x = element_blank())

ggplot2::ggsave(plot = marginal_accuracy_publ, 
                filename = paste0("../plots/publication/marginal_accuracy", 
                                  undshort, ".pdf"), 
                width = 280, 
                height = 100, 
                units = "mm")
ggplot2::ggsave(plot = marginal_accuracy_publ, 
                filename = paste0("../plots/publication/marginal_accuracy", 
                                  undshort, ".png"), 
                width = 280, 
                height = 100, 
                units = "mm")
### FOUR PANEL ACCURACY PLOT 

fourpanel_df <- 
  bind_rows(
  accuracy_df |>
    filter(measure == "marginal", set != "all") |>
    select(-denom_cond) |>
    mutate(set = as.character(set)) |>
    mutate(set = map_chr(set, ~ str_to_sentence(paste(.x, "taxa")))) |>
    mutate(accuracy = accuracy * 100),
  accuracies |> 
    filter(set != "All") |>
    mutate(set = as.character(set)) |>
    mutate(set = paste(set, "species")) |>
    mutate(measure = "accuracy"))

fourpanel_df$set <- factor(fourpanel_df$set, levels = c("Observed species",
                                                        "Novel species", 
                                                        "Observed taxa",
                                                        "Novel taxa"))
ggplot2::ggplot() + 
  ggplot2::theme_bw() + 
  ggplot2::theme(aspect.ratio = 1) + 
  ggplot2::geom_line(data = fourpanel_df, 
                     mapping = ggplot2::aes(x = rank, 
                                            y = accuracy, 
                                            color = model, 
                                            group = model)) + 
  ggplot2::geom_point(data = fourpanel_df, 
                      mapping = ggplot2::aes(x = rank, 
                                             y = accuracy, 
                                             color = model)) + 
  ggplot2::facet_wrap(~set) + 
  ggplot2::ylab("Accuracy (%)") + 
  ggplot2::labs(color = "Model") + 
  theme(panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1))


### True vs false novelty 
predicted_novelty <- function(results, data_true, id_novel){
  ranks <- colnames(data_true)[-1]
  pred_df <- data.frame()
  
  for(i in 1:length(results)){
    
    sum_novel <- 0 
    cumid_novel <- c()
    pred_novel <- 0 
    
    for(r in ranks){
      
      sum_novel <- sum_novel + length(id_novel[[r]])
      cumid_novel <- c(cumid_novel, id_novel[[r]])
      
      pred_novel <- pred_novel + 
        sum(stringr::str_ends(results[[i]][,r], "_new"))
      
      df <- data.frame(model = names(results)[i], 
                       rank = r,
                       sum_novel = sum_novel, 
                       pred_novel = pred_novel
      )
      
      pred_df <- rbind(pred_df, df)
    }
    
  }

  return(pred_df)
}

pred_novel <- predicted_novelty(results, data_true, id_novel)

pred_novel$rank <- factor(pred_novel$rank, levels = names(data_true)[-1])

#sum_novel <- 
  ggplot() + 
  geom_line(data = pred_novel |> filter(model == "BayesANT"), 
            aes(x = rank, y = sum_novel, group = model), 
            color = "black", 
            linetype = "dashed") + 
  geom_line(data = pred_novel, 
            aes(x = rank, y = pred_novel, color = model, group = model)) + 
  geom_point(data = pred_novel, 
            aes(x = rank, y = pred_novel, color = model, group = model)) + 
  theme_bw() + 
  scale_y_log10()

ggsave(plot = sum_novel, 
       filename = "../plots/sum_novel.pdf", 
       height = 150, 
       width = 200, 
       units = "mm")

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
  ggplot2::theme(axis.title.x = ggplot2::element_blank()) + 
  ggplot2::ylab("Accuracy %")

# Proportion novel taxa correctly predicted novel (on any taxonomic level, i.e. correct or above)
true_novel <- p +
  ggplot2::geom_line(data = novel_accuracy,
                     mapping = ggplot2::aes(x = rank,
                                            y = novel_accuracy,
                                            color = model,
                                            group = model)) +
  ggplot2::geom_point(data = novel_accuracy,
                      mapping = ggplot2::aes(x = rank,
                                             y = novel_accuracy,
                                             color = model)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::labs(color = "Model") +
  ggplot2::facet_wrap(~paste("True novel"))

# How many were predicted new (on any level) that were actually observed?
false_novel <- p +
  ggplot2::geom_line(data = novel_accuracy,
                     mapping = ggplot2::aes(x = rank,
                                            y = false_novelty,
                                            color = model,
                                            group = model)) +
  ggplot2::geom_point(data = novel_accuracy,
                      mapping = ggplot2::aes(x = rank,
                                             y = false_novelty,
                                             color = model)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::labs(color = "Model") +
  ggforce::facet_zoom(ylim = c(0,500))  

# Partition the correctly predicted novel taxa 
partitioned_novelty <- 
   partitioned_novelty |> 
   dplyr::filter(novelty_rank != "None") #|> 
#   dplyr::mutate(novelty_rank = dplyr::if_else(novelty_rank == rank, 
#                                        "Correct rank", novelty_rank))

# partitioned_novelty$novelty_rank <- factor(partitioned_novelty$novelty_rank, 
#                                            levels = c(colnames(data_true)[-1], "Correct rank"))
partitioned_novelty$rank <- factor(partitioned_novelty$rank, 
                                   levels = colnames(data_true)[-1])

partitioned_novelty$novelty_rank <- factor(partitioned_novelty$novelty_rank,
                                           levels = colnames(data_true)[-1])

color_alternatives <- c("#EFF3FF", "#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C")


#part_novelty <- 
  p + 
  ggplot2::geom_bar(data = partitioned_novelty, 
                    mapping = ggplot2::aes(x = rank, 
                                           y = perc, 
                                           fill = novelty_rank), 
                    position = "stack", 
                    stat = "identity") + 
  ggplot2::facet_wrap(~model) + 
  ggplot2::scale_fill_manual(limits = c(colnames(data_true)[-1], "Correct rank"), 
                             values = c(color_alternatives[1:length(colnames(data_true)[-1])], "grey"))

  
part_w_neg <- 
  partitioned_novelty  |> 
  mutate(rank_no = map_dbl(rank, ~which(ranks == .x)), 
         novelty_rank_no = map_dbl(novelty_rank, ~which(ranks == .x))) |> 
  mutate(pred_lower = (rank_no-novelty_rank_no)<0) |> 
  mutate(perc = map2_dbl(pred_lower, perc, 
                         ~if_else(.x, -.y, .y))) |> 
  mutate(novelty_rank = if_else(novelty_rank==rank, "Correct", novelty_rank))
  
part_w_neg$novelty_rank <- factor(part_w_neg$novelty_rank,
                                           levels = c(colnames(data_true)[-1], "Correct"))

part_w_neg_plot <- 
  p + 
  geom_bar(data = part_w_neg, 
           mapping = aes(x = rank, 
                         y = perc, 
                         fill = novelty_rank), 
           position = "stack", 
           stat = "identity") + 
  facet_wrap(~model) + 
  scale_fill_manual(limits = c(colnames(data_true)[-1], "Correct"), 
                    values = c(color_alternatives[1:length(colnames(data_true)[-1])], "grey")) + 
  ylab("Proportion of predictions (%)") + 
  labs(fill = "Predicted novel rank") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot = part_w_neg_plot, 
       filename = "../plots/part_novelty.pdf", 
       width = 220, 
       height = 130, 
       units = "mm")

novel_accuracy_plot <-
  ggpubr::ggarrange(plotlist = list(ggpubr::ggarrange(plotlist = list(true_novel, false_novel), 
                                                    ncol = 2, widths = c(1,2), 
                                                    common.legend = T, 
                                                    labels = c("", "False novel")),
                                  part_novelty), 
                  nrow = 2)

ggplot2::ggsave(filename = paste0("../plots/novel_accuracy", undshort, "_", level, ".pdf"), 
                plot = novel_accuracy_plot, 
                width = 13, 
                height = 9.5)

