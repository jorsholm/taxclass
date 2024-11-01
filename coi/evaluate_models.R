
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
keep_na <- T
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

# CORRECT COLUMN NAMES AND ORDER -----------------------------------------------

ranks <- colnames(data_true)[-1]
correct_cols <- c("ID",
                  ranks,
                  sapply(ranks, function(x) paste0("Prob_", x), USE.NAMES = F))

# READ RESULTS -----------------------------------------------------------------

results <- readRDS(paste0("result_list", natxt, ".rds"))

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

# data_true has a subset of sequences 
results <- lapply(results, function(x) x[which(x[,1] %in% data_true[,1]),])


# Some quick checks 
if(!length(unique(lapply(results, function(x) colnames(x)))) == 1) print("One of the results data frames have wrong column names.")
if(!all(sapply(results, function(x) all(x[,1] == data_true[,1])))) print("One of the results data frames is not sorted by ID.")
if(!all(colnames(data_true) == colnames(results[[1]])[1:ncol(data_true)])) print("Column names containing taxonomic information not identical to data_true")


# PLOT CALIBRATION CURVES ------------------------------------------------------

calibrations <- get_calibration(results, data_true, observed_everywhere)
calibrations$rank <- factor(calibrations$rank, 
                            levels = colnames(data_true)[-1])
calibrations$set <- factor(calibrations$set, 
                           levels = c("All", "Observed", "Novel"))

#### Plot all calibration curves (all ranks, taxa sets and models) #####
p_cal <- plot_calibration(calibrations) 

ggplot2::ggsave(plot = p_cal, 
                filename = paste0("../plots/calibration_COI_all", undshort, ".pdf"), 
                width = 210, 
                height = 297, 
                units = "mm")

#### Subset plot for publication: observed and novel sets only #### 
obs_perc <- round(sum(observed_everywhere)/nrow(data_true) * 100, 1)
novel_perc <- round(sum(!observed_everywhere)/nrow(data_true) * 100, 1)

p_cal_subset <- 
  plot_calibration(calibrations |> 
                     dplyr::filter(set != "All")) +
  ggplot2::facet_grid(set~model, 
                      labeller = ggplot2::labeller(set = c(Observed = paste0("Observed species (", obs_perc, "%)"),
                                                           Novel = paste0("Novel species (", novel_perc, "%)")))) + 
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 legend.position = "bottom") + 
  ggplot2::guides(color = ggplot2::guide_legend(nrow = 1))

ggplot2::ggsave(plot = p_cal_subset, 
                filename = paste0("../plots/calibration_COI_subset", undshort, ".pdf"), 
                width = 235, 
                height = 130, 
                units = "mm")

# PLOT BINNED CALIBRATIONS -----------------------------------------------------
binned_calibrations <- get_calibration(results, data_true, observed_everywhere, 
                                       binned = T) 

binned_calibrations$rank <- factor(binned_calibrations$rank, 
                                   levels = colnames(data_true)[-1])
binned_calibrations$set <- factor(binned_calibrations$set, 
                                   levels = c("All", "Observed", "Novel"))
binned_calibrations <- 
  binned_calibrations |> 
  dplyr::mutate(bin = purrr::map_dbl(bin,
                                     function(x)
                                       stringr::str_remove(x, "\\(") |>
                                       stringr::str_remove("\\[") |>
                                       stringr::str_split(",") |>
                                       unlist() |>
                                       head(1) |>
                                       as.numeric())) |> 
  dplyr::mutate(bin = bin + 0.05)

#### Genus-level plot #### 
p_cal_binned_genus <- 
  binned_calibrations |> 
  dplyr::filter(set != "All") |> 
  dplyr::filter(rank == "Genus") |> 
  ggplot2::ggplot() + 
  ggplot2::ylim(0, 1) +
  ggplot2::xlim(0, 1.05) +
  ggplot2::theme_bw() +
  ggplot2::theme(aspect.ratio = 1) +
  ggplot2::geom_abline(intercept = 0, slope = 1, color = "grey") +
  ggplot2::geom_bar(ggplot2::aes(x = bin, 
                                 y = correct, 
                                 fill = count), 
                    stat = "identity", 
                    alpha = 0.5) + 
  ggplot2::facet_grid(set~model, 
                      labeller = ggplot2::labeller(set = 
                                                     c(Observed = paste0("Observed species (", obs_perc, "%)"),
                                                       Novel = paste0("Novel species (", novel_perc, "%)")))) + 
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) + 
  ggplot2::ylab("Prop. correct") + 
  ggplot2::xlab("Pred. probability") +
  ggplot2::scale_fill_gradient(name = "# predictions",
                               trans = "log",
                               breaks = c(10, 100, 1000),
                               labels = c(10, 100, 1000))

ggsave(plot = binned_genus, 
       filename = paste0("../plots/binned_calibration_genus_COI", undshort, ".pdf"), 
       width = 200, 
       height = 90, 
       units = "mm")

#### All ranks for all models #### 
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

pdf(paste0("../plots/binned_calibration_COI", undshort, ".pdf"), 
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
                filename = paste0("../plots/accuracy_COI", undshort, ".pdf"), 
                width = 210, 
                height = 150, 
                units = "mm")

# MARGINAL AND CONDITIONAL ACCURACY --------------------------------------------

marg_accuracies <- marginal_accuracy(results, data_true, 
                                     id_list = list(All = id_all, 
                                                    Observed = id_observed, 
                                                    Novel = id_novel))

cond_accuracies <- conditional_accuracy(results, data_true,
                                        id_list = list(All = id_all,
                                                       Observed = id_observed,
                                                       Novel = id_novel))

accuracy_df <- 
  marg_accuracies |> 
  dplyr::left_join(cond_accuracies) |>
  dplyr::rename(marginal = marg_accuracy, 
                conditional = cond_accuracy) |> 
  tidyr::pivot_longer(c("marginal", "conditional"), 
                      names_to = "measure", 
                      values_to = "accuracy") |> 
  dplyr::mutate(denom_cond = dplyr::if_else(measure == "marginal", 0, denom_cond)) |> 
  dplyr::mutate(denom_cond = dplyr::na_if(denom_cond, 0))

accuracy_df$rank <- factor(accuracy_df$rank, 
                           levels = colnames(data_true)[-1])
accuracy_df$set <- factor(accuracy_df$set, 
                          levels = c("All", "Observed", "Novel"))
accuracy_df$measure <- factor(accuracy_df$measure, 
                              levels = c("marginal", "conditional"))

#### Plot marginal and conditional accuracy #### 

fourpanel_accuracy <-
  accuracy_df |> 
  dplyr::filter(measure == "marginal", set != "All") |> 
  dplyr::select(-denom_cond) |> 
  dplyr::mutate(set = purrr::map_chr(set, ~stringr::str_to_sentence(paste(.x, "taxa"))),
                accuracy = accuracy * 100) |> 
  dplyr::bind_rows(accuracies |>
                     dplyr::filter(set != "All") |>
                     dplyr::mutate(set = paste(set, "species"),
                                   measure = "accuracy")) |> 
  dplyr::mutate(set = factor(set, levels = c("Observed species",
                                      "Novel species",
                                      "Observed taxa",
                                      "Novel taxa"))) |> 
  ggplot2::ggplot() + 
  ggplot2::theme_bw() + 
  ggplot2::theme(aspect.ratio = 1, 
                 panel.grid.minor = ggplot2::element_blank(), 
                 axis.title.x = ggplot2::element_blank(), 
                 axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) + 
  ggplot2::aes(x = rank, 
               y = accuracy, 
               color = model) + 
  ggplot2::geom_line(mapping = ggplot2::aes(group = model)) + 
  ggplot2::geom_point() + 
  ggplot2::facet_wrap(~set) + 
  ggplot2::labs(color = "Model", y = "Accuracy (%)") 

ggplot2::ggsave(filename = paste0("../plots/fourpanel_accuracy_COI", 
                                  undshort, ".pdf"), 
                plot = fourpanel_accuracy, 
                height = 150, 
                width = 200, 
                units = "mm")

#### Single panel conditional accuracy ####
# TODO: Fix this for new models 
cond_accuracy_novel_taxa <- 
  accuracy_df |> 
  filter(measure == "conditional" & set == "Novel") |> 
  # filter(model != "RDP" & model != "SINTAX") |> 
  ggplot() + 
  geom_line(aes(x = rank, y = accuracy * 100, 
                group = model, color = model)) + 
  geom_point(aes(x = rank, y = accuracy * 100, 
                 group = model, color = model)) + 
  geom_text(aes(x = rank, y = accuracy * 100, 
                color = model, group = model, 
                label = denom_cond), 
            nudge_y = 5, size = 3) + 
  # scale_color_manual(limits = c("BayesANT", "EPA-ng taxtree", "PROTAX"), 
  #                    values = c("#f8766d", "#a3a500", "#00bf7d")) + 
  theme_bw() + 
  ylab("Conditional recall (%)") + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.minor.y = element_blank()) + 
  labs(color = "Model")

# ggsave(plot = cond_accuracy_novel_taxa, 
#        filename = "../plots/publication/cond_recall.pdf", 
#        width = 130, 
#        height = 90, 
#        units = "mm")

# NOVELTY: TOTAL NUMBER OF PREDICTION ------------------------------------------

pred_novel <- count_novel(results, data_true, id_novel)
pred_novel$rank <- factor(pred_novel$rank, levels = names(data_true)[-1])

p_tot_novel <-
  pred_novel |> 
  ggplot2::ggplot() + 
  ggplot2::geom_line(data = pred_novel |> filter(model == model[1]),
                     ggplot2::aes(x = rank, y = sum_novel, group = model),
                     color = "black",
                     linetype = "dashed") + 
  ggplot2::aes(x = rank, y = pred_novel,
               color = model, group = model) +
  ggplot2::geom_line() + 
  ggplot2::geom_point() + 
  ggplot2::theme_bw() + 
  ggplot2::scale_y_log10()

ggsave(plot = p_tot_novel, 
       filename = paste0("../plots/sum_novel_COI", undshort, ".pdf"), 
       height = 150, 
       width = 200, 
       units = "mm")

# NOVELTY: PARTITIONED ---------------------------------------------------------

partitioned_novelty <- 
  get_partitioned_novelty(data_true, results, observed_everywhere) |> 
  dplyr::filter(pred_rank != "None")

partitioned_novelty <- 
  partitioned_novelty |> 
  dplyr::mutate(pred_lower = purrr::map2_lgl(true_rank, pred_rank,
                                             ~which(ranks == .x) < which(ranks == .y))) |> 
  dplyr::mutate(perc = purrr::map2_dbl(pred_lower, perc,
                                       ~dplyr::if_else(.x, -.y, .y))) |> 
  dplyr::mutate(pred_rank = dplyr::if_else(pred_rank == true_rank, 
                                           "Correct", pred_rank))
  
partitioned_novelty$true_rank <- factor(partitioned_novelty$true_rank,
                                        levels = colnames(data_true)[-1])
partitioned_novelty$pred_rank <- factor(partitioned_novelty$pred_rank,
                                        levels = c(colnames(data_true)[-1], 
                                                   "Correct"))

plot_colors <- c("#EFF3FF", "#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C")

p_part_novelty <- 
  ggplot2::ggplot() + 
  ggplot2::geom_bar(data = partitioned_novelty, 
                    ggplot2::aes(x = true_rank, 
                                 y = perc, 
                                 fill = pred_rank), 
                    position = "stack", 
                    stat = "identity") + 
  ggplot2::facet_wrap(~model) + 
  ggplot2::scale_fill_manual(limits = c(colnames(data_true)[-1],
                                        "Correct"), 
                             values = c(color_alternatives[1:length(colnames(data_true)[-1])],
                                        "grey")) + 
  ggplot2::labs(fill = "Predicted novel rank",
                y= "Proportion of predictions (%)") + 
  ggplot2::theme_bw() + 
  ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                 axis.title.x = element_blank()) + 
  ggplot2::geom_hline(yintercept = 0, 
                      color = "grey20", linetype = "dashed")

ggsave(plot = p_part_novelty, 
       filename = paste0("../plots/part_novelty_COI", undshort, ".pdf"), 
       width = 220, 
       height = 130, 
       units = "mm")

# % CLASSIFIED ~ % CORRECT -----------------------------------------------------


# TODO: BLAST can be included with % identity 

threshold_curve <- function(results, data_true, thresholds = seq(0, 1, 0.01)){
  
  out <- data.frame()
  ranks <- colnames(data_true)[-1]
  n <- nrow(data_true)
  
  for(i in 1:length(results)){
    
    for(r in ranks){
      
      probcol <- paste0("Prob_", r)
      
      correct <- sapply(thresholds, function(x)
        sum(results[[i]][which(results[[i]][, probcol] >= x), r] == data_true[which(results[[i]][, probcol] >= x), r]) /
          length(which(results[[i]][, probcol] >= x)) * 100)
      classified <- sapply(thresholds, 
                           function(x) length(which(results[[i]][,probcol] >= x))/n * 100)
      
      out <- rbind(out, 
            data.frame(model = names(results)[i], 
                       rank = r, 
                       threshold = thresholds, 
                       correct = correct, 
                       classified = classified))
      
    }
  }
  return(out)
}

test <- threshold_curve(results, data_true)

test |> 
  ggplot2::ggplot() + 
  ggplot2::geom_line(ggplot2::aes(x = classified, 
                                  y = correct, 
                                  color = model)) + 
  ggplot2::facet_wrap(~rank, 
                      scales = "free_x") + 
  ggplot2::theme_bw() + 
  ggplot2::theme(aspect.ratio = 1) + 
  ggplot2::labs(x = "% classified", 
                y = "% correct", 
                color = "Model")


