
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

# LOAD DATA --------------------------------------------------------------------

# Load test and train data 
data <- load_FinBOL_GBOL()

# Replace labels of taxa unique to test
# Sort test data in alphabetical order
data_true <- get_data_true(data$test, data$train) |>
  dplyr::arrange(ID)

# CORRECT COLUMN NAMES AND ORDER -----------------------------------------------

ranks <- colnames(data_true)[-1]
correct_cols <- c("ID", 
                  ranks, 
                  sapply(ranks, function(x) paste0("Prob_", x), USE.NAMES = F))

# READ RESULTS -----------------------------------------------------------------

# BayesANT
# result_BayesANT <-
#   read.table(paste0("models/results_bayesant_finbol_gbol/bayesant_test", 
#                     shorttxt, "_finbol-gbol_", level, ".txt"),
#              header = T) |>
#   dplyr::arrange(ID)

# RDP
# result_RDP <-
#   read.table(paste0("models/results_rdp_finbol_gbol/rdp_test", 
#                     shorttxt, "_finbol-gbol_", level, ".txt"),
#              header = T) |>
#   dplyr::arrange(ID)

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
  read.table(paste0("results/sintax/sintax_test", shorttxt, "_nt_4.tsv"),
                            header = TRUE)
result_SINTAX <- dplyr::arrange(result_SINTAX, ID)
result_SINTAX <- arrange_columns(result_SINTAX, correct_cols)

# EPA-ng phylogenetic tree 
result_epaphyl <- 
  read.table(paste0("results/epang_phyltree/epang_phyltree_test",
                  shorttxt, "_finbol-gbol.txt"),
           header = T, sep = "\t")  |>
  dplyr::arrange(ID)
result_epaphyl$species <- sapply(result_epaphyl$species, 
                                 function(x) stringr::str_replace(x, " ", "_"),
                                 USE.NAMES = F)
result_epaphyl <- rename_unk_output(result_epaphyl)
result_epaphyl <- arrange_columns(result_epaphyl, correct_cols)

# BLAST 
# TODO: what is the species column separated by ";"? 
# result_blast <- 
#   read.table(paste0("results/blast/blast_top_hit_test", shorttxt, "_nt_16.tsv"), 
#            header = T)

# DNA-barcoder
# result_dnabarcoder <- 
#   read.table(paste0("results/dnabarcoder/test", shorttxt, "_nt_4.tsv"), 
#            header = T)
# result_dnabarcoder <- arrange_columns(result_dnabarcoder, correct_cols)

# IDTAXA 
# TODO: Some are missing an identification 
#result_idtaxa <- 
#  read.table(paste0("results/idtaxa/idtaxa_test", shorttxt, "_nt_4.tsv"))

# MycoAI-BERT
# TODO: deal with missing class 
# result_aibert <- 
#   read.table(paste0("results/mycoai_bert/mycoai_bert_test", shorttxt, 
#                     "_nt_gpu.tsv"), 
#              header = T) |> 
#   dplyr::arrange(ID)
# result_aibert <- arrange_columns(result_aibert, correct_cols)

# MycoAI-CNN 
# TODO: deal with missing class 
# result_aicnn <-
#   read.table(paste0("results/mycoai_cnn/mycoai_cnn_test", shorttxt,
#                     "_nt_gpu.tsv"),
#              header = T) |>
#   dplyr::arrange(ID)
# result_aicnn <- arrange_columns(result_aicnn, correct_cols)

# Brendan's mystery model 
result_mystery <- read.table("models/mystery/mystery_model_finbol_gbol.txt",
                             header = T) |> 
  dplyr::arrange(ID)
result_mystery[is.na(result_mystery)] <- 0
result_mystery <- arrange_columns(result_mystery, correct_cols)


# PREPARE AND CHECK RESULT -----------------------------------------------------

results <- list(
#  "BayesANT" = result_BayesANT,
#  "PROTAX" = result_PROTAX,
  "EPA-ng taxtree" = result_epatax,
#  "RDP" = result_RDP,
  "SINTAX" = result_SINTAX,
  "Mystery" = result_mystery,
  "EPA-ng phyltree" = result_epaphyl
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


# PLOT CALIBRATION CURVES ------------------------------------------------------

calibrations <- get_calibration(results, data_true, observed_everywhere)
calibrations$rank <- factor(calibrations$rank, 
                            levels = colnames(data_true)[-1])
calibrations$set <- factor(calibrations$set, 
                           levels = c("All", "Observed", "Novel"))

#### Plot all calibration curves (all ranks, taxa sets and models) #####
p_cal <- plot_calibration(calibrations) 

ggplot2::ggsave(plot = p_cal, 
                filename = paste0("../plots/calibration_COI_all", undshort,
                                  "_", level, ".pdf"), 
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
                filename = paste0("../plots/accuracy_COI", undshort, 
                                  "_", level, ".pdf"), 
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
  filter(measure == "marginal", set != "All") |> 
  select(-denom_cond) |> 
  mutate(set = map_chr(set, ~str_to_sentence(paste(.x, "taxa"))), 
         accuracy = accuracy * 100) |> 
  bind_rows(accuracies |> 
              filter(set != "All") |> 
              mutate(set = paste(set, "species"), 
                     measure = "accuracy")) |> 
  mutate(set = factor(set, levels = c("Observed species",
                                      "Novel species",
                                      "Observed taxa",
                                      "Novel taxa"))) |> 
  ggplot2::ggplot() + 
  ggplot2::theme_bw() + 
  ggplot2::theme(aspect.ratio = 1, 
                 panel.grid.minor = element_blank(), 
                 axis.title.x = element_blank(), 
                 axis.text.x = element_text(angle = 45, hjust = 1)) + 
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
