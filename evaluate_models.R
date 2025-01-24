
rm(list = ls())

# LOAD FUNCTIONS AND PACKAGES --------------------------------------------------

library(tidyverse)
source("functions.R")

# case <- "its" 
case <- "coi"

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

if(case == "its"){
  ranks <- c("Kingdom",
             "Phylum",
             "Class", 
             "Order", 
             "Family", 
             "Genus", 
             "Species")
}else{
  ranks <- c("Class", 
             "Order",
             "Family",
             "Subfamily",
             "Tribe",
             "Genus",
             "Species")
}


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

results <- readRDS(paste0(case, "/result_list", natxt, ".rds"))

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
results <- lapply(results, function(x) x[which(x[,1] %in% data_true[,1]),])

# Some quick checks of results structure 
if(!length(unique(lapply(results, function(x) colnames(x)))) == 1) print("One of the results data frames have wrong column names.")
if(!all(sapply(results, function(x) all(x[,1] == data_true[,1])))) print("One of the results data frames is not sorted by ID.")
if(!all(colnames(data_true) == colnames(results[[1]])[1:ncol(data_true)])) print("Column names containing taxonomic information not identical to data_true")


# PLOT CALIBRATION CURVES ------------------------------------------------------

calibrations <- 
  get_calibration(results, data_true, observed_everywhere) |> 
  mutate(rank = factor(rank, levels = ranks), 
         set = factor(set, levels = c("All", "Observed", "Novel")), 
         model = factor(model, levels = algs_sorted))

#### Plot all calibration curves (all ranks, taxa sets and models) #####
p_cal <- plot_calibration(calibrations) 

ggsave(plot = p_cal, 
       filename = paste0("plots/calibration_all_", case, natxt, undshort, ".pdf"), 
       width = 500, 
       height = 300, 
       units = "mm")

#### Subset plot for publication: observed and novel sets only #### 
obs_perc <- round(sum(observed_everywhere)/nrow(data_true) * 100, 1)
novel_perc <- round(sum(!observed_everywhere)/nrow(data_true) * 100, 1)

p_cal_subset <- 
  plot_calibration(calibrations |> 
                     filter(set != "All")) +
  facet_grid(set~model, 
             labeller = labeller(set = c(Observed = paste0("Observed species (", obs_perc, "%)"),
                                         Novel = paste0("Novel species (", novel_perc, "%)")))) + 
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  guides(color = guide_legend(nrow = 1))

plotlist_cal <- list() 
for(i in 1:2){
  plotsets <- c("Observed", "Novel")
  
  plotlist_cal[[i]] <- 
    plot_calibration(calibrations |>
                       filter(set == plotsets[i])) +
    facet_wrap(~model, ncol = 6) + 
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom") + 
    guides(color = guide_legend(nrow = 1))
}

p_cal_subset_divided <- 
  ggpubr::ggarrange(plotlist = plotlist_cal, 
                    nrow = round(length(results)/6), 
                    common.legend = T, 
                    legend = "bottom", 
                    labels = c("a", "b"))

ggsave(plot = p_cal_subset_divided, 
       filename = paste0("plots/calibration_subset_", case, natxt, undshort, ".pdf"), 
       width = 200, 
       height = 180, 
       units = "mm")

# PLOT BINNED CALIBRATIONS -----------------------------------------------------

binned_calibrations <- get_calibration(results, data_true, observed_everywhere, 
                                       binned = T) |> 
  mutate(rank = factor(rank, levels = ranks), 
         set = factor(set, levels = c("All", "Observed", "Novel")))

binned_calibrations <- 
  binned_calibrations |> 
  mutate(bin = purrr::map_dbl(bin,
                              function(x)
                                str_remove(x, "\\(") |>
                                str_remove("\\[") |>
                                str_split(",") |>
                                unlist() |>
                                head(1) |>
                                as.numeric())) |> 
  mutate(bin = bin + 0.05)

#### Genus-level plot #### 

plotlist_cal <- list() 
plotsets <- c("Observed", "Novel")

for(i in 1:2){
  plotlist_cal[[i]] <- 
    binned_calibrations |> 
    filter(set == plotsets[i], 
           rank == "Genus") |>
    ggplot() + 
    ylim(0, 1) +
    xlim(0, 1.05) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    geom_abline(intercept = 0, slope = 1, color = "grey") +
    geom_bar(aes(x = bin, 
                 y = correct, 
                 fill = count), 
             stat = "identity", 
             alpha = 0.5, 
             width = 0.05) + 
    facet_wrap(~model, ncol = 6)  + 
    theme(axis.text.x = element_text(angle = 90)) + 
    ylab("Prop. correct") + 
    xlab("Pred. probability") +
    scale_fill_gradient(name = "# predictions",
                        trans = "log",
                        breaks = c(10, 100, 1000),
                        labels = c(10, 100, 1000))
}

p_cal_binned_genus <- 
  ggpubr::ggarrange(plotlist = plotlist_cal, 
                    nrow = round(length(results)/6), 
                    common.legend = T, 
                    legend = "right", 
                    labels = c("a", "b"))

ggsave(plot = p_cal_binned_genus, 
       filename = paste0("plots/binned_calibration_genus_", case, natxt, undshort, ".pdf"), 
       width = 250, 
       height = 200, 
       units = "mm")

#### All ranks for all models #### 
plotlist <- list()

for(m in names(results)){
  
  plotlist[[m]] <- 
    ggplot() + 
    theme_bw() +
    theme(aspect.ratio = 1) +
    geom_abline(intercept = 0, slope = 1, color = "grey") +
    geom_bar(data = binned_calibrations |> filter(model == m),
             aes(x = bin,
                 y = correct,
                 fill = log(count)),
             stat = "identity",
             alpha = 0.5,
             width = 0.05)  + 
    ylim(0, 1) +
    xlim(0, 1.05) +
    facet_grid(set~rank) + 
    theme(axis.text.x = element_text(angle = 90)) + 
    ylab("Prop. correct") + 
    xlab("Pred. probability") + 
    ggtitle(m)
}

pdf(paste0("plots/binned_calibration_all_", case, natxt, undshort, ".pdf"), 
    width = 10)
lapply(plotlist, print)
dev.off()

# PLOT ACCURACIES --------------------------------------------------------------

### If a sequence belongs to a new taxon (on any rank), how well is it predicted on higher ranks? 
accuracies <- get_accuracy_from_cal(calibrations) |> 
  mutate(rank = factor(rank, levels = ranks), 
         set = factor(set, levels = c("All", "Observed", "Novel")))

#### Plot across all models and species sets ####
p_acc <- 
  plot_accuracies(accuracies) |> 
  change_plot_colors()

ggsave(plot = p_acc, 
       filename = paste0("plots/accuracy_", case, natxt, undshort, ".pdf"), 
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
  left_join(cond_accuracies) |>
  rename(marginal = marg_accuracy,
         conditional = cond_accuracy) |> 
  pivot_longer(c("marginal", "conditional"),
               names_to = "measure",
               values_to = "accuracy") |> 
  mutate(denom_cond = if_else(measure == "marginal", 0, denom_cond)) |> 
  mutate(denom_cond = na_if(denom_cond, 0)) |> 
  mutate(rank = factor(rank, levels = ranks), 
         set = factor(set, levels = c("All", "Observed", "Novel")), 
         measure = factor(measure, levels = c("marginal", "conditional")))

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
                                      "Novel taxa")), 
         model = factor(model, levels = algs_sorted)) |> 
  ggplot() + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  aes(x = rank, 
      y = accuracy, 
      color = model, 
      shape = model) + 
  geom_line(mapping = aes(group = model)) + 
  geom_point() + 
  facet_wrap(~set) + 
  labs(color = "Model", y = "Accuracy (%)") 
fourpanel_accuracy <- change_plot_colors(fourpanel_accuracy)
  
ggsave(filename = paste0("plots/fourpanel_accuracy_", case,
                         natxt, undshort, ".pdf"), 
       plot = fourpanel_accuracy, 
       height = 150, 
       width = 200, 
       units = "mm")

#### Single panel conditional accuracy ####
no_novel <- accuracy_df |>
  filter(measure == "conditional" & set == "Novel") |>
  filter(all(accuracy == 0), .by = model) |> 
  distinct(model) |> 
  pull(model)

# TODO: Fix this plot 
#cond_accuracy_novel_taxa <- 
  accuracy_df |> 
  filter(measure == "conditional" & set == "Novel") |> 
  filter(!(model %in% no_novel)) |> 
  mutate(model_category = if_else(model %in% c("Crest4", "BLAST threshold", "DNABarcoder"), 
                                  "Similarity", 
                                  if_else(model %in% c("BayesANT", "PROTAX"), 
                                          "Probabilistic", 
                                          if_else(str_starts(model, "EPA-"), 
                                                  "Phylogenetic placement", 
                                                  "Composition")))) |> 
  mutate(rank = factor(rank, levels = ranks),
         model_category = factor(model_category, levels = c("Similarity",
                                                            "Composition", 
                                                            "Probabilistic", 
                                                            "Phylogenetic placement"))) |> 
  mutate(accuracy = if_else(is.na(denom_cond), 0, accuracy)) |> 
  mutate(denom_cond = if_else(is.na(denom_cond), 0, denom_cond)) |> 
  ggplot() + 
  aes(x = rank,
      y = accuracy * 100, 
      group = model, 
      color = model, 
      shape = model, 
      label = denom_cond) + 
  geom_line() + 
  geom_point() + 
  geom_text(nudge_y = 5, 
            size = 3, 
            show.legend = F) + 
  scale_color_manual(name = "Model",
                     labels = algs_sorted[which(!(algs_sorted %in% no_novel))],
                     limits = algs_sorted[which(!(algs_sorted %in% no_novel))],
                     values = plot_colors[which(!(algs_sorted %in% no_novel))]) + 
  scale_shape_manual(name = "Model",
                     labels = algs_sorted[which(!(algs_sorted %in% no_novel))],
                     limits = algs_sorted[which(!(algs_sorted %in% no_novel))], 
                     values = point_shapes[which(!(algs_sorted %in% no_novel))]) + 
  theme_bw() + 
  ylab("Conditional recall (%)") + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.minor.y = element_blank()) + 
  labs(color = "Model") + 
  facet_wrap(~model_category)

ggsave(plot = cond_accuracy_novel_taxa,
       filename = paste0("../plots/cond_recall", natxt, undshort, "_", case, ".pdf"),
       width = 200,
       height = 130,
       units = "mm")

# NOVELTY: TOTAL NUMBER OF PREDICTION ------------------------------------------

pred_novel <- count_novel(results, data_true, id_novel)
pred_novel$rank <- factor(pred_novel$rank, levels = names(data_true)[-1])

p_tot_novel <-
  pred_novel |> 
  filter(!all(pred_novel == 0), .by = model) |> 
  ggplot() + 
  geom_line(data = pred_novel |> filter(model == model[1]),
            aes(x = rank, 
                y = sum_novel, 
                group = model),
            color = "black",
            linetype = "dashed") + 
  aes(x = rank, 
      y = pred_novel,
      color = model, 
      group = model, 
      shape = model) +
  geom_line() + 
  geom_point() + 
  theme_bw() + 
  scale_y_log10(labels = c("10", "1000", "100000"), 
                breaks = c(10, 1000, 100000)) +
  labs(y = "# predicted novel") + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
p_tot_novel <- change_plot_colors(p_tot_novel)

ggsave(plot = p_tot_novel, 
       filename = paste0("plots/sum_novel_", case, natxt, undshort, ".pdf"), 
       height = 100, 
       width = 150, 
       units = "mm")

# NOVELTY: PARTITIONED ---------------------------------------------------------

partitioned_novelty <- 
  get_partitioned_novelty(data_true, results, observed_everywhere) |> 
  filter(pred_rank != "None")

partitioned_novelty <- 
  partitioned_novelty |> 
  mutate(pred_lower = map2_lgl(true_rank, pred_rank,
                               ~which(ranks == .x) < which(ranks == .y))) |> 
  mutate(perc = map2_dbl(pred_lower, perc,
                         ~if_else(.x, -.y, .y))) |> 
  mutate(pred_rank = if_else(pred_rank == true_rank, 
                             "Correct", pred_rank)) |> 
  mutate(true_rank = factor(true_rank, levels = colnames(data_true)[-1]), 
         pred_rank = factor(pred_rank, levels = c(colnames(data_true)[-1], 
                                                  "Correct")))

plot_colors_blue <- c(
  "#eff3ff",
  "#c6dbef",
  "#9ecae1",
  "#6baed6",
  "#4292c6",
  "#2171b5",
  "#084594"
  )

p_part_novelty <- 
  ggplot() + 
  geom_bar(data = partitioned_novelty, 
           aes(x = true_rank, 
               y = perc, 
               fill = pred_rank), 
           position = "stack", 
           stat = "identity") + 
  facet_wrap(~model) + 
  scale_fill_manual(limits = c(colnames(data_true)[-1],
                               "Correct"),
                    values = c(plot_colors_blue[1:length(colnames(data_true)[-1])],
                               "grey")) + 
  labs(fill = "Predicted novel rank",
       y = "Proportion of predictions (%)") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank()) + 
  geom_hline(yintercept = 0, 
             color = "grey20", 
             linetype = "dashed")

ggsave(plot = p_part_novelty, 
       filename = paste0("plots/part_novelty_", case, natxt, undshort, ".pdf"), 
       width = 220, 
       height = 130, 
       units = "mm")

# % CLASSIFIED ~ % CORRECT -----------------------------------------------------

# TODO: Fix this 

# Add similarity as "probability" for BLAST top hit 
result_blast_top_similarity <-
  read.table(paste0("results/blast/blast_top_hit_test", shorttxt, "_nt_16.tsv"),
             header = T) |> 
  arrange(ID) 
result_blast_top_similarity <- arrange_columns(result_blast_top_similarity, correct_cols)
result_blast_top_similarity <- 
  result_blast_top_similarity |> 
  select(ID, all_of(ranks)) |> 
  separate(Species, into = c("Species", "Similarity"), 
           sep = ";") 
result_blast_top_similarity[,correct_cols[which(str_starts(correct_cols, "Prob_"))]] <- as.numeric(result_blast_top_similarity$Similarity)/100
result_blast_top_similarity <- result_blast_top_similarity |> select(-Similarity)

add_blast_top_similarity <- data.frame(ID = data_true$ID[which(!(data_true$ID %in% result_blast_top_similarity$ID))])
add_blast_top_similarity[,correct_cols[-1]] <- NA

result_blast_top_similarity <- 
  rbind(result_blast_top_similarity, 
        add_blast_top_similarity) |> 
  arrange(ID)

result_blast_top_similarity <- 
  result_blast_top_similarity |> 
  filter(ID %in% data_true$ID)

# When NA = skip 

results_skipNA <- readRDS(paste0("result_list", "_keepNA", ".rds"))
has_nas <- names(which(sapply(results_skipNA, function(x) any(is.na(x$Species)))))
has_nas <- has_nas[which(has_nas != "BLAST top hit")]

results_similarity <- results
results_similarity[paste0(has_nas, "_NA")] <- results_skipNA[has_nas]
results_similarity$`BLAST top hit` <- result_blast_top_similarity

results_similarity <- lapply(results_similarity, function(x) x[which(x[,1] %in% data_true[,1]),])

classified_correct <- 
  threshold_curve(results_similarity, data_true) |>   
  mutate(rank = factor(rank, levels = ranks))

point_models <- 
  classified_correct |> 
    group_by(model) |> 
    distinct(correct, classified) |> 
    summarise(count = n()) |> 
    filter(count == length(ranks)) |> 
    pull(model)

classified_correct <- 
  classified_correct |>
  mutate(rank = factor(rank, levels = ranks)) |> 
  separate(model, into = c("model", "has_na"), sep = "_") |> 
  mutate(has_na = !is.na(has_na))

#p_class_correct <- 
  ggplot() +  
  geom_line(data = classified_correct |> 
              filter(!model %in% point_models), 
              aes(x = classified, 
                y = correct, 
                color = model,
                linetype = has_na),
            alpha = 0.7) +
  geom_point(data = classified_correct |> 
               filter(model %in% point_models), 
             aes(x = classified, 
                 y = correct, 
                 color = model,
                 shape = has_na)) +  
  facet_wrap(~rank, 
             nrow = 2) + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid.minor = element_blank()) + 
  labs(x = "% classified", 
       y = "% correct", 
       color = "Model") + 
  ylim(c(0, 100)) + 
  scale_color_manual(name = "Model",
                     labels = algs_sorted, 
                     values = plot_colors, 
                     breaks = algs_sorted) 

ggsave(plot = p_class_correct, 
       filename = paste0("../plots/classified_correct_COI", natxt, undshort, ".pdf"), 
       width = 250, 
       height = 130, 
       units = "mm")

# NA COUNTS --------------------------------------------------------------------

# TODO: fix this 
if(keep_na){
  p_na_count <- 
    bind_rows(results, .id = "model") |> 
    mutate(across(all_of(ranks), is.na)) |> 
    summarise(across(all_of(ranks), sum), 
              .by = model) |> 
    pivot_longer(cols = all_of(ranks), 
                 names_to = "rank", 
                 values_to = "na_count") |> 
    filter(!all(na_count == 0), .by = model) |>
    mutate(na_perc = na_count/nrow(data_true) * 100) |> 
    mutate(rank = factor(rank, levels = ranks)) |> 
    ggplot() + 
    geom_bar(aes(x = rank, 
                 y = na_perc), 
             stat = "identity") + 
    facet_wrap(~model) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ylim(c(0, 100)) + 
    ylab("Missing predictions (%)")
  
  ggsave(plot = p_na_count,
         filename = paste0("../plots/na_count", undshort, ".pdf"), 
         width = 6, 
         height = 5, 
         units = "in")
}

# PREDICTION PROBABILITY -------------------------------------------------------

p_pred_prob <- 
  bind_rows(results, .id = "model") |> 
  pivot_longer(cols = all_of(paste0("Prob_", ranks)), 
               names_to = "rank", 
               values_to = "probability", 
               names_prefix = "Prob_") |> 
  mutate(rank = factor(rank, levels = ranks)) |> 
  filter(!all(probability == 1), .by = model) |> 
  ggplot() + 
  geom_histogram(aes(x = probability), 
                 fill = "steelblue", 
                 alpha = 0.7) + 
  facet_grid(rank~model, scales = "free") + 
  ylab("") + 
  xlab("Probability") + 
  theme_bw() +
  scale_y_log10()

ggsave(filename = paste0("plots/pred_prob_", case, undshort, natxt, ".pdf"), 
       plot = p_pred_prob, 
       width = 350, 
       height = 150, 
       units = "mm")

# MIS-, OVER-, AND UNDERCLASSIFICATION -----------------------------------------

oclass_df <- overclass_rate(results, id_novel)
uclass_df <- underclass_rate(results, id_observed)
misclass_df <- misclass_rate(results, data_true, id_observed)

errorrates <- 
  oclass_df |> 
  full_join(uclass_df) |> 
  full_join(misclass_df) |> 
  mutate(rank = factor(rank, levels = ranks)) |> 
  pivot_longer(c("overclassification", "underclassification", "misclassification"), 
               names_to = "error_type",
               values_to = "error_rate") |> 
  mutate(error_type = str_to_sentence(error_type)) |> 
  ggplot() + 
  geom_line(aes(x = rank, 
                y = error_rate, 
                color = model, 
                group = model)) +
  geom_point(aes(x = rank, 
                 y = error_rate, 
                 color = model, 
                 shape = model)) + 
  facet_wrap(~error_type) + 
  theme_bw() + 
  ylab("Error rate (%)") + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1), 
        aspect.ratio = 1)
errorrates <- change_plot_colors(errorrates)

ggsave(filename = paste0("plots/errorrates_", case, undshort, natxt, ".pdf"), 
       plot = errorrates, 
       width = 10,
       height = 4, 
       units = "in")

# REPRESENTATIVE SEQ -----------------------------------------------------------

p_repseq <- 
  data$train |> 
  as_tibble() |> 
  pivot_longer(2:8, 
               names_to = "rank", 
               values_to = "taxon") |> 
  group_by(rank, taxon) |> 
  summarise(n_rep = n()) |>
  filter(n_rep <= 5) |> 
  mutate(rank = factor(rank, levels = ranks)) |> 
  ggplot() + 
  aes(x = n_rep) +
  geom_histogram(binwidth = 1, 
                 fill = "steelblue", 
                 color = "black") +
  facet_wrap(~rank, scales = "free_y", nrow = 2) + 
  theme_bw() + 
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) + 
  ylab("# of taxa") + 
  xlab("# representative sequences")

ggsave(filename = paste0("plots/repseq_", case, undshort, natxt, ".pdf"), 
       plot = p_repseq, 
       width = 6,
       height = 4, 
       units = "in")


# FOR THE TEXT -----------------------------------------------------------------
# TODO: fix this 
accuracies |> 
  filter(set == "Observed") |> 
  ggplot() + 
  aes(x = model, 
      y = accuracy, 
      color = model) + 
  geom_point() + 
  facet_wrap(~rank, scales = "free_y") + 
  theme(axis.text.x = element_text(angle = 90))
  
accuracies |> 
  group_by(rank) |> 
  filter(#model != "Crest4",
         set == "Observed") |> 
  summarise(min(accuracy))

accuracy_df |> 
  filter(set == "Observed", 
         measure == "marginal") |> 
  mutate(model = factor(model, levels = algs_sorted)) |> 
  ggplot() + 
  aes(x = model, 
      y = accuracy*100, 
      color = model) + 
  geom_point() + 
  facet_wrap(~rank, scales = "free_y") + 
  theme(axis.text.x = element_text(angle = 90))

accuracies |> 
  filter(model %in% c("EPA-ng taxtree", "EPA-ng phyltree", "BayesANT", "PROTAX"),
         set == "Novel") |> 
  group_by(rank) |>
  summarise(max = max(accuracy), 
            min = min(accuracy)) 

accuracies |> 
  filter(set == "Novel") |> 
  ggplot() + 
  aes(x = model, 
      y = accuracy, 
      color = model) + 
  geom_point() + 
  facet_wrap(~rank, scales = "free_y") + 
  theme(axis.text.x = element_text(angle = 90))
