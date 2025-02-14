
rm(list = ls())

# LOAD FUNCTIONS AND PACKAGES --------------------------------------------------

library(tidyverse)
source("functions.R")

# LOAD DATA --------------------------------------------------------------------

ranks_its <- 
  c("Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species")

ranks_coi <- 
  c("Class", 
    "Order",
    "Family",
    "Subfamily",
    "Tribe",
    "Genus",
    "Species")

# Load test and train data
data_its <- 
  load_train_test(train_file = "its/data/train_tax.tsv", 
                  test_file = "its/data/test_tax.tsv",
                  ranks = ranks_its)

data_coi <- 
  load_train_test(train_file = "coi/data/train_tax.tsv", 
                  test_file = "coi/data/test_tax.tsv",
                  ranks = ranks_coi)

# Replace labels of taxa unique to test
# Sort test data in alphabetical order
data_true_its <- 
  get_data_true(data_its$test,
                data_its$train,
                ranks = ranks_its) |>
  dplyr::arrange(ID)

data_true_coi <- 
  get_data_true(data_coi$test,
                data_coi$train,
                ranks = ranks_coi) |>
  dplyr::arrange(ID)

# CORRECT COLUMN NAMES AND ORDER -----------------------------------------------

correct_cols_its <- 
  c("ID",
    ranks_its,
    sapply(ranks_its, function(x) paste0("Prob_", x), USE.NAMES = F))

correct_cols_coi <- 
  c("ID",
    ranks_coi,
    sapply(ranks_coi, function(x) paste0("Prob_", x), USE.NAMES = F))

# READ RESULTS -----------------------------------------------------------------

# Without any missing predictions
results_its <- readRDS("its/result_list.rds")
results_coi <- readRDS("coi/result_list_nt.rds")

# Allowing missing predictions (and using author-recommended thresholds)
results_its_mp <- readRDS("its/result_list_keepNA.rds")
results_coi_mp <- readRDS("coi/result_list_keepNA_nt.rds")

# IDENTIFY TAXA SETS -----------------------------------------------------------

# Observed on all ranks ("Observed species")
observed_everywhere_its <- !grepl("_new$", data_true_its[, ncol(data_true_its)])
observed_everywhere_coi <- !grepl("_new$", data_true_coi[, ncol(data_true_coi)])

id_all_its <- lapply(data_true_its, function(x) 1:length(x))
id_all_coi <- lapply(data_true_coi, function(x) 1:length(x))

# ITS 
# Observed and novel at each rank ("Novel taxa")
id_novel_its <- list()
id_observed_its <- list()
for(rank in colnames(data_true_its)[-1]){
  id_novel_its[[rank]] <- which(stringr::str_ends(data_true_its[,rank], paste0(rank, "_new")))
  id_observed_its[[rank]] <- which(!stringr::str_ends(data_true_its[,rank], "_new"))
}

# COI 
# Observed and novel at each rank ("Novel taxa")
id_novel_coi <- list()
id_observed_coi <- list()
for(rank in colnames(data_true_coi)[-1]){
  id_novel_coi[[rank]] <- which(stringr::str_ends(data_true_coi[,rank], paste0(rank, "_new")))
  id_observed_coi[[rank]] <- which(!stringr::str_ends(data_true_coi[,rank], "_new"))
}

# Check if length are the same for all 
lapply(results_its, function(x) nrow(x)) |> unique()
lapply(results_its_mp, function(x) nrow(x)) |> unique()

lapply(results_coi, function(x) nrow(x)) |> unique()
lapply(results_coi_mp, function(x) nrow(x)) |> unique()

# TODO: Some of our results have more sequences than data_true
# Here, I remove them as a temporary solution 
results_coi <- lapply(results_coi, function(x) x[which(x[,1] %in% data_true_coi[,1]),])
results_coi_mp <- lapply(results_coi_mp, function(x) x[which(x[,1] %in% data_true_coi[,1]),])

# PLOT CALIBRATION CURVES ------------------------------------------------------

# Without allowing missing predictions 
calibrations_coi <- 
  get_calibration(results_coi, data_true_coi, observed_everywhere_coi) |> 
  mutate(rank = factor(rank, levels = ranks_coi), 
         set = factor(set, levels = c("All", "Observed", "Novel")), 
         model = factor(model, levels = algs_sorted))
calibrations_its <- 
  get_calibration(results_its, data_true_its, observed_everywhere_its) |> 
  mutate(rank = factor(rank, levels = ranks_its), 
         set = factor(set, levels = c("All", "Observed", "Novel")), 
         model = factor(model, levels = algs_sorted))

#### Novel species #### 

p_cal_novel <- list()

p_cal_novel[["coi"]] <- 
  plot_calibration(calibrations_coi |> 
                     filter(set == "Novel"), 
                   data_true = data_true_coi) +
  facet_wrap(~model, 
             nrow = 2) + 
  theme(panel.grid.minor = element_blank(),
        legend.position = "right") + 
  guides(color = guide_legend(ncol = 1))

p_cal_novel[["its"]] <- 
  plot_calibration(calibrations_its |> 
                     filter(set == "Novel"), 
                   data_true = data_true_its) +
  facet_wrap(~model, 
             nrow = 2) + 
  theme(panel.grid.minor = element_blank(),
        legend.position = "right") + 
  guides(color = guide_legend(ncol = 1))

p_cal_novel_joint <- 
  ggpubr::ggarrange(plotlist = p_cal_novel, 
                    ncol = 1, 
                    labels = str_to_upper(names(p_cal)))

ggsave(plot = p_cal_novel_joint, 
       filename = "plots/calibration_novel_joint.pdf", 
       width = 10, 
       height = 8, 
       units = "in")

#### Observed species #### 

p_cal_obs <- list()

p_cal_obs[["coi"]] <- 
  plot_calibration(calibrations_coi |> 
                     filter(set == "Observed"), 
                   data_true = data_true_coi) +
  facet_wrap(~model, 
             nrow = 2) + 
  theme(panel.grid.minor = element_blank(),
        legend.position = "right") + 
  guides(color = guide_legend(ncol = 1))

p_cal_obs[["its"]] <- 
  plot_calibration(calibrations_its |> 
                     filter(set == "Observed"), 
                   data_true = data_true_its) +
  facet_wrap(~model, 
             nrow = 2) + 
  theme(panel.grid.minor = element_blank(),
        legend.position = "right") + 
  guides(color = guide_legend(ncol = 1))

p_cal_obs_joint <- 
  ggpubr::ggarrange(plotlist = p_cal_obs, 
                    ncol = 1, 
                    labels = str_to_upper(names(p_cal)))

ggsave(plot = p_cal_obs_joint, 
       filename = "plots/calibration_obs_joint.pdf", 
       width = 10, 
       height = 8, 
       units = "in")

# PLOT ACCURACIES --------------------------------------------------------------

accuracies_coi <-
  get_accuracy_from_cal(calibrations_coi) |> 
  mutate(rank = factor(rank, levels = ranks_coi), 
         set = factor(set, levels = c("All", "Observed", "Novel")))

accuracies_its <-
  get_accuracy_from_cal(calibrations_its) |> 
  mutate(rank = factor(rank, levels = ranks_its), 
         set = factor(set, levels = c("All", "Observed", "Novel")))

marg_accuracies_coi <- 
  marginal_accuracy(results_coi, 
                    data_true_coi,
                    id_list = list(All = id_all_coi,
                                   Observed = id_observed_coi,
                                   Novel = id_novel_coi))
marg_accuracies_its <- 
  marginal_accuracy(results_its, 
                    data_true_its,
                    id_list = list(All = id_all_its,
                                   Observed = id_observed_its,
                                   Novel = id_novel_its))

p_acc <- list()

p_acc[["coi"]] <- 
  marg_accuracies_coi |> 
  filter(set != "All") |> 
  rename(accuracy = marg_accuracy) |> 
  mutate(set = paste(set, "taxa"), 
         accuracy = accuracy * 100) |> 
  bind_rows(accuracies_coi |> 
              filter(set != "All") |> 
              mutate(set = paste(set, "species"))) |> 
    mutate(set = factor(set, levels = c("Observed species",
                                        "Novel species",
                                        "Observed taxa",
                                        "Novel taxa")), 
           model = factor(model, levels = algs_sorted), 
           rank = factor(rank, levels = ranks_coi)) |> 
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
p_acc[["coi"]] <- change_plot_colors(p_acc[["coi"]])

p_acc[["its"]] <- 
  marg_accuracies_its |> 
  filter(set != "All") |> 
  rename(accuracy = marg_accuracy) |> 
  mutate(set = paste(set, "taxa"), 
         accuracy = accuracy * 100) |> 
  bind_rows(accuracies_its |> 
              filter(set != "All") |> 
              mutate(set = paste(set, "species"))) |> 
  mutate(set = factor(set, levels = c("Observed species",
                                      "Novel species",
                                      "Observed taxa",
                                      "Novel taxa")), 
         model = factor(model, levels = algs_sorted), 
         rank = factor(rank, levels = ranks_its)) |> 
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
p_acc[["its"]] <- change_plot_colors(p_acc[["its"]])

p_acc_joint <- 
  ggpubr::ggarrange(plotlist = p_acc, 
                  ncol = 2, 
                  common.legend = T, 
                  labels = str_to_upper(names(p_acc)), 
                  legend = "bottom")
ggsave(plot = p_acc_joint, 
       filename = "plots/fourpanel_joint.pdf", 
       width = 9, 
       height = 6, 
       units = "in")


# % CLASSIFIED ~ % CORRECT -----------------------------------------------------

# Use results which allow missing predictions

# Add similarity as "probability" for BLAST top hit 
result_blast_top_similarity_coi <-
  read.table("coi/results/blast/blast_top_hit_test_nt_16.tsv",
             header = T) |> 
  arrange_columns(correct_cols = correct_cols_coi) |> 
  separate(Species, into = c("Species", "Similarity"), 
           sep = ";") %>%
  mutate(Similarity = as.numeric(Similarity)/100) |> 
  mutate(across(starts_with("Prob_"), ~Similarity)) |> 
  select(-Similarity) %>%
  bind_rows(data.frame(ID = data_true_coi$ID[which(!(data_true_coi$ID %in% .$ID))])) |> 
  arrange(ID) 
result_blast_top_similarity_its <-
  read.table("its/results/blast/blast_top_hit_test_nt_16.tsv",
             header = T) |> 
  arrange_columns(correct_cols = correct_cols_its) |> 
  separate(Species, into = c("Species", "Similarity"), 
           sep = ";") %>%
  mutate(Similarity = as.numeric(Similarity)/100) |> 
  mutate(across(starts_with("Prob_"), ~Similarity)) |> 
  select(-Similarity) %>%
  bind_rows(data.frame(ID = data_true_its$ID[which(!(data_true_its$ID %in% .$ID))])) |> 
  arrange(ID) 

# Replace BLAST result with similarity data frame 
results_coi_mp_similarity <- results_coi_mp
results_coi_mp_similarity$`BLAST top hit` <- result_blast_top_similarity_coi
results_its_mp_similarity <- results_its_mp
results_its_mp_similarity$`BLAST top hit` <- result_blast_top_similarity_its

classified_correct_coi <- 
  threshold_curve(results_coi_mp_similarity, data_true_coi) |>   
  mutate(rank = factor(rank, levels = ranks_coi))
classified_correct_its <- 
  threshold_curve(results_its_mp_similarity, data_true_its) |>   
  mutate(rank = factor(rank, levels = ranks_its))

point_models_coi <- 
  classified_correct_coi |> 
  group_by(model) |> 
  distinct(correct, classified) |> 
  summarise(count = n()) |> 
  filter(count == length(ranks_coi)) |> 
  pull(model)
point_models_its <- 
  classified_correct_its |> 
  group_by(model) |> 
  distinct(correct, classified) |> 
  summarise(count = n()) |> 
  filter(count == length(ranks_its)) |> 
  pull(model)

p_class <- list()

p_class[["coi"]] <- 
  ggplot() +  
  geom_line(data = classified_correct_coi |> 
              filter(!model %in% point_models_coi), 
            aes(x = classified, 
                y = correct, 
                color = model),
            alpha = 0.7) +
  geom_point(data = classified_correct_coi |> 
               filter(model %in% point_models_coi), 
             aes(x = classified, 
                 y = correct, 
                 color = model)) +  
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

p_class[["its"]] <- 
  ggplot() +  
  geom_line(data = classified_correct_its |> 
              filter(!model %in% point_models_its), 
            aes(x = classified, 
                y = correct, 
                color = model),
            alpha = 0.7) +
  geom_point(data = classified_correct_its |> 
               filter(model %in% point_models_its), 
             aes(x = classified, 
                 y = correct, 
                 color = model)) +  
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

p_class_joint <- 
  ggpubr::ggarrange(plotlist = p_class, 
                  ncol = 1, 
                  labels = str_to_upper(names(p_class)), 
                  common.legend = T, 
                  legend = "right")

ggsave(plot = p_class_joint, 
       filename = "plots/classified_correct_joint.pdf", 
       width = 8, 
       height = 7, 
       units = "in")

# MIS-, OVER-, AND UNDERCLASSIFICATION -----------------------------------------

# Here we use thresholds for IDTAXA, RDP, SINTAX
# (novel and NA will have the same interpretation here)
oclass_df_coi <- overclass_rate(results_coi_mp, id_novel_coi, ranks = ranks_coi)
uclass_df_coi <- underclass_rate(results_coi_mp, id_observed_coi, ranks = ranks_coi)
misclass_df_coi <- misclass_rate(results_coi_mp, data_true_coi, id_observed_coi, ranks = ranks_coi)

oclass_df_its <- overclass_rate(results_its_mp, id_novel_its, ranks = ranks_its)
uclass_df_its <- underclass_rate(results_its_mp, id_observed_its, ranks = ranks_its)
misclass_df_its <- misclass_rate(results_its_mp, data_true_its, id_observed_coi, ranks = ranks_its)

p_error <- list()

p_error[["coi"]] <- 
  oclass_df_coi |> 
  full_join(uclass_df_coi) |> 
  full_join(misclass_df_coi) |> 
  mutate(rank = factor(rank, levels = ranks_coi)) |> 
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
p_error[["coi"]] <- change_plot_colors(p_error[["coi"]])

p_error[["its"]] <- 
  oclass_df_its |> 
  full_join(uclass_df_its) |> 
  full_join(misclass_df_its) |> 
  mutate(rank = factor(rank, levels = ranks_its)) |> 
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
p_error[["its"]] <- change_plot_colors(p_error[["its"]])

p_error_joint <- 
  ggpubr::ggarrange(plotlist = p_error, 
                  ncol = 1, 
                  common.legend = T, 
                  legend = "right", 
                  labels = str_to_upper(names(p_error)))

ggsave(plot = p_error_joint, 
       filename = "plots/errorrates_joint.pdf", 
       height = 6, 
       width = 8, 
       units = "in")
