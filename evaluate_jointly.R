
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

allranks <- 
  c("Kingdom",
    "Phylum",
    "Class", 
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

caldf <- 
  calibrations_coi |> 
  filter(set != "All") |> 
  mutate(gene = "COI") |> 
  bind_rows(calibrations_its |> 
              filter(set != "All") |> 
              mutate(gene = "ITS")) |> 
  filter(rank == "Genus") 
  
p_cal <- 
  ggplot() +
  xlim(0, 100) +
  ylim(0, 100) +
  theme_bw() +
  theme(aspect.ratio = 1, 
        panel.spacing = unit(0, "lines")) +
  geom_abline(intercept = 0, slope = 1, color = "grey") +
  xlab("Cumulative probability %") +
  ylab("Cumulative correct %") +
  facet_grid(set~gene) + 
  #grids(linetype = "dashed") + 
  geom_line(data = caldf, 
            aes(x = cumprob, 
                y = cumcorr, 
                color = model)) +
  geom_point(data = caldf |>
               group_by(across(all_of(c("model", "gene", "set")))) |>
               summarise(cumprob = max(cumprob),
                         cumcorr = max(cumcorr)), 
             aes(x = cumprob, 
                 y = cumcorr, 
                 color = model, 
                 shape = model)) 

p_cal <- change_plot_colors(p_cal)

ggsave(plot = p_cal, 
       filename = "plots/calibration_joint.pdf", 
       width = 6, 
       height = 5, 
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


p_acc <- 
  marg_accuracies_coi |> 
  filter(set != "All") |> 
  rename(accuracy = marg_accuracy) |> 
  mutate(set = paste(set, "Taxa"), 
         accuracy = accuracy * 100) |> 
  bind_rows(accuracies_coi |> 
              filter(set != "All") |> 
              mutate(set = paste(set, "Species"))) |> 
    mutate(gene = "COI") |> 
    bind_rows(marg_accuracies_its |> 
                filter(set != "All") |> 
                rename(accuracy = marg_accuracy) |> 
                mutate(set = paste(set, "Taxa"), 
                       accuracy = accuracy * 100) |> 
                bind_rows(accuracies_its |> 
                            filter(set != "All") |> 
                            mutate(set = paste(set, "Species"))) |> 
                mutate(gene = "ITS")) |> 
    separate(set, into = c("obs", "set"), sep = " ") |> 
    mutate(rank = factor(rank, levels = allranks), 
           obs = factor(obs, levels = c("Observed", "Novel")), 
           model = factor(model, levels = algs_sorted)) |> 
    ggplot() + 
    aes(x = rank, 
        y = accuracy, 
        color = model, 
        shape = model) +
    geom_line(mapping = aes(group = model)) + 
    geom_point() +
    theme_bw() + 
    theme(aspect.ratio = 1, 
          panel.grid.minor = element_blank(), 
          axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.position = "bottom",
          panel.spacing = unit(0, "lines")) + 
    ggh4x::facet_nested(set ~ gene + obs, scales = "free_x") + 
    labs(color = "Model", y = "Marginal recall (%)                  Accuracy (%)") # Stupid double axis label fix 
p_acc <- change_plot_colors(p_acc) 

ggsave(plot = p_acc, 
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

point_models <- c(
  "BLAST threshold",
  "Crest4", 
  "DNABarcoder"
)

classified_correct <- 
  classified_correct_coi |> 
  mutate(gene = "COI") |> 
  bind_rows(classified_correct_its |> 
              mutate(gene = "ITS")) |> 
  mutate(rank = factor(rank, levels = allranks)) |> 
  filter(rank %in% c("Genus", "Species"))

p_class <- 
  ggplot() +  
  geom_line(data = classified_correct |> 
              filter(!model %in% point_models), 
            aes(x = classified, 
                y = correct, 
                color = model, 
                shape = model),
            alpha = 0.7) +
  geom_point(data = classified_correct |> 
               filter(model %in% point_models), 
             aes(x = classified, 
                 y = correct, 
                 color = model, 
                 shape = model)) + 
  facet_grid(rank~gene) + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid.minor = element_blank()) + 
  labs(x = "% classified", 
       y = "% correct", 
       color = "Model") + 
  ylim(c(0, 100)) 
p_class <- change_plot_colors(p_class)


ggsave(plot = p_class, 
       filename = "plots/classified_correct_joint.pdf", 
       width = 6, 
       height = 5, 
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


p_error <- 
  oclass_df_coi |> 
  full_join(uclass_df_coi) |> 
  full_join(misclass_df_coi) |> 
  mutate(gene = "COI") |> 
    bind_rows(oclass_df_its |> 
                full_join(uclass_df_its) |> 
                full_join(misclass_df_its) |> 
                mutate(gene = "ITS")) |> 
  mutate(rank = factor(rank, levels = allranks)) |> 
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
  facet_grid(error_type~gene, scales = "free_x") +
  theme_bw() + 
  ylab("Error rate (%)") + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1), 
        aspect.ratio = 1, 
        panel.spacing = unit(0, "lines"))
p_error <- change_plot_colors(p_error)

ggsave(plot = p_error, 
       filename = "plots/errorrates_joint.pdf", 
       height = 6.5, 
       width = 6.5, 
       units = "in")
