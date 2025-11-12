
rm(list = ls())

# LOAD FUNCTIONS AND PACKAGES --------------------------------------------------

library(tidyverse)
source("functions.R")

# LOAD DATA --------------------------------------------------------------------

ranks_its <- 
  c(#"Kingdom",
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
  c(#"Kingdom",
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
                  ranks = c("Kingdom", ranks_its))

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

# Remove Kingdom from ITS data 
results_its <- lapply(results_its, function(df) df[,correct_cols_its])
results_its_mp <- lapply(results_its_mp, function(df) df[,correct_cols_its])

# # Do not evaluate incertae sedis 
# sorted_test_coi <- data_coi$test |> arrange(id)
# incertae_coi <- lapply(as.list(sorted_test_coi[2:ncol(sorted_test_coi)]), function(x) which(str_ends(x, "_incertae_sedis")))
# 
# results_coi <- 
#   lapply(results_coi, function(df) {
#     for (colname in names(incertae_coi)) {
#       if (colname %in% names(df)) {
#         df[incertae_coi[[colname]], colname] <- NA
#       }
#       }
#     df
#     })
# results_coi_mp <- 
#   lapply(results_coi_mp, function(df) {
#     for (colname in names(incertae_coi)) {
#       if (colname %in% names(df)) {
#         df[incertae_coi[[colname]], colname] <- NA
#       }
#     }
#     df
#   })

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

# # TODO: Some of our results have more sequences than data_true
# # Here, I remove them as a temporary solution 
# results_coi <- lapply(results_coi, function(x) x[which(x[,1] %in% data_true_coi[,1]),])
# results_coi_mp <- lapply(results_coi_mp, function(x) x[which(x[,1] %in% data_true_coi[,1]),])

# FOR PLOTTING -----------------------------------------------------------------

model_details <- create_model_details(unique(c(names(results_coi), 
                                               names(results_its))))

models <- split(model_details$model, model_details$type, drop = TRUE)
models <- models[unique(model_details$type)]
      

# PLOT CALIBRATION CURVES ------------------------------------------------------

# Without allowing missing predictions 
calibrations_coi <- 
  get_calibration(results_coi, data_true_coi, observed_everywhere_coi) |> 
  mutate(rank = factor(rank, levels = ranks_coi), 
         set = factor(set, levels = c("All", "Observed", "Novel")), 
         model = factor(model, levels = model_details$model))
calibrations_its <- 
  get_calibration(results_its, data_true_its, observed_everywhere_its) |> 
  mutate(rank = factor(rank, levels = ranks_its), 
         set = factor(set, levels = c("All", "Observed", "Novel")), 
         model = factor(model, levels = model_details$model))

# With missing predictions
calibrations_coi_mp <- 
  get_calibration(results_coi_mp, data_true_coi, observed_everywhere_coi) |> 
  mutate(rank = factor(rank, levels = ranks_coi), 
         set = factor(set, levels = c("All", "Observed", "Novel")), 
         model = factor(model, levels = model_details$model))
calibrations_its_mp <- 
  get_calibration(results_its_mp, data_true_its, observed_everywhere_its) |> 
  mutate(rank = factor(rank, levels = ranks_its), 
         set = factor(set, levels = c("All", "Observed", "Novel")), 
         model = factor(model, levels = model_details$model))

#### Genus-level calibrations #### 

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
  facet_grid(set~gene) 

for(type in names(models)){
  p_cal <- 
    p_cal + 
    plot_line(df = caldf |> 
                filter(model %in% models[[type]]), 
              x = cumprob, 
              y = cumcorr, 
              model_type = type, 
              plot_details = model_details) + 
    plot_point(df = caldf |>
                 filter(model %in% models[[type]]) |>
                 group_by(across(all_of(c("model", "gene", "set")))) |>
                 summarise(cumprob = max(cumprob),
                           cumcorr = max(cumcorr)),
               x = cumprob,
               y = cumcorr,
               model_type = type,
               plot_details = model_details) +
    ggnewscale::new_scale_color() + ggnewscale::new_scale("shape") 
}

p_cal <- 
  p_cal + 
  theme(legend.spacing.x = unit(0.1, "lines"), 
        legend.title = element_text(size = 10), 
        legend.position = "bottom", 
        legend.box = "horizontal", 
        legend.title.position = "top")

ggsave(plot = p_cal, 
       filename = "plots/calibration_joint.pdf", 
       width = 6.5, 
       height = 6.5, 
       units = "in")

#### All ranks #### 

p_cal_all <- list()

p_cal_all[["coi_obs"]] <- 
  plot_calibration(calibrations_coi |> 
                     filter(set == "Observed"), 
                   data_true = data_true_coi) +
  ggplot2::scale_color_manual(name = "Rank", 
                              limits = ranks_coi, 
                              breaks = ranks_coi, 
                              values = RColorBrewer::brewer.pal(name = "Blues",
                                                                n = length(ranks_coi) + 1)[2:(length(ranks_coi) + 1)])

p_cal_all[["its_obs"]] <- 
  plot_calibration(calibrations_its |> 
                     filter(set == "Observed"), 
                   data_true = data_true_its) +
  ggplot2::scale_color_manual(name = "Rank", 
                              limits = ranks_its, 
                              breaks = ranks_its, 
                              values = RColorBrewer::brewer.pal(name = "Blues",
                                                                n = length(ranks_its) + 1)[2:(length(ranks_its) + 1)])

p_cal_all[["coi_nov"]] <- 
  plot_calibration(calibrations_coi |> 
                     filter(set == "Novel"), 
                   data_true = data_true_coi) +
  ggplot2::scale_color_manual(name = "Rank", 
                              limits = ranks_coi, 
                              breaks = ranks_coi, 
                              values = RColorBrewer::brewer.pal(name = "Blues",
                                                                n = length(ranks_coi) + 1)[2:(length(ranks_coi) + 1)])

p_cal_all[["its_nov"]] <- 
  plot_calibration(calibrations_its |> 
                     filter(set == "Novel"), 
                   data_true = data_true_its) +
  ggplot2::scale_color_manual(name = "Rank", 
                              limits = ranks_its, 
                              breaks = ranks_its, 
                              values = RColorBrewer::brewer.pal(name = "Blues",
                                                                n = length(ranks_its) + 1)[2:(length(ranks_its) + 1)])
p_cal_all <- 
  lapply(p_cal_all, function(p){
  p <- 
    p + 
    facet_wrap(~model, 
               ncol = 5) + 
    theme(panel.grid.minor = element_blank(),
          legend.position = "right") 
  p
})

p_cal_all_coi <- 
  ggpubr::ggarrange(plotlist = p_cal_all[c(1,3)], 
                    ncol = 1, 
                    hjust = 0,
                    labels = c("COI_obs", "COI_novel"), 
                    common.legend = T, 
                    legend = "right")

ggsave(plot = p_cal_all_coi, 
       filename = "plots/calibration_coi_all.pdf", 
       width = 8, 
       height = 9, 
       units = "in")

p_cal_all_its <- 
  ggpubr::ggarrange(plotlist = p_cal_all[c(2,4)], 
                    ncol = 1, 
                    hjust = 0,
                    labels = c("ITS_obs", "ITS_novel"), 
                    common.legend = T, 
                    legend = "right")

ggsave(plot = p_cal_all_its, 
       filename = "plots/calibration_its_all.pdf", 
       width = 8, 
       height = 9, 
       units = "in")

# PLOT ACCURACIES --------------------------------------------------------------

##### Forced predictions #####
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

p_acc_df <- 
  marg_accuracies_coi |> 
  filter(set != "All") |> 
  rename(accuracy = marg_accuracy) |> 
  mutate(set = paste(set, "Taxa"), 
         accuracy = accuracy * 100) |> 
  bind_rows(accuracies_coi |> 
              filter(set != "All") |> 
              mutate(set = map_chr(set, ~if_else(.x == "Observed", 
                                                 paste(.x, "everywhere"),
                                                 paste(.x, "anywhere"))))) |> 
  mutate(gene = "COI") |> 
  bind_rows(marg_accuracies_its |> 
              filter(set != "All") |> 
              rename(accuracy = marg_accuracy) |> 
              mutate(set = paste(set, "Taxa"), 
                     accuracy = accuracy * 100) |> 
              bind_rows(accuracies_its |> 
                          filter(set != "All") |> 
                          mutate(set = map_chr(set, ~if_else(.x == "Observed", 
                                                             paste(.x, "everywhere"),
                                                             paste(.x, "anywhere"))))) |> 
              mutate(gene = "ITS"))  |> 
  separate(set, into = c("obs", NA), sep = " ", remove = F) |> 
  mutate(rank = factor(rank, levels = allranks), 
         obs = factor(obs, levels = c("Observed", "Novel")), 
         model = factor(model, levels = model_details$model)) 

p_acc <- list(Species = NA, 
              Taxa = NA)

p_acc <- 
  lapply(p_acc, function(x){
  p <- ggplot() +
    theme_bw() + 
    theme(aspect.ratio = 1, 
          panel.grid.minor = element_blank(), 
          axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.spacing = unit(0, "lines"),
          strip.text.y = element_blank(),legend.spacing.x = unit(0.1, "lines"), 
          legend.title = element_text(size = 10), 
          legend.position = "bottom", 
          legend.box = "horizontal", 
          legend.title.position = "top") +
    ggh4x::facet_nested(~ gene + set, scales = "free_x") + 
    labs(color = "Model") 
  p
})

for(s in names(p_acc)){
  
  if(s == "Taxa"){
    subdf <- 
      p_acc_df |> 
      filter(str_ends(set, s))
  }else{
    subdf <- 
      p_acc_df |> 
      filter(str_ends(set, "where"))
  }
  
  for(type in names(models)){
    p_acc[[s]] <- 
      p_acc[[s]] + 
      plot_line(df = subdf |> filter(model %in% models[[type]]), 
                x = rank, 
                y = accuracy, 
                model_type = type, 
                plot_details = model_details, 
                group = model) +
      plot_point(df = subdf |> filter(model %in% models[[type]]), 
                 x = rank, 
                 y = accuracy, 
                 model_type = type, 
                 plot_details = model_details) + 
      ggnewscale::new_scale_color() + ggnewscale::new_scale("shape") 
  }
  if(s == "Species"){
    p_acc[[s]] <- 
      p_acc[[s]] + 
      labs(y = "Accuracy (%)")
  }else{
    p_acc[[s]] <- 
      p_acc[[s]] + 
      labs(y = "Marginal recall (%)")
  }
}

p_acc_joint <- 
  ggpubr::ggarrange(plotlist = p_acc, 
                  ncol = 1, 
                  labels = c("a", "b"), 
                  common.legend = T, 
                  legend = "bottom")

ggsave(plot = p_acc_joint, 
       filename = "plots/fourpanel_joint.pdf", 
       width = 7, 
       height = 7, 
       units = "in")

##### Missing predictions #####
accuracies_coi_mp <-
  get_accuracy_from_cal(calibrations_coi_mp) |> 
  mutate(rank = factor(rank, levels = ranks_coi), 
         set = factor(set, levels = c("All", "Observed", "Novel")))

accuracies_its_mp <-
  get_accuracy_from_cal(calibrations_its_mp) |> 
  mutate(rank = factor(rank, levels = ranks_its), 
         set = factor(set, levels = c("All", "Observed", "Novel")))

marg_accuracies_coi_mp <- 
  marginal_accuracy(results_coi_mp, 
                    data_true_coi,
                    id_list = list(All = id_all_coi,
                                   Observed = id_observed_coi,
                                   Novel = id_novel_coi))
marg_accuracies_its_mp <- 
  marginal_accuracy(results_its_mp, 
                    data_true_its,
                    id_list = list(All = id_all_its,
                                   Observed = id_observed_its,
                                   Novel = id_novel_its))

p_acc_mp_df <- 
  marg_accuracies_coi_mp |> 
  filter(set != "All") |> 
  rename(accuracy = marg_accuracy) |> 
  mutate(set = paste(set, "Taxa"), 
         accuracy = accuracy * 100) |> 
  bind_rows(accuracies_coi_mp |> 
              filter(set != "All") |> 
              mutate(set = map_chr(set, ~if_else(.x == "Observed", 
                                                 paste(.x, "everywhere"),
                                                 paste(.x, "anywhere"))))) |> 
  mutate(gene = "COI") |> 
  bind_rows(marg_accuracies_its_mp |> 
              filter(set != "All") |> 
              rename(accuracy = marg_accuracy) |> 
              mutate(set = paste(set, "Taxa"), 
                     accuracy = accuracy * 100) |> 
              bind_rows(accuracies_its_mp |> 
                          filter(set != "All") |> 
                          mutate(set = map_chr(set, ~if_else(.x == "Observed", 
                                                             paste(.x, "everywhere"),
                                                             paste(.x, "anywhere"))))) |> 
              mutate(gene = "ITS")) |> 
  separate(set, into = c("obs", NA), sep = " ", remove = F) |> 
  mutate(rank = factor(rank, levels = allranks), 
         obs = factor(obs, levels = c("Observed", "Novel")), 
         model = factor(model, levels = model_details$model))

p_acc_mp <- list(Species = NA, 
              Taxa = NA)

p_acc_mp <- 
  lapply(p_acc_mp, function(x){
    p <- ggplot() +
      theme_bw() + 
      theme(aspect.ratio = 1, 
            panel.grid.minor = element_blank(), 
            axis.title.x = element_blank(), 
            axis.text.x = element_text(angle = 45, hjust = 1),
            panel.spacing = unit(0, "lines"),
            strip.text.y = element_blank(),legend.spacing.x = unit(0.1, "lines"), 
            legend.title = element_text(size = 10), 
            legend.position = "bottom", 
            legend.box = "horizontal", 
            legend.title.position = "top") +
      ggh4x::facet_nested(~ gene + set, scales = "free_x") + 
      labs(color = "Model") 
    p
  })

for(s in names(p_acc_mp)){
  if(s == "Taxa"){
    subdf <- 
      p_acc_mp_df |> 
      filter(str_ends(set, s))
  }else{
    subdf <- 
      p_acc_mp_df |> 
      filter(str_ends(set, "where"))
  }
  
  for(type in names(models)){
    p_acc_mp[[s]] <- 
      p_acc_mp[[s]] + 
      plot_line(df = subdf |> filter(model %in% models[[type]]), 
                x = rank, 
                y = accuracy, 
                model_type = type, 
                plot_details = model_details, 
                group = model) +
      plot_point(df = subdf |> filter(model %in% models[[type]]), 
                 x = rank, 
                 y = accuracy, 
                 model_type = type, 
                 plot_details = model_details) + 
      ggnewscale::new_scale_color() + ggnewscale::new_scale("shape") 
  }
  if(s == "Species"){
    p_acc_mp[[s]] <- 
      p_acc_mp[[s]] + 
      labs(y = "Accuracy (%)")
  }else{
    p_acc_mp[[s]] <- 
      p_acc_mp[[s]] + 
      labs(y = "Marginal recall (%)")
  }
}

p_acc_mp_joint <- 
  ggpubr::ggarrange(plotlist = p_acc_mp, 
                    ncol = 1, 
                    labels = c("a", "b"), 
                    common.legend = T, 
                    legend = "bottom")

ggsave(plot = p_acc_mp_joint, 
       filename = "plots/fourpanel_joint_mp.pdf", 
       width = 7, 
       height = 7, 
       units = "in")

# CONDITIONAL ACCURACY ---------------------------------------------------------

cond_accuracies_coi <- 
  conditional_accuracy(results_coi,
                       data_true_coi, 
                       id_list = list(All = id_all_coi,
                                      Observed = id_observed_coi,
                                      Novel = id_novel_coi))

cond_accuracies_its <- 
  conditional_accuracy(results_its,
                       data_true_its, 
                       id_list = list(All = id_all_its,
                                      Observed = id_observed_its,
                                      Novel = id_novel_its))

cond_df <- 
  bind_rows(cond_accuracies_coi |> 
              mutate(gene = "COI"), 
            cond_accuracies_its |> 
              mutate(gene = "ITS"))

no_novel <- 
  cond_df |>
  filter(set == "Novel") |>
  filter(all(cond_accuracy == 0), .by = c(model, gene)) |> 
  distinct(model) |>
  pull(model)

p_cond <- 
  cond_df |> 
  filter(set == "Novel" & !(model %in% no_novel)) |>
  mutate(model_category = 
           case_when(
             model == "Crest4" ~ "Similarity", 
             model == "BLAST threshold" ~ "Similarity",
             model == "BLAST top hit" ~ "Similarity",
             model == "dnabarcoder" ~ "Similarity",
             model == "BayesANT" ~ "Probabilistic",
             model == "PROTAX" ~ "Probabilistic",
             model == "PROTAX aug." ~ "Probabilistic", 
             str_starts(model, "EPA-") ~ "Phylogenetic placement",
             str_starts(model, "MycoAI-") ~ "Neural network",
             TRUE ~ "Composition"
           )) |> 
  mutate(rank = factor(rank, levels = allranks),
         model_category = factor(model_category,
                                 levels = c("Similarity",
                                            "Composition", 
                                            "Probabilistic", 
                                            "Phylogenetic placement", 
                                            "Neural network"))) |> 
  mutate(cond_accuracy = if_else(is.na(cond_accuracy), 0, cond_accuracy), 
         denom_cond = if_else(is.na(denom_cond), 0, denom_cond)) |>  
  ggplot() + 
  aes(x = rank,
      y = cond_accuracy * 100, 
      group = model, 
      color = model, 
      shape = model, 
      label = denom_cond) + 
  geom_line() + 
  geom_point() + 
  geom_text(nudge_y = 5, 
            size = 3, 
            show.legend = F) +
  facet_grid(gene~model_category) + 
  theme_bw() + 
  ylab("Conditional recall (%)") + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(0, "lines"), 
        legend.position = "bottom") + 
  labs(color = "Model", 
       shape = "Model") 

p_cond <- change_plot_colors(p_cond, model_details = model_details |>
                               filter(!(model %in% no_novel)))
  
ggsave(plot = p_cond,
       filename = "plots/cond_recall.pdf",
       width = 250,
       height = 130,
       units = "mm")

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
  select(-contains("kingdom")) |> 
  arrange_columns(correct_cols = correct_cols_its) |> 
  separate(Species, into = c("Species", "Similarity"), 
           sep = ";") %>%
  mutate(Similarity = as.numeric(Similarity)/100) |> 
  mutate(across(starts_with("Prob_"), ~Similarity)) |> 
  select(-Similarity) %>%
  bind_rows(data.frame(ID = data_true_its$ID[which(!(data_true_its$ID %in% .$ID))])) |> 
  arrange(ID) 

# Replace BLAST result with similarity data frame 
results_coi_coverage <- results_coi
results_coi_coverage$`BLAST top hit` <- result_blast_top_similarity_coi
results_coi_coverage$`BLAST threshold` <- results_coi_mp$`BLAST threshold`
results_coi_coverage$Crest4 <- results_coi_mp$Crest4
results_coi_coverage$DNABarcoder <- results_coi_mp$DNABarcoder

results_coi_coverage$RDP_point <-
  results_coi_mp$RDP |> 
  pivot_longer(cols = starts_with("Prob_"), 
               names_to = "rank", 
               values_to = "prob") |> 
  mutate(prob = if_else(is.na(prob), NA, 1)) |> 
  pivot_wider(names_from = rank, 
              values_from = prob)

results_coi_coverage$SINTAX_point <-
  results_coi_mp$SINTAX |> 
  pivot_longer(cols = starts_with("Prob_"), 
               names_to = "rank", 
               values_to = "prob") |> 
  mutate(prob = if_else(is.na(prob), NA, 1)) |> 
  pivot_wider(names_from = rank, 
              values_from = prob)

results_coi_coverage$IDTAXA_point <-
  results_coi_mp$IDTAXA |> 
  pivot_longer(cols = starts_with("Prob_"), 
               names_to = "rank", 
               values_to = "prob") |> 
  mutate(prob = if_else(is.na(prob), NA, 1)) |> 
  pivot_wider(names_from = rank, 
              values_from = prob)

results_its_coverage <- results_its
results_its_coverage$`BLAST top hit` <- result_blast_top_similarity_its
results_its_coverage$`BLAST threshold` <- results_its_mp$`BLAST threshold`
results_its_coverage$Crest4 <- results_its_mp$Crest4
results_its_coverage$DNABarcoder <- results_its_mp$DNABarcoder

results_its_coverage$RDP_point <-
  results_its_mp$RDP |> 
  pivot_longer(cols = starts_with("Prob_"), 
               names_to = "rank", 
               values_to = "prob") |> 
  mutate(prob = if_else(is.na(prob), NA, 1)) |> 
  pivot_wider(names_from = rank, 
              values_from = prob)

results_its_coverage$SINTAX_point <-
  results_its_mp$SINTAX |> 
  pivot_longer(cols = starts_with("Prob_"), 
               names_to = "rank", 
               values_to = "prob") |> 
  mutate(prob = if_else(is.na(prob), NA, 1)) |> 
  pivot_wider(names_from = rank, 
              values_from = prob)

results_its_coverage$IDTAXA_point <-
  results_its_mp$IDTAXA |> 
  pivot_longer(cols = starts_with("Prob_"), 
               names_to = "rank", 
               values_to = "prob") |> 
  mutate(prob = if_else(is.na(prob), NA, 1)) |> 
  pivot_wider(names_from = rank, 
              values_from = prob)

classified_correct_coi <- 
  threshold_curve(results_coi_coverage, data_true_coi) |>   
  mutate(rank = factor(rank, levels = ranks_coi))
classified_correct_its <- 
  threshold_curve(results_its_coverage, data_true_its) |>   
  mutate(rank = factor(rank, levels = ranks_its))

point_models <- c(
  "BLAST threshold",
  "Crest4", 
  "dnabarcoder", 
  "RDP_point", 
  "SINTAX_point", 
  "IDTAXA_point"
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
  facet_grid(rank~gene) + 
  theme_bw() + 
  theme(aspect.ratio = 1, 
        panel.grid.minor = element_blank(), 
        panel.spacing = unit(0, "lines"), 
        legend.position = "bottom", 
        legend.title.position = "top") + 
  labs(x = "% classified", 
       y = "% correct", 
       color = "Model") + 
  ylim(c(0, 100)) 

for(type in names(models)){
  if(any(str_detect(point_models, paste0(models[[type]], collapse = "|")))){
    p_class <- 
      p_class + 
      plot_line(df = classified_correct |> 
                  filter(model %in% models[[type]] & !model %in% point_models), 
                x = classified, 
                y = correct, 
                model_type = type, 
                plot_details = model_details, 
                shape = model) + 
      plot_point(df = classified_correct |> 
                   filter(model %in% point_models) |> 
                   mutate(model = map_chr(model, ~str_remove(.x, "_point"))) |> 
                   filter(model %in% models[[type]]), 
                 x = classified, 
                 y = correct, 
                 model_type = type, 
                 plot_details = model_details)
  }else{
    p_class <- 
      p_class + 
      plot_line(df = classified_correct |> 
                  filter(model %in% models[[type]] & !model %in% point_models), 
                x = classified, 
                y = correct, 
                model_type = type, 
                plot_details = model_details, 
                shape = model)
  }
  p_class <- 
    p_class + 
    ggnewscale::new_scale_color() + ggnewscale::new_scale("shape")
}

p_class

ggsave(plot = p_class, 
       filename = "plots/classified_correct_joint.pdf", 
       width = 7.2, 
       height = 7.2, 
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


p_error_df <- 
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
  filter(!is.na(error_rate)) 

p_error <- 
  ggplot() + 
  facet_grid(gene~error_type, scales = "free_x") +
  theme_bw() + 
  ylab("Error rate (%)") + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1), 
        aspect.ratio = 1, 
        panel.spacing = unit(0, "lines"), 
        legend.position = "bottom", 
        legend.title.position = "top")

for(type in names(models)){
  p_error <- 
    p_error + 
    plot_line(df = p_error_df |> filter(model %in% models[[type]]), 
              x = rank, 
              y = error_rate, 
              model_type = type, 
              plot_details = model_details, 
              group = model) + 
    plot_point(df = p_error_df |> filter(model %in% models[[type]]), 
               x = rank, 
               y = error_rate, 
               model_type = type, 
               plot_details = model_details) + 
    ggnewscale::new_scale_color() + ggnewscale::new_scale("shape")
}

p_error

ggsave(plot = p_error, 
       filename = "plots/errorrates_joint.pdf", 
       height = 7, 
       width = 7.5, 
       units = "in")

# TOTAL NUMBER OF NOVELTY PREDICTIONS ------------------------------------------

pred_novel_coi <- count_novel(results_coi, data_true_coi, id_novel_coi)
pred_novel_its <- count_novel(results_its, data_true_its, id_novel_its)

pred_novel_df <- 
  bind_rows(pred_novel_coi |> 
              mutate(gene = "COI"), 
            pred_novel_its |> 
              mutate(gene = "ITS")) |> 
  mutate(rank = factor(rank, levels = allranks))

p_tot_novel <-
  pred_novel_df |> 
  filter(!all(pred_novel == 0), .by = model) |> 
  ggplot() + 
  geom_line(data = pred_novel_df |> filter(model == model[1]),
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
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  facet_wrap(~gene,
             scales = "free_x")
p_tot_novel <- change_plot_colors(p_tot_novel, model_details)

ggsave(plot = p_tot_novel, 
       filename = "plots/sum_novel.pdf", 
       height = 100, 
       width = 200, 
       units = "mm")

# TESTSHORT DATA ---------------------------------------------------------------

# Read short data 
results_coi_short <- readRDS("coi/result_list_short_nt.rds")
results_its_short <- readRDS("its/result_list_short.rds")

# Remove Kingdom from ITS data 
results_its_short <- lapply(results_its_short, function(df) df[,correct_cols_its])

# Get calibrations data 
calibrations_coi_short <-  
  get_calibration(results_coi_short, 
                  data_true_coi, 
                  observed_everywhere_coi) 
calibrations_its_short <- 
  get_calibration(results_its_short,
                  data_true_its, 
                  observed_everywhere_its) 

acc_long_short <- 
  get_accuracy_from_cal(calibrations_coi) |> 
  left_join(get_accuracy_from_cal(calibrations_coi_short), 
            by = c("model", 
                   "rank", 
                   "set")) |> 
  rename(accuracy_long = accuracy.x, 
         accuracy_short = accuracy.y) |> 
  mutate(gene = "COI") |> 
  bind_rows(get_accuracy_from_cal(calibrations_its) |> 
              left_join(get_accuracy_from_cal(calibrations_its_short), 
                        by = c("model", 
                               "rank", 
                               "set")) |> 
              rename(accuracy_long = accuracy.x, 
                     accuracy_short = accuracy.y) |> 
              mutate(gene = "ITS")) |> 
  mutate(diff = accuracy_short- accuracy_long,
         rank = factor(rank, levels = allranks), 
         set = factor(set, levels = c("All", "Observed", "Novel"))) |> 
  filter(set != "All")

p_acc_short <- 
  ggplot() + 
  facet_grid(set~gene, 
             scales = "free") + 
  geom_hline(yintercept = 0, 
             linetype = "dashed") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.spacing = unit(0, "lines"), 
        axis.title.x = element_blank(), 
        legend.position = "bottom", 
        legend.title.position = "top") + 
  labs(y = "Accuracy difference (%pt.)")

for(type in names(models)){
  p_acc_short <- 
    p_acc_short + 
    plot_line(df = acc_long_short |> filter(model %in% models[[type]]), 
              x = rank, 
              y = diff, 
              model_type = type, 
              plot_details = model_details, 
              group = model) + 
    plot_point(df = acc_long_short |> filter(model %in% models[[type]]), 
               x = rank, 
               y = diff, 
               model_type = type, 
               plot_details = model_details) + 
    ggnewscale::new_scale_color() + ggnewscale::new_scale("shape")
}

p_acc_short

ggsave(plot = p_acc_short, 
       filename = "plots/accuracy_diff_short.pdf", 
       height = 7, 
       width = 7, 
       units = "in")

# TEXT INFO --------------------------------------------------------------------

# Accuracy 
acc_df <- 
  accuracies_coi |> 
  mutate(gene = "COI") |> 
  bind_rows(accuracies_its |> 
              mutate(gene = "ITS")) 
# at genus and species level, observed species   
acc_df |> 
  filter(set == "Observed", rank %in% c("Genus", "Species")) |> 
  pivot_wider(names_from = rank, 
              values_from = accuracy) |> 
  arrange(Genus) 

# novel, genus
acc_df |> 
  filter(set == "Novel", rank == "Genus") |> 
  pivot_wider(names_from = gene, values_from = accuracy) |> 
  arrange(COI)

# novel, species 
acc_df |> 
  filter(set == "Novel", rank == "Species", str_detect(model, "EPA-ng|PROTAX|BayesANT")) |> 
  pivot_wider(names_from = gene, values_from = accuracy) |> 
  arrange(COI)

# misclassification 
misclass_df_coi |> 
  mutate(gene = "COI") |> 
  bind_rows(misclass_df_its |> 
              mutate(gene = "ITS")) |>
  pivot_wider(names_from = gene, values_from = misclassification) |> 
  arrange(-COI) |> View()

# overclassification, 
oclass_df_coi |> 
  mutate(gene = "COI") |> 
  bind_rows(oclass_df_its |> 
              mutate(gene = "ITS")) |> 
  pivot_wider(names_from = gene, values_from = overclassification) |> 
  filter(model == "RDP" | model == "SINTAX")

