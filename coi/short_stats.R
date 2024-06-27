# Short and snappy statistics about the raw data 

rm(list = ls())
setwd("coi")

# LOAD FUNCTIONS ---------------------------------------------------------------

library(tidyverse)
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

ranks <- colnames(data_true |> select(-ID))

# SEQUENCE COUNT ---------------------------------------------------------------

print(
  tibble(dataset = c("Train", "Test"), 
       sequences = sapply(data, nrow), 
       unique_sequences = sapply(data, function(x) length(unique(x$DNA))))
  )

# TAXA COUNT -------------------------------------------------------------------

bind_rows(data$train |> mutate(dataset = "Train"), 
          data$test |> mutate(dataset = "Test")) |> 
  select(-DNA, -ID) |> 
  relocate(dataset, .before = "Class") %>% 
  pivot_longer(2:ncol(.), names_to = "rank", values_to = "taxon") |> 
  distinct() |> 
  group_by(dataset, rank) |> 
  summarise(count = n()) |> 
  pivot_wider(values_from = count, names_from = rank) |> 
  select(dataset, all_of(ranks)) |> 
  arrange(desc(dataset))

# NOVEL ANYWHERE ---------------------------------------------------------------

tibble(novel_count = sum(str_ends(data_true$Species, "_new")), 
       novel_perc = sum(str_ends(data_true$Species, "_new"))/nrow(data_true))

# RANK-SPECIFIC NOVELTY --------------------------------------------------------

data_true %>%
  pivot_longer(2:ncol(.), 
               names_to = "rank", 
               values_to = "taxon") |> 
  mutate(new = map2_lgl(rank, taxon, ~str_ends(.y, paste0(.x, "_new")))) |> 
  group_by(rank) |> 
  summarise(novel_count = sum(new)) |> 
  mutate(novel_perc = novel_count/nrow(data_true)) |> 
  arrange(match(rank, ranks))

# TAXONOMIC OVERLAP ------------------------------------------------------------

data$train |> 
  select(all_of(ranks)) |> 
  mutate(train = T) |> 
  relocate(train, .before = "Class") %>%
  pivot_longer(2:ncol(.), 
               names_to = "rank", 
               values_to = "taxon") |> 
  distinct() |> 
  full_join(data$test |> 
              select(all_of(ranks)) |> 
              mutate(test = T) |>
              relocate(test, .before = "Class") %>%
              pivot_longer(2:ncol(.), 
                           names_to = "rank", 
                           values_to = "taxon") |> 
              distinct()) |> 
  mutate(set = if_else(!is.na(test) & !is.na(train), "both", 
                       dplyr::if_else(is.na(test), "train", "test"))) |> 
  group_by(rank, set) |> 
  summarise(count = n()) |> 
  pivot_wider(names_from = set, values_from = count) |> 
  mutate(train = map_dbl(train, ~if_else(is.na(.x), 0, .x))) |> 
  mutate(total_count = test+train+both) |> 
  mutate(perc_common = both/total_count * 100) |> 
  arrange(match(rank, ranks))


