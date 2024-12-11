
separate_labels <- function(df, taxgroups = c("Class", 
                                              "Order", 
                                              "Family", 
                                              "Subfamily", 
                                              "Tribe", 
                                              "Genus", 
                                              "Species")){
  
  return(df |> 
    separate(label, into = c("ID", "label"), sep = " ") |> 
    separate(label, into = taxgroups, sep = "\\|"))
}


# Alluvial plot of overlap between FinBOL and GBOL 

#devtools::install_github("davidsjoberg/ggsankey")

library(ggsankey)
library(tidyverse)

setwd("coi")

# READ DATA ----------------

gbol_raw <- ape::read.FASTA("data/test_nt_aln_label.fasta")
finbol_raw <- ape::read.FASTA("data/train_nt_aln_label.fasta")

# LABELS --------------------

common_lab <- "b_common"
gbol_lab <- "a_uniq_gbol"
finbol_lab <- "c_uniq_finbol"

# TIBBLES OF SEQUENCES -----------------------

finbol_seq <- tibble(seq = sapply(as.character(finbol_raw),
                                  function(x) paste(x, collapse = "")))
finbol_seq <- bind_cols(finbol_seq, tibble(label = names(finbol_raw))) |> 
  separate_labels() |> 
  select(seq, ID)

gbol_seq <- tibble(seq = sapply(as.character(gbol_raw),
                                function(x) paste(x, collapse = "")))
gbol_seq <- bind_cols(gbol_seq, tibble(label = names(gbol_raw))) |> 
  separate_labels() |> 
  select(seq, ID)

# SEPARATE LABELS --------------------

gbol_taxa <- tibble(label = names(gbol_raw)) |> 
  separate_labels() |> 
  left_join(gbol_seq)

finbol_taxa <- tibble(label = names(finbol_raw)) |> 
  separate_labels() |> 
  left_join(finbol_seq)

# FIND COMMON TAXA AND SEQUENCES -------------------

common <- 
  gbol_taxa %>%
  pivot_longer(2:ncol(.), names_to = "rank", values_to = "taxon") |> 
  distinct(rank, taxon) |> 
  inner_join(finbol_taxa %>%
               pivot_longer(2:ncol(.), names_to = "rank", values_to = "taxon") |> 
               distinct(rank, taxon))

# DIVIDE DATASETS INTO UNIQUE AND SHARED ---------------------

sets <- 
  # gbol shared 
  gbol_taxa %>% 
  pivot_longer(2:ncol(.), 
               names_to = "rank", 
               values_to = "taxon") |> 
  inner_join(common, by = c("rank", "taxon")) |> 
  mutate(set = common_lab) |> 
  # gbol unique 
  bind_rows(gbol_taxa %>% 
              pivot_longer(2:ncol(.), 
                           names_to = "rank", 
                           values_to = "taxon") |> 
              anti_join(common, by = c("rank", "taxon")) |> 
              mutate(set = gbol_lab)) |> 
  # finbol shared 
  bind_rows(finbol_taxa %>%
              pivot_longer(2:ncol(.), 
                           names_to = "rank", 
                           values_to = "taxon") |> 
              inner_join(common, by = c("rank", "taxon")) |> 
              mutate(set = common_lab)) |> 
  # finbol unique 
  bind_rows(finbol_taxa %>%
              pivot_longer(2:ncol(.), 
                           names_to = "rank", 
                           values_to = "taxon") |> 
              anti_join(common, by = c("rank", "taxon")) |> 
              mutate(set = finbol_lab)) |> 
  select(-taxon) |> 
  pivot_wider(names_from = rank, values_from = set) |> 
  rename(Sequence = seq)

# COUNT OF CLASSES --------------------------

counts <- 
  common |>
  group_by(rank) |> 
  summarise(count = n()) |> 
  mutate(node = common_lab) |> 
  bind_rows(gbol_taxa %>% 
              pivot_longer(2:ncol(.), names_to = "rank", values_to = "taxon") |> 
              anti_join(common) |> 
              distinct(rank, taxon) |> 
              group_by(rank) |> 
              summarise(count = n()) |> 
              mutate(node = gbol_lab)) |> 
  bind_rows(finbol_taxa %>% 
              pivot_longer(2:ncol(.), names_to = "rank", values_to = "taxon") |> 
              anti_join(common) |> 
              distinct(rank, taxon) |> 
              group_by(rank) |> 
              summarise(count = n()) |> 
              mutate(node = finbol_lab)) |> 
  rename(x = rank) |> 
  mutate(x = if_else(x == "seq", "Sequence", x)) |> 
  mutate(x = factor(x)) 
  

# SHAPE FOR PLOT ---------------------

sets_plot <- 
  sets |> 
  make_long(Class, Order, Family, Subfamily, Tribe, Genus, Species, Sequence) |> 
  left_join(counts)

# PLOT ----------------------------

p <- 
  sets_plot |> 
  ggplot() +
  aes(x = x, 
      next_x = next_x, 
      node = node, 
      next_node = next_node, 
      fill = node, 
      label = count) + 
  geom_sankey(flow.alpha = 0.5, width = 0.35) + 
  geom_sankey_label() + 
  theme_minimal() + 
  scale_y_continuous(breaks = c(-40000, -20000, 0, 20000, 40000),
                     labels = c("0", "20000", "40000", "60000", "80000"), 
                     position = "right") + 
  scale_fill_discrete(limits = c(finbol_lab, common_lab, gbol_lab), 
                      labels = c("Unique to FinBOL", "Shared", "Unique to GBOL")) + 
  ylab("Number of sequences") + 
  theme(axis.title.x = element_blank(), 
        panel.grid.major.x = element_blank(), 
        legend.title = element_blank(), 
        axis.text = element_text(size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.text = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        legend.position = "top") 
p

ggsave(plot = p, 
       filename = "../plots/overlap.pdf", 
       width = 1900, 
       height = 1600, 
       units = "px", 
       dpi = 300)
