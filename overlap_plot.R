# Alluvial plot of overlap between train and test sets 

#devtools::install_github("davidsjoberg/ggsankey")
library(ggsankey)
library(tidyverse)

# LABELS --------------------

common_lab <- "b_common"
test_lab <- "a_uniq_test"
train_lab <- "c_uniq_train"

# READ DATA ----------------
test_raw_coi <- ape::read.FASTA("coi/data/test_nt.fasta")
train_raw_coi <- ape::read.FASTA("coi/data/train_nt.fasta")
test_raw_its <- ape::read.FASTA("its/data/test_nt.fasta")
train_raw_its <- ape::read.FASTA("its/data/train_nt.fasta")
  
correctcols_coi <- c("Id", "Class", "Order", "Family", "Subfamily", "Tribe", "Genus", "Species", "Sequence")
correctcols_its <- c("Id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Sequence")

# TIBBLES OF SEQUENCES AND TAXONOMY -----------------------
  
  # Train data 
train_seq_coi <- tibble(Sequence = sapply(as.character(train_raw_coi),
                                          function(x) paste(x, collapse = "")))
train_seq_its <- tibble(Sequence = sapply(as.character(train_raw_its),
                                          function(x) paste(x, collapse = "")))
  
  train_taxa_coi <- 
    bind_cols(train_seq_coi, tibble(id = names(train_raw_coi))) |> 
    left_join(read_delim("coi/data/train_tax.tsv",
                         delim = "\t"))
  train_taxa_coi <- train_taxa_coi[,c(2:ncol(train_taxa_coi), 1)]
  colnames(train_taxa_coi) <- correctcols_coi
  
  train_taxa_its <- 
    bind_cols(train_seq_its, tibble(id = names(train_raw_its))) |> 
    left_join(read_delim("its/data/train_tax.tsv",
                         delim = "\t"))
  train_taxa_its <- train_taxa_its[,c(2:ncol(train_taxa_its), 1)]
  colnames(train_taxa_its) <- correctcols_its
  
  # Test data 
  test_seq_coi <- tibble(Sequence = sapply(as.character(test_raw_coi),
                                       function(x) paste(x, collapse = "")))
  test_seq_its <- tibble(Sequence = sapply(as.character(test_raw_its),
                                           function(x) paste(x, collapse = "")))
  
  test_taxa_coi <- 
    bind_cols(test_seq_coi, tibble(id = names(test_raw_coi))) |> 
    left_join(read_delim("coi/data/test_tax.tsv",
                         delim = "\t"))
  test_taxa_coi <- test_taxa_coi[,c(2:ncol(test_taxa_coi), 1)]
  colnames(test_taxa_coi) <- correctcols_coi
  
  test_taxa_its <- 
    bind_cols(test_seq_its, tibble(id = names(test_raw_its))) |> 
    left_join(read_delim("its/data/test_tax.tsv",
                         delim = "\t"))
  test_taxa_its <- test_taxa_its[,c(2:ncol(test_taxa_its), 1)]
  colnames(test_taxa_its) <- correctcols_its
  
  # FIND COMMON TAXA AND SEQUENCES -------------------
  
  common_coi <- 
    test_taxa_coi %>%
    pivot_longer(2:ncol(.), names_to = "rank", values_to = "taxon") |> 
    distinct(rank, taxon) |> 
    inner_join(train_taxa_coi %>%
                 pivot_longer(2:ncol(.), names_to = "rank", values_to = "taxon") |> 
                 distinct(rank, taxon))
  common_its <- 
    test_taxa_its %>%
    pivot_longer(2:ncol(.), names_to = "rank", values_to = "taxon") |> 
    distinct(rank, taxon) |> 
    inner_join(train_taxa_its %>%
                 pivot_longer(2:ncol(.), names_to = "rank", values_to = "taxon") |> 
                 distinct(rank, taxon))
  
  # DIVIDE DATASETS INTO UNIQUE AND SHARED ---------------------
  
  sets_coi <- 
    # test shared 
    test_taxa_coi %>% 
    pivot_longer(2:ncol(.), 
                 names_to = "rank", 
                 values_to = "taxon") |> 
    inner_join(common_coi, by = c("rank", "taxon")) |> 
    mutate(set = common_lab) |> 
    # test unique 
    bind_rows(test_taxa_coi %>% 
                pivot_longer(2:ncol(.), 
                             names_to = "rank", 
                             values_to = "taxon") |> 
                anti_join(common_coi, by = c("rank", "taxon")) |> 
                mutate(set = test_lab)) |> 
    # train shared 
    bind_rows(train_taxa_coi %>%
                pivot_longer(2:ncol(.), 
                             names_to = "rank", 
                             values_to = "taxon") |> 
                inner_join(common_coi, by = c("rank", "taxon")) |> 
                mutate(set = common_lab)) |> 
    # train unique 
    bind_rows(train_taxa_coi %>%
                pivot_longer(2:ncol(.), 
                             names_to = "rank", 
                             values_to = "taxon") |> 
                anti_join(common_coi, by = c("rank", "taxon")) |> 
                mutate(set = train_lab)) |> 
    select(-taxon) |> 
    pivot_wider(names_from = rank, values_from = set) 
  
  sets_its <- 
    # test shared 
    test_taxa_its %>% 
    pivot_longer(2:ncol(.), 
                 names_to = "rank", 
                 values_to = "taxon") |> 
    inner_join(common_its, by = c("rank", "taxon")) |> 
    mutate(set = common_lab) |> 
    # test unique 
    bind_rows(test_taxa_its %>% 
                pivot_longer(2:ncol(.), 
                             names_to = "rank", 
                             values_to = "taxon") |> 
                anti_join(common_its, by = c("rank", "taxon")) |> 
                mutate(set = test_lab)) |> 
    # train shared 
    bind_rows(train_taxa_its %>%
                pivot_longer(2:ncol(.), 
                             names_to = "rank", 
                             values_to = "taxon") |> 
                inner_join(common_its, by = c("rank", "taxon")) |> 
                mutate(set = common_lab)) |> 
    # train unique 
    bind_rows(train_taxa_its %>%
                pivot_longer(2:ncol(.), 
                             names_to = "rank", 
                             values_to = "taxon") |> 
                anti_join(common_its, by = c("rank", "taxon")) |> 
                mutate(set = train_lab)) |> 
    select(-taxon) |> 
    pivot_wider(names_from = rank, values_from = set) 
  
  # COUNT OF CLASSES --------------------------
  
  counts_coi <- 
    common_coi |>
    group_by(rank) |> 
    summarise(count = n()) |> 
    mutate(node = common_lab) |> 
    bind_rows(test_taxa_coi %>% 
                pivot_longer(2:ncol(.), names_to = "rank", values_to = "taxon") |> 
                anti_join(common_coi) |> 
                distinct(rank, taxon) |> 
                group_by(rank) |> 
                summarise(count = n()) |> 
                mutate(node = test_lab)) |> 
    bind_rows(train_taxa_coi %>% 
                pivot_longer(2:ncol(.), names_to = "rank", values_to = "taxon") |> 
                anti_join(common_coi) |> 
                distinct(rank, taxon) |> 
                group_by(rank) |> 
                summarise(count = n()) |> 
                mutate(node = train_lab)) |> 
    rename(x = rank) |> 
    mutate(x = factor(x)) 
  
  counts_its <- 
    common_its |>
    group_by(rank) |> 
    summarise(count = n()) |> 
    mutate(node = common_lab) |> 
    bind_rows(test_taxa_its %>% 
                pivot_longer(2:ncol(.), names_to = "rank", values_to = "taxon") |> 
                anti_join(common_its) |> 
                distinct(rank, taxon) |> 
                group_by(rank) |> 
                summarise(count = n()) |> 
                mutate(node = test_lab)) |> 
    bind_rows(train_taxa_its %>% 
                pivot_longer(2:ncol(.), names_to = "rank", values_to = "taxon") |> 
                anti_join(common_its) |> 
                distinct(rank, taxon) |> 
                group_by(rank) |> 
                summarise(count = n()) |> 
                mutate(node = train_lab)) |> 
    rename(x = rank) |> 
    mutate(x = factor(x))
  
  # SHAPE FOR PLOT ---------------------
  
 sets_plot <- 
    sets_coi |> 
    rename(rank1 = Class, 
           rank2 = Order, 
           rank3 = Family, 
           rank4 = Subfamily, 
           rank5 = Tribe, 
           rank6 = Genus, 
           rank7 = Species, 
           rank8 = Sequence) |> 
    make_long(all_of(paste0("rank", 1:8))) |> 
    left_join(counts_coi |> 
                mutate(x = if_else(x == "Class", "rank1", 
                                   if_else(x == "Order", "rank2", 
                                           if_else(x == "Family", "rank3", 
                                                   if_else(x == "Subfamily", "rank4", 
                                                           if_else(x == "Tribe", "rank5", 
                                                                   if_else(x == "Genus", "rank6",
                                                                           if_else(x == "Species", "rank7", "rank8"))))))))) |> 
    mutate(set = "COI") |> 
    bind_rows(sets_its |> 
                rename(rank1 = Kingdom, 
                       rank2 = Phylum, 
                       rank3 = Class, 
                       rank4 = Order, 
                       rank5 = Family,
                       rank6 = Genus, 
                       rank7 = Species, 
                       rank8 = Sequence) |> 
                make_long(all_of(paste0("rank", 1:8))) |> 
                left_join(counts_its |> 
                            mutate(x = if_else(x == "Kingdom", "rank1", 
                                               if_else(x == "Phylum", "rank2", 
                                                       if_else(x == "Class", "rank3", 
                                                               if_else(x == "Order", "rank4", 
                                                                       if_else(x == "Family", "rank5", 
                                                                               if_else(x == "Genus", "rank6",
                                                                                       if_else(x == "Species", "rank7", "rank8"))))))))) |> 
                mutate(set = "ITS")) |> 
    mutate(x = factor(x, levels = paste0("rank", 1:8)), 
           next_x = factor(next_x, levels = paste0("rank", 1:8)))
  
  # PLOT ----------------------------
  
  # allranks <- c("Kingdom", 
  #               "Phylum", 
  #               "Class", 
  #               "Order", 
  #               "Family", 
  #               "Subfamily", 
  #               "Tribe", 
  #               "Genus", 
  #               "Species", 
  #               "Sequence")
  
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
    geom_sankey_label(show.legend = F, size = 3) + 
    theme_minimal() + 
    scale_fill_manual(limits = c(train_lab, common_lab, test_lab), 
                        labels = c("Unique to train", "Shared", "Unique to test"), 
                        values = c("#faab5c", "#40a373", "#2a94d1")) + 
    ylab("Number of sequences") + 
    facet_wrap(~set) + 
    theme(axis.title.x = element_blank(), 
          panel.grid.major.x = element_blank(), 
          legend.title = element_blank(), 
          axis.text = element_text(size = 11),
          axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.text = element_text(size = 11),
          axis.title.y = element_text(size = 11),
          legend.position = "top",
          panel.spacing = unit(0, "lines")) + 
    scale_y_continuous(breaks = c(-40000, -20000, 0, 20000, 40000),
                       labels = c("0", "20000", "40000", "60000", "80000"))
    
p
    
ggsave(plot = p, 
       filename = "plots/overlap.pdf", 
       width = 2500, 
       height = 1500, 
       units = "px", 
       dpi = 300)
