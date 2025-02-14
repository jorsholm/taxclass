# Alluvial plot of overlap between train and test sets 

#devtools::install_github("davidsjoberg/ggsankey")
library(ggsankey)
library(tidyverse)

cases <- c("coi", "its")

p <- list()

# LABELS --------------------

common_lab <- "b_common"
test_lab <- "a_uniq_test"
train_lab <- "c_uniq_train"

for(case in cases){
  
  # READ DATA ----------------
  test_raw <- ape::read.FASTA(paste0(case, "/data/test_nt.fasta"))
  train_raw <- ape::read.FASTA(paste0(case, "/data/train_nt.fasta"))
  
  if(case == "coi"){
    correctcols <- c("Id", "Class", "Order", "Family", "Subfamily", "Tribe", "Genus", "Species", "Sequence")
  }else{
    correctcols <- c("Id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Sequence")
  }
  
  # TIBBLES OF SEQUENCES AND TAXONOMY -----------------------
  
  # Train data 
  train_seq <- tibble(Sequence = sapply(as.character(train_raw),
                                        function(x) paste(x, collapse = "")))
  
  train_taxa <- 
    bind_cols(train_seq, tibble(id = names(train_raw))) |> 
    left_join(read_delim(paste0(case, "/data/train_tax.tsv"),
                         delim = "\t"))
  train_taxa <- train_taxa[,c(2:ncol(train_taxa), 1)]
  colnames(train_taxa) <- correctcols
  
  # Test data 
  test_seq <- tibble(Sequence = sapply(as.character(test_raw),
                                       function(x) paste(x, collapse = "")))
  
  test_taxa <- 
    bind_cols(test_seq, tibble(id = names(test_raw))) |> 
    left_join(read_delim(paste0(case, "/data/test_tax.tsv"),
                         delim = "\t"))
  test_taxa <- test_taxa[,c(2:ncol(test_taxa), 1)]
  colnames(test_taxa) <- correctcols
  
  # FIND COMMON TAXA AND SEQUENCES -------------------
  
  common <- 
    test_taxa %>%
    pivot_longer(2:ncol(.), names_to = "rank", values_to = "taxon") |> 
    distinct(rank, taxon) |> 
    inner_join(train_taxa %>%
                 pivot_longer(2:ncol(.), names_to = "rank", values_to = "taxon") |> 
                 distinct(rank, taxon))
  
  # DIVIDE DATASETS INTO UNIQUE AND SHARED ---------------------
  
  sets <- 
    # test shared 
    test_taxa %>% 
    pivot_longer(2:ncol(.), 
                 names_to = "rank", 
                 values_to = "taxon") |> 
    inner_join(common, by = c("rank", "taxon")) |> 
    mutate(set = common_lab) |> 
    # test unique 
    bind_rows(test_taxa %>% 
                pivot_longer(2:ncol(.), 
                             names_to = "rank", 
                             values_to = "taxon") |> 
                anti_join(common, by = c("rank", "taxon")) |> 
                mutate(set = test_lab)) |> 
    # train shared 
    bind_rows(train_taxa %>%
                pivot_longer(2:ncol(.), 
                             names_to = "rank", 
                             values_to = "taxon") |> 
                inner_join(common, by = c("rank", "taxon")) |> 
                mutate(set = common_lab)) |> 
    # train unique 
    bind_rows(train_taxa %>%
                pivot_longer(2:ncol(.), 
                             names_to = "rank", 
                             values_to = "taxon") |> 
                anti_join(common, by = c("rank", "taxon")) |> 
                mutate(set = train_lab)) |> 
    select(-taxon) |> 
    pivot_wider(names_from = rank, values_from = set) 
  
  # COUNT OF CLASSES --------------------------
  
  counts <- 
    common |>
    group_by(rank) |> 
    summarise(count = n()) |> 
    mutate(node = common_lab) |> 
    bind_rows(test_taxa %>% 
                pivot_longer(2:ncol(.), names_to = "rank", values_to = "taxon") |> 
                anti_join(common) |> 
                distinct(rank, taxon) |> 
                group_by(rank) |> 
                summarise(count = n()) |> 
                mutate(node = test_lab)) |> 
    bind_rows(train_taxa %>% 
                pivot_longer(2:ncol(.), names_to = "rank", values_to = "taxon") |> 
                anti_join(common) |> 
                distinct(rank, taxon) |> 
                group_by(rank) |> 
                summarise(count = n()) |> 
                mutate(node = train_lab)) |> 
    rename(x = rank) |> 
    mutate(x = factor(x)) 
  
  # SHAPE FOR PLOT ---------------------
  
  sets_plot <- 
    sets |> 
    make_long(all_of(correctcols[-1])) |> 
    left_join(counts)
  
  # PLOT ----------------------------
  
  p[[case]] <- 
    sets_plot |> 
    ggplot() +
    aes(x = x, 
        next_x = next_x, 
        node = node, 
        next_node = next_node, 
        fill = node, 
        label = count) + 
    geom_sankey(flow.alpha = 0.5, width = 0.35) + 
    geom_sankey_label(show.legend = F) + 
    theme_minimal() + 
    scale_fill_discrete(limits = c(train_lab, common_lab, test_lab), 
                        labels = c("Unique to train", "Shared", "Unique to test")) + 
    ylab("Number of sequences") + 
    theme(axis.title.x = element_blank(), 
          panel.grid.major.x = element_blank(), 
          legend.title = element_blank(), 
          axis.text = element_text(size = 11),
          axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.text = element_text(size = 11),
          axis.title.y = element_text(size = 11),
          legend.position = "top") 
  
  if(case == "its"){
    p[[case]] <- p[[case]] + 
      scale_y_continuous(breaks = c(-20000, 0, 20000),
                         labels = c("0", "20000", "40000"), 
                         position = "right")
  }else{
    p[[case]] <- p[[case]] + 
      scale_y_continuous(breaks = c(-40000, -20000, 0, 20000, 40000),
                         labels = c("0", "20000", "40000", "60000", "80000"), 
                         position = "right")
  }
  
  ggsave(plot = p[[case]], 
         filename = paste0("plots/overlap_", case, ".pdf"), 
         width = 1900, 
         height = 1600, 
         units = "px", 
         dpi = 300)
  
}

overlap_plot <- 
  ggpubr::ggarrange(plotlist = p,
                    ncol = 1,
                    common.legend = T,
                    legend = "bottom", 
                    labels = c("COI", "ITS"))

ggsave(plot = overlap_plot, 
       filename = "plots/overlap.pdf", 
       width = 1900, 
       height = 2000, 
       units = "px", 
       dpi = 300)
