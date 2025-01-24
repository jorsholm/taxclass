################################################################################
# Functions for analysing performance of taxonomic classification algorithms.
################################################################################

# FOR PLOTTING -----------------------------------------------------------------

algs_sorted <- c(
  # Similarity
  "BLAST top hit", 
  "BLAST threshold",
  "DNABarcoder",
  "Crest4", 
  # K-mer 
  "IDTAXA",
  "RDP", 
  "SINTAX",
  # Probabilistic
  "BayesANT", 
  #"PROTAX", 
  #Neural networks
  "MycoAI-BERT", 
  "MycoAI-CNN", 
  # phylogenetic 
  "EPA-ng phyltree", 
  "EPA-ng taxtree")

plot_colors <- c(
  # Similarity
  "#bdd7e7",
  "#6baed6",
  "#3182bd",
  "#08519c", 
  # K-mer
  "#bae4b3",
  "#74c476",
  "#238b45",
  # Probabilistic
  "#bcbddc",
  #"#756bb1",
  #Neural networks
  "#fdbe85",
  "#fd8d3c",
  # phylogenetic
  "#fbb4b9",
  "#f768a1")

point_shapes <- c( 
  #similarity 
  20,
  2,
  4,
  5, 
  # k-mer
  20,
  4,
  5,
  # probabilistic 
  2,
  # 5,
  # neural
  20,
  4,
  # phylogenetic
  2,
  5
)

change_plot_colors <- function(x){
  return(x + 
    scale_color_manual(name = "Model",
                     labels = algs_sorted, 
                     values = plot_colors, 
                     breaks = algs_sorted) + 
    scale_shape_manual(name = "Model",
                       labels = algs_sorted,
                       values = point_shapes, 
                       breaks = algs_sorted)
  )
} 

# OTHER FUNCTIONS --------------------------------------------------------------

#' Loads FinBOL data (train) and GBOL data (test) as listed dataframes. 
#' @param level lowest taxonomic level (genus/species)
load_train_test <- function(train_file = "data/train_tax.tsv",
                             test_file = "data/test_tax.tsv", 
                             ranks = c("Class", 
                                       "Order",
                                       "Family",
                                       "Subfamily",
                                       "Tribe",
                                       "Genus",
                                       "Species")){
  
  data_train <- read.delim(train_file)
  data_test <- read.delim(test_file)
  
  colnames(data_train) <- colnames(data_test) <- c(colnames(data_train)[1], ranks)
  
  return(list("train" = data_train, "test"= data_test))
}

#' Replaces labels of classes that are unique to test dataset. 
#' Example of new labels: Araneae_Family_new
#' Or, if first level is new: Class_new
#' @param test 
#' @param train 
get_data_true <- function(test, train, 
                          ranks = c("Class", 
                                    "Order",
                                    "Family",
                                    "Subfamily",
                                    "Tribe",
                                    "Genus",
                                    "Species")){

  # Keep only ranks columns 
  train <- train[,ranks]
  data_true <- test[,ranks]

  # Get level names
  n_ranks <- length(ranks)

  for(l in 1:n_ranks){
    # Find taxa that are unique to test
    id <- !(data_true[,l] %in% unique(train[,l]))
    n_new <- sum(id)

    # Replace new classes with "_new"-name
    if(n_new > 0){
      if(l == 1){
        name_clust <- paste0(ranks[1], "_new")
        data_true[id, l:n_ranks] <- matrix(rep(name_clust, n_new*n_ranks),
                                            nrow = n_new, ncol = n_ranks)
      }else{

        for(i in which(id == TRUE)){
          if(!grepl("_new$", data_true[i, l-1])){
            name_clust <- paste0(data_true[i, l-1], "_", ranks[l], "_new")
            data_true[i, l:n_ranks] <- rep(name_clust, n_ranks-l + 1)
          }
        }
      }
    }
  }
  return(cbind("ID" = test[,1], data_true))
}

#' Renames PROTAX new taxa ("unk") to BayesANT standard. 
#' Only deals with taxonomy-columns (not prob-columns). 
#' @param x one row data frame with taxonomy columns of protax output 
rename_unk <- function(x, unktxt = "unk"){
  # Remove NAs 
  cols <- names(x)
  na_cols <- names(x[which(is.na(x))])
  if(length(na_cols) > 0) x <- x[-which(names(x) %in% na_cols)]
  
  is_unk <- x == unktxt
  
  if (sum(is_unk) > 0) {
    if (all(x == unktxt)) {
      x[1:length(x)] <- "Class_new"
    } else{
      d <- min(which(is_unk))
      x[d:length(x)] <- paste0(x[d - 1], "_", stringr::str_to_title(names(x)[d]), "_new")
    }
  }
  # Replace NAs
  x[na_cols] <- NA
  # Return sorted to original order 
  return(x[cols])
}

#' Renames 'unk' to BayesANT standard + sort columns 
#' @param out data frame with PROTAX results. Cannot contain NA values  
rename_unk_output <- function(out, unktxt = "unk"){
  # Rename unk 
  renamed <- cbind("ID" = out[,which(colnames(out) == "ID")],
                   t(apply(out[, which(!stringr::str_detect(colnames(out), "Prob|ID"))],
                           1, function(x) rename_unk(x, unktxt = unktxt))), 
                   out[,which(stringr::str_detect(colnames(out), "Prob_"))])
  
  return(renamed)
}

#' Arrange and rename columns in a specified order
#' @param df The dataframe to rearrange 
#' @param correct_cols The order and name of columns wanted
arrange_columns <- function(df, correct_cols){
  if(!all(sort(stringr::str_to_lower(colnames(df))) ==
          sort(stringr::str_to_lower(correct_cols)))) print("Mismatching columns.")
  
  df <- df[, match(stringr::str_to_lower(correct_cols),
                   stringr::str_to_lower(colnames(df)))]
  colnames(df) <- correct_cols
  return(df)
}

#' Build data frame with calibration information for each model and set of taxa 
#' (all/obs/novel)
#' If binned = TRUE, will return values for binned intervals 
get_calibration <- function(results, data_true, observed_everywhere, 
                            binned = F, bins = seq(0, 1, 0.05)){
  
  ranks <- colnames(data_true)[-1]
  calibration <- data.frame()
  
  sets <- c("All", "Observed", "Novel")
  
  for(set in sets){
    
    for(rank in ranks){
      
      col <- which(colnames(data_true) == rank)
      col_prob <- which(colnames(results[[1]]) == paste0("Prob_", rank))
      
      for(i in 1:length(results)){
        
        if(set == "All"){
          correct <- results[[i]][,col] == data_true[,col]
          prob <- results[[i]][,col_prob]
        }else if(set == "Observed"){
          correct <- results[[i]][observed_everywhere,col] == data_true[observed_everywhere,col]
          prob <- results[[i]][observed_everywhere,col_prob]
        }else if(set == "Novel"){
          correct <- results[[i]][!observed_everywhere,col] == data_true[!observed_everywhere,col]
          prob <- results[[i]][!observed_everywhere,col_prob]
        }
        
        # Skip over NAs
        na_pos <- which(is.na(correct))
        if(length(na_pos) > 0){
          correct <- correct[-na_pos]
          prob <- prob[-na_pos]
        }
        
        if(binned){
          cal <- calc_calibration_binned(prob, correct, bins)
          
          cal <- cal |>
            dplyr::mutate(model = names(results)[i],
                          rank = rank,
                          set = set)
          
          calibration <- rbind(calibration,
                               cal)
        }else{
          cal <- calc_calibration(prob, correct)
          
          calibration <- rbind(calibration, 
                               data.frame(model = names(results)[i], 
                                          set = set, 
                                          rank = rank, 
                                          cumcorr = cal$cumcorr, 
                                          cumprob = cal$cumprob)) 
        }
      }
    }
  }
  return(calibration)
}

#' Calculate calibration curve for a single model, rank and taxa set 
calc_calibration <- function(prob, correct){
  n <- length(prob)
  s <- sort(prob, dec = F, index.return = T)
  cumprob <- cumsum(prob[s$ix])/n * 100 
  cumcorr <- cumsum(correct[s$ix])/n * 100 
  return(list(cumprob = cumprob, 
              cumcorr = cumcorr))
}

calc_calibration_binned <- function(prob, correct, bins){
  
  n <- length(prob)
  
  df <- 
    data.frame(prob = prob,
               correct = correct) |> 
    dplyr::mutate(bin = cut(prob, bins, include.lowest = T)) |>
    dplyr::group_by(bin) |> 
    dplyr::summarise(correct = sum(correct)/dplyr::n(), 
                     count = dplyr::n()) 
  
  return(df)
}

# My attempt for flexible plotting of calibrations 
# TODO: When plotting a single model + set, include accuracies in legend 
# TODO: When plotting a single rank + set, include accuracies in legend + info about sample size 
# TODO: sample size in All Observed and Novel headings? 
plot_calibration <- function(calibrations){
  
  # Make plot base 
  p <- ggplot2::ggplot() +
    ggplot2::xlim(0, 100) +
    ggplot2::ylim(0, 100) +
    ggplot2::theme_bw() +
    ggplot2::theme(aspect.ratio = 1) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "grey") +
    ggplot2::xlab("Cumulative probability %") +
    ggplot2::ylab("Cumulative correct %") +
    ggpubr::grids(linetype = "dashed")
  
  # Which of model rank and set has more than one unique value? 
  groups <- c()
  
  for(col in c("model", "rank", "set")){
    if(length(unique(calibrations[,col])) > 1) groups <- c(groups, col)
  }
  
  if("rank" %in% groups){
    p <- p + 
      ggplot2::geom_line(data = calibrations, 
                         mapping = ggplot2::aes(x = cumprob, 
                                                y = cumcorr, 
                                                color = rank)) + 
      ggplot2::geom_point(data = calibrations |>
                            dplyr::group_by(across(all_of(groups))) |> 
                            dplyr::summarise(cumprob = max(cumprob),
                                      cumcorr = max(cumcorr)), 
                          mapping = ggplot2::aes(x = cumprob,
                                                 y = cumcorr,
                                                 color = rank)) +
      ggplot2::scale_color_manual(name = "Rank", 
                                  limits = colnames(data_true)[-1], 
                                  values = RColorBrewer::brewer.pal(name = "Blues",
                                                                    n = (length(unique(calibrations$rank)) + 2))[3:(length(unique(calibrations$rank)) + 2)])
  }else{
    p <- p + 
      ggplot2::geom_line(data = calibrations, 
                         mapping = ggplot2::aes(x = cumprob, 
                                                y = cumcorr, 
                                                color = model)) +
      ggplot2::geom_point(data = calibrations |>
                            group_by(across(all_of(groups))) |> 
                                       summarise(cumprob = max(cumprob),
                                                 cumcorr = max(cumcorr)), 
                                     mapping = ggplot2::aes(x = cumprob,
                                                            y = cumcorr,
                                                            color = model))
  }

  if(length(groups) == 3){
    p <- p + 
      ggplot2::facet_grid(set~model)
  }else if(all(groups == c("model", "rank"))){
    p <- p + 
      ggplot2::facet_wrap(~model)
  }else if(all(groups == c("rank", "set"))){
    p <- p + 
      ggplot2::facet_wrap(~set)
  }else if(all(groups == c("model", "set"))){
    p <- p + 
      ggplot2::facet_wrap(~set)
  }

  return(p)
}

#' Calculate accuracy from calibrations data frame 
get_accuracy_from_cal <- function(calibrations, groups = c("model", "rank", "set")){
  
  accuracies <- calibrations |> 
    dplyr::group_by(across(all_of(groups))) |> 
    dplyr::summarise(accuracy = max(cumcorr)) 
  
  return(accuracies)
}

get_partitioned_novelty <- function(data_true, results, observed_everywhere){
  
  ranks <- colnames(data_true)[-1]
  
  true_novel_rank <- 
    sapply(data_true[!observed_everywhere, tail(ranks, 1)],
           function(x) stringr::str_split(x, "_")[[1]] |>
             tail(2) |> head(1),
           USE.NAMES = F)
  
  df <- data.frame() 
  
  for(i in 1:length(results)){
    
    sub_df <- 
      dplyr::tibble(true_rank = true_novel_rank,
                    pred = results[[i]][!observed_everywhere, tail(ranks, 1)]) |> 
      dplyr::mutate(pred_rank = purrr::map_chr(pred, 
                                        ~if(!is.na(.x) & stringr::str_ends(.x, "_new")) 
                                          stringr::str_split(.x, "_")[[1]] |> 
                                          tail(2) |> head(1) 
                                        else "None")) |> 
      dplyr::mutate(true_count = dplyr::n(), 
                    .by = true_rank) |> 
      dplyr::summarise(pred_count = dplyr::n(), 
                       .by = c(true_rank, pred_rank, true_count)) |> 
      dplyr::mutate(model = names(results)[i], 
                    perc = (pred_count/true_count)*100) 

    df <- rbind(df, sub_df)
  }
  return(df)
}

plot_accuracies <- function(accuracies){
  
  groups <- c()
  
  for(col in c("model", "rank", "set")){
    if(!(col %in% colnames(accuracies))) next
    if(nrow(unique(accuracies[,col])) > 1) groups <- c(groups, col)
  }
  
  p <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::theme(aspect.ratio = 1) + 
    ggplot2::ylim(0, 100) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank()) + 
    ggplot2::ylab("Accuracy %")
  
  p <- p + 
    ggplot2::geom_line(data = accuracies, 
                       mapping = ggplot2::aes(x = rank, 
                                              y = accuracy, 
                                              color = model, 
                                              group = model)) + 
    ggplot2::geom_point(data = accuracies, 
                        mapping = ggplot2::aes(x = rank,
                                               y = accuracy,
                                               color = model, 
                                               shape = model)) + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) + 
    ggplot2::labs(color = "Model")
  
  if(length(groups) == 3){
    p <- p +
      ggplot2::facet_wrap(~set)
  }
  
  return(p)
}

#' What proportion of sequences was correctly predicted as novel *within the correct higher rank*? 
marginal_accuracy <- function(results, data_true, id_list){
  
  ranks <- colnames(data_true)[-1]
  marg_accuracies <- data.frame()
  
  for(set in names(id_list)){
    
    for(r in ranks){
      
      for(i in 1:length(results)){
        
        positions <- id_list[[set]][[r]]
        na_pos <- which(is.na(results[[i]][positions,r]))
        
        # Skip NAs in the calculation
        if(length(na_pos) > 0) positions <- positions[-na_pos]
        
        correct <- sum(results[[i]][positions,r] == data_true[positions,r])
        
        marg_accuracies <- rbind(marg_accuracies, 
                                 data.frame(model = names(results[i]),
                                            set = set,
                                            rank = r,
                                            marg_accuracy = correct/length(positions)))
      }
    }
  }
  return(marg_accuracies)
}

conditional_accuracy <- function(results, data_true, id_list){
  
  ranks <- colnames(data_true)[-1]
  cond_accuracies <- data.frame()
  
  for(set in names(id_list)){
    
    ids <- id_list[[set]]
    
    for(r in 1:length(ranks)){
      
      for(i in 1:length(results)){
        
        positions <- ids[[ranks[r]]]
        
        # Both the target rank and the rank above 
        if(r == 1){
          na_pos <- which(is.na(results[[i]][positions, ranks[r]]))
        }else{
          na_pos <- union(which(is.na(results[[i]][positions, ranks[r]])),
                          which(is.na(results[[i]][positions, ranks[r-1]])))
        }
        
        if(length(na_pos) > 0) positions <- positions[-na_pos]
        
        # Denominator: novel taxa correctly predicted on the rank above
        if(r == 1){
          denom <- length(positions)
        }else{
          denom <-
            sum(
              results[[i]][positions, ranks[r-1]] == data_true[positions, ranks[r-1]]
            )
        }
        
        correct <- 
          sum(
            results[[i]][positions, ranks[r]] == data_true[positions, ranks[r]]
          )
        
        if(correct > denom) print("Error: number of correct predictions larger than the denominator.")
        
        cond_accuracies <- rbind(cond_accuracies, 
                                 data.frame(model = names(results[i]),
                                            set = set, 
                                            rank = ranks[r],
                                            cond_accuracy = correct/denom, 
                                            denom_cond = denom))
      }
    }
  }
  return(cond_accuracies)
}

count_novel <- function(results, data_true, id_novel){
  
  ranks <- colnames(data_true)[-1]
  pred_df <- data.frame()
  
  for(i in 1:length(results)){
    
    sum_novel <- 0 
    pred_novel <- 0 
    
    for(r in ranks){
      
      sum_novel <- sum_novel + length(id_novel[[r]])
      pred_novel <- pred_novel + 
        sum(stringr::str_ends(results[[i]][,r], "_new"))
      
      df <- data.frame(model = names(results)[i], 
                       rank = r,
                       sum_novel = sum_novel, 
                       pred_novel = pred_novel)
      
      pred_df <- rbind(pred_df, df)
    }
  }
  return(pred_df)
}

#' @param taxdf data frame with allowed taxonomy. Must only contain rank columns (ordered). 
check_taxonomy <- function(taxdf, results, ranks){
  
  # Build a list of all allowed taxonomic combinations 
  taxlist <- unique(taxdf[,1])
  for(r in 2:length(ranks)){
    taxlist <- c(taxlist, 
                 apply(taxdf[1:r] |> dplyr::distinct(), 1, function(row) paste(row, collapse = "|")))
  }
  
  mismatch <- data.frame() 
  
  for(i in 1:length(results)){
    
    print(names(results)[i])
    
    df <- results[[i]][,ranks] |> dplyr::distinct()
    
    for(j in 1:nrow(df)){
      
      if(is.na(df[j,1])) next
      
      df_row <- df[j,which(!is.na(df[j,]))]
      
      nonew_ranks <- which(!stringr::str_ends(df_row, "_new"))
      
      # Skip if Class_new
      if(length(nonew_ranks) == 0) next 
      
      taxstr <- paste(df[j, 1:max(nonew_ranks)], collapse = "|")
      match <- taxstr %in% taxlist
      
      if(!match) mismatch <- rbind(mismatch, 
                                   df[j,] |> dplyr::mutate(model = names(results)[i]))
    }
  }
  return(mismatch)
}

threshold_curve <- function(results, data_true, thresholds = seq(0, 1, 0.01)){
  
  out <- data.frame()
  ranks <- colnames(data_true)[-1]
  n <- nrow(data_true)
  
  # Find which have single value 1 for all classification probabilities 
  point_models <- 
    names(which(sapply(results, 
                       function(x) all(x$Prob_Species == 1 | is.na(x$Prob_Species)))))
  
  for(i in 1:length(results)){
    
    for(r in ranks){
      
      probcol <- paste0("Prob_", r)
      
      not_na_pos <- which(!is.na(results[[i]][,r]))
      
      if(names(results)[i] %in% point_models){
        
        correct <- sum(results[[i]][not_na_pos, r] == data_true[not_na_pos, r])/nrow(results[[i]][not_na_pos,]) * 100
        classified <- nrow(results[[i]][not_na_pos,])/n * 100
        
        out <- rbind(out,
                     data.frame(model = names(results)[i],
                                rank = r,
                                threshold = unique(results[[i]][not_na_pos, probcol]),
                                correct = correct,
                                classified = classified))
        
      }else{

        results_without_na <- results[[i]][not_na_pos,]
        data_true_without_na <- data_true[not_na_pos,]

        correct <- sapply(thresholds, function(x)
          sum(results_without_na[which(results_without_na[, probcol] >= x), r] == data_true_without_na[which(results_without_na[, probcol] >= x), r]) /
            length(which(results_without_na[, probcol] >= x)) * 100)
        classified <- sapply(thresholds,
                             function(x) length(which(results_without_na[,probcol] >= x))/n * 100)
        
        out <- rbind(out,
                     data.frame(model = names(results)[i],
                                rank = r,
                                threshold = thresholds,
                                correct = correct,
                                classified = classified))
      }
    }
  }
  return(out)
}

underclass_rate <- function(results, id_observed){
  uclass_df <- data.frame() 
  
  for(i in 1:length(results)){
    
    uclass_df <- 
      rbind(uclass_df, 
            data.frame(underclassification = sapply(ranks, function(x)
              sum(
                stringr::str_ends(results[[i]][id_observed[[x]], x], "_new") |
                  is.na(results[[i]][id_observed[[x]], x])
              )/length(id_observed[[x]]) * 100)) |> 
              tibble::rownames_to_column("rank") |> 
              dplyr::mutate(model = names(results)[i]))
  }
  return(uclass_df)
}

overclass_rate <- function(results, id_novel){
  oclass_df <- data.frame() 
  
  for(i in 1:length(results)){
    
    oclass_df <- 
      rbind(oclass_df, 
            data.frame(overclassification = sapply(ranks, function(x)
              sum(
                !stringr::str_ends(results[[i]][id_novel[[x]], x], "_new") &
                  !is.na(results[[i]][id_novel[[x]], x])
              )/length(id_novel[[x]]) * 100)) |> 
              tibble::rownames_to_column("rank") |> 
              dplyr::mutate(model = names(results)[i]))
  }
  return(oclass_df)
}

misclass_rate <- function(results, data_true, id_observed){
  misclass_df <- data.frame() 
  
  for(i in 1:length(results)){
    
    misclass_df <- 
      rbind(misclass_df, 
            data.frame(misclassification = sapply(ranks, function(x)
              sum(
                results[[i]][id_observed[[x]], x] != data_true[id_observed[[x]], x] &
                  !is.na(results[[i]][id_observed[[x]], x]) & 
                  !stringr::str_ends(results[[i]][id_observed[[x]], x], "_new")
              )/length(id_observed[[x]]) * 100)) |> 
              tibble::rownames_to_column("rank") |> 
              dplyr::mutate(model = names(results)[i]))
  }
  return(misclass_df)
}
