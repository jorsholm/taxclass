################################################################################
# Functions for analysing performance of taxonomic classification algorithms.
################################################################################
#library(tidyverse)

#' Replaces labels of classes that are unique to test dataset. 
#' Example of new labels: Araneae_Family_new
#' Or, if first level is new: Class_new
#' @param test 
#' @param train 
get_data_true <- function(test, train){

  # Remove ID and DNA columns
  train <- train[, -c(1, ncol(train))]
  data_true <- test[, -c(1, ncol(test))]

  # Get level names
  level_names <- names(train)
  n_levels <- length(level_names)

  for(l in 1:n_levels){
    # Find taxa that are unique to test
    id <- !(data_true[,l] %in% unique(train[,l]))
    n_new <- sum(id)

    # Replace new classes with "_new"-name
    if(n_new > 0){
      if(l == 1){
        name_clust <- paste0(level_names[1], "_new")
        data_true[id, l:n_levels] <- matrix(rep(name_clust, n_new*n_levels),
                                            nrow = n_new, ncol = n_levels)
      }else{

        for(i in which(id == TRUE)){
          if(!grepl("_new$", data_true[i, l-1])){
            name_clust <- paste0(data_true[i, l-1], "_", level_names[l], "_new")
            data_true[i, l:n_levels] <- rep(name_clust, n_levels-l + 1)
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
rename_unk <- function(x){
  is_unk <- x == "unk"
  
  if(sum(is_unk) > 0){
    if(all(x == "unk")){
      x[1:length(x)] <- "Class_new"
    }else{
      d <- min(which(is_unk))
      x[d:length(x)] <- paste0(x[d-1], "_", stringr::str_to_title(names(x)[d]), "_new")
    }
  }
  return(x)
}

#' Renames PROTAX unk to BayesANT standard + sort columns 
#' @param out data frame with PROTAX results. Cannot contain NA values  
rename_PROTAX_output <- function(out, level = "species"){
  # Rename unk 
  renamed <- cbind("ID" = out[,which(colnames(out) == "ID")],
                   t(apply(out[, which(!stringr::str_detect(colnames(out), "Prob|ID"))],
                           1, function(x) rename_unk(x))), 
                   out[,which(stringr::str_detect(colnames(out), "Prob_"))])
  
  return(renamed)
}

#' Build data frame with calibration information for each model and set of taxa 
#' (all/obs/novel)
get_calibration <- function(results, data_true, id_observed){
  
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
          correct <- results[[i]][id_observed,col] == data_true[id_observed,col]
          prob <- results[[i]][id_observed,col_prob]
        }else if(set == "Novel"){
          correct <- results[[i]][!id_observed,col] == data_true[!id_observed,col]
          prob <- results[[i]][!id_observed,col_prob]
        }
        
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

get_binned_calibration <- function(results, data_true, observed_everywhere, 
                                   bins = seq(0, 1, 0.05)){
  
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
        
        cal <- calc_calibration_binned(prob, correct, bins)
        
        cal <- cal |> 
          dplyr::mutate(model = names(results)[i], 
                        rank = rank, 
                        set = set)
        
        calibration <- rbind(calibration, 
                             cal)
      }
    }
  }
  return(calibration)
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
      ggplot2::facet_grid(model~set)
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

get_novelty_accuracy <- function(results, id_novel){
  ranks <- colnames(data_true)[-1]
  n <- list() 
  
  acc_df <- data.frame()
  part_novelty <- data.frame()
  
  for(r in ranks){
    n[[r]] <- length(id_novel[[r]])
    
    for(i in 1:length(results)){
      novel_acc <- sum(stringr::str_ends(results[[i]][id_novel[[r]], r], "_new"))/n[[r]] * 100
      
      # False novelty is calc as proportion of falsely predicted novel species on that rank / true novel on that rank  
      false_novelty <- sum(stringr::str_ends(results[[i]][-id_novel[[r]], r],  paste0(r, "_new")))/n[[r]] * 100
      
      acc_df <- rbind(acc_df, 
                data.frame(model = names(results)[i], 
                           rank = r, 
                           novel_accuracy = novel_acc, 
                           false_novelty = false_novelty)) 

      suppressMessages(part <- 
        data.frame(rank = r, 
                 model = names(results)[i], 
                 novelty_rank = 
                   sapply(results[[i]][id_novel[[r]],r],
                          function(x) dplyr::if_else(stringr::str_detect(x, "_new"),
                                              stringr::str_split(x, "_")[[1]] |> tail(2) |> head(1),
                                              "None"),
                          USE.NAMES = F)) |> 
        dplyr::group_by(rank, model, novelty_rank) |> 
        dplyr::summarise(count = dplyr::n()) |> 
        dplyr::mutate(perc = count / n[[r]] * 100))
      
      part_novelty <- rbind(part_novelty, part)
      
    }
  }
  return(list(count = n, 
              acc = acc_df, 
              part_novelty = part_novelty))
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
                                               color = model)) + 
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
  
  for(r in ranks){
    taxgroups <- unique(data_true[id_list[[r]],r])
    
    denom <- length(data_true[id_list[[r]], r])
    
    for(i in 1:length(results)){
      
      correct <- 0
      
      for(tax in taxgroups){
        ids <- which(data_true[,r] == tax)
        
        correct <- correct + sum(results[[i]][ids,r] == tax)
      }
      
      marg_accuracies <- rbind(marg_accuracies, 
                               data.frame(model = names(results[i]),
                                          rank = r,
                                          marg_accuracy = correct/denom))
      
    }
  }
  return(marg_accuracies)
}

conditional_accuracy <- function(results, data_true, id_list){
  
  ranks <- colnames(data_true)[-1]
  cond_accuracies <- data.frame()
  
  for(r in ranks){
    taxgroups <- unique(data_true[id_list[[r]],r])
    
    for(i in 1:length(results)){
      
      correct <- 0
      denom <- 0 
      
      for(tax in taxgroups){
        # Which are actually this taxon? 
        ids_sugg <- which(data_true[,r] == tax)
        
        if(r == "Class"){
          correct <- correct + sum(results[[i]][ids_sugg,r] == tax)
          denom <- length(id_list[[r]])
        }else{
          
          higher_rank <- ranks[which(ranks==r)-1]
          
          higher_tax <- data_true[ids_sugg, higher_rank][1]
          
          # Which are this taxon + have been correctly predicted to the higher rank 
          correct_higher <- ids_sugg[which(results[[i]][ids_sugg, higher_rank] == higher_tax)]
          
          denom <- denom + length(correct_higher)
          correct <- correct + sum(results[[i]][correct_higher,r] == tax)
        }
      }
      cond_accuracies <- rbind(cond_accuracies, 
                               data.frame(model = names(results[i]),
                                          rank = r,
                                          cond_accuracy = correct/denom, 
                                          denom_cond = denom))
    }
  }
  return(cond_accuracies)
}
# 