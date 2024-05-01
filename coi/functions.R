################################################################################
# Functions for analysing performance of taxonomic classification algorithms.
################################################################################
library(tidyverse)

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

make_plot_base <- function(){
  # Base for the plot 
  p <- ggplot2::ggplot() +
    ggplot2::xlim(0, 100) +
    ggplot2::ylim(0, 100) +
    ggplot2::theme_bw() +
    ggplot2::theme(aspect.ratio = 1) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "grey") +
    ggplot2::xlab("Cumulative probability %") +
    ggplot2::ylab("Cumulative correct %") +
    ggpubr::grids(linetype = "dashed")
  
  return(p)
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

# My attempt for flexible plotting of calibrations 
# TODO: When plotting a single model + set, include accuracies in legend 
# TODO: When plotting a single rank + set, include accuracies in legend + info about sample size 
# TODO: sample size in All Observed and Novel headings? 
plot_calibration <- function(calibrations){
  p <- make_plot_base()
  
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
                            group_by(across(all_of(groups))) |> 
                            summarise(cumprob = max(cumprob),
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

#' Calculate accuracy on a single result data frame 
#' @param x Result data frame 
calc_accuracy <- function(x, data_true){
  return(colMeans(x[,2:ncol(data_true)] == data_true[,-1]))
}

#' Create a single accuracy facet 
plot_acc_facet <- function(df, data_true){
  accuracies <- t(sapply(df, function(x) calc_accuracy(x, data_true))) 
  
  plot_acc <- cbind("Model" = rownames(accuracies),
                    data.frame(accuracies, row.names = NULL)) |>
    pivot_longer(2:ncol(data_true), names_to = "Rank", values_to = "Accuracy")
  
  plot_acc$Rank <- factor(plot_acc$Rank, levels = colnames(data_true)[-1])
  
  p <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::theme(aspect.ratio = 1) + 
    ggplot2::ylim(0, 1)
  
  p <- p + 
    ggplot2::geom_line(data = plot_acc, 
                       aes(x = Rank, y = Accuracy, 
                           color = Model, group = Model)) + 
    ggplot2::geom_point(data = plot_acc, 
                        aes(x = Rank, y = Accuracy, color = Model)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

plot_accuracies <- function(results, data_true, id_observed = NULL){
  
  # bind column with "all taxa", "known", "unknown"
  # Combine into the same data frame 
  t(sapply(results, function(x) calc_accuracy(x, data_true)))
  
  if(is.null(id_observed)){
    return(plot_acc_facet(df = results, data_true))
  }else{
    plotlist <- list()
    
    plotlist[[1]] <- plot_acc_facet(df = results, data_true)
    plotlist[[2]] <- plot_acc_facet(df = lapply(results, function(x) x[id_observed, ]), 
                                    data_true[id_observed, ])
    plotlist[[3]] <- plot_acc_facet(df = lapply(results, function(x) x[!id_observed, ]), 
                                    data_true[!id_observed, ])
    
    ggpubr::ggarrange(plotlist = plotlist, 
                      ncol = 3, nrow = 1, 
                      labels = c("All taxa", "Observed taxa", "New taxa"), 
                      common.legend = T, legend = "right")
    
  }
}





