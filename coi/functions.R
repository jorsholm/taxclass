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

#' Data frames in result must have identical colnames with each other, 
#' and colnames of taxonomic ranks (Species, Genus, etc.) must be identical to 
#' data_true (including order of columns). 
#' @param results Named list of data frames with tax class results for each model.
#' @param data_true Data frame with the true taxonomic classifications 
#' @param level Taxonomic rank to plot 
#' @return Calibration plot 
plot_comparison <- function(results, data_true, level = "Species"){
  
  # Find columns that contain taxonomic classifications + probabilities 
  col <- which(colnames(data_true) == level)
  col_prob <- which(colnames(results[[1]]) == paste0("Prob_", level))
  
  # Base for the plot 
  p <- make_plot_base()
  
  # Add calibration line + point for each model 
  for(i in 1:length(results)){
    algorithm <- names(results)[i]
    prob <- results[[i]][,col_prob]
    correct <- results[[i]][,col] == data_true[,col]
    plotelements <- plot_calibration(prob, correct, algorithm)
    p <- p + plotelements[[1]] + plotelements[[2]]
  }
  
  # Add header with sample size 
  p <- p + 
    ggplot2::facet_wrap(~paste0(level, " - sample size: ", nrow(data_true)))
  
  return(p)
}

#' Plot calibration line + point for a single model 
#' @param prob vector of classification probabilities 
#' @param correct vector of T/F taxonomic assignment results 
#' @param algorithm name of model 
#' @return list of plot elements; one line and one point 
plot_calibration <- function(prob, correct, algorithm){
  s <- sort(prob, dec = F, index.return = T)
  n <- length(prob)
  x <- cumsum(prob[s$ix])/n * 100 
  y <- cumsum(correct[s$ix])/n * 100 
  df <- data.frame(x, y)
  
  Model <- paste0(algorithm, " - ", 
                     as.character(round(mean(correct) * 100, 1)), "%")
  l <- ggplot2::geom_line(data = df, 
                         aes(x = x, y = y, color = Model))
  m <- ggplot2::geom_point(aes(x = x[length(x)], 
                               y = y[length(y)], 
                               color = Model))
  
  return(list(l, m))
}

#' Plots calibration of a single model, across all ranks 
plot_ranks <- function(model_result, data_true, model_name){
  ranks <- colnames(data_true)[-1]
  
  p <- make_plot_base()
  limits <- c() 
  
  for(r in ranks){
    
    col <- which(colnames(model_result) == r)
    col_prob <- which(colnames(model_result) == paste0("Prob_", r))
    
    prob <- model_result[,col_prob]
    correct <- model_result[,col] == data_true[,col]
    
    limits <- c(limits,
                paste0(r, " - ", 
                       as.character(round(mean(correct) * 100, 1)), "%"))
    
    p <- p + 
      plot_calibration(prob, correct, r)

  }
  
  blue_palette <- RColorBrewer::brewer.pal(name = "Blues",
                        n = (length(ranks) + 3))[3:(length(ranks) + 3)]
  
  p <- p + 
    ggplot2::scale_color_manual(name = "Rank", 
                                limits = limits, 
                                values = blue_palette) + 
    ggplot2::ggtitle(model_name) + 
    ggplot2::theme(plot.title = element_text(hjust = 0.5))
    
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

plot_cal_unknown <- function(model_result, data_true, model_name, id_observed, legend = "bottom"){
  
  # All species 
  all <- plot_ranks(model_result, data_true, model_name = paste(model_name, "- All taxa"))
  obs <- plot_ranks(model_result[id_observed,], data_true[id_observed,], model_name = paste(model_name, "- Observed taxa"))
  unobs <- plot_ranks(model_result[!id_observed,], data_true[!id_observed,], model_name = paste(model_name, "- New taxa"))
  
  plotlist <- list(all, obs, unobs)
  
  p <- ggpubr::ggarrange(plotlist = plotlist, ncol = 3, legend = legend)
  
  return(p)
  
}



