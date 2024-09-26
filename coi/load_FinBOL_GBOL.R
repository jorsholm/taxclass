
#' Loads FinBOL data (train) and GBOL data (test) as listed dataframes. 
#' @param level lowest taxonomic level (genus/species)
load_FinBOL_GBOL <- function(train_file = "data/train_tax.tsv",
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
