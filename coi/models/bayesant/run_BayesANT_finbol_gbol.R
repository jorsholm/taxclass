# This file runs BayesANT for the finbol-gbol split (new version)
library(tidyverse)
library(BayesANT)
library(seqinr)
library(ggpubr)
library(foreach)
library(doParallel)

readBayesANT <- function(fasta.file, rank = NULL, rank_names = NULL, sep = ";") {
  data_fasta <- seqinr::read.fasta(file = fasta.file, forceDNAtolower = F, as.string = T)
  annot <- lapply(data_fasta, function(x) attr(x, "Annot"))
  annot <- gsub(">", "", annot)
  IDs <- names(data_fasta)
  taxa <- stringr::str_split(annot, pattern = sep, simplify = T)
  taxa <- taxa[, -1]
  if (is.null(rank)) {
    lv <- ncol(taxa)
  } else {
    lv <- rank
    if (rank > ncol(taxa)) {
      stop(cat(
        "Value for rank = ", rank, "exceeds the length for the taxonomic annotations in the library. Specify a value for rank lower or equal to ",
        ncol(taxa)
      ))
    }
  }
  taxa <- taxa[, 1:lv]
  if (is.null(rank_names)) {
    rank_names <- paste0("Level", c(1:lv))
  }
  colnames(taxa) <- rank_names[1:lv]
  data <- data.frame(cbind(taxa, DNA = unlist(lapply(
    data_fasta,
    function(x) {
      stringr::str_to_upper(x)
    }
  ))))
  class(data) <- c("data.frame", "BayesANT.data")
  return(data)
}

# Replace with directory where the files are
main_dir <- "~/TaxonomicClassifier/FinBOL/data/finbol-gbol/"
out_dir <- "~/TaxonomicClassifier/FinBOL/models/results_bayesant_finbol_gbol/"
registerDoParallel(22)
################################################################################
# Part 1 - train and test using the .nt file
################################################################################

#----- Training

# Load data
file.train <- paste0(main_dir, "train_finbol-gbol_nt.fasta")
rank_names <- c("Class", "Order", "Family", "Subfamily", "Genus", "Tribe", "Species")
data_train <- readBayesANT(file.train, rank_names = rank_names, sep = "\\s+|\\|")

# Train BayesANT
time_train <- Sys.time()
classifier <- BayesANT(
  data = data_train,
  typeseq = "aligned",
  type_location = "single",
  newtaxa = TRUE,
  verbose = TRUE
)
time_train <- as.numeric(Sys.time() - time_train)
size_classifier <- as.numeric(object.size(classifier)) / 1e6

#-----  Predict out of sample
# Full sequences
file.test <- paste0(main_dir, "test_finbol-gbol_nt.fasta")
data_test <- read.BayesANT.testDNA(file.test, rank_names = rank_names, sep = "\\s+|\\|")
# Predict nt
time_pred <- Sys.time()
output <- predict(classifier, DNA = data_test, rho = 0.1, cores = 22)
time_pred <- Sys.time() - time_pred

# Testshort sequences
file.testshort <- paste0(main_dir, "testshort_finbol-gbol_nt.fasta")
data_testshort <- read.BayesANT.testDNA(file.testshort, rank_names = rank_names, sep = "\\s+|\\|")
# Predict testshort nt
time_pred_short <- Sys.time()
output_short <- predict(classifier, DNA = data_testshort, rho = 0.1, cores = 22)
time_pred_short <- Sys.time() - time_pred_short

# Summary of prediction
summary_prediction <- as.matrix(c(
  "size_classifier_Mb" = size_classifier,
  "time_training" = time_train,
  "time_prediction" = time_pred,
  "time_prediction_short" = time_pred_short
))

# Save the outputs
write.table(output, file = paste0(out_dir, "bayesant_", "test_finbol-gbol_nt.txt"))
write.table(output_short, file = paste0(out_dir, "bayesant_", "testshort_finbol-gbol_nt.txt"))
write.table(summary_prediction, file = paste0(out_dir, "bayesant_", "summary_finbol-gbol_nt.txt"))

################################################################################
# Part 2 - train and test using the .translated files
################################################################################

# Note: there is no Tribe here

#----- Training
# "FinBOL/data/finbol-gbol/train_finbol-gbol_species.translated.fasta"
# "FinBOL/data/finbol-gbol/testshort_finbol-gbol_species.translated.fasta"
# "FinBOL/data/finbol-gbol/test_finbol-gbol_species.translated.fasta"

# Load data
file.train <- paste0(main_dir, "train_finbol-gbol_species.translated.fasta")
rank_names <- c("Class", "Order", "Family", "Subfamily", "Genus", "Species")
data_train <- readBayesANT(file.train, rank_names = rank_names, sep = "\\s+|\\|")

data_fasta <- seqinr::read.fasta(file = file.train, forceDNAtolower = F, as.string = T)
head(data_fasta)
data_train$Tribe
# Train BayesANT
time_train <- Sys.time()
classifier <- BayesANT(
  data = data_train,
  typeseq = "aligned",
  type_location = "single",
  newtaxa = TRUE,
  verbose = TRUE
)
time_train <- as.numeric(Sys.time() - time_train)
size_classifier <- as.numeric(object.size(classifier)) / 1e6

#-----  Predict out of sample
# Full sequences
file.test <- paste0(main_dir, "test_finbol-gbol_species.translated.fasta")
data_test <- read.BayesANT.testDNA(file.test, rank_names = rank_names, sep = "\\s+|\\|")
# Predict nt
time_pred <- Sys.time()
output <- predict(classifier, DNA = data_test, rho = 0.1, cores = 22)
time_pred <- Sys.time() - time_pred

# Testshort sequences
file.testshort <- paste0(main_dir, "testshort_finbol-gbol_species.translated.fasta")
data_testshort <- read.BayesANT.testDNA(file.testshort, rank_names = rank_names, sep = "\\s+|\\|")
# Predict testshort nt
time_pred_short <- Sys.time()
output_short <- predict(classifier, DNA = data_testshort, rho = 0.1, cores = 22)
time_pred_short <- Sys.time() - time_pred_short

# Summary of prediction
summary_prediction <- as.matrix(c(
  "size_classifier_Mb" = size_classifier,
  "time_training" = time_train,
  "time_prediction" = time_pred,
  "time_prediction_short" = time_pred_short
))

# Save the outputs
write.table(output, file = paste0(out_dir, "bayesant_", "test_finbol-gbol_species.translated.txt"))
write.table(output_short, file = paste0(out_dir, "bayesant_", "testshort_finbol-gbol_species.translated.txt"))
write.table(summary_prediction, file = paste0(out_dir, "bayesant_", "summary_finbol-gbol_species.translated.txt"))


################################################################################
# Part 3 - train and test using the genus translated files
################################################################################

# Note: there is no species here, but I detected the Tribe

#----- Training
"FinBOL/data/finbol-gbol/train_finbol-gbol_genus.translated.fasta"
"FinBOL/data/finbol-gbol/test_finbol-gbol_genus.translated.fasta"
"FinBOL/data/finbol-gbol/testshort_finbol-gbol_genus.translated.fasta"

# Load data
file.train <- paste0(main_dir, "train_finbol-gbol_genus.translated.fasta")
rank_names <- c("Class", "Order", "Family", "Subfamily", "Genus", "Tribe")
data_train <- readBayesANT(file.train, rank_names = rank_names, sep = "\\s+|\\|")

data_fasta <- seqinr::read.fasta(file = file.train, forceDNAtolower = F, as.string = T)
head(data_fasta)
data_train$Tribe
# Train BayesANT
time_train <- Sys.time()
classifier <- BayesANT(
  data = data_train,
  typeseq = "aligned",
  type_location = "single",
  newtaxa = TRUE,
  verbose = TRUE
)
time_train <- as.numeric(Sys.time() - time_train)
size_classifier <- as.numeric(object.size(classifier)) / 1e6

#-----  Predict out of sample
# Full sequences
file.test <- paste0(main_dir, "test_finbol-gbol_genus.translated.fasta")
data_test <- read.BayesANT.testDNA(file.test, rank_names = rank_names, sep = "\\s+|\\|")
# Predict nt
time_pred <- Sys.time()
output <- predict(classifier, DNA = data_test, rho = 0.1, cores = 22)
time_pred <- Sys.time() - time_pred

# Testshort sequences
file.testshort <- paste0(main_dir, "testshort_finbol-gbol_genus.translated.fasta")
data_testshort <- read.BayesANT.testDNA(file.testshort, rank_names = rank_names, sep = "\\s+|\\|")
# Predict testshort nt
time_pred_short <- Sys.time()
output_short <- predict(classifier, DNA = data_testshort, rho = 0.1, cores = 22)
time_pred_short <- Sys.time() - time_pred_short

# Summary of prediction
summary_prediction <- as.matrix(c(
  "size_classifier_Mb" = size_classifier,
  "time_training" = time_train,
  "time_prediction" = time_pred,
  "time_prediction_short" = time_pred_short
))

# Save the outputs
write.table(output, file = paste0(out_dir, "bayesant_", "test_finbol-gbol_genus.translated.txt"))
write.table(output_short, file = paste0(out_dir, "bayesant_", "testshort_finbol-gbol_genus.translated.txt"))
write.table(summary_prediction, file = paste0(out_dir, "bayesant_", "summary_finbol-gbol_genus.txt"))

################################################################################
# Part 4 - aminoacids
################################################################################

# Note: there is no species here, but I detected the Tribe

#----- Training
"FinBOL/data/finbol-gbol/train_finbol-gbol_aa.fasta"
"FinBOL/data/finbol-gbol/test_finbol-gbol_aa.fasta"
"FinBOL/data/finbol-gbol/testshort_finbol-gbol_aa.fasta"

# Load data
file.train <- paste0(main_dir, "train_finbol-gbol_aa.fasta")
rank_names <- c("Class", "Order", "Family", "Subfamily", "Genus", "Tribe", "Species")
data_train <- readBayesANT(file.train, rank_names = rank_names, sep = "\\s+|\\|")
data_train$DNA

data_fasta <- seqinr::read.fasta(file = file.train, forceDNAtolower = F, as.string = T)

# Train BayesANT
time_train <- Sys.time()
classifier <- BayesANT_amino(
  data = data_train,
  typeseq = "aligned",
  type_location = "single",
  newtaxa = TRUE,
  verbose = TRUE, amino = TRUE,
)
time_train <- as.numeric(Sys.time() - time_train)
size_classifier <- as.numeric(object.size(classifier)) / 1e6

#-----  Predict out of sample
# Full sequences
file.test <- paste0(main_dir, "test_finbol-gbol_aa.fasta")
data_test <- read.BayesANT.testDNA(file.test, rank_names = rank_names, sep = "\\s+|\\|")
# Predict nt
time_pred <- Sys.time()
output <- predict(classifier, DNA = data_test, rho = 0.3, cores = 22)
time_pred <- Sys.time() - time_pred

# Testshort sequences
file.testshort <- paste0(main_dir, "testshort_finbol-gbol_aa.fasta")
data_testshort <- read.BayesANT.testDNA(file.testshort, rank_names = rank_names, sep = "\\s+|\\|")
# Predict testshort nt
time_pred_short <- Sys.time()
output_short <- predict(classifier, DNA = data_testshort, rho = 0.3, cores = 22)
time_pred_short <- Sys.time() - time_pred_short

# Summary of prediction
summary_prediction <- as.matrix(c(
  "size_classifier_Mb" = size_classifier,
  "time_training" = time_train,
  "time_prediction" = time_pred,
  "time_prediction_short" = time_pred_short
))

# Save the outputs
write.table(output, file = paste0(out_dir, "bayesant_", "test_finbol-gbol_aa.txt"))
write.table(output_short, file = paste0(out_dir, "bayesant_", "testshort_finbol-gbol_aa.txt"))
write.table(summary_prediction, file = paste0(out_dir, "bayesant_", "summary_finbol-gbol_aa.txt"))

