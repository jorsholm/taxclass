# This file runs BayesANT for the finbol-gbol split (new version)
library(tidyverse)
library(BayesANT)
library(seqinr)
library(ggpubr)
library(foreach)
library(doParallel)

# Replace with directory where the files are
data_dir <- "~/taxclass/coi/data/"
out_dir <- "~/taxclass/coi/results/bayesant/"
registerDoParallel(22)

################################################################################
# Part 1 - train and test using the nt file with labels
################################################################################

#----- Training
"~/taxclass/coi/data/train_nt_aln_label.fasta"
"~/taxclass/coi/data/test_nt_aln_label.fasta"
"~/taxclass/coi/data/testshort_nt_aln_label.fasta"

# Load data
file.train <- paste0(data_dir, "train_nt_aln_label.fasta")
rank_names <- c("Class", "Order", "Family", "Subfamily", "Tribe", "Genus", "Species")
data_train <- read.BayesANT.data(file.train, rank_names = rank_names, sep = "\\s+|\\|")

# Train BayesANT
classifier <- BayesANT(
  data = data_train,
  typeseq = "aligned",
  type_location = "single",
  newtaxa = TRUE,
  verbose = TRUE
)

#-----  Predict out of sample
# Full sequences
file.test <- paste0(data_dir, "test_nt_aln_label.fasta")
data_test <- read.BayesANT.testDNA(file.test, rank_names = rank_names, sep = "\\s+|\\|")
# Predict nt
output <- predict(classifier, DNA = data_test, rho = 0.1, cores = 22)

# Testshort sequences
file.testshort <- paste0(data_dir, "testshort_nt_aln_label.fasta")
data_testshort <- read.BayesANT.testDNA(file.testshort, rank_names = rank_names, sep = "\\s+|\\|")
# Predict testshort nt
output_short <- predict(classifier, DNA = data_testshort, rho = 0.1, cores = 22)

# Save the outputs
write.table(output, file = paste0(out_dir, "bayesant_", "test_nt_aln_label.txt"))
write.table(output_short, file = paste0(out_dir, "bayesant_", "testshort_nt_aln_label.txt"))

#--------------------------------- Show accuracies 
output <- read.table(paste0(out_dir, "bayesant_", "test_nt_aln_label.txt"))
data_test_lab <- read.BayesANT.data(file.test, rank_names = rank_names, sep = "\\s+|\\|")
data_test_lab <- add_novelty_to_test_data(data_train, data_test_lab)
plot_accuracies(output, data_test_lab)

# Show accuracies in full case
output_short <- read.table(paste0(out_dir, "bayesant_", "testshort_nt_aln_label.txt"))
plot_accuracies(output_short, data_test_lab)
#---------------------------------




################################################################################
# Part 2 - aminoacids
################################################################################

# Note: there is no species here, but I detected the Tribe

#----- Training
"~/taxclass/coi/data/train_aa_aln_label.fasta"
"~/taxclass/coi/data/test_aa_aln_label.fasta"
"~/taxclass/coi/data/testshort_aa_aln_label.fasta"

# Load data
file.train <- paste0(data_dir, "train_aa_aln_label.fasta")
rank_names <- c("Class", "Order", "Family", "Subfamily", "Tribe", "Genus", "Species")
data_train <- read.BayesANT.data(file.train, rank_names = rank_names, sep = "\\s+|\\|")

# Train BayesANT
classifier <- BayesANT_amino(
  data = data_train,
  typeseq = "aligned",
  type_location = "single",
  newtaxa = TRUE,
  verbose = TRUE, amino = TRUE,
)

#-----  Predict out of sample
# Full sequences
file.test <- paste0(data_dir, "test_aa_aln_label.fasta")
data_test <- read.BayesANT.testDNA(file.test, rank_names = rank_names, sep = "\\s+|\\|")
# Predict nt
output <- predict(classifier, DNA = data_test, rho = 0.3, cores = 22)

# Testshort sequences
file.testshort <- paste0(data_dir, "testshort_aa_aln_label.fasta")
data_testshort <- read.BayesANT.testDNA(file.testshort, rank_names = rank_names, sep = "\\s+|\\|")
# Predict testshort nt
output_short <- predict(classifier, DNA = data_testshort, rho = 0.3, cores = 22)

# Save the outputs
write.table(output, file = paste0(out_dir, "bayesant_", "test_aa_aln_label.txt"))
write.table(output_short, file = paste0(out_dir, "bayesant_", "testshort_aa_aln_label.txt"))

# Show accuracies 
output <- read.table(paste0(out_dir, "bayesant_", "test_aa_aln_label.txt"))
data_testshort_lab <- read.BayesANT.data(file.testshort, rank_names = rank_names, sep = "\\s+|\\|")
data_testshort_lab <- add_novelty_to_test_data(data_train, data_test_lab)
plot_accuracies(output, data_testshort_lab)

# Show accuracies in full case
output_short <- read.table(paste0(out_dir, "bayesant_", "testshort_aa_aln_label.txt"))
plot_accuracies(output_short, data_testshort_lab)


