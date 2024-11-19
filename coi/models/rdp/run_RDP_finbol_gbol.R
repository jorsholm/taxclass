# This file runs the  RDP classifier on the finbol-gbol coi sequences
library(tidyverse)
library(seqinr)

# Source the file to run the RDP classifier
source("~/taxclass/coi/models/rdp/rRDP.R")

# Path to the rdp classifier in java code
path_to_rdp_jar <- "~/rdp_classifier/rdp_classifier_2.14/dist/classifier.jar"

################################################################################
# Part 1 - train and test using the nt file with labels
################################################################################

RDP_dir <- "~/taxclass/coi/models/rdp/RDP_nt_aln_label"

file.train <- "~/taxclass/coi/data/train_nt_aln_label.fasta"
  
dataRDP_train <- readDNAStringSet(file = file.train)
# Modify the names to make the annotation coherent with the model
names_old <- names(dataRDP_train)
names_new <- gsub("[|]", ";", sub(" ", " Root;", names_old))
names(dataRDP_train) <- names_new

# Run the RDP classifier
model_RDP <- train_RDP(dataRDP_train,
                       path_to_rdp_jar = path_to_rdp_jar, 
                       dir = RDP_dir, 
                       java_args = "-Xmx5g",
                       rank = "species")


# Predict RDP 

#--- test
file.test <- "~/taxclass/coi/data/test_nt_aln_label.fasta"
data_test <- readDNAStringSet(file = file.test)
output <- predict_RDP_parallel(model_RDP = model_RDP, 
                               path_to_rdp_jar = path_to_rdp_jar,
                               dataRDP_test = data_test,
                               cores = 22)

#--- testshort
file.testshort <- "~/taxclass/coi/data/testshort_nt_aln_label.fasta"
data_testshort <- readDNAStringSet(file = file.test)
output_short <- predict_RDP_parallel(model_RDP = model_RDP, 
                                     path_to_rdp_jar = path_to_rdp_jar,
                                     dataRDP_test = data_testshort,
                                     cores = 22)

# Save output
write.table(output, 
            file = "~/taxclass/coi/results/rdp/rdp_test_nt_aln_label.txt")
write.table(output_short, 
            file = "~/taxclass/coi/results/rdp/rdp_testshort_nt_aln_label.txt")



