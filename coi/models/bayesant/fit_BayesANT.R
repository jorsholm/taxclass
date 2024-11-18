# This file runs BayesANT for the finbol-gbol split (new version)

# Load BayesANT (install it if not present)
if (!"BayesANT" %in% rownames(installed.packages())) {
  devtools::install_github("alessandrozito/BayesANT")
}

# Replace with directory where the files are
# Load the environmental variables
train_file <- Sys.getenv("TRAIN_FILE")
model_file <- Sys.getenv("MODEL_FILE")

# train_file <- "coi/data/train_nt_aln_label.fasta"
# model_file <- "coi/models/bayesant/all_models/BayesANT_nt_aln_label.rds.gzip"
data_train <- BayesANT::read.BayesANT.data(train_file, sep = "\\s+|\\|")

# Run the model
if (grepl("aa", train_file)) {     # <---- Check for aminoacid sequences
  model <- BayesANT::BayesANT_amino(
    data = data_train,
    typeseq = "aligned",
    type_location = "single",
    newtaxa = TRUE,
    verbose = TRUE, 
    amino = TRUE,
  )
} else {           # <---- run the standard nucleotide kernel
  model <- BayesANT::BayesANT(
    data = data_train,
    typeseq = "aligned",
    type_location = "single",
    newtaxa = TRUE,
    verbose = TRUE
  )
}

# Write the output
saveRDS(model, model_file) # <----- suggestion: do `compress = "gzip"` to save space
