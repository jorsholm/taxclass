# input and output file names are set in main script using environmental
# variables
nthread <- as.integer(Sys.getenv("OMP_NUM_THREADS"))
model_file <- Sys.getenv("MODEL_FILE")
test_file <- Sys.getenv("TEST_FILE")
result_file <- Sys.getenv("RESULT_FILE")

# Read the model
model <- readRDS(model_file)

#test_file <- "coi/data/test_nt_aln_label.fasta"
data_test <- BayesANT::read.BayesANT.testDNA(test_file, sep = "\\s+|\\|")

# Predict nt
if (grepl("aa", test_file)) {
  rho <- 0.3
} else {
  rho <- 0.1
}

# Make prediction
result <- predict(model,
  DNA = data_test,
  rho = rho,
  cores = nthread
)

# Save the output. Need relabel according to the order needed for the other files
#result_file <- "coi/results/bayesant/bayesant_test_nt_aln_label.txt"
write.table(result, file = result_file)



