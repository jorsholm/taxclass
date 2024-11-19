# Environment variables
nthread <- as.integer(Sys.getenv("OMP_NUM_THREADS"))
test_file <- Sys.getenv("TEST_FILE")
result_file <- Sys.getenv("RESULT_FILE")
path_to_rdp_jar <- Sys.getenv("PATH_TO_CLASSIFIER")
model_dir <- Sys.getenv("MODEL_DIR")

# Point to the RDP classifier
model_RDP <- list(dir = model_dir)
attr(model_RDP, "class") <- "RDPClassifier"

# Run the prediction
data_test <- readDNAStringSet(file = test_file)

result <- predict_RDP_parallel(model_RDP = model_RDP, 
                               path_to_rdp_jar = path_to_rdp_jar,
                               dataRDP_test = data_test,
                               cores = nthread)

# Save the results
write.table(result, file = result_file)
