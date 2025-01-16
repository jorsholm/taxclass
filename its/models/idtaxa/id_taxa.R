# input and output file names are set in main script using environmental
# variables
nthread <- as.integer(Sys.getenv("OMP_NUM_THREADS"))
model_file <- Sys.getenv("MODEL_FILE")
test_file <- Sys.getenv("TEST_FILE")
result_file <- Sys.getenv("RESULT_FILE")

# Read the model
model <- readRDS(model_file)

# Read the test data.  We will be doing AA sometimes and DNA sometimes,
# so read first as BString, then try converting to DNA (which has a more
# restricted alphabet than AA) and if that fails convert to AA.
test_data <- Biostrings::readBStringSet(test_file)
test_data <- tryCatch(
    Biostrings::DNAStringSet(test_data),
    error = \(e) Biostrings::AAStringSet(test_data)
)

#### Classify the test data ####
result <- DECIPHER::IdTaxa(
    test = test_data,
    trainingSet = model,
    type = "collapsed", # this output is faster to output and parse
    strand = "top", # test sequences are oriented properly
    threshold = 1, # Also report low-confidence results
    processors = nthread
)

#### format results ####

# remove "Root"
result <- sub("Root[^;]+; *", "", result)

# parse 2-digit percentage bootstraps
result <- gsub(
    pattern = " *\\[([^.\\d]*)(\\d{2})\\.(\\d*)%\\];? *",
    replacement = "\t0.\\2\\3\t",
    result,
    perl = TRUE
)

# parse 1-digit percentage bootstraps
result <- gsub(
    pattern = " *\\[([^.\\d]*)(\\d)\\.(\\d*)%\\];? *",
    replacement = "\t0.0\\2\\3\t",
    result,
    perl = TRUE
)

# parse 100% bootstraps
result <- gsub(
    pattern = " *\\[([^.\\d]*)100\\.0%\\];? *",
    replacement = "\t1.0\t",
    result,
    perl = TRUE
)

# remove trailing white space
result <- trimws(result)

# add IDs
result <- paste(names(result), result, sep = "\t")
header <- paste("ID", "class", "Prob_class", "order", "Prob_order", "family",
               "Prob_family", "subfamily", "Prob_subfamily", "tribe",
               "Prob_tribe", "genus", "Prob_genus", "species", "Prob_species",
               sep = "\t")

writeLines(c(header, result), result_file)
