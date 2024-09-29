# input and output file names are set in main script using environmental
# variables
train_file <- Sys.getenv("TRAIN_FILE")
tax_file <- Sys.getenv("TAX_FILE")
model_file <- Sys.getenv("MODEL_FILE")

# Read the training data.  We will be doing AA sometimes and DNA sometimes,
# so read first as BString, then try converting to DNA (which has a more
# restricted alphabet than AA) and if that fails convert to AA.
train_data <- Biostrings::readBStringSet(train_file)
train_data <- tryCatch(
    Biostrings::DNAStringSet(train_data),
    error = \(e) Biostrings::AAStringSet(train_data)
)

# Read the training taxonomy in tabular format.
train_tax <- read.delim(
    Sys.getenv("TAX_FILE"),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

# DECIPHER wants the taxonomy as single strings, semicolon delimited, beginning
# with "Root"
train_tax$id <- NULL
train_tax$sep <- ";"
train_tax <- do.call(paste, train_tax)
train_tax <- paste("Root", train_tax, sep = ";")

# Train the model.
model <- DECIPHER::LearnTaxa(train = train_data, taxonomy = train_tax)

# Write the output
saveRDS(model, model_file)
