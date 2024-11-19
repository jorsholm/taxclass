# Load BayesANT (install it if not present)

# Replace with directory where the files are
# Load the environmental variables
train_file <- Sys.getenv("TRAIN_FILE")
model_dir <- Sys.getenv("MODEL_DIR")
path_to_rdp_jar <- Sys.getenv("PATH_TO_CLASSIFIER")

#path_to_rdp_jar <- "~/rdp_classifier/rdp_classifier_2.14/dist/classifier.jar"
#RDP_dir <- "~/taxclass/coi/models/rdp/RDP_nt_aln_label"

# Source useful functions
source("~/taxclass/coi/models/rdp/rRDP.R")

# Load the data
dataRDP_train <- readDNAStringSet(file = train_file)
# Modify the names to make the annotation coherent with the model
names_old <- names(dataRDP_train)
names_new <- gsub("[|]", ";", sub(" ", " Root;", names_old))
names(dataRDP_train) <- names_new

# Run the RDP classifier. This code automatically creates a repository with the 
# model, no need to save the output in RDS.
model_RDP <- train_RDP(dataRDP_train,
                       path_to_rdp_jar = path_to_rdp_jar, 
                       dir = RDP_dir, 
                       java_args = "-Xmx5g",
                       rank = "species")


