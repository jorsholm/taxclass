###############################################################
# This file runs the RDP classifier from R. The structure and some functions are
# similar to the package "rRDP". However, this package uses an old version of the
# classifier (RDP version 2.05). Here we create
# functions that work with the latest version of RDP (version 2.13).
# In the folder "java", we include the latest version of the RDP classifier,
# downloaded from https://sourceforge.net/projects/rdp-classifier/
###############################################################

## Create the rdp attribute to the folder
rdp <- function(dir = NULL) {
  structure(list(dir = dir), class = "RDPClassifier")
}

## Call the java
.javaExecutable <- function() Sys.which("java")

## Call the RDP classifier
.get_rdp <- function(path_to_rdp_jar) {
  file.path(path_to_rdp_jar) # <- Place here the RDP classifier directory
}

## Train the RDP classifier
train_RDP <- function(x, path_to_rdp_jar, dir = "classifier", rank = "genus", java_args = "-Xmx5g") {
  # if (file.exists(dir)) stop("Classifier directory already exists! Choose a different directory")

  taxonomyNames <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  rankNumber <- which(tolower(rank) == tolower(taxonomyNames))
  dir.create(dir, showWarnings = FALSE)
  l <- strsplit(names(x), "Root;")
  # annot is the hierarchy starting from kingdom down to genus (or desired rank)
  annot <- sapply(l, FUN = function(x) x[2])
  # Make names of x so that it only has hierarchy info upto the desired rank
  names(x) <- sapply(strsplit(names(x), ";"), FUN = function(y) {
    paste(y[1:(rankNumber + 1)], collapse = ";")
  })
  # get the indices of sequences that are to be removed because they have incomplete hierarchy information
  removeIdx <- as.integer(sapply(annot, FUN = function(y) {
    if (length(unlist(strsplit(y, ";"))) < rankNumber || grepl(";;", y) ||
      grepl("unknown", y)) {
      1
    } else {
      0
    }
  }))

  if (any(removeIdx == 1)) {
    removeIdx <- which(removeIdx == 1)
    # get greengenes ids to be removed
    idsRemoved <- as.character(sapply(names(x),
      FUN = function(y) as.character(unlist(strsplit(y, " "))[1])
    ))[removeIdx]
    x <- x[-removeIdx]
    annot <- annot[-removeIdx]
    cat("Warning! The following sequences did not contain complete hierarchy information and have been ignored:", idsRemoved, "\n")
  }
  if (length(x) <= 0) stop("No sequences with complete information found")
  # writeXStringSet(x,file.path(dir,"train.fasta"))
  # h<-matrix(ncol=rankNumber,nrow=0)
  # colnames(h) <-c("Kingdom","Phylum","Class","Order","Family","Genus")
  # colnames(h) <- taxonomyNames[1:rankNumber]
  # for(i in 1:length(annot)) {
  #  h<-rbind(h,unlist(strsplit(annot[i],";"))[1:rankNumber])
  # }
  h <- str_split(annot, ";", simplify = TRUE)[, 1:rankNumber]
  colnames(h) <- taxonomyNames[1:rankNumber]
  m <- matrix(ncol = 5, nrow = 0)
  # first row of the file
  f <- "0*Root*-1*0*rootrank"
  m <- rbind(m, unlist(strsplit(f, split = "\\*")))
  taxNames <- colnames(h)
  badHierarchy <- vector()
  cat("Building the hierarchy \n")
  pb <- txtProgressBar(style = 3)
  for (i in 1:nrow(h)) {
    # if(i%%100 == 0) print(i)
    setTxtProgressBar(pb, i / nrow(h))
    for (j in 1:ncol(h)) {
      taxId <- nrow(m)
      taxonName <- h[i, j]
      if (j == 1) {
        parentTaxId <- 0
      } else {
        parentTaxId <- previousTaxId
      }
      depth <- j
      rank <- colnames(h)[j]

      # search if already there
      if (length(which(m[, 2] == taxonName & m[, 5] == rank)) == 0) {
        str <- paste(taxId, taxonName, parentTaxId, depth, rank, sep = "*")
        m <- rbind(m, unlist(strsplit(str, split = "\\*")))
        previousTaxId <- taxId
      } else if (length(which(m[, 2] == taxonName & m[, 5] == rank)) > 0) {
        row <- which(m[, 2] == taxonName & m[, 5] == rank)
        # error check -> is parent tax id same for both or not?
        # if not, remove the sequence
        if (parentTaxId != m[row, 3]) {
          # x<-x[-i]
          badHierarchy <- c(badHierarchy, i)
        }
        previousTaxId <- m[row, 1]
      }
      # end seach
    }
  }

  out <- apply(m, MARGIN = 1, FUN = function(x) paste(x, collapse = "*"))
  write(out, file = file.path(dir, "train.txt"))

  # create parsed training files
  if (length(badHierarchy) > 0) {
    warning(
      "Following sequences had bad sequence hierarchy information, so removing: ",
      names(x)[badHierarchy], "\n"
    )
    x <- x[-badHierarchy]
  }
  writeXStringSet(x, file.path(dir, "train.fasta"))

  ## Construct the RDP classifier
  cat(" \nRunning RDP via Java \n")
  java_script <- paste(
    .javaExecutable(), java_args, "-jar", .get_rdp(path_to_rdp_jar), "train",
    "-o", dir, "-t", file.path(dir, "train.txt"),
    "-s", file.path(dir, "train.fasta")
  )
  system(java_script, ignore.stderr = F, ignore.stdout = F, intern = T)

  # file.copy(system.file("java/rRNAClassifier.properties",
  #                      package="rRDP"), dir)
  file.copy("~/BayesANT/Files_for_Zenodo_submission/RDP/java/rRNAClassifier.properties", dir) # <- Place here the RDP classifier directory


  rdp(dir)
}


## Predict!
predict_RDP <- function(object,
                        newdata,
                        path_to_rdp_jar,
                        confidence = .8,
                        rdp_args = "",
                        java_args = "-Xmx5g") {
  # RDPmodel: an object of class RDP.
  # newdata: DNA sequence set

  classifier <- object$dir
  x <- newdata

  # get temp files and change working directory
  wd <- tempdir()
  dir <- getwd()
  temp_file <- tempfile(tmpdir = wd)
  #  temp_file <- "newdata"
  on.exit({
    file.remove(Sys.glob(paste(temp_file, "*", sep = "")))
    setwd(dir)
  })

  ## Files o classify
  infile <- paste(temp_file, ".fasta", sep = "")
  outfile <- paste(temp_file, "_tax_assignments.txt", sep = "")

  # Property
  property <- paste("-t", file.path(classifier, "rRNAClassifier.properties"))

  ## RUN THE RDP model calling Java from R. Notice that the format that we want is "-f allrank".
  ## In the package rRDP instead we have "-f fixrank", which causes problems and is not rigid.
  javascript_code <- paste(
    .javaExecutable(), java_args, "-jar",
    .get_rdp(path_to_rdp_jar), "classify",
    property, "-o", outfile, "-c",
    confidence, "-format allrank", rdp_args, infile
  )
  writeXStringSet(x, infile, append = FALSE)
  system(javascript_code, ignore.stdout = TRUE, ignore.stderr = TRUE, intern = FALSE)

  ## read and parse rdp output
  cl_tab <- read.table(outfile, sep = "\t")

  seq_names <- cl_tab[, 1] ## sequence names are in first column

  i <- seq(3, ncol(cl_tab), by = 3) ## start at col. 3 and use 3 columns for each tax. level

  ## get classification
  cl <- cl_tab[, i]
  dimnames(cl) <- list(seq_names, as.matrix(cl_tab[1, i + 1])[1, ])

  ## get confidence
  conf <- as.matrix(cl_tab[, i + 2])
  dimnames(conf) <- list(seq_names, paste0("Prob_", as.matrix(cl_tab[1, i + 1])[1, ]))

  ### don't show assignments with too low confidence
  # if(confidence>0) cl[conf < confidence] <- NA

  cbind(cl[, -1], conf[, -1])
}


predict_RDP_parallel <- function(model_RDP, dataRDP_test, path_to_rdp_jar,  cores = 1) {
  doParallel::registerDoParallel(cores)
  d <- c(1:length(dataRDP_test))
  aa <- base::split(d, sort(d %% cores))
  out_RDP <- foreach(i = 1:length(aa), .combine = "rbind") %dopar% {
    predict_RDP(model_RDP, 
                path_to_rdp_jar = path_to_rdp_jar, 
                newdata = dataRDP_test[aa[[i]]],
                java_args = "-Xmx2g")
  }
  rownames(out_RDP) <- gsub(" .*", "", names(dataRDP_test))
  return(out_RDP)
}
