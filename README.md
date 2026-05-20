# taxclass

## Repository structure

This repository is organized around two marker datasets: COI and ITS. Each marker has a parallel directory structure containing raw data preparation notes, formatted datasets, model training/running files, and classifier outputs.

### Top-level files

- `functions.R`  
  Shared R helper functions used by the result-reading, evaluation, and plotting scripts.

- `read_result_coi.R` and `read_result_its.R`  
  Scripts for reading classifier-specific output files and converting them into a common format for COI and ITS, respectively.

- `read_results.qmd`  
  Quarto file that runs the result-reading scripts for the main combinations of marker, sequence length, datatype, and NA-handling settings.

- `evaluate_jointly.R`  
  Main downstream evaluation script for comparing classification results across methods and datasets.

- `aa_comparison.R`  
  Script for comparing amino-acid and nucleotide-based analyses for COI.

- `overlap_plot.R`  
  Script for plotting train/test overlap-related summaries.

- `taxclass.Rproj`  
  RStudio project file.

### Marker-specific directories

#### `coi/`

Contains the COI analysis.

- `coi/data_raw/`  
  Raw COI input files and documentation for downloading and processing the original data sources. The README in this folder documents the steps used to generate the COI training and test datasets from BOLD, GBOL, and related reference files.

- `coi/data/`  
  Formatted COI training and test datasets used by the classifiers. This includes nucleotide and amino-acid FASTA files, aligned and unaligned versions, labelled and unlabelled files, shortened test sequences, and taxonomy mapping files in method-specific formats.

- `coi/models/`  
  Method-specific scripts and files for training or preparing COI classifiers. Subdirectories correspond to the different classification approaches, including BayesANT, BLAST, CREST4, Dnabarcoder, EPA-ng-based workflows, IDTAXA, MycoAI, PROTAX-A, RDP NBC, and SINTAX.

- `coi/results/`  
  Method-specific COI classifier outputs. The subdirectory names mirror the classifier names used in `coi/models/`.

- `coi/short_stats.R`  
  Script for summarizing properties of the shortened COI test sequences.

#### `its/`

Contains the ITS analysis.

- `its/data_raw/`  
  Raw ITS input data and documentation for downloading and processing the original data sources. The README in this folder documents the steps used to generate the ITS training and test datasets from UNITE and Westerdijk/CBG datasets.

- `its/data/`  
  Formatted ITS training and test datasets used by the classifiers. This includes FASTA files, labelled and unlabelled files, shortened ITS2-only test sequences, and taxonomy mapping files in method-specific formats.

- `its/models/`  
  Method-specific scripts and files for training or preparing ITS classifiers. Subdirectories correspond to the different classification approaches, including BayesANT, BLAST, CREST4, Dnabarcoder, IDTAXA, MycoAI, PROTAX, RDP NBC, and SINTAX.

- `its/results/`  
  Method-specific ITS classifier outputs. The subdirectory names mirror the classifier names used in `its/models/`.

### Shared output and utility directories

- `results/`  
  Repository-level summary outputs, including model size summaries.

- `scripts/`  
  Utility scripts for extracting, parsing, and summarizing runtime/resource statistics and model sizes.
