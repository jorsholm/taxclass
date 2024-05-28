
#' Loads FinBOL data (train) and GBOL data (test) as listed dataframes.
#' @param level lowest taxonomic level (genus/species)
load_FinBOL_GBOL <- function(level = "species"){

  data_train <- seqinr::read.fasta(paste0("data/finbol-gbol/train_finbol-gbol_",
                                          level, ".fasta"),
                                   whole.header = T)
  data_test <- seqinr::read.fasta(paste0("data/finbol-gbol/test_finbol-gbol_",
                                         level, ".fasta"),
                                  whole.header = T)

  ref_seq <- data.frame(do.call("rbind", stringr::str_split(names(data_train), " |[|]")),
                        unlist(unname(lapply(data_train,
                                             function(x) stringr::str_to_upper(
                                               stringr::str_c(x, sep = "", collapse = ""))))))
  ref_seq_test <- data.frame(do.call("rbind", stringr::str_split(names(data_test), " |[|]")),
                             unlist(unname(lapply(data_test,
                                                  function(x) stringr::str_to_upper(
                                                    stringr::str_c(x, sep = "", collapse = ""))))))

  colnames(ref_seq) = colnames(ref_seq_test) = c("ID", "Class", "Order", "Family", "Subfamily", "Tribe", "Genus", "Species", "DNA")

  return(list("train" = ref_seq, "test"= ref_seq_test))
}
