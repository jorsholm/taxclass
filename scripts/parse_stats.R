#### convert *_time.tsv files to formatted table ####

# function to print a Period object in a condensed format
print_period <- function(x) {
  dplyr::case_when(
    x@day > 0 ~ sprintf("%dd %02d:%02d:%02d", x@day, x@hour, x@minute, lubridate::second(x)),
    x@hour > 0 ~ sprintf("%d:%02d:%02d", x@hour, x@minute, lubridate::second(x)),
    TRUE ~ sprintf("%d:%02d", x@minute, lubridate::second(x))
  )
}


# find all the time files
stats_files <- list.files(
  path = here::here(),
  pattern = "_time.tsv$",
  full.names = TRUE,
  recursive = TRUE
)

stats <- readr::read_tsv(stats_files, id = "file", col_types = "cccc") |>
  tidyr::extract(
    file,
    into = c("marker", "model", "ncpu"),
    regex = paste0(here::here(), "/(its|coi)/results/([^/]+)/.+_(\\d+|gpu)_time.tsv"),
    remove = FALSE
  ) |>
  dplyr::mutate(
    ncpu = ordered(ncpu, levels = c("1", "4", "16", "40", "gpu")),
    marker = toupper(marker)
  ) |>
  dplyr::mutate(
    
    # the time is given as MM:SS.sss for times less than 1 hour, or
    # HH:MM:SS for times greater than 1 hour
    time = dplyr::case_when(
      grepl("^\\d+:\\d+:\\d+$", time) ~ lubridate::hms(time),
      grepl("^\\d+:\\d+\\.\\d+$", time) ~ lubridate::ms(time),
    ),
    
    # determine the stage based on the command line
    stage = dplyr::case_when(
      # in most cases, for test we should be referring to one of the test or
      # testshort files; failing that we should refer to a train file
      grepl("testshort_(nt|aa)", command) ~ "testshort",
      grepl("test_(nt|aa)", command) ~ "test",
      grepl("train_(nt|aa)", command) ~ "train",
      
      # train_protax does not refer to the train file
      grepl("train_protax.sh", command) ~ "train",
      
      # BayesANT, IDTAXA, and Protax are run as R or BASH scripts, and the file 
      # is passed as an environmental variable rather than as part of the
      # command line. In these cases we just rely on the order.
      model %in% c("bayesant", "idtaxa", "protax") ~
        rep_len(c("train", "test", "testshort"), dplyr::n()),
      grepl("gappa prepare taxonomy-tree", command) ~ "train"
    ) |>
      ordered(levels = c("train", "test", "testshort")),
    
    # determine the sequence type based on the command line
    seq_type = dplyr::case_when(
      grepl("(train|test|testshort)_aa", command) ~ "aa",
      grepl("(train|test|testshort)_nt", command) ~ "nt",
      
      # train_protax does not refer to the train file, but it only works on nt
      model %in% c("protax", "protax-a") ~ "nt",
      
      # BayesANT and IDTAXA rely on ordering
      model %in% c("bayesant", "idtaxa") ~ rep(c("nt", "aa"), each = 3, length.out = dplyr::n()),
      
      # Building the taxonomic constraint tree is done just once for
      # epang-taxtree, but is used for both nt and aa.  Apply its time to both.
      grepl("gappa prepare taxonomy-tree", command) ~ "nt,aa",
      TRUE ~ NA_character_
    ),
    .by = c(marker, model, ncpu)
  ) |>
  # duplicate the rows which were used for both nt and aa
  tidyr::separate_rows(seq_type, sep = ",") |>
  dplyr::mutate(
    seq_type = ordered(seq_type, levels = c("nt", "aa"), labels = c("NT", "AA"))
  )

summary_stats <- dplyr::summarise(
  stats,
  # for cases where the same stage used more than one command, weight the CPU
  # usage of each command by the time it took to run
  CPU = (sum(readr::parse_number(CPU) * lubridate::period_to_seconds(time)) /
    sum(lubridate::period_to_seconds(time))) |>
    round() |>
    paste0("%"),
  # sum the time; just doing a sum on the period object does not seem to work.
  time = round(lubridate::seconds_to_period(sum(lubridate::period_to_seconds(time)))),
  # max memory is the only thing that is relevant
  mem = max(as.integer(mem)),
  .by = c(marker, model, seq_type, stage, ncpu)
) |>
  dplyr::arrange(marker, model, seq_type, stage, ncpu)

# load the model size data
# these have already been collected into a single file by `scripts/model_size.sh`

size_stats <- readr::read_delim(
  here::here("results/model_sizes.tsv"),
  delim = " ",
  col_names = c("model_size", "file"),
  col_types = "ic"
) |>
  tidyr::extract(
    file,
    into = c("marker", "model", "seq_type", "ncpu"),
    regex = paste0("(its|coi)/models/([^/]+)/.+_(nt|aa)_[^\\d]*(\\d+|gpu).*"),
    remove = FALSE
  ) |>
  dplyr::filter(!is.na(model)) |>
  dplyr::mutate(
    marker = toupper(marker),
    seq_type = toupper(seq_type),
    ncpu = ordered(ncpu, levels = c("1", "4", "16", "40", "gpu"))
  ) |>
  dplyr::summarize(
    model_size = sum(model_size),
    .by = c(marker, model, seq_type, ncpu)
  ) |>
  dplyr::arrange(marker, model, seq_type, ncpu)

# join the two tables
all_stats <- dplyr::full_join(
  summary_stats,
  size_stats,
  by = c("marker", "model", "seq_type", "ncpu")
)

# format the results as LaTeX

knitr::kable(
  all_stats,
  format = "latex",
  booktabs = TRUE,
  linesep = with(
    summary_stats,
    ifelse(
      marker == dplyr::lead(marker, default = NA) &
        model == dplyr::lead(model, default = NA) &
        seq_type == dplyr::lead(seq_type, default = NA),
      "",
      "\\addlinespace"
    )
  ),
  format.args = list(big.mark = "§")
) |>
  gsub("§", "\\,", x = _, fixed = TRUE)

# condensed for main text

maintext_stats <- 
  all_stats |>
  dplyr::mutate(
    time = print_period(time),
    marker = paste(marker, seq_type),
    model = dplyr::case_match(
      model,
      "bayesant" ~ "BayesANT",
      "blast" ~ "BLAST",
      "crest4" ~ "CREST4",
      "idtaxa" ~ "IDTAXA",
      "epang-freetree" ~ "EPA-ng free",
      "epang-taxtree" ~ "EPA-ng taxa",
      "epang-phyltree" ~ "EPA-ng phyl",
      "mycoai_bert" ~ "MycoAI-BERT",
      "mycoai_cnn" ~ "MycoAI-CNN",
      "protax" ~ "Protax",
      "protax-a" ~ "ProtaxA",
      "rdp_nbc" ~ "RDP-NBC",
      "sintax" ~ "SINTAX",
      .default = model
    ),
    ncpu = toupper(ncpu)
  ) |>
  tidyr::pivot_wider(
    names_from = stage,
    values_from = c(CPU, time, mem),
    names_vary = "slowest"
  ) |>
  tidyr::pivot_longer(
    cols = c(ends_with("test"), ends_with("testshort")),
    names_to = c(".value", "stage"),
    names_sep = "_"
  ) |>
  dplyr::mutate(
    dplyr::across(
      c(marker, model, ncpu, time_train, CPU_train, mem_train, model_size),
      ~ if (dplyr::cur_group()$stage == "testshort") ""
      else paste0("\\multirow{2}{*}[0.3em]{", if (is.numeric(.)) format(., big.mark = "§") else ., "}")
    ),
    mem = format(mem, big.mark = "§"),
    dplyr::across(c(time, CPU, mem), ~ if (dplyr::cur_group()$stage == "testshort") paste0("(", trimws(.), ")") else .),
    .by = stage
  ) |>
  dplyr::select(marker, model, ncpu, time_train, CPU_train, mem_train,
                time, CPU, mem, model_size)

knitr::kable(
  maintext_stats,
  format = "latex",
  booktabs = TRUE,
  linesep = with(
    maintext_stats,
    ifelse(
      marker == "",
      ifelse(
        dplyr::lag(
          marker == dplyr::lead(marker, default = NA, n = 2) &
            model == dplyr::lead(model, default = NA, n = 2)
        ),
        "",
        "\\addlinespace"
      ),
      "\\addlinespace[-0.6em]"
    )
  ),
  format.args = list(big.mark = "§")
) |>
  gsub("§", "\\,", x = _, fixed = TRUE) |>
  gsub("\\textbackslash{}", "\\", x = _, fixed = TRUE) |>
  gsub("\\{", "{", x = _, fixed = TRUE) |>
  gsub("\\}", "}", x = _, fixed = TRUE) |>
  sub(".*\\\\midrule\n", "", x = _) |>
  sub("\n\\\\bottomrule.*", "", x = _) |>
  clipr::write_clip()
