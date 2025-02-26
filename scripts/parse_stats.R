#### convert *_time.tsv files to formatted table ####

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
    ncpu = ordered(ncpu, levels = c("1", "4", "16", "40", "gpu"))
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
      
      # BayesANT and IDTAXA are run as R scripts, and the file is passed as
      # an environmental variable rather than as a script. In these cases
      # we just rely on the order.
      model %in% c("bayesant", "idtaxa") ~ rep_len(c("train", "test", "testshort"), dplyr::n()),
    ) |>
      ordered(levels = c("train", "test", "testshort")),
    # determine the sequence type based on the command line
    seq_type = dplyr::case_when(
      grepl("(train|test|testshort)_aa", command) ~ "aa",
      grepl("(train|test|testshort)_nt", command) ~ "nt",
      
      # train_protax does not refer to the train file, but it only works on nt
      model == "protax-a" ~ "nt",
      
      # BayesANT and IDTAXA rely on ordering
      model %in% c("bayesant", "idtaxa") ~ rep(c("nt", "aa"), each = 3, length.out = dplyr::n()),
      TRUE ~ NA_character_
    ) |>
      ordered(levels = c("nt", "aa")),
    .by = c(marker, model, ncpu)
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
  time = round(lubridate::seconds_to_period(sum(lubridate::period_to_seconds(time))), 2),
  # max memory is the only thing that is relevant
  mem = max(as.integer(mem)),
  .by = c(marker, model, seq_type, stage, ncpu)
) |>
  dplyr::arrange(marker, model, seq_type, stage, ncpu)

knitr::kable(
  summary_stats,
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
  )
)
