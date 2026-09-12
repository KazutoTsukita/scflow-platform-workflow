#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})


# コマンドライン引数の取得
args <- commandArgs(trailingOnly = TRUE)

# 引数が指定されているかチェック
if (length(args) == 0) {
  stop("No file path provided", call. = FALSE)
}

# ファイルパス
file_path <- args[1]
output_tsv <- file_path
output_csv <- NA_character_

#フィルタリング条件
filter_conditions <- args[-1]
output_arg <- stringr::str_detect(filter_conditions, "^output_tsv=")
if (any(output_arg)) {
  output_tsv <- stringr::str_replace(filter_conditions[which(output_arg)[1]], "^output_tsv=", "")
  filter_conditions <- filter_conditions[!output_arg]
}
output_csv_arg <- stringr::str_detect(filter_conditions, "^output_csv=")
if (any(output_csv_arg)) {
  output_csv <- stringr::str_replace(filter_conditions[which(output_csv_arg)[1]], "^output_csv=", "")
  filter_conditions <- filter_conditions[!output_csv_arg]
}


# データの読み込み
a <- read_tsv(file_path, col_types = cols())

empty_to_na <- function(x) {
  x <- as.character(x)
  x[stringr::str_trim(x) == ""] <- NA_character_
  x
}

get_col_or_empty <- function(df, col) {
  if (col %in% names(df)) {
    as.character(df[[col]])
  } else {
    rep("", nrow(df))
  }
}

extract_gsm <- function(x) {
  stringr::str_to_upper(stringr::str_extract(as.character(x), regex("GSM[0-9]+", ignore_case = TRUE)))
}

resolve_sample_alias <- function(df) {
  sample_alias <- empty_to_na(get_col_or_empty(df, "sample_alias"))
  sample_alias_gsm <- empty_to_na(extract_gsm(sample_alias))
  experiment_alias_gsm <- empty_to_na(extract_gsm(get_col_or_empty(df, "experiment_alias")))
  sample_title_gsm <- empty_to_na(extract_gsm(get_col_or_empty(df, "sample_title")))
  sample_accession <- empty_to_na(get_col_or_empty(df, "sample_accession"))

  # Prefer GEO sample IDs when available. ENA sample_alias can contain SAMN/SRS IDs,
  # while the GEO GSM can appear in experiment_alias or sample_title instead.
  dplyr::coalesce(sample_alias_gsm, experiment_alias_gsm, sample_title_gsm, sample_alias, sample_accession)
}

a <- a %>% dplyr::mutate(.uniscflow_resolved_sample_alias = resolve_sample_alias(.))

is_accession_key <- function(key) {
  stringr::str_detect(key, "(^|_)accession$") || key %in% c("study_alias", "experiment_alias")
}

exact_match <- function(x, value) {
  stringr::str_to_upper(stringr::str_trim(dplyr::coalesce(as.character(x), ""))) ==
    stringr::str_to_upper(stringr::str_trim(value))
}

selected <- rep(TRUE, nrow(a))
for (condition in filter_conditions) {
  key_value <- str_split(condition, "=", n = 2, simplify = TRUE)
  key <- key_value[1]
  if (ncol(key_value) < 2 || key == "" || key_value[2] == "") {
    stop(sprintf("Invalid filter condition '%s'; expected key=value", condition), call. = FALSE)
  }
  values <- str_split(key_value[2], ",", simplify = FALSE)[[1]] %>%
    stringr::str_trim()
  values <- values[values != ""]
  if (length(values) == 0) {
    stop(sprintf("Filter condition '%s' has no non-empty values", condition), call. = FALSE)
  }

  if (key == "sample_alias") {
    candidate_columns <- intersect(
      c("sample_alias", ".uniscflow_resolved_sample_alias", "experiment_alias", "sample_title", "sample_accession", "secondary_sample_accession"),
      names(a)
    )
    condition_matched <- rep(FALSE, nrow(a))
    for(iter in values){
      matched <- rep(FALSE, nrow(a))
      for (column in candidate_columns) {
        column_values <- dplyr::coalesce(as.character(a[[column]]), "")
        matched <- matched |
          exact_match(column_values, iter) |
          exact_match(dplyr::coalesce(extract_gsm(column_values), ""), iter)
      }
      if (!any(matched)) {
        warning(sprintf("No rows matched sample_alias=%s across sample_alias/experiment_alias/sample metadata", iter))
      }
      condition_matched <- condition_matched | matched
    }
  } else if (key %in% names(a)) {
    condition_matched <- rep(FALSE, nrow(a))
    for(iter in values){
      column_values <- dplyr::coalesce(as.character(a[[key]]), "")
      if (is_accession_key(key)) {
        matched <- exact_match(column_values, iter)
      } else {
        matched <- stringr::str_detect(column_values, fixed(iter, ignore_case = TRUE))
      }
      if (!any(matched)) {
        warning(sprintf("No rows matched %s=%s", key, iter))
      }
      condition_matched <- condition_matched | matched
    }
  } else {
    stop(sprintf("Column name '%s' does not exist in the data frame.", key), call. = FALSE)
  }
  selected <- selected & condition_matched
}

if(length(filter_conditions) >= 1) {
  a <- a[selected, , drop = FALSE]
}
raw <- a


# 必要な列の選択
required_columns <- c("study_alias", "sample_accession", "sample_alias", ".uniscflow_resolved_sample_alias", "run_accession", "experiment_title", "study_title")
missing_columns <- setdiff(required_columns,colnames(a))

if(length(missing_columns) >=1){
  missing_columns_dataframe <- c()
  for(s in 1:length(missing_columns)){
    tibble("V1" = rep("none",nrow(a))) -> temp1
    colnames(temp1) <- missing_columns[s]
    missing_columns_dataframe <- missing_columns_dataframe %>% dplyr::bind_cols(temp1)
  }
  a <- a %>% dplyr::bind_cols(missing_columns_dataframe) %>%
  dplyr::select(all_of(required_columns))
}

a <- a %>% dplyr::select(all_of(required_columns))

  # データの集約
a <- a %>%
  dplyr::mutate(sample_alias = dplyr::coalesce(empty_to_na(.uniscflow_resolved_sample_alias), empty_to_na(sample_alias), empty_to_na(sample_accession))) %>%
  dplyr::select(-.uniscflow_resolved_sample_alias) %>%
  group_by(study_alias,sample_accession,sample_alias) %>%
  summarise(
    run_accessions = paste(run_accession, collapse = " "),
    experiment_title = first(experiment_title),
    study_title = first(study_title),
    .groups = "drop"
  ) %>%
  unique()


a$sequencing<- tibble(title=a$experiment_title) %>% apply(MARGIN=1,FUN=function(vector){
  vector %>% stringr::str_split(pattern = ": ") -> temp1
  temp1[[1]][1] %>% return()
})

a$sampleinfo<- tibble(title=a$experiment_title) %>% apply(MARGIN=1,FUN=function(vector){
  vector %>% stringr::str_split(pattern = ": ") -> temp1
  temp1[[1]][3] %>% return()
})

a %>% dplyr::select(-experiment_title) -> a

dir_path <- dirname(file_path)
pattern <- "PRJNA[0-9]+"
filename <- stringr::str_extract(file_path, pattern=pattern)
if (is.na(output_csv)) {
  output_csv <- paste0(dir_path,"/",filename,".csv")
}

a <- a %>%
  dplyr::mutate(PRJNA=filename) %>%
  dplyr::select(PRJNA, study_alias, sample_accession, sample_alias, sampleinfo, everything()) %>%
  dplyr::mutate(sample_alias = ifelse(is.na(sample_alias), sample_accession, sample_alias)) %>% ungroup()


write_csv(a, output_csv)
write_tsv(raw, output_tsv)
