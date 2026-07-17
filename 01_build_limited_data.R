rm(list = ls())
options(scipen = 999)

library(tidyverse)
library(haven)
library(psych)
library(janitor)

# PATHS --------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grepl(file_arg, args)])
script_path <- gsub("~\\+~", " ", script_path)

if (length(script_path) == 0) {
  script_dir <- getwd()
} else {
  script_dir <- dirname(normalizePath(script_path))
}

data_path <- file.path(
  script_dir,
  "..",
  "..",
  "AU",
  "DATA",
  "2025_backslide_AU_November 5, 2025_12.34.sav"
)

csv_out <- file.path(script_dir, "au_replication_limited_data.csv")
rds_out <- file.path(script_dir, "au_replication_limited_data.rds")

# HELPERS ------------------------------------------------------------------

add_factor_scores <- function(data, score_name) {
  cm_items <- c("CM_1", "CM_2", "CM_3", "CM_4", "CM_5")

  fa_data <- data %>%
    select(all_of(cm_items)) %>%
    drop_na()

  fa_result <- psych::fa(fa_data, nfactors = 1, fm = "ml", scores = "regression")
  scores <- as.numeric(fa_result$scores)

  # Keep the factor direction intuitive: higher values mean stronger
  # conspiracy mentality.
  if (mean(as.numeric(fa_result$loadings[, 1]), na.rm = TRUE) < 0) {
    scores <- -scores
  }

  data_indexed <- data %>%
    mutate(.row_id = row_number())

  scored_rows <- tibble(
    .row_id = data_indexed$.row_id[complete.cases(data_indexed[, cm_items])],
    score = scores
  )

  data_indexed %>%
    left_join(scored_rows, by = ".row_id") %>%
    rename(!!score_name := score) %>%
    select(-.row_id)
}

# READ DATA ----------------------------------------------------------------

raw_data <- read_sav(data_path)

numeric_vars <- c(
  "Q2", "Q3", "Q4", "Q5", "Q8", "Q9_partyID",
  "Q13", "Q14", "Q15", "Q16",
  "Q19_1", "Q19_2", "Q19_3", "Q19_4", "Q19_5",
  "Q20", "Q21", "Q23", "Q24",
  "Q25_1", "Q25_3", "Q25_4",
  "Q26_1", "Q26_2", "Q26_3", "Q26_5",
  "Q27", "Q28", "CM_1", "CM_2", "CM_3", "CM_4", "CM_5"
)

# RECODE -------------------------------------------------------------------

data_au <- raw_data %>%
  filter(T30_Page_Submit != "") %>%
  mutate(across(any_of(numeric_vars), ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(
    state = State,

    sex = case_when(
      Q4 == 1 ~ "male",
      Q4 == 2 ~ "female",
      TRUE ~ NA_character_
    ),
    sex = factor(sex, levels = c("female", "male")),
    sex_num = ifelse(sex == "male", 1, 0),

    age = Q5,

    edu_raw = case_when(
      Q27 == 1 ~ 1,
      Q27 == 2 ~ 2,
      Q27 == 4 ~ 3,
      Q27 == 3 ~ 4,
      TRUE ~ NA_real_
    ),
    edu = (edu_raw - 1) / 2,

    hincome = (Q28 - 1) / 14,

    ideology = case_when(
      Q20 == 99 ~ NA_real_,
      TRUE ~ (Q20 - 1) / 6
    ),

    party_bloc = case_when(
      Q9_partyID == 1 ~ "ruling",
      Q9_partyID %in% c(2, 3, 4) ~ "opposition",
      TRUE ~ "Independent"
    ),
    party_bloc = factor(party_bloc, levels = c("ruling", "opposition", "Independent")),

    party_id_detail = case_when(
      Q9_partyID == 1 ~ "Labor",
      Q9_partyID %in% c(2, 3, 4) ~ "Coalition",
      Q9_partyID == 7 ~ "Greens",
      Q9_partyID == 8 ~ "One Nation",
      Q9_partyID %in% c(9, 99, 100) ~ "Independent",
      TRUE ~ "No party ID"
    ),
    party_id_detail = factor(
      party_id_detail,
      levels = c("Coalition", "Labor", "Greens", "One Nation", "Independent", "No party ID")
    ),

    major_party_identifier = ifelse(Q9_partyID %in% c(1, 2, 3, 4), 1, 0),

    candidate_party = case_when(
      ruling == "0" ~ "opposition",
      ruling == "1" ~ "ruling",
      TRUE ~ NA_character_
    ),
    candidate_party = factor(candidate_party, levels = c("ruling", "opposition")),

    ingroup3 = case_when(
      party_bloc == "ruling" & candidate_party == "ruling" ~ "In-group",
      party_bloc == "opposition" & candidate_party == "opposition" ~ "In-group",
      party_bloc %in% c("ruling", "opposition") ~ "Out-group",
      TRUE ~ "Independent"
    ),
    ingroup3 = factor(ingroup3, levels = c("In-group", "Out-group", "Independent")),

    ft_rul = Q19_1 / 100,
    ft_opp = rowMeans(select(., Q19_2, Q19_3, Q19_4, Q19_5), na.rm = TRUE) / 100,

    vignette = case_when(
      vignette == 0 ~ "Control",
      vignette == 1 ~ "Election",
      vignette == 2 ~ "Protests",
      vignette == 3 ~ "Court",
      TRUE ~ NA_character_
    ),
    vignette = factor(vignette, levels = c("Control", "Election", "Protests", "Court")),

    attention_pass = ifelse(Q2 == 1, 1, 0),
    issue_check_pass = case_when(
      vignette == "Control" & Q23 == 99 ~ 1,
      vignette == "Election" & Q23 == 1 ~ 1,
      vignette == "Protests" & Q23 == 2 ~ 1,
      vignette == "Court" & Q23 == 3 ~ 1,
      TRUE ~ 0
    ),
    main_analytic = issue_check_pass,

    intent = Q24 - 1,

    enthusiasm = (Q25_4 - 1) / 4,
    negative = ((Q25_1 - 1) / 4 + (Q25_3 - 1) / 4) / 2,

    competent = (5 - Q26_1) / 4,
    sincere = (Q26_2 - 1) / 4,
    trustworthy = (5 - Q26_3) / 4,
    undemocratic_perception = (Q26_5 - 1) / 4,

    CM_1_bin = case_when(is.na(CM_1) ~ NA_integer_, CM_1 >= 16 ~ 1L, TRUE ~ 0L),
    CM_2_bin = case_when(is.na(CM_2) ~ NA_integer_, CM_2 >= 4 ~ 1L, TRUE ~ 0L),
    CM_3_bin = case_when(is.na(CM_3) ~ NA_integer_, CM_3 >= 4 ~ 1L, TRUE ~ 0L),
    CM_4_bin = case_when(is.na(CM_4) ~ NA_integer_, CM_4 >= 4 ~ 1L, TRUE ~ 0L),
    CM_5_bin = case_when(is.na(CM_5) ~ NA_integer_, CM_5 >= 4 ~ 1L, TRUE ~ 0L),
    conspiracy_additive = CM_1_bin + CM_2_bin + CM_3_bin + CM_4_bin + CM_5_bin
  )

# CONSPIRACY SCALE ---------------------------------------------------------

data_au <- data_au %>%
  add_factor_scores("conspiracy_full") %>%
  mutate(case_id = row_number())

analytic_scores <- data_au %>%
  filter(main_analytic == 1) %>%
  add_factor_scores("conspiracy") %>%
  select(case_id, conspiracy)

data_au <- data_au %>%
  left_join(analytic_scores, by = "case_id") %>%
  mutate(conspiracy_z = as.numeric(scale(conspiracy)))

# KEEP LIMITED VARIABLES ---------------------------------------------------

limited_data <- data_au %>%
  transmute(
    case_id,
    main_analytic,
    issue_check_pass,
    attention_pass,
    major_party_identifier,
    intent,
    vignette,
    candidate_party,
    party_bloc,
    party_id_detail,
    ingroup3,
    state,
    sex,
    sex_num,
    age,
    edu,
    hincome,
    ideology,
    ft_rul,
    ft_opp,
    CM_1,
    CM_2,
    CM_3,
    CM_4,
    CM_5,
    CM_1_bin,
    CM_2_bin,
    CM_3_bin,
    CM_4_bin,
    CM_5_bin,
    conspiracy,
    conspiracy_z,
    conspiracy_full,
    conspiracy_additive,
    negative,
    enthusiasm,
    competent,
    sincere,
    trustworthy,
    undemocratic_perception
  ) %>%
  clean_names()

# SAVE ---------------------------------------------------------------------

write_csv(limited_data, csv_out, na = "")
saveRDS(limited_data, rds_out)

message("Wrote limited data to:")
message(csv_out)
message(rds_out)
message(paste("Rows in limited submitted-sample dataset:", nrow(limited_data)))
message(paste("Rows in main analytic sample:", sum(limited_data$main_analytic == 1, na.rm = TRUE)))
