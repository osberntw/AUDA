rm(list = ls())
options(scipen = 999)

library(tidyverse)
library(emmeans)
library(flextable)

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

data_path <- file.path(script_dir, "au_replication_limited_data.csv")
output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# READ DATA ----------------------------------------------------------------

data_au <- read_csv(data_path, show_col_types = FALSE) %>%
  filter(main_analytic == 1)

data_au <- data_au %>%
  mutate(
    vignette = factor(vignette, levels = c("Control", "Election", "Protests", "Court")),
    candidate_party = factor(candidate_party, levels = c("ruling", "opposition")),
    party_id_detail = factor(
      party_id_detail,
      levels = c("Coalition", "Labor", "Greens", "One Nation", "Independent", "No party ID")
    )
  )

# DERIVE VARIABLES ---------------------------------------------------------

data_au <- data_au %>%
  mutate(
    ingroup_ft = case_when(
      ft_opp > ft_rul & candidate_party == "opposition" ~ "In-group",
      ft_opp > ft_rul & candidate_party == "ruling" ~ "Out-group",
      ft_opp < ft_rul & candidate_party == "ruling" ~ "In-group",
      ft_opp < ft_rul & candidate_party == "opposition" ~ "Out-group",
      TRUE ~ NA_character_
    ),
    ingroup_ft = factor(ingroup_ft, levels = c("In-group", "Out-group"))
  )

# MODELS -------------------------------------------------------------------

model_party_detail <- lm(intent ~ vignette * party_id_detail, data = data_au)
model_party_ft <- lm(intent ~ vignette * ingroup_ft, data = data_au)

# OUTPUT TABLES ------------------------------------------------------------

party_counts <- data_au %>%
  count(party_id_detail, name = "N")

party_detail_long <- contrast(
  emmeans(model_party_detail, ~ vignette | party_id_detail),
  method = "trt.vs.ctrl",
  ref = 1
) %>%
  summary(infer = TRUE) %>%
  as_tibble() %>%
  left_join(
    party_counts %>% mutate(party_id_detail = as.character(party_id_detail)),
    by = "party_id_detail"
  ) %>%
  mutate(
    treatment = recode(
      contrast,
      "Election - Control" = "Election",
      "Protests - Control" = "Protests",
      "Court - Control" = "Court"
    ),
    stars = case_when(
      p.value < 0.01 ~ "***",
      p.value < 0.05 ~ "**",
      p.value < 0.1 ~ "*",
      TRUE ~ ""
    ),
    estimate_se = sprintf("%.2f%s (%.2f)", estimate, stars, SE)
  ) %>%
  select(
    party_id_detail,
    N,
    treatment,
    estimate,
    SE,
    lower.CL,
    upper.CL,
    p.value,
    estimate_se
  )

party_detail_table <- party_detail_long %>%
  select(party_id_detail, N, treatment, estimate_se) %>%
  pivot_wider(names_from = treatment, values_from = estimate_se) %>%
  rename(
    `Party identification` = party_id_detail,
    `Election - Control` = Election,
    `Protests - Control` = Protests,
    `Court - Control` = Court
  )

write_csv(
  party_detail_long,
  file.path(output_dir, "table_a7_party_alignment_detail_long.csv")
)

write_csv(
  party_detail_table,
  file.path(output_dir, "table_a7_party_alignment_detail.csv")
)

ft_counts <- data_au %>%
  filter(!is.na(ingroup_ft)) %>%
  count(ingroup_ft, name = "N")

party_ft_long <- contrast(
  emmeans(model_party_ft, ~ vignette | ingroup_ft),
  method = "trt.vs.ctrl",
  ref = 1
) %>%
  summary(infer = TRUE) %>%
  as_tibble() %>%
  mutate(
    treatment = recode(
      contrast,
      "Election - Control" = "Election",
      "Protests - Control" = "Protests",
      "Court - Control" = "Court"
    ),
    stars = case_when(
      p.value < 0.01 ~ "***",
      p.value < 0.05 ~ "**",
      p.value < 0.1 ~ "*",
      TRUE ~ ""
    ),
    estimate_se = sprintf("%.2f%s (%.2f)", estimate, stars, SE)
  ) %>%
  left_join(ft_counts %>% mutate(ingroup_ft = as.character(ingroup_ft)), by = "ingroup_ft") %>%
  select(
    ingroup_ft,
    N,
    treatment,
    estimate,
    SE,
    lower.CL,
    upper.CL,
    p.value,
    estimate_se
  )

party_ft_table <- party_ft_long %>%
  select(ingroup_ft, N, treatment, estimate_se) %>%
  pivot_wider(names_from = treatment, values_from = estimate_se) %>%
  rename(
    `Feeling-thermometer alignment` = ingroup_ft,
    `Election - Control` = Election,
    `Protests - Control` = Protests,
    `Court - Control` = Court
  )

write_csv(
  party_ft_long,
  file.path(output_dir, "table_a8_party_alignment_feeling_thermometer_long.csv")
)

write_csv(
  party_ft_table,
  file.path(output_dir, "table_a8_party_alignment_feeling_thermometer.csv")
)

ft_party_detail <- regulartable(party_detail_table) %>%
  autofit() %>%
  theme_booktabs() %>%
  align(align = "left", j = 1, part = "all") %>%
  align(align = "center", j = 2:5, part = "all")

save_as_docx(
  "Table A7. Treatment contrasts by detailed party identification" = ft_party_detail,
  path = file.path(output_dir, "table_a7_party_alignment_detail.docx")
)

ft_party_ft <- regulartable(party_ft_table) %>%
  autofit() %>%
  theme_booktabs() %>%
  align(align = "left", j = 1, part = "all") %>%
  align(align = "center", j = 2:5, part = "all")

save_as_docx(
  "Table A8. Treatment contrasts by feeling-thermometer alignment" = ft_party_ft,
  path = file.path(output_dir, "table_a8_party_alignment_feeling_thermometer.docx")
)
