rm(list = ls())
options(scipen = 999)

library(tidyverse)
library(broom)
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

# HELPERS ------------------------------------------------------------------

star_code <- function(p_value) {
  case_when(
    is.na(p_value) ~ "",
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ ""
  )
}

format_estimate <- function(estimate, p_value) {
  ifelse(is.na(estimate), "", sprintf("%.3f%s", estimate, star_code(p_value)))
}

format_se <- function(std_error) {
  ifelse(is.na(std_error), "", sprintf("(%.3f)", std_error))
}

format_pvalue <- function(p_value) {
  ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
}

save_docx_table <- function(table_df, title, filename) {
  ft <- regulartable(table_df) %>%
    autofit() %>%
    theme_booktabs() %>%
    align(align = "left", j = 1, part = "all") %>%
    align(align = "center", j = 2:ncol(table_df), part = "all")

  ft_list <- setNames(list(ft), title)
  do.call(save_as_docx, c(ft_list, list(path = file.path(output_dir, filename))))
}

build_regression_table <- function(models, column_names, term_map) {
  term_order <- names(term_map)

  model_tables <- map2(models, column_names, function(model, column_name) {
    tidy(model) %>%
      filter(term %in% term_order) %>%
      transmute(
        term,
        estimate_row = format_estimate(estimate, p.value),
        se_row = format_se(std.error)
      ) %>%
      pivot_longer(
        cols = c(estimate_row, se_row),
        names_to = "row_type",
        values_to = column_name
      )
  })

  table_body <- tibble(term = term_order) %>%
    crossing(row_type = c("estimate_row", "se_row")) %>%
    mutate(row_type = factor(row_type, levels = c("estimate_row", "se_row"))) %>%
    arrange(match(term, term_order), row_type) %>%
    reduce(
      model_tables,
      .init = .,
      .f = ~ left_join(.x, .y, by = c("term", "row_type"))
    ) %>%
    mutate(Predictor = ifelse(row_type == "estimate_row", unname(term_map[term]), "")) %>%
    select(Predictor, all_of(column_names))

  n_row <- tibble(Predictor = "N")
  for (i in seq_along(models)) {
    n_row[[column_names[i]]] <- format(nobs(models[[i]]), big.mark = ",")
  }

  bind_rows(table_body, n_row)
}

build_balance_table <- function(data, table_title, csv_name, docx_name) {
  balance_vars <- c("sex_num", "age", "edu", "hincome", "ideology")
  balance_labels <- c(
    sex_num = "Male",
    age = "Age",
    edu = "Educ.",
    hincome = "Income",
    ideology = "Ideology"
  )

  balance_table <- map_dfr(balance_vars, function(variable_name) {
    values <- data[[variable_name]]

    tibble(
      Variable = balance_labels[[variable_name]],
      Control = sprintf(
        "%.2f (%.2f)",
        mean(values[data$vignette == "Control"], na.rm = TRUE),
        sd(values[data$vignette == "Control"], na.rm = TRUE)
      ),
      Election = sprintf(
        "%.2f (%.2f)",
        mean(values[data$vignette == "Election"], na.rm = TRUE),
        sd(values[data$vignette == "Election"], na.rm = TRUE)
      ),
      Protest = sprintf(
        "%.2f (%.2f)",
        mean(values[data$vignette == "Protests"], na.rm = TRUE),
        sd(values[data$vignette == "Protests"], na.rm = TRUE)
      ),
      Court = sprintf(
        "%.2f (%.2f)",
        mean(values[data$vignette == "Court"], na.rm = TRUE),
        sd(values[data$vignette == "Court"], na.rm = TRUE)
      ),
      `p-value` = format_pvalue(
        kruskal.test(reformulate("vignette", response = variable_name), data = data)$p.value
      )
    )
  })

  write_csv(balance_table, file.path(output_dir, csv_name), na = "")
  save_docx_table(balance_table, table_title, docx_name)
}

# READ DATA ----------------------------------------------------------------

data_au <- read_csv(data_path, show_col_types = FALSE) %>%
  mutate(
    vignette = factor(vignette, levels = c("Control", "Election", "Protests", "Court")),
    candidate_party = factor(candidate_party, levels = c("ruling", "opposition")),
    party_bloc = factor(party_bloc, levels = c("ruling", "opposition", "Independent")),
    party_id_detail = factor(
      party_id_detail,
      levels = c("Coalition", "Labor", "Greens", "One Nation", "Independent", "No party ID")
    ),
    ingroup3 = factor(ingroup3, levels = c("In-group", "Out-group", "Independent")),
    sex = factor(sex, levels = c("female", "male"))
  )

data_main <- data_au %>%
  filter(main_analytic == 1)

# APPENDIX 1: CMQ DISTRIBUTIONS -------------------------------------------

cm_labels <- c(
  cm_1_bin = "Many very important things happen\nin the world, which the public\nis never informed about",
  cm_2_bin = "Politicians usually do not tell us\nthe true motives for their decisions",
  cm_3_bin = "Government agencies closely\nmonitor all citizens",
  cm_4_bin = "Events which superficially seem\nto lack a connection are often\nthe result of secret activities",
  cm_5_bin = "There are secret organisations\nthat greatly influence\npolitical decisions"
)

figure_a1 <- data_au %>%
  select(all_of(names(cm_labels))) %>%
  pivot_longer(cols = everything(), names_to = "item", values_to = "agree_top2") %>%
  group_by(item) %>%
  summarise(percent = 100 * mean(agree_top2, na.rm = TRUE), .groups = "drop") %>%
  mutate(label = factor(cm_labels[item], levels = rev(cm_labels)))

write_csv(figure_a1, file.path(output_dir, "figure_a1_cmq_raw_agreement_broader_sample.csv"))

plot_a1 <- ggplot(figure_a1, aes(x = percent, y = label)) +
  geom_col(fill = "grey40") +
  geom_text(
    aes(label = sprintf("%.2f%%", percent), x = percent + 1),
    hjust = 0,
    size = 3.8
  ) +
  scale_x_continuous(limits = c(0, max(figure_a1$percent) + 6)) +
  labs(
    x = "Percentage somewhat or strongly agreeing",
    y = NULL
  ) +
  theme_classic()

ggsave(
  file.path(output_dir, "figure_a1_cmq_raw_agreement_broader_sample.png"),
  plot_a1,
  width = 8,
  height = 6,
  dpi = 600
)

figure_a2_summary <- data_au %>%
  summarise(
    n = sum(!is.na(conspiracy_full)),
    mean = mean(conspiracy_full, na.rm = TRUE),
    sd = sd(conspiracy_full, na.rm = TRUE),
    min = min(conspiracy_full, na.rm = TRUE),
    max = max(conspiracy_full, na.rm = TRUE)
  )

write_csv(figure_a2_summary, file.path(output_dir, "figure_a2_cmq_factor_score_summary_broader_sample.csv"))

plot_a2 <- ggplot(data_au, aes(x = conspiracy_full)) +
  geom_histogram(bins = 35, fill = "grey40", colour = "white", linewidth = 0.2) +
  geom_vline(
    xintercept = figure_a2_summary$mean,
    linewidth = 0.6,
    colour = "grey15"
  ) +
  geom_vline(
    xintercept = c(
      figure_a2_summary$mean - figure_a2_summary$sd,
      figure_a2_summary$mean + figure_a2_summary$sd
    ),
    linewidth = 0.5,
    linetype = "dashed",
    colour = "grey35"
  ) +
  labs(
    x = "Conspiracy-mentality factor score",
    y = "Number of respondents"
  ) +
  theme_classic()

ggsave(
  file.path(output_dir, "figure_a2_cmq_factor_score_distribution_broader_sample.png"),
  plot_a2,
  width = 8,
  height = 5,
  dpi = 600
)

# APPENDIX 1: DESCRIPTIVE TABLES ------------------------------------------

desc_vars <- c("intent", "sex_num", "age", "edu", "hincome", "ideology", "conspiracy")
desc_labels <- c(
  intent = "Vote intent.",
  sex_num = "Male",
  age = "Age",
  edu = "Educ.",
  hincome = "Income",
  ideology = "Conservative Ideology (0-1)",
  conspiracy = "Conspiracy score"
)

desc_table <- data_main %>%
  select(all_of(desc_vars)) %>%
  summarise(across(
    everything(),
    list(
      n = ~ sum(!is.na(.)),
      mean = ~ mean(., na.rm = TRUE),
      sd = ~ sd(., na.rm = TRUE),
      min = ~ min(., na.rm = TRUE),
      max = ~ max(., na.rm = TRUE)
    ),
    .names = "{.col}_{.fn}"
  )) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("variable", ".value"),
    names_pattern = "^(.*)_(n|mean|sd|min|max)$"
  ) %>%
  mutate(
    variable = desc_labels[variable],
    across(c(mean, sd, min, max), ~ round(., 3))
  ) %>%
  rename(Variable = variable, N = n, Mean = mean, SD = sd, Min = min, Max = max)

write_csv(desc_table, file.path(output_dir, "table_a1_descriptive_statistics.csv"), na = "")
save_docx_table(desc_table, "Table A1. Descriptive statistics for the analytic sample", "table_a1_descriptive_statistics.docx")

categorical_table <- bind_rows(
  data_main %>%
    count(Variable = "State", Category = state, name = "N"),
  data_main %>%
    count(Variable = "Detailed party ID", Category = as.character(party_id_detail), name = "N"),
  data_main %>%
    count(Variable = "Partisan alignment", Category = as.character(ingroup3), name = "N"),
  data_main %>%
    mutate(
      Category = recode(
        as.character(vignette),
        Control = "Control",
        Election = "Election denial",
        Protests = "Protest restriction",
        Court = "Judicial attack"
      )
    ) %>%
    count(Variable = "Experimental vignette", Category, name = "N"),
  data_main %>%
    mutate(
      Category = recode(
        as.character(candidate_party),
        opposition = "Opposition",
        ruling = "Ruling"
      )
    ) %>%
    count(Variable = "MP party in vignette", Category, name = "N")
) %>%
  group_by(Variable) %>%
  mutate(Percent = round(100 * N / sum(N), 2)) %>%
  ungroup()

write_csv(categorical_table, file.path(output_dir, "table_a2_categorical_descriptives.csv"), na = "")
save_docx_table(categorical_table, "Table A2. Descriptive statistics for categorical variables", "table_a2_categorical_descriptives.docx")

# APPENDIX 2: BALANCE TABLES ----------------------------------------------

build_balance_table(
  data = data_au,
  table_title = "Table A3a. Balance table by vignette condition, broader sample",
  csv_name = "table_a3a_balance_broader_sample.csv",
  docx_name = "table_a3a_balance_broader_sample.docx"
)

build_balance_table(
  data = data_main,
  table_title = "Table A3b. Balance table by vignette condition, analytic sample",
  csv_name = "table_a3b_balance_analytic_sample.csv",
  docx_name = "table_a3b_balance_analytic_sample.docx"
)

# MAIN MODELS --------------------------------------------------------------

model_h1 <- lm(intent ~ vignette, data = data_main)
model_h2 <- lm(intent ~ vignette * ingroup3, data = data_main)
model_h3 <- lm(intent ~ vignette * conspiracy, data = data_main)

term_map <- c(
  "vignetteElection" = "Election denial",
  "vignetteProtests" = "Protest restriction",
  "vignetteCourt" = "Judicial attack",
  "ingroup3Out-group" = "Outgroup",
  "ingroup3Independent" = "Independent/non-aligned",
  "vignetteElection:ingroup3Out-group" = "Election denial*outgroup",
  "vignetteProtests:ingroup3Out-group" = "Protest restriction*outgroup",
  "vignetteCourt:ingroup3Out-group" = "Judicial attack*outgroup",
  "vignetteElection:ingroup3Independent" = "Election denial*non-aligned",
  "vignetteProtests:ingroup3Independent" = "Protest restriction*non-aligned",
  "vignetteCourt:ingroup3Independent" = "Judicial attack*non-aligned",
  "conspiracy" = "Conspiracy mentality",
  "vignetteElection:conspiracy" = "Election denial*conspiracy mentality",
  "vignetteProtests:conspiracy" = "Protest restriction*conspiracy mentality",
  "vignetteCourt:conspiracy" = "Judicial attack*conspiracy mentality",
  "(Intercept)" = "Constant"
)

table_1 <- build_regression_table(
  models = list(model_h1, model_h2, model_h3),
  column_names = c("H1: Baseline", "H2: Partisan alignment", "H3: Conspiracy moderation"),
  term_map = term_map
)

write_csv(table_1, file.path(output_dir, "table_1_vote_intention_models.csv"), na = "")
save_docx_table(table_1, "Table 1. Vote intention models", "table_1_vote_intention_models.docx")

# MAIN FIGURE 1 ------------------------------------------------------------

conspiracy_mean <- mean(data_main$conspiracy, na.rm = TRUE)
conspiracy_sd <- sd(data_main$conspiracy, na.rm = TRUE)
conspiracy_low <- conspiracy_mean - conspiracy_sd
conspiracy_high <- conspiracy_mean + conspiracy_sd

figure_1 <- emmeans(
  model_h3,
  ~ vignette | conspiracy,
  at = list(conspiracy = c(conspiracy_low, conspiracy_mean, conspiracy_high))
) %>%
  contrast(method = "trt.vs.ctrl", ref = 1) %>%
  confint() %>%
  as.data.frame() %>%
  mutate(
    conspiracy_group = case_when(
      abs(conspiracy - conspiracy_low) < 1e-8 ~ "Low conspiracy (Mean - 1 SD)",
      abs(conspiracy - conspiracy_mean) < 1e-8 ~ "Mean conspiracy",
      abs(conspiracy - conspiracy_high) < 1e-8 ~ "High conspiracy (Mean + 1 SD)",
      TRUE ~ NA_character_
    ),
    treatment = recode(
      contrast,
      "Election - Control" = "Election denial",
      "Protests - Control" = "Protest restriction",
      "Court - Control" = "Judicial attack"
    ),
    treatment = factor(
      treatment,
      levels = c("Judicial attack", "Protest restriction", "Election denial")
    ),
    conspiracy_group = factor(
      conspiracy_group,
      levels = c("Low conspiracy (Mean - 1 SD)", "Mean conspiracy", "High conspiracy (Mean + 1 SD)")
    ),
    coefficient_label = sprintf("%4.2f", estimate),
    treatment_y = as.numeric(treatment),
    group_offset = case_when(
      conspiracy_group == "Low conspiracy (Mean - 1 SD)" ~ -0.23,
      conspiracy_group == "Mean conspiracy" ~ 0,
      conspiracy_group == "High conspiracy (Mean + 1 SD)" ~ 0.23,
      TRUE ~ 0
    ),
    plot_y = treatment_y + group_offset,
    label_y = plot_y + 0.11
  )

write_csv(figure_1, file.path(output_dir, "figure_1_conditional_treatment_penalties.csv"), na = "")

plot_1 <- ggplot(figure_1, aes(x = estimate, y = plot_y, colour = conspiracy_group)) +
  geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed", colour = "grey45") +
  geom_errorbar(
    aes(xmin = lower.CL, xmax = upper.CL),
    orientation = "y",
    width = 0.12,
    linewidth = 0.55
  ) +
  geom_point(size = 2.4) +
  geom_label(
    aes(y = label_y, label = coefficient_label),
    size = 3.1,
    fill = "white",
    linewidth = 0,
    label.padding = grid::unit(0.08, "lines"),
    show.legend = FALSE
  ) +
  scale_colour_manual(values = c("#007A87", "#666666", "#C7432B")) +
  scale_y_continuous(
    breaks = seq_along(levels(figure_1$treatment)),
    labels = levels(figure_1$treatment),
    limits = c(0.6, 3.6),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Treatment effect on vote intention",
    y = NULL,
    colour = NULL
  ) +
  theme_minimal(base_family = "Helvetica", base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(8, 18, 8, 8)
  )

ggsave(
  file.path(output_dir, "figure_1_conditional_treatment_penalties.png"),
  plot_1,
  width = 8,
  height = 4.8,
  dpi = 600
)

# SAVE MODEL LOG -----------------------------------------------------------

sink(file.path(output_dir, "model_output.txt"))
cat("Broader submitted sample:", nrow(data_au), "\n")
cat("Main analytic sample:", nrow(data_main), "\n\n")
cat("Table 1 model: H1\n")
print(summary(model_h1))
cat("\nTable 1 model: H2\n")
print(summary(model_h2))
cat("\nTable 1 model: H3\n")
print(summary(model_h3))
sink()

message("Main-text and descriptive appendix outputs written to: ", output_dir)
