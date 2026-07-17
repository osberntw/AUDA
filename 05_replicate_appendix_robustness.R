rm(list = ls())
options(scipen = 999)

library(tidyverse)
library(broom)
library(MASS)
library(psych)
library(flextable)
library(emmeans)

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
    p_value < 0.01 ~ "***",
    p_value < 0.05 ~ "**",
    p_value < 0.1 ~ "*",
    TRUE ~ ""
  )
}

format_entry <- function(estimate, std_error, p_value) {
  ifelse(
    is.na(estimate),
    "",
    sprintf("%.3f%s (%.3f)", estimate, star_code(p_value), std_error)
  )
}

format_pvalue <- function(p_value) {
  ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
}

tidy_polr <- function(model) {
  tidy(model) %>%
    mutate(
      statistic = estimate / std.error,
      p.value = 2 * pnorm(abs(statistic), lower.tail = FALSE)
    )
}

build_model_table <- function(models, column_names, coef_map, first_col = "Term", gof_rows = NULL) {
  term_order <- names(coef_map)

  coef_tables <- map2(models, column_names, function(model, column_name) {
    tidy_fn <- if (inherits(model, "polr")) tidy_polr else broom::tidy

    tidy_fn(model) %>%
      filter(term %in% term_order) %>%
      transmute(
        term,
        !!column_name := format_entry(estimate, std.error, p.value)
      )
  })

  merged <- tibble(term = term_order) %>%
    reduce(
      coef_tables,
      .init = .,
      .f = ~ left_join(.x, .y, by = "term")
    ) %>%
    mutate(!!first_col := unname(coef_map[term])) %>%
    dplyr::select(all_of(first_col), all_of(column_names))

  if (!is.null(gof_rows)) {
    gof_table <- map_dfr(names(gof_rows), function(row_name) {
      values <- as.list(gof_rows[[row_name]])
      tibble(!!first_col := row_name, !!!values)
    }) %>%
      dplyr::select(all_of(first_col), all_of(column_names))

    merged <- bind_rows(merged, gof_table)
  }

  merged
}

save_table_outputs <- function(table_df, table_title, csv_name, docx_name) {
  write_csv(table_df, file.path(output_dir, csv_name), na = "")

  ft <- regulartable(table_df) %>%
    autofit() %>%
    theme_booktabs() %>%
    align(align = "left", j = 1, part = "all") %>%
    align(align = "center", j = 2:ncol(table_df), part = "all")

  ft_list <- setNames(list(ft), table_title)
  do.call(save_as_docx, c(ft_list, list(path = file.path(output_dir, docx_name))))
}

build_conditional_effect <- function(model, moderator_values, moderator_name, focal_term) {
  vc <- vcov(model)
  beta_focal <- coef(model)[[focal_term]]
  interaction_term <- paste0(focal_term, ":", moderator_name)
  beta_interaction <- coef(model)[[interaction_term]]

  if (is.na(beta_interaction)) {
    interaction_term <- paste0(moderator_name, ":", focal_term)
    beta_interaction <- coef(model)[[interaction_term]]
  }

  var_focal <- vc[focal_term, focal_term]
  var_interaction <- vc[interaction_term, interaction_term]
  covar <- vc[focal_term, interaction_term]

  tibble(conspiracy = moderator_values) %>%
    mutate(
      estimate = beta_focal + beta_interaction * conspiracy,
      std_error = sqrt(var_focal + (conspiracy ^ 2) * var_interaction + 2 * conspiracy * covar),
      lower = estimate - 1.96 * std_error,
      upper = estimate + 1.96 * std_error
    )
}

# READ DATA ----------------------------------------------------------------

data_au <- read_csv(data_path, show_col_types = FALSE) %>%
  mutate(
    vignette = factor(vignette, levels = c("Control", "Election", "Protests", "Court")),
    ingroup3 = factor(ingroup3, levels = c("In-group", "Out-group", "Independent")),
    sex = factor(sex, levels = c("female", "male")),
    party_id_detail = factor(
      party_id_detail,
      levels = c("Coalition", "Labor", "Greens", "One Nation", "Independent", "No party ID")
    )
  )

data_main <- data_au %>%
  filter(main_analytic == 1) %>%
  mutate(intent_ord = ordered(intent, levels = 0:4))

# APPENDIX 3: SAMPLE FLOW --------------------------------------------------

screening_counts <- tibble(
  Sample = c(
    "Submitted cases with experimental fields",
    "Passed post-vignette issue-identification check"
  ),
  N = c(nrow(data_au), sum(data_au$main_analytic == 1, na.rm = TRUE)),
  `Share of submitted cases (%)` = round(100 * N / nrow(data_au), 1),
  Use = c(
    "Broadest submitted sample",
    "Main analytic sample used in the manuscript"
  )
)

save_table_outputs(
  screening_counts,
  "Screening counts and analytic samples",
  "table_a3c_screening_counts.csv",
  "table_a3c_screening_counts.docx"
)

profile_data <- data_au %>%
  mutate(main_sample = ifelse(main_analytic == 1, "Main analytic sample", "Excluded from main sample"))

profile_specs <- tribble(
  ~variable, ~label, ~type,
  "sex_num", "Male (%)", "binary",
  "age", "Age", "continuous",
  "edu", "Education", "continuous",
  "hincome", "Income", "continuous",
  "ideology", "Ideology", "continuous",
  "major_party_identifier", "Major-party identifier (%)", "binary",
  "attention_pass", "Passed attention check (%)", "binary"
)

profile_rows <- map_dfr(seq_len(nrow(profile_specs)), function(i) {
  variable_name <- profile_specs$variable[i]
  variable_label <- profile_specs$label[i]
  variable_type <- profile_specs$type[i]
  x <- profile_data[[variable_name]]
  g <- profile_data$main_sample

  included_mean <- mean(x[g == "Main analytic sample"], na.rm = TRUE)
  excluded_mean <- mean(x[g == "Excluded from main sample"], na.rm = TRUE)

  if (variable_type == "binary") {
    p_value <- suppressWarnings(chisq.test(table(g, x))$p.value)
    included_value <- sprintf("%.1f", 100 * included_mean)
    excluded_value <- sprintf("%.1f", 100 * excluded_mean)
  } else {
    p_value <- t.test(x ~ g)$p.value
    included_value <- sprintf("%.2f", included_mean)
    excluded_value <- sprintf("%.2f", excluded_mean)
  }

  tibble(
    Variable = variable_label,
    `Main analytic sample (N = 1227)` = included_value,
    `Excluded from main sample (N = 786)` = excluded_value,
    `p-value` = format_pvalue(p_value)
  )
})

save_table_outputs(
  profile_rows,
  "Profile of respondents excluded from the main analytic sample",
  "table_a3d_excluded_profile.csv",
  "table_a3d_excluded_profile.docx"
)

# APPENDIX 3: BROADER-SAMPLE MODELS ---------------------------------------

full_h1 <- lm(intent ~ vignette, data = data_au)
full_h2 <- lm(intent ~ vignette * ingroup3, data = data_au)
full_h3 <- lm(intent ~ vignette * conspiracy_full, data = data_au)

coef_map_full <- c(
  "vignetteElection" = "Election",
  "vignetteProtests" = "Protests",
  "vignetteCourt" = "Court",
  "ingroup3Out-group" = "Out-group",
  "ingroup3Independent" = "Independent",
  "vignetteElection:ingroup3Out-group" = "Election*Out-group",
  "vignetteProtests:ingroup3Out-group" = "Protests* Out-group",
  "vignetteCourt:ingroup3Out-group" = "Court*Out-group",
  "vignetteElection:ingroup3Independent" = "Election*Independent",
  "vignetteProtests:ingroup3Independent" = "Protests*Independent",
  "vignetteCourt:ingroup3Independent" = "Court*Independent",
  "conspiracy_full" = "Conspiracy mentality",
  "vignetteElection:conspiracy_full" = "Election*Conspiracy mentality",
  "vignetteProtests:conspiracy_full" = "Protests*Conspiracy mentality",
  "vignetteCourt:conspiracy_full" = "Court*Conspiracy mentality"
)

table_a4 <- build_model_table(
  models = list(full_h1, full_h2, full_h3),
  column_names = c("H1: Baseline", "H2: Partisan alignment", "H3: Conspiracy moderation"),
  coef_map = coef_map_full,
  first_col = "Variable",
  gof_rows = list(
    "N" = c(
      "H1: Baseline" = as.character(nobs(full_h1)),
      "H2: Partisan alignment" = as.character(nobs(full_h2)),
      "H3: Conspiracy moderation" = as.character(nobs(full_h3))
    )
  )
)

save_table_outputs(
  table_a4,
  "Table A4. Broader submitted-sample models",
  "table_a4_broader_submitted_sample.csv",
  "table_a4_broader_submitted_sample.docx"
)

conspiracy_mean <- mean(data_au$conspiracy_full, na.rm = TRUE)
conspiracy_sd <- sd(data_au$conspiracy_full, na.rm = TRUE)
conspiracy_low <- conspiracy_mean - conspiracy_sd
conspiracy_high <- conspiracy_mean + conspiracy_sd

figure_a3 <- emmeans(
  full_h3,
  ~ vignette | conspiracy_full,
  at = list(conspiracy_full = c(conspiracy_low, conspiracy_mean, conspiracy_high))
) %>%
  contrast(method = "trt.vs.ctrl", ref = 1) %>%
  confint() %>%
  as.data.frame() %>%
  mutate(
    conspiracy_group = case_when(
      abs(conspiracy_full - conspiracy_low) < 1e-8 ~ "Low conspiracy (Mean - 1 SD)",
      abs(conspiracy_full - conspiracy_mean) < 1e-8 ~ "Mean conspiracy",
      abs(conspiracy_full - conspiracy_high) < 1e-8 ~ "High conspiracy (Mean + 1 SD)",
      TRUE ~ NA_character_
    ),
    treatment = recode(
      contrast,
      "Election - Control" = "Election denial",
      "Protests - Control" = "Protest restriction",
      "Court - Control" = "Judicial attack"
    ),
    treatment = factor(treatment, levels = c("Judicial attack", "Protest restriction", "Election denial")),
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

write_csv(figure_a3, file.path(output_dir, "figure_a3_broader_sample_penalties.csv"), na = "")

plot_a3 <- ggplot(figure_a3, aes(x = estimate, y = plot_y, colour = conspiracy_group)) +
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
    breaks = seq_along(levels(figure_a3$treatment)),
    labels = levels(figure_a3$treatment),
    limits = c(0.6, 3.6),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = "Treatment effect on vote intention", y = NULL, colour = NULL) +
  theme_minimal(base_family = "Helvetica", base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(8, 18, 8, 8)
  )

ggsave(
  file.path(output_dir, "figure_a3_broader_sample_penalties.png"),
  plot_a3,
  width = 8,
  height = 4.8,
  dpi = 600
)

conspiracy_range <- seq(
  min(data_au$conspiracy_full, na.rm = TRUE),
  max(data_au$conspiracy_full, na.rm = TRUE),
  length.out = 300
)

figure_a4 <- build_conditional_effect(
  model = full_h3,
  moderator_values = conspiracy_range,
  moderator_name = "conspiracy_full",
  focal_term = "vignetteElection"
) %>%
  mutate(significance = ifelse(lower > 0 | upper < 0, "Significant", "Not significant"))

write_csv(figure_a4, file.path(output_dir, "figure_a4_election_denial_conditional_effect_broader_sample.csv"), na = "")

plot_a4 <- ggplot(figure_a4, aes(x = conspiracy, y = estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#BFD7DC", alpha = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed", colour = "grey45") +
  geom_vline(
    xintercept = c(conspiracy_low, conspiracy_mean, conspiracy_high),
    linewidth = 0.35,
    linetype = "dotted",
    colour = "grey55"
  ) +
  geom_line(linewidth = 0.8, colour = "#007A87") +
  labs(
    x = "Conspiracy mentality",
    y = "Conditional effect of election denial"
  ) +
  theme_minimal(base_family = "Helvetica", base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 18, 8, 8)
  )

ggsave(
  file.path(output_dir, "figure_a4_election_denial_conditional_effect_broader_sample.png"),
  plot_a4,
  width = 8,
  height = 4.8,
  dpi = 600
)

# TABLE A5: ORDERED LOGIT --------------------------------------------------

ologit_h1 <- polr(intent_ord ~ vignette, data = data_main, Hess = TRUE)
ologit_h2 <- polr(intent_ord ~ vignette * ingroup3, data = data_main, Hess = TRUE)
ologit_h3 <- polr(intent_ord ~ vignette * conspiracy, data = data_main, Hess = TRUE)

coef_map_core <- c(
  "vignetteElection" = "Election",
  "vignetteProtests" = "Protests",
  "vignetteCourt" = "Court",
  "ingroup3Out-group" = "Out-group",
  "ingroup3Independent" = "Independent",
  "vignetteElection:ingroup3Out-group" = "Election*Out-group",
  "vignetteProtests:ingroup3Out-group" = "Protests*Out-group",
  "vignetteCourt:ingroup3Out-group" = "Court*Out-group",
  "vignetteElection:ingroup3Independent" = "Election*Independent",
  "vignetteProtests:ingroup3Independent" = "Protests*Independent",
  "vignetteCourt:ingroup3Independent" = "Court*Independent",
  "conspiracy" = "Conspiracy mentality",
  "vignetteElection:conspiracy" = "Election*Conspiracy mentality",
  "vignetteProtests:conspiracy" = "Protests*Conspiracy mentality",
  "vignetteCourt:conspiracy" = "Court*Conspiracy mentality"
)

table_a5 <- build_model_table(
  models = list(ologit_h1, ologit_h2, ologit_h3),
  column_names = c("H1: Baseline", "H2: Partisan alignment", "H3: Conspiracy moderation"),
  coef_map = coef_map_core,
  first_col = "Term",
  gof_rows = list(
    "N" = c(
      "H1: Baseline" = as.character(nobs(ologit_h1)),
      "H2: Partisan alignment" = as.character(nobs(ologit_h2)),
      "H3: Conspiracy moderation" = as.character(nobs(ologit_h3))
    )
  )
)

save_table_outputs(
  table_a5,
  "Table A5. Ordered logit models",
  "table_a5_ordered_logit.csv",
  "table_a5_ordered_logit.docx"
)

# TABLE A6: COVARIATE-ADJUSTED MODELS --------------------------------------

adj_h1 <- lm(intent ~ vignette + sex + edu + hincome + ideology, data = data_main)
adj_h2 <- lm(intent ~ vignette * ingroup3 + sex + edu + hincome + ideology, data = data_main)
adj_h3 <- lm(intent ~ vignette * conspiracy + sex + edu + hincome + ideology, data = data_main)

coef_map_adjusted <- c(
  coef_map_core,
  "sexmale" = "Male",
  "edu" = "Education",
  "hincome" = "Household income",
  "ideology" = "Ideology"
)

table_a6 <- build_model_table(
  models = list(adj_h1, adj_h2, adj_h3),
  column_names = c("H1: Baseline", "H2: Partisan alignment", "H3: Conspiracy moderation"),
  coef_map = coef_map_adjusted,
  first_col = "Term",
  gof_rows = list(
    "N" = c(
      "H1: Baseline" = as.character(nobs(adj_h1)),
      "H2: Partisan alignment" = as.character(nobs(adj_h2)),
      "H3: Conspiracy moderation" = as.character(nobs(adj_h3))
    )
  )
)

save_table_outputs(
  table_a6,
  "Table A6. Covariate-adjusted models",
  "table_a6_covariate_adjustment.csv",
  "table_a6_covariate_adjustment.docx"
)

# TABLE A9: FACTOR LOADINGS ------------------------------------------------

cm_items <- c("cm_1", "cm_2", "cm_3", "cm_4", "cm_5")

fa_data <- data_main %>%
  dplyr::select(all_of(cm_items)) %>%
  drop_na()

fa_result <- psych::fa(fa_data, nfactors = 1, fm = "ml")
loadings <- as.numeric(fa_result$loadings[, 1])

if (mean(loadings, na.rm = TRUE) < 0) {
  loadings <- -loadings
}

table_a9 <- tibble(
  Item = c(
    "Many very important things happen in the world, which the public is never informed about",
    "Politicians usually do not tell us the true motives for their decisions",
    "Government agencies closely monitor all citizens",
    "Events that seem unconnected are often the result of secret activities",
    "Secret organisations greatly influence political decisions"
  ),
  Loading = sprintf("%.3f", loadings),
  Uniqueness = sprintf("%.3f", 1 - fa_result$communality)
) %>%
  bind_rows(
    tibble(
      Item = "N",
      Loading = as.character(nrow(fa_data)),
      Uniqueness = ""
    )
  )

save_table_outputs(
  table_a9,
  "Table A9. One-factor loadings for the conspiracy mentality scale",
  "table_a9_factor_loadings.csv",
  "table_a9_factor_loadings.docx"
)

# TABLE A10: ADDITIVE INDEX ROBUSTNESS -------------------------------------

additive_model <- lm(intent ~ vignette * conspiracy_additive, data = data_main)

coef_map_additive <- c(
  "vignetteElection" = "Election",
  "vignetteProtests" = "Protests",
  "vignetteCourt" = "Court",
  "conspiracy_additive" = "Additive conspiracy index",
  "vignetteElection:conspiracy_additive" = "Election*Additive conspiracy index",
  "vignetteProtests:conspiracy_additive" = "Protests*Additive conspiracy index",
  "vignetteCourt:conspiracy_additive" = "Court*Additive conspiracy index"
)

table_a10 <- build_model_table(
  models = list(additive_model),
  column_names = c("Model 1"),
  coef_map = coef_map_additive,
  first_col = "Term",
  gof_rows = list(
    "N" = c("Model 1" = as.character(nobs(additive_model)))
  )
)

save_table_outputs(
  table_a10,
  "Table A10. Additive conspiracy-index robustness check",
  "table_a10_conspiracy_additive.csv",
  "table_a10_conspiracy_additive.docx"
)

message("Appendix robustness outputs written to: ", output_dir)
