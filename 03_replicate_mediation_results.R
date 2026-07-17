rm(list = ls())
options(scipen = 999)

library(tidyverse)
library(broom)
library(mediation)
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

set.seed(12234)
mediation_sims <- 500

# HELPERS ------------------------------------------------------------------

format_pvalue <- function(x) {
  ifelse(x < 0.001, "<.001", sprintf("%.3f", x))
}

format_point <- function(x) {
  ifelse(abs(x) < 0.005, "0.00", sprintf("%.2f", x))
}

format_bound <- function(x) {
  ifelse(x < 0 & abs(x) < 0.005, sprintf("%.3f", x), format_point(x))
}

format_interval <- function(estimate, lower, upper) {
  sprintf("%s [%s, %s]", format_point(estimate), format_bound(lower), format_bound(upper))
}

get_term <- function(tidy_df, term_name, column_name) {
  candidate_terms <- term_name

  if (grepl(":", term_name, fixed = TRUE)) {
    term_parts <- strsplit(term_name, ":", fixed = TRUE)[[1]]
    candidate_terms <- unique(c(term_name, paste(rev(term_parts), collapse = ":")))
  }

  out <- tidy_df %>%
    filter(term %in% candidate_terms) %>%
    slice(1) %>%
    pull(all_of(column_name))

  if (length(out) == 0) {
    return(NA_real_)
  }

  out[[1]]
}

save_docx_table <- function(table_df, title, filename) {
  ft <- regulartable(table_df) %>%
    autofit() %>%
    theme_booktabs() %>%
    align(align = "left", j = 1:2, part = "all") %>%
    align(align = "center", j = 3:ncol(table_df), part = "all")

  ft_list <- setNames(list(ft), title)
  do.call(save_as_docx, c(ft_list, list(path = file.path(output_dir, filename))))
}

# READ DATA ----------------------------------------------------------------

data_au <- read_csv(data_path, show_col_types = FALSE) %>%
  filter(main_analytic == 1) %>%
  mutate(
    vignette = factor(vignette, levels = c("Control", "Election", "Protests", "Court")),
    election = ifelse(vignette == "Election", 1, 0),
    protests = ifelse(vignette == "Protests", 1, 0),
    court = ifelse(vignette == "Court", 1, 0)
  )

# CONDITIONAL MEDIATION ----------------------------------------------------

pretty_mediator <- c(
  undemocratic_perception = "Undemocraticness",
  negative = "Negative emotion",
  enthusiasm = "Enthusiasm",
  trustworthy = "Trustworthiness",
  competent = "Competence",
  sincere = "Sincerity"
)

core_mediators <- names(pretty_mediator)
treatments <- c("Election", "Protests", "Court")
treat_vars <- c("election", "protests", "court")
interaction_vars <- c("election:conspiracy_z", "protests:conspiracy_z", "court:conspiracy_z")

run_pooled_mediation_one <- function(data, mediator_name, sims = 500) {
  analysis_data <- data %>%
    filter(
      complete.cases(
        intent,
        .data[[mediator_name]],
        election,
        protests,
        court,
        conspiracy_z
      )
    )

  mediator_formula <- as.formula(
    paste0(
      mediator_name,
      " ~ election * conspiracy_z + protests * conspiracy_z + court * conspiracy_z"
    )
  )

  outcome_formula <- as.formula(
    paste0(
      "intent ~ ",
      mediator_name,
      " * conspiracy_z + election * conspiracy_z + protests * conspiracy_z + court * conspiracy_z"
    )
  )

  mediator_model <- lm(mediator_formula, data = analysis_data)
  outcome_model <- lm(outcome_formula, data = analysis_data)
  total_model <- lm(
    intent ~ election * conspiracy_z + protests * conspiracy_z + court * conspiracy_z,
    data = analysis_data
  )

  mediator_tidy <- broom::tidy(mediator_model)
  total_tidy <- broom::tidy(total_model)

  result_list <- vector("list", length(treatments))

  for (i in seq_along(treatments)) {
    med_low <- mediate(
      mediator_model,
      outcome_model,
      treat = treat_vars[i],
      mediator = mediator_name,
      covariates = list(conspiracy_z = -1),
      robustSE = TRUE,
      sims = sims
    )

    med_high <- mediate(
      mediator_model,
      outcome_model,
      treat = treat_vars[i],
      mediator = mediator_name,
      covariates = list(conspiracy_z = 1),
      robustSE = TRUE,
      sims = sims
    )

    b_med_treat <- get_term(mediator_tidy, treat_vars[i], "estimate")
    p_med_treat <- get_term(mediator_tidy, treat_vars[i], "p.value")
    b_med_int <- get_term(mediator_tidy, interaction_vars[i], "estimate")
    p_med_int <- get_term(mediator_tidy, interaction_vars[i], "p.value")

    b_total_treat <- get_term(total_tidy, treat_vars[i], "estimate")
    p_total_treat <- get_term(total_tidy, treat_vars[i], "p.value")
    b_total_int <- get_term(total_tidy, interaction_vars[i], "estimate")
    p_total_int <- get_term(total_tidy, interaction_vars[i], "p.value")

    result_list[[i]] <- tibble(
      treatment = treatments[i],
      mediator = mediator_name,
      mediator_label = pretty_mediator[[mediator_name]],
      n = nrow(analysis_data),
      mediator_treat_coef = round(b_med_treat, 3),
      mediator_treat_p = format_pvalue(p_med_treat),
      mediator_treat_consp_coef = round(b_med_int, 3),
      mediator_treat_consp_p = format_pvalue(p_med_int),
      mediator_effect_low = round(b_med_treat - b_med_int, 3),
      mediator_effect_high = round(b_med_treat + b_med_int, 3),
      total_treat_coef = round(b_total_treat, 3),
      total_treat_p = format_pvalue(p_total_treat),
      total_treat_consp_coef = round(b_total_int, 3),
      total_treat_consp_p = format_pvalue(p_total_int),
      total_effect_low = round(b_total_treat - b_total_int, 3),
      total_effect_high = round(b_total_treat + b_total_int, 3),
      acme_low = med_low$d.avg,
      acme_low_lower = med_low$d.avg.ci[1],
      acme_low_upper = med_low$d.avg.ci[2],
      acme_low_p = format_pvalue(med_low$d.avg.p),
      ade_low = round(med_low$z.avg, 3),
      ade_low_p = format_pvalue(med_low$z.avg.p),
      prop_med_low = round(med_low$n.avg, 3),
      prop_med_low_p = format_pvalue(med_low$n.avg.p),
      acme_high = med_high$d.avg,
      acme_high_lower = med_high$d.avg.ci[1],
      acme_high_upper = med_high$d.avg.ci[2],
      acme_high_p = format_pvalue(med_high$d.avg.p),
      ade_high = round(med_high$z.avg, 3),
      ade_high_p = format_pvalue(med_high$z.avg.p),
      prop_med_high = round(med_high$n.avg, 3),
      prop_med_high_p = format_pvalue(med_high$n.avg.p)
    )
  }

  bind_rows(result_list)
}

core_results <- bind_rows(
  lapply(core_mediators, run_pooled_mediation_one, data = data_au, sims = mediation_sims)
)

write_csv(core_results, file.path(output_dir, "mediation_core_results_long.csv"), na = "")

table_2 <- core_results %>%
  filter(treatment %in% c("Election", "Court")) %>%
  mutate(
    treatment = factor(treatment, levels = c("Election", "Court")),
    mediator_label = factor(mediator_label, levels = unname(pretty_mediator)),
    `Low-conspiracy ACME` = format_interval(acme_low, acme_low_lower, acme_low_upper),
    `High-conspiracy ACME` = format_interval(acme_high, acme_high_lower, acme_high_upper)
  ) %>%
  arrange(treatment, mediator_label) %>%
  transmute(
    Treatment = as.character(treatment),
    Mediator = as.character(mediator_label),
    `Low-conspiracy ACME`,
    `High-conspiracy ACME`
  )

write_csv(table_2, file.path(output_dir, "table_2_conditional_mediation_effects.csv"), na = "")
save_docx_table(table_2, "Table 2. Conditional mediation effects for core mechanisms", "table_2_conditional_mediation_effects.docx")

# MAIN FIGURE 2 ------------------------------------------------------------

figure_2 <- bind_rows(
  core_results %>%
    filter(treatment %in% c("Election", "Court")) %>%
    transmute(
      treatment,
      mediator_label,
      conspiracy_group = "Low conspiracy (Mean - 1 SD)",
      estimate = acme_low,
      lower = acme_low_lower,
      upper = acme_low_upper
    ),
  core_results %>%
    filter(treatment %in% c("Election", "Court")) %>%
    transmute(
      treatment,
      mediator_label,
      conspiracy_group = "High conspiracy (Mean + 1 SD)",
      estimate = acme_high,
      lower = acme_high_lower,
      upper = acme_high_upper
    )
) %>%
  mutate(
    treatment = factor(treatment, levels = c("Election", "Court")),
    mediator_label = factor(mediator_label, levels = rev(unname(pretty_mediator))),
    conspiracy_group = factor(
      conspiracy_group,
      levels = c("Low conspiracy (Mean - 1 SD)", "High conspiracy (Mean + 1 SD)")
    )
  )

write_csv(figure_2, file.path(output_dir, "figure_2_conditional_mediation_effects.csv"), na = "")

plot_2 <- ggplot(
  figure_2,
  aes(x = estimate, y = mediator_label, colour = conspiracy_group)
) +
  geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed", colour = "grey45") +
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    orientation = "y",
    position = position_dodge(width = 0.55),
    width = 0.2,
    linewidth = 0.55
  ) +
  geom_point(position = position_dodge(width = 0.55), size = 2.2) +
  facet_wrap(~ treatment) +
  scale_colour_manual(values = c("#007A87", "#C7432B")) +
  labs(
    x = "Average causal mediation effect",
    y = NULL,
    colour = NULL
  ) +
  theme_minimal(base_family = "Helvetica", base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  file.path(output_dir, "figure_2_conditional_mediation_effects.png"),
  plot_2,
  width = 9,
  height = 5.8,
  dpi = 600
)

message("Mediation replication outputs written to: ", output_dir)
