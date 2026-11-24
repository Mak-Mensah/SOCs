######################################################################
## LAND-USE META-ANALYSIS + PUBLICATION PIPELINE
## Land-use transitions -> Soil Organic Carbon (SOC) stock change
##
## PRIMARY ANALYSIS:
##   Raw mean difference (MD)
##   Multilevel random-effects REML
##   Shared-control sampling covariance
##   Publication-clustered CR2 robust inference
##
## SHARED-CONTROL MATCHING:
##   1. Explicit control_group_id, when available
##   2. Otherwise: publication + soil depth + control mean + SD + n
##
## IMPORTANT:
##   Summary-statistic matching is an inferred control identity.
##   Inspect the shared-control audit before interpreting results.
######################################################################


## ================================================================
## 0. PACKAGES
## ================================================================

required_packages <- c(
  "tidyverse", "metafor", "clubSandwich", "janitor",
  "readr", "fs", "glue", "stringr", "patchwork",
  "scales", "forcats", "flextable", "officer"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install these packages before running the script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(metafor)
  library(clubSandwich)
  library(janitor)
  library(readr)
  library(fs)
  library(glue)
  library(stringr)
  library(patchwork)
  library(scales)
  library(forcats)
  library(flextable)
  library(officer)
})


## ================================================================
## 1. PATHS AND OUTPUT DIRECTORIES
## ================================================================

input_csv <- "Data_112326.csv"

output_dir <- "meta_outputs"
pub_dir <- file.path(output_dir, "publication_ready")

dir_create(output_dir)
dir_create(pub_dir)

if (!file.exists(input_csv)) {
  stop(
    "Cannot find ", input_csv,
    ". Put the CSV in your working directory."
  )
}


## ================================================================
## 2. GENERAL HELPERS
## ================================================================

safe_file <- function(x) {
  x |>
    as.character() |>
    str_replace_all("[/\\\\]", "-") |>
    str_replace_all("[:*?\"<>|]", "") |>
    str_replace_all("\\s+", "_") |>
    str_replace_all("_+", "_") |>
    str_trim()
}

theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = rel(0.9)),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

z_na <- function(x) {
  if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) {
    return(rep(NA_real_, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

fmt_num <- function(x, digits = 2) {
  ifelse(
    is.na(x) | !is.finite(x),
    "",
    formatC(as.numeric(x), format = "f", digits = digits)
  )
}

fmt_p <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "<0.001",
    TRUE ~ fmt_num(p, 3)
  )
}

fmt_ci <- function(est, lo, hi, digits = 2) {
  paste0(
    fmt_num(est, digits), " (",
    fmt_num(lo, digits), " to ",
    fmt_num(hi, digits), ")"
  )
}

safe_range <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(NA_real_, NA_real_))
  range(x)
}

extract_publication_id <- function(x) {
  str_remove(as.character(x), "_[0-9]+$")
}

canon_landuse <- function(x) {
  x0 <- str_to_lower(str_squish(as.character(x)))
  
  case_when(
    str_detect(
      x0,
      "cropland|crop|arable|cultivat|agric|farmland|tillage|row\\s*crop|paddy|rice"
    ) ~ "cropland",
    
    str_detect(
      x0,
      "agroforest|forest|woodland|woody|tree|timber|plantation|afforest|reforest|shrub"
    ) ~ "forest",
    
    str_detect(
      x0,
      "grassland|pasture|rangeland|meadow|steppe|savann|prairie|graz"
    ) ~ "grassland",
    
    str_detect(
      x0,
      "degrad|bare|erod|desert|mine|spoils?|tailings?|wasteland|abandon|salin|gull"
    ) ~ "degraded",
    
    str_detect(
      x0,
      "restor|reclam|rehab|reveg|remediat|improv|recover|regenerat"
    ) ~ "restored",
    
    TRUE ~ "other"
  )
}

parse_depth_interval <- function(x) {
  x <- as.character(x)
  
  x <- str_replace_all(
    x, "[\u2012–\u2013\u2014\u2212]", "-"
  )
  x <- str_replace_all(x, "_", "-")
  
  x <- str_remove_all(
    x,
    regex("\\s*(cm|centimeters?|centimetres?)\\s*$", ignore_case = TRUE)
  )
  
  x <- str_replace_all(
    x,
    regex("\\s*to\\s*", ignore_case = TRUE),
    "-"
  )
  
  m <- str_match(
    x,
    "^\\s*([0-9.]+)\\s*-\\s*([0-9.]+)\\s*$"
  )
  
  tibble(
    depth_top_cm = suppressWarnings(as.numeric(m[, 2])),
    depth_bottom_cm = suppressWarnings(as.numeric(m[, 3]))
  )
}

write_updated_csv <- function(x, path, ...) {
  obsolete <- intersect(
    c("transition_family", "transition_family_directional"),
    names(x)
  )
  
  if (length(obsolete) > 0) {
    x <- x |> select(-all_of(obsolete))
  }
  
  readr::write_csv(x, path, ...)
  invisible(x)
}

assert_updated_family <- function(
    df,
    object_name = "data",
    allow_empty = FALSE) {
  
  if (nrow(df) == 0 && allow_empty) return(invisible(TRUE))
  
  if (!"transition_family_k10" %in% names(df)) {
    stop(object_name, " is missing transition_family_k10.")
  }
  
  if (any(
    c("transition_family", "transition_family_directional") %in%
    names(df)
  )) {
    stop(object_name, " contains an obsolete family column.")
  }
  
  invisible(TRUE)
}

optional_numeric <- function(data, variable) {
  if (variable %in% names(data)) {
    suppressWarnings(as.numeric(data[[variable]]))
  } else {
    rep(NA_real_, nrow(data))
  }
}

optional_character <- function(data, variable) {
  if (variable %in% names(data)) {
    as.character(data[[variable]])
  } else {
    rep(NA_character_, nrow(data))
  }
}


## ================================================================
## 3. MODEL AND CR2 FUNCTIONS
## ================================================================

## Shared-control covariance for raw mean differences:
##
## yi = mean_treat - mean_control
##
## If two effects share the same control group:
##
## Cov(yi, yj) = SD_control^2 / n_control
##
## Diagonal elements remain:
##
## Var(yi) = SD_treat^2/n_treat + SD_control^2/n_control
##
## The control key is already restricted to publication and depth.
## Each shared-control block is constructed by setting its
## off-diagonal elements to the common control sampling variance.

build_V <- function(df, shared_control = TRUE) {
  
  n <- nrow(df)
  
  if (n == 0) stop("Cannot construct V for an empty dataset.")
  
  if (any(!is.finite(df$vi)) || any(df$vi <= 0)) {
    stop("All sampling variances must be finite and strictly positive.")
  }
  
  V <- diag(df$vi, nrow = n)
  
  if (!shared_control) return(V)
  
  if (!all(c("control_key", "var_ctrl") %in% names(df))) {
    stop("df must contain control_key and var_ctrl.")
  }
  
  ## Group by the unique control key. The key itself includes
  ## publication ID and depth, preventing cross-publication and
  ## cross-depth matching.
  
  control_groups <- split(
    seq_len(n),
    df$control_key
  )
  
  for (idx in control_groups) {
    
    if (length(idx) < 2) next
    
    control_variances <- unique(df$var_ctrl[idx])
    
    if (
      length(control_variances) != 1 ||
      !is.finite(control_variances[1])
    ) {
      stop(
        "Inconsistent control sampling variances within control key: ",
        df$control_key[idx[1]]
      )
    }
    
    ## Add shared-control covariance to off-diagonal cells only.
    ## The diagonal already contains treatment + control variance.
    
    V[idx, idx] <- V[idx, idx] +
      control_variances[1] *
      (matrix(1, length(idx), length(idx)) -
         diag(length(idx)))
  }
  
  V
}


fit_mv <- function(df, mods = NULL, shared_control = TRUE) {
  
  df <- df |>
    filter(
      is.finite(yi),
      is.finite(vi),
      vi > 0,
      !is.na(study_id),
      !is.na(comparison_id)
    )
  
  if (n_distinct(df$study_id) < 2) {
    stop("At least two publication clusters are required.")
  }
  
  V <- build_V(df, shared_control)
  
  if (is.null(mods)) {
    
    metafor::rma.mv(
      yi,
      V = V,
      random = ~ 1 | study_id/comparison_id,
      data = df,
      method = "REML"
    )
    
  } else {
    
    metafor::rma.mv(
      yi,
      V = V,
      mods = mods,
      random = ~ 1 | study_id/comparison_id,
      data = df,
      method = "REML"
    )
  }
}


cr2_results <- function(m, cluster, conf_level = 0.95) {
  
  rob <- clubSandwich::coef_test(
    m,
    cluster = cluster,
    vcov = "CR2"
  )
  
  df_col <- if ("df_Satt" %in% names(rob)) {
    as.numeric(rob$df_Satt)
  } else {
    rep(Inf, nrow(rob))
  }
  
  alpha <- 1 - conf_level
  
  critical <- ifelse(
    is.finite(df_col) & df_col > 0,
    qt(1 - alpha / 2, df = df_col),
    qnorm(1 - alpha / 2)
  )
  
  tibble(
    term = rownames(rob),
    estimate = as.numeric(rob$beta),
    se = as.numeric(rob$SE),
    df = df_col,
    p = as.numeric(rob$p_Satt),
    lo = as.numeric(rob$beta) - critical * as.numeric(rob$SE),
    hi = as.numeric(rob$beta) + critical * as.numeric(rob$SE)
  )
}


cr2_omnibus <- function(m, cluster, coefficient_names) {
  
  coefficient_names <- intersect(
    coefficient_names,
    names(coef(m))
  )
  
  if (length(coefficient_names) == 0) return(NULL)
  
  constraints <- clubSandwich::constrain_zero(coefficient_names)
  
  tryCatch(
    clubSandwich::Wald_test(
      m,
      constraints = constraints,
      vcov = "CR2",
      cluster = cluster,
      test = "HTZ"
    ),
    error = function(e) {
      message("Omnibus CR2 test failed: ", conditionMessage(e))
      NULL
    }
  )
}


tau2_total <- function(m) {
  sum(m$sigma2, na.rm = TRUE)
}


i2_mv_approx <- function(m, df) {
  
  tau2 <- tau2_total(m)
  mean_vi <- mean(df$vi, na.rm = TRUE)
  
  if (
    !is.finite(tau2) ||
    !is.finite(mean_vi) ||
    tau2 + mean_vi <= 0
  ) {
    return(NA_real_)
  }
  
  100 * tau2 / (tau2 + mean_vi)
}


report_mv <- function(m, df, label = "Section") {
  
  rob <- tryCatch(
    cr2_results(m, df$study_id),
    error = function(e) structure(
      list(message = conditionMessage(e)),
      class = "cr2_failure"
    )
  )
  
  est <- as.numeric(coef(m)[1])
  
  cr2_ok <- !inherits(rob, "cr2_failure") &&
    !is.null(rob) &&
    nrow(rob) >= 1 &&
    is.finite(rob$se[1]) &&
    is.finite(rob$lo[1]) &&
    is.finite(rob$hi[1])
  
  se_model <- sqrt(as.numeric(vcov(m)[1, 1]))
  
  ci_low_model <- est - qnorm(0.975) * se_model
  ci_high_model <- est + qnorm(0.975) * se_model
  
  se_cr2 <- if (cr2_ok) rob$se[1] else NA_real_
  df_cr2 <- if (cr2_ok) rob$df[1] else NA_real_
  ci_low_cr2 <- if (cr2_ok) rob$lo[1] else NA_real_
  ci_high_cr2 <- if (cr2_ok) rob$hi[1] else NA_real_
  p_cr2 <- if (cr2_ok) rob$p[1] else NA_real_
  
  QE <- suppressWarnings(
    tryCatch(as.numeric(m$QE), error = function(e) NA_real_)
  )
  
  QEp <- suppressWarnings(
    tryCatch(as.numeric(m$QEp), error = function(e) NA_real_)
  )
  
  tau2 <- tau2_total(m)
  I2 <- i2_mv_approx(m, df)
  
  out_tbl <- tibble(
    Section = label,
    pooled_SOC_change = est,
    SE_CR2 = se_cr2,
    df_CR2 = df_cr2,
    CI_low_CR2 = ci_low_cr2,
    CI_high_CR2 = ci_high_cr2,
    CI_low = ci_low_cr2,
    CI_high = ci_high_cr2,
    p_CR2 = p_cr2,
    CR2_status = if (cr2_ok) "Successful" else "Failed",
    CR2_failure_reason = if (cr2_ok) NA_character_ else
      if (inherits(rob, "cr2_failure")) rob$message else
        "CR2 returned no usable coefficient result",
    SE_model_based = se_model,
    CI_low_model_based = ci_low_model,
    CI_high_model_based = ci_high_model,
    Q = QE,
    Q_p = QEp,
    tau2_total = tau2,
    I2_approx = I2,
    k = nrow(df),
    studies = n_distinct(df$study_id)
  )
  
  direction <- case_when(
    est > 0 ~ "increases",
    est < 0 ~ "decreases",
    TRUE ~ "does not change"
  )
  
  finding <- if (cr2_ok) {
    glue(
      "{label}: Land-use change {direction} SOC stocks by ",
      "{round(est, 2)} Mg C ha-1 ",
      "(95% CR2 CI {round(ci_low_cr2, 2)} to {round(ci_high_cr2, 2)})."
    )
  } else {
    glue(
      "{label}: model estimate {round(est, 2)} Mg C ha-1; ",
      "CR2 inference unavailable. Model-based 95% CI: ",
      "{round(ci_low_model, 2)} to {round(ci_high_model, 2)}."
    )
  }
  
  variability <- glue(
    "Heterogeneity: approximate I2={round(I2, 1)}%, ",
    "tau2={round(tau2, 3)}, Q={round(QE, 2)}, p={QEp}."
  )
  
  list(
    table = out_tbl,
    finding = finding,
    variability = variability,
    cr2 = rob
  )
}


compare_moderator_ml <- function(df, moderator_formula, moderator_label) {
  
  V <- build_V(df)
  
  fit0 <- tryCatch(
    metafor::rma.mv(
      yi, V = V,
      random = ~ 1 | study_id/comparison_id,
      data = df, method = "ML"
    ),
    error = function(e) e
  )
  
  fit1 <- tryCatch(
    metafor::rma.mv(
      yi, V = V,
      mods = moderator_formula,
      random = ~ 1 | study_id/comparison_id,
      data = df, method = "ML"
    ),
    error = function(e) e
  )
  
  if (inherits(fit0, "error") || inherits(fit1, "error")) {
    return(tibble(
      moderator = moderator_label,
      k = nrow(df),
      studies = n_distinct(df$study_id),
      AIC_null = NA_real_,
      AIC_full = NA_real_,
      delta_AIC = NA_real_,
      BIC_null = NA_real_,
      BIC_full = NA_real_,
      delta_BIC = NA_real_,
      LRT = NA_real_,
      df_LRT = NA_real_,
      p_LRT = NA_real_,
      status = "ML fit failed"
    ))
  }
  
  ll0 <- logLik(fit0)
  ll1 <- logLik(fit1)
  
  ddf <- as.numeric(attr(ll1, "df") - attr(ll0, "df"))
  lrt_raw <- 2 * (as.numeric(ll1) - as.numeric(ll0))
  
  valid <- is.finite(ddf) &&
    ddf > 0 &&
    is.finite(lrt_raw) &&
    lrt_raw >= -1e-7
  
  lrt <- if (valid) max(0, lrt_raw) else NA_real_
  
  tibble(
    moderator = moderator_label,
    k = nrow(df),
    studies = n_distinct(df$study_id),
    AIC_null = AIC(fit0),
    AIC_full = AIC(fit1),
    delta_AIC = AIC(fit1) - AIC(fit0),
    BIC_null = BIC(fit0),
    BIC_full = BIC(fit1),
    delta_BIC = BIC(fit1) - BIC(fit0),
    LRT = lrt,
    df_LRT = if (valid) ddf else NA_real_,
    p_LRT = if (valid) pchisq(lrt, ddf, lower.tail = FALSE) else NA_real_,
    status = if (valid) "Successful ML comparison" else
      "AIC/BIC available; LRT unavailable or invalid"
  )
}


write_section_outputs <- function(rep_obj, out_prefix) {
  
  write_updated_csv(
    rep_obj$table,
    file.path(output_dir, paste0(out_prefix, "_table.csv"))
  )
  
  writeLines(
    c(rep_obj$finding, rep_obj$variability, ""),
    file.path(output_dir, paste0(out_prefix, "_text.txt"))
  )
}


empty_model_table <- function() {
  tibble(
    Section = character(),
    pooled_SOC_change = double(),
    SE_CR2 = double(),
    df_CR2 = double(),
    CI_low = double(),
    CI_high = double(),
    p_CR2 = double(),
    Q = double(),
    Q_p = double(),
    tau2_total = double(),
    I2_approx = double(),
    k = integer(),
    studies = integer()
  )
}


## ================================================================
## 4. READ AND VALIDATE INPUT
## ================================================================

dat_raw <- readr::read_csv(
  input_csv,
  show_col_types = FALSE,
  na = c("", "NA", "N/A")
) |>
  janitor::clean_names() |>
  mutate(
    across(
      where(is.character),
      ~ na_if(str_trim(.x), "")
    )
  )

required_columns <- c(
  "study_id",
  "land_use_before",
  "land_use_after",
  "soil_depth_cm",
  "treatment_mean",
  "treatment_sd",
  "treatment_n",
  "control_mean",
  "control_sd",
  "control_n",
  "outcome_group",
  "unit"
)

missing_columns <- setdiff(required_columns, names(dat_raw))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

unit_audit <- dat_raw |> count(unit)

write_updated_csv(
  unit_audit,
  file.path(output_dir, "diag_units.csv")
)

if (n_distinct(na.omit(dat_raw$unit)) > 1) {
  stop("Multiple measurement units detected.")
}

outcome_audit <- dat_raw |> count(outcome_group)

write_updated_csv(
  outcome_audit,
  file.path(output_dir, "diag_outcome_groups.csv")
)

dat_raw <- dat_raw |>
  filter(outcome_group == "SOC stock")

if (nrow(dat_raw) == 0) {
  stop("No SOC stock observations found.")
}


## ================================================================
## 5. CONSTRUCT ANALYSIS DATA
## ================================================================

dat <- tibble(
  source_row = seq_len(nrow(dat_raw)),
  Study_ID = as.character(dat_raw$study_id),
  study_id = extract_publication_id(dat_raw$study_id),
  comparison_id = as.character(dat_raw$study_id),
  
  country = optional_character(dat_raw, "country"),
  land_use_before = optional_character(dat_raw, "land_use_before"),
  land_use_after = optional_character(dat_raw, "land_use_after"),
  duration_years = optional_numeric(dat_raw, "duration_years"),
  soil_depth_cm = optional_character(dat_raw, "soil_depth_cm"),
  climate_zone = optional_character(dat_raw, "climate_zone"),
  soil_type = optional_character(dat_raw, "soil_type"),
  soil_texture = optional_character(dat_raw, "soil_texture"),
  outcome_group = optional_character(dat_raw, "outcome_group"),
  unit = optional_character(dat_raw, "unit"),
  
  management_type = optional_character(dat_raw, "management_type"),
  vegetation_type = optional_character(dat_raw, "vegetation_type"),
  land_use_intensity = optional_character(dat_raw, "land_use_intensity"),
  measurement_method = optional_character(dat_raw, "measurement_method"),
  depth_method = optional_character(dat_raw, "depth_method"),
  study_design = optional_character(dat_raw, "study_design"),
  
  control_group_id = optional_character(dat_raw, "control_group_id"),
  site_id = optional_character(dat_raw, "site_id"),
  experiment_id = optional_character(dat_raw, "experiment_id"),
  
  mean_treat = optional_numeric(dat_raw, "treatment_mean"),
  sd_treat = optional_numeric(dat_raw, "treatment_sd"),
  n_treat = optional_numeric(dat_raw, "treatment_n"),
  
  mean_ctrl = optional_numeric(dat_raw, "control_mean"),
  sd_ctrl = optional_numeric(dat_raw, "control_sd"),
  n_ctrl = optional_numeric(dat_raw, "control_n")
) |>
  mutate(
    transition = paste(
      land_use_before, "\u2192", land_use_after
    )
  )

class_audit <- tibble(
  variable = names(dat_raw),
  class = map_chr(dat_raw, ~ paste(class(.x), collapse = "/"))
)

write_updated_csv(
  class_audit,
  file.path(output_dir, "diag_input_variable_classes.csv")
)

missing_audit <- tibble(
  variable = names(dat_raw),
  missing_n = map_int(dat_raw, ~ sum(is.na(.x))),
  missing_percent = 100 * missing_n / nrow(dat_raw)
) |>
  arrange(desc(missing_n))

write_updated_csv(
  missing_audit,
  file.path(output_dir, "diag_missingness.csv")
)

quality_audit <- dat |>
  summarise(
    n_input = n(),
    missing_study_id = sum(is.na(study_id)),
    missing_treatment_mean = sum(!is.finite(mean_treat)),
    missing_treatment_sd = sum(!is.finite(sd_treat)),
    missing_treatment_n = sum(!is.finite(n_treat)),
    missing_control_mean = sum(!is.finite(mean_ctrl)),
    missing_control_sd = sum(!is.finite(sd_ctrl)),
    missing_control_n = sum(!is.finite(n_ctrl)),
    missing_soil_depth = sum(is.na(soil_depth_cm)),
    noninteger_sample_sizes =
      sum(is.finite(n_treat) & n_treat != floor(n_treat)) +
      sum(is.finite(n_ctrl) & n_ctrl != floor(n_ctrl)),
    zero_treatment_sd = sum(sd_treat == 0, na.rm = TRUE),
    zero_control_sd = sum(sd_ctrl == 0, na.rm = TRUE),
    negative_treatment_sd = sum(sd_treat < 0, na.rm = TRUE),
    negative_control_sd = sum(sd_ctrl < 0, na.rm = TRUE),
    treatment_n_below_2 = sum(n_treat < 2, na.rm = TRUE),
    control_n_below_2 = sum(n_ctrl < 2, na.rm = TRUE)
  )

write_updated_csv(
  quality_audit,
  file.path(output_dir, "diag_data_quality_audit.csv")
)

dat <- dat |>
  filter(
    !is.na(study_id),
    !is.na(comparison_id),
    !is.na(soil_depth_cm),
    is.finite(mean_treat),
    is.finite(mean_ctrl),
    is.finite(sd_treat),
    is.finite(sd_ctrl),
    is.finite(n_treat),
    is.finite(n_ctrl),
    n_treat >= 2,
    n_ctrl >= 2,
    n_treat == floor(n_treat),
    n_ctrl == floor(n_ctrl),
    sd_treat >= 0,
    sd_ctrl >= 0
  )

if (nrow(dat) == 0) {
  stop("No valid observations remain after cleaning.")
}

publication_audit <- dat |>
  distinct(Study_ID, study_id) |>
  arrange(study_id, Study_ID)

write_updated_csv(
  publication_audit,
  file.path(output_dir, "diag_publication_identifier_mapping.csv")
)


## ================================================================
## 6. EFFECT SIZE CALCULATION AND CONTROL IDENTIFICATION
## ================================================================

## Raw mean difference:
## yi = treatment mean - control mean
##
## Sampling variance:
## vi = SDt^2 / nt + SDc^2 / nc
##
## Control variance:
## var_ctrl = SDc^2 / nc
##
## Matching:
## - Explicit control_group_id takes precedence.
## - Otherwise, identical control mean, SD and n within the same
##   publication and depth are treated as a shared-control pattern.
##
## Numeric values are formatted consistently before key creation.

format_key_number <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    format(
      as.numeric(x),
      scientific = FALSE,
      trim = TRUE,
      digits = 15
    )
  )
}

dat_es <- dat |>
  mutate(
    yi = mean_treat - mean_ctrl,
    
    vi = sd_treat^2 / n_treat +
      sd_ctrl^2 / n_ctrl,
    
    sei = sqrt(vi),
    
    var_ctrl = sd_ctrl^2 / n_ctrl,
    
    control_key_source = if_else(
      !is.na(control_group_id) & control_group_id != "",
      "explicit_control_group_id",
      "inferred_from_identical_control_statistics"
    ),
    
    control_key = if_else(
      control_key_source == "explicit_control_group_id",
      
      paste(
        study_id,
        soil_depth_cm,
        "ID",
        control_group_id,
        sep = "|"
      ),
      
      paste(
        study_id,
        soil_depth_cm,
        "SUMMARY",
        format_key_number(mean_ctrl),
        format_key_number(sd_ctrl),
        format_key_number(n_ctrl),
        sep = "|"
      )
    )
  ) |>
  filter(
    is.finite(yi),
    is.finite(vi),
    vi > 0,
    is.finite(var_ctrl)
  ) |>
  mutate(
    duration_years_z = z_na(duration_years)
  )

if (nrow(dat_es) == 0) {
  stop("No valid effect sizes could be calculated.")
}


## ================================================================
## 7. DEPENDENCE AND SHARED-CONTROL AUDITS
## ================================================================

dependence_audit <- dat_es |>
  count(study_id, name = "k_effects") |>
  summarise(
    studies = n(),
    effects = sum(k_effects),
    min_effects_per_study = min(k_effects),
    median_effects_per_study = median(k_effects),
    max_effects_per_study = max(k_effects),
    studies_with_multiple_effects = sum(k_effects > 1)
  )

write_updated_csv(
  dependence_audit,
  file.path(output_dir, "diag_within_study_dependence.csv")
)


## Summary-statistic repetitions, including depth.

control_summary_audit <- dat_es |>
  group_by(
    study_id,
    soil_depth_cm,
    mean_ctrl,
    sd_ctrl,
    n_ctrl
  ) |>
  summarise(
    k_rows = n(),
    comparisons = paste(comparison_id, collapse = "; "),
    transitions = paste(unique(na.omit(transition)), collapse = "; "),
    source_rows = paste(source_row, collapse = "; "),
    .groups = "drop"
  ) |>
  filter(k_rows > 1) |>
  arrange(desc(k_rows))

write_updated_csv(
  control_summary_audit,
  file.path(output_dir, "diag_repeated_control_summary_statistics.csv")
)


## Row-level mapping: inspect exactly which effects share each key.

shared_control_mapping <- dat_es |>
  group_by(control_key) |>
  mutate(
    control_group_k = n(),
    shared_control_group = control_group_k > 1
  ) |>
  ungroup() |>
  select(
    source_row,
    Study_ID,
    study_id,
    comparison_id,
    soil_depth_cm,
    transition,
    control_group_id,
    control_key_source,
    control_key,
    mean_ctrl,
    sd_ctrl,
    n_ctrl,
    var_ctrl,
    control_group_k,
    shared_control_group
  ) |>
  arrange(study_id, control_key, source_row)

write_updated_csv(
  shared_control_mapping,
  file.path(output_dir, "diag_shared_control_mapping_row_level.csv")
)


## Group-level summary.

shared_control_key_audit <- dat_es |>
  group_by(
    study_id,
    soil_depth_cm,
    control_key,
    control_key_source,
    mean_ctrl,
    sd_ctrl,
    n_ctrl,
    var_ctrl
  ) |>
  summarise(
    k_effects = n(),
    comparisons = paste(comparison_id, collapse = "; "),
    transitions = paste(unique(na.omit(transition)), collapse = "; "),
    source_rows = paste(source_row, collapse = "; "),
    .groups = "drop"
  ) |>
  mutate(
    shared_control = k_effects > 1
  ) |>
  arrange(desc(k_effects), study_id)

write_updated_csv(
  shared_control_key_audit,
  file.path(output_dir, "diag_shared_control_key_audit.csv")
)


## Keep the original audit object name for downstream compatibility.

shared_control_audit <- shared_control_key_audit |>
  filter(shared_control) |>
  transmute(
    study_id,
    control_group_id = control_key,
    k = k_effects,
    soil_depth_cm,
    control_key_source,
    comparisons,
    transitions
  )

write_updated_csv(
  shared_control_audit,
  file.path(output_dir, "diag_potential_shared_controls.csv")
)


## ================================================================
## 8. TRANSITION FAMILY CLASSIFICATION
## ================================================================

dat_es <- dat_es |>
  mutate(
    lu_before_class = canon_landuse(land_use_before),
    lu_after_class = canon_landuse(land_use_after),
    
    transition_family_directional = paste(
      str_to_sentence(str_squish(land_use_before)),
      "\u2192",
      str_to_sentence(str_squish(land_use_after))
    )
  )

classification_audit <- dat_es |>
  distinct(
    land_use_before,
    lu_before_class,
    land_use_after,
    lu_after_class,
    transition_family_directional
  ) |>
  arrange(transition_family_directional)

write_updated_csv(
  classification_audit,
  file.path(output_dir, "diag_landuse_classification_audit.csv")
)

family_support <- dat_es |>
  group_by(transition_family_directional) |>
  summarise(
    k = n(),
    studies = n_distinct(study_id),
    .groups = "drop"
  ) |>
  arrange(desc(k))

write_updated_csv(
  family_support,
  file.path(output_dir, "diag_transition_family_support.csv")
)

large_families <- family_support |>
  filter(k > 10) |>
  pull(transition_family_directional)

other_family <- "Other / low-k transitions"

family_levels <- c(large_families, other_family)

dat_es <- dat_es |>
  mutate(
    transition_family_k10 = if_else(
      transition_family_directional %in% large_families,
      transition_family_directional,
      other_family
    ),
    transition_family_k10 = factor(
      transition_family_k10,
      levels = family_levels
    )
  )

write_updated_csv(
  dat_es |> count(transition_family_k10, sort = TRUE),
  file.path(output_dir, "diag_transition_family_k_gt_10_counts.csv")
)

write_updated_csv(
  dat_es |>
    filter(transition_family_k10 == other_family) |>
    count(transition, transition_family_k10, sort = TRUE),
  file.path(output_dir, "diag_transitions_in_other_low_k.csv")
)

dat_es_export <- dat_es |>
  select(-transition_family_directional)

assert_updated_family(dat_es_export, "effect-size export")

write_updated_csv(
  dat_es_export,
  file.path(output_dir, "effect_sizes_MD_soc_change.csv")
)


## ================================================================
## 9. OVERALL MODEL
## ================================================================

m_overall <- fit_mv(dat_es)

rep_overall <- report_mv(
  m_overall,
  dat_es,
  "Overall effect: all land-use transitions"
)

write_section_outputs(rep_overall, "report_overall")

rob_overall <- cr2_results(m_overall, dat_es$study_id)

m_overall_diag <- fit_mv(
  dat_es,
  shared_control = FALSE
)

rep_overall_diag <- report_mv(
  m_overall_diag,
  dat_es,
  "Overall effect: diagonal V (shared controls ignored)"
)

write_updated_csv(
  bind_rows(
    rep_overall$table,
    rep_overall_diag$table
  ),
  file.path(output_dir, "table_overall_shared_control_vs_diagonal_V.csv")
)

write_updated_csv(
  rob_overall,
  file.path(output_dir, "table_overall_model_CR2_coefficients.csv")
)


## ================================================================
## 10. ASSUMED WITHIN-PUBLICATION CORRELATION SENSITIVITY
## ================================================================

## These are alternative working covariance structures.
## They are NOT added to the known shared-control covariance.

make_working_V <- function(vi, cluster, rho) {
  
  stopifnot(
    length(vi) == length(cluster),
    rho >= 0,
    rho < 1
  )
  
  V <- diag(vi, nrow = length(vi))
  
  for (ii in split(seq_along(vi), cluster)) {
    if (length(ii) > 1L) {
      block <- rho * outer(sqrt(vi[ii]), sqrt(vi[ii]))
      diag(block) <- vi[ii]
      V[ii, ii] <- block
    }
  }
  
  V
}

rho_values <- c(0, 0.25, 0.50, 0.75)

rho_sensitivity_models <- purrr::map(
  rho_values,
  function(rho) {
    
    V_rho <- make_working_V(
      dat_es$vi,
      dat_es$study_id,
      rho
    )
    
    fit <- tryCatch(
      metafor::rma.mv(
        yi,
        V = V_rho,
        random = ~ 1 | study_id/comparison_id,
        data = dat_es,
        method = "REML"
      ),
      error = function(e) e
    )
    
    list(rho = rho, V = V_rho, fit = fit)
  }
)

names(rho_sensitivity_models) <- paste0("rho_", rho_values)

rho_sensitivity_results <- purrr::map_dfr(
  rho_sensitivity_models,
  function(x) {
    
    if (inherits(x$fit, "error")) {
      return(tibble(
        rho = x$rho,
        estimate = NA_real_,
        SE_CR2 = NA_real_,
        df_CR2 = NA_real_,
        CI_low_CR2 = NA_real_,
        CI_high_CR2 = NA_real_,
        p_CR2 = NA_real_,
        CR2_status = "Model fitting failed",
        failure_reason = conditionMessage(x$fit),
        k = nrow(dat_es),
        studies = n_distinct(dat_es$study_id)
      ))
    }
    
    cr <- tryCatch(
      cr2_results(x$fit, dat_es$study_id),
      error = function(e) e
    )
    
    ok <- !inherits(cr, "error") &&
      nrow(cr) >= 1 &&
      is.finite(cr$se[1])
    
    tibble(
      rho = x$rho,
      estimate = as.numeric(coef(x$fit)[1]),
      SE_CR2 = if (ok) cr$se[1] else NA_real_,
      df_CR2 = if (ok) cr$df[1] else NA_real_,
      CI_low_CR2 = if (ok) cr$lo[1] else NA_real_,
      CI_high_CR2 = if (ok) cr$hi[1] else NA_real_,
      p_CR2 = if (ok) cr$p[1] else NA_real_,
      CR2_status = if (ok) "Successful" else "CR2 failed",
      failure_reason = if (ok) NA_character_ else
        if (inherits(cr, "error")) conditionMessage(cr) else
          "No usable CR2 result",
      tau2_total = tau2_total(x$fit),
      k = nrow(dat_es),
      studies = n_distinct(dat_es$study_id)
    )
  }
)

write_updated_csv(
  rho_sensitivity_results,
  file.path(output_dir, "table_assumed_within_publication_rho_sensitivity.csv")
)

saveRDS(
  rho_sensitivity_models,
  file.path(output_dir, "models_assumed_within_publication_rho_sensitivity.rds")
)


## ================================================================
## 11. TRANSITION-FAMILY MODELS
## ================================================================

family_results <- empty_model_table()
family_text <- character()

for (fam in family_levels) {
  
  dsub <- dat_es |>
    filter(transition_family_k10 == fam)
  
  if (
    nrow(dsub) < 3 ||
    n_distinct(dsub$study_id) < 2
  ) next
  
  m_fam <- tryCatch(
    fit_mv(dsub),
    error = function(e) {
      message("Family model failed for ", fam, ": ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(m_fam)) next
  
  rep_fam <- report_mv(
    m_fam,
    dsub,
    paste("Transition family:", fam)
  )
  
  family_results <- bind_rows(
    family_results,
    rep_fam$table |>
      mutate(transition_family_k10 = as.character(fam))
  )
  
  family_text <- c(
    family_text,
    rep_fam$finding,
    rep_fam$variability,
    ""
  )
}

assert_updated_family(
  family_results,
  "family_results",
  allow_empty = TRUE
)

write_updated_csv(
  family_results,
  file.path(output_dir, "report_transition_family_sections_table.csv")
)

writeLines(
  family_text,
  file.path(output_dir, "report_transition_family_sections_text.txt")
)


## ================================================================
## 12. FACTOR MODERATOR MODELS
## ================================================================

fit_factor_moderator <- function(
    df,
    fvar,
    min_k = 10,
    min_studies = 5,
    min_k_level = 3,
    max_levels = 25) {
  
  d <- df |>
    filter(
      !is.na(.data[[fvar]]),
      is.finite(yi),
      is.finite(vi)
    )
  
  if (
    nrow(d) < min_k ||
    n_distinct(d$study_id) < min_studies
  ) return(NULL)
  
  level_counts <- d |>
    count(.data[[fvar]], name = "k_level") |>
    mutate(level_chr = as.character(.data[[fvar]])) |>
    filter(
      !is.na(level_chr),
      level_chr != "",
      k_level >= min_k_level
    ) |>
    arrange(desc(k_level))
  
  if (nrow(level_counts) < 2) return(NULL)
  
  if (nrow(level_counts) > max_levels) {
    level_counts <- level_counts |> slice_head(n = max_levels)
  }
  
  d <- d |>
    filter(as.character(.data[[fvar]]) %in% level_counts$level_chr) |>
    mutate(
      moderator_value = droplevels(factor(.data[[fvar]]))
    )
  
  fml <- ~ moderator_value - 1
  
  m <- tryCatch(
    metafor::rma.mv(
      yi,
      V = build_V(d),
      mods = fml,
      random = ~ 1 | study_id/comparison_id,
      data = d,
      method = "REML"
    ),
    error = function(e) {
      message("Factor model failed for ", fvar, ": ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(m)) return(NULL)
  
  tt <- tryCatch(
    cr2_results(m, d$study_id),
    error = function(e) NULL
  )
  
  if (is.null(tt)) return(NULL)
  
  tt <- tt |>
    mutate(
      level = str_remove(term, "^moderator_value"),
      level = str_squish(level)
    )
  
  omnibus <- cr2_omnibus(m, d$study_id, names(coef(m)))
  
  p_omnibus <- if (
    !is.null(omnibus) &&
    "p_val" %in% names(omnibus)
  ) {
    as.numeric(omnibus$p_val[1])
  } else {
    NA_real_
  }
  
  tt <- tt |>
    mutate(
      p_omnibus = p_omnibus,
      moderator = fvar,
      k = nrow(d),
      studies = n_distinct(d$study_id)
    )
  
  p <- ggplot(
    tt,
    aes(
      x = estimate,
      y = fct_reorder(level, estimate)
    )
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbar(
      aes(xmin = lo, xmax = hi),
      orientation = "y",
      width = 0.2
    ) +
    geom_point(aes(size = 1 / se), shape = 21, fill = "white") +
    scale_size_continuous(
      range = c(2, 6),
      name = "Precision\n(1/SE)"
    ) +
    labs(
      title = paste("SOC change by", fvar),
      subtitle = paste("CR2 omnibus p =", fmt_p(p_omnibus)),
      x = "SOC change (Mg C ha-1)",
      y = NULL
    ) +
    theme_pub(11)
  
  ml_comparison <- compare_moderator_ml(d, fml, fvar)
  
  list(
    model = m,
    table = tt,
    plot = p,
    ml_comparison = ml_comparison,
    analysis_data = d
  )
}

factor_vars <- c(
  "transition_family_k10",
  "climate_zone",
  "soil_type",
  "soil_texture",
  "management_type",
  "vegetation_type",
  "land_use_intensity",
  "measurement_method",
  "depth_method",
  "study_design",
  "country"
)

factor_results <- tibble()
factor_ml_comparisons <- tibble()

for (fvar in factor_vars) {
  
  out <- fit_factor_moderator(dat_es, fvar)
  
  if (is.null(out)) next
  
  fig_name <- paste0("figure_factor_", safe_file(fvar))
  tab_name <- paste0("table_factor_", safe_file(fvar), "_effects_CR2")
  
  ggsave(
    file.path(output_dir, paste0(fig_name, ".png")),
    out$plot,
    width = 7.8,
    height = 5.4,
    dpi = 600,
    bg = "white"
  )
  
  saveRDS(
    out$plot,
    file.path(output_dir, paste0(fig_name, ".rds"))
  )
  
  write_updated_csv(
    out$table,
    file.path(output_dir, paste0(tab_name, ".csv"))
  )
  
  factor_results <- bind_rows(factor_results, out$table)
  factor_ml_comparisons <- bind_rows(
    factor_ml_comparisons,
    out$ml_comparison
  )
}

write_updated_csv(
  factor_results,
  file.path(output_dir, "table_factor_all_effects_CR2.csv")
)

factor_ml_comparisons <- factor_ml_comparisons |>
  mutate(p_LRT_holm = p.adjust(p_LRT, method = "holm"))

write_updated_csv(
  factor_ml_comparisons,
  file.path(output_dir, "table_factor_model_comparisons_ML.csv")
)


## ================================================================
## 13. DEPTH INTERVAL AUDIT AND ANALYSIS
## ================================================================

if (
  all(c("soil_depth_top_cm", "soil_depth_bottom_cm") %in% names(dat_raw))
) {
  
  depth_top <- suppressWarnings(as.numeric(dat_raw$soil_depth_top_cm))
  depth_bottom <- suppressWarnings(as.numeric(dat_raw$soil_depth_bottom_cm))
  depth_source <- "explicit top/bottom columns"
  
} else if ("soil_depth_cm" %in% names(dat_raw)) {
  
  parsed_depth <- parse_depth_interval(dat_raw$soil_depth_cm)
  depth_top <- parsed_depth$depth_top_cm
  depth_bottom <- parsed_depth$depth_bottom_cm
  depth_source <- "parsed soil_depth_cm"
  
} else {
  
  depth_top <- rep(NA_real_, nrow(dat_raw))
  depth_bottom <- rep(NA_real_, nrow(dat_raw))
  depth_source <- "no depth data available"
}

depth_lookup <- tibble(
  source_row = seq_len(nrow(dat_raw)),
  depth_top_cm = depth_top,
  depth_bottom_cm = depth_bottom
)

dat_es_depth <- dat_es |>
  left_join(depth_lookup, by = "source_row") |>
  mutate(
    depth_valid = is.finite(depth_top_cm) &
      is.finite(depth_bottom_cm) &
      depth_top_cm >= 0 &
      depth_bottom_cm > depth_top_cm,
    
    depth_interval = if_else(
      depth_valid,
      paste0(depth_top_cm, "-", depth_bottom_cm, " cm"),
      NA_character_
    ),
    
    depth_width_cm = if_else(
      depth_valid,
      depth_bottom_cm - depth_top_cm,
      NA_real_
    ),
    
    depth_bin = depth_interval,
    
    is_cumulative_profile = depth_valid &
      depth_top_cm == 0 &
      depth_bottom_cm > 20
  )

write_updated_csv(
  dat_es_depth |>
    distinct(
      soil_depth_cm,
      depth_top_cm,
      depth_bottom_cm,
      depth_interval,
      depth_bin,
      is_cumulative_profile
    ),
  file.path(output_dir, "diag_depth_interval_audit.csv")
)

overlap_depth_audit <- dat_es_depth |>
  filter(depth_valid) |>
  group_by(study_id) |>
  summarise(
    has_cumulative_profile = any(is_cumulative_profile),
    n_depth_intervals = n_distinct(depth_interval),
    depth_intervals = paste(
      sort(unique(depth_interval)),
      collapse = "; "
    ),
    .groups = "drop"
  ) |>
  filter(has_cumulative_profile, n_depth_intervals > 1)

write_updated_csv(
  overlap_depth_audit,
  file.path(output_dir, "diag_cumulative_and_multiple_depths.csv")
)

message("Depth classification source: ", depth_source)

depth_results <- empty_model_table() |>
  mutate(
    depth_bin = character(),
    transition_family_k10 = character()
  )

depth_text <- character()

depth_levels <- dat_es_depth |>
  filter(!is.na(depth_bin)) |>
  distinct(depth_bin, depth_top_cm, depth_bottom_cm) |>
  arrange(depth_top_cm, depth_bottom_cm) |>
  pull(depth_bin)

for (db in depth_levels) {
  
  for (fam in family_levels) {
    
    dsub <- dat_es_depth |>
      filter(
        depth_bin == db,
        transition_family_k10 == fam
      )
    
    if (
      nrow(dsub) < 3 ||
      n_distinct(dsub$study_id) < 2
    ) next
    
    m_depth <- tryCatch(
      fit_mv(dsub),
      error = function(e) NULL
    )
    
    if (is.null(m_depth)) next
    
    rep_depth <- report_mv(
      m_depth,
      dsub,
      paste("Transition family:", fam, "|", db)
    )
    
    depth_results <- bind_rows(
      depth_results,
      rep_depth$table |>
        mutate(
          depth_bin = db,
          transition_family_k10 = as.character(fam)
        )
    )
    
    depth_text <- c(
      depth_text,
      rep_depth$finding,
      rep_depth$variability,
      ""
    )
  }
}

assert_updated_family(
  depth_results,
  "depth_results",
  allow_empty = TRUE
)

write_updated_csv(
  depth_results,
  file.path(output_dir, "report_transition_family_by_depth_table.csv")
)

writeLines(
  depth_text,
  file.path(output_dir, "report_transition_family_by_depth_text.txt")
)


## ================================================================
## 14. NUMERIC MODERATOR META-REGRESSION
## ================================================================

numeric_meta_results <- tibble()
numeric_ml_comparisons <- tibble()

numeric_vars <- c("duration_years_z")

for (xvar in numeric_vars) {
  
  d <- dat_es |>
    filter(
      is.finite(yi),
      is.finite(vi),
      is.finite(.data[[xvar]])
    )
  
  if (
    nrow(d) < 10 ||
    n_distinct(d$study_id) < 5
  ) next
  
  m_num <- tryCatch(
    fit_mv(d, mods = as.formula(paste("~", xvar))),
    error = function(e) {
      message("Meta-regression failed for ", xvar, ": ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(m_num)) next
  
  rob <- tryCatch(
    cr2_results(m_num, d$study_id),
    error = function(e) NULL
  )
  
  if (is.null(rob)) next
  
  numeric_meta_results <- bind_rows(
    numeric_meta_results,
    rob |>
      mutate(
        moderator = xvar,
        k = nrow(d),
        studies = n_distinct(d$study_id)
      )
  )
  
  p <- ggplot(d, aes(x = .data[[xvar]], y = yi)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(
      aes(size = 1 / vi),
      alpha = 0.6,
      shape = 21,
      fill = "white"
    ) +
    geom_smooth(method = "lm", se = TRUE, linetype = "dotted") +
    scale_size_continuous(
      range = c(1, 6),
      name = "Precision\n(1/vi)"
    ) +
    labs(
      title = paste("Observed effects by", xvar),
      subtitle = "Descriptive line only; meta-regression is reported separately",
      x = xvar,
      y = "SOC change (Mg C ha-1)"
    ) +
    theme_pub()
  
  ggsave(
    file.path(output_dir, paste0("figure_bubble_", safe_file(xvar), ".png")),
    p,
    width = 7,
    height = 5,
    dpi = 600,
    bg = "white"
  )
  
  numeric_ml_comparisons <- bind_rows(
    numeric_ml_comparisons,
    compare_moderator_ml(
      d,
      as.formula(paste("~", xvar)),
      xvar
    )
  )
}

write_updated_csv(
  numeric_meta_results,
  file.path(output_dir, "table_numeric_meta_regression_CR2.csv")
)

numeric_ml_comparisons <- numeric_ml_comparisons |>
  mutate(p_LRT_holm = p.adjust(p_LRT, method = "holm"))

write_updated_csv(
  numeric_ml_comparisons,
  file.path(output_dir, "table_numeric_model_comparisons_ML.csv")
)

all_moderator_ml_comparisons <- bind_rows(
  factor_ml_comparisons,
  numeric_ml_comparisons
) |>
  mutate(
    p_LRT_holm_all_moderators = p.adjust(p_LRT, method = "holm")
  )

write_updated_csv(
  all_moderator_ml_comparisons,
  file.path(output_dir, "table_all_moderator_model_comparisons_ML.csv")
)


## ================================================================
## 15. LEAVE-ONE-PUBLICATION-OUT SENSITIVITY
## ================================================================

loco_results <- tibble(
  left_out = character(),
  est = double(),
  se_CR2 = double(),
  df_CR2 = double(),
  lo = double(),
  hi = double(),
  p_CR2 = double(),
  studies_remaining = integer(),
  k_remaining = integer()
)

for (sid in unique(dat_es$study_id)) {
  
  dsub <- dat_es |> filter(study_id != sid)
  
  if (
    nrow(dsub) < 3 ||
    n_distinct(dsub$study_id) < 2
  ) next
  
  msub <- tryCatch(
    fit_mv(dsub),
    error = function(e) NULL
  )
  
  if (is.null(msub)) next
  
  rob <- tryCatch(
    cr2_results(msub, dsub$study_id),
    error = function(e) NULL
  )
  
  if (
    is.null(rob) ||
    nrow(rob) == 0 ||
    !is.finite(rob$estimate[1])
  ) next
  
  loco_results <- bind_rows(
    loco_results,
    tibble(
      left_out = as.character(sid),
      est = rob$estimate[1],
      se_CR2 = rob$se[1],
      df_CR2 = rob$df[1],
      lo = rob$lo[1],
      hi = rob$hi[1],
      p_CR2 = rob$p[1],
      studies_remaining = n_distinct(dsub$study_id),
      k_remaining = nrow(dsub)
    )
  )
}

overall_est <- rep_overall$table$pooled_SOC_change
overall_lo <- rep_overall$table$CI_low
overall_hi <- rep_overall$table$CI_high

if (nrow(loco_results) > 0) {
  
  p_loco <- ggplot(
    loco_results,
    aes(x = reorder(left_out, est), y = est)
  ) +
    geom_hline(yintercept = overall_est, linetype = "dashed") +
    geom_point() +
    coord_flip() +
    labs(
      title = "Leave-one-publication-out sensitivity",
      subtitle = "CR2 estimates after omitting each publication",
      x = "Publication omitted",
      y = "SOC change (Mg C ha-1)"
    ) +
    theme_pub()
  
  ggsave(
    file.path(output_dir, "figure_leave_one_out.png"),
    p_loco,
    width = 8,
    height = 6,
    dpi = 600,
    bg = "white"
  )
  
  loo_summary <- tibble(
    pooled = overall_est,
    pooled_lo = overall_lo,
    pooled_hi = overall_hi,
    loo_min = min(loco_results$est, na.rm = TRUE),
    loo_max = max(loco_results$est, na.rm = TRUE),
    loo_range = diff(range(loco_results$est, na.rm = TRUE)),
    studies_analysed = n_distinct(dat_es$study_id),
    successful_loo_models = nrow(loco_results)
  )
  
} else {
  
  loo_summary <- tibble(
    pooled = overall_est,
    pooled_lo = overall_lo,
    pooled_hi = overall_hi,
    loo_min = NA_real_,
    loo_max = NA_real_,
    loo_range = NA_real_,
    studies_analysed = n_distinct(dat_es$study_id),
    successful_loo_models = 0L
  )
}

write_updated_csv(
  loo_summary,
  file.path(output_dir, "table_leave_one_out_summary.csv")
)

write_updated_csv(
  loco_results,
  file.path(output_dir, "table_leave_one_out_estimates.csv")
)


## ================================================================
## 16. EXPLORATORY FUNNEL PLOTS
## ================================================================

funnel_results <- tibble(
  transition_family_k10 = character(),
  k_effect_sizes = integer(),
  k_studies = integer(),
  diagnostic = character(),
  formal_test_performed = logical(),
  note = character()
)

for (fam in setdiff(family_levels, other_family)) {
  
  d <- dat_es |>
    filter(
      transition_family_k10 == fam,
      is.finite(yi),
      is.finite(vi),
      vi > 0
    )
  
  n_study <- n_distinct(d$study_id)
  
  if (nrow(d) < 5 || n_study < 5) next
  
  m <- tryCatch(fit_mv(d), error = function(e) NULL)
  if (is.null(m)) next
  
  png(
    file.path(
      output_dir,
      paste0("funnel_exploratory_", safe_file(fam), ".png")
    ),
    width = 1400,
    height = 1050,
    res = 220
  )
  
  metafor::funnel(
    x = d$yi,
    sei = sqrt(d$vi),
    main = paste("Exploratory funnel plot:", fam),
    xlab = "SOC change (Mg C ha-1)"
  )
  
  abline(v = as.numeric(coef(m)[1]), lty = 2)
  dev.off()
  
  funnel_results <- bind_rows(
    funnel_results,
    tibble(
      transition_family_k10 = as.character(fam),
      k_effect_sizes = nrow(d),
      k_studies = n_study,
      diagnostic = "Exploratory funnel plot only",
      formal_test_performed = FALSE,
      note = paste(
        "Potentially dependent effects are present.",
        "Visual asymmetry is not proof of publication bias."
      )
    )
  )
}

write_updated_csv(
  funnel_results,
  file.path(output_dir, "table_funnel_diagnostics_by_transition_family.csv")
)


## ================================================================
## 17. OVERALL RESULTS BY TRANSITION FAMILY
## ================================================================

write_updated_csv(
  family_results,
  file.path(output_dir, "table_overall_by_transition_family.csv")
)


## ================================================================
## 18. SOFTWARE VERSIONS
## ================================================================

package_versions <- tibble(
  Package = required_packages,
  Version = map_chr(
    required_packages,
    ~ as.character(packageVersion(.x))
  )
)

write_updated_csv(
  package_versions,
  file.path(output_dir, "table_software_versions.csv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)


## ================================================================
## 19. PUBLICATION TABLE HELPERS
## ================================================================

make_ft <- function(df, caption, font_size = 9) {
  
  flextable(df) |>
    set_caption(caption) |>
    theme_booktabs() |>
    fontsize(size = font_size, part = "all") |>
    font(fontname = "Arial", part = "all") |>
    align(align = "center", part = "all") |>
    align(j = 1, align = "left", part = "all") |>
    autofit()
}


## ================================================================
## 20. PUBLICATION TABLE 1: EVIDENCE BASE
## ================================================================

n_effects <- nrow(dat_es)
n_studies <- n_distinct(dat_es$study_id)

duration_range <- safe_range(dat_es$duration_years)
depth_range <- safe_range(dat_es_depth$depth_bottom_cm)

climate_n <- sum(!is.na(dat_es$climate_zone))
soil_n <- sum(!is.na(dat_es$soil_type))
texture_n <- sum(!is.na(dat_es$soil_texture))

climate_top <- dat_es |>
  filter(!is.na(climate_zone)) |>
  count(climate_zone, sort = TRUE) |>
  slice_head(n = 3) |>
  transmute(txt = paste0(climate_zone, " (n=", n, ")")) |>
  pull(txt) |>
  paste(collapse = "; ")

soil_top <- dat_es |>
  filter(!is.na(soil_type)) |>
  count(soil_type, sort = TRUE) |>
  slice_head(n = 3) |>
  transmute(txt = paste0(soil_type, " (n=", n, ")")) |>
  pull(txt) |>
  paste(collapse = "; ")

texture_top <- dat_es |>
  filter(!is.na(soil_texture)) |>
  count(soil_texture, sort = TRUE) |>
  slice_head(n = 3) |>
  transmute(txt = paste0(soil_texture, " (n=", n, ")")) |>
  pull(txt) |>
  paste(collapse = "; ")

table1_summary <- tibble(
  Component = c(
    "Effect sizes",
    "Publication clusters",
    "Effect-size metric",
    "Duration range",
    "Missing duration",
    "Sampling-depth bottom range",
    "Climate-zone information",
    "Main climatic zones",
    "Soil-type information",
    "Main soil types",
    "Soil-texture information",
    "Main soil textures"
  ),
  Summary = c(
    as.character(n_effects),
    as.character(n_studies),
    "Raw mean difference in SOC stock (Mg C ha-1)",
    ifelse(
      all(is.na(duration_range)),
      "Not available",
      paste(duration_range[1], "to", duration_range[2], "years")
    ),
    as.character(sum(is.na(dat_es$duration_years))),
    ifelse(
      all(is.na(depth_range)),
      "Not available",
      paste(depth_range[1], "to", depth_range[2], "cm")
    ),
    paste0(
      climate_n, " (",
      fmt_num(100 * climate_n / n_effects, 1), "%)"
    ),
    climate_top,
    paste0(soil_n, " (", fmt_num(100 * soil_n / n_effects, 1), "%)"),
    soil_top,
    paste0(
      texture_n, " (",
      fmt_num(100 * texture_n / n_effects, 1), "%)"
    ),
    texture_top
  )
)

table1_families <- dat_es |>
  count(transition_family_k10, name = "n") |>
  mutate(
    Component = paste("Transition family:", transition_family_k10),
    Summary = paste0(
      n, " (", fmt_num(100 * n / sum(n), 1), "%)"
    )
  ) |>
  select(Component, Summary)

table1 <- bind_rows(table1_summary, table1_families)

write_updated_csv(
  table1,
  file.path(pub_dir, "Table_1_evidence_base.csv")
)

save_as_docx(
  make_ft(
    table1,
    "Table 1. Evidence base and transition-family structure."
  ),
  path = file.path(pub_dir, "Table_1_evidence_base.docx")
)


## ================================================================
## 21. PUBLICATION TABLE 2: OVERALL AND FAMILY EFFECTS
## ================================================================

overall_pub <- rep_overall$table |>
  mutate(Group = "Overall")

family_pub <- family_results |>
  mutate(Group = as.character(transition_family_k10))

table2 <- bind_rows(overall_pub, family_pub) |>
  transmute(
    Group,
    `SOC change, Mg C ha-1 (95% CI)` =
      fmt_ci(pooled_SOC_change, CI_low, CI_high),
    k,
    Studies = studies,
    Q = fmt_num(Q, 2),
    `p(Q)` = fmt_p(Q_p),
    `I2 (%)` = fmt_num(I2_approx, 1)
  )

write_updated_csv(
  table2,
  file.path(pub_dir, "Table_2_transition_family_effects.csv")
)

save_as_docx(
  make_ft(
    table2,
    "Table 2. Overall and transition-family effects on SOC stock.",
    8.5
  ),
  path = file.path(pub_dir, "Table_2_transition_family_effects.docx")
)


## ================================================================
## 22. PUBLICATION TABLE 3: DEPTH-STRATIFIED EFFECTS
## ================================================================

if (nrow(depth_results) > 0) {
  
  table3 <- depth_results |>
    arrange(
      match(depth_bin, depth_levels),
      transition_family_k10
    ) |>
    transmute(
      Depth = depth_bin,
      `Transition family` = transition_family_k10,
      `SOC change, Mg C ha-1 (95% CI)` =
        fmt_ci(pooled_SOC_change, CI_low, CI_high),
      k,
      Studies = studies,
      Q = fmt_num(Q, 2),
      `p(Q)` = fmt_p(Q_p),
      `I2 (%)` = fmt_num(I2_approx, 1)
    )
  
} else {
  
  table3 <- tibble(
    Depth = character(),
    `Transition family` = character(),
    `SOC change, Mg C ha-1 (95% CI)` = character(),
    k = integer(),
    Studies = integer(),
    Q = character(),
    `p(Q)` = character(),
    `I2 (%)` = character()
  )
}

write_updated_csv(
  table3,
  file.path(pub_dir, "Table_3_depth_stratified_effects.csv")
)

if (nrow(table3) > 0) {
  save_as_docx(
    make_ft(
      table3,
      "Table 3. Depth-interval-stratified SOC effects.",
      8
    ),
    path = file.path(pub_dir, "Table_3_depth_stratified_effects.docx")
  )
}


## ================================================================
## 23. PUBLICATION TABLE 4: FUNNEL DIAGNOSTICS
## ================================================================

table4 <- funnel_results |>
  transmute(
    `Transition family` = transition_family_k10,
    `Effect sizes` = k_effect_sizes,
    `Publications` = k_studies,
    Diagnostic = diagnostic,
    `Formal test performed` = formal_test_performed,
    Interpretation = note
  )

write_updated_csv(
  table4,
  file.path(pub_dir, "Table_4_funnel_diagnostics.csv")
)

if (nrow(table4) > 0) {
  save_as_docx(
    make_ft(
      table4,
      "Table 4. Exploratory funnel-plot diagnostics.",
      8
    ),
    path = file.path(pub_dir, "Table_4_funnel_diagnostics.docx")
  )
}


## ================================================================
## 24. PUBLICATION TABLE S1: LEAVE-ONE-PUBLICATION-OUT
## ================================================================

table_s1 <- loo_summary |>
  transmute(
    `Full pooled estimate` = fmt_num(pooled, 2),
    `95% CI lower` = fmt_num(pooled_lo, 2),
    `95% CI upper` = fmt_num(pooled_hi, 2),
    `Minimum leave-one-out estimate` = fmt_num(loo_min, 2),
    `Maximum leave-one-out estimate` = fmt_num(loo_max, 2),
    `Leave-one-out range` = fmt_num(loo_range, 2)
  )

write_updated_csv(
  table_s1,
  file.path(pub_dir, "Table_S1_leave_one_out_summary.csv")
)

save_as_docx(
  make_ft(
    table_s1,
    "Table S1. Leave-one-publication-out sensitivity analysis."
  ),
  path = file.path(pub_dir, "Table_S1_leave_one_out_summary.docx")
)


## ================================================================
## 25. FIGURE 1: TRANSITION-FAMILY EFFECTS
## ================================================================

family_plot_results <- family_results |>
  filter(as.character(transition_family_k10) != other_family) |>
  mutate(
    transition_family_k10 = droplevels(
      factor(transition_family_k10)
    )
  )

if (nrow(family_plot_results) > 0) {
  
  fig1 <- ggplot(
    family_plot_results,
    aes(
      x = pooled_SOC_change,
      y = fct_reorder(transition_family_k10, pooled_SOC_change)
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey40"
    ) +
    geom_errorbar(
      aes(xmin = CI_low, xmax = CI_high),
      orientation = "y",
      width = 0.2
    ) +
    geom_point(aes(size = 1 / SE_CR2), shape = 21, fill = "white") +
    scale_size_continuous(
      range = c(2, 6),
      name = "Precision\n(1/SE)"
    ) +
    labs(
      title = "SOC response by transition family",
      x = "SOC change (Mg C ha-1)",
      y = NULL
    ) +
    theme_pub(11)
  
  ggsave(
    file.path(pub_dir, "Fig_1_transition_family.tiff"),
    fig1,
    width = 180,
    height = 130,
    units = "mm",
    dpi = 600,
    compression = "lzw"
  )
  
  ggsave(
    file.path(pub_dir, "Fig_1_transition_family.png"),
    fig1,
    width = 180,
    height = 130,
    units = "mm",
    dpi = 600,
    bg = "white"
  )
}


## ================================================================
## 26. FIGURE 2: DEPTH-STRATIFIED EFFECTS
## ================================================================

depth_plot_results <- depth_results |>
  filter(as.character(transition_family_k10) != other_family) |>
  mutate(
    transition_family_k10 = droplevels(
      factor(transition_family_k10)
    ),
    depth_bin = factor(depth_bin, levels = depth_levels)
  )

if (nrow(depth_plot_results) > 0) {
  
  fig2 <- ggplot(
    depth_plot_results,
    aes(
      x = pooled_SOC_change,
      y = fct_reorder(transition_family_k10, pooled_SOC_change)
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey40"
    ) +
    geom_errorbar(
      aes(xmin = CI_low, xmax = CI_high),
      orientation = "y",
      width = 0.2
    ) +
    geom_point(aes(size = 1 / SE_CR2), shape = 21, fill = "white") +
    facet_wrap(~ depth_bin, scales = "free_y", ncol = 1) +
    scale_size_continuous(
      range = c(2, 5),
      name = "Precision\n(1/SE)"
    ) +
    labs(
      title = "Depth-interval-stratified SOC responses",
      x = "SOC change (Mg C ha-1)",
      y = NULL
    ) +
    theme_pub(10)
  
  ggsave(
    file.path(pub_dir, "Fig_2_depth_stratified.tiff"),
    fig2,
    width = 180,
    height = 210,
    units = "mm",
    dpi = 600,
    compression = "lzw"
  )
  
  ggsave(
    file.path(pub_dir, "Fig_2_depth_stratified.png"),
    fig2,
    width = 180,
    height = 210,
    units = "mm",
    dpi = 600,
    bg = "white"
  )
}


## ================================================================
## 27. FIGURE 3: CLIMATE AND SOIL MODERATORS
## ================================================================

climate_file <- file.path(
  output_dir,
  "table_factor_climate_zone_effects_CR2.csv"
)

soil_file <- file.path(
  output_dir,
  "table_factor_soil_type_effects_CR2.csv"
)

if (file.exists(climate_file) && file.exists(soil_file)) {
  
  climate_plot_data <- read_csv(
    climate_file,
    show_col_types = FALSE
  )
  
  soil_plot_data <- read_csv(
    soil_file,
    show_col_types = FALSE
  )
  
  moderator_plot <- function(df, plot_title) {
    
    if (
      nrow(df) == 0 ||
      !all(c("estimate", "lo", "hi", "se", "level") %in% names(df))
    ) {
      return(NULL)
    }
    
    df <- df |>
      filter(
        is.finite(estimate),
        is.finite(lo),
        is.finite(hi),
        is.finite(se),
        se > 0
      ) |>
      mutate(
        Level = level,
        precision = 1 / se
      )
    
    if (nrow(df) == 0) return(NULL)
    
    ggplot(
      df,
      aes(x = estimate, y = fct_reorder(Level, estimate))
    ) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_errorbar(
        aes(xmin = lo, xmax = hi),
        orientation = "y",
        width = 0.2
      ) +
      geom_point(aes(size = precision), shape = 21, fill = "white") +
      scale_size_continuous(
        range = c(2, 5),
        name = "Precision\n(1/SE)"
      ) +
      labs(
        title = plot_title,
        x = "SOC change (Mg C ha-1)",
        y = NULL
      ) +
      theme_pub(10)
  }
  
  p_climate <- moderator_plot(climate_plot_data, "Climate zone")
  p_soil <- moderator_plot(soil_plot_data, "Soil type")
  
  if (!is.null(p_climate) && !is.null(p_soil)) {
    
    fig3 <- p_climate / p_soil +
      plot_annotation(tag_levels = "a")
    
    ggsave(
      file.path(pub_dir, "Fig_3_climate_soil.tiff"),
      fig3,
      width = 180,
      height = 210,
      units = "mm",
      dpi = 600,
      compression = "lzw"
    )
    
    ggsave(
      file.path(pub_dir, "Fig_3_climate_soil.png"),
      fig3,
      width = 180,
      height = 210,
      units = "mm",
      dpi = 600,
      bg = "white"
    )
  }
}


## ================================================================
## 28. FIGURE 4: LEAVE-ONE-PUBLICATION-OUT
## ================================================================

if (nrow(loco_results) > 0) {
  
  fig4 <- ggplot(loco_results, aes(x = est, y = 1)) +
    geom_rect(
      data = tibble(
        xmin = overall_lo,
        xmax = overall_hi,
        ymin = -Inf,
        ymax = Inf
      ),
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      inherit.aes = FALSE,
      alpha = 0.08
    ) +
    geom_vline(xintercept = overall_est, linetype = "dashed") +
    geom_point(alpha = 0.6) +
    scale_y_continuous(NULL, breaks = NULL) +
    labs(
      title = "Leave-one-publication-out sensitivity",
      x = "SOC change after omitting one publication (Mg C ha-1)"
    ) +
    theme_pub(11)
  
  ggsave(
    file.path(pub_dir, "Fig_4_leave_one_out.tiff"),
    fig4,
    width = 180,
    height = 80,
    units = "mm",
    dpi = 600,
    compression = "lzw"
  )
  
  ggsave(
    file.path(pub_dir, "Fig_4_leave_one_out.png"),
    fig4,
    width = 180,
    height = 80,
    units = "mm",
    dpi = 600,
    bg = "white"
  )
}


## ================================================================
## 29. COPY SUPPORTING FIGURES
## ================================================================

supporting_figures <- list.files(
  output_dir,
  pattern = "^figure_.*\\.(png|tiff|pdf)$",
  full.names = TRUE,
  recursive = FALSE
)

for (f in supporting_figures) {
  file.copy(
    f,
    file.path(pub_dir, paste0("supporting_", basename(f))),
    overwrite = TRUE
  )
}


## ================================================================
## 30. FINAL CSV VALIDATION
## ================================================================

all_csv_files <- list.files(
  output_dir,
  pattern = "\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

legacy_csv_files <- character()

for (csv_file in all_csv_files) {
  
  csv_names <- tryCatch(
    names(read_csv(
      csv_file,
      n_max = 0,
      show_col_types = FALSE
    )),
    error = function(e) character()
  )
  
  if (any(
    c("transition_family", "transition_family_directional") %in%
    csv_names
  )) {
    legacy_csv_files <- c(legacy_csv_files, csv_file)
  }
}

if (length(legacy_csv_files) > 0) {
  stop(
    "Legacy family columns detected in CSV exports:\n",
    paste(legacy_csv_files, collapse = "\n")
  )
}

effect_export <- read_csv(
  file.path(output_dir, "effect_sizes_MD_soc_change.csv"),
  n_max = 0,
  show_col_types = FALSE
)

if (!"transition_family_k10" %in% names(effect_export)) {
  stop("Effect-size export is missing transition_family_k10.")
}


## ================================================================
## 31. FINAL MODEL AND DATA DIAGNOSTICS
## ================================================================

shared_key_counts <- dat_es |>
  count(control_key, name = "k") |>
  filter(k > 1)

model_diagnostics <- tibble(
  input_rows = nrow(dat_raw),
  
  valid_effect_sizes = nrow(dat_es),
  
  publication_clusters = n_distinct(dat_es$study_id),
  
  transition_families_k_gt_10 = length(large_families),
  
  repeated_control_summary_patterns = nrow(control_summary_audit),
  
  identified_shared_control_groups = nrow(shared_control_audit),
  
  explicit_control_id_rows = sum(
    dat_es$control_key_source == "explicit_control_group_id"
  ),
  
  inferred_control_key_rows = sum(
    dat_es$control_key_source ==
      "inferred_from_identical_control_statistics"
  ),
  
  effects_in_shared_control_groups = sum(
    dat_es$control_key %in% shared_key_counts$control_key
  ),
  
  shared_control_groups_in_model = nrow(shared_key_counts),
  
  rows_with_valid_depth_intervals = sum(dat_es_depth$depth_valid),
  
  rows_with_unparsed_depth_intervals = sum(!dat_es_depth$depth_valid),
  
  cumulative_depth_rows = sum(
    dat_es_depth$is_cumulative_profile,
    na.rm = TRUE
  ),
  
  sampling_covariance_off_diagonal_specified = TRUE
)

write_updated_csv(
  model_diagnostics,
  file.path(output_dir, "diag_final_model_and_data_summary.csv")
)


## ================================================================
## 32. REPRODUCIBILITY BUNDLE AND RUN RECORD
## ================================================================

bundle_names <- c(
  "dat_raw",
  "dat_es",
  "dat_es_depth",
  "m_overall",
  "m_overall_diag",
  "rep_overall",
  "rep_overall_diag",
  "family_results",
  "factor_results",
  "factor_ml_comparisons",
  "numeric_meta_results",
  "numeric_ml_comparisons",
  "all_moderator_ml_comparisons",
  "rho_sensitivity_results",
  "rho_sensitivity_models",
  "loco_results",
  "funnel_results",
  "depth_results",
  "control_summary_audit",
  "shared_control_audit",
  "shared_control_mapping",
  "shared_control_key_audit",
  "overlap_depth_audit",
  "model_diagnostics"
)

analysis_bundle <- mget(
  intersect(bundle_names, ls(envir = .GlobalEnv)),
  envir = .GlobalEnv
)

analysis_bundle$session_info <- utils::sessionInfo()

analysis_bundle$input_file <- normalizePath(
  input_csv,
  mustWork = TRUE
)

analysis_bundle$created_at <- Sys.time()

analysis_bundle$method_notes <- c(
  "Primary effect size: raw mean difference in Mg C ha-1.",
  "Primary model: multilevel REML with shared-control sampling covariance.",
  "Primary robust inference: CR2 clustered by publication.",
  paste(
    "Shared controls matched using explicit control_group_id when available;",
    "otherwise by publication, soil depth, control mean, control SD and n."
  ),
  paste(
    "Summary-statistic-based control matching is inferred;",
    "inspect shared-control audit files."
  ),
  paste(
    "Assumed-rho models are alternative working covariance sensitivity analyses,",
    "not additive to the shared-control covariance."
  ),
  "ML nested-model comparisons are complementary to CR2 inference.",
  "Small-study-effect diagnostics remain exploratory because effect sizes are dependent.",
  "Publication identifiers extracted by removing the final numeric suffix.",
  "Repeated-measure covariance is not explicitly modelled.",
  "Cross-depth covariance for overlapping soil profiles is not explicitly modelled."
)

saveRDS(
  analysis_bundle,
  file.path(output_dir, "soc_meta_analysis_reproducibility_bundle.rds")
)

capture.output(
  utils::sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)

writeLines(
  c(
    "SOC meta-analysis pipeline completed.",
    paste("Input:", normalizePath(input_csv)),
    paste("Effect sizes:", nrow(dat_es)),
    paste("Publication clusters:", n_distinct(dat_es$study_id)),
    "Primary inference: publication-clustered CR2.",
    "Primary sampling covariance: shared-control covariance.",
    "Sensitivity: diagonal V and assumed within-publication rho working structures.",
    "Moderator comparisons: ML; CR2 remains the robust inferential analysis.",
    "Review shared-control identity, publication mapping and repeated-measure dependence before publication."
  ),
  file.path(output_dir, "RUN_NOTES.txt")
)


## ================================================================
## 33. FINAL OUTPUT SUMMARY
## ================================================================

message("")
message("====================================================")
message("META-ANALYSIS PIPELINE COMPLETED")
message("====================================================")

message("Input rows: ", nrow(dat_raw))
message("Valid effect sizes: ", nrow(dat_es))
message("Publication clusters: ", n_distinct(dat_es$study_id))

message(
  "Transition families retained with k > 10: ",
  length(large_families)
)

message(
  "Transition-family levels: ",
  paste(family_levels, collapse = "; ")
)

message("")
message("Authoritative exported family variable: transition_family_k10")

message(
  "Sampling covariance: shared-control covariance included. ",
  "Diagonal-V sensitivity model saved. ",
  "Repeated-measure and cross-depth covariance are not modelled."
)

message("CR2 inference: clustered at publication level.")

message(
  "Shared-control groups inferred from identical control statistics: ",
  sum(
    shared_control_key_audit$shared_control &
      shared_control_key_audit$control_key_source ==
      "inferred_from_identical_control_statistics"
  )
)

message("Legacy family columns in CSV exports: NONE")
message("CSV validation passed.")
message("")
message("Main outputs: ", output_dir)
message("Publication-ready outputs: ", pub_dir)
message("")
message(
  "IMPORTANT: Review publication identifier mapping, shared-control ",
  "mapping, repeated-control audit and depth-overlap audit before publication."
)

message("====================================================")