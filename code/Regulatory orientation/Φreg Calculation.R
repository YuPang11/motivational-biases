## =========================================================
## HC + OCD: extract subject-level β, compute PCI / QCI / Φreg,
## and save results tables only
## HC: use M16 + M19 only
## OCD: use M16 + M21 only
## No plotting
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(rstan)
})

## ------------------------- 1. Common functions -------------------------

# Extract subject-level β (x6[*]) from stanfit .Rdata
get_beta_draws <- function(rda_path, par_regex = "^x6\\[") {
  e <- new.env()
  load(rda_path, envir = e)
  stopifnot(exists("f2", e))
  
  M <- as.matrix(e$f2)
  sel <- grep(par_regex, colnames(M), value = TRUE)
  
  if (!length(sel)) {
    stop("No x6[*] (β) columns found in ", basename(rda_path))
  }
  
  draws <- M[, sel, drop = FALSE]
  
  list(
    draws = draws,
    mean  = colMeans(draws),
    sd    = apply(draws, 2, sd),
    p_neg = colMeans(draws < 0),
    p_pos = colMeans(draws > 0)
  )
}

# z-standardization
zfun <- function(x) as.numeric(scale(x))

# Build metrics table
build_metrics_df <- function(
    beta_pi, sd_pi, p_pi_neg, p_pi_pos,
    beta_q,  sd_q,  p_q_neg,  p_q_pos,
    qci_sign = -1,
    q_col_name = "betaQ",
    model_label = NA_character_,
    type1_use_q = c("neg"),
    type2_use_q = c("pos")
) {
  stopifnot(length(beta_pi) == length(beta_q))
  nSub <- length(beta_pi)
  
  q_use1 <- match.arg(type1_use_q[1], c("neg", "pos"))
  q_use2 <- match.arg(type2_use_q[1], c("neg", "pos"))
  
  q_prob1 <- if (q_use1 == "neg") p_q_neg else p_q_pos
  q_prob2 <- if (q_use2 == "neg") p_q_neg else p_q_pos
  
  df <- tibble::tibble(
    sub        = 1:nSub,
    beta16     = as.numeric(beta_pi),
    sd_beta16  = as.numeric(sd_pi),
    p16_neg    = as.numeric(p_pi_neg),
    p16_pos    = as.numeric(p_pi_pos),
    betaQ      = as.numeric(beta_q),
    sd_betaQ   = as.numeric(sd_q),
    pQ_neg     = as.numeric(p_q_neg),
    pQ_pos     = as.numeric(p_q_pos),
    model_used = model_label
  )
  
  names(df)[names(df) == "betaQ"]    <- q_col_name
  names(df)[names(df) == "sd_betaQ"] <- paste0("sd_", q_col_name)
  names(df)[names(df) == "pQ_neg"]   <- paste0("p_", q_col_name, "_neg")
  names(df)[names(df) == "pQ_pos"]   <- paste0("p_", q_col_name, "_pos")
  
  beta_q_col   <- q_col_name
  sd_q_col     <- paste0("sd_", q_col_name)
  
  df <- df %>%
    mutate(
      PCI       = -beta16,
      QCI       = qci_sign * .data[[beta_q_col]],
      PCI_z     = zfun(PCI),
      QCI_z     = zfun(QCI),
      H_NGtW    = PCI_z + QCI_z,
      H_GtA     = PCI_z - QCI_z,
      Balance   = PCI_z - QCI_z,
      Opponency = -beta16 * .data[[beta_q_col]],
      R_mag     = sqrt(PCI_z^2 + QCI_z^2),
      `Φreg`    = atan2(QCI_z, PCI_z) * 180 / pi,
      w         = 1 / (sd_beta16^2 + .data[[sd_q_col]]^2),
      Type1_prob   = p16_neg * q_prob1,
      Type2_prob   = p16_pos * q_prob2,
      TypeContrast = Type1_prob - Type2_prob
    )
  
  df
}

## ------------------------- 2. HC group: M16 + M19  -------------------------

hc_fit16 <- "F:/PIT/model/EEGmodelHC/HC/results/stan/M16fits_obsData.Rdata"
hc_fit19 <- "F:/PIT/model/EEGmodelHC/HC/results/stan/M19fits_obsData.Rdata"

hc_out_dir <- "F:/PIT/model/EEGmodelHC/HC/plots/βR"
if (!dir.exists(hc_out_dir)) dir.create(hc_out_dir, recursive = TRUE)

stopifnot(file.exists(hc_fit16), file.exists(hc_fit19))

# Extract β
hc_ext16 <- get_beta_draws(hc_fit16)
hc_ext19 <- get_beta_draws(hc_fit19)

# Compute metrics
# Original HC logic:
# QCI = -beta19
# Type1 = p16_neg * p19_neg
# Type2 = p16_pos * p19_pos
hc_df <- build_metrics_df(
  beta_pi    = hc_ext16$mean,
  sd_pi      = hc_ext16$sd,
  p_pi_neg   = hc_ext16$p_neg,
  p_pi_pos   = hc_ext16$p_pos,
  beta_q     = hc_ext19$mean,
  sd_q       = hc_ext19$sd,
  p_q_neg    = hc_ext19$p_neg,
  p_q_pos    = hc_ext19$p_pos,
  qci_sign   = -1,
  q_col_name = "beta19",
  model_label = "M16 + M19",
  type1_use_q = "neg",
  type2_use_q = "pos"
)

hc_out_csv <- file.path(hc_out_dir, "subject_metrics_HC_M16_M19_Phireg.csv")
readr::write_csv(hc_df, hc_out_csv)

## ------------------------- 3. OCD group: M16 + M21 -------------------------

ocd_fit16 <- "F:/PIT/model/EEGmodelHC/OCD/results/stan/M16fits_obsData.Rdata"
ocd_fit21 <- "F:/PIT/model/EEGmodelHC/OCD/results/stan/M21fits_obsData.Rdata"

ocd_out_dir <- "F:/PIT/model/EEGmodelHC/OCD/plots/βR"
if (!dir.exists(ocd_out_dir)) dir.create(ocd_out_dir, recursive = TRUE)

stopifnot(file.exists(ocd_fit16), file.exists(ocd_fit21))

# Extract β
ocd_ext16 <- get_beta_draws(ocd_fit16)
ocd_ext21 <- get_beta_draws(ocd_fit21)

# Compute metrics
# Original OCD logic:
# QCI = -beta21
# Type1 = p16_neg * p21_neg
# Type2 = p16_pos * p21_pos
ocd_df <- build_metrics_df(
  beta_pi    = ocd_ext16$mean,
  sd_pi      = ocd_ext16$sd,
  p_pi_neg   = ocd_ext16$p_neg,
  p_pi_pos   = ocd_ext16$p_pos,
  beta_q     = ocd_ext21$mean,
  sd_q       = ocd_ext21$sd,
  p_q_neg    = ocd_ext21$p_neg,
  p_q_pos    = ocd_ext21$p_pos,
  qci_sign   = -1,
  q_col_name = "beta21",
  model_label = "M16 + M21",
  type1_use_q = "neg",
  type2_use_q = "pos"
)

ocd_out_csv <- file.path(ocd_out_dir, "subject_metrics_OCD_M16_M21_Phireg.csv")
readr::write_csv(ocd_df, ocd_out_csv)

## ------------------------- 4. Completion message -------------------------

message("✅ Done. Results tables saved only; no plots generated.")
message("HC output: ", hc_out_csv)
message("OCD output: ", ocd_out_csv)