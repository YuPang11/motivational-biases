# =================== Top-level posteriors & summary (M9/10/12/13/14) ===================
suppressPackageStartupMessages({
  library(rstan)
})

# ---- Paths ----
rootdir    <- "F:/PIT/model/EEGmodelHC/HC"
datadir    <- file.path(rootdir, "results")
resultsdir <- file.path(datadir, "stan")
outdir     <- file.path(resultsdir, "posteriors")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Models to process ----
model_list <- c(9, 10, 12, 13, 14, 15, 16, 19, 21)

# ---- Helper: safe inv_logit (avoids dependency on boot) ----
inv_logit <- function(x) 1/(1+exp(-x))

# ---- Helper: parameter names (based on your current model-family convention) ----
# K=5:  rho, epsilon, go_bias, pav_bias, kappa
# K=6:  the above 5 plus betaEEG (in the 6th dimension)
get_param_names <- function(K) {
  if (K == 5) c("rho","epsilon","go_bias","pav_bias","kappa")
  else if (K == 6) c("rho","epsilon","go_bias","pav_bias","kappa","betaEEG")
  else paste0("X", seq_len(K))
}

# ---- Main loop ----
for (m in model_list) {
  cat("\n================  M", m, "  ================\n", sep = "")
  stanfitsfile   <- file.path(resultsdir, sprintf("M%dfits_obsData.Rdata", m))
  posteriorsfile <- file.path(outdir,     sprintf("M%d_posteriors_obsData.Rdata", m))
  
  if (!file.exists(stanfitsfile)) {
    warning("Missing fit file: ", stanfitsfile)
    next
  }
  
  # Load fit objects
  load(stanfitsfile)  # -> f2, dList, mstr (based on your save format)
  
  # Extract posterior samples from vector parameter X: iterations x K
  ext <- rstan::extract(f2, pars = "X")
  if (is.list(ext) && !is.null(ext$X)) {
    posteriors <- ext$X
  } else {
    # Fallback: use as.matrix() and take the first K columns
    # (not recommended, but kept as a backup)
    draws <- as.matrix(f2)
    posteriors <- draws[, 1:dList$K, drop = FALSE]
  }
  
  K <- ncol(posteriors)
  pnames <- get_param_names(K)
  colnames(posteriors) <- pnames
  
  # Save posterior samples for later reuse
  save(posteriors, file = posteriorsfile)
  
  # ---- Basic summary statistics ----
  means <- apply(posteriors, 2, mean)
  sds   <- apply(posteriors, 2, sd)
  ppos  <- apply(posteriors, 2, function(x) mean(x > 0) * 100)
  
  summary_tab <- data.frame(param = pnames,
                            mean  = as.numeric(means),
                            sd    = as.numeric(sds),
                            pct_pos = as.numeric(ppos))
  print(summary_tab, row.names = FALSE)
  
  # ---- Common parameter transforms (computed if the relevant dimensions exist) ----
  if ("rho" %in% pnames) {
    rho_exp_mean <- mean(exp(posteriors[,"rho"]))
    rho_exp_sd   <- sd(exp(posteriors[,"rho"]))
    cat(sprintf("exp(rho): mean=%.4f, sd=%.4f\n", rho_exp_mean, rho_exp_sd))
  }
  
  if ("epsilon" %in% pnames) {
    eps_logit_mean <- mean(inv_logit(posteriors[,"epsilon"]))
    eps_logit_sd   <- sd(inv_logit(posteriors[,"epsilon"]))
    cat(sprintf("logit(epsilon): mean=%.4f, sd=%.4f\n", eps_logit_mean, eps_logit_sd))
  }
  
  if (all(c("epsilon","pav_bias") %in% pnames)) {
    e_plus_pi  <- inv_logit(posteriors[,"epsilon"] + posteriors[,"pav_bias"])
    e_minus_pi <- inv_logit(posteriors[,"epsilon"] - posteriors[,"pav_bias"])
    cat(sprintf("logit(epsilon + pav_bias): mean=%.4f, sd=%.4f\n",
                mean(e_plus_pi), sd(e_plus_pi)))
    cat(sprintf("logit(epsilon - pav_bias): mean=%.4f, sd=%.4f\n",
                mean(e_minus_pi), sd(e_minus_pi)))
  }
  
  # Optional: if kappa is treated as a rate in your model, report its logistic transform
  if ("kappa" %in% pnames) {
    kap_logit_mean <- mean(inv_logit(posteriors[,"kappa"]))
    kap_logit_sd   <- sd(inv_logit(posteriors[,"kappa"]))
    cat(sprintf("logit(kappa) [if needed]: mean=%.4f, sd=%.4f\n", kap_logit_mean, kap_logit_sd))
  }
}