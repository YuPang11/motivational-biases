

# ---- 0) Setup ----
suppressPackageStartupMessages(library(rstan))

# Your paths
rootdir    <- "F:/PIT/model/EEGmodelHC/HC"
datadir    <- file.path(rootdir, "results")
resultsdir <- file.path(datadir, "stan")
r2dir      <- file.path(resultsdir, "R2")
dir.create(r2dir, showWarnings = FALSE, recursive = TRUE)

# Evaluate only these models
model_list <- c(1, 2, 3, 4, 5)

# ---- 1) Compute and save R2 ----
get_R2_for_model <- function(m) {
  stanfitsfile <- file.path(resultsdir, sprintf("M%dfits_obsData.Rdata", m))
  if (!file.exists(stanfitsfile)) {
    warning("Missing file: ", stanfitsfile)
    return(data.frame(model = sprintf("M%d", m), R2 = NA_real_, R2_percent = NA_real_))
  }
  
  # Load saved fit results (should contain dList and f2)
  load(stanfitsfile)  # -> f2, dList, mstr
  
  # log_lik is available only if calcLogLik = TRUE was used during fitting
  ext <- try(rstan::extract(f2, "log_lik"), silent = TRUE)
  if (inherits(ext, "try-error") || is.null(ext$log_lik)) {
    warning("log_lik not found in file: ", stanfitsfile, " (calcLogLik=TRUE is required during fitting)")
    return(data.frame(model = sprintf("M%d", m), R2 = NA_real_, R2_percent = NA_real_))
  }
  
  # Expected likelihood (mean hit probability across trials)
  likelihood <- colMeans(exp(ext$log_lik))  # dimension: nIter x nData -> column means
  nData <- length(likelihood)
  
  # Total variance
  # Original author: baseline error per trial = 2/3; equivalent to 1 - 1/nResp
  base   <- 1 - 1 / dList$nResp      # for your data, nResp=3 -> base=2/3
  SS_tot <- base^2 * nData
  SS_res <- sum((1 - likelihood)^2)
  
  R2  <- 1 - (SS_res / SS_tot)
  R2p <- 100 * R2
  
  # Save
  save(R2p, file = file.path(r2dir, sprintf("M%d_R2.Rdata", m)))
  
  data.frame(model = sprintf("M%d", m), R2 = R2, R2_percent = R2p)
}

r2_table <- do.call(rbind, lapply(model_list, get_R2_for_model))

# ---- 2) Print summary ----
print(r2_table)




# ---- 0) Setup ----
suppressPackageStartupMessages(library(rstan))

# Your paths
rootdir    <- "F:/PIT/model/EEGmodelHC/OCD"
datadir    <- file.path(rootdir, "results")
resultsdir <- file.path(datadir, "stan", "OCD1")
r2dir      <- file.path(resultsdir, "R2")
dir.create(r2dir, showWarnings = FALSE, recursive = TRUE)

# Evaluate only these models
model_list <- c(1, 2, 3, 4, 5)

# ---- 1) Compute and save R2 ----
get_R2_for_model <- function(m) {
  stanfitsfile <- file.path(resultsdir, sprintf("M%dfits_obsData.Rdata", m))
  if (!file.exists(stanfitsfile)) {
    warning("Missing file: ", stanfitsfile)
    return(data.frame(model = sprintf("M%d", m), R2 = NA_real_, R2_percent = NA_real_))
  }
  
  # Load saved fit results (should contain dList and f2)
  load(stanfitsfile)  # -> f2, dList, mstr
  
  # log_lik is available only if calcLogLik = TRUE was used during fitting
  ext <- try(rstan::extract(f2, "log_lik"), silent = TRUE)
  if (inherits(ext, "try-error") || is.null(ext$log_lik)) {
    warning("log_lik not found in file: ", stanfitsfile, " (calcLogLik=TRUE is required during fitting)")
    return(data.frame(model = sprintf("M%d", m), R2 = NA_real_, R2_percent = NA_real_))
  }
  
  # Expected likelihood (mean hit probability across trials)
  likelihood <- colMeans(exp(ext$log_lik))  # dimension: nIter x nData -> column means
  nData <- length(likelihood)
  
  # Total variance
  # Original author: baseline error per trial = 2/3; equivalent to 1 - 1/nResp
  base   <- 1 - 1 / dList$nResp      # for your data, nResp=3 -> base=2/3
  SS_tot <- base^2 * nData
  SS_res <- sum((1 - likelihood)^2)
  
  R2  <- 1 - (SS_res / SS_tot)
  R2p <- 100 * R2
  
  # Save
  save(R2p, file = file.path(r2dir, sprintf("M%d_R2.Rdata", m)))
  
  data.frame(model = sprintf("M%d", m), R2 = R2, R2_percent = R2p)
}

r2_table <- do.call(rbind, lapply(model_list, get_R2_for_model))

# ---- 2) Print summary ----
print(r2_table)


