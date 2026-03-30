rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(rstan)
  library(matrixStats)
  library(dplyr)
})

# =========================================================
# 0) Path settings: modify these paths as needed
# =========================================================

# ---- M5 group-level posterior export ----
path_m5_hc_fit  <- "F:/PIT/model/EEGmodelHC/HC/results/stan/M5/M5fits_obsData.Rdata"
path_m5_ocd_fit <- "F:/PIT/model/EEGmodelHC/OCD/results/stan/M5/M5fits_obsData.Rdata"

out_plot_dir <- "F:/PIT/model/EEGmodelHC/HC/results/stan/plot"
if (!dir.exists(out_plot_dir)) dir.create(out_plot_dir, recursive = TRUE)

# ---- HC: M5 posterior predictive ----
hc_dir      <- "F:/PIT/model/EEGmodelHC/HC/results/stan/M5"
hc_fit_file <- file.path(hc_dir, "M5fits_obsData.Rdata")
hc_obs_file <- file.path(hc_dir, "obs4stanhc3.Rdata")   # Change here if the filename is different

# ---- OCD: M5 posterior predictive ----
ocd_dir      <- "F:/PIT/model/EEGmodelHC/OCD/results/stan/M5"
ocd_fit_file <- file.path(ocd_dir, "M5fits_obsData.Rdata")
ocd_obs_file <- file.path(ocd_dir, "obs4stanocd3.Rdata") # Change here if the filename is different

# =========================================================
# 1) Utility functions
# =========================================================

inv_logit <- function(z) 1 / (1 + exp(-z))

# -------- Extract group-level X[1:5] for M5_posteriors_wide.csv --------
load_group_X <- function(rdata_path) {
  if (!file.exists(rdata_path)) stop("File not found: ", rdata_path)
  
  e <- new.env()
  load(rdata_path, envir = e)
  
  if (!exists("f2", envir = e)) {
    stop("Object 'f2' not found in: ", rdata_path)
  }
  
  S <- as.data.frame(as.matrix(e$f2))
  Xcols <- grep("^X\\[[1-5]\\]$", colnames(S), perl = TRUE)
  
  if (length(Xcols) != 5) {
    stop("Did not find exactly 5 group-level X[1:5] columns in: ", rdata_path)
  }
  
  Sx <- S[, Xcols, drop = FALSE]
  names(Sx) <- c("rho", "epsilon", "go_bias", "pav_bias", "kappa")
  Sx
}

transform_m5_group <- function(df) {
  df %>%
    mutate(
      rho     = exp(rho),
      epsilon = inv_logit(epsilon)
    )
}

# -------- Extract subject-level x1..x5 from f2 by parameter name --------
extract_subject_draws_M5 <- function(f2, nSub) {
  mat <- as.matrix(f2)
  
  pull_block <- function(prefix) {
    idx <- grep(paste0("^", prefix, "\\["), colnames(mat))
    if (length(idx) == 0) stop("Parameter columns not found: ", prefix)
    nm  <- colnames(mat)[idx]
    ord <- order(as.integer(sub(".*\\[(\\d+)\\]$", "\\1", nm)))
    idx <- idx[ord]
    out <- mat[, idx, drop = FALSE]
    if (ncol(out) != nSub) {
      stop(prefix, " returned ", ncol(out), " subject columns, but expected ", nSub)
    }
    out
  }
  
  list(
    rho     = pull_block("x1"),
    epsilon = pull_block("x2"),
    gobias  = pull_block("x3"),
    pibias  = pull_block("x4"),
    kappa   = pull_block("x5")
  )
}

# -------- M5 simulation function: behavioral model (shared by HC/OCD) --------
M5_sim <- function(fbsens, lr, gb, pi, kap,
                   stimSeq, respSeq, fbSeq,
                   nTrial, nResp, nStim, nPerm) {
  
  nPerm <- min(nPerm, nrow(as.matrix(fbsens)))
  
  allp  <- array(NA_real_, dim = c(nPerm, nTrial, nResp))
  confl <- array(NA_real_, dim = c(nPerm, nTrial))
  
  for (iIter in 1:nPerm) {
    rho_i <- exp(fbsens[iIter])
    eps_i <- inv_logit(lr[iIter])
    
    if (eps_i < 0.5) {
      bias_lr    <- inv_logit(c(NA, lr[iIter] - kap[iIter]))  # (Go, NoGo)
      bias_lr[1] <- 2 * eps_i - bias_lr[2]
    } else {
      bias_lr    <- inv_logit(c(lr[iIter] + kap[iIter], NA))
      bias_lr[2] <- 2 * eps_i - bias_lr[1]
    }
    
    # Initialize Q and V
    Q <- matrix(rho_i * c(.5, .5, -.5, -.5), nrow = nResp, ncol = nStim, byrow = TRUE)
    V <- rep(c(.5, .5, -.5, -.5), nStim / 4)
    
    valenced <- rep(0, nStim)
    q <- matrix(0, nrow = nTrial, ncol = nResp)
    
    for (t in 1:nTrial) {
      s <- stimSeq[t]
      
      # Action weights for the current trial
      q[t, ]    <- valenced[s] * Q[, s]
      q[t, 1:2] <- q[t, 1:2] + gb[iIter] + pi[iIter] * valenced[s] * V[s]
      
      # Model-based conflict (not standardized)
      confl[iIter, t] <- -valenced[s] * (V[s] * (mean(Q[1:2, s]) - Q[3, s]))
      
      # Softmax
      p <- exp(q[t, ]) / sum(exp(q[t, ]))
      p[is.nan(p)] <- 1
      allp[iIter, t, ] <- p
      
      # Q update
      resp <- respSeq[t]
      out  <- fbSeq[t]
      iEps <- ifelse(resp %in% c(1, 2), 1, 2)  # 1=Go, 2=NoGo
      
      if ((iEps == 1 && out == 1) || (iEps == 2 && out == -1)) {
        Q[resp, s] <- Q[resp, s] + bias_lr[iEps] * (rho_i * out - Q[resp, s])
      } else {
        Q[resp, s] <- Q[resp, s] + eps_i * (rho_i * out - Q[resp, s])
      }
      
      if (out != 0) valenced[s] <- 1
    }
  }
  
  list(
    meanp = colMeans(allp, na.rm = TRUE),                 # nTrial × nResp
    confl = matrixStats::colMedians(confl, na.rm = TRUE) # Median across posterior samples
  )
}

# -------- Export M5 pGo / pCorrect for one group --------
export_M5_predictions <- function(fit_file, obs_file, out_dir, group_name = "HC") {
  cat("\n=== Running ", group_name, " M5 posterior predictive export ===\n", sep = "")
  
  if (!file.exists(fit_file)) stop("Missing fit file: ", fit_file)
  if (!file.exists(obs_file)) stop("Missing obs file: ", obs_file)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Load fitted model
  e_fit <- new.env()
  load(fit_file, envir = e_fit)
  if (!exists("f2", envir = e_fit)) stop("Object 'f2' not found in: ", fit_file)
  f2 <- e_fit$f2
  
  # Load dList
  e_obs <- new.env()
  load(obs_file, envir = e_obs)
  if (!exists("dList", envir = e_obs)) stop("Object 'dList' not found in: ", obs_file)
  dList <- e_obs$dList
  
  nSub   <- if (!is.null(dList$nSub))   as.integer(dList$nSub)   else nrow(dList$s)
  nStim  <- if (!is.null(dList$nStim))  as.integer(dList$nStim)  else length(unique(as.vector(dList$s)))
  nResp  <- if (!is.null(dList$nResp))  as.integer(dList$nResp)  else 3L
  nTrial <- if (!is.null(dList$nTrial)) as.integer(dList$nTrial) else ncol(dList$s)
  
  stopifnot(nResp == 3L)
  
  pars <- extract_subject_draws_M5(f2, nSub)
  rho     <- pars$rho
  epsilon <- pars$epsilon
  gobias  <- pars$gobias
  pibias  <- pars$pibias
  kappa   <- pars$kappa
  
  # Number of posterior samples
  nPerm <- nrow(rho)
  
  meanp <- array(NA_real_, dim = c(nSub, nTrial, nResp))
  confl <- array(NA_real_, dim = c(nSub, nTrial))
  
  for (iSub in 1:nSub) {
    cat(group_name, " subject ", iSub, "/", nSub, "\n", sep = "")
    out <- M5_sim(
      fbsens  = rho[, iSub, drop = FALSE],
      lr      = epsilon[, iSub, drop = FALSE],
      gb      = gobias[, iSub, drop = FALSE],
      pi      = pibias[, iSub, drop = FALSE],
      kap     = kappa[, iSub, drop = FALSE],
      stimSeq = dList$s[iSub, ],
      respSeq = dList$ya[iSub, ],
      fbSeq   = dList$r[iSub, ],
      nTrial  = nTrial,
      nResp   = nResp,
      nStim   = nStim,
      nPerm   = nPerm
    )
    meanp[iSub, , ] <- out$meanp
    confl[iSub, ]   <- out$confl
  }
  
  # Correct response code for each stimulus
  correctResp <- c(1, 2, 1, 2, 3, 3, 3, 3)
  
  nRep <- nTrial / nStim
  if (nRep != floor(nRep)) {
    stop("nTrial / nStim is not an integer, so trials cannot be reshaped by stimulus repetitions.")
  }
  
  pGo      <- array(NA_real_, dim = c(nSub, nStim, nRep))
  pCorrect <- array(NA_real_, dim = c(nSub, nStim, nRep))
  
  for (iSub in 1:nSub) {
    for (iStim in 1:nStim) {
      idx <- dList$s[iSub, ] %in% iStim
      pGo[iSub, iStim, ]      <- meanp[iSub, idx, 1] + meanp[iSub, idx, 2]
      pCorrect[iSub, iStim, ] <- meanp[iSub, idx, correctResp[iStim]]
    }
  }
  
  # Convert to 2D: subject × (stimulus × repetition)
  pGo_2d      <- matrix(aperm(pGo,      c(1, 2, 3)), nrow = nSub)
  pCorrect_2d <- matrix(aperm(pCorrect, c(1, 2, 3)), nrow = nSub)
  
  file_pgo <- file.path(out_dir, "M5_osap_samples_pGo.csv")
  file_pco <- file.path(out_dir, "M5_osap_samples_pCorrect.csv")
  
  write.csv(pGo_2d,      file = file_pgo, row.names = FALSE)
  write.csv(pCorrect_2d, file = file_pco, row.names = FALSE)
  
  cat("✅ Saved: ", file_pgo, "\n", sep = "")
  cat("✅ Saved: ", file_pco, "\n", sep = "")
  
  invisible(list(
    pgo = file_pgo,
    pco = file_pco
  ))
}

# =========================================================
# 2) Export M5_posteriors_wide.csv
# =========================================================

cat("Loading M5 HC group-level posteriors...\n")
df_hc <- load_group_X(path_m5_hc_fit) %>%
  transform_m5_group() %>%
  mutate(Group = "HC")

cat("Loading M5 OCD group-level posteriors...\n")
df_ocd <- load_group_X(path_m5_ocd_fit) %>%
  transform_m5_group() %>%
  mutate(Group = "OCD")

df_all <- bind_rows(df_hc, df_ocd)

out_m5_csv <- file.path(out_plot_dir, "M5_posteriors_wide.csv")
write.csv(df_all, out_m5_csv, row.names = FALSE)
cat("✅ Saved: ", out_m5_csv, "\n", sep = "")

# =========================================================
# 3) Export posterior predictive outputs for both HC and OCD using M5
# =========================================================

hc_out <- export_M5_predictions(
  fit_file   = hc_fit_file,
  obs_file   = hc_obs_file,
  out_dir    = hc_dir,
  group_name = "HC"
)

ocd_out <- export_M5_predictions(
  fit_file   = ocd_fit_file,
  obs_file   = ocd_obs_file,
  out_dir    = ocd_dir,
  group_name = "OCD"
)

# =========================================================
# 4) Final message
# =========================================================

cat("\n==============================\n")
cat("Done. A total of 5 core files were generated:\n")
cat("1) ", out_m5_csv, "\n", sep = "")
cat("2) ", hc_out$pgo, "\n", sep = "")
cat("3) ", hc_out$pco, "\n", sep = "")
cat("4) ", ocd_out$pgo, "\n", sep = "")
cat("5) ", ocd_out$pco, "\n", sep = "")
cat("==============================\n")