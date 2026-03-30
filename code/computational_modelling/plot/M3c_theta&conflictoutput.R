# ============================================================================
# Rebuild EEGestimates4stan.Rdata from dList
# Also export: thetapower.csv, conflicttrials.csv
# Additional step: generate conflict_raw.csv / conflict_z.csv from the M3c fit
# Adapted for: 37 subjects, 8 stimuli, 480 trials, conflict stimuli {3,4,5,6}
# ============================================================================

rm(list=ls()); graphics.off(); cat("\014")

suppressPackageStartupMessages({
  library(R.matlab)     # For fallback loading of tPow from .mat
  library(boot)         # inv.logit
  library(matrixStats)  # colMedians
})

# ---------------- Paths ----------------
rootdir    <- "F:/PIT/model/EEGmodelHC/HC"
datadir    <- file.path(rootdir, "results")
resultsdir <- file.path(datadir, "stan")
dir.create(datadir,    showWarnings = FALSE, recursive = TRUE)
dir.create(resultsdir, showWarnings = FALSE, recursive = TRUE)
setwd(datadir)

# ---------------- 1) Load dList ----------------
obs_candidates <- c(
  file.path(datadir, "obs4stanhc.Rdata"),
  file.path(datadir, "obs4stan.Rdata"),
  file.path(datadir, "EEG4stan.Rdata")
)
obsfile <- obs_candidates[file.exists(obs_candidates)][1]
if (length(obsfile) == 0) stop("Cannot find obs4stan*.Rdata / EEG4stan.Rdata under ", datadir, ".")
message("Load dList from: ", obsfile)
load(obsfile)  # -> dList
stopifnot(exists("dList"))

# Basic dimensions
stopifnot(is.matrix(dList$s))
nSub   <- nrow(dList$s)
nTrial <- ncol(dList$s)
nStim  <- if (!is.null(dList$nStim)) dList$nStim else 8L
message(sprintf("Detected nSub=%d, nTrial=%d, nStim=%d", nSub, nTrial, nStim))

# ---------------- 2) Extract / fill EEG covariates from dList ----------------
# Prefer tPow from dList; if missing, fall back to .mat
tPow <- NULL
if (!is.null(dList$tPow)) {
  tPow <- as.matrix(dList$tPow)
  message("tPow found in dList.")
} else {
  # Fallback: try reading from .mat in the same directory
  mf_candidates <- c("MFpower4Rhc.mat", "MFpower4hc.mat", "MFpower.mat", "thetapower.mat")
  mf_candidates <- file.path(datadir, mf_candidates)
  mf_file <- mf_candidates[file.exists(mf_candidates)][1]
  if (length(mf_file) == 0) stop("tPow is missing in dList, and no usable *.mat file was found (e.g., MFpower4Rhc.mat).")
  message("tPow not in dList, fallback to MAT: ", mf_file)
  m <- readMat(mf_file)
  cand <- intersect(c("tPow","MFpower","theta","thetaPower","MF_power"), names(m))
  if (length(cand) == 0) stop("Cannot find tPow/MFpower/theta variable in ", basename(mf_file), ".")
  tPow <- as.matrix(m[[cand[1]]])
}
stopifnot(all(dim(tPow) == c(nSub, nTrial)))

# Other covariates: collect if present; otherwise skip
lICPC    <- if (!is.null(dList$lICPC))   as.matrix(dList$lICPC)   else NULL
rICPC    <- if (!is.null(dList$rICPC))   as.matrix(dList$rICPC)   else NULL
pfcICPC  <- if (!is.null(dList$pfcICPC)) as.matrix(dList$pfcICPC) else NULL

# ---------------- 3) Build trial-wise motivconflict (nSub × nTrial) ----------------
# Rule: conflict stimuli = {3,4,5,6} (Go×Avoid and NoGo×Win)
winStim      <- c(1L,2L,5L,6L)
avoidStim    <- c(3L,4L,7L,8L)
conflictStim <- c(3L,4L,5L,6L)

motivconflict <- NULL
if (!is.null(dList$motivconflict)) {
  if (is.matrix(dList$motivconflict) &&
      all(dim(dList$motivconflict) == c(nSub, nTrial))) {
    motivconflict <- dList$motivconflict
    message("Use trial-wise motivconflict from dList (matrix).")
  } else if (is.vector(dList$motivconflict) &&
             length(dList$motivconflict) %in% c(nStim, 2L)) {
    message("Map stimulus-level motivconflict (vector) to trial-wise matrix.")
    mc_vec <- as.integer(dList$motivconflict)  # may be {1=conflict, 2=non-conflict} or {0/1}
    mc_vec[mc_vec == 2L] <- 0L                 # unify to 0/1
    motivconflict <- matrix(0L, nrow = nSub, ncol = nTrial)
    for (i in 1:nSub) motivconflict[i, ] <- mc_vec[dList$s[i, ]]
  }
}
if (is.null(motivconflict)) {
  message("Build trial-wise motivconflict by rule {3,4,5,6}.")
  motivconflict <- matrix(as.integer(dList$s %in% conflictStim), nrow = nSub, ncol = nTrial)
}
stopifnot(all(dim(motivconflict) == c(nSub, nTrial)))

# ---------------- 4) Save EEGestimates4stan.Rdata ----------------
eeg_rdata_out <- file.path(datadir, "EEGestimates4stan.Rdata")

vars_to_save <- c("tPow", "motivconflict")
if (exists("lICPC")   && !is.null(lICPC))   vars_to_save <- c(vars_to_save, "lICPC")
if (exists("rICPC")   && !is.null(rICPC))   vars_to_save <- c(vars_to_save, "rICPC")
if (exists("pfcICPC") && !is.null(pfcICPC)) vars_to_save <- c(vars_to_save, "pfcICPC")

save(list = vars_to_save, file = eeg_rdata_out)
message("Saved: ", eeg_rdata_out,
        "\n   - ", paste(vars_to_save, collapse = ", "),
        "\n   dimensions: tPow = ", paste(dim(tPow), collapse = "×"),
        if (!is.null(lICPC))   paste0("; lICPC = ",   paste(dim(lICPC), collapse="×")) else "",
        if (!is.null(rICPC))   paste0("; rICPC = ",   paste(dim(rICPC), collapse="×")) else "",
        if (!is.null(pfcICPC)) paste0("; pfcICPC = ", paste(dim(pfcICPC), collapse="×")) else "",
        "\n   motivconflict = ", paste(dim(motivconflict), collapse = "×"),
        " (1=conflict, 0=non-conflict)")

# ---------------- 5) Export thetapower.csv / conflicttrials.csv ----------------
write.csv(tPow,               file = file.path(resultsdir, "thetapower.csv"),     row.names = FALSE)
write.csv(motivconflict,      file = file.path(resultsdir, "conflicttrials.csv"), row.names = FALSE)
message("Also wrote CSVs to: ", resultsdir,
        "\n - thetapower.csv",
        "\n - conflicttrials.csv")

# ---------------- 6) Additional step: generate conflict_raw.csv / conflict_z.csv from M3c ----------------
# Required: M3c fit file (commonly named M9fits_obsData.Rdata; old naming M5fits_obsData.Rdata)
fit_candidates <- c(file.path(resultsdir, "M9fits_obsData.Rdata"),
                    file.path(resultsdir, "M5fits_obsData.Rdata"))
fitfile <- fit_candidates[file.exists(fit_candidates)][1]
if (length(fitfile) == 0) {
  stop("M3c fit file not found: ", paste(fit_candidates, collapse=" or "))
}
message("Load M3c fit: ", fitfile)
load(fitfile)  # -> f2 (stanfit)

draws <- as.matrix(f2)
pick_cols <- function(nm){
  cols <- grep(paste0("^", nm, "\\["), colnames(draws))
  if (length(cols)==0) stop("Cannot find parameter in posterior draws: ", nm)
  draws[, cols, drop=FALSE]
}
# M3c parameters: x1..x5 = rho, epsilon, go_bias, pi, kappa
x1 <- pick_cols("x1"); x2 <- pick_cols("x2"); x3 <- pick_cols("x3"); x4 <- pick_cols("x4"); x5 <- pick_cols("x5")
stopifnot(ncol(x1) == nSub)
nPerm <- nrow(x1)

# Compute modeled conflict trial-wise using the paper definition
# and take the posterior median
M3c_conflict <- function(fbsens, lr, gobias, pibias, kappa,
                         stimSeq, respSeq, fbSeq,
                         nTrial, nStim, nPerm,
                         winStim = c(1,2,5,6), avoidStim = c(3,4,7,8)) {
  V <- numeric(nStim); V[winStim] <- +0.5; V[avoidStim] <- -0.5
  confl <- matrix(NA_real_, nrow=nPerm, ncol=nTrial)
  for (i in 1:nPerm) {
    rho     <- exp(fbsens[i])
    epsilon <- boot::inv.logit(lr[i])
    if (epsilon < .5) {
      biasEps <- boot::inv.logit(c(NA, lr[i]-kappa[i])); biasEps[1] <- 2*epsilon - biasEps[2]
    } else {
      biasEps <- boot::inv.logit(c(lr[i]+kappa[i], NA));  biasEps[2] <- 2*epsilon - biasEps[1]
    }
    # Initialize Q: 3(left/right/NoGo) × nStim
    Q <- matrix(rep(rho * c(+.5, +.5, -.5, -.5, +.5, +.5, -.5, -.5), each=3),
                nrow=3, ncol=nStim, byrow=FALSE)
    valenced <- rep(0, nStim)
    for (t in 1:nTrial) {
      s <- stimSeq[t]
      # Conflict definition (same as in the paper)
      confl[i,t] <- - valenced[s] * ( V[s] * ( mean(Q[1:2,s]) - Q[3,s] ) )
      # Update Q using observed response and outcome
      resp <- respSeq[t]; out <- fbSeq[t]; ie <- ifelse(resp %in% c(1,2), 1, 2)
      if ((ie==1 && out== 1) || (ie==2 && out==-1)) {
        Q[resp,s] <- Q[resp,s] + biasEps[ie] * (rho*out - Q[resp,s])
      } else {
        Q[resp,s] <- Q[resp,s] + epsilon     * (rho*out - Q[resp,s])
      }
      if (out != 0) valenced[s] <- 1
    }
  }
  matrixStats::colMedians(confl, na.rm=TRUE)
}

# Compute trial-wise conflict for each subject
conflict_raw <- matrix(NA_real_, nrow=nSub, ncol=nTrial)
for (iSub in 1:nSub) {
  message("Compute conflict: subject ", iSub, "/", nSub)
  conflict_raw[iSub, ] <- M3c_conflict(
    fbsens = x1[,iSub], lr = x2[,iSub], gobias = x3[,iSub],
    pibias = x4[,iSub], kappa = x5[,iSub],
    stimSeq = dList$s[iSub,], respSeq = dList$ya[iSub,], fbSeq = dList$r[iSub,],
    nTrial = nTrial, nStim = nStim, nPerm = nPerm,
    winStim = winStim, avoidStim = avoidStim
  )
}

# If trial2use==0 exists, set to NA before z-scoring
if (!is.null(dList$trial2use)) conflict_raw[dList$trial2use == 0] <- NA_real_

# Within-subject z-score (row-wise z)
conflict_z <- t(scale(t(conflict_raw)))  # [nSub × nTrial]

# Export
write.csv(conflict_raw, file = file.path(resultsdir, "conflict_raw.csv"), row.names = FALSE)
write.csv(conflict_z,   file = file.path(resultsdir, "conflict_z.csv"),   row.names = FALSE)
message("Also wrote: ",
        "\n - ", file.path(resultsdir, "conflict_raw.csv"),
        "\n - ", file.path(resultsdir, "conflict_z.csv"))






####### Generate M5_osap_conflict.csv from the model ####################
## ================== M5: posterior predictive conflict to CSV ==================
rm(list = ls()); graphics.off(); cat("\014")

suppressPackageStartupMessages({
  library(boot)          # inv.logit
  library(matrixStats)   # colMedians
  library(rstan)
})

## -------- Paths (input=H drive; output=F drive) -----------------------------------------
input_fit_dir  <- "H:/mt_model/HC"                 # Your M5 fit directory
stanfitfile    <- file.path(input_fit_dir, "M5fits_obsData.Rdata")

output_dir     <- "F:/PIT/model/EEGmodelHC/HC/results"      # Target output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## -------- Load fit (must contain f2 and dList) ------------------------------------
stopifnot(file.exists(stanfitfile))
load(stanfitfile)   # -> f2, dList

# Basic dimensions (use dList as the reference)
if (is.null(dList$nSub) || is.null(dList$nTrial)) {
  dList$nSub   <- nrow(dList$s)
  dList$nTrial <- ncol(dList$s)
}
nSub   <- as.integer(dList$nSub)
nStim  <- as.integer(dList$nStim)
nResp  <- as.integer(dList$nResp)
nTrial <- as.integer(dList$nTrial)
stopifnot(nResp == 3L)

cat(sprintf("Loaded fit: nSub=%d, nStim=%d, nTrial=%d\n", nSub, nStim, nTrial))

## ================== Your specified column-slice extraction block (kept unchanged) ==================
estimates <- as.matrix(f2)[ , 11:(10 + 5L * dList$nSub) ]
nPerm <- 4000L
rho     <- estimates[ , (0*dList$nSub+1):(1*dList$nSub), drop=FALSE]
epsilon <- estimates[ , (1*dList$nSub+1):(2*dList$nSub), drop=FALSE]
gobias  <- estimates[ , (2*dList$nSub+1):(3*dList$nSub), drop=FALSE]
pibias  <- estimates[ , (3*dList$nSub+1):(4*dList$nSub), drop=FALSE]
kappa   <- estimates[ , (4*dList$nSub+1):(5*dList$nSub), drop=FALSE]
## ============================================================================

# For robustness: reduce nPerm automatically if the actual number of draws is smaller
nPerm <- min(nPerm, nrow(estimates))

## -------- M5 simulation function (same as original author; median across posterior draws) -------------------
M5_sim <- function(fbsens, lr, gb, pi, kap, stimSeq, respSeq, fbSeq,
                   nTrial, nResp, nStim, nPerm){
  # Use the smaller value between requested nPerm and available draw rows
  nPerm <- min(nPerm, nrow(as.matrix(fbsens)))
  
  allp  <- array(NA_real_, dim = c(nPerm, nTrial, nResp))
  confl <- array(NA_real_, dim = c(nPerm, nTrial))
  
  for(iIter in 1:nPerm){
    rho_i <- exp(fbsens[iIter])       # feedback sensitivity
    eps_i <- inv.logit(lr[iIter])     # baseline learning rate
    
    if (eps_i < .5){
      bias_lr <- inv.logit(c(NA, lr[iIter]-kap[iIter]))  # (go, nogo)
      bias_lr[1] <- 2*eps_i - bias_lr[2]
    } else {
      bias_lr <- inv.logit(c(lr[iIter]+kap[iIter], NA))
      bias_lr[2] <- 2*eps_i - bias_lr[1]
    }
    
    # Initial Q/V (same as original)
    Q <- matrix(rho_i * c(.5,.5,-.5,-.5), nrow=nResp, ncol=nStim, byrow=TRUE)
    V <- rep(c(.5,.5,-.5,-.5), nStim/4)
    valenced <- rep(0, nStim)
    q <- matrix(0, nrow=nTrial, ncol=nResp)
    
    for (t in 1:nTrial){
      s <- stimSeq[t]
      
      # Action weights for the current trial
      q[t,] <- valenced[s] * Q[,s]
      q[t,1:2] <- q[t,1:2] + gb[iIter] + pi[iIter] * valenced[s] * V[s]
      
      # Unstandardized Pavlovian–instrumental conflict
      # (same formula as in the original author code)
      confl[iIter,t] <- -valenced[s] * ( V[s] * ( mean(Q[1:2,s]) - Q[3,s] ) )
      
      # Softmax probabilities
      # (kept for full consistency with the original output)
      p <- exp(q[t,]) / sum(exp(q[t,])); p[is.nan(p)] <- 1
      allp[iIter,t,] <- p
      
      # Q update
      resp <- respSeq[t]; out <- fbSeq[t]
      iEps <- ifelse(resp %in% c(1,2), 1, 2)  # 1=Go, 2=NoGo
      if ((iEps==1 && out== 1) || (iEps==2 && out==-1)){
        Q[resp,s] <- Q[resp,s] + bias_lr[iEps] * (rho_i*out - Q[resp,s])
      } else {
        Q[resp,s] <- Q[resp,s] + eps_i * (rho_i*out - Q[resp,s])
      }
      if (out != 0) valenced[s] <- 1
    }
  }
  
  list(
    meanp = colMeans(allp, na.rm = TRUE),                 # nTrial × nResp
    confl = matrixStats::colMedians(confl, na.rm = TRUE)  # median across posterior draws
  )
}

## -------- Simulate by subject and collect conflict medians --------------------------------------
meanp <- array(NA_real_, dim = c(nSub, nTrial, nResp))
confl <- array(NA_real_, dim = c(nSub, nTrial))

for (iSub in 1:nSub){
  cat("Simulating subject", iSub, "/", nSub, "...\n")
  out <- M5_sim(
    fbsens = rho[1:nPerm, iSub, drop=FALSE],
    lr     = epsilon[1:nPerm, iSub, drop=FALSE],
    gb     = gobias[1:nPerm, iSub, drop=FALSE],
    pi     = pibias[1:nPerm, iSub, drop=FALSE],
    kap    = kappa[1:nPerm, iSub, drop=FALSE],
    stimSeq= dList$s[iSub, ],
    respSeq= dList$ya[iSub, ],
    fbSeq  = dList$r[iSub, ],
    nTrial = nTrial, nResp = nResp, nStim = nStim, nPerm = nPerm
  )
  meanp[iSub,,] <- out$meanp
  confl[iSub, ] <- out$confl
}

## -------- Write the unstandardized conflict matrix with posterior medians -----------------------
out_csv <- file.path(output_dir, "M5_osap_conflict.csv")
write.csv(confl, file = out_csv, row.names = FALSE, col.names = FALSE)

cat("\nOutput written: ", out_csv,
    "\n   shape = ", nSub, " × ", nTrial,
    "\n   meaning = modeled conflict strength for each subject × trial",
    " (unstandardized; posterior median across draws).\n", sep = "")