## 0) Libs & Stan options ------------------------------------------------------
suppressPackageStartupMessages({
  require(rstan)
  require(data.table); require(parallel); require(boot)
  require(R.matlab);   require(devtools); require(corrgram)
  require(ggplot2);    require(psych)
})
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## 1) Paths --------------------------------------------------------------------
rootdir    <- "F:/PIT/model/EEGmodelHC/HC"
scriptdir  <- file.path(rootdir, "scripts")
datadir    <- file.path(rootdir, "results")
plotdir    <- file.path(rootdir, "plots")
resultsdir <- file.path(datadir, "stan")

dir.create(plotdir,    recursive = TRUE, showWarnings = FALSE)
dir.create(resultsdir, recursive = TRUE, showWarnings = FALSE)

## 2) Data files ---------------------------------------------------------------
obs4stanfile <- file.path(datadir, "obs4stanhc.Rdata")
raw_matfile  <- file.path(datadir, "MFpower4Rhc.mat")

## 3) Load scripts -------------------------------------------------------------
source(file.path(scriptdir, "models4stan.R"))
source(file.path(scriptdir, "fitmodel.R"))
source(file.path(scriptdir, "extractEstimates.R"))
source(file.path(scriptdir, "plotFits.R"))
source(file.path(scriptdir, "getWAIC.R"))
source(file.path(scriptdir, "waic.R"))
# Enable this only if you need to rebuild RData from the .mat file:
# source(file.path(scriptdir, "prepEEG4stan.R"))

## 4) Prepare obs data once (only if missing) ---------------------------------
if (!file.exists(obs4stanfile)) {
  if (!file.exists(raw_matfile))
    stop("Missing raw .mat file: ", raw_matfile, "; cannot create obs4stanhc.Rdata")
  message("Preparing obs4stanhc.Rdata from MFpower4Rhc.mat ...")
  prepEEG4stan(obs4stanfile, raw_matfile)
}

## 5) Quick sanity check on dList ----------------------------------------------
load(obs4stanfile)  # -> dList
stopifnot(exists("dList"))
stopifnot(all(c("ya","nSub","nTrial","nResp","nStim") %in% names(dList)))
stopifnot(nrow(dList$ya) == dList$nSub, ncol(dList$ya) == dList$nTrial)
if (dList$nResp != 3L)  stop("nResp is not 3; it does not match your model family")
if (dList$nStim != 8L)  stop("nStim is not 8; please check how dList was constructed")

## 6) motivconflict (used by EEG models) ---------------------------------------
if (!("motivconflict" %in% names(dList))) {
  motivconflict <- as.integer(c(2,2,1,1,1,1,2,2))  # 1 = conflict, 2 = non-conflict
  assign("motivconflict", motivconflict, envir = .GlobalEnv)  # fitmodel will read this
}

## 7) K (number of parameters) per model ---------------------------------------
# Based on your current setup, elements 21 and 22 are added
# (M21/M22 each has 1 EEG weight -> K = 6)
Kvec <- c(2,3,4,4,5,2,2,4,5,6,6,6,6,6,6,6,7,7,6,6,6,6,7,7)
# Index:   1 2 3 4 5 6 7 8  9.............................24
stopifnot(length(Kvec) >= 24)

## 8) Define initQV (original rule; placed outside the loop for fitmodel) ------
initQV <- function(nResp, nStim) {
  # Same as the original rule:
  # initial Q values set Win/Go to +.5 and Avoid/NoGo to -.5 (repeating pattern)
  Qi <- matrix(c(.5, .5, -.5, -.5), nrow = nResp, ncol = nStim, byrow = TRUE)
  Vi <- rep(c(.5, .5, -.5, -.5), nStim / 4)
  list(Qi = Qi, Vi = Vi)
}
# Original special case for M7
# (not used in the current 9..14 run, but harmless to keep)
initQV_M7 <- function(nResp, nStim) {
  Qi <- matrix(0, nrow = nResp, ncol = nStim, byrow = TRUE)
  Vi <- rep(c(.5, .5, -.5, -.5), nStim / 4)
  list(Qi = Qi, Vi = Vi)
}

## 9) Helper: safe execution with success flag ---------------------------------
safe_do <- function(expr, on_error_msg) {
  ok <- tryCatch({ force(expr); TRUE }, error = function(e) {
    message(on_error_msg, ": ", conditionMessage(e)); FALSE
  })
  ok
}

## 10) Fit loop (models 9..14) -------------------------------------------------
for (model2fit in 14) {
  
  ##for (model2fit in 1:24) {
  message("\n================= MODEL ", model2fit, " =================")
  K <- Kvec[model2fit]
  
  stanfitsfile <- file.path(resultsdir, paste0("M", model2fit, "fits_obsData.Rdata"))
  estimatefile <- file.path(resultsdir, paste0("M", model2fit, "estimates_obsData.Rdata"))
  waicfile     <- file.path(resultsdir, paste0("M", model2fit, "waic_obsData.Rdata"))
  plotfitsfile <- file.path(plotdir,    paste0("M", model2fit, "fits_obsData"))
  
  ## -- Fit --
  if (!file.exists(stanfitsfile)) {
    message("Fitting model ", model2fit, " ...")
    ok_fit <- safe_do(
      fitmodel(model2fit,
               dat4stanfile = obs4stanfile,
               stanfitsfile = stanfitsfile,
               calcLogLik   = TRUE),
      on_error_msg = paste("!! Fitting failed for model", model2fit)
    )
    if (!ok_fit) next
  } else {
    message("Fit exists: ", stanfitsfile)
  }
  
  ## -- Extract parameters --
  if (!file.exists(estimatefile)) {
    message("Extracting estimates for model ", model2fit, " ...")
    ok_est <- safe_do(
      { library(rstan); extractEstimates(stanfitsfile, estimatefile, model2fit) },
      on_error_msg = paste("!! extractEstimates failed for model", model2fit)
    )
    if (!ok_est) next
  } else {
    message("Estimates exist: ", estimatefile)
  }
  
  ## -- Plot fits --
  message("Plotting fits for model ", model2fit, " ...")
  ok_plot <- safe_do(
    plotFits(estimatefile, plotfitsfile, model2fit),
    on_error_msg = paste("!! plotFits failed for model", model2fit)
  )
  if (!ok_plot) next
  
  ## -- Calculate WAIC --
  message("Calculating WAIC for model ", model2fit, " ...")
  ok_waic <- safe_do(
    { library(rstan); getWAIC(stanfitsfile, waicfile, model2fit) },
    on_error_msg = paste("!! getWAIC failed for model", model2fit)
  )
  if (!ok_waic) next
  
  message("END OF MODEL ", model2fit)
}