## ============== Step 4: Extract subject-level parameters (M9/10/12/13/14) ==============
rm(list = ls())
suppressPackageStartupMessages({
  library(rstan)
  library(data.table)
})

# ---- Paths (updated to your new directory) ----
rootdir    <- "F:/PIT/model/EEGmodelHC/HC"
resultsdir <- file.path(rootdir, "results", "stan")
dir.create(resultsdir, recursive = TRUE, showWarnings = FALSE)

# ---- Function: extract subject-level posterior means from M?fits_obsData.Rdata ----
extractEstimates <- function(stanfitsfile, estimatefile) {
  # 1) Load Stan fit results (should contain f2 and dList)
  if (!file.exists(stanfitsfile)) {
    stop("File not found: ", stanfitsfile)
  }
  load(stanfitsfile)  # -> f2, dList, mstr (based on your save format)
  if (!exists("f2"))    stop("Object f2 does not exist: ", stanfitsfile)
  if (!exists("dList")) stop("Object dList does not exist: ", stanfitsfile)
  
  nSub <- as.integer(dList$nSub)
  K    <- as.integer(dList$K)
  if (!is.finite(nSub) || !is.finite(K) || K < 1L) stop("Invalid dList$nSub or dList$K")
  
  # 2) Subject level: posterior means of x1..xK (one value per subject)
  # First try extract(); if it fails, fall back to as.matrix() + regex
  xbay <- data.table(sub = seq_len(nSub))
  if ("ID" %in% names(dList)) xbay[, ID := dList$ID]
  
  # Backup objects (created only once if fallback is triggered)
  fmeans <- NULL
  cols   <- NULL
  
  for (iK in seq_len(K)) {
    parname <- paste0("x", iK)
    
    # Preferred: extract() -> posterior draws × nSub -> column means by subject
    ok <- TRUE
    vals <- tryCatch({
      ext <- rstan::extract(f2, pars = parname)
      if (is.list(ext) && !is.null(ext[[parname]])) {
        colMeans(ext[[parname]])
      } else stop("extract() did not return this parameter")
    }, error = function(e) { ok <<- FALSE; NULL })
    
    # Fallback: as.matrix() + anchored regex to avoid matching x1 to x10
    if (!ok) {
      if (is.null(fmeans)) {
        mtx    <- as.matrix(f2)
        fmeans <- colMeans(mtx)
        cols   <- names(fmeans)
      }
      idx <- grep(paste0("^", parname, "\\["), cols)
      if (length(idx) == 0) stop("Parameter column not found: ", parname, " (fallback mode)")
      # Sort by subject index inside brackets
      ord <- order(as.integer(sub(".*\\[(\\d+)\\]$", "\\1", cols[idx])))
      idx <- idx[ord]
      vals <- as.numeric(fmeans[idx])
    }
    
    if (length(vals) != nSub) {
      warning(sprintf("%s returned %d subject values (expected %d)", parname, length(vals), nSub))
      # Truncate or pad to nSub to avoid save failure
      length(vals) <- nSub
    }
    xbay[[parname]] <- vals
  }
  
  # 3) Save output (same as your old script: at least save xbay and K)
  save(xbay, K, file = estimatefile)
  cat("Saved subject-level parameters to: ", estimatefile, "\n", sep = "")
}

# ---- Batch-generate estimate files (run only 9/10/12/13/14) ----
models <- c(9, 10, 12, 13, 14, 15, 16, 19, 21)
for (m in models) {
  cat("\n=== Extracting subject-level parameters for M", m, " ===\n", sep = "")
  stanfitsfile <- file.path(resultsdir, sprintf("M%dfits_obsData.Rdata", m))
  estimatefile <- file.path(resultsdir, sprintf("M%destimates_obsData.Rdata", m))
  try(extractEstimates(stanfitsfile, estimatefile), silent = FALSE)
}

# ---- Optional: quick preview of one model (e.g., M14) ----
m_demo <- 14
ef_demo <- file.path(resultsdir, sprintf("M%destimates_obsData.Rdata", m_demo))
if (file.exists(ef_demo)) {
  load(ef_demo)  # -> xbay, K
  cat("\nPreview of xbay for M", m_demo, " (first 6 rows):\n", sep = "")
  print(head(xbay))
  cat("K = ", K, " (this model contains x1..xK, with ", K, " subject-level parameters)\n", sep = "")
}
## ======================================================================