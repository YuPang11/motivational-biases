rm(list = ls())

suppressPackageStartupMessages(require(rstan))

rootdir    <- "F:/PIT/model/EEGmodelHC/HC"
datadir    <- file.path(rootdir, "results")
resultsdir <- file.path(datadir, "stan")

models_to_check <- c(9, 10, 12, 13, 14, 15, 16, 19, 21) # Specify the models to inspect, e.g., c(9,10,12,14)

summ_list <- list()

for (m in models_to_check) {
  stanfitsfile <- file.path(resultsdir, paste0("M", m, "fits_obsData.Rdata"))
  if (!file.exists(stanfitsfile)) {
    message("Skip M", m, ": file not found.")
    next
  }
  load(stanfitsfile)                 # -> dList, f1, f2, mstr
  f2@.MISC <- new.env(parent = globalenv())
  
  nParam <- dList$K * (dList$nSub + 2)
  ss <- summary(f2)$summary
  Rhat_vec <- ss[1:nParam, "Rhat"]
  neff_vec <- ss[1:nParam, "n_eff"]
  
  summ_list[[paste0("M", m)]] <- data.frame(
    model   = paste0("M", m),
    K       = dList$K,
    nSub    = dList$nSub,
    nParam  = nParam,
    max_Rhat= max(Rhat_vec, na.rm = TRUE),
    p_Rhat_gt_1.1 = mean(Rhat_vec > 1.1, na.rm = TRUE),
    min_n_eff     = min(neff_vec, na.rm = TRUE),
    median_n_eff  = median(neff_vec, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

diag_table <- do.call(rbind, summ_list)
print(diag_table, row.names = FALSE)

## Quickly inspect the Rhat distribution for one model:
if (nrow(diag_table) > 0) {
  pick <- diag_table$model[1]
  m <- as.integer(sub("M", "", pick))
  load(file.path(resultsdir, paste0("M", m, "fits_obsData.Rdata")))
  f2@.MISC <- new.env(parent = globalenv())
  nParam <- dList$K * (dList$nSub + 2)
  Rhat_vec <- summary(f2)$summary[1:nParam, "Rhat"]
  plot(density(Rhat_vec, na.rm = TRUE), main = paste0(pick, " — Rhat density"))
}
