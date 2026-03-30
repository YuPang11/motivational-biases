rm(list = ls())
suppressPackageStartupMessages(library(rstan))

# ==========================================
# Standalone Script 2: Export raw Stan f2 draws to CSV
# ==========================================

# ---------- 1) Configuration ----------
# Define the export job list here
# Each row represents one export job: input Rdata -> output CSV
jobs <- data.frame(
  fit_rda = c(
    "F:/PIT/model/EEGmodelHC/HC/results/stan/M12fits_obsData.Rdata",
    "F:/PIT/model/EEGmodelHC/HC/results/stan/M16fits_obsData.Rdata",
    "F:/PIT/model/EEGmodelHC/HC/results/stan/M19fits_obsData.Rdata",
    "F:/PIT/model/EEGmodelHC/OCD/results/stan/M12fits_obsData.Rdata",
    "F:/PIT/model/EEGmodelHC/OCD/results/stan/M16fits_obsData.Rdata",
    "F:/PIT/model/EEGmodelHC/OCD/results/stan/M21fits_obsData.Rdata"
  ),
  out_csv = c(
    "F:/PIT/ISPSplot/M12_f2_drawshc.csv",
    "F:/PIT/ISPSplot/M16_f2_drawshc.csv",
    "F:/PIT/ISPSplot/M19_f2_drawshc.csv",
    "F:/PIT/ISPSplot/M12_f2_drawsocd.csv",
    "F:/PIT/ISPSplot/M16_f2_drawsocd.csv",
    "F:/PIT/ISPSplot/M21_f2_drawsocd.csv"
  ),
  stringsAsFactors = FALSE
)

# ---------- 2) Export function ----------
export_f2_to_csv <- function(fit_rda, out_csv) {
  if (!file.exists(fit_rda)) {
    stop("File not found: ", fit_rda)
  }
  
  out_dir <- dirname(out_csv)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  e <- new.env()
  load(fit_rda, envir = e)
  
  if (!exists("f2", envir = e)) {
    stop("Object f2 was not found in ", basename(fit_rda), ".")
  }
  
  f2mat <- as.matrix(e$f2)
  write.csv(f2mat, out_csv, row.names = FALSE)
  
  message("Saved: ", out_csv,
          " | draws = ", nrow(f2mat),
          " | params = ", ncol(f2mat))
}

# ---------- 3) Run all jobs ----------
for (i in seq_len(nrow(jobs))) {
  cat("\n[", i, "/", nrow(jobs), "] Exporting:\n", sep = "")
  cat("  input : ", jobs$fit_rda[i], "\n", sep = "")
  cat("  output: ", jobs$out_csv[i], "\n", sep = "")
  export_f2_to_csv(jobs$fit_rda[i], jobs$out_csv[i])
}

cat("\n✅ Raw Stan draw CSV export completed.\n")