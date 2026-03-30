rm(list = ls())
options(stringsAsFactors = FALSE)

# ================================
# Standalone Script 1: Export WAIC to CSV
# ================================

# ---------- 1) Settings ----------
# Set the directory containing the WAIC files
base_dir <- "F:/PIT/model/EEGmodelHC/OCD/results/stan"

# Set the output directory for CSV files
out_dir  <- "F:/PIT/model/EEGmodelHC/OCD/results/stan/plot"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Define the model sets to export
# Modify these as needed
sets <- list(
  Set0 = c(1, 2, 3, 4, 5),
  Set1 = c(5, 10, 12, 13, 14),
  Set2 = c(5, 10, 15, 16),
  Set3 = c(5, 12, 19, 21)
)

# Output filename prefix
group_tag <- "OCD"   # e.g., "HC" or "OCD"

# ---------- 2) Utility functions ----------
get_model <- function(p) {
  sub("^.*(M[0-9]+).*?$", "\\1", p)
}

load_waic <- function(path) {
  stopifnot(file.exists(path))
  e <- new.env()
  load(path, envir = e)
  objs <- mget(ls(e, all.names = TRUE), envir = e, inherits = FALSE)
  
  extract_waic <- function(x) {
    # 1) Single numeric value
    if (is.numeric(x) && length(x) == 1L) return(as.numeric(x))
    
    # 2) matrix / data.frame
    if (is.matrix(x) || is.data.frame(x)) {
      rn <- rownames(x)
      cn <- colnames(x)
      
      if (!is.null(rn) && "waic" %in% rn) {
        if (!is.null(cn) && "Estimate" %in% cn) return(as.numeric(x["waic", "Estimate"]))
        if (ncol(x) >= 1L) return(as.numeric(x["waic", 1]))
      }
      
      if (!is.null(rn) && "elpd_waic" %in% rn) {
        val <- if (!is.null(cn) && "Estimate" %in% cn) x["elpd_waic", "Estimate"] else x["elpd_waic", 1]
        return(as.numeric(-2 * val))
      }
      
      if (!is.null(cn) && "waic" %in% cn && nrow(x) >= 1L) {
        return(as.numeric(x[1, "waic"]))
      }
      
      if (!is.null(cn) && "elpd_waic" %in% cn && nrow(x) >= 1L) {
        return(as.numeric(-2 * x[1, "elpd_waic"]))
      }
    }
    
    # 3) list
    if (is.list(x)) {
      if (!is.null(x$waic) && is.numeric(x$waic) && length(x$waic) == 1L) {
        return(as.numeric(x$waic))
      }
      if (!is.null(x$WAIC) && is.numeric(x$WAIC) && length(x$WAIC) == 1L) {
        return(as.numeric(x$WAIC))
      }
      if (!is.null(x$elpd_waic) && is.numeric(x$elpd_waic) && length(x$elpd_waic) == 1L) {
        return(as.numeric(-2 * x$elpd_waic))
      }
      if (!is.null(x$estimates)) {
        val <- try(extract_waic(x$estimates), silent = TRUE)
        if (!inherits(val, "try-error") && !is.null(val)) return(as.numeric(val))
      }
      for (nm in names(x)) {
        val <- try(extract_waic(x[[nm]]), silent = TRUE)
        if (!inherits(val, "try-error") && !is.null(val)) return(as.numeric(val))
      }
    }
    
    NULL
  }
  
  for (obj in objs) {
    val <- extract_waic(obj)
    if (!is.null(val)) return(as.numeric(val))
  }
  
  stop("No parsable WAIC was found in ", basename(path), ".")
}

save_csv_for_matlab <- function(model_ids, filename, base_dir) {
  paths <- file.path(base_dir, paste0("M", model_ids, "waic_obsData.Rdata"))
  
  miss <- paths[!file.exists(paths)]
  if (length(miss) > 0) {
    stop("The following WAIC files do not exist:\n", paste(miss, collapse = "\n"))
  }
  
  df <- data.frame(
    model = vapply(paths, get_model, character(1)),
    WAIC  = vapply(paths, load_waic, numeric(1)),
    stringsAsFactors = FALSE
  )
  
  write.csv(df, filename, row.names = FALSE)
  message("Saved: ", filename)
  invisible(df)
}

# ---------- 3) Batch export ----------
for (nm in names(sets)) {
  out_file <- file.path(out_dir, paste0(nm, "_data", group_tag, ".csv"))
  save_csv_for_matlab(
    model_ids = sets[[nm]],
    filename  = out_file,
    base_dir  = base_dir
  )
}

cat("\n✅ WAIC CSV export completed.\n")