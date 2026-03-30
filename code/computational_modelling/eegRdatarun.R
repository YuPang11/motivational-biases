########## EEG model fitting data preparation (R) ##########

# 0) Load function ------------------------------------------------------
root <- "F:/PIT/model/EEGmodelHC"     # unified root path
src_file <- file.path(root, "prepEEG4stan.R")
stopifnot(file.exists(src_file))
source(src_file)
stopifnot(exists("prepEEG4stan"))

# 1) Paths --------------------------------------------------------------
rawfile      <- file.path(root, "MFpower4Rocd.mat")     # .mat exported from MATLAB
obs4stanfile <- file.path(root, "obs4stanocd.Rdata")    # output .Rdata for Stan
stopifnot(file.exists(rawfile))
dir.create(dirname(obs4stanfile), showWarnings = FALSE, recursive = TRUE)

# 2) Generate RData for Stan -------------------------------------------
prepEEG4stan(obs4stanfile, rawfile)

# 3) Quick output check -------------------------------------------------
load(obs4stanfile)   # loads dList
stopifnot(exists("dList"))

# 3.1 Shape and basic info
cat("\n== Basic shape ==\n")
print(dim(dList$ya))          # expected 37 x 480
print(dim(dList$tPow))        # expected 37 x 480
cat("unique s: ", sort(unique(as.vector(dList$s))), "\n")          # expected 1..8
cat("unique rew: ", sort(unique(as.vector(dList$rew))), "\n")      # expected {1,2}
cat("nTrial =", dList$nTrial, " nSub =", dList$nSub,
    " nResp =", dList$nResp, " nStim =", dList$nStim, "\n")

# 3.2 Type check (all should be integer)
for (nm in c("ya","a","s","r","rew","yax","nTrial","nSub","nResp","nStim","nData")) {
  stopifnot(is.integer(dList[[nm]]))
}
stopifnot(length(dList$motivconflict) == dList$nStim)

# 3.3 Value range check
stopifnot(all(dList$ya  %in% 1:3))     # 1/2/3 = left/right/NoGo
stopifnot(all(dList$a   %in% 1:2))     # Go/NoGo binary
stopifnot(all(dList$rew %in% 1:2))     # Win=1, Avoid=2
stopifnot(all(sort(unique(as.vector(dList$s))) %in% 1:dList$nStim))

# 3.4 Check zero-masking for trial2use (if available)
if (!is.null(dList$trial2use)) {
  mask0 <- dList$trial2use == 0
  cat("Number of entries with trial2use==0:", sum(mask0), "\n")
  if (!is.null(dList$tPow))     stopifnot(all(dList$tPow    [mask0] == 0))
  if (!is.null(dList$pfcICPC))  stopifnot(all(dList$pfcICPC [mask0] == 0))
  if (!is.null(dList$lICPC))    stopifnot(all(dList$lICPC   [mask0] == 0))
  if (!is.null(dList$rICPC))    stopifnot(all(dList$rICPC   [mask0] == 0))
}

cat("\n✅ Check passed: obs4stanocd.Rdata is ready.\n")





# Data check
## Show all variables for subID == 201 (trial by trial)
file <- "F:/PIT/model/EEGmodelHC/obs4stanocd.Rdata" 
stopifnot(file.exists(file))
load(file)                                  # -> dList
stopifnot(exists("dList"))

# Find the row index for subID == 201
if (!("subID" %in% names(dList))) stop("subID is not available in dList; cannot match by subject ID.")
idx <- which(as.integer(dList$subID) == 201) 
if (length(idx) == 0) stop("subID == 201 was not found. Available IDs: ", paste(dList$subID, collapse=", "))

# Helper: extract the full row for this subject from a matrix
# return an NA vector if the variable is missing
pick_row <- function(M) {
  if (is.null(M)) return(rep(NA, dList$nTrial))
  as.vector(M[idx, ])
}

# Conflict label (if motivconflict exists, map from s to 1=conflict / 2=non-conflict)
conflict_label <- if (!is.null(dList$motivconflict)) {
  as.integer(dList$motivconflict[ as.vector(dList$s[idx, ]) ])
} else {
  rep(NA_integer_, dList$nTrial)
}

# Build a data.frame
# variable names kept consistent with the current dataset
df201 <- data.frame(
  sub        = dList$subID[idx],
  trial      = seq_len(dList$nTrial),
  s          = as.integer(pick_row(dList$s)),
  rew        = as.integer(pick_row(dList$rew)),
  r          = as.integer(pick_row(dList$r)),
  a          = as.integer(pick_row(dList$a)),
  ya         = as.integer(pick_row(dList$ya)),
  tPow       = as.numeric(pick_row(dList$tPow)),
  pfcICPC    = as.numeric(pick_row(dList$pfcICPC)),
  lICPC      = as.numeric(pick_row(dList$lICPC)),
  rICPC      = as.numeric(pick_row(dList$rICPC)),
  IVleft     = as.integer(pick_row(dList$IVleft)),
  IVright    = as.integer(pick_row(dList$IVright)),
  trial2use  = as.integer(pick_row(dList$trial2use)),
  conflict   = conflict_label,
  stringsAsFactors = FALSE
)

# Print all rows
# use View(df201) if needed
print(df201, row.names = FALSE)

# Quick preview:
# head(df201, 20)

# Show selected columns only (example: behavior + conflict):
# df201[, c("trial","s","rew","r","a","ya","conflict")]