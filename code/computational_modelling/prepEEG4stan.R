# ------------------------------------------------------------
# Convert MATLAB MFpower4Rhc*.mat into RData containing
# "behavior + EEG (transformed)" for Stan
# ------------------------------------------------------------
library(R.matlab)

prepEEG4stan <- function(obs4stanfile, rawfile, zero_mask = TRUE, verbose = TRUE) {
  ## 1) Key: keep original .mat variable names (do not auto-fix names),
  ##    and add alias patches when needed
  E <- R.matlab::readMat(rawfile, fixNames = FALSE)
  if (!"l_ICPCmotor" %in% names(E) && "l.ICPCmotor" %in% names(E)) E[["l_ICPCmotor"]] <- E[["l.ICPCmotor"]]
  if (!"r_ICPCmotor" %in% names(E) && "r.ICPCmotor" %in% names(E)) E[["r_ICPCmotor"]] <- E[["r.ICPCmotor"]]
  
  # Required keys
  req <- c("MFpower","ICPCpfc","l_ICPCmotor","r_ICPCmotor",
           "IVleft","IVright","trial2use","N","Nsub","ya","s","r","rew")
  stopifnot(all(req %in% names(E)))
  
  ## ---------- Basic dimensions ----------
  nTrial <- as.integer(E$N[1,1])
  nSub   <- as.integer(E$Nsub[1,1])
  nResp  <- 3L
  nStim  <- 8L
  nData  <- nSub * nTrial
  
  ## ---------- Behavior ----------
  ya <- as.matrix(E$ya)            # expected values: {0,1,2}
  ya[is.na(ya)] <- 0
  ya <- ifelse(ya == 0, 3L, ya)    # 0 (NoGo) -> 3 (action index: 1/2/3)
  
  # a: Go/NoGo binary coding (Go=1, NoGo=2)
  ab <- ifelse(ya == 3L, 2L, 1L)
  
  sb   <- as.matrix(E$s)           # stimulus index: 1..8
  rb   <- as.matrix(E$r)           # feedback (-1/0/1)
  rewb <- as.matrix(E$rew)         # 0/1: Avoid/Win
  rb[is.na(rb)] <- 0
  rewb <- ifelse(rewb == 1, 1L, 2L)  # -> 1/2 (Win=1, Avoid=2)
  
  # Vectorized actions (concatenate trials by subject)
  yax <- as.integer(as.vector(t(ya)))
  
  ## ---------- Read EEG matrices ----------
  MFpower      <- as.matrix(E$MFpower)
  ICPCpfc      <- as.matrix(E$ICPCpfc)
  l_ICPCmotor  <- as.matrix(E$l_ICPCmotor)
  r_ICPCmotor  <- as.matrix(E$r_ICPCmotor)
  IVleft       <- as.matrix(E$IVleft)
  IVright      <- as.matrix(E$IVright)
  trial2use    <- as.matrix(E$trial2use)
  
  # If shape is nTrial x nSub, transpose to nSub x nTrial
  fix_shape <- function(M) {
    if (is.null(dim(M))) return(M)
    if (identical(dim(M), c(nTrial, nSub))) return(t(M))
    M
  }
  MFpower     <- fix_shape(MFpower)
  ICPCpfc     <- fix_shape(ICPCpfc)
  l_ICPCmotor <- fix_shape(l_ICPCmotor)
  r_ICPCmotor <- fix_shape(r_ICPCmotor)
  IVleft      <- fix_shape(IVleft)
  IVright     <- fix_shape(IVright)
  trial2use   <- fix_shape(trial2use)
  
  ## ---------- Monotonic inverse transforms ----------
  inv_pow_transform <- function(M) {
    M <- as.matrix(M)
    idx <- is.finite(M)
    if (any(idx)) {
      minv <- min(M[idx])
      M[idx] <- 1 - 1/(M[idx] + (1 - minv))
    }
    M
  }
  inv_icpc_transform <- function(M) {
    M <- as.matrix(M)
    idx <- is.finite(M)
    if (any(idx)) {
      maxv <- max(M[idx])
      M[idx] <- 1/(-M[idx] + maxv + 1)
    }
    M
  }
  
  # Apply inverse transforms first and keep raw copies for checking
  tPow_raw  <- inv_pow_transform(MFpower)
  pfc_raw   <- inv_icpc_transform(ICPCpfc)
  lICPC_raw <- inv_icpc_transform(l_ICPCmotor)
  rICPC_raw <- inv_icpc_transform(r_ICPCmotor)
  
  # Replace non-finite values with zero
  tPow   <- tPow_raw;   tPow   [!is.finite(tPow)]    <- 0
  pfcICPC<- pfc_raw;    pfcICPC[!is.finite(pfcICPC)] <- 0
  lICPC  <- lICPC_raw;  lICPC  [!is.finite(lICPC)]   <- 0
  rICPC  <- rICPC_raw;  rICPC  [!is.finite(rICPC)]   <- 0
  
  # trial2use -> 0/1 mask; zero-out invalid trials if requested
  trial2use01 <- ifelse(is.finite(trial2use) & trial2use != 0, 1L, 0L)
  if (zero_mask) {
    tPow   [trial2use01 == 0] <- 0
    pfcICPC[trial2use01 == 0] <- 0
    lICPC  [trial2use01 == 0] <- 0
    rICPC  [trial2use01 == 0] <- 0
  }
  
  ## ---------- Conflict labels for the 8 stimuli ----------
  # 1 = conflict, 2 = non-conflict
  motivconflict <- as.integer(c(2,2,1,1,1,1,2,2))
  
  ## ---------- Force integer storage ----------
  storage.mode(ya)   <- "integer"
  storage.mode(ab)   <- "integer"
  storage.mode(sb)   <- "integer"
  storage.mode(rb)   <- "integer"
  storage.mode(rewb) <- "integer"
  storage.mode(yax)  <- "integer"
  nTrial <- as.integer(nTrial); nSub <- as.integer(nSub)
  nResp  <- as.integer(nResp);  nStim <- as.integer(nStim)
  nData  <- as.integer(nData)
  
  ## ---------- Pack output ----------
  dList <- list(
    # Behavior
    ya = ya, yax = yax, a = ab, s = sb, r = rb, rew = rewb,
    nTrial = nTrial, nSub = nSub, nData = nData, nResp = nResp, nStim = nStim,
    # EEG (transformed matrices)
    tPow = tPow, pfcICPC = pfcICPC, lICPC = lICPC, rICPC = rICPC,
    # Additional helpers
    trial2use = trial2use01, IVleft = IVleft, IVright = IVright,
    motivconflict = motivconflict
  )
  if ("groupVec" %in% names(E)) dList$group <- as.integer(as.vector(E$groupVec))
  if ("subIDs"   %in% names(E)) dList$subID <- as.integer(as.vector(E$subIDs))
  
  save(dList, file = obs4stanfile)
  
  if (verbose) {
    cat("\n== Save completed: ", obs4stanfile, " ==\n", sep = "")
    cat("lICPC_raw range (finite): ", range(lICPC_raw[is.finite(lICPC_raw)]), "\n")
    cat("rICPC_raw range (finite): ", range(rICPC_raw[is.finite(rICPC_raw)]), "\n")
    cat("Proportion of non-zero lICPC at unmasked positions: ",
        mean(lICPC_raw[trial2use01==1] != 0, na.rm = TRUE), "\n")
    cat("Proportion of non-zero rICPC at unmasked positions: ",
        mean(rICPC_raw[trial2use01==1] != 0, na.rm = TRUE), "\n")
    cat("Number of trial2use==0: ", sum(trial2use01==0), " / ", length(trial2use01), "\n", sep="")
  }
}