fitmodel <- function(model2fit, dat4stanfile, stanfitsfile, calcLogLik = TRUE){
  # Load prepared data for Rstan.
  load(dat4stanfile)  # -> dList
  stopifnot(exists("dList"))
  
  # 1) More robustly check whether K comes from the parent environment
  #    (avoids interference from same-named local variables)
  if (!exists("K", inherits = TRUE)) {
    if (!is.null(dList$K)) {
      K <- as.integer(dList$K)
    } else {
      K <- if (model2fit > 8) 6L else 5L
    }
  }
  # 2) Force integer type regardless of the source (critical)
  K <- as.integer(K)
  stopifnot(is.integer(K), length(K)==1, K >= 1, K <= 8)
  
  # Specify which parameters Stan should output.
  paramOut  <- c("X","sdX","x1","x2","x3","x4","x5","x6","x7","x8","log_lik")
  paramOut  <- paramOut[c(1:(2+K), length(paramOut))]     # same trimming logic as original
  paramInit <- paramOut[1:(length(paramOut)-1)]           # do not include log_lik for f1
  if(!calcLogLik){ paramOut <- paramInit }
  
  # Retrieve the model string to use.
  mstr <- get(paste("mstr", model2fit, sep=""))
  
  # Print some checks to make sure everything is OK.
  print(paste("model2fit =", model2fit))
  print(paste("Data file =", dat4stanfile))
  print(paste("K =", K))
  print(paste("paramOut =", paste(paramOut, collapse=", ")))
  print("STAN model string:\n")
  print(mstr)
  
  # MODEL FITTING + PARAMETER ESTIMATION.
  # ============================================================================.
  dList$K  <- K
  tmp      <- initQV(dList$nResp, dList$nStim)
  dList$Qi <- tmp$Qi
  dList$Vi <- tmp$Vi
  
  ## Sanity checks for behavioral dimensions and values (recommended to keep)
  stopifnot(
    is.matrix(dList$ya),  dim(dList$ya)  [1]==dList$nSub, dim(dList$ya)  [2]==dList$nTrial,
    is.matrix(dList$a),   dim(dList$a)   [1]==dList$nSub, dim(dList$a)   [2]==dList$nTrial,
    is.matrix(dList$s),   dim(dList$s)   [1]==dList$nSub, dim(dList$s)   [2]==dList$nTrial,
    is.matrix(dList$r),   dim(dList$r)   [1]==dList$nSub, dim(dList$r)   [2]==dList$nTrial,
    is.matrix(dList$rew), dim(dList$rew) [1]==dList$nSub, dim(dList$rew) [2]==dList$nTrial
  )
  # Optional extra safeguard for value ranges
  stopifnot(
    all(dList$ya  %in% 1:3),          # 1/2/3 = left/right/NoGo
    all(dList$a   %in% 1:2),          # Go/NoGo binary
    all(dList$rew %in% 1:2),          # Win=1, Avoid=2
    all(dList$r   %in% c(-1,0,1)),    # feedback
    all(dList$s   %in% 1:dList$nStim) # stimuli 1..8
  )
  
  ## --- Different from the original here:
  ##     no separate EEG loading, no outlier exclusion, subjects remain complete
  if (model2fit > 8) {
    need <- c("tPow","pfcICPC","lICPC","rICPC")
    miss <- setdiff(need, names(dList))
    if (length(miss)) stop("dList is missing EEG matrices: ", paste(miss, collapse=", "))
    # Add motivconflict if not already included
    # (kept consistent with your prepEEG4stan)
    if (!("motivconflict" %in% names(dList))) {
      dList$motivconflict <- as.integer(c(2,2,1,1,1,1,2,2))
    }
    stopifnot(length(dList$motivconflict) == dList$nStim)
    
    # Dimension consistency (same meaning as original)
    stopifnot(
      nrow(dList$tPow)    == dList$nSub,  ncol(dList$tPow)    == dList$nTrial,
      nrow(dList$pfcICPC) == dList$nSub,  ncol(dList$pfcICPC) == dList$nTrial,
      nrow(dList$lICPC)   == dList$nSub,  ncol(dList$lICPC)   == dList$nTrial,
      nrow(dList$rICPC)   == dList$nSub,  ncol(dList$rICPC)   == dList$nTrial
    )
    
    # Optional: make sure EEG values are zeroed where trial2use == 0
    # (if trial2use exists)
    if ("trial2use" %in% names(dList)) {
      mask0 <- dList$trial2use == 0
      stopifnot(all(dList$tPow   [mask0] == 0),
                all(dList$pfcICPC[mask0] == 0),
                all(dList$lICPC  [mask0] == 0),
                all(dList$rICPC  [mask0] == 0))
    }
  }
  
  # Equivalent to the original: recompute nData / yax
  # (keep the author's original row-major flattening behavior)
  dList$nData <- dList$nSub * dList$nTrial
  dList$yax   <- as.integer(as.vector(t(dList$ya)))
  
  # Specify seed for reproducibility.
  seed <- 69
  
  # Run a few iterations for compilation...
  f1 <- rstan::stan(model_code = mstr, data = dList, iter = 10, init = 0,
                    chains = 1, seed = seed, chain_id = 1, pars = paramInit)
  options(max.print = 2000)
  print(f1)
  
  # Main sampling
  f2 <- rstan::stan(fit = f1, data = dList,
                    warmup = 200, iter = 1200, chains = 4,
                    seed = seed + 1, pars = paramOut, init = 0)
  print(f2)
  
  if (!file.exists(stanfitsfile)) {
    f1@.MISC <- emptyenv()
    f2@.MISC <- emptyenv()
    save(dList, mstr, f1, f2, file = stanfitsfile)
  }
}
