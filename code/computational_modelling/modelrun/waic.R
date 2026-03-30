# WAIC...
# From Aki & Gelman 2014
# Little function to calculate posterior variances from simulation
colVars <- function (a){
  diff <- a - matrix (colMeans(a), nrow(a), ncol(a), byrow=TRUE)
  vars <- colMeans (diff^2)*nrow(a)/(nrow(a)-1)
  return (vars)
}
log_sum_exp <- function(x){
  x_max <- max(x)
  x_max + log(sum(exp(x-x_max)))
}

# The calculation of Waic!  Returns lppd, p_waic_1, p_waic_2, and waic, which we define
# as 2*(lppd - p_waic_2), as recommmended in BDA
waic <- function(stanfit){
  log_lik <- extract (stanfit, "log_lik")$log_lik
  dim(log_lik) <- if (length(dim(log_lik))==1) c(length(log_lik),1) else
    c(dim(log_lik)[1], prod(dim(log_lik)[2:length(dim(log_lik))]))
  S <- nrow(log_lik)
  n <- ncol(log_lik)
  
  summands <- apply(log_lik, 2, log_sum_exp)
  correc   <- -ncol(log_lik) * log(nrow(log_lik))
  lppd     <- sum(summands) + correc
  p_waic_1 <- 2*(sum(summands - colMeans(log_lik)) + correc)
  p_waic_2 <- sum(colVars(log_lik))
  waic_2   <- -2*lppd + 2*p_waic_2
  
  lpd       <- log(colMeans(exp(log_lik)))
  p_waic    <- colVars(log_lik)
  elpd_waic <- lpd - p_waic
  waic      <- -2*elpd_waic
  loo_weights_raw <- 1/exp(log_lik-max(log_lik))
  loo_weights_normalized <- loo_weights_raw/
    matrix(colMeans(loo_weights_raw),nrow=S,ncol=n,byrow=TRUE)
  loo_weights_regularized <- pmin (loo_weights_normalized, sqrt(S))
  elpd_loo <- log(colMeans(exp(log_lik)*loo_weights_regularized)/
                    colMeans(loo_weights_regularized))
  p_loo <- lpd - elpd_loo
  pointwise <- cbind(waic,lpd,p_waic,elpd_waic,p_loo,elpd_loo)
  total <- c(colSums(pointwise), waic_2=waic_2, p_waic_1=p_waic_1, p_waic_2=p_waic_2, lppd=lppd)
  se <- sqrt(n*colVars(pointwise))
  return(list(waic=total["waic"], elpd_waic=total["elpd_waic"],
              p_waic=total["p_waic"], elpd_loo=total["elpd_loo"], p_loo=total["p_loo"],
              waic_2=total["waic_2"], p_waic_1=total["p_waic_1"], p_waic_2=total["p_waic_2"], 
              lppd=total["lppd"],
              pointwise=pointwise, total=total, se=se))
}
