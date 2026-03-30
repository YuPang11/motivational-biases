extractEstimates <- function(stanfitsfile,estimatefile,model2fit) {
  
  # retrieve parameter estimates.
  # ===================================================================.
  # load parameter estimates.
    load(stanfitsfile)
    f2@.MISC<-new.env(parent=globalenv())
    
    # retrieve estimates per parameter.
    fmeans  <- colMeans(as.matrix(f2))
    rm(f2)
    ix1     <- grep("x1", names(fmeans))
    xbay    <- data.table(x1=fmeans[ix1])
    
    # loop over parameters.
    for(iK in 2:dList$K){ 
      # parameter and index name.
      xname     = paste("x",iK,sep="")
      ixname    = paste("i",xname,sep="")
      
      # get index of parameter values.
      assign(ixname, grep(xname,names(fmeans))) # examp: iK=2 -> all indexes of 'x2' in fmeans are stored in variable 'ix2'.
      
      # get parameter values.
      xbay[, paste("x",iK,sep="") := data.table(fmeans[get(ixname)])] # examp: iK=2 -> all x2 values are stored in xbay$x2.
      
    } # end for iK-loop.
    
    K = dList$K
    
    # save parameter estimates.
    save(xbay,K,file=estimatefile)
    
} # end extractEstimates.