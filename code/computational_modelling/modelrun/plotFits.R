plotFits <- function(estimatefile,plotfitsfile,model2fit) {
  
  # load parameter estimates.
  load(estimatefile)
  
  # print relevant stats/diagnostics.
  # summary statistics.
  print("means:")
  print(apply(xbay, 2, mean)) # mean 
  print("SD:")
  print(apply(xbay, 2, sd))# standard deviation 
  
  
  # PLOTTING.
  # ===================================================================.
  col.corrgram <-function(ncol){colorRampPalette(c("darkgoldenrod4","burlywood1","darkkhaki","darkgreen"))(ncol)}
  
  # correlation of parameters.
  filename = paste(plotfitsfile,"_corrgram.png",sep="")
  if (!file.exists(filename)){
    X11(); corrgram(xbay,order=F,upper.panel = panel.pts,lower.panel=panel.shade,text.panel=panel.txt,main="Correlogram of fitted parameters",col.regions = col.corrgram)
    savePlot(filename = filename)}
  
  # posterior distributions.
  filename = paste(plotfitsfile,"_posterior.png",sep="")
  if (!file.exists(filename)){
    load(stanfitsfile)
    X11(); par(mfrow=c(4,2), mar=c(4,4,2,1))
    for(iK in 1:K){
      parname = paste("X[",iK,"]",sep="")
      post <- unlist(extract(f2, parname), use.names=FALSE)
      plot(density(post),xlab=parname, col=grey(0, 0.8),main=paste("Model",model2fit),ylab="posterior density")
    } # end iK-loop.
    savePlot(filename = filename)}    
  
} # end plot fits.
