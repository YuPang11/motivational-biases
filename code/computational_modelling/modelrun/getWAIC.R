getWAIC <- function (stanfitsfile,waicfile,model2fit){
  
  waicname = paste("M",model2fit,"WAIC",sep="")
  
  if(!file.exists(waicfile)){
    
    # calculate WAIC.
    load(stanfitsfile)
    f2@.MISC<-new.env(parent=globalenv())
    modelWAIC <- list(waicname=waicname,waic=waic(f2),mstr=mstr)
    rm(f2) # clear from memory.
    
    # store waic with model.
    save(modelWAIC,file=waicfile)
    
  } else {
    load(waicfile)
  } # end if waicfile exists.
  
  # print relevat information.
  print(modelWAIC$waic$waic)
  
} # end getWAIC.