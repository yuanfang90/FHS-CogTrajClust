gridsearch3 <- function(nrep,niter,minit,G,mod.dat){
  fmla <- as.character(minit$call) # parse the formula of the model
  loglik <- numeric(nrep)
  gridsearch.B.list <- vector(mode <- "list", length <- nrep)
  cat(paste0(G,"-component model."),"Gridsearch initialization with",nrep,"random starts",niter,"iterations each is running ... ")
  binit = minit$best
  for(rep in 1:nrep){
    mod1 <- hlme(fixed=as.formula(fmla[2]),random=as.formula(fmla[3]),mixture=as.formula(paste("~",sub(".*~ ", "", fmla[2]))),subject=fmla[4],ng=G,data=mod.dat,maxiter=niter,B=random(minit), verbose=FALSE)
    loglik[rep] <- mod1$loglik
    gridsearch.B.list[[rep]] <- mod1$best
    if(rep==nrep){cat("*","\n")}else{cat("*")}
  }
  cat("Initialization done. ")
  binit.up = gridsearch.B.list[[which.max(loglik)]]
  return(binit.up)
}


