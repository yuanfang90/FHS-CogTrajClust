multcp.step <- function(mod,nrep,niter,n.cp.in){# n.cp = number of change points; n.cov = number of other covariates including intercept
  sum <- data.frame(summary(mod)) # grab the summary table
  fmla <- as.character(mod$call) # parse the formula of the model
  G <- ifelse(nchar(fmla[4])>10,as.numeric(fmla[6]),1) ##fmla[4] will store "random" if is multiple-component model and fmla[6] will store number of classes in that case
  
  parnames <- row.names(sum) # stores names of all parameters estimated in mod, will not change unless one or more random effect change points were removed for all classes
  cp.coef.ind <- which(startsWith(parnames,prefix="age.3")) # specific for my way of setting up change pt variables. Once parnames changes, this need to be changed too
  pars.cp <- parnames[cp.coef.ind] # parameter names for the change points, may be multiple classes. Need to update this if parnames changed. If G=1, the same as cp.names defined below; if G > 1, stores as the format "cpname classg".
  sum1 <- sum[cp.coef.ind,] # extract the summary table for change points only
  pval <- sum1$p.value # extract pvalues for change points. This and the sum1 table need to be updated each iteration
  
  n.cp <- n.cp.in # number of change points
  cpclasses <- vector(mode = "list", length = n.cp.in) # a list with the same length as number of change points, each element is a vector contains the G parameter for that change point
  if(G!=1){
    for(i in 1:n.cp.in){
      cpclasses[[i]] <- pars.cp[((i-1)*G+1):(i*G)]
    }
    cp.names <- unique(sub(" class.*", "", pars.cp)) # grab the names of the change points
  }else{
    cp.names <- unique(pars.cp)
  }
  cp.names.in <- cp.names # define initial cp names
  
  binit <- mod$best # B for the model
  forforposfix <- rep(0,length(binit))
  
  m <- nrep
  s <- niter
  
  mod.up <- mod
  steps <- 1
  while(steps!=0){
    cat("round ",steps,"\n")
    mod.old <- mod.up

    ### extract p-values for the change points
    p <- pval
    rm.ind <- which.max(p) # notice that pval have the same length as pars.cp but not necessary cp.names, not the same as parnames
    
    ### make the initial value for the change point with biggest p-value 0
    rm <- pars.cp[rm.ind]
    binit[rm] <- 0
    fix.ind <- which(names(binit)==rm)
    forforposfix[fix.ind] <- fix.ind
    forposfix <<- as.numeric(forforposfix[!(forforposfix==0)])
    cat("removing",rm,"with p-value",p[rm.ind], "\n")
    
    if(nchar(fmla[4])>10){ # for G > 1. fmla[4] will store "random" if is multiple-component model, fmla[3] is "mixture"
      for(i in 1:n.cp.in){ # this list should be of the length of initial number of change points
        for(g in 1:G){
          if(cpclasses[[i]][g]==rm) cpclasses[[i]][g] <- paste0(cpclasses[[i]][g],"-removing")
        }
      }
      ##### check if if a change point for all classes will be removed,
      ##### if so, refit model with removing random effect with this change point too
      for(i in 1:n.cp.in){
        if(all(endsWith(cpclasses[[i]],suffix="removing"))){
          rm.random <- cp.names.in[i]
          rm.fixeff <- sub("-removing", "", cpclasses[[i]])
          cpclasses[[i]] <- paste0(cpclasses[[i]],"-removed")
        }
      }
      if(exists("rm.random")){ # only proceed the following if any cp for all classes were removed at this step
        #### update formula######## update formula
        cp.names.after.rem <- cp.names[!(cp.names==rm.random)]
        fmla[4] <- paste0("~age.1 + ",stringr::str_c(cp.names.after.rem, collapse = " + ")) # random effect formula
        fmla[2] <- paste0(sub(" ~.*", "", fmla[2])," ~ age.1 + ",stringr::str_c(cp.names.after.rem, collapse = " + ")," + SEX + EDUcomp") # fixed effect formula
        fmla[3] <- paste0("~age.1 + ",stringr::str_c(cp.names.after.rem, collapse = " + ")," + SEX + EDUcomp") # mixture formula
        
        cp.ind <- which(cp.names == rm.random)
        ######## update random effects initial values
        if(cp.ind == n.cp){
          random.ind <- c((((cp.ind+2)*(cp.ind+1)*0.5+1):((cp.ind+2)*(cp.ind+3)*0.5)))
        }else{
          random.ind <- c((((cp.ind+2)*(cp.ind+1)*0.5+1):((cp.ind+2)*(cp.ind+3)*0.5)),((cp.ind+2)*(cp.ind+3)*0.5)+cumsum(c((cp.ind+2):(n.cp+1))))
        }
        varcov.rm <- paste("varcov",random.ind)
        binit <- binit[which(!names(binit)%in%varcov.rm)]
        change.varcov.name.ind <- which(names(binit)== paste("varcov",(cp.ind+2)*(cp.ind+1)*0.5))
        names(binit)[(change.varcov.name.ind+1):(length(binit)-1)] <- paste("varcov",((cp.ind+2)*(cp.ind+1)*0.5+1):((cp.ind+2)*(cp.ind+1)*0.5+(length(binit)-1)-change.varcov.name.ind))
        
        ######## update fix effect initial values
        # fixed.rm <- which(names(binit)%in%rm.fixeff)
        binit <- binit[which(!names(binit)%in%rm.fixeff)]
        forforposfix <- rep(0,length(binit))
        fixed.rm <- which(binit==0)
        forforposfix[fixed.rm] <- fixed.rm
        forposfix <<- as.numeric(forforposfix[!(forforposfix==0)])
        # names(forposfix) <- NULL
        rm("rm.random","cp.ind","rm.fixeff","fixed.rm")
      }
      
      loglik <- numeric(m)
      gridsearch.B.list <- vector(mode <- "list", length <- m)
      cat(paste0(fmla[6],"-component model."),"Gridsearch initialization with",m,"random starts",s,"iterations each is running ... ")
      
      for(rep in 1:m){
        mod1 <- hlme(fixed=as.formula(fmla[2]),random=as.formula(fmla[4]),mixture=as.formula(fmla[3]),subject=fmla[5],ng=G,data=model.data,maxiter=s,B=binit, posfix=forposfix, verbose=FALSE)
        loglik[rep] <- mod1$loglik
        gridsearch.B.list[[rep]] <- mod1$best
        if(rep==m){cat("*","\n")}else{cat("*")}
      }
      cat("Initialization done. ")
      binit.up =gridsearch.B.list[[which.max(loglik)]]
      mod.up <- hlme(fixed=as.formula(fmla[2]),random=as.formula(fmla[4]),mixture=as.formula(fmla[3]),subject=fmla[5],ng=G,data=model.data,B=binit.up,maxiter=1000,posfix=forposfix, verbose=TRUE)
    }else{ #fmla[4] will store "subject" if is 1-component model, "random" is in fmla[3]
      # mod.old <- mod.up
      cat("1-component model. ")
      cp.ind <- which(cp.names == rm)
      cp.names.after.rem <- cp.names[!(cp.names==rm)]
      fmla[3] <- paste0("~age.1 + ",stringr::str_c(cp.names.after.rem, collapse=" + "))
      # fmla[3] <- sub(rm.random, "", fmla[3])
      
      if(cp.ind == n.cp){
        random.ind <- c((((cp.ind+2)*(cp.ind+1)*0.5+1):((cp.ind+2)*(cp.ind+3)*0.5)))
      }else{
        random.ind <- c((((cp.ind+2)*(cp.ind+1)*0.5+1):((cp.ind+2)*(cp.ind+3)*0.5)),((cp.ind+2)*(cp.ind+3)*0.5)+cumsum(c((cp.ind+2):(n.cp+1))))
      }
      varcov.rm <- paste("varcov",random.ind)
      binit <- binit[which(!names(binit)%in%varcov.rm)]
      change.varcov.name.ind <- which(names(binit)== paste("varcov",(cp.ind+2)*(cp.ind+1)*0.5))
      names(binit)[(change.varcov.name.ind+1):(length(binit)-1)] <- paste("varcov",((cp.ind+2)*(cp.ind+1)*0.5+1):((cp.ind+2)*(cp.ind+1)*0.5+(length(binit)-1)-change.varcov.name.ind))
      
      mod.up <- hlme(fixed=as.formula(fmla[2]),random=as.formula(fmla[3]),subject=fmla[4],data=model.data,B=binit,maxiter=1000,posfix=forposfix,verbose=TRUE)
    } 
    
    ##### update everything after fitting the new model
    sum <- data.frame(summary(mod.up))
    parnames <- row.names(sum)
    cp.coef.ind <- which(startsWith(parnames,prefix="age.3"))
    pars.cp <- parnames[cp.coef.ind]
    sum1 <- sum[cp.coef.ind,]
    pval <- sum1$p.value
    if(G!=1){
      cp.names <- unique(sub(" class.*", "", pars.cp)) # grab the names of the change points
    }else{
      cp.names <- unique(cp.names.after.rem)
    }
    n.cp <- length(cp.names)
    if(all(is.na(pval))|all(na.omit(pval)<0.05)){
      steps <- 0
      break
    }else{
      steps <- steps+1
      binit <- mod.up$best # B for the model
    }
  } # end of while loop
  if(all(is.na(pval))){return(mod.old)}else{return(mod.up)}
  # return(mod.up)
} # end of function