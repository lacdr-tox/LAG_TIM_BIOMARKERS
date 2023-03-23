parnames <- c("k_LAG","k_TIM","k_PD",
              "k_inf","k_exp","de",
              "d_lag","d_tim","d_pd1","d_pdl1", "d_ifn",
              "k_ifn","ke")

b16 <- function(t,state,parms){
  with(as.list(c(state,parms)), {
    
    k_A <- exp(1) #fixed
    k_exh <- exp(1) #fixed
    
    SIZE <- (G1+SG2M)/(10^5) ## if we say TC density is 10^5 TC/s mm^-3, then this is  the tumour volume in mm^3
    
    ## exhaustion status of CTLs, for simplicity this is a weighted sum of all the checkpoints
    E_exh <- k_LAG*LAG3/E + k_TIM*TIM3/E + k_PD*PD1/E*(PDL1)/SIZE 
    
    ## activity is inversely proportionaly to exhaistion
    E_activ <- 1-E_exh/(E_exh+1)
    if(is.na(E_activ)) E_activ <- 1 ##necessary because E starts at 0.
    
    dE <- 0
    if(t>0) dE <- k_inf*SIZE+k_exp*E_activ*E-de*E  ## activated phenotype (IFNg secretion)
    
    ##  IFNg increases basal antigenicity (CONCENTRATION DEPENDENT)
    A_total <- 1+IFN/SIZE
    if(E>0){ ## IFN produced by t cells, concentration of checkpoint molecules on t cells
      ## slight problem in that the basal level of checkpoint expression might not be zero...
      dIFN <- E*E_activ - d_ifn*IFN
      dLAG3 <- E*A_total-LAG3*d_lag
      dTIM3 <- E*A_total-TIM3*d_tim 
      dPD1 <- E*A_total-PD1*d_pd1 
      dPDL1 <- E*A_total-PDL1*d_pdl1
    }
    if(E<=0){
      dIFN <- 0
      dPD1 <- 0
      dPDL1 <- 0
      dTIM3 <- 0
      dLAG3 <- 0
    }
    
    
    
    dSG2M <- krg*G1/(1+k_ifn*IFN/SIZE) - kgr*SG2M - ke * E *E_activ*SG2M/(SG2M+G1)
    dG1 <- -krg*G1/(1+k_ifn*IFN/SIZE) +2*kgr*SG2M - ke *E_activ* E *G1/(SG2M+G1)
    
    return(list(c(dE, dLAG3,dTIM3,dPD1,dPDL1, dIFN, dSG2M, dG1)))
    
  })
}

run_model <- function(pars,g0=0.4,GR0=1.1,TCdensity=10^5,V0=25.56,times=seq(-1,15,0.1),expon=TRUE,fixpars=NULL){
  
  TC0 <- V0*TCdensity
  SG2M <- TC0*GR0/(1+GR0)
  G1 <- TC0-SG2M
  state <- c(E=0,LAG3=0,TIM3=0,PD1=0,PDL1=0,IFN=0,SG2M=SG2M,G1=G1)  
  if(expon==TRUE) pars <- exp(pars)
  
  parnames <- c("k_LAG","k_TIM","k_PD",
    "k_inf","k_exp","de",
    "d_lag","d_tim","d_pd1","d_pdl1", "d_ifn",
    "k_ifn","ke")
  if(!is.null(fixpars)) parnames <- parnames[!parnames%in%names(fixpars)]
  names(pars) <- parnames
  pars <- c(pars,fixpars)
  
  pars <- c(pars,get.growth(g0,GR0))
  out <- data.frame(ode(y=state,times=times,func = b16,parms = pars))
  return(out)
}

transform_output <- function(df,dat_times){
  ids <- unlist(lapply(1:length(dat_times), function(i) rep(names(dat_times)[i],length(dat_times[[i]]))))[-1]
  yctl <- df$E[df$time%in%dat_times$ctl_dat]
  ylag <- df$LAG3[df$time%in%dat_times$lag_dat]
  ytim <- df$TIM3[df$time%in%dat_times$tim_dat]
  ypd1 <- df$PD1[df$time%in%dat_times$pd1_dat]
  ypdl1 <- df$PDL1[df$time%in%dat_times$pdl1_dat]
  yifn <- df$IFN[df$time%in%dat_times$ifn_dat]
  ygr <- (df$SG2M/df$G1)[df$time%in%dat_times$GR_dat]
  vg <- rowSums(df[df$time%in%c(-1,1,3,5,7,10,14),c("SG2M","G1")])
  yg <- log(vg[2:7]/vg[1:6])/diff(c(-1,1,3,5,7,10,14))## estimate growth rate between treatment points
  #note we have zero sd for the data at the first ctl timepoint, so it is removed
  y <- list(yctl=yctl[-1],ylag=ylag,ytim=ytim,ypd1=ypd1,ypdl1=ypdl1,yifn=yifn,ygr=ygr,yg=yg)
  y <- unlist(y)  
  y <- data.frame(y=as.numeric(y),ids=ids)
  rownames(y) <- NULL
  
  return(y)
}

transform_dat <- function(dat){
  dat <- dat[c("ctl_dat","lag_dat","tim_dat","pd1_dat","pdl1_dat","ifn_dat","GR_dat","v_dat")]
  
  values <- lapply(dat,function(di) di$value)
  sd <- lapply(dat,function(di) di$sd)
  day <- lapply(dat,function(di) di$day)
  
  ids <- lapply(1:length(dat), function(i) rep(names(dat)[i],nrow(dat[[i]])))
  
  dfd <- data.frame(day=unlist(day),mean=unlist(values),sd=unlist(sd),ids=unlist(ids))
  dfd$sd[is.na(dfd$sd)]<-0.054+0.29*dfd$mean[is.na(dfd$sd)]
  dfd <- dfd[!dfd$sd==0,]#note this removes the data at the first ctl timepoint
  rownames(dfd) <- NULL
  
  rna_ids <- c("tim_dat","lag_dat","ifn_dat","pd1_dat","pdl1_dat")
  
  for(id in rna_ids){
    dfd$mean[dfd$ids==id] <- dfd$mean[dfd$ids==id]-min(dfd$mean[dfd$ids==id])
  }
  
  return(dfd)
}


get_activity <- function(x, pars, fixpars=NULL){
  pars <- exp(as.numeric(pars))
  parnames <- c("k_LAG","k_TIM","k_PD",
                "k_inf","k_exp","de",
                "d_lag","d_tim","d_pd1","d_pdl1", "d_ifn",
                "k_ifn","ke")
  if(!is.null(fixpars)) parnames <- parnames[!parnames%in%names(fixpars)]
  names(pars) <- parnames
  pars <- c(pars,fixpars)
  
  df <- data.frame(pd=pars['k_PD']*x$PD1/x$E*(x$PDL1)*10^5/(x$G1+x$SG2M),
             lag=pars['k_LAG']*x$LAG3/x$E,
             tim=pars['k_TIM']*x$TIM3/x$E)
  df$E_exh <- rowSums(df)
  df$E_activ <- 1-1/(1+(1/df$E_exh))
  df$time <- x$time
  df
}

get_ll <- function(pars,dat,dat_times,SIZE0=10^5,
                   reject.waves=TRUE,fixpars=NULL){

  x <- tryCatch(run_model(pars,expon = T,fixpars = fixpars),error=function(e) return(10^9))
  if(is.character(x)) return(10^5)
  if(is.null(nrow(x))) return(10^5)
  activity <- get_activity(x,pars,fixpars)$E_activ
  if(max(diff(activity[!is.na(activity)]))>0&reject.waves) return(10^3)
  X <- transform_output(x,dat = dat_times)
  
  rna_ids <- c("tim_dat","lag_dat","ifn_dat","pd1_dat","pdl1_dat")
  
  rescalr <- function(scalr,inpt,dt){
    sum((inpt*exp(scalr)-dt$mean)^2/dt$sd^2)/2
  }
  
  for(id in rna_ids){
    inpt <- X$y[X$ids==id]
    dt <- dat[dat$ids==id,]
    opt <- optimise(f=rescalr,interval = c(-20,20),inpt=inpt,dt=dt)$minimum
    X$y[X$ids==id]<-X$y[X$ids==id]*exp(opt)
  }
  
  
  errs <- (X$y-dat$mean)^2/dat$sd^2
  errs <- errs[-length(errs)] ## remove last growth point from fit
  err <- sum(errs)/2
  if(is.na(err)) return(10^5)
  return(err)
}

get.growth <- function(g,GR){ 
  kgr <- g*(1+GR)/GR
  krg <- g*(1+2*GR)
  c(krg=krg,kgr=kgr)
}

setup_opt <- function(path=""){
  minval <- 0.001
  maxval <- 1000
  lower <-  c("k_LAG"=minval,"k_TIM"=minval,"k_PD"=minval,
              "k_inf"=minval,"k_exp"=minval,"de"=0.01,
              "d_lag"=0.01,"d_tim"=0.01,"d_pd1"=0.01,"d_pdl1"=0.01, "d_ifn"=0.01,
              "k_ifn"=minval,"ke"=0.01)
  
  upper <- c("k_LAG"=maxval,"k_TIM"=maxval,"k_PD"=maxval,
             "k_inf"=maxval,"k_exp"=maxval,"de"=10,
             "d_lag"=100,"d_tim"=100,"d_pd1"=100,"d_pdl1"=100, "d_ifn"=100,
             "k_ifn"=maxval,"ke"=100)
  
  df_p0 <- data.frame(parname=names(upper),lower,upper)
  
  dat <- readRDS(paste0(path,"fit_dat.rds"))
  dat <- dat[c("ctl_dat","lag_dat","tim_dat","pd1_dat","pdl1_dat","ifn_dat","GR_dat","v_dat")]
  
  dat_times <- lapply(dat,function(di) di$day)
  
  dat <- transform_dat(dat)
  
  list(df_p0=df_p0,dat=dat,dat_times=dat_times)
  
}
