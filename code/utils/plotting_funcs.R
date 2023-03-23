parnames <- c("k_LAG","k_TIM","k_PD",
              "k_inf","k_exp","de",
              "d_lag","d_tim","d_pd1","d_pdl1", "d_ifn",
              "k_ifn","ke")
run_and_rescale <- function(pars,d0,SIZE0=10^5,fixpars=NULL){
  dat <- d0$dat
  dat_times <- d0$dat_times
  x <- run_model(pars,expon = T,fixpars=fixpars)
  X <- transform_output(x,dat = dat_times)
  
  rna_ids <- c("tim_dat","lag_dat","ifn_dat","pd1_dat","pdl1_dat")
  
  rescalr <- function(scalr,inpt,dt){
    sum((inpt*exp(scalr)-dt$mean)^2/dt$sd^2)/2
  }
  
  par_scalrs <- sapply(rna_ids, function(id){
    inpt <- X$y[X$ids==id]
    dt <- dat[dat$ids==id,]
    opt <- optimise(f=rescalr,interval = c(-20,20),inpt=inpt,dt=dt)$minimum
    return(opt)
  })
  
  x$TIM3 <- x$TIM3*exp(par_scalrs["tim_dat"])
  x$LAG3 <- x$LAG3*exp(par_scalrs["lag_dat"])
  x$PD1 <- x$PD1*exp(par_scalrs["pd1_dat"])
  x$PDL1 <- x$PDL1*exp(par_scalrs["pdl1_dat"])
  x$IFN <- x$IFN*exp(par_scalrs["ifn_dat"])
  
  list(x=x,par_scalrs=par_scalrs)
}

reshape_model <- function(x,SIZE0=10^5){
  modTC <- x$SG2M+x$G1
  modV <- modTC/SIZE0
  gr<-x$SG2M/x$G1
  
  y <- log(modV[2:length(modV)]/modV[1:(length(modV)-1)])/diff(x$time)## estimate growth rate between treatment points
  mod_g <- c(y,NA)
  
  x <- x[,c("time","E","LAG3","TIM3","PD1","PDL1","IFN")]
  x$g <- mod_g
  x$gr <- gr
  x <- reshape2::melt(x,id.vars="time")
  return(x)
}

reshape_data <- function(dat){
  renamr <- c("ctl_dat" ="E", "lag_dat" = "LAG3", "tim_dat"="TIM3",
    "pd1_dat"="PD1", "pdl1_dat" = "PDL1", "ifn_dat"="IFN",
    "GR_dat" ="gr", "v_dat" = "g")
  dat$variable <- renamr[dat$ids]
  return(dat)
}

load_fit <- function(id, pb,fixpars=NULL,path="data/PAR_EST/"){
  parnames <- c("k_LAG","k_TIM","k_PD",
                "k_inf","k_exp","de",
                "d_lag","d_tim","d_pd1","d_pdl1", "d_ifn",
                "k_ifn","ke")
  if(!is.null(fixpars)) parnames <- parnames[!parnames%in%names(fixpars)]
  names(pb) <- parnames
  #pb <- c(pb,fixpars)
  
  ### TESTS IF ANY OF THESE PARAMETERS HITTING THE LOWER BOUNDS MATTER (THEY DON'T).
  #if(id=="TIM\ndominant"){
  # pb["k_LAG"] <- -30
  #pb["d_tim"] <- -30
  #}
  # if(id=="LAG\ndominant"){
  #pb["k_TIM"] <- -30
  #pb["d_lag"] <- -30
  #}
  
  #pb["k_inf"] <- -12
  
  d0 <- setup_opt(path = path)
  ll <- get_ll(pb,dat = d0$dat,dat_times=d0$dat_times,fixpars = fixpars)
  
  pbko <- pb
  pbko["k_PD"] <- -30
  llko <- get_ll(pbko,dat = d0$dat,dat_times=d0$dat_times,fixpars=fixpars)
  
  df <- run_and_rescale(pb,d0,fixpars=fixpars)$x
  df <- reshape_model(df)
  df$id <- id
  
  df_inh <- df[df$variable%in%c("LAG3","TIM3","PD1","PDL1"),]
  df <- df[!df$variable%in%c("LAG3","TIM3","PD1","PDL1"),]
  
  list(df=df,df_inh=df_inh,dfll=data.frame(ll,llko))
}

load_activity <- function(id, pars,fixpars=NULL){
  x <- run_model(pars,expon = T,fixpars=fixpars)
  activity <- get_activity(x,pars,fixpars)
  for(inh in c("pd","lag","tim")){
    activity[,inh] <- activity[,inh]/activity[,"E_exh"]
  }
  activity <- reshape2::melt(activity,id.vars=c("E_exh","E_activ","time"))
  activity$id <- id
  activity
}

## function sweeps a parameter and then calculates the likelihood
## note there is no re-optimisation.
sweep_par<- function(pb,sweep="d_pd1",range=seq(-7,7,0.1),path="data/PAR_EST/"){
  parnames <- c("k_LAG","k_TIM","k_PD",
                "k_inf","k_exp","de",
                "d_lag","d_tim","d_pd1","d_pdl1", "d_ifn",
                "k_ifn","ke")

  names(pb) <- parnames
  #p0 <- pb[sweep]
  d0 <- setup_opt(path = path)
  ll <- sapply(range, function(ri){
    pb[sweep] <- ri
    get_ll(pb,dat = d0$dat,dat_times=d0$dat_times)
  })
  
  data.frame(swept=sweep,parval = range,ll=ll)
}
