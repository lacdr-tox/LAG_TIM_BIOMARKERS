require(deSolve)
library(pbapply)
#setwd("/home/4473331/projects/000_phd/rev_01b/optim/")
setwd("C:/Users/4473331/Documents/05_phd/02_CTLcomms/")
getE <- function(t){
  Edat <- data.frame(day=c(-1.01,1,3,5,7,10,16),value=c(0,0,10718,9667,4487,897,0)) ## boilerplate for convenience 
  x0 <- tail(Edat[Edat$day<t,],1)
  x1 <- head(Edat[Edat$day>=t,],1)
  
  dx <- (t-x0$day)/(x1$day-x0$day)
  y <- x0$value + dx*(x1$value - x0$value)
}

getIFN <- function(t){
  IFNdat <- data.frame(day=c(-1.01,1,3,5,7,10,16),value=c(0,0.0,1.5,0.3,0.1,0,0)) ## boilerplate for convenience 
  x0 <- tail(IFNdat[IFNdat$day<t,],1)
  x1 <- head(IFNdat[IFNdat$day>=t,],1)
  
  dx <- (t-x0$day)/(x1$day-x0$day)
  y <- x0$value + dx*(x1$value - x0$value)
}

b16 <- function(t,state,parms){
  with(as.list(c(state,parms)), {
    
    Tt <- G1+SG2M
    #SIZE <- Tt/TC0
    SIZE0 <- 10^5 ## if we say TC density is 10^5 TC/s mm^-3, then this is 1mm^3 
    SIZE <- Tt/SIZE0
    
    E <- getE(t)*SIZE
    IFN <- getIFN(t)*SIZE
    
    dSG2M <- krg*G1/(1+k_ifn*IFN/SIZE) - kgr*SG2M - ke*E*SG2M/(SG2M+G1)
    dG1 <- -krg*G1/(1+k_ifn*IFN/SIZE) +2*kgr*SG2M - ke*E *G1/(SG2M+G1)
    
    return(list(c(dSG2M, dG1)))
    
  })
}

get.growth <- function(g,GR){ 
  kgr <- g*(1+GR)/GR
  krg <- g*(1+2*GR)
  c(krg=krg,kgr=kgr)
}

run.model <- function(pars,g0=0.4,GR0=1.1,TCdensity=10^5,V0=25.56,times=seq(-1,15,0.1)){
  TC0 <- V0*TCdensity
  SG2M <- TC0*GR0/(1+GR0)
  G1 <- TC0-SG2M
  state <- c(SG2M=SG2M,G1=G1)  
  
  names(pars) <- c("k_ifn","ke")
  pars <- c(pars,get.growth(g0,GR0))
  out <- data.frame(ode(y=state,times=times,func = b16,parms = pars))
  return(out)
}

ssq <- function(df,dat,SIZE0=10^5){
  
  ygr <- (df$SG2M/df$G1)[df$time%in%dat$GR_dat$day]
  vg <- rowSums(df[df$time%in%c(-1,1,3,5,7,10,14),c("SG2M","G1")])
  yg <- log(vg[2:7]/vg[1:6])/diff(c(-1,1,3,5,7,10,14))
  
  yg <- yg[-length(yg)]
  dat$v_dat <- dat$v_dat[-nrow(dat$v_dat),]
  
  logl <- 1/2*sum((c(ygr,yg)-c(dat$GR_dat$value,dat$v_dat$value))^2/c(dat$GR_dat$sd,dat$v_dat$sd)^2)
  
  return(logl)
}

get_err <- function(pars){
  pdf <- readRDS("fit_dat.rds")
  x <- run.model(pars)
  err <- ssq(x,pdf)
  out <- data.frame(k_ifn=pars[1],k_e=pars[2],err=err)
  out
}

if(eval.sweep==FALSE){
  print("NOT EVALUATING SWEEP")
}else{
  pars <- expand.grid(seq(0,20,.1),seq(0,20,.1))
  
  cl <- makeCluster(getOption("cl.cores",75))
  clusterExport(cl=cl,varlist=c("ssq","run.model","get.growth","b16","getIFN","getE"))
  clusterCall(cl,function(dummyvar) library("deSolve"))
  
  pars <- lapply(1:nrow(pars), function(i) as.numeric(pars[i,]))
  
  out <- do.call(rbind,parLapplyLB(cl=cl,X = pars,get_err))
  
  saveRDS(out,file="sweep_ke_ki.RDS")
}

