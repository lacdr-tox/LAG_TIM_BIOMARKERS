## this script performs adaptive metropolis parallel hierarchical sampling 
## as detailed here: Hug, S., et al. "High-dimensional Bayesian parameter 
## estimation: Case study for a model of JAK2/STAT5 signaling." 
## Mathematical biosciences 246.2 (2013): 293-304.
setwd("path/to/LAG_TIM_BIOMARKERS")
set.seed(42)
ncores <- 3

source("code/utils/optim_funcs.R")
library(mvtnorm)
library(parallel)
library(deSolve)


p0 <- readRDS("data/RESULTS/opt_out_wider_bounds.Rds")
v0 <- sapply(p0, function(p0i) p0i$value)

p0 <- p0[order(v0)]
p0 <- p0[1:21]

d0 <- setup_opt()
upper <- d0$df_p0$upper
lower <- d0$df_p0$lower




chains <- lapply(p0, function(p0i) {
  m <- matrix(p0i$par,nrow=1)
  colnames(m) <- parnames
  list(par=m,ll=p0i$value)})

chainswap <- function(chains){
  swapr <- sample(2:length(chains),1)
  tmp <- tail(chains[[1]]$par,1)
  chains[[1]]$par <- rbind(chains[[1]]$par,tail(chains[[swapr]]$par,1))
  chains[[swapr]]$par <- rbind(chains[[swapr]]$par,tmp)
  tmp <- tail(chains[[swapr]]$ll,1)
  chains[[swapr]]$ll <- c(chains[[swapr]]$ll,tail(chains[[1]]$ll,1))
  chains[[1]]$ll <- c(chains[[1]]$ll,tmp)
  list(chains=chains,swapr=swapr)
}



chains <- chainswap(chains)
swapr <- chains$swapr
chains <- chains$chains



ids <- (1:length(chains))[-c(1,swapr)]

px <- lapply(ids, function(i){
  xi <- chains[[i]]
  lli <- Inf
  accept <- FALSE
  ntrials <- 0
  while(!accept){
    p1 <- as.numeric(rmvnorm(1,mean=as.numeric(xi$par),sigma = diag(length(xi$par))*(10^-7)))
    if(prod(log(upper)>=p1)==1 & prod(log(lower)<=p1)==1){ ## bounds check
      #print("TRUE")
      #print(p1)
      #print(lli)
      lli <- get_ll(p1,dat = d0$dat,dat_times=d0$dat_times)
     accept <- runif(1)<exp(xi$ll-lli)
    }
    ntrials <- ntrials+1
    print(ntrials)
  }
  xi$par <- rbind(xi$par,p1)
  xi$ll <- c(xi$ll,lli)
  return(xi)
})

chains <- c(list(chains[[1]]),px,list(chains[[swapr]]))


ma_step <- function(xi,d0 = d0){
  eps <- 0.01
  lli <- Inf
  accept <- FALSE
  ntrials <- 0
  scf <- 0.2
  sigma_mat <- cov(xi$par)*scf
  while(!accept){
    p1 <- as.numeric(rmvnorm(1,mean=as.numeric(xi$par[nrow(xi$par),]),sigma = sigma_mat))
    if(prod((log(d0$df_p0$upper)+eps)>p1)==1 & prod((log(d0$df_p0$lower)-eps)<p1)==1){ ## bounds check
    lli <- get_ll(p1,dat = d0$dat,dat_times=d0$dat_times)
    accept <- runif(1)<exp(tail(xi$ll,1)-lli)
    }
    ntrials <- ntrials+1
    print(ntrials)
  }
  xi$par <- rbind(xi$par,p1)
  xi$ll <- c(xi$ll,lli)
  return(xi)
  
}


ncores <- min(ncores,length(chains)-2)



cl <- makeCluster(getOption("cl.cores", ncores))
clusterCall(cl, function() {
  source("code/utils/optim_funcs.R")
  library(mvtnorm)
  library(deSolve)
})

for(i in 1:100000){
  t0=Sys.time()
  chains <- chainswap(chains)
  swapr <- chains$swapr
  chains <- chains$chains
  ids <- (1:length(chains))[-c(1,swapr)]
  xi <- chains[ids]
  px <- parLapply(cl=cl,X=xi,fun=ma_step,d0=d0)
  #px <- lapply(xi, ma_step,d0=d0 )
  chains <- c(list(chains[[1]]),px,list(chains[[swapr]]))
  t1 <- Sys.time()
  print(paste0("step ",i,": time diff of ",round(t1-t0,digits=3)," seconds"))
  if(i%%250==0) saveRDS(chains,paste0("data/RESULTS/all_pars_wider_chain.Rds"))
}


