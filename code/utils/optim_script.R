## This script performs multi-start optimisation of the model
## with N0 starts using Nc cores

N0 <- 20000
Nc <- 70
setwd("path/to/LAG_TIM_BIOMARKERS")

library(parallel)
source("code/utils/optim_funcs.R")
datapath <- "data/PAR_EST/fit_dat.rds"
outputpath <- "data/RESULTS/opt_out_wider_bounds.Rds"





d0 <- setup_opt(datapath)

wrap_opt <- function(p0,d0){

  dat <- d0$dat
  dat_times <- d0$dat_times
  
  upper=log(d0$df_p0$upper)
  lower=log(d0$df_p0$lower)
  tryCatch(optimr(p0,get_ll,"grcentral", method="L-BFGS-B",control=list(kkt=FALSE,trace=0),
                  upper=upper,lower=lower,dat=dat,dat_times=dat_times), 
           error=function(e) return(list(par=NaN,value=NaN)))
}

cl <- makeCluster(getOption("cl.cores", Nc))
clusterCall(cl, function(nothing_here) {
  library(deSolve)
  library(optimx)
  source("code/utils/optim_funcs.R")
})

eval_p0 <- function(i,d0){
  upper=log(d0$df_p0$upper)
  lower=log(d0$df_p0$lower)
  p0 <- runif(length(upper),min=lower,max=upper)
  ll <- get_ll(p0,dat = d0$dat,dat_times=d0$dat_times)
  list(p0=p0,ll0=ll)
}

r0 <- 1:(10*N0)

p0 <- clusterApplyLB(cl,x=r0,fun=eval_p0,d0=d0)
ll0 <- sapply(p0, function(ppi) ppi$ll0)
p0 <- p0[order(ll0)]
p0 <- p0[1:N0]
ll0 <- sapply(p0, function(ppi) ppi$ll0)
print(ll0)
p0 <- lapply(p0, function(ppi) ppi$p0)


  #opt <- lapply(x, function(xi) wrap_opt(xi,dat,fitmat))
opt <- clusterApplyLB(cl,x=p0,fun=wrap_opt,d0=d0)

saveRDS(opt,outputpath)






