setwd("GENE_EXPRESSION/GEOdata/")
require(data.table)
## note: there are some stats at the head of each text file (now removed - in headers). After these individual rows are prefixed by DATA - colnames row is prefixed by FEATURES...

load_ff <- function(file,treatment,day){
  df <- fread(file)
  df$treatment <- treatment
  df$day <- day
  df
}


keys <- read.csv("key.csv")
setwd("feature_tables")
require(parallel)
  
data <- mclapply(1:nrow(keys), function(i) load_ff(file=keys$filename[i],treatment=keys$treatment[i],day=keys$day[i]),mc.cores=8) 



rename <- function(b){
  x <- data.frame(x=b$gProcessedSignal)
  names(x) <- paste(b$treatment[1],b$day[1],sep="_")
  x
}




bg_check <- do.call(cbind,lapply(data,function(df) df$gIsWellAboveBG))
to_remove <- rowSums(bg_check)==0
print(paste("removing ", sum(to_remove)," probes that were never above BG... "))
bndf <- lapply(data, function(df) df[!to_remove,])
#bndf <- data
df <- bndf[[1]][,c( "FeatureNum", "Row","Col", "ProbeName", "GeneName", "SystematicName")]

dat <- do.call(cbind,lapply(bndf,rename))

df <- cbind(df,dat)


write.csv(df,"../PROCESSED/all_gProcessedSignal_BGfiltered.csv",row.names = FALSE)
