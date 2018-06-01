

source("Functions.R")

# Files TEs.RData and tfnet.RData has to be requested from robin.lindstrom@gmail.com
load("genedata/TEs.RData")
write.csv(TEs, file = "TEs.csv")
load("genedata/tfnet.RData")
sapply(c("huge", "expm", "corrplot", "coop", "igraph", "dbscan", "rowr", "pracma", "tictoc", "fields", "xtable"),
       require, character.only = TRUE)

genedat_cov = cov(TEs)
write.csv(rownames(genedat_cov), file = "rownames_genedat_cov.csv")

nlinks = nrow(tfnet)
genedat_prec = diag(1,nrow(genedat_cov))

colnames(genedat_prec) = colnames(genedat_cov)
rownames(genedat_prec) = rownames(genedat_cov)
nlinks_actual = 0

for(i in 1:nlinks){
  a = as.character(tfnet[i,1])
  b = as.character(tfnet[i,2])
  if(a == b){
  
  }else if(tryCatch(genedat_prec[a,], error=function(e) return(TRUE)) || tryCatch(genedat_prec[,b], error=function(e) return(TRUE))){
    
  }else{
    genedat_prec[a,b] = 1
    genedat_prec[b,a] = 1
    nlinks_actual = nlinks_actual + 1
  }
}
print(nlinks_actual)

write.csv(genedat_prec, file = "genedata_prec.csv")
adj_plot(genedat_prec)

# Test to plot the written data
#prec = as.matrix(read.table(paste("Genedata", "/", "genedata_prec",".csv",sep = ""),sep=',',row.names = 1, header= TRUE))
#adj_plot(prec)

# Create a subsample

set.seed(1234)
ii = sample(ncol(TEs),500)
TEs2 = TEs[,ii]
write.csv(TEs2, file = "TEs_subsample500.csv")
genedat_cov = cov(TEs2)


nlinks = nrow(tfnet)
genedat_prec = diag(1,nrow(genedat_cov))

colnames(genedat_prec) = colnames(genedat_cov)
rownames(genedat_prec) = rownames(genedat_cov)

for(i in 1:nlinks){
  a = as.character(tfnet[i,1])
  b = as.character(tfnet[i,2])
  
  if(tryCatch(genedat_prec[a,], error=function(e) return(TRUE)) || tryCatch(genedat_prec[,b], error=function(e) return(TRUE))){
    print(i)
  }else{
    genedat_prec[a,b] = 1
    genedat_prec[b,a] = 1
  }

}
adj_plot(genedat_prec)
write.csv(genedat_prec, file = "genedata_prec_subsample500.csv")


