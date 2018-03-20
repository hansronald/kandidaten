library(corrplot)
library(ppcor)
library(huge)
library(sna)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(glasso)


#TODO 
# - fixed settings possible to parametrize
#     - number of blocks
#     - size of blocks
#     - n, p, relationship between p and n
#     - add correlation in zero blocks
#     - type of graph for every block
#     - v for every block
#     - zero threshold
# - choose k (in huge.select, or loop k's - choose lambda for every k and then choose k)
# what is the conditional nr for?



n = 50
p = 50
tol = 1e-8
k = 2
lambda_seq = seq(1 / k, 0.05 / k, by = - 0.05 / k)
par(mfrow = c(2,2))
v_par = c(0.01, 5, 10, 20)

#m = 100
#q = 100

# g number of clusters, v parameter to change spread in eigenvalues (larger v larger spread)



k_glasso_eval = matrix(integer(20), 4, 5)
glasso_eval = matrix(integer(20), 4, 5)
rownames(k_glasso_eval) = c("v = 0.1", "v = 5", "v = 10", "v = 20")
rownames(glasso_eval) = c("v = 0.1", "v = 5", "v = 10", "v = 20")
colnames(k_glasso_eval) = c("cond nr", "MCC", "FDR", "lambda", "k")
colnames(glasso_eval) = c("cond nr", "MCC", "FDR", "lambda", "k")


for(ind in 1:length(v_par)){
  #M = matrix(, nrow = 0, ncol = 0)
  prec = matrix(, nrow = 0, ncol = 0)
  for(nblocks in 1:3){
    #n = round(runif(1,50,200))
    #p = round(runif(1,50,200))
    #v = runif(1,0.01,10)
    #g = round(runif(1,1,5))
    
    # g doesnt work on scale-free
    # v gives a higher correlation in the block
    v = nblocks*v_par[ind]
    data = huge.generator(n=n,d=p,graph="scale-free", g = 1, v = v)
    #corrmat = data$sigmahat
    precmat = data$omega
    #mat = matrix(c(1,2,3,4), 2,2)
    #M = as.matrix(bdiag(M,corrmat))
    prec = as.matrix(bdiag(prec,precmat))
  }
  M=solve(prec) #M is covariance matrix to draw data from
  #kappa_M = kappa(M)
  #scree(M)

  
    
    #pheatmap(M, cluster_row = FALSE, cluster_col = FALSE, color=brewer.pal(9,"RdGy"))
    
    
    
    
    
    datagen = mvrnorm(n = n, mu = integer(dim(M)[1]), Sigma = M, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    
    
    sigmahat = cor(datagen)
    scree_obj = scree(sigmahat)
    cond_nr = round(scree_obj[1])
    title(paste("v = ", v_par[ind], ", cond nr = ", cond_nr))

    prec[which(abs(prec) < tol)] = 0 
    
    huge_obj_k = huge(datagen, method = "glasso",
                    lambda = lambda_seq, k = k)
    huge_select_k = huge.select(huge_obj_k)
    huge_prec_k = as.matrix(huge_select_k$opt.icov)
    lambda_sel_k = huge_select_k$opt.lambda
    
    #TODO own lambda?
    huge_obj_g = huge(datagen, method = "glasso",
                      lambda = lambda_seq, k = 1)
    huge_select_g = huge.select(huge_obj_g)
    huge_prec_g = as.matrix(huge_select_g$opt.icov)
    lambda_sel_g = huge_select_g$opt.lambda
    
    glasso_eval[ind,] = c(cond_nr, get_eval(huge_prec_g, prec), lambda_sel_g, 1)
    k_glasso_eval[ind,] = c(cond_nr, get_eval(huge_prec_k, prec), lambda_sel_k, k)

}

glasso_eval
k_glasso_eval







##-------------------------------------------
scree = function(data){
  eig = eigen(data) 
  eig$values[which(abs(eig$values) < 1e-8)] = 0
  #par(mfrow = c(1,1))
  plot(eig$values, type = 'b')
  
  #Only plot non-zero eigenvalues
  #plot(eig$values[which(eig$values > 0)], type = 'b')
  cond = range(eig$values[which(eig$values>0)])[2]/range(eig$values[which(eig$values>0)])[1]
  print("Conditional number ")
  print(cond)
}

get_eval = function(est_prec, prec){
  result = integer(2)
  TP = length(est_prec[est_prec != 0 & prec != 0]) # icke-nollor som ska vara icke-nollor
  TN = length(est_prec[est_prec == 0 & prec == 0]) # nollor som ska vara nollor
  FP = length(est_prec[est_prec != 0 & prec == 0]) # icke-nollor som ska vara nollor
  FN = length(est_prec[est_prec == 0 & prec != 0]) # nollor som ska vara icke-nollor
  SPEC = TN / (TN + FP) # alla rätt gissade nollor/ alla som ska vara nollor
  SENS = TP / (TP + FN) # alla rätt gissade icke-nollor/ alla som ska vara icke-nollor
  result[1] = (TP * TN - FP * FN) / (sqrt(TP + FP) * sqrt(TP + FN) * sqrt(TN + FP) * sqrt(TN + FN)) # matthews correlation coeff.
  result[2] = FP / (FP + TP) # 
  return(result)
}














##---------------------------------------------
## Testplottar
par(mfrow = c(2,2))
colfunc <- colorRampPalette(c("white", "red"))
M[which(abs(M) < 10^-8)] = 0




data2 = huge.generator(n=m,d=q,graph="scale-free", g = 4, v = 5, vis = TRUE)

# Set all close to 0 values to 0 in the precision matrix
data$omega[which(data$omega < 10^-8)] = 0

#m = matrix(c(1,2,3,4),2,2)
#n = matrix(c(5,6,7,8),2,2)
#bdiag(m,n)

scree(cov(data$data))

M = as.matrix(bdiag(data$sigma, data2$sigma))


plot.sociomatrix(M,drawlab=FALSE,diaglab=FALSE)
plot.sociomatrix(data$sigma,drawlab=FALSE,diaglab=FALSE)

plotcor(M,data$sigma)

# Scree plot
# scree(data$sigmahat)

# Correlation and partial correlation
corr = cor(data$data)
omega = data$omega
parcorr = -cov2cor(omega)

# Constructed block matrix
corr = cov2cor(M)
omega = solve(M)
parcorr = -cov2cor(omega)

plotcor(corr,parcorr)

# Pair plot
pairs(data$data)





plotcor = function(C,P){
  par(mfrow = c(1,2))
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  corrplot(C, method="color", col=col(200), main = "Correlation",  
           order="hclust", 
           #addCoef.col = "black", # Add coefficient of correlation
           #tl.col="black", tl.srt=45, #Text label color and rotation
           # Combine with significance
           #p.mat = p.mat, sig.level = 0.01, insig = "blank", 
           # hide correlation coefficient on the principal diagonal
           diag=FALSE,
           mar=c(0,0,1,0),
           tl.pos = "n"
  )
  corrplot(P, method="color", col=col(200), main = "Partial correlation",
           order="hclust", 
           #addCoef.col = "black", # Add coefficient of correlation
           #tl.col="black", tl.srt=45, #Text label color and rotation
           # Combine with significance
           #p.mat = p.mat, sig.level = 0.01, insig = "blank", 
           # hide correlation coefficient on the principal diagonal
           diag=FALSE,
           mar=c(0,0,1,0),
           tl.pos = "n"
  )
}





