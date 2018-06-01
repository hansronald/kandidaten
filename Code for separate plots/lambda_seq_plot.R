# ----------
# for every model: lambda sequence should cover true sparsity +- x %
#   for every k, lambda[i]_k should respond to lambda[i]_k for all k

library(coop)
library(MASS)

k_seq = c(1, 1.5, 2)
tol = 1e-10

sparsity_list = vector("list", length(k_seq))

lambda_list = vector("list", length(k_seq))

# model = c("model_p200") 
# model = c("model_p1000")
 model = c("Model-4-200")
# model = c("model_4-200-offlinks")


# n = 80
# n = 400
 n = 1000
# n = 5000

prec = as.matrix(read.table(paste(model,".csv",sep = ""),sep=','))
p = dim(prec)[1]
prec_s = sparsity(prec)

sigmagen = solve(prec)
datagen = mvrnorm(n = n, mu = integer(dim(sigmagen)[1]), Sigma = sigmagen, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
sigmahat = cov(datagen)



# p200n80 OK
# ls = seq(0.187, 0.045, length.out = 48)
# lambda_seq = matrix(integer(3*length(ls)), nrow = 3)
# lambda_seq[1,] = ls
# lambda_seq[2,] = ls
# lambda_seq[3,] = (ls)^(1/1.4)/1.5


# p200n1000 OK
# ls = seq(0.105, 0.014, length.out = 48)
# lambda_seq = matrix(integer(3*length(ls)), nrow = 3)
# lambda_seq[1,] = ls
# lambda_seq[2,] = ls
# lambda_seq[3,] = (ls)^(1/1.2)*1.5

# # p1000n400
# ls = seq(0.02, 0.004, length.out = 48)
# lambda_seq = matrix(integer(3*length(ls)), nrow = 3)
# lambda_seq[1,] = ls
# lambda_seq[2,] = ls^(1/1.1)*1.1
# lambda_seq[3,] = (ls)*4.5

# # p1000n5000 kolla!
#ls = seq(0.007, 0.0025, length.out = 48)
#lambda_seq = matrix(integer(3*length(ls)), nrow = 3)
#lambda_seq[1,] = ls
#lambda_seq[2,] = (ls)^(1/1.1)*1.05
#lambda_seq[3,] = (ls)

#avagyann80 OK
# ls = seq(0.635, 0.40, length.out = 48)
# lambda_seq = matrix(integer(3*length(ls)), nrow = 3)
# lambda_seq[1,] = ls*3.7
# lambda_seq[2,] = ls
# lambda_seq[3,] = (ls)^(1/2.4)/1.5

#avagyann80 OK
# ls = seq(0.635, 0.40, length.out = 48)
# lambda_seq = matrix(integer(3*length(ls)), nrow = 3)
# lambda_seq[1,] = ls*3.7
# lambda_seq[2,] = ls
# lambda_seq[3,] = (ls)^(1/2.4)/1.5


# # avagyann1000  OK
# ls = seq(0.2, 0.13, length.out = 48)
# lambda_seq = matrix(integer(3*length(ls)), nrow = 3)
# lambda_seq[1,] = (ls)*3.7
# lambda_seq[2,] = ls
# lambda_seq[3,] = ((ls)^(1/1.5))/1.8


# Avagyans k
lambda_seq[1,] =  seq(1/1,0.05/1,length.out = 48)
lambda_seq[2,] =  seq(1/(1.5),0.05/(1.5),length.out = 48)
lambda_seq[3,] =  seq(1/2,0.05/2,length.out = 48)



# #avagyanoffn80 OK
# ls = seq(0.117, 0.046, length.out = 48)
# lambda_seq = matrix(integer(3*length(ls)), nrow = 3)
# lambda_seq[1,] = ls
# lambda_seq[2,] = ls/1.1
# lambda_seq[3,] = (ls)^(1/3.1)/5



# avagyanoffn1000 OK

# ls = seq(0.023, 0.016, length.out = 48)
# lambda_seq = matrix(integer(3*length(ls)), nrow = 3)
# lambda_seq[1,] = (ls)
# lambda_seq[2,] = ls/1.1
# lambda_seq[3,] = ((ls)^(1/2.3))/7.52



tic()
i = 1
for (k in k_seq){
  E = eigen(sigmahat) 
  E_values = E$values
  E_values[which(abs(E_values) < tol)] = 0 
  D = diag(E_values)^(1/k) 
  P = E$vectors
  sigmahat_k = P %*% D %*% t(P)


  huge_obj =  huge(sigmahat_k, method = "glasso",
                       lambda = lambda_seq[i,], verbose = FALSE)
  
  lambda_list[[i]] = huge_obj$lambda
  sparsity_list[[i]] = get_sparsity_seq(huge_obj = huge_obj, nlambdas = length(lambda_seq[i,]), k = k, tol = tol)
  
  
  i = i + 1
}

#par(mfrow = c(1,4))
par(mfrow = c(1,1))

par(mar = c(5.1, 4.5, 4.1, 1.7))
plot(seq(1, length(lambda_list[[1]])), sparsity_list[[1]], ylim = c(0,1), type = "b", col = "black", 
     xlab = "lambda", ylab = "sparsity",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
lines(seq(1, length(lambda_list[[2]])), sparsity_list[[2]], type = "b", col = "red")
lines(seq(1, length(lambda_list[[3]])), sparsity_list[[3]], type = "b", col = "green")
legend('bottomleft', legend=c("k = 1", "k = 1.5", 'k = 2'), lty=c(1,1,1),
       col = c('black', 'red', 'green'),
       cex=1.5, bty = 'n')
par(mar = c(5.1, 4.1, 4.1, 2.1))


print(range(sparsity_list[[1]]))
print(range(sparsity_list[[2]]))
print(range(sparsity_list[[3]]))


toc()

# --------- functions ----------
get_sparsity_seq = function(huge_obj, nlambdas, k, tol){
  
  sparse_mat = integer(nlambdas)
  
  for(i in seq(1, nlambdas)){
    est_prec = as.matrix(huge_obj$icov[[i]]) %^% k
    est_prec[which(abs(est_prec) < tol)] = 0
    
    sparse_mat[i] = sparsity(est_prec)
  }
  return(sparse_mat)
}

make_log_scale = function(mi, ma, steps){
  b = log(mi/ma)/(mi-ma)
  a = mi/exp(b*mi)
  return(rev(a*exp(b*seq(mi, ma, by = steps))))
}
