setwd("/Users/rebecka/Documents/Dokument/MSG900/Kod/precisionsmatriser")

prec = as.matrix(read.table("model_p200.csv", sep=','))
cov200 = solve(prec)
n = 80
datagen = mvrnorm(n = n, mu = integer(dim(cov200)[1]), Sigma = cov200, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
cov200 = cov(datagen)

prec = as.matrix(read.table("model_p500.csv",sep=','))
cov500 = solve(prec)
n = 200
datagen = mvrnorm(n = n, mu = integer(dim(cov500)[1]), Sigma = cov500, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
cov500 = cov(datagen)

prec = as.matrix(read.table("Model-4-200.csv",sep=','))
cov4 = solve(prec)
n = 80
datagen = mvrnorm(n = n, mu = integer(dim(cov4)[1]), Sigma = cov4, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
cov4 = cov(datagen)

prec = as.matrix(read.table("model_4-200-offlinks_new.csv",sep=','))
cov4plus = solve(prec)
n = 80
datagen = mvrnorm(n = n, mu = integer(dim(cov4plus)[1]), Sigma = cov4plus, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
cov4plus = cov(datagen)

tol = 1e-10

par(mar=c(5.1, 5.1, 4.1, 2.1))
E = eigen(cov200) 
E_values = E$values
E_values[which(abs(E_values) < tol)] = 0 
D1 = diag(E_values)
D15 = diag(E_values)^(1/1.5)
D2 = diag(E_values)^(1/2)
P = E$vectors
cov200 = P %*% D1 %*% t(P)
cov200_15 = P %*% D15 %*% t(P)
cov200_2 = P %*% D2 %*% t(P)
range(E_values)


# kappa200 = signif(kappa(cov200),2)
# kappa200_15 = signif(kappa(cov200_15),2)
# kappa200_2 = signif(kappa(cov200_2),2)

kappa200 = signif(max(E_values)/min(E_values[which(E_values != 0)]),2)
kappa200_15 = signif((max(E_values)^(1/1.5))/(min(E_values[which(E_values != 0)])^(1/1.5)),2)
kappa200_2 = signif((max(E_values)^(1/2))/(min(E_values[which(E_values != 0)])^(1/2)),2)

plot(E$values[seq(1,5)], ylab = "Egenvärden", xlab = "Index för egenvärden", ylim=c(0, max(E$values)), type = 'b', cex.lab=1.8, cex.axis=1.8)
lines(diag(D15), type="b", lty = 2)
lines(diag(D2),type="b", lty = 3)
legend("topright", c(paste("k = 1 (",kappa200,")",sep=""), paste("k = 1.5 (", kappa200_15, ")",sep=""), paste("k = 2 (", kappa200_2, ")",sep="")), lty = c(1,2,3), cex=1.8)



E = eigen(cov500) 
E_values = E$values
E_values[which(abs(E_values) < tol)] = 0 
D1 = diag(E_values)
D15 = diag(E_values)^(1/1.5)
D2 = diag(E_values)^(1/2)
P = E$vectors
cov500 = P %*% D1 %*% t(P)
cov500_15 = P %*% D15 %*% t(P)
cov500_2 = P %*% D2 %*% t(P)
range(E_values)


 kappa500 = signif(kappa(cov500),2)
 kappa500_15 = signif(kappa(cov500_15),2)
 kappa500_2 = signif(kappa(cov500_2),2)

kappa500 = signif(max(E_values)/min(E_values[which(E_values != 0)]),2)
kappa500_15 = signif((max(E_values)^(1/1.5))/(min(E_values[which(E_values != 0)])^(1/1.5)),2)
kappa500_2 = signif((max(E_values)^(1/2))/(min(E_values[which(E_values != 0)])^(1/2)),2)

plot(E$values[seq(1,5)], ylab = "Egenvärden", xlab = "Index för egenvärden", ylim=c(0, max(E$values)), type = 'b', cex.lab=1.8, cex.axis=1.8)
lines(diag(D15), type="b", lty = 2)
lines(diag(D2), type="b",lty = 3)
legend("topright", c(paste("k = 1 (", kappa500, ")",sep=""), paste("k = 1.5 (", kappa500_15, ")",sep=""), paste("k = 2 (", kappa500_2, ")",sep="")), lty = c(1,2,3), cex=1.8)




E = eigen(cov4) 
E_values = E$values
E_values[which(abs(E_values) < tol)] = 0 
D1 = diag(E_values)
D15 = diag(E_values)^(1/1.5)
D2 = diag(E_values)^(1/2)
P = E$vectors
cov4 = P %*% D1 %*% t(P)
cov4_15 = P %*% D15 %*% t(P)
cov4_2 = P %*% D2 %*% t(P)
range(E_values)


# kappa4 = signif(kappa(cov4),2)
# kappa4_15 = signif(kappa(cov4_15),2)
# kappa4_2 = signif(kappa(cov4_2),2)

kappa4 = signif(max(E_values)/min(E_values[which(E_values != 0)]),2)
kappa4_15 = signif((max(E_values)^(1/1.5))/(min(E_values[which(E_values != 0)])^(1/1.5)),2)
kappa4_2 = signif((max(E_values)^(1/2))/(min(E_values[which(E_values != 0)])^(1/2)),2)

plot(E$values[seq(1,5)], ylab = "Egenvärden", xlab = "Index för egenvärden", ylim=c(0, max(E$values)), type = 'b', cex.lab=1.8, cex.axis=1.8)
lines(diag(D15), type="b", lty = 2)
lines(diag(D2), type="b",lty = 3)
legend("topright", c(paste("k = 1 (", kappa4, ")",sep=""), paste("k = 1.5 (", kappa4_15, ")",sep=""), paste("k = 2 (", kappa4_2, ")",sep="")), lty = c(1,2,3), cex=1.8)




E = eigen(cov4plus) 
E_values = E$values
E_values[which(abs(E_values) < tol)] = 0 
D1 = diag(E_values)
D15 = diag(E_values)^(1/1.5)
D2 = diag(E_values)^(1/2)
P = E$vectors
cov4plus = P %*% D1 %*% t(P)
cov4plus_15 = P %*% D15 %*% t(P)
cov4plus_2 = P %*% D2 %*% t(P)
range(E_values)


# kappa4plus = signif(kappa(cov4plus),2)
# kappa4plus_15 = signif(kappa(cov4plus_15),2)
# kappa4plus_2 = signif(kappa(cov4plus_2),2)

kappa4plus = signif(max(E_values)/min(E_values[which(E_values != 0)]),2)
kappa4plus_15 = signif((max(E_values)^(1/1.5))/(min(E_values[which(E_values != 0)])^(1/1.5)),2)
kappa4plus_2 = signif((max(E_values)^(1/2))/(min(E_values[which(E_values != 0)])^(1/2)),2)

plot(E$values[seq(1,5)], ylab = "Egenvärden", xlab = "Index för egenvärden", ylim=c(0, max(E$values)), type = 'b', cex.lab=1.8, cex.axis=1.8)
lines(diag(D15), type="b", lty = 2)
lines(diag(D2), type="b",lty = 3)
legend("topright", c(paste("k = 1 (", kappa4plus, ")",sep=""), paste("k = 1.5 (", kappa4plus_15, ")",sep=""), paste("k = 2 (", kappa4plus_2, ")",sep="")), lty = c(1,2,3), cex=1.8)


#----------
# par(mar=c(5.1, 5.1, 4.1, 2.1))
# E = eigen(sigmahat) 
# E_values = E$values
# E_values[which(abs(E_values) < tol)] = 0 
# D1 = diag(E_values)
# D15 = diag(E_values)^(1/1.5)
# D2 = diag(E_values)^(1/2)
# P = E$vectors
# sigmahat = P %*% D1 %*% t(P)
# sigmahat_15 = P %*% D15 %*% t(P)
# sigmahat_2 = P %*% D2 %*% t(P)
# range(E_values)
# 
# 
# # kappa200 = signif(kappa(cov200),2)
# # kappa200_15 = signif(kappa(cov200_15),2)
# # kappa200_2 = signif(kappa(cov200_2),2)
# 
# kappasigmahat = signif(max(E_values)/min(E_values[which(E_values != 0)]),2)
# kappasigmahat_15 = signif((max(E_values)^(1/1.5))/(min(E_values[which(E_values != 0)])^(1/1.5)),2)
# kappasigmahat_2 = signif((max(E_values)^(1/2))/(min(E_values[which(E_values != 0)])^(1/2)),2)
# 
plot(E$values[seq(1,5)], ylab = "Egenvärden", xlab = "Index för egenvärden", ylim=c(0, max(E$values)), type = 'b', cex.lab=1.8, cex.axis=1.8)
lines(diag(D15),type="b", lty = 2)
lines(diag(D2),type="b", lty = 3)
legend("topright", c(paste("k = 1 (", kappasigmahat, ")",sep=""), paste("k = 1.5 (", kappasigmahat_15, ")",sep=""), paste("k = 2 (", kappasigmahat_2, ")",sep="")), lty = c(1,2,3), cex=1.8)
# 
# 
