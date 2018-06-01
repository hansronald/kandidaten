##-------

setwd("~/Google Drive/Skola/Kandidatarbete drive/Kod/Simuleringskod")
for(l in 1:4){
  models = c("model_p200", "model_p500", "Model-4-200", "model_4-200-offlinks_new" )
  prec = as.matrix(read.table(paste("precisionsmatriser", "/", models[l],".csv",sep = ""),sep=','))
  
  p = ncol(prec)
  ns = c(0.4*p, 5*p)
  n = ns[1]
  sigmagen = solve(prec)
  datagen = mvrnorm(n = n, mu = integer(dim(sigmagen)[1]), Sigma = sigmagen, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  sigmahat = cov(datagen)
  
  E = eigen(sigmahat)
  E_values = E$values
  E_values[which(abs(E_values) < tol)] = 0
  D = diag(E_values)^(1/2)
  D2 = diag(E_values)^(1/(1.5))
  P = E$vectors
  
  file_name = paste("bilder/scree_plot_", models[l], ".png", sep = "")
  png(filename=file_name)
  plot_scree2(E = E, Ds = list(D2,D), k = k)
  dev.off()
}

##------

plot_scree2 = function(E, Ds, k){
  par(mfrow = c(1,1))
  par(mar = c(3.1, 4.8, 2.1, 1.5))
  plot(E$values, ylab = "Egenvärde", xlab = "", ylim=c(0, max(E$values)), xlim = c(1,10),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, type = 'l')
  lines(diag(Ds[[1]]), lty = 2)
  lines(diag(Ds[[2]]), lty = 3)
  legend("topright", c("k = 1", "k = 1.5", "k = 2"), lty = c(1,2,3),
         cex=1.5, bty = 'n')
  
  par(mar = c(5.1, 4.1, 4.1, 2.1))
}

