setwd("~/Documents/msg900/Rcode/csvs")

main_dir = getwd()
result_dir <- file.path("resultat", Sys.Date())
dir.create(file.path(main_dir, result_dir), showWarnings = FALSE)

p200n80 = read.csv("p200n80.csv", sep =";", header = FALSE)
p200n80 = edit_mse_fdr(p200n80)
p500n200 = read.csv("p500n200.csv", sep =";", header = FALSE)
p500n200 = edit_mse_fdr(p500n200)

p200n1000= read.csv("p200n1000.csv", sep =";", header = FALSE)
p200n1000 = edit_mse_fdr(p200n1000)
p500n2500 = read.csv("p500n2500.csv", sep =";", header = FALSE)
p500n2500 = edit_mse_fdr(p500n2500)

plot_evals(p200n80, p500n200, file = paste(file.path(main_dir, result_dir), "/", "p200p500lown_spars_evals.png", sep = ""))
plot_evals(p200n1000, p500n2500, file = paste(file.path(main_dir, result_dir), "/", "p200p500highn_spars_evals.png", sep = ""))



p200n80 = read.csv("p200n80BIC.csv", sep =";", header = FALSE)
p200n80 = edit_mse_fdr(p200n80)
p500n200 = read.csv("p500n200BIC.csv", sep =";", header = FALSE)
p500n200 = edit_mse_fdr(p500n200)

p200n1000= read.csv("p200n1000BIC.csv", sep =";", header = FALSE)
p200n1000 = edit_mse_fdr(p200n1000)
p500n2500 = read.csv("p500n2500BIC.csv", sep =";", header = FALSE)
p500n2500 = edit_mse_fdr(p500n2500)

plot_evals(p200n80, p500n200, file = paste(file.path(main_dir, result_dir), "/", "p200p500lown_bic_evals.png", sep = ""))
plot_evals(p200n1000, p500n2500, file = paste(file.path(main_dir, result_dir), "/", "p200p500highn_bic_evals.png", sep = ""))



ava1n80 = read.csv("avagyan200n80.csv", sep =";", header = FALSE)
ava1n80 = edit_mse_fdr(ava1n80)
ava2n80 = read.csv("avagyanofflinksp200n80.csv", sep =";", header = FALSE)
ava2n80 = edit_mse_fdr(ava2n80)

ava1n1000= read.csv("avagyanp200n1000.csv", sep =";", header = FALSE)
ava1n1000 = edit_mse_fdr(ava1n1000)
ava2n1000 = read.csv("avagyanofflinksp200n1000.csv", sep =";", header = FALSE)
ava2n1000 = edit_mse_fdr(ava2n1000)

plot_evals(ava1n80, ava2n80, file = paste(file.path(main_dir, result_dir), "/", "ava1ava2lown_spars_evals.png", sep = ""))
plot_evals(ava1n1000, ava2n1000, file = paste(file.path(main_dir, result_dir), "/", "ava1ava2highn_spars_evals.png", sep = ""))



ava1n80 = read.csv("avagyanp200n80BIC.csv", sep =";", header = FALSE)
ava1n80 = edit_mse_fdr(ava1n80)
ava2n80 = read.csv("avagyanofflinksp200n80BIC.csv", sep =";", header = FALSE)
ava2n80 = edit_mse_fdr(ava2n80)

ava1n1000= read.csv("avagyanp200n1000BIC.csv", sep =";", header = FALSE)
ava1n1000 = edit_mse_fdr(ava1n1000)
ava2n1000 = read.csv("avagyanofflinksp200n1000BIC.csv", sep =";", header = FALSE)
ava2n1000 = edit_mse_fdr(ava2n1000)

plot_evals(ava1n80, ava2n80, file = paste(file.path(main_dir, result_dir), "/", "ava1ava2lown_bic_evals.png", sep = ""))
plot_evals(ava1n1000, ava2n1000, file = paste(file.path(main_dir, result_dir), "/", "ava1ava2highn_bic_evals.png", sep = ""))


# cancerbic = read.csv("cancerBIC.csv", sep =";", header = FALSE)
# cancerbic = edit_mse_fdr(cancerbic)
# cancerspars = read.csv("cancer_sparse.csv", sep =";", header = FALSE)
# cancerspars = edit_mse_fdr(cancerspars)
# 
# plot_evals2(cancerbic, file = paste(file.path(main_dir, result_dir), "/", "cancer_bic_evals.png", sep = ""))
# plot_evals2(cancerspars, file = paste(file.path(main_dir, result_dir), "/", "cancer_spars_evals.png", sep = ""))




plot_evals = function(model1, model2, file, save = TRUE){
  
  if(save == TRUE){
    png(filename=file)
  }
  
  methods = c("SPEC", "SENS", 'PPV','g-centr','s-centr')

  par(mfrow = c(1,1))
  matplot(t(model1), type = 'l', xaxt = 'n', ylim=c(0, 1),
          ylab = "", lty = c(1,1,1),
          col = c('black', 'darkorange', 'deepskyblue2'),
          main = "", cex.lab=2, cex.axis=1.7, cex.main=1.5, cex.sub=1.7)
  matlines (t(model2), type = "l", lty = c(2,2,2),
            col =c('black', 'darkorange', 'deepskyblue2'))
  legend('bottomleft', legend=c('k = 1','k = 1.5',
                         'k = 2'), col=c('black', 'darkorange', 'deepskyblue2'),
         lty=c(1,1,1), cex=1.5, bty = 'n')
  axis(1, at=1:5, labels=methods, cex.axis = 1.4)
  
  if(save == TRUE){
    dev.off()
  }
}


plot_evals2 = function(model1, file, save = TRUE){
  
  if(save == TRUE){
    png(filename=file)
  }
  
  methods = c("SPEC", "SENS", 'FDR','g-centr','s-centr')
  
  par(mfrow = c(1,1))
  matplot(t(model1), type = 'l',  yaxt = 'n', xaxt = 'n', ylim=c(0, 1),
          ylab = "", lty = c(1,1,1),
          col = c('black', 'green', 'red'),
          main = "", cex.lab=2, cex.axis=1.7, cex.main=1.5, cex.sub=1.7)
  legend('bottomleft', legend=c('k = 1','k = 1.5',
                                'k = 2'), col=c('black', 'green', 'red'),
         lty=c(1,1,1), cex=1.5, bty = 'n')
  axis(1, at=1:5, labels=methods, cex.axis = 1.4)
  
  if(save == TRUE){
    dev.off()
  }
}


edit_mse_fdr = function(x){
  # x[,4] = (x[,4] - min(x[,4]))/(max(x[,4]) - min(x[,4]))
  # x[,4] = 1 - x[,4]
  x[,3] = 1 - x[,3]
  x = x[,-c(4,5)]
  return(x)
}
