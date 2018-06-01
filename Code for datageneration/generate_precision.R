library(huge)
library(igraph)
library(coop)

# fixed parameters
set.seed(1234)
tol = 1e-10
offblock_links_level = 0.0025
nblocks = 4
block_size = 0.25

# loop parameters
p_seq = c(200, 500)
list_i = 1

for(p in p_seq){
    prec_temp = matrix(numeric(0), nrow = 0, ncol = 0)
    for(i in 1:nblocks){
      n = 100
      prec_block = barabasi_me(p*block_size)
      prec_temp = as.matrix(bdiag(prec_temp,prec_block))
    }
    
    # add offblocks links
    add_links = ifelse(matrix(runif(length(prec_temp)),dim(prec_temp)[1], dim(prec_temp)[2]) < offblock_links_level, rnorm(length(prec_temp)), 0)
    add_links[lower.tri(add_links)] = t(add_links)[lower.tri(add_links)]
    
    # -- Set every entry inside a block to zero --------

    start_i = 1
    end_i = p*block_size
    
    for(i in 1:nblocks){
      add_links[seq(start_i, end_i), seq(start_i, end_i)] = 0
      start_i = start_i + p*block_size
      end_i = end_i + p*block_size
    }

    prec_temp = prec_temp + add_links
    prec_temp = prec_temp - (min(eigen(prec_temp)$values)-.00001) * diag(1,p)
    prec_temp[which(abs(prec_temp) < tol)] = 0
    
    prec[[list_i]] = prec_temp 
    list_i = list_i + 1
}

model_p200 = prec[[1]]
model_p500 = prec[[2]]

adj_p200 = adjacency_me(model_p200)
adj_p500 = adjacency_me(model_p500)

par(mfrow = c(1,1))
adj_plot(adj_p200, title = "model_p200")
adj_plot(adj_p500, title = "model_p500")

sparsity(model_p200)
sparsity(model_p500)

setwd("~/Google Drive/Skola/Kandidatarbete drive/Kod/Simuleringskod/precisionsmatriser")
write.table(model_p200, file = "model_p200.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")
write.table(model_p500, file = "model_p500.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")

# Avagyan model with added off diagonal links
set.seed(1234)
data = as.matrix(read.table(paste("Model-4-200",".csv",sep = ""),sep=','))
adj_plot(adjacency_me(data))
prec_temp = data
p = 50
nblocks = 4
block_size = 0.25


add_links = ifelse(   matrix(runif(length(prec_temp)),dim(prec_temp)[1],dim(prec_temp)[2])  # 200*200 matrix with U dist between 0 and 1 
                      < offblock_links_level, sample(data[which(data != 0)], length(prec_temp), replace = TRUE), 0) # if entrys is lower than 0.0025, add_links are rnorm otherwise 0
add_links[lower.tri(add_links)] = t(add_links)[lower.tri(add_links)]

start_i = 1
end_i = p*block_size

for(i in 1:nblocks){
  add_links[seq(start_i, end_i), seq(start_i, end_i)] = 0
  start_i = start_i + p*block_size
  end_i = end_i + p*block_size
}

prec_temp = prec_temp + add_links
prec_temp = prec_temp - (min(eigen(prec_temp)$values)-.001) * eye(p*nblocks)

prec_temp[which(abs(prec_temp) < tol)] = 0
model_avagyan_offlinks = prec_temp
par(mfrow = c(1,1))
adj_plot(adjacency_me(prec))

setwd("~/Google Drive/Skola/Kandidatarbete drive/Kod/precisionsmatriser")
write.table(model_avagyan_offlinks, file = "model_4-200-offlinks.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")

adjacency_me = function(adjacency_matrix){
  adj = adjacency_matrix
  adj[which(adj != 0)]=1
  return(adj)
}

adj_plot = function(adj, title = "Adjacency matrix"){
  par(mar=c(2,0.5,5,0.5))
  image(t(adj)[,nrow(adj):1],
        col=c("white", "midnightblue"),
        xaxt= "n", yaxt= "n", frame.plot=F, main = title)
  par(mar=c(5.1,4.1,4.1,2.1))
}

barabasi_me = function(p){
  prec = matrix(0, nrow = 0, ncol = 0)

  if(p == 50 | p == 125){
    power = runif(1,1,1.4)
    m = runif(1,1,8)
  }else if(p == 1){
    power = runif(1,1,1.5)
    m = runif(1,1,20)
  }
  
  precmat = barabasi.game(n = p, power = power, m = m, out.dist = NULL, out.seq = NULL, out.pref = FALSE, 
                          directed=TRUE)
  precmat = as.matrix(as_adj(precmat))
  t = sum(precmat)
  precmat = ifelse(precmat == 0, 0, rnorm(t, 0, 1))
  precmat[upper.tri(precmat)] = t(precmat)[upper.tri(precmat)]
  prec = as.matrix(bdiag(prec,precmat))

  return(prec)  

}
