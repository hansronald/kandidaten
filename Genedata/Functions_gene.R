##########################################################################
# Description of the functions                                           #
# Copyright Jenny Andersson, Rebecka Bertilsson, Helena Foogde,          #
# Lovisa Köllerström,  Robin Lindström. Gothenburg University 2018 .     #
##########################################################################

if (exists("find_opt_sparsity") ) rm(find_opt_sparsity)

# --------------------- Functions ---------------------

plot_lambda_spars = function(lambda_list, sparsity_list, file, save = FALSE){
  if(save == TRUE){
    png(filename=file)
  }
  
  par(mfrow = c(1,1))
  
  #plot(seq(1, length(lambda_list[[1]])), sparsity_list[[1]], ylim = c(0.7,1), type = "b", col = "black", 
  #     xlab = "lambda", ylab = "sparsity", main = paste("p", p, "n", n))
  plot(seq(1, length(lambda_list[[1]])), sparsity_list[[1]], ylim = c(0.7,1), type = "b", col = "black", 
       xlab = "lambda index", ylab = "gleshet")
  lines(seq(1, length(lambda_list[[2]])), sparsity_list[[2]], type = "b", col = "red")
  lines(seq(1, length(lambda_list[[3]])), sparsity_list[[3]], type = "b", col = "green")
  legend('bottomright', legend=c("k = 1", "k = 1.5", 'k = 2'), lty=c(1,1,1),
         col = c('black', 'red', 'green'))
  
  if(save == TRUE){
    dev.off()
  }
}



# make_invertible = function(origMat){
#   #cholStatus <- try(u <- chol(origMat), silent = FALSE)
#   detError <- ifelse(det(origMat) == 0, TRUE, FALSE)
#   print(detError)
#   iter <- 0
#   if(detError == FALSE){
#     return(origMat)
#   }
#   newMat <- origMat
#   
#   while (detError) {
#     
#     iter <- iter + 1
#     cat("iteration ", iter, "\n")
#     
#     # replace -ve eigen values with small +ve number
#     newEig <- eigen(newMat)
#     newEig2 <- ifelse(newEig$values < 0, 0, newEig$values)
#     
#     # create modified matrix eqn 5 from Brissette et al 2007, inv = transp for
#     # eig vectors
#     newMat <- newEig$vectors %*% diag(newEig2) %*% t(newEig$vectors)
#     
#     # normalize modified matrix eqn 6 from Brissette et al 2007
#     newMat <- newMat/sqrt(diag(newMat) %*% t(diag(newMat)))
#     
#     # try chol again
#     #cholStatus <- try(u <- chol(newMat), silent = TRUE)
#     detError <- ifelse(det(newMat == 0), TRUE, FALSE)
#   }
#   return(newMat)
# }

get_sparsity_seq = function(huge_obj, nlambdas, k, tol){
  
  sparse_mat = integer(nlambdas)
  
  for(i in seq(1, nlambdas)){
    est_prec = as.matrix(huge_obj$icov[[i]]) %^% k
    est_prec[which(abs(est_prec) < tol)] = 0
    
    sparse_mat[i] = sparsity(est_prec)
  }
  return(sparse_mat)
}


plot_bics_standardize = function(model_bics, save = TRUE, show_title = FALSE, file){
  
  if(save == TRUE){
    png(filename=file)
  }
  
  ymin = min(unlist(lapply(lapply(model_bics,FUN=function(x)apply(x,MARGIN=1,FUN=min)), `[[`, 2)))
  ymax = max(unlist(lapply(lapply(model_bics,FUN=function(x)apply(x,MARGIN=1,FUN=max)), `[[`, 2)))
  
  # normalization
  model_bics[[1]][1,] = (model_bics[[1]][1,] - mean(model_bics[[1]][1,]))/sd(model_bics[[1]][1,])
  model_bics[[1]][2,] = (model_bics[[1]][2,] - mean(model_bics[[1]][2,]))/sd(model_bics[[1]][2,]) 
  model_bics[[1]][3,] = (model_bics[[1]][3,] - mean(model_bics[[1]][3,]))/sd(model_bics[[1]][3,]) 
  
  model_bics[[2]][1,] = (model_bics[[2]][1,] - mean(model_bics[[2]][1,]))/sd(model_bics[[2]][1,])
  model_bics[[2]][2,] = (model_bics[[2]][2,] - mean(model_bics[[2]][2,]))/sd(model_bics[[2]][2,]) 
  model_bics[[2]][3,] = (model_bics[[2]][3,] - mean(model_bics[[2]][3,]))/sd(model_bics[[2]][3,]) 
  
  
  if(is.finite(ymin) && is.finite(ymax)){
    par(mfrow = c(1, dim(model_bics[[1]])[1]))
    for(k_i in 1:dim(model_bics[[1]])[1]){
      if(show_title == TRUE){
        title_text = rownames(model_bics[[1]])[k_i]
      }else{
        title_text = ""
      }
      
      #ymin = min(model_bics[[1]][k_i,],model_bics[[2]][k_i,])
      #ymax = max(model_bics[[1]][k_i,],model_bics[[2]][k_i,])
      
      par(mar = c(5.1,5.1,4.1,2.1))
      
      plot(seq(1, dim(model_bics[[1]])[2]),model_bics[[1]][k_i,], type = 'l',
           xlab = "lambda index", ylab = "BIC",
           main = title_text, 
           cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
      lines(seq(1, dim(model_bics[[2]])[2]), model_bics[[2]][k_i,], type = 'l', col = "red")
      legend('topright', legend=c("p > n", "p < n"), lty=c(1,1), cex=1,3,
             col = c('black', 'red'))
      
      par(mar = c(5.1,4.1,4.1,2.1))
      
    }
    if(save == TRUE){
      dev.off()
    }
  } else {
    print(paste("couldn't print ", file))
  }
}




plot_bics_normalize = function(model_bics, save = TRUE, show_title = FALSE, file){
  
  if(save == TRUE){
    png(filename=file)
  }
  
  ymin = min(unlist(lapply(lapply(model_bics,FUN=function(x)apply(x,MARGIN=1,FUN=min)), `[[`, 2)))
  ymax = max(unlist(lapply(lapply(model_bics,FUN=function(x)apply(x,MARGIN=1,FUN=max)), `[[`, 2)))
  
  # normalization
  model_bics[[1]][1,] = (model_bics[[1]][1,] - min(model_bics[[1]][1,]))/(max(model_bics[[1]][1,]) - min(model_bics[[1]][1,]))
  model_bics[[1]][2,] = (model_bics[[1]][2,] - min(model_bics[[1]][2,]))/(max(model_bics[[1]][2,]) - min(model_bics[[1]][2,]))
  model_bics[[1]][3,] = (model_bics[[1]][3,] - min(model_bics[[1]][3,]))/(max(model_bics[[1]][3,]) - min(model_bics[[1]][3,]))
  
  model_bics[[2]][1,] = (model_bics[[2]][1,] - min(model_bics[[2]][1,]))/(max(model_bics[[2]][1,]) - min(model_bics[[2]][1,]))
  model_bics[[2]][2,] = (model_bics[[2]][2,] - min(model_bics[[2]][2,]))/(max(model_bics[[2]][2,]) - min(model_bics[[2]][2,]))
  model_bics[[2]][3,] = (model_bics[[2]][3,] - min(model_bics[[2]][3,]))/(max(model_bics[[2]][3,]) - min(model_bics[[2]][3,]))
  
  
  if(is.finite(ymin) && is.finite(ymax)){
    par(mfrow = c(1, dim(model_bics[[1]])[1]))
    for(k_i in 1:dim(model_bics[[1]])[1]){
      if(show_title == TRUE){
        title_text = rownames(model_bics[[1]])[k_i]
      }else{
        title_text = ""
      }
     
      #ymin = min(model_bics[[1]][k_i,],model_bics[[2]][k_i,])
      #ymax = max(model_bics[[1]][k_i,],model_bics[[2]][k_i,])
    
      par(mar = c(5.1,5.1,4.1,2.1))
      
      plot(seq(1, dim(model_bics[[1]])[2]),model_bics[[1]][k_i,], type = 'l',
           xlab = "lambda index", ylab = "BIC",
           ylim = c(0,1),
           main = title_text, 
           cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
      lines(seq(1, dim(model_bics[[2]])[2]), model_bics[[2]][k_i,], type = 'l', col = "red")
      legend('topright', legend=c("p > n", "p < n"), lty=c(1,1), cex=1,3,
             col = c('black', 'red'))
      
      par(mar = c(5.1,4.1,4.1,2.1))

    }
    if(save == TRUE){
      dev.off()
    }
  } else {
    print(paste("couldn't print ", file))
  }
}


plot_bics = function(model_bics, save = TRUE, show_title = FALSE, file){

  if(save == TRUE){
    png(filename=file)
  }

  ymin = min(unlist(lapply(lapply(model_bics,FUN=function(x)apply(x,MARGIN=1,FUN=min)), `[[`, 2)))
  ymax = max(unlist(lapply(lapply(model_bics,FUN=function(x)apply(x,MARGIN=1,FUN=max)), `[[`, 2)))


  if(is.finite(ymin) && is.finite(ymax)){
    par(mfrow = c(1, dim(model_bics[[1]])[1]))
    for(k_i in 1:dim(model_bics[[1]])[1]){
      if(show_title == TRUE){
        title_text = rownames(model_bics[[1]])[k_i]
      }else{
        title_text = ""
      }

      ymin = min(model_bics[[1]][k_i,],model_bics[[2]][k_i,])
      ymax = max(model_bics[[1]][k_i,],model_bics[[2]][k_i,])

      par(mar = c(5.1,5.1,4.1,2.1))

      plot(seq(1, dim(model_bics[[1]])[2]),model_bics[[1]][k_i,], type = 'l',
           xlab = "lambda index", ylab = "BIC",
           ylim = c(ymin,ymax),
           main = title_text,
           cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
      lines(seq(1, dim(model_bics[[2]])[2]), model_bics[[2]][k_i,], type = 'l', col = "red")
      legend('topright', legend=c("p > n", "p < n"), lty=c(1,1), cex=1,3,
             col = c('black', 'red'))

      par(mar = c(5.1,4.1,4.1,2.1))

    }
    if(save == TRUE){
      dev.off()
    }
  } else {
    print(paste("couldn't print ", file))
  }
}




get_bic = function(huge_obj, nlambdas, n, p, sigmahat, k, tol){
  
  bic_mat = integer(nlambdas)
  
  for(i in seq(1, nlambdas)){
    est_prec = get_k_glasso(huge_obj = huge_obj, lambda_ind = i, k = k, tol = tol)
    oo = est_prec[lower.tri(est_prec, diag = TRUE)]
    nz = oo[oo != 0]
    nz = length(nz)
    mult = est_prec %*% sigmahat
    trace = sum(diag(mult))
    L = chol(est_prec)
    logdet_estprec = 2*sum(log(diag(L)))
    
    # linear algebra below!
    # det(k*A) = k^p*det(A) (A is p*p)
    # log(det(k*A)) = log(k^p*det(A)) = p*log(k) + log(det(A))
    if(is.finite(det(est_prec))){
      bic_mat[i] = n * (   -log(det(est_prec * 2)) + p * log(2) + trace)    + nz * log(n)
    } else {
      bic_mat[i] = n * (-p*log(2) - logdet_estprec + p * log(2) + trace) + nz * log(n)
    }
    
    #bic_mat[i] = n * (   -log(det(est_prec * 2)) + p * log(2) + trace)    + nz * log(n)
    #bic_mat[i] = n * (-p*log(2) - logdet_estprec + p * log(2) + trace) + nz * log(n)
  }
  
  return(bic_mat)
}


# get_bic = function(huge_obj, nlambdas, n, p, sigmahat, k, tol){
#   
#   bic_mat = integer(nlambdas)
#   
#   for(i in seq(1, nlambdas)){
#     est_prec = get_k_glasso(huge_obj = huge_obj, lambda_ind = i, k = k, tol = tol)
#     oo = est_prec[lower.tri(est_prec, diag = TRUE)]
#     nz = oo[oo != 0]
#     nz = length(nz)
#     mult = est_prec %*% sigmahat
#     trace = sum(diag(mult))
#     bic_mat[i] = n * (   -log(det(est_prec * 2)) + p * log(2) + trace)    + nz * log(n)
#   }
#   
#   return(bic_mat)
# }



get_target_sparsity_index = function(huge_obj, k, tol, nlambdas, target_sparsity, dev_tol){
  sparsitys = integer(nlambdas)
  for(i in seq(1, nlambdas)){
    #est_prec = get_k_glasso(huge_obj = huge_obj, lambda_ind = i, k = k, tol = tol)
    est_prec = as.matrix(huge_obj$icov[[i]]) %^% k
    est_prec[which(abs(est_prec) < tol)] = 0
    sparsitys[i] = sparsity(est_prec)
  }
  if(abs(target_sparsity - sparsitys[which.min(abs(sparsitys - target_sparsity))]) > dev_tol){
    target_index = 0
  } else{
    target_index = which.min(abs(sparsitys - target_sparsity))
  }
  
  return(target_index)
}


get_target_fdr_index = function(huge_obj, prec, k, tol, nlambdas, target_fdr, dev_tol){
  fdr = integer(nlambdas)
  for(i in seq(1, nlambdas)){
    est_prec = get_k_glasso(huge_obj = huge_obj, lambda_ind = i, k = k, tol = tol)
    fdr[i] = get_eval(est_prec = est_prec, prec = prec, RUN_ENT = FALSE)[3]
  }
  if(abs(target_fdr - fdr[which.min(abs(fdr - target_fdr))]) > dev_tol){
    target_index = 0
  } else{
    target_index = which.min(abs(fdr - target_fdr))
  }
  #print(fdr)
  return(target_index)
}

sd_spec = function(vals, N_divide){
  vals_mean = mean(vals)
  #print(vals)
  #print(vals_mean)
  sd_eval = sqrt(sum( (vals - vals_mean)^2 ) / (N_divide - 1))
  #print(sd_eval)
  #print(N_divide)
  ifelse(N_divide <= 1, return(0), return(sd_eval))
}

mean_spec = function(vals, N_divide){
  ifelse(N_divide == 0, 0, sum(vals)/N_divide)
}

get_sd_evals = function(k_seq, evals, N, N_divide = rep(N, length(k_seq)), t = 1){
  
  nks = length(k_seq)
  nmeasures = length(evals[[1]][1,])
  
  sd_eval = matrix(integer(nks*nmeasures), nrow = nks, ncol = nmeasures)
  #list_i = 1
  
  for(j in 1:nks){
    sd_eval[j,] = apply(evals[[j]], 2, sd_spec, N_divide = N_divide[j])
    #for(i in 1:nmeasures){
    #  sd_eval[j,i]=sd(unlist(lapply(evals[list_i :(list_i + N - 1)], "[[", i)))
    #}
    #list_i = list_i + N
    if(N_divide[j]/N < t){
      sd_eval[j,] = rep(NA, nmeasures)
    }
  }
  
  round(sd_eval,3)
  rownames(sd_eval) = paste("k = ", k_seq)
  colnames(sd_eval) = names(evals[[1]])
  
  return(sd_eval)
}


# get_sd_evals = function(k_seq, evals, N, N_divide = rep(N, length(k_seq))){
#   nks = length(k_seq)
#   nmeasures = length(evals[[1,1]])
#   
#   sd_eval = matrix(integer(nks*nmeasures), nrow = nks, ncol = nmeasures)
#   list_i = 1
#   for(j in 1:nks){
#     for(i in 1:nmeasures){
#       sd_eval[j,i]=sd(unlist(lapply(evals[list_i :(list_i + N - 1)], "[[", i)))
#     }
#     list_i = list_i + N
#   }
#   
#   round(sd_eval,3)
#   rownames(sd_eval) = paste("k = ", k_seq)
#   colnames(sd_eval) = names(evals[[1,1]])
#   
#   return(sd_eval)
# }

# #Centrality function and plots without zeros
# print_centrality_2 = function(c_degree, c_degree_g, c_degree_k){
#   topnodes=0.1*p
#   if(is.null(c_degree$res)){
#     ii = which(c_degree$vector >= sort(c_degree$vector, decreasing=T)[topnodes], arr.ind=TRUE)
#     ii_g = which(c_degree_g$vector >= sort(c_degree_g$vector, decreasing=T)[topnodes], arr.ind=TRUE)
#     ii_k = which(c_degree_k$vector >= sort(c_degree_k$vector, decreasing=T)[topnodes], arr.ind=TRUE)
#   } else {
#     ii = which(c_degree$res >= sort(c_degree$res, decreasing=T)[topnodes], arr.ind=TRUE)
#     unordered_values=c_degree$res[ii] #innehåller de värden vi ska jämföra
#     ordered_values = sort(unordered_values, decreasing=TRUE)
#     ordered_indexes=order(unordered_values, decreasing=TRUE)   # innehåller sorterade index på de 
#     # värden vi vill jämföra, index för 
#     # maxvärdet först (index av u_v)
#     
#     ii_g = which(c_degree_g$res >= sort(c_degree_g$res, decreasing=T)[topnodes], arr.ind=TRUE)
#     unordered_values_g=c_degree_g$res[ii_g]
#     ordered_values_g = sort(unordered_values_g, decreasing=TRUE)
#     ordered_indexes_g=order(unordered_values_g, decreasing=TRUE)
#     
#     ii_k = which(c_degree_k$res >= sort(c_degree_k$res, decreasing=T)[topnodes], arr.ind=TRUE)
#     unordered_values_k=c_degree_k$res[ii_k] 
#     ordered_values_k = sort(unordered_values_k, decreasing=TRUE)
#     ordered_indexes_k=order(unordered_values_k, decreasing=TRUE) 
#   }
#   c_degree_max = cbind.fill(ii,ii_g,ii_k, fill = "")
#   colnames(c_degree_max) = c("True omega", "glasso", "k-glasso")
#   print(c_degree_max, row.names = FALSE) # nodes with 10% largest degree centrality values
#   centr = c(c_degree$centralization, c_degree_g$centralization, c_degree_k$centralization) # degree centrality measure
#   names(centr) = c("True omega", "glasso", "k-glasso")
#   print(centr)
#   
#   # How well does the topnodes match?
#   matching_g = na.omit(match(unordered_values_g,unordered_values)) 
#   matching_k = na.omit(match(unordered_values_k,unordered_values))
#   match_vector = c((length(ii)),round(length(matching_g)/length(ii_g),2),round(length(matching_k)/length(ii_k),2))
#   names(match_vector) = c("Number of top 10% nodes","glasso match %","k-glasso match %")
#   print(match_vector)
#   
#   # Matching with one topvalue at a time
#   W_g = matrix()
#   for(index in 1:length(which(ordered_values_g>0))){
#     matching = na.omit(match(unordered_values_g[ordered_indexes_g[1:index]],unordered_values))
#     W_g[index] = length(matching)*100/index
#   }
#   
#   W_k = matrix()
#   for(index in 1:length(which(ordered_values_k>0))){
#     matching = na.omit(match(unordered_values_k[ordered_indexes_k[1:index]],unordered_values))
#     W_k[index] = length(matching)*100/index
#   }
#   
#   plot(W_g,type = "l", xlab = "N.o. topnodes added",ylab = "Percent", xlim = c(0, max(length(ii_g), length(ii_k))), ylim = c(0,100))
#   title("Matching topnodes")
#   lines(W_k,lty = 2)
#   legend('topright', legend = c("Glasso", "k-Glasso"), lty = 1:2)
# }


############ This function is used in the replicate loop
plot_avg_path = function(paths, nks, N, save, show_title = FALSE, file){
  
  if(save == TRUE){
    png(filename=file)
  }
  
  avg_paths = get_avg_path(paths, nks, N)
  
  xmin = min(unlist(lapply(lapply(avg_paths,FUN=function(x)apply(x,MARGIN=1,FUN=min)), `[[`, 1)))
  xmax = max(unlist(lapply(lapply(avg_paths,FUN=function(x)apply(x,MARGIN=1,FUN=max)), `[[`, 1)))
  ymin = min(unlist(lapply(lapply(avg_paths,FUN=function(x)apply(x,MARGIN=1,FUN=min)), `[[`, 2)))
  ymax = max(unlist(lapply(lapply(avg_paths,FUN=function(x)apply(x,MARGIN=1,FUN=max)), `[[`, 2)))
  
  if(show_title == TRUE){
    title_text = paste("N = ", N, "\n", "p = ", p, "\n", "n = ", n, sep = "")
  }else{
    title_text = ""
  }
  
  par(mar=c(5.1,5.1,4.1,2.1))
  plot(avg_paths[[1]][1,], avg_paths[[1]][2,], type = 'l',
       xlab = "FDR", ylab = "SENS",
       #xlim = c(0,1), ylim = c(0,1))
       xlim = c(xmin,xmax), ylim = c(ymin,ymax),
       cex.lab=2, cex.axis=1.7, cex.main=1.7, cex.sub=1.7)
  
  par(mar=c(5.1,4.1,4.1,2.1))
  
  legend_text = "k = 1"
  for(m in 2:nks){
    lines(avg_paths[[m]][1,], avg_paths[[m]][2,], type = 'l', lty = m)
    legend_text = c(legend_text, paste("k = ", k_seq[m], sep = ""))
  }
  
  legend('topleft', legend=legend_text,
         lty=1:nks, cex=1.7, bty = 'n', lwd = rep(1,nks))
  
  string = paste("N = ", N, "\n", "p = ", p, "\n", "n = ", n, sep = "")
  #text(xmax*0.9, xmin*0.9, string, cex = 1)
  
  if(save == TRUE){
    dev.off()
  }
}

get_avg_path = function(paths, nks, N){
  for(nk in 1:nks){
    paths[[nk]] = paths[[nk]] / N
  }
  return(paths)
}

get_scaled_MSE = function(est_prec,prec){
  aa=as.matrix(est_prec-prec)
  bb=as.matrix(0-prec)
  MSE2=norm(aa,"2")^2/norm(bb,"2")
}


get_eval = function(est_prec, prec, RUN_ENT = FALSE){
  result = integer(10)
  names(result) = c("SPEC", "SENS", "FDR", "MCC", "ENT", "MSE","scaled-MSE", "Sparsity","Overlap degree","Overlap between")
  
  TP = length(est_prec[est_prec != 0 & prec != 0])
  TN = length(est_prec[est_prec == 0 & prec == 0])
  FP = length(est_prec[est_prec != 0 & prec == 0])
  FN = length(est_prec[est_prec == 0 & prec != 0])
  spars = sparsity(est_prec)
  
  SPEC = TN / (TN + FP) # alla rätt gissade nollor/ alla som ska vara nollor
  SENS = TP / (TP + FN) # alla rätt gissade icke-nollor/ alla som ska vara icke-nollor
  FDR = FP / (FP + TP)
  MCC = (TP * TN - FP * FN) / (sqrt(TP + FP) * sqrt(TP + FN) * sqrt(TN + FP) * sqrt(TN + FN))
  
  ENT = NaN
  
  if(RUN_ENT){
    # prec needs to be square
    sigma = solve(prec)
    mult = est_prec %*% sigma
    trace = sum(diag(mult))
    ENT = trace - log(det(prec %*% sigma)) - p
  }
  #MSE = (norm(est_prec - prec, "F")) ^ 2
  MSE = get_mse(est_prec,prec)
  scaled_MSE=get_scaled_MSE(est_prec,prec)
  
  
    overlap_degree = centrality_overlap(prec,est_prec)[[1]]
    overlap_between = centrality_overlap(prec,est_prec)[[2]]
 
  
  result[] = c(SPEC, SENS, FDR, MCC, ENT, MSE, scaled_MSE, spars,overlap_degree,overlap_between)
  return(result)
}

get_TPetc = function(est_prec, prec){
  result = integer(4)
  names(result) = c("TP", "TN", "FP", "FN")
  
  TP = length(est_prec[est_prec != 0 & prec != 0])
  TN = length(est_prec[est_prec == 0 & prec == 0])
  FP = length(est_prec[est_prec != 0 & prec == 0])
  FN = length(est_prec[est_prec == 0 & prec != 0])
  
  result[] = c(TP, TN, FP, FN)
  return(result)
}

get_FS_paths = function(huge_obj, k, tol, opt_lambda, plot = FALSE){
  
  lambdas = huge_obj$lambda
  #spar_g = integer(length(lambdas))
  eval_g = matrix(integer(4*length(lambdas)),length(lambdas),4)
  
  for(i in 1:length(lambdas)){
    mat_g = get_k_glasso(huge_obj = huge_obj, lambda_ind = i, k = k, tol = tol)
    eval_g[i,] = get_eval(mat_g, prec, RUN_ENT = FALSE)[1:4]
  }
  
  if(plot == TRUE){
    par(mfrow = c(1,1))
    
    #x_min = min(eval_g[,3], eval_k[,3])
    #x_max = max(eval_g[,3], eval_k[,3])
    #y_min = min(eval_g[,2], eval_k[,2])
    #y_max = max(eval_g[,2], eval_k[,2])
    
    par(mar=c(5.1, 5.1, 4.1, 2.1))
    plot(eval_g[,3], eval_g[,2], type = 'l',
         main = paste("FDR vs SENS (lambda = [", min(lambdas), ",", max(lambdas),"])"),
         xlab = "FDR", ylab = "SENS",
         xlim = c(0,1), ylim = c(0,1), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2
)
    
    points(eval_g[which(huge_obj$lambda == opt_lambda),3],
           eval_g[which(huge_obj$lambda == opt_lambda),2],
           col = 'blue', lwd = 2)
    
    legend('topleft', legend=c("Glasso", "k-Glasso", 'Glasso lambda', "k-Glasso lambda"),
           lty=c(1,2,NA, NA), cex=1.5, bty = 'n', pch = c(NA,NA,1,1), lwd = c(1,1,2,2),
           col = c('black', 'black', 'steelblue', 'tomato'))
  }
  return(rbind(eval_g[,3], eval_g[,2]))
}


plot_scree = function(E, D, k){
  par(mfrow = c(1,1))
  plot(E$values[seq(1,10)], ylab = "Egenvärden", xlab = "", ylim=c(0, max(E$values)), type = 'l')
  lines(diag(D), lty = 2)
  title(main="Egenvärden")
  legend("topright", c("Otransformerade", paste(k, "-roten ur")), lty = c(1,2))
}

adjacency_me = function(adjacency_matrix){
  adj = adjacency_matrix
  adj[which(adj != 0)]=1
  #diag(adj) = 0
  return(adj)
}

get_optimal = function(huge_obj, nlambdas, n, p, sigmahat, k, tol){
  
  bic_mat = get_bic(huge_obj = huge_obj, nlambdas = nlambdas, n = n, p = p, sigmahat =  sigmahat, k = k, tol = tol)
  
  # for(i in seq(1, nlambdas)){
  #   est_prec = get_k_glasso(huge_obj = huge_obj, lambda_ind = i, k = k, tol = tol)
  #   oo = est_prec[lower.tri(est_prec, diag = TRUE)]
  #   nz = oo[oo != 0]
  #   nz = length(nz)
  #   mult = est_prec %*% sigmahat
  #   trace = sum(diag(mult))
  #   bic_mat[i] = n * (-log(det(est_prec * 2)) + p * log(2) + trace) + nz * log(n)
  #   #print(det(as.matrix(huge_obj$cov[[i]])))
  # }
  
  m = min(na.omit(bic_mat))
  m_ind = which(bic_mat == m, arr.ind = TRUE)
  m_ind = m_ind[1] # TODO first index if several matches same min
  #m_ind = m_ind[length(m_ind)] # TODO first index if several matches same min
  
    
  opt_prec = as.matrix(huge_obj$icov[[m_ind]]) %^% k
  opt_prec[which(abs(opt_prec) < tol)] = 0
  opt_lambda = huge_obj$lambda[[m_ind]]
  opts = list(opt_prec, opt_lambda)
  return(opts)
}

# TODO - divide by zero
# get_eval = function(est_prec, prec, RUN_ENT = FALSE){
#   result = integer(7)
#   names(result) = c("SPEC", "SENS", "FDR", "MCC", "ENT", "MSE", "Sparsity")
#   
#   TP = length(est_prec[est_prec != 0 & prec != 0])
#   TN = length(est_prec[est_prec == 0 & prec == 0])
#   FP = length(est_prec[est_prec != 0 & prec == 0])
#   FN = length(est_prec[est_prec == 0 & prec != 0])
#   spars = sparsity(est_prec)
#   
#   SPEC = TN / (TN + FP) # alla rätt gissade nollor/ alla som ska vara nollor
#   SENS = TP / (TP + FN) # alla rätt gissade icke-nollor/ alla som ska vara icke-nollor
#   FDR = FP / (FP + TP)
#   MCC = (TP * TN - FP * FN) / (sqrt(TP + FP) * sqrt(TP + FN) * sqrt(TN + FP) * sqrt(TN + FN))
#   
#   
#   
#   ENT = NaN
#   
#   if(RUN_ENT){
#     # prec needs to be square
#     sigma = solve(prec)
#     mult = est_prec %*% sigma
#     trace = sum(diag(mult))
#     ENT = trace - log(det(prec %*% sigma)) - p
#   }
#   #MSE = (norm(est_prec - prec, "F")) ^ 2
#   MSE = get_mse(est_prec,prec)
#   
#   
#   result[] = c(SPEC, SENS, FDR, MCC, ENT, MSE, spars)
#   return(result)
# }
# 
# get_mse = function(est_prec, prec){
#   return(sum((est_prec - prec) ^2)/(dim(prec)[1]^2))
# }

# divide_N is number of replicates to divide with
# t is threshold for how many valid replicates out of total replicates must exist to produce an average
# leave diide_n and t as default for original get_aver_eval


get_aver_eval = function(k_seq, evals, N, N_divide = rep(N, length(k_seq)), t = 1){
  
  nks = length(k_seq)
  nmeasures = length(evals[[1]][1,])
  
  aver_eval = matrix(integer(nks*nmeasures), nrow = nks, ncol = nmeasures)
  #list_i = 1
  
  for(j in 1:nks){
    aver_eval[j,] = apply(evals[[j]], 2, mean_spec, N_divide = N_divide[j])
    if(N_divide[j]/N < t){
      aver_eval[j,] = rep(NA, nmeasures)
    }
    #for(i in 1:nmeasures){
    #  sd_eval[j,i]=sd(unlist(lapply(evals[list_i :(list_i + N - 1)], "[[", i)))
    #}
    #list_i = list_i + N
  }
  
  round(aver_eval,3)
  rownames(aver_eval) = paste("k = ", k_seq)
  colnames(aver_eval) = names(evals[[1]])
  return(aver_eval)
}

get_aver_bics = function(k_seq, bics, N){
  nks = length(k_seq)
  nbics = length(bics[[1,1]])

  aver_bics = matrix(integer(nks*nbics), nrow = nks, ncol = nbics)
  list_i = 1
  for(j in 1:nks){
    for(i in 1:nbics){
      aver_bics[j,i] = sum(unlist(lapply(bics[list_i :(list_i + N - 1)], "[[", i)))/N
    }
    list_i = list_i + N
  }

  round(aver_bics,3)
  rownames(aver_bics) = paste("k = ", k_seq)
  colnames(aver_bics) = names(bics[[1,1]])
  return(aver_bics)
}

get_matrix_from_blocks = function(M, p, nblocks){
  diag_blocks = list(0)
  off_diag_blocks = list(0)
  
  start_i = 1
  end_i = p/nblocks
  
  for(block in 1:nblocks){
    diag_blocks[[block]] = M[seq(start_i, end_i), seq(start_i, end_i)]
    off_diag_blocks[[block]] = M[seq(start_i, end_i), -seq(start_i, end_i)]      
    
    start_i = start_i + p/nblocks
    end_i = end_i + p/nblocks
  }
  diag_blocks = do.call(cbind, diag_blocks)
  off_diag_blocks = do.call(cbind, off_diag_blocks)
  return(list(diag_blocks, off_diag_blocks))
}

adj_plot = function(adj, title = "Adjacency matrix", save_image = FALSE, file_name){
  
  if(save_image == TRUE){
    png(filename=file_name)
  }
  
  par(mar=c(2,0.5,5,0.5))
  image(t(adj)[,nrow(adj):1],
        col=c("white", "midnightblue"),
        xaxt= "n", yaxt= "n", frame.plot=T, main = title, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
)
  par(mar=c(5.1,4.1,4.1,2.1))
  
  if(save_image == TRUE){
    dev.off()
  }

}

prop_of_FP = function(prec, est_prec, diags, est_diags, off_diags, est_off_diags){
  FP = length(est_prec[est_prec != 0 & prec == 0])
  FP_diags = length(est_diags[est_diags != 0 & diags == 0])
  FP_off_diags = length(est_off_diags[est_off_diags != 0 & off_diags == 0])
  prop_FP_diags = FP_diags/FP
  prop_FP_off_diags = FP_off_diags/FP
  return(c(prop_FP_diags, prop_FP_off_diags))
}


get_k_glasso = function(huge_obj, lambda = NA, lambda_ind = NA, k, tol){
  if(is.na(lambda_ind)){
    lambda_seq_huge = huge_obj$lambda
    lambda_ind = which(lambda_seq_huge==lambda)
  }
  
  prec_k = as.matrix(huge_obj$icov[[lambda_ind]]) %^% k
  prec_k[which(abs(prec_k) < tol)] = 0
  return(prec_k)
}

adj_plot_with_evals = function(true_mat, est_mat, model = "", file_name,FN = FALSE, save_image = FALSE, show_title = TRUE){
  
  col = c("white","olivedrab", "red")
  mat = matrix(0, nrow(true_mat), ncol(true_mat))
  mat_TP = matrix(0, nrow(true_mat), ncol(true_mat))
  mat_FP = matrix(0, nrow(true_mat), ncol(true_mat))
  
  mat_TP = apply(true_mat == 1 & est_mat == 1, c(1,2), as.numeric)  # TP
  TP = sum(mat_TP)

  mat_FP = apply(true_mat == 0 & est_mat == 1, c(1,2), as.numeric)*2   # FP
  FP = sum(mat_FP/2)

  mat = mat_TP + mat_FP

  if (FN == TRUE){
    mat_FN = apply(true_mat == 1 & est_mat == 0, c(1,2), as.numeric)*3      # FN
    mat = mat + mat_FN
    col = c(col,'khaki')
  }
  
  rot_mat = t(mat)[,nrow(mat):1]

  par(mar=c(2,0.5,5,0.5))
  
  if(save_image == TRUE){
    png(filename=file_name)
  }
  
  title_text = ifelse(show_title == FALSE, paste(model, ":", "TP/FP", "(", TP, "/", FP, ")", sep = ""), "")
  image(rot_mat, col = col,
        xaxt= "n", yaxt= "n", frame.plot=T,
        main = title_text)
  
  if(save_image == TRUE){
    dev.off()
  }
  
  par(mar=c(5.1,4.1,4.1,2.1))

}

# find_opt_sparsity = function(huge_obj, huge_obj_k, k, opt_lambda_g, opt_lambda_k, plot_sparse = FALSE){
# 
#   #opt_lambda_g = huge_lambda
#   #opt_lambda_k = huge_lambda_k
#   
#   lambdas_g = huge_obj$lambda
#   lambdas_k = huge_obj_k$lambda
#   
#   spar_g = integer(length(lambdas_g))
#   spar_k = integer(length(lambdas_k))
# 
#   for(i in 1:length(lambda_seq)){
#     mat_g = as.matrix(huge_obj$icov[[i]])
#     mat_k = get_k_glasso(huge_objs_k = huge_obj_k, k_ind = i)
#     #mat_k = as.matrix(huge_obj_k$icov[[i]])%
# 
#     spar_g[i] = sparsity(mat_g)
#     spar_k[i] = sparsity(mat_k)
#   }
#   
#   
#   # Find optimal (BIC) points
#   
#   #lambdas_g = huge_obj$lambda
#   #lambdas_k = huge_obj_k$lambda
#   
#   opt_g_ind = which(lambdas_g==opt_lambda_g)
#   opt_k_ind = which(lambdas_k==opt_lambda_k)
#   
#   mat_opt_g = as.matrix(huge_obj$icov[[opt_g_ind]])
#   mat_opt_k = get_k_glasso(huge_objs_k = huge_obj_k, k_ind = opt_k_ind)
#   #mat_opt_k = as.matrix(huge_obj_k$icov[[opt_k_ind]])%
#   spar_opt_g = sparsity(mat_opt_g)
#   spar_opt_k = sparsity(mat_opt_k)
# 
#   
#   if(plot_sparse == TRUE){
#     #y_min = min(spar_g, spar_k)
#     #y_max = max(spar_g, spar_k)
#     x_min = min(lambdas_g, lambdas_k)
#     plot(lambdas_g, spar_g, type = 'l', xlab = "Lambdas", ylab = "Sparsity",
#          ylim = c(0.8,1), main = "Sparsity for each lambda")
#     lines(lambdas_k, spar_k, lty = 2)
#     points(opt_lambda_g, spar_opt_g, col = "blue")
#     points(opt_lambda_k, spar_opt_k, col = "red")
#     
#     legend('bottomright', legend=c("Glasso", "k-Glasso"),
#            lty=c(1,2), cex=0.8, bty = 'n')
#   }
# 
#   ind_g = which(abs(spar_g-prec_sparsity)==min(abs(spar_g-prec_sparsity)))
#   ind_g = ind_g[1]
#   ind_k = which(abs(spar_k-prec_sparsity)==min(abs(spar_k-prec_sparsity)))
#   ind_k = ind_k[1]
#   
#   
#   sparsity_lambda_g = lambdas_g[ind_g]
#   sparsity_lambda_k = lambdas_k[ind_k]
# 
#   mat = t(matrix(c(ind_g, ind_k,sparsity_lambda_g,sparsity_lambda_k),2,2))
#   colnames(mat) = c("Glasso", "k-Glasso")
#   rownames(mat) = c("Index", "Lambdas")
#   return(mat)
# }

# plot_path = function(huge_obj, huge_obj_k, k, opt_lambda, opt_lambda_k){
# 
#   lambdas = huge_obj$lambda
#   lambdas_k = huge_obj_k$lambda
# 
#   spar_g = integer(length(lambdas))
#   spar_k = integer(length(lambdas_k))
# 
#   eval_g = matrix(integer(4*length(lambda_seq)),length(lambda_seq),4)
#   eval_k = matrix(integer(4*length(lambdas_k)),length(lambdas_k),4)
# 
#   for(i in 1:length(lambda_seq)){
# 
#     mat_g = as.matrix(huge_obj$icov[[i]])
#     mat_k = get_k_glasso(huge_objs_k = huge_obj_k, k_ind = i)
#     #mat_k = as.matrix(huge_obj_k$icov[[i]])
# 
#     eval_g[i,] = get_eval(mat_g, prec, RUN_ENT = FALSE)[1:4]
#     eval_k[i,] = get_eval(mat_k, prec, RUN_ENT = FALSE)[1:4]
# 
#     #spar_k[i] = sparsity(mat_k)
#     #spar_g[i] = sparsity(mat_g)
#   }
# 
#   par(mfrow = c(1,1))
# 
#   x_min = min(eval_g[,3], eval_k[,3])
#   x_max = max(eval_g[,3], eval_k[,3])
#   y_min = min(eval_g[,2], eval_k[,2])
#   y_max = max(eval_g[,2], eval_k[,2])
# 
#   plot(eval_g[,3], eval_g[,2], type = 'l',
#        main = paste("FDR vs SENS (lambda = [", min(lambda_seq), ",", max(lambda_seq),"])"),
#        xlab = "FDR", ylab = "SENS",
#        xlim = c(x_min,x_max), ylim = c(y_min,y_max))
#   lines(eval_k[,3], eval_k[,2], type = 'l', lty = 2)
# 
#   points(eval_g[which(huge_obj$lambda == opt_lambda),3],
#          eval_g[which(huge_obj$lambda == opt_lambda),2],
#          col = 'blue', lwd = 2)
# 
#   points(eval_k[which(huge_obj_k$lambda == opt_lambda_k),3],
#          eval_k[which(huge_obj_k$lambda == opt_lambda_k),2],
#          col = 'red', lwd = 2)
# 
#   legend('topleft', legend=c("Glasso", "k-Glasso", 'Glasso lambda', "k-Glasso lambda"),
#          lty=c(1,2,NA, NA), cex=0.8, bty = 'n', pch = c(NA,NA,1,1), lwd = c(1,1,2,2),
#          col = c('black', 'black', 'steelblue', 'tomato'))
# 
#   #plot(eval_g[,3], eval_g[,4], type = 'l',
#   #     main = paste("FDR vs SPEC (lambda = [", min(lambda_seq), ",", max(lambda_seq),"])"),
#   #     xlab = "FDR", ylab = "SPEC")
#   #lines(eval_k[,3], eval_k[,4], type = 'l', lty = 2)
#   #legend('topright', legend=c("Glasso", "k-Glasso"), lty=1:2, cex=0.8, bty = 'n')
# 
# }

# plot_confusion_sep = function(true_mat,est_mat,lambda){
#   
#   #true_mat = adj_org
#   #est_mat = adj_k
#   
#   par(mfrow = c(2,4))
#   plot.new()
#   par(mar=c(0.1,0,6,0))
#   string = paste("\n\n", "p = ", p, "\n", "n = ", n, "\n", "Sparsity = ",
#                  round(prec_sparsity,4), "\n", "Noise = " ,noise_level)
#   text(0.5, 0.5, string, cex = 1.5)
#   adj_plot(true_mat, title = "Original")
#   par(mar=c(2,0.5,5,0.5))
# 
#   mat = apply(true_mat == 1 & est_mat == 1, c(1,2), as.numeric)  # TP
#   rot_mat = t(mat)[,nrow(mat):1]
#   image(rot_mat, col = c("white","springgreen4"),
#         xaxt= "n", yaxt= "n", frame.plot=T,
#         main = "True positive")
# 
#   mat = apply(true_mat == 0 & est_mat == 1, c(1,2), as.numeric)    # FP
#   rot_mat = t(mat)[,nrow(mat):1]
#   image(rot_mat, col = c("white","mediumorchid"),
#         xaxt= "n", yaxt= "n", frame.plot=T,
#         main = "False positive")
# 
#   plot.new()
#   par(mar=c(0.1,0,6,0))
#   string = paste("Lambda = ", lambda, "\n Sparsity = ", round(sparsity(huge_prec), 4))
#   text(0.5, 0.5, string, cex = 1.5)
#   adj_plot(est_mat, title = "Estimate")
#   par(mar=c(2,0.5,5,0.5))
# 
#   mat = apply(true_mat == 1 & est_mat == 0, c(1,2), as.numeric)      # FN
#   rot_mat = t(mat)[,nrow(mat):1]
#   image(rot_mat,col = c("white","tomato1"),
#         xaxt= "n", yaxt= "n", frame.plot=T,
#         main = "False negative")
#   rot_mat = t(mat)[,nrow(mat):1]
# 
#   mat = apply(true_mat == 0 & est_mat == 0, c(1,2), as.numeric)    # TN
#   rot_mat = t(mat)[,nrow(mat):1]
#   image(rot_mat, col = c("white","steelblue4"),
#         xaxt= "n", yaxt= "n", frame.plot=T,
#         main = "True negative")
# 
# 
#   par(mar=c(5.1,4.1,4.1,2.1))
# 
# }

# plot_icovs = function(adj_org, huge_objs, opt_lambda, from, k){
#   par(mfrow = c(6,9))
#   par(mar=c(0.1,0.1,1,0.1))
#   lambdas = huge_objs$lambda
#   opt = which(lambdas==opt_lambda)
#   #print(opt_lambda)
# 
#   image(t(adj_org)[,nrow(adj_org):1],
#         col=c("white", "midnightblue"),
#         xaxt= "n", yaxt= "n", frame.plot=F, main = "Original")
# 
#   for(i in from:(from+52)){
#     if(k == 1){
#       prec_temp = as.matrix(huge_objs$icov[[i]])
#     }else{
#       #prec_temp = prec_temp
#       prec_temp = get_k_glasso(huge_objs = huge_obj_k, k_ind = i)
#     }
#     adj_temp = adjacency_me(prec_temp)
#     #title = paste("lambda = ", huge_obj_k$lambda[i])
#     col = ifelse(lambdas[i] == opt_lambda, "red", "black")
#     image(t(adj_temp)[,nrow(adj_temp):1],
#           col=c("white", "midnightblue"),
#           xaxt= "n", yaxt= "n", frame.plot=F, main = lambdas[i],
#           col.main=col)
# 
#   }
#   par(mar=c(5.1,4.1,4.1,2.1))
#   par(mfrow = c(1,1))
# }

plot_eval_zoom = function(huge_obj, huge_obj_k, k, opt_lambda_g, opt_lambda_k, zoom = TRUE){
  
  lambda_seq_k = huge_obj_k$lambda
  lambda_seq = huge_obj$lambda
  
  spar_g = integer(length(lambda_seq))
  spar_k = integer(length(lambda_seq_k))
  
  eval_g = matrix(integer(4*length(lambda_seq)),length(lambda_seq),4)
  eval_k = matrix(integer(4*length(lambda_seq_k)),length(lambda_seq_k),4)
  
  for(ind in 1:length(lambda_seq)){
    mat_g = as.matrix(huge_obj$icov[[i]])
    
    #prec_temp = get_k_glasso(huge_objs = huge_obj_k, k_ind = i)
    
    mat_k = get_k_glasso(huge_objs = huge_obj_k, lambda_ind = ind, k = k, tol = tol)
    #mat_k = as.matrix(huge_obj_k$icov[[i]])%^%k
    
    eval_g[ind,] = get_eval(mat_g, prec, RUN_ENT = FALSE)[1:4]
    spar_g[ind] = sparsity(mat_g)
    
    # mat_k = as.matrix(huge_obj_k$icov[[ind]])
    eval_k[ind,] = get_eval(mat_k, prec, RUN_ENT = FALSE)[1:4]
    spar_k[ind] = sparsity(mat_k)
  }
  
  
  # Find optimal (BIC) points
  
  lambdas_g = huge_obj$lambda
  lambdas_k = huge_obj_k$lambda
  
  opt_lambda_g=huge_lambda
  opt_lambda_k=huge_lambda_k
  opt_g_ind = which(lambdas_g==opt_lambda_g)
  opt_k_ind = which(lambdas_k==opt_lambda_k)
  
  mat_opt_g = as.matrix(huge_obj$icov[[opt_g_ind]])
  
  mat_opt_k = get_k_glasso(huge_objs = huge_obj_k, k_ind = opt_k_ind)
  #mat_opt_k = as.matrix(huge_obj_k$icov[[opt_k_ind]])%^%k
  
  eval_opt_g = get_eval(mat_opt_g, prec, RUN_ENT = FALSE)[1:4]
  spar_opt_g = sparsity(mat_opt_g)
  
  # mat_k = as.matrix(huge_obj_k$icov[[i]])
  eval_opt_k = get_eval(mat_opt_k, prec, RUN_ENT = FALSE)[1:4]
  spar_opt_k = sparsity(mat_opt_k)
  
  title = c("SPEC", "SENS", "FDR", "MCC")
  par(mfrow = c(2,2))
  
  for(j in 1:4){
    #y_min = min(eval_g[,j], eval_k[,j])
    #y_max = max(eval_g[,j], eval_k[,j])
    
    plot_range = round(length(lambda_seq)/2):length(lambda_seq)
    plot_range = 1:length(lambda_seq)
    plot_range = min(opt_g_ind, opt_k_ind):length(lambda_seq)
    legend_pos = c("topleft", "bottomleft", "bottomleft", "topleft")
    if(zoom == TRUE){
      xlim =  c(spar_opt_g, spar_opt_k)
      ylim = c(min(eval_opt_g, eval_opt_k), c(max(eval_opt_g, eval_opt_k)))
      
      plot(spar_g[plot_range], eval_g[plot_range,j], type = 'l', main = title[j], xlab = "sparsity", ylab = "",
           xlim = xlim , ylim = ylim)
      lines(spar_k, eval_k[,j], type = 'l', lty = 2)
      
    } else {
      xlim = c(0.95, max(spar_g,spar_k))
      ylim = c(min(na.omit(eval_g[,j]), na.omit(eval_k[,j])), c(max(na.omit(eval_g[,j]), na.omit(eval_k[,j])))) 
      xlim = c(0,1)
      plot(spar_g[plot_range], eval_g[plot_range,j], type = 'l', main = title[j], xlab = "sparsity", ylab = "",
           xlim = xlim, ylim = ylim)
      lines(spar_k, eval_k[,j], type = 'l', lty = 2)
    }
    
    
    
    points(spar_opt_g, eval_opt_g[j], col = 'blue')
    points(spar_opt_k, eval_opt_k[j], col = 'red')
    legend(legend_pos[j], legend=c("Glasso", "k-Glasso"), lty=1:2, cex=0.8, bty = 'n')
  }
  
}


# plot_eval = function(huge_obj, huge_obj_k, k){
# 
#   lambda_seq = huge_obj$lambda
#   lambda_seq_k = huge_obj_k$lambda
#   spar_g = integer(length(lambda_seq))
#   spar_k = integer(length(lambda_seq_k))
# 
#   eval_g = matrix(integer(4*length(lambda_seq)),length(lambda_seq),4)
#   eval_k = matrix(integer(4*length(lambda_seq_k)),length(lambda_seq_k),4)
# 
#   for(i in 1:length(lambda_seq)){
#     mat_g = as.matrix(huge_obj$icov[[i]])
#     mat_k = as.matrix(huge_obj_k$icov[[i]])%^%k
#     eval_g[i,] = get_eval(mat_g, prec, RUN_ENT = FALSE)[1:4]
#     spar_g[i] = sparsity(mat_g)
# 
#     # mat_k = as.matrix(huge_obj_k$icov[[i]])
#     eval_k[i,] = get_eval(mat_k, prec, RUN_ENT = FALSE)[1:4]
#     spar_k[i] = sparsity(mat_k)
#   }
# 
#   title = c("SPEC", "SENS", "FDR", "MCC")
# 
#   par(mfrow = c(1,4))
# 
#   for(j in 1:4){
#     plot(spar_g, eval_g[,j], type = 'l', main = title[j], xlab = "sparsity", ylab = "")
#     lines(spar_k, eval_k[,j], type = 'l', lty = 2)
#     legend('topright', legend=c("Glasso", "k-Glasso"), lty=1:2, cex=0.8, bty = 'n')
#   }
# 
# }

# get_optimal = function(huge_obj, nlambdas, n, p, sigmahat, k){
# 
#   bic_mat = integer(length(lambda_seq))
# 
#   for(i in seq(1, nlambdas)){
#     est_prec = as.matrix(huge_obj$icov[[i]]) %^% k
#     oo = est_prec[lower.tri(est_prec, diag = TRUE)]
#     nz = oo[oo != 0]
#     nz = length(nz)
#     mult = est_prec %*% sigmahat
#     trace = sum(diag(mult))
#     bic_mat[i] = n * (-log(det(est_prec * 2)) + p * log(2) + trace) + nz * log(n)
#   }
# 
#   m = min(na.omit(bic_mat))
#   m_ind = which(bic_mat == m, arr.ind = TRUE)
#   m_ind = m_ind[length(m_ind)] # TODO last index if several matches same min
#   # Which lambda is picked?
# 
#   opt_prec = as.matrix(huge_obj$icov[[m_ind]]) %^% k
#   opt_prec[which(abs(opt_prec) < tol)] = 0
#   opt_lambda = huge_obj$lambda[[m_ind]]
#   opts = list(opt_prec, opt_lambda)
#   return(opts)
# }

# adjacency_me = function(adjacency_matrix){
#   adj = adjacency_matrix
#   adj[which(adj != 0)]=1
#   #diag(adj) = 0
#   return(adj)
# }

# print_eval = function(est_prec, est_prec_k, prec, name_matrix, RUNENT = FALSE){
#   #eval = rbind(get_eval(est_prec = est_prec, prec = prec, RUN_ENT = RUNENT), get_eval(est_prec = est_prec_k, prec = prec, RUN_ENT = RUNENT))
#   eval_g = get_eval(est_prec = est_prec, prec = prec, RUN_ENT = RUNENT)
#   eval_k = get_eval(est_prec = est_prec_k, prec = prec, RUN_ENT = RUNENT)
#   eval = rbind(eval_g, eval_k)
#   eval = round(eval,3)
#   org_sparsity = sparsity(prec)
#   rownames(eval) = c("k = 1", "k = 2")
#   names(dimnames(eval)) <- list("", paste(paste("n = ", n), paste(", p = ", p), paste(", sparsity in org = ", format(round(org_sparsity, 3), nsmall = 3)), paste(", nblocks = ", nblocks), paste(", ", name_matrix)))
#   print(eval)
#   FDR=rbind(round(eval_g[3],3),round(eval_k[3],3))
#   rownames(FDR)=c("k=1", "k=2")
#   return(FDR)
# }

get_mse = function(est_prec, prec){
  return(sum((est_prec - prec) ^2)/(dim(prec)[1]^2))
}

# get_centrality_degree = function(c_obj){
#   adj_obj = adjacency_me(c_obj)
#   g_obj = graph_from_adjacency_matrix(adj_obj, mode = "undirected")
#   c_degree_obj = centr_degree(g_obj)
#   return(c_degree_obj$centralization)
# }
# 
# get_centrality_between = function(c_obj){
#   adj_obj = adjacency_me(c_obj)
#   g_obj = graph_from_adjacency_matrix(adj_obj, mode = "undirected")
#   c_between_obj = centr_betw(g_obj)
#   return(c_between_obj$centralization)
# }
centrality_top = function(c_obj){
  adj_obj = adjacency_me(c_obj)
  g_obj = graph_from_adjacency_matrix(adj_obj, mode = "undirected")
  degree_obj = centr_degree(g_obj)
  between_obj = centr_betw(g_obj)
  
  topnodes=0.1*p
  if(is.null(degree_obj$res)){
    ii_d = which(degree_obj$vector >= sort(degree_obj$vector, decreasing=T)[topnodes], arr.ind=TRUE)
  } else {
    ii_d = which(degree_obj$res >= sort(degree_obj$res, decreasing=T)[topnodes], arr.ind=TRUE)
    #unordered_values_d = degree_obj$res[ii_d] #innehåller de värden vi ska jämföra
    #ordered_values_d = sort(unordered_values_d, decreasing=TRUE)
    #ordered_indexes_d = order(unordered_values_d, decreasing=TRUE)   # innehåller sorterade index på de 
  }
  if(is.null(between_obj$res)){
    ii_b = which(between_obj$vector >= sort(between_obj$vector, decreasing=T)[topnodes], arr.ind=TRUE)
  } else {
    ii_b = which(between_obj$res >= sort(between_obj$res, decreasing=T)[topnodes], arr.ind=TRUE)
    #unordered_values_b = between_obj$res[ii_b] #innehåller de värden vi ska jämföra
    #ordered_values_b = sort(unordered_values_b, decreasing=TRUE)
    #ordered_indexes_b = order(unordered_values_b, decreasing=TRUE)   # innehåller sorterade index på de 
  } 
  centrality_list = list(ii_d,ii_b)
  return(centrality_list)
}
centrality_overlap = function(prec,est){
  matching_degree = na.omit(match(centrality_top(est)[[1]],centrality_top(prec)[[1]])) 
  matching_between = na.omit(match(centrality_top(est)[[2]],centrality_top(prec)[[2]])) 
  match_vector = c(round(length(matching_degree)/length(centrality_top(est)[[1]]),2),round(length(matching_between)/length(centrality_top(est)[[2]]),2))
  names(match_vector) = c("match % degree", "match % between")
  return(match_vector)
}



