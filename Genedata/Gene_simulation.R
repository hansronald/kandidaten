# Load packages

#setwd("/Users/rebecka/Documents/Dokument/MSG900/Kod/")
setwd("~/Documents/msg900/Rcode/gene")
#setwd("~/Google Drive/Skola/Kandidatarbete drive/Kod/Simuleringskod")

source("Functions_gene.R")
sapply(c("huge", "expm", "corrplot", "coop", "igraph", "dbscan", "rowr", "pracma", "tictoc", "fields", "xtable"),
       require, character.only = TRUE)

# ---------------------DATA GENERATION/SETUP ---------------------
load("TEs.RData") # 508 observations of 2371 variables
# --------------------- fixed parameters ---------------------

tol = 1e-10
k_seq = c(1,1.5,2)
nmeasures = 10
N = 1 # replicates. for testing use lower N!
j = 1
n_index = 1
N_FDR = integer(length(k_seq))
N_sparsity = integer(length(k_seq))
# load("lambda_mat_gene.R")

ls = seq(1.3,0.6, length.out = 42)
lambda_seq = list(ls^(1.2)*1.3, ls^(1/1.65)/4, (ls)^(1/3)/6.58)


nlambdas = length(lambda_seq[[1]]) # needed!
target_fdr = 0.2
dev_tol = 0.005
t = 0.8

# Create todays results folder
main_dir = getwd()
result_dir <- file.path("resultat", Sys.Date())
dir.create(file.path(main_dir, result_dir), showWarnings = FALSE)

models = c("genedata_prec_subsample500", "genedata_prec", "genedata_prec_subsample1000")
#l = 1
l = 2
#l = 3

  # --------------------- true omega ---------------------
prec = as.matrix(read.table(paste(models[l],".csv",sep = ""),sep=',',row.names = 1, header= TRUE))

adj_plot(adj = adjacency_me(prec), title = "", save_image = TRUE,
         file_name = paste(file.path(main_dir, result_dir), "/", models[l], "_adj_plot_original.png", sep = ""))
  p = ncol(prec)
  #datagen = as.matrix(read.table(paste("TEs_subsample500",".csv",sep = ""),sep=',',row.names = 1, header= TRUE))
  datagen = TEs 
  #datagen = as.matrix(read.table(paste("TEs_subsample1000",".csv",sep = ""),sep=',',row.names = 1, header= TRUE))
  
  sigmahat = cov(datagen)
  n = dim(datagen)[1]
  
    N_FDR = integer(length(k_seq)) 
    N_sparsity = integer(length(k_seq))
    
    prec_sparsity = sparsity(prec)
    adj_org = adjacency_me(prec)
    
  
    
    # -------------------- allocate space --------------
    
    glasso_opt_lambda = matrix(integer(N*length(k_seq)), N, length(k_seq))
    evals = rep(list(matrix(0,N,nmeasures)), length(k_seq))
    TP_evals = rep(list(matrix(0,N,4)), length(k_seq))

    TP_evals_spars = rep(list(matrix(0,N,4)), length(k_seq))
    TP_evals_fdr = rep(list(matrix(0,N,4)), length(k_seq))
    

    fdr_evals = rep(list(matrix(0,N,nmeasures)), length(k_seq))
    sparsity_evals = rep(list(matrix(0,N,nmeasures)), length(k_seq))
    
    sparsity_list = vector("list", length(k_seq))
    lambda_list = vector("list", length(k_seq))
    
    paths = rep(list(matrix(0,2,length(lambda_seq[[1]]))), length(k_seq))
    

      
      i = 1
      for(k in k_seq){
        
        
        # --------------------- empirical k-covariance matrix ---------------------
        E = eigen(sigmahat) 
        E_values = E$values
        E_values[which(abs(E_values) < tol)] = 0 
        D = diag(E_values)^(1/k) 
        P = E$vectors
        sigmahat_k = P %*% D %*% t(P)
        
        plot_scree(E = E, D = D, k = 2)
        # --------------------- glasso ---------------------
        
        tic(paste("huge_obj, k =", k))
        huge_obj = huge(sigmahat_k, method = "glasso",
                       lambda = lambda_seq[[i]], verbose = TRUE)
        toc()
        # --------------------- selecting opt lambda and opt prec ---------------------
        tic("sel_opts")
        
        lambda_list[[i]] = huge_obj$lambda
        sparsity_list[[i]] = get_sparsity_seq(huge_obj = huge_obj, nlambdas = length(ls), k = k, tol = tol)
        
        sel_opts = get_optimal(huge_obj = huge_obj, nlambdas = length(lambda_seq[[i]]), n = n,
                               p = p, sigmahat = sigmahat_k, k = k, tol = tol)
        toc()
        
        glasso_opt_lambda[j,i] = sel_opts[[2]]
        
        # evaluation measures
        evals[[i]][j,] = get_eval(est_prec = as.matrix(sel_opts[[1]]), prec = prec, RUN_ENT = FALSE)
        TP_evals[[i]][j,] = get_TPetc(est_prec = as.matrix(sel_opts[[1]]), prec = prec)
        
        # target sparsity evals
        sparsity_i = get_target_sparsity_index(huge_obj = huge_obj, k = k, tol = tol, nlambdas = length(huge_obj$lambda), target_sparsity = prec_sparsity, dev_tol = dev_tol)
        if(sparsity_i == 0){
          sparsity_evals[[i]][j,] = integer(nmeasures) 
          TP_evals_spars[[i]][j,] = integer(4)
        
        } else{
          sparsity_prec = get_k_glasso(huge_obj = huge_obj, lambda_ind = sparsity_i, k = k, tol = tol)
          sparsity_evals[[i]][j,] = get_eval(est_prec = sparsity_prec,
                                           prec = prec, RUN_ENT = FALSE)
          TP_evals_spars[[i]][j,] = get_TPetc(est_prec = sparsity_prec, prec = prec)
          
          N_sparsity[i] = N_sparsity[i] + 1
        }
        
        # target fdr evals
        fdr_i = get_target_fdr_index(huge_obj = huge_obj, prec = prec, k = k, tol = tol, nlambdas = length(huge_obj$lambda), target_fdr = target_fdr, dev_tol = dev_tol)
        if(fdr_i == 0){
          # fdr_evals[[j,i]] = integer(nmeasures) #TODO change to variable
          fdr_evals[[i]][j,] = integer(nmeasures)
          TP_evals_fdr[[i]][j,] = integer(4)
          
        } else{
          fdr_prec = get_k_glasso(huge_obj = huge_obj, lambda_ind = fdr_i, k = k, tol = tol)
          fdr_evals[[i]][j,] = get_eval(est_prec = fdr_prec,
                                      prec = prec, RUN_ENT = FALSE)
          TP_evals_fdr[[i]][j,] = get_TPetc(est_prec = sparsity_prec, prec = prec)
          
          N_FDR[i] = N_FDR[i] + 1
        }
        
        paths[[i]] = paths[[i]] + get_FS_paths(huge_obj, k, tol, opt_lambda = sel_opts[[2]], plot = FALSE)
        
      
          adj_plot_with_evals(true_mat = adjacency_me(prec), est_mat = adjacency_me(sel_opts[[1]]),
                              model = models[l],
                              file_name = paste(file.path(main_dir, result_dir), "/", models[l],"_n", n, "_k", k, "_adj_plot_evals.png", sep = ""),
                              FN = FALSE, save_image = TRUE, show_title = TRUE)
          if(fdr_i != 0){
            fdr_prec = get_k_glasso(huge_obj = huge_obj, lambda_ind = fdr_i, k = k, tol = tol)
            adj_plot_with_evals(true_mat = adjacency_me(prec), est_mat = adjacency_me(fdr_prec),
                                model = models[l],
                                file_name = paste(file.path(main_dir, result_dir), "/", models[l],"_n", n, "_k", k, "_adj_plot_fdr.png", sep = ""),
                                FN = FALSE, save_image = TRUE, show_title = TRUE)
          }
          if(sparsity_i != 0){
            spars_prec = get_k_glasso(huge_obj = huge_obj, lambda_ind = sparsity_i, k = k, tol = tol)
            adj_plot_with_evals(true_mat = adjacency_me(prec), est_mat = adjacency_me(spars_prec),
                                model = models[l],
                                file_name = paste(file.path(main_dir, result_dir), "/", models[l],"_n", n, "_k", k, "_adj_plot_spars.png", sep = ""),
                                FN = FALSE, save_image = TRUE, show_title = TRUE)
          }
        i = i + 1
      }
    
    
    
    # --------------------- RESULTS ---------------------
    # Plots -------------
      
      plot_lambda_spars(lambda_list, sparsity_list, file = paste(file.path(main_dir, result_dir), "/", models[l],"_n", n, "_lambda_spars.png", sep = ""), save = TRUE)
    
    plot_avg_path(paths, nks = length(k_seq), N = N, save = TRUE, show_title = TRUE,
                  file = paste(file.path(main_dir, result_dir), "/", models[l],"_n", n, "_FDR_paths.png", sep = ""))
    # --------------------- tables --------------------
    sd_lambda = apply(glasso_opt_lambda, 2, sd)
    sd_lambda = matrix(sd_lambda, 1, byrow=TRUE)
    colnames(sd_lambda) = c("k = 1", "k = 1.5", "k = 2")
    rownames(sd_lambda) = "sd_lambda"
    
    aver_eval = get_aver_eval(k_seq = k_seq, evals = evals, N = N)
    aver_eval = aver_eval[,-c(4,5)]
    sd_eval = get_sd_evals(k_seq = k_seq, evals = evals, N = N)
    
    sparsity_aver = get_aver_eval(k_seq = k_seq, evals = sparsity_evals, N = N, N_divide = N_sparsity, t = t)
    sd_sparsity = get_sd_evals(k_seq = k_seq, evals = sparsity_evals, N = N, N_divide = N_sparsity)
    
    fdr_aver = get_aver_eval(k_seq = k_seq, evals = fdr_evals, N = N, N_divide = N_FDR, t = t)
    sd_fdr = get_sd_evals(k_seq = k_seq, evals = fdr_evals, N = N, N_divide = N_FDR, t = t) # 
    
    aver_TP_evals = get_aver_eval(k_seq = k_seq, evals = TP_evals, N = N)

    aver_TP_evals_spars = get_aver_eval(k_seq = k_seq, evals = TP_evals_spars, N_divide = N_sparsity, N = N, t = t)
    aver_TP_evals_fdr = get_aver_eval(k_seq = k_seq, evals = TP_evals_fdr, N = N, N_divide = N_FDR, t = t)
    
    
    txtBIC = xtable(aver_eval,digits = 3)
    txtBIC_sd = xtable(sd_eval,digits = 3)
    txtFDR = xtable(fdr_aver,digits = 3)
    txtspars = xtable(sparsity_aver,digits = 3)
    txtFDR_sd = xtable(sd_fdr,digits = 3)
    txtspars_sd = xtable(sd_sparsity,digits = 3)
    txtlambda_sd = xtable(sd_lambda, digits = 3)
    txtTPetc = xtable(aver_TP_evals, digits = 3)
    txtTPetc_spars = xtable(aver_TP_evals_spars, digits = 3)
    txtTPetc_fdr = xtable(aver_TP_evals_fdr, digits = 3)
    
    
    
    print(txtBIC, file=paste(result_dir, "/", models[l],"_n", n, "-tableBIC.txt", sep = ""))
    print(txtBIC_sd, file=paste(result_dir, "/", models[l],"_n", n, "-tableBIC_sd.txt", sep = ""))
    print(txtFDR, file=paste(result_dir, "/", models[l],"_n", n, "-tableFDR.txt", sep = ""))
    print(txtspars, file=paste(result_dir, "/", models[l],"_n", n, "-tablespars.txt", sep = ""))
    print(txtFDR_sd, file=paste(result_dir, "/", models[l],"_n", n, "-tableFDR_sd.txt", sep = ""))
    print(txtspars_sd, file=paste(result_dir, "/", models[l],"_n", n, "-tablespars_sd.txt", sep = ""))
    print(txtlambda_sd, file=paste(result_dir, "/", models[l],"_n", n, "-tablelambda_sd.txt", sep = ""))
    print(txtTPetc, file=paste(result_dir, "/", models[l],"_n", n, "-tableTPetc.txt", sep = ""))
    print(txtTPetc_spars, file=paste(result_dir, "/", models[l],"_n", n, "-tableTPetc_spars.txt", sep = ""))
    print(txtTPetc_fdr, file=paste(result_dir, "/", models[l],"_n", n, "-tableTPetc_fdr.txt", sep = ""))
    
    print(length(prec[prec != 0]))
    
    

