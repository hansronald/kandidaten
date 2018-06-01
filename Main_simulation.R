# Load packages

source("Functions.R")
sapply(c("huge", "expm", "corrplot", "coop", "igraph", "dbscan", "rowr", "pracma", "tictoc", "fields", "xtable"),
       require, character.only = TRUE)

# ---------------------DATA GENERATION/SETUP ---------------------
# --------------------- fixed parameters ---------------------

tol = 1e-10
nblocks = 4
k_seq = c(1,1.5,2)
nmeasures = 10
N = 2 # replicates. for testing use lower N!
N_FDR = integer(length(k_seq))
N_sparsity = integer(length(k_seq))
load("lambda_mat.R")
nlambdas = length(lambda_mat[[1,1,1]]) # needed!
target_fdr = 0.2
dev_tol = 0.02
t = 0.8

# Create todays results folder
main_dir = getwd()
result_dir <- file.path("resultat", Sys.Date())
dir.create(file.path(main_dir, result_dir), showWarnings = FALSE)

models = c("model_p200", "model_p500", "Model-4-200", "model_4-200-offlinks_new" ) # do not change order!!!!

tic("Total time")
for(l in 1:length(models)){

  tic(models[l])
  # --------------------- true omega ---------------------
  prec = as.matrix(read.table(paste("precisionsmatriser", "/", models[l],".csv",sep = ""),sep=','))

  par(mfrow = c(1,1))
  adj_plot(adj = adjacency_me(prec), title = "", save_image = FALSE,
           file_name = paste(file.path(main_dir, result_dir), "/", models[l], "_adj_plot_original.png", sep = ""))
  
  adj_plot_with_evals(true_mat = adjacency_me(prec), est_mat = adjacency_me(sel_opts[[1]]),
                      model = models[l],
                      file_name = paste(file.path(main_dir, result_dir), "/", models[l],"_n", n, "_k", k, "_adj_plot_evals.png", sep = ""),
                      FN = FALSE, save_image = FALSE, show_title = TRUE)

  
  p = ncol(prec)
  ns = c(0.4*p, 5*p)
  model_bics = lapply(1:length(ns), matrix, data = NA, nrow = length(k_seq), ncol = nlambdas)
  
  n_index = 1
  for(n in ns){
    N_FDR = integer(length(k_seq)) 
    N_sparsity = integer(length(k_seq))
    
    prec_sparsity = sparsity(prec)
    adj_org = adjacency_me(prec)
    
    d_temp = get_matrix_from_blocks(M = prec, p = p, nblocks = nblocks) 
    diags = d_temp[[1]]
    off_diags = d_temp[[2]]
    
    # -------------------- allocate space --------------
    
    glasso_opt_lambda = matrix(integer(N*length(k_seq)), N, length(k_seq))
  
    evals = rep(list(matrix(0,N,nmeasures)), length(k_seq))
    
    diag_evals = rep(list(matrix(0,N,nmeasures)), length(k_seq))
    off_diag_evals = rep(list(matrix(0,N,nmeasures)), length(k_seq))
    
    FP_evals =  rep(list(matrix(0,N,2)), length(k_seq))
    FP_spars =  rep(list(matrix(0,N,2)), length(k_seq))
    FP_fdr =  rep(list(matrix(0,N,2)), length(k_seq))
    
  
    fdr_evals = rep(list(matrix(0,N,nmeasures)), length(k_seq))
    sparsity_evals = rep(list(matrix(0,N,nmeasures)), length(k_seq))
    
    bics = matrix(list(), N, length(k_seq))
    paths = rep(list(matrix(0,2,length(lambda_mat[[l,n_index, 1]]))), length(k_seq))
    
    experiment_sample = sample(N,1)
    
    for(j in 1:N){
      # --------------------- generate data ---------------------
      
      sigmagen = solve(prec)
      datagen = mvrnorm(n = n, mu = integer(dim(sigmagen)[1]), Sigma = sigmagen, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
      sigmahat = cov(datagen)
      
      i = 1
      for(k in k_seq){
        
        
        # --------------------- empirical k-covariance matrix ---------------------
        E = eigen(sigmahat)
        E_values = E$values
        E_values[which(abs(E_values) < tol)] = 0
        D = diag(E_values)^(1/k)
        P = E$vectors
        sigmahat_k = P %*% D %*% t(P)
        
        plot_scree(E = E, D = D, k = k)
        # --------------------- glasso ---------------------
        
        tic(paste("huge_obj, k =", k))
        huge_obj = huge(sigmahat_k, method = "glasso",
                        lambda = lambda_mat[[l, n_index, i]], cov.output = TRUE, verbose = FALSE)
        
        toc()
        # --------------------- selecting opt lambda and opt prec ---------------------
        tic("sel_opts")
        sel_opts = get_optimal(huge_obj = huge_obj, nlambdas = length(lambda_mat[[l,n_index,i]]), n = n,
                               p = p, sigmahat = sigmahat_k, k = k, tol = tol)
        toc()
  
        # lambdas
        glasso_opt_lambda[j,i] = sel_opts[[2]]
  
        # evaluation measures
        evals[[i]][j,] = get_eval(est_prec = as.matrix(sel_opts[[1]]), prec = prec, RUN_ENT = TRUE)
  
        est_diags_temp = get_matrix_from_blocks(M = as.matrix(sel_opts[[1]]), p = p, nblocks = nblocks)
        diag_evals[[i]][j,] = get_eval(est_prec = est_diags_temp[[1]], prec = diags, RUN_ENT = FALSE)
        off_diag_evals[[i]][j,] = get_eval(est_prec = est_diags_temp[[2]], prec = off_diags, RUN_ENT = FALSE)
        
        FP_evals[[i]][j,] = prop_of_FP(prec = prec, est_prec = as.matrix(sel_opts[[1]]), diags = diags, est_diags = est_diags_temp[[1]],
                                     off_diags = off_diags, est_off_diags = est_diags_temp[[2]])
        
        # target sparsity evals
        sparsity_i = get_target_sparsity_index(huge_obj = huge_obj, k = k, tol = tol, nlambdas = length(huge_obj$lambda), target_sparsity = prec_sparsity, dev_tol = dev_tol)
        if(sparsity_i == 0){
          sparsity_evals[[i]][j,] = integer(nmeasures) 
          FP_spars[[i]][j,] = integer(2)
        } else{
          sparsity_prec = get_k_glasso(huge_obj = huge_obj, lambda_ind = sparsity_i, k = k, tol = tol)
          est_diags_temp = get_matrix_from_blocks(M = sparsity_prec, p = p, nblocks = nblocks)
          
          sparsity_evals[[i]][j,] = get_eval(est_prec = sparsity_prec,
                                           prec = prec, RUN_ENT = TRUE)
          FP_spars[[i]][j,] = prop_of_FP(prec = prec, est_prec = sparsity_prec, diags = diags, est_diags = est_diags_temp[[1]],
                                       off_diags = off_diags, est_off_diags = est_diags_temp[[2]])
          N_sparsity[i] = N_sparsity[i] + 1
        }
        
        # target fdr evals
        fdr_i = get_target_fdr_index(huge_obj = huge_obj, prec = prec, k = k, tol = tol, nlambdas = length(huge_obj$lambda), target_fdr = target_fdr, dev_tol = dev_tol)
        if(fdr_i == 0){
          fdr_evals[[i]][j,] = integer(nmeasures)
          FP_fdr[[i]][j,] = integer(2)
        } else{
          fdr_prec = get_k_glasso(huge_obj = huge_obj, lambda_ind = fdr_i, k = k, tol = tol)
          est_diags_temp = get_matrix_from_blocks(M = fdr_prec, p = p, nblocks = nblocks)
          fdr_evals[[i]][j,] = get_eval(est_prec = fdr_prec,
                                      prec = prec, RUN_ENT = TRUE)
          FP_fdr[[i]][j,] = prop_of_FP(prec = prec, est_prec = fdr_prec, diags = diags, est_diags = est_diags_temp[[1]],
                                       off_diags = off_diags, est_off_diags = est_diags_temp[[2]])
          N_FDR[i] = N_FDR[i] + 1
        }
        
        paths[[i]] = paths[[i]] + get_FS_paths(huge_obj, k, tol, opt_lambda = sel_opts[[2]], plot = FALSE)
        
        
        bics[[j,i]] = get_bic(huge_obj = huge_obj, nlambdas = length(huge_obj$lambda), n = n, p = p, sigmahat = sigmahat_k, k = k, tol = tol)

        # Here we can print, plot and save a sample from 1 random replicate 
        # (runs for all p, n and k for that specific random N)
        if(j == experiment_sample){
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
        }
        i = i + 1
      }
      print(j)
    }
    
    
    # --------------------- RESULTS ---------------------
    # Plots -------------
    
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
    
    aver_diag_eval = get_aver_eval(k_seq = k_seq, evals = diag_evals, N = N)
    aver_diag_eval = aver_diag_eval[,-c(4,5,9,10)] # without MCC and entropy loss
    sd_diag = get_sd_evals(k_seq = k_seq, evals = diag_evals, N = N)
    
    
    aver_off_diag_eval = get_aver_eval(k_seq = k_seq, evals = off_diag_evals, N = N)
    aver_off_diag_eval = aver_off_diag_eval[,-c(4,5,9,10)] # without MCC and entropy loss
    sd_off_diag = get_sd_evals(k_seq = k_seq, evals = off_diag_evals, N = N)
    
    sparsity_aver = get_aver_eval(k_seq = k_seq, evals = sparsity_evals, N = N, N_divide = N_sparsity, t = t)
    sd_sparsity = get_sd_evals(k_seq = k_seq, evals = sparsity_evals, N = N, N_divide = N_sparsity)
    
    fdr_aver = get_aver_eval(k_seq = k_seq, evals = fdr_evals, N = N, N_divide = N_FDR, t = t)
    sd_fdr = get_sd_evals(k_seq = k_seq, evals = fdr_evals, N = N, N_divide = N_FDR, t = t) # 
    
    aver_FP_eval = get_aver_eval(k_seq = k_seq, evals = FP_evals, N = N)
    sd_FP_eval = get_sd_evals(k_seq = k_seq, evals = FP_evals, N = N)
    
    aver_FP_spars= get_aver_eval(k_seq = k_seq, evals = FP_spars, N = N, N_divide = N_sparsity, t = t)
    sd_FP_spars = get_sd_evals(k_seq = k_seq, evals = FP_spars, N = N, N_divide = N_sparsity, t = t)
    
    aver_FP_fdr = get_aver_eval(k_seq = k_seq, evals = FP_fdr, N = N, N_divide = N_FDR, t = t)
    sd_FP_fdr = get_sd_evals(k_seq = k_seq, evals = FP_fdr, N = N, N_divide = N_FDR, t = t)
    
    
    txtBIC = xtable(aver_eval,digits = 3)
    txtBIC_sd = xtable(sd_eval,digits = 3)
    txtFDR = xtable(fdr_aver,digits = 3)
    txtspars = xtable(sparsity_aver,digits = 3)
    txtdiag = xtable(aver_diag_eval,digits = 3)
    txtoff_diag = xtable(aver_off_diag_eval,digits = 3)
    txtFDR_sd = xtable(sd_fdr,digits = 3)
    txtspars_sd = xtable(sd_sparsity,digits = 3)
    txtoff_diag_sd = xtable(sd_off_diag,digits = 3)
    txtdiag_sd = xtable(sd_diag,digits = 3)
    txtlambda_sd = xtable(sd_lambda, digits = 3)
    
    txtFP_eval = xtable(aver_FP_eval, digits = 3)
    txtFP_eval_sd = xtable(sd_FP_eval, digits = 3)
    txtFP_spars = xtable(aver_FP_spars, digits = 3)
    txtFP_spars_sd = xtable(sd_FP_spars, digits = 3)
    txtFP_fdr = xtable(aver_FP_fdr, digits = 3)
    txtFP_fdr_sd = xtable(sd_FP_fdr, digits = 3)
    
    print(txtBIC, file=paste(result_dir, "/", models[l],"_n", n, "-tableBIC.txt", sep = ""))
    print(txtBIC_sd, file=paste(result_dir, "/", models[l],"_n", n, "-tableBIC_sd.txt", sep = ""))
    print(txtFDR, file=paste(result_dir, "/", models[l],"_n", n, "-tableFDR.txt", sep = ""))
    print(txtspars, file=paste(result_dir, "/", models[l],"_n", n, "-tablespars.txt", sep = ""))
    print(txtdiag, file=paste(result_dir, "/", models[l],"_n", n, "-tablediag.txt", sep = ""))
    print(txtoff_diag, file=paste(result_dir, "/", models[l],"_n", n, "-tableoff_diag.txt", sep = ""))
    print(txtFDR_sd, file=paste(result_dir, "/", models[l],"_n", n, "-tableFDR_sd.txt", sep = ""))
    print(txtspars_sd, file=paste(result_dir, "/", models[l],"_n", n, "-tablespars_sd.txt", sep = ""))
    print(txtoff_diag_sd, file=paste(result_dir, "/", models[l],"_n", n, "-tableoff_diag_sd.txt", sep = ""))
    print(txtdiag_sd, file=paste(result_dir, "/", models[l],"_n", n, "-tablediag_sd.txt", sep = ""))
    print(txtlambda_sd, file=paste(result_dir, "/", models[l],"_n", n, "-tablelambda_sd.txt", sep = ""))
    
    print(txtFP_eval, file=paste(result_dir, "/", models[l],"_n", n, "-tableFP_eval.txt", sep = ""))
    print(txtFP_eval_sd, file=paste(result_dir, "/", models[l],"_n", n, "-tableFP_eval_sd.txt", sep = ""))
    print(txtFP_spars, file=paste(result_dir, "/", models[l],"_n", n, "-tableFP_spars.txt", sep = ""))
    print(txtFP_spars_sd, file=paste(result_dir, "/", models[l],"_n", n, "-tableFP_spars_sd.txt", sep = ""))
    print(txtFP_fdr, file=paste(result_dir, "/", models[l],"_n", n, "-tableFP_fdr.txt", sep = ""))
    print(txtFP_fdr_sd, file=paste(result_dir, "/", models[l],"_n", n, "-tableFP_fdr_sd.txt", sep = ""))
    

    model_bics[[n_index]] = get_aver_bics(k_seq = k_seq, bics = bics, N = N)
    n_index = n_index + 1
  }
  plot_bics_normalize(model_bics = model_bics, save = TRUE, show_title = TRUE,
            file = paste(file.path(main_dir, result_dir), "/", models[l], "_bics_normalized.png", sep = ""))
  
  toc()
}
toc()

