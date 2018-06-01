#models = c("model_p200", "model_p500", "Model-4-200", "model_4-200-offlinks")
models = "gene"
# ns = c(0.4, 5)
ns = 500
k_seq = c(1,1.5,2)

lambda_mat = array(list(), c(length(models), length(ns), length(k_seq)))

# genedata
ls = seq(0.7,0.1, length.out = 48)
lambda_mat[1,1,1] = list(ls^(1.2)*1.3) 
lambda_mat[1,1,2] = list(ls^(1/1.7)/4) 
lambda_mat[1,1,3] = list((ls)^(1/6.8)/8) 



ls = seq(0.187, 0.045, length.out = 48)
lambda_mat[1,1,1] = list(ls) # p200_n80_k1
lambda_mat[1,1,2] = list(ls) # p200_n80_k15
lambda_mat[1,1,3] = list((ls)^(1/1.4)/1.5) #p200_n80_k2

ls = seq(0.105, 0.014, length.out = 48)
lambda_mat[1,2,1] = list(ls) # p200_n1000_k1
lambda_mat[1,2,2] = list(ls) # p200_n1000_k15
lambda_mat[1,2,3] = list((ls)^(1/1.2)*1.5) # p200_n1000_k2

# ls = seq(0.02, 0.004, length.out = 48)
# lambda_mat[2,1,1] = list(ls) # p1000_n400_k1
# lambda_mat[2,1,2] = list(ls^(1/1.1)*1.1) # p1000_n400_k15
# lambda_mat[2,1,3] = list((ls)*4.5) #p1000_n400_k2

# ls = seq(0.16, 0.05, length.out = 48)
# lambda_mat[2,1,1] = list(ls) # p200old_n80_k1
# lambda_mat[2,1,2] = list(ls/1.1) # p200old_n80_k15
# lambda_mat[2,1,3] = list((ls)^(1/2.3)/3.9) #p200old_n80_k2
# 
# ls = seq(0.05, 0.018, length.out = 48)
# lambda_mat[2,2,1] = list(ls) # p200old_n1000_k1
# lambda_mat[2,2,2] = list(ls/1.05) # p200old_n1000_k15
# lambda_mat[2,2,3] = list(ls^(1/1.5)/3.1) # p200old_n1000_k2

ls = seq(0.013, 0.004, length.out = 48)
lambda_mat[2,1,1] = list((ls*3)) # p500_n200_k1
lambda_mat[2,1,2] = list((ls)^(1/1.1)*2.5) # p500o_n200_k15
lambda_mat[2,1,3] = list((ls)^(1/1.17)*5) #p500_n200_k2


ls = seq(0.01, 0.0009, length.out = 48)
lambda_mat[2,2,1] = list((ls^(1/1.2))*1.1) # p500_n2500_k1
lambda_mat[2,2,2] = list((ls^(1/1.25))*1.3) # p500_n2500_k15
lambda_mat[2,2,3] = list((ls^(1/1.43))*3.8) # p500_n2500_k2




ls = seq(0.635, 0.40, length.out = 48)
lambda_mat[3,1,1] = list(ls*3.7) # p200ava_n80_k1
lambda_mat[3,1,2] = list(ls) # p200ava_n80_k15
lambda_mat[3,1,3] = list((ls)^(1/2.4)/1.5) # p200ava_n80_k2

ls = seq(0.2, 0.13, length.out = 48)
lambda_mat[3,2,1] = list((ls)*3.7) # p200ava_n1000_k1
lambda_mat[3,2,2] = list(ls) # p200ava_n1000_k15
lambda_mat[3,2,3] = list(((ls)^(1/1.5))/1.8) # p200ava_n1000_k2

ls = seq(0.117, 0.046, length.out = 48)
lambda_mat[4,1,1] = list(ls) # p200avaoff_n80_k1
lambda_mat[4,1,2] = list(ls/1.1) # p200avaoff_n80_k15
lambda_mat[4,1,3] = list((ls)^(1/3.1)/5) # p200avaoff_n80_k2

ls = seq(0.023, 0.016, length.out = 48)
lambda_mat[4,2,1] = list(ls) # p200avaoff_n1000_k1
lambda_mat[4,2,2] = list(ls/1.1) # p200avaoff_n1000_k15
lambda_mat[4,2,3] = list(((ls)^(1/2.3))/7.52) # p200avaoff_n1000_k2

save(lambda_mat,
     file = "lambda_mat_gene.R")

# save(lambda_mat,
     file = "lambda_mat.R")
#load("lambda_mat.R")

