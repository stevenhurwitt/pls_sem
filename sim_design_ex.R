library(seminr)
library(mnormt)
library(rlist)

getwd()
setwd("D:\\gerson")
load("D:/gerson/sim_design_ex.RData")
set.seed(421)

n_input = list(50, 100, 500)

var_names = c('opt1','opt2', 'opt3', 'pt1', 'pt2', 'pt3', 'pt4', 'pt5', 'fse1','fse2','fse3', 'teamperf1', 'teamperf2', 'teamperf3', 'ict1','ict2','ict3', 'ITexp1', 'ITexp2', 'ITexp3', 'infint1', 'infint2', 'infint3', 'lmx1', 'lmx2', 'lmx3', 'lmx4', 'lmx5', 'lmx6', 'lmx7')

## param specification
opt_lambda = c(.86, .88, .8)
pt_lambda = c(.78, .8, .64, .8, .75)
fse_lambda = c(.85, .89, .7)
teamperf_lambda = c(.85, .82, .75)
ict_lambda = c(.83, .68, .67)
ITexp_lambda = c(.84, .74, .75)
infint_lambda = c(.73, .82, .78)
lmx_lambda = c(.74, .68, .75, .78, .7, .62, .79)

lambda_input = list(opt_lambda, pt_lambda, fse_lambda, teamperf_lambda, ict_lambda, ITexp_lambda, infint_lambda, lmx_lambda)

opt_alpha = c(2, 4, 6)
pt_alpha = c(4, 3, 5, 3, 4)
fse_alpha = c(3, 5, 4)
teamperf_alpha = c(3, 4, 5)
ict_alpha =  c(4, 3, 4)
ITexp_alpha = c(2, 5, 4)
infint_alpha = c(5, 4, 3)
lmx_alpha = c(3, 4, 3, 5, 2, 6, 3)

alpha_input = list(opt_alpha, pt_alpha, fse_alpha, teamperf_alpha, ict_alpha, ITexp_alpha, infint_alpha, lmx_alpha)

opt_phi = matrix(c(1, .2, .3, 
                   .2, 1, .3, 
                   .3, .3, 1), nrow = 3)

pt_phi = matrix(c(1, .2, .3, .4, .5, 
                  .2, 1, .4, .5, .4,
                  .3, .4, 1, .5, .2,
                  .4, .5, .5, 1, .2,
                  .5, .4, .2, .2, 1), nrow = 5)

fse_phi = matrix(c(1, .2, .3, 
                   .2, 1, .2, 
                   .3, .2, 1), nrow = 3)

teamperf_phi = matrix(c(1, .3, .2, 
                        .3, 1, .2,
                        .2, .2, 1), nrow = 3)

ict_phi =matrix(c(1, .3, .5, 
                  .3, 1, .4,
                  .5, .4, 1), nrow = 3)

ITexp_phi = matrix(c(1, .5, .3, 
                     .5, 1, .3,
                     .3, .3, 1), nrow = 3)

infint_phi = matrix(c(1, .2, .3, 
                      .2, 1, .2,
                      .3, .2, 1), nrow = 3)

lmx_phi = matrix(c(1, .3, .4, .2, .5, .4, .3,
                   .3, 1, .2, .4, .5, .3, .4,
                   .4, .2, 1, .4, .3, .2, .4,
                   .2, .4, .4, 1, .3, .5, .2,
                   .5, .5, .3, .3, 1, .4, .4,
                   .4, .3, .2, .5, .4, 1, .3,
                   .3, .4, .4, .2, .4, .3, 1), nrow = 7)

phi_input = list(opt_phi, pt_phi, fse_phi, teamperf_phi, ict_phi, ITexp_phi, infint_phi, lmx_phi)

lambda = c(opt_lambda, pt_lambda, fse_lambda, teamperf_lambda, ict_lambda, ITexp_lambda, infint_lambda, lmx_lambda)

CR = c(.5, .7, .9)

sum_err_var = sum(lambda)^2/CR - sum(lambda)^2
err_var = sum_err_var / 30
err_var = as.list(err_var)



## data gen
master_data_gen = function(N, err, alpha, phi, lambda, lambda_input){

eta_gen = function(N, alpha, phi){
  eta = rmnorm(N, mean = alpha, varcov = phi)
  return(eta)
}


raw_data_list = mapply(function(alpha, phi){eta_gen(N, alpha, phi)},
                       alpha = alpha, phi = phi)


# error var
error_gen = function(error_var, alpha, N){
  theta = diag(error_var, nrow = length(alpha))
  error = rmnorm(N, varcov = theta)
  return(error)
}


errors_input = lapply(alpha, function(alpha){error_gen(err, alpha, N)})

# gen lambda * alpha + error
lambda_gen = function(lambda_vec){
  Lambda = matrix(lambda_vec, length(lambda_vec), length(lambda_vec))
  return(t(Lambda))
}

Lambda_input = lapply(lambda_input, lambda_gen)


data = function(eta, lambda, error){
  return(tcrossprod(eta, lambda) + error)}

final_data_list = mapply(function(eta, lambda, error){data(eta, lambda, error)}, eta = raw_data_list, lambda = Lambda_input, error = errors_input)

final_data = data.frame(matrix(unlist(final_data_list), ncol = 30), stringsAsFactors = F)

colnames(final_data) = var_names
final_data = ceiling(final_data)
return(final_data)
}


err_var
n_input
df_list = list()

for (i in 1:length(n_input)){
  for (j in 1:length(err_var)){
    df = master_data_gen(n_input[[i]], err_var[[j]], alpha_input, phi_input, lambda, lambda_input)
    df_list = list.append(df_list, df)
  }
}

PLS_SEM = function(data){
## pls sem model
df_mm <- constructs(
  composite("Optimism",         multi_items("opt", 1:3), weights = mode_A),
  composite("Perspective",   multi_items("pt", 1:5), weights = mode_A),
  composite("Self",       multi_items("fse", 1:3), weights = mode_A),
  composite("Team Performance",         multi_items("teamperf", 1:3), weights = mode_A),
  composite("ICT",       multi_items("ict", 1:3), weights = mode_A),
  composite("IT",       multi_items("ITexp", 1:3), weights = mode_A),
  composite("Info",       multi_items("infint", 1:3), weights = mode_A),
  composite("Leader",       multi_items("lmx", 1:7), weights = mode_A)
)

#create structural model
df_sm <- relationships(
  paths(from = "Optimism",        to = c("Leader", "Self")),
  paths(from = "Perspective",  to = c("Leader", "Self", "ICT", "IT")),
  paths(from = "Self",      to = "Team Performance"),
  paths(from = "ICT",      to = "IT"),
  paths(from = "IT",      to = "Team Performance"),
  paths(from = "Info",      to = c("Leader", "Self", "ICT", "IT")),
  paths(from = "Leader",      to = "Team Performance")
)

#model estimation
df_pls <- estimate_pls(data = data,
                        measurement_model = df_mm,
                        structural_model = df_sm,
                        inner_weights = path_weighting)

return(df_pls)}

pls_results = lapply(df_list, PLS_SEM)

#bootstrap
bootstrap_pls = function(pls){
boot_df_pls <- bootstrap_model(seminr_model = pls,
                                nboot = 5000,
                                cores = 6)

return(boot_df_pls)}

bootstrap_results = lapply(pls_results, bootstrap_pls)
