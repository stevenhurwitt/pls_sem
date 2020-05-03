library(seminr)
library(mnormt) 


getwd()
setwd("/media/steven/samsung_t5/gerson")
set.seed(421)
N = 100

#generate data
opt_lambda = c(.86, .88, .8)
opt_alpha = c(2, 4, 6)
opt_phi = matrix(c(1, .2, .3, .2, .3, 1, .3, .2, 1), nrow = 3)
opt_eta = rmnorm(N, mean = opt_alpha, varcov = opt_phi)
#opt_eta = ceiling(opt_eta)
colnames(opt_eta)= c('opt1', 'opt2', 'opt3')


pt_lambda = c(.78, .8, .64, .8, .75)
pt_alpha = c(4, 3, 5, 3, 4)
pt_phi = matrix(c(1, .2, .3, .4, .5, 
                  .2, 1, .4, .5, .3,
                  .3, .4, .5, 1, .2,
                  .4, .5, .3, 1, .2,
                  .5, .4, .3, .2, 1), nrow = 5)
pt_eta = rmnorm(N, mean = pt_alpha, varcov = pt_phi)
#pt_eta = ceiling(pt_eta)
colnames(pt_eta)= c('pt1', 'pt2', 'pt3', 'pt4', 'pt5')


fse_lambda = c(.85, .89, .7)
fse_alpha = c(3, 5, 4)
fse_phi = matrix(c(1, .2, .3, .2, .3, 1, .3, .2, 1), nrow = 3)
fse_eta = rmnorm(N, mean = fse_alpha, varcov = fse_phi)
#fse_eta = ceiling(fse_eta)
colnames(fse_eta)= c('fse1', 'fse2', 'fse3')


teamperf_lambda = c(.85, .82, .75)
teamperf_alpha = c(3, 5, 4)
teamperf_phi = matrix(c(1, .2, .3, .2, .3, 1, .3, .2, 1), nrow = 3)
teamperf_eta = rmnorm(N, mean = teamperf_alpha, varcov = teamperf_phi)
#teamperf_eta = ceiling(teamperf_eta)
colnames(teamperf_eta)= c('teamperf1', 'teamperf2', 'teamperf3')

#calculate variance based on CR
lambda = c(opt_lambda, pt_lambda, fse_lambda, teamperf_lambda)
CR = c(.4, .6, .8)

cr = sum(lambda)^2/(sum(lambda)^2 + sum_err_var)

sum_err_var = sum(lambda)^2/CR - sum(lambda)^2
err_var = sum_err_var / 14

#generate errors
theta1 = diag(err_var[1], nrow = 3)
theta1_pt = diag(err_var[1], nrow = 5)
e1_opt = rmnorm(N, varcov = theta1)
e1_pt = rmnorm(N, varcov = theta1_pt)
e1_fse = rmnorm(N, varcov = theta1)
e1_teamperf = rmnorm(N, varcov = theta1)

theta2 = diag(err_var[2], nrow = 3)
theta2_pt = diag(err_var[2], nrow = 5)
e2_opt = rmnorm(N, varcov = theta1)
e2_pt = rmnorm(N, varcov = theta1_pt)
e2_fse = rmnorm(N, varcov = theta1)
e2_teamperf = rmnorm(N, varcov = theta1)

theta3 = diag(err_var[3], nrow = 3)
theta3_pt = diag(err_var[3], nrow = 5)
e3_opt = rmnorm(N, varcov = theta1)
e3_pt = rmnorm(N, varcov = theta1_pt)
e3_fse = rmnorm(N, varcov = theta1)
e3_teamperf = rmnorm(N, varcov = theta1)

#create final variables factor loadings * normal data + error
opt_Lambda = cbind(opt_lambda, opt_lambda, opt_lambda)
pt_Lambda = cbind(pt_lambda, pt_lambda, pt_lambda, pt_lambda, pt_lambda)
fse_Lambda = cbind(fse_lambda, fse_lambda, fse_lambda)
teamperf_Lambda = cbind(teamperf_lambda, teamperf_lambda, teamperf_lambda)

final_opt = tcrossprod(opt_eta, opt_Lambda) + e1_opt
final_opt2 = tcrossprod(opt_eta, opt_Lambda) + e2_opt
final_opt3 = tcrossprod(opt_eta, opt_Lambda) + e3_opt
final_opt = ceiling(final_opt)
final_opt2 = ceiling(final_opt2)
final_opt3 = ceiling(final_opt3)

final_pt = tcrossprod(pt_eta, pt_Lambda) + e1_pt
final_pt2 = tcrossprod(pt_eta, pt_Lambda) + e2_pt
final_pt3 = tcrossprod(pt_eta, pt_Lambda) + e3_pt
final_pt = ceiling(final_pt)
final_pt2 = ceiling(final_pt2)
final_pt3 = ceiling(final_pt3)

final_fse = tcrossprod(fse_eta, fse_Lambda) + e1_fse
final_fse2 = tcrossprod(fse_eta, fse_Lambda) + e2_fse
final_fse3 = tcrossprod(fse_eta, fse_Lambda) + e3_fse
final_fse = ceiling(final_fse)
final_fse2 = ceiling(final_fse2)
final_fse3 = ceiling(final_fse3)

final_teamperf = tcrossprod(teamperf_eta, teamperf_Lambda) + e1_teamperf
final_teamperf2 = tcrossprod(teamperf_eta, teamperf_Lambda) + e2_teamperf
final_teamperf3 = tcrossprod(teamperf_eta, teamperf_Lambda) + e3_teamperf
final_teamperf = ceiling(final_teamperf)
final_teamperf2 = ceiling(final_teamperf2)
final_teamperf3 = ceiling(final_teamperf3)

#final dataset
df1 = cbind(final_opt, final_pt, final_fse, final_teamperf)
df2 = cbind(final_opt2, final_pt2, final_fse2, final_teamperf2)
df3 = cbind(final_opt3, final_pt3, final_fse3, final_teamperf3)
colnames(df1) = c('opt1', 'opt2', 'opt3', 'pt1', 'pt2', 'pt3', 'pt4', 'pt5', 'fse1', 'fse2', 'fse3', 'teamperf1', 'teamperf2', 'teamperf3')
colnames(df2) = c('opt1', 'opt2', 'opt3', 'pt1', 'pt2', 'pt3', 'pt4', 'pt5', 'fse1', 'fse2', 'fse3', 'teamperf1', 'teamperf2', 'teamperf3')
colnames(df3) = c('opt1', 'opt2', 'opt3', 'pt1', 'pt2', 'pt3', 'pt4', 'pt5', 'fse1', 'fse2', 'fse3', 'teamperf1', 'teamperf2', 'teamperf3')

#write data
write.csv(df1, 'df1.csv', row.names = F)
write.csv(df2, 'df2.csv', row.names = F)
write.csv(df3, 'df3.csv', row.names = F)

### PLS SEM ###
#create measurement model matrix
df1_mm <- constructs(
  composite("COMP",         multi_items("opt", 1:3), weights = mode_A),
  composite("LIKE",   multi_items("pt", 1:5), weights = mode_A),
  composite("CUSA",       multi_items("fse", 1:3), weights = mode_A),
  composite("CUSL",         multi_items("teamperf", 1:3), weights = mode_A)
)

#create structural model
df1_sm <- relationships(
  paths(from = "COMP",        to = c("CUSA", "CUSL")),
  paths(from = "LIKE",  to = c("CUSA", "CUSL")),
  paths(from = "CUSA",      to = "CUSL")
)

#model estimation
df1_pls <- estimate_pls(data = df1,
                         measurement_model = df1_mm,
                         structural_model = df1_sm,
                         inner_weights = path_weighting)

summary(df1_pls)
df1_pls$path_coef
df1_pls$outer_loadings

#bootstrap
boot_df1_pls <- bootstrap_model(seminr_model = df1_pls,
                                 nboot = 5000,
                                 cores = 2)

summary(boot_df1_pls)

#create measurement model matrix
df2_mm <- constructs(
  composite("COMP",         multi_items("opt", 1:3), weights = mode_A),
  composite("LIKE",   multi_items("pt", 1:5), weights = mode_A),
  composite("CUSA",       multi_items("fse", 1:3), weights = mode_A),
  composite("CUSL",         multi_items("teamperf", 1:3), weights = mode_A)
)

#create structural model
df2_sm <- relationships(
  paths(from = "COMP",        to = c("CUSA", "CUSL")),
  paths(from = "LIKE",  to = c("CUSA", "CUSL")),
  paths(from = "CUSA",      to = "CUSL")
)

#model estimation
df2_pls <- estimate_pls(data = df2,
                        measurement_model = df2_mm,
                        structural_model = df2_sm,
                        inner_weights = path_weighting)

summary(df2_pls)

#bootstrap
boot_df2_pls <- bootstrap_model(seminr_model = df2_pls,
                                nboot = 5000,
                                cores = 2)

summary(boot_df2_pls)

#create measurement model matrix
df3_mm <- constructs(
  composite("COMP",         multi_items("opt", 1:3), weights = mode_A),
  composite("LIKE",   multi_items("pt", 1:5), weights = mode_A),
  composite("CUSA",       multi_items("fse", 1:3), weights = mode_A),
  composite("CUSL",         multi_items("teamperf", 1:3), weights = mode_A)
)

#create structural model
df3_sm <- relationships(
  paths(from = "COMP",        to = c("CUSA", "CUSL")),
  paths(from = "LIKE",  to = c("CUSA", "CUSL")),
  paths(from = "CUSA",      to = "CUSL")
)

#model estimation
df3_pls <- estimate_pls(data = df3,
                        measurement_model = df3_mm,
                        structural_model = df3_sm,
                        inner_weights = path_weighting)

summary(df3_pls)

#bootstrap
boot_df3_pls <- bootstrap_model(seminr_model = df3_pls,
                                nboot = 5000,
                                cores = 2)


summary(boot_df3_pls)
