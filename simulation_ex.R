library(seminr)
library(mnormt)

getwd()
setwd("D:/gerson")
set.seed(421)

opt_lambda = c(.86, .88, .8)
opt_alpha = c(2, 4, 6)
opt1
opt2
opt3

pt_lambda = c(.78, .8, .64, .8, .75)
pt_alpha = c(4, 3, 5, 3, 4)
pt1
pt2
pt3
pt4
pt5

fse_lambda = c(.85, .89, .7)
fse_alpha = c(3, 5, 4)
fse1
fse2
fse3

teamperf_lambda = c(.85, .82, .75)
teamperf_alpha = c(3, 5, 4)
teamperf1
teamperf2
teamperf3

lambda = c(opt_lambda, pt_lambda, fse_lambda, teamperf_lambda)
CR = c(.4, .6, .8)

cr = sum(lambda)^2/(sum(lambda)^2 + sum_err_var)

sum_err_var = sum(lambda)^2/CR - sum(lambda)^2
err_var = sum_err_var / 14

theta1 = diag(err_var[1], nrow = 3)
theta1_pt = diag(err_var[1], nrow = 5)

theta2 = diag(err_var[2], nrow = 3)
theta2_pt = diag(err_var[2], nrow = 5)

theta3 = diag(err_var[3], nrow = 3)
theta3_pt = diag(err_var[3], nrow = 5)
