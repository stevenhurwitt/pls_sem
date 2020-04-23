library(lavaan)

getwd()
setwd("C:/Users/steve/Documents/gerson")
pol_dem = read.csv("Political_Democracy.csv", sep = ",", header = TRUE)

# measurement model
model = "ind60 =~ x1 + x2 + x3
dem60 =~ y1 + y2 + y3 + y4
dem65 =~ y5 + y6 + y7 + y8
# regressions
dem60 ~ ind60
dem65 ~ ind60 + dem60
# residual correlations
y1 ~~ y5
y2 ~~ y4 + y6
y3 ~~ y7
y4 ~~ y8
y6 ~~ y8"

fit = sem(model, data=PoliticalDemocracy)
summary(fit, standardized=TRUE)

library(seminr)
pol_dem$ind60 = pol_dem$x1 + pol_dem$x2 + pol_dem$x3
pol_dem$dem60 = pol_dem$y1 + pol_dem$y2 + pol_dem$y3 + pol_dem$y4
pol_dem$dem65 = pol_dem$y5 + pol_dem$y6 + pol_dem$y7 + pol_dem$y8

model.mm = constructs(
  composite('y1'. single_item('y1'), weights = mode_A),
  composite('y1'. single_item('y1'), weights = mode_A),
  composite('y1'. single_item('y1'), weights = mode_A),
  composite('y1'. single_item('y1'), weights = mode_A),
  composite('y1'. single_item('y1'), weights = mode_A),
  composite('y1'. single_item('y1'), weights = mode_A),
  composite('y1'. single_item('y1'), weights = mode_A),
  
)