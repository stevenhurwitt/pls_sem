library(seminr)
library(lavaan)

#make dataframe
set.seed(123)
# Define sample size
N <- 100
# Define the Fixed Parameters
alpha <- c(1, 0.5)  # latent means
Phi <- matrix(c(1, 0.1, 
                0.1, 0.2), nrow = 2)  # latent variances/covariances
Lambda <- cbind(c(1, 1, 1, 1), 
                c(0, 1, 2, 3))  # factor loadings
Theta <- diag(0.5, nrow = 4)
# Generate latent factor scores
eta <- rmnorm(N, mean = alpha, varcov = Phi)
# Generate residuals:
e <- rmnorm(N, varcov = Theta)
# Compute outcome scores: y_i = t(Lambda %*% eta_i) + e
y <- tcrossprod(eta, Lambda) + e
# Make it a data frame
colnames(y) <- paste0("y", 1:4)
df <- as.data.frame(y)

#run lavaan model
m1 <- 'i =~ 1*y1 + 1*y2 + 1*y3 + 1*y4
       s =~ 0*y1 + 1*y2 + 2*y3 + 3*y4
       i ~~ s'

m1_fit <- growth(m1, data = df)

fitMeasures(m1_fit)
coef(m1_fit)
vcov(m1_fit)

bootstrapLavaan(m1_fit, R = 1000L, type = "bollen.stine",
                FUN = fitMeasures, fit.measures = c("cfi"))

#create latent slope and intercept
fl = as.matrix(c(0, 1, 2, 3))
df$s = as.matrix(df) %*% fl

df$i = rowSums(df[,1:4])

#create measurement model
df_mm = constructs(
  composite('y1', single_item('y1'), weights = mode_A),
  composite('y2', single_item('y2'), weights = mode_A),
  composite('y3', single_item('y3'), weights = mode_A),
  composite('y4', single_item('y4'), weights = mode_A),
  composite('i', single_item('i'), weights = mode_A),
  composite('s', single_item('s'), weights = mode_A)
)


#create structural model
df_sm = relationships(
  paths(from = 'i', to = c('y1', 'y2', 'y3', 'y4', 's')),
  paths(from = 's', to = c('y1', 'y2', 'y3', 'y4', 'i'))
)


#parameter estimates
df_pls = estimate_pls(data = df,
                      measurement_model = df_mm,
                      structural_model = df_sm,
                      inner_weights = path_weighting)

summary(df_pls)

#bootstrap
boot_df_pls <- bootstrap_model(seminr_model = df_pls,
                                 nboot = 1000,
                                 cores = 2)

summary(boot_df_pls)
