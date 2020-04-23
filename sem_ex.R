library(tidyverse)
library(mnormt)
library(lavaan)
library(SimDesign)

theme_set(theme_classic() + theme(panel.grid.major.y = element_line(color = "grey92")))

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

#generate test data
N_test <- 1e5
eta_test <- rmnorm(N_test, mean = alpha, varcov = Phi)
colMeans(eta_test)

cov(eta_test)

e_test <- rmnorm(N_test, varcov = Theta)
colMeans(e_test)

cov(e_test)

## function to gen data
gen_lgm_data <- function(N, alpha, Phi, Lambda, Theta) {
  # Generate latent factor scores
  eta <- rmnorm(N, mean = alpha, varcov = Phi)
  # Generate residuals:
  e <- rmnorm(N, varcov = Theta)
  # Compute outcome scores
  y <- tcrossprod(eta, Lambda) + e
  colnames(y) <- paste0("y", 1:4)
  # Make it a data frame
  as.data.frame(y)
}
# Test it:
set.seed(123)
gen_lgm_data(100, 
             alpha = c(1, 0.5), 
             Phi = matrix(c(1, 0.1, 
                            0.1, 0.2), nrow = 2), 
             Lambda = cbind(c(1, 1, 1, 1), 
                            c(0, 1, 2, 3)), 
             Theta = diag(0.5, nrow = 4)) %>% 
head() # shows only the first six cases

##test w/ lavaan
# True model
m1 <- 'i =~ 1*y1 + 1*y2 + 1*y3 + 1*y4
       s =~ 0*y1 + 1*y2 + 2*y3 + 3*y4
       i ~~ s'
# Model without random slopes
m2 <- 'i =~ 1*y1 + 1*y2 + 1*y3 + 1*y4
       s =~ 0*y1 + 1*y2 + 2*y3 + 3*y4
       s ~~ 0*i + 0*s'
m1_fit <- growth(m1, data = df)
m2_fit <- growth(m2, data = df)

fitMeasures(m1_fit)  # return more than 30 fit indices
# Specific fit indices
fitMeasures(m1_fit, c("chisq", "df", "pvalue", "cfi", "rmsea"))
coef(m1_fit)  # parameter estimates
vcov(m1_fit)  # asymptotic covariance matrix of parameters
# Asymptotic standard errors of latent mean of s
sqrt(diag(vcov(m1_fit))["s~1"])
parameterEstimates(m1_fit, standardized = FALSE)  # coefficients, SE, CI
modificationIndices(m1_fit, minimum.value = 3.84)  # modification indices
# Bootstrapped fit index
bootstrapLavaan(m1_fit, R = 1000L, type = "bollen.stine",
                FUN = fitMeasures, fit.measures = c("cfi"))

##small simulation
set.seed(515)  # set the seed for reproducibility
NREP <- 500  # number of replications
# Fixed parameters (all caps)
ALPHA1 <- 1  # latent mean of intercepts
PHI11 <- 1  # intercept variance
LAMBDA <- cbind(c(1, 1, 1, 1), 
                c(0, 1, 2, 3))  # factor loadings
THETA <- diag(0.5, nrow = 4)  # residual variances

# Function for generating data:
gen_lgm_data <- function(N, alpha, Phi, Lambda, Theta) {
  # Generate latent factor scores
  eta <- rmnorm(N, mean = alpha, varcov = Phi)
  # Generate residuals:
  e <- rmnorm(N, varcov = Theta)
  # Compute outcome scores
  y <- tcrossprod(eta, Lambda) + e
  colnames(y) <- paste0("y", 1:4)
  # Make it a data frame
  as.data.frame(y)
}
# lavaan syntax
# True model
m1 <- 'i =~ 1*y1 + 1*y2 + 1*y3 + 1*y4
       s =~ 0*y1 + 1*y2 + 2*y3 + 3*y4
       i ~~ s'
# Model without random slopes
m2 <- 'i =~ 1*y1 + 1*y2 + 1*y3 + 1*y4
       s =~ 0*y1 + 1*y2 + 2*y3 + 3*y4
       s ~~ 0*i + 0*s'

# Design factors:
DESIGNFACTOR <- expand.grid(
  N = c(50, 100, 200), 
  phi22 = c(0.1, 0.5), 
  alpha2 = c(1, 0.5)
)
# Add condition number:
DESIGNFACTOR <- rowid_to_column(DESIGNFACTOR, "cond")
DESIGNFACTOR

##function
runsim <- function(to_run,  # conditions to run
                   nrep,  # number of replications
                   alpha1 = ALPHA1,  # latent mean of intercepts
                   phi11 = PHI11,  # intercept variance
                   Lambda = LAMBDA,  # factor loadings
                   Theta = THETA,  # residual variances
                   designfactors = DESIGNFACTOR) {
  # Extract design parameters for the given condition
  N <- designfactors[to_run, "N"]
  phi22 <- designfactors[to_run, "phi22"]
  alpha2 <- designfactors[to_run, "alpha2"]
  # Put the values back to the matrix
  alpha <- c(alpha1, alpha2)
  Phi <- matrix(c(phi11, phi22 / 2,
                  phi22 / 2, phi22), nrow = 2)
  # Initialize place holders for the results
  rep_result <- vector("list", nrep)
  for (i in seq_len(nrep)) {
    # Generate data
    df <- gen_lgm_data(N, alpha, Phi, Lambda, Theta)
    # Run model 1
    m1_fit <- growth(m1, data = df)
    # Run model 2
    m2_fit <- growth(m2, data = df)
    # Save results 
    rep_result[[i]] <- list(m1_fit = m1_fit, 
                            m2_fit = m2_fit)
  }
  # Return results
  return(rep_result)
}

sim_results <- runsim(1, 2)

sim_results[[1]]

sim_results <- map(seq_len(nrow(DESIGNFACTOR)),  ~ runsim(.x, nrep = NREP))

with(sim_results[[1]][[1]], 
     tibble(coef(m1_fit)["s~1"], 
            sqrt(vcov(m1_fit)["s~1", "s~1"]), 
            coef(m2_fit)["s~1"], 
            sqrt(vcov(m2_fit)["s~1", "s~1"])))
# Now, wrap it as a function:
extract_coef <- function(res) {
  out <- with(res, 
              tibble(coef(m1_fit)["s~1"], 
                     sqrt(vcov(m1_fit)["s~1", "s~1"]), 
                     coef(m2_fit)["s~1"], 
                     sqrt(vcov(m2_fit)["s~1", "s~1"])))
  # Add names:
  names(out) <- c("m1_est", "m1_se", "m2_est", "m2_se")
  out
}
extract_coef(sim_results[[1]][[1]])

sim_outdata <- sim_results %>% 
  # Iterate across conditions
  map_dfr(
    # Iterate across replications
    ~ map_dfr(.x, extract_coef, .id = "rep"), 
    .id = "cond") %>% 
  # Make `cond` a integer variable
  mutate(rep = as.integer(rep), 
         cond = as.integer(cond)) %>% 
  # Add design factors
  left_join(DESIGNFACTOR)
# Save the results

getwd()
setwd("C:/Users/steve/Documents/gerson")
write_csv(sim_outdata, "example_sem_results.csv")

#long version
sim_outdata <- sim_outdata %>% 
  # Make it long format
  gather("var", "val", m1_est:m2_se) %>% 
  separate(var, c("model", "var")) %>% 
  spread(var, val) %>% 
  # (Optional) Rename the levels of design factors to make it better looking in 
  # the plots
  mutate(model = factor(model, labels = c("Correctly specified", 
                                          "Misspecified")), 
         alpha2_lab = as_factor(paste0("alpha[2] == ", alpha2)), 
         phi22_lab = as_factor(paste0("phi[22] == ", phi22)), 
         N_lab = as_factor(paste0("italic(N) == ", N)))
sim_sum <- sim_outdata %>% 
  # Summarize results by conditions
  group_by(model, alpha2, phi22, N, alpha2_lab, phi22_lab, N_lab) %>% 
  summarise(ave_est = mean(est), 
            emp_sd = sd(est), 
            ave_se = mean(se)) %>% 
  # Compute standardized bias and relative SE bias
  mutate(bias = ave_est - alpha2, 
         bias_mcse = emp_sd / sqrt(NREP), 
         std_bias = bias / emp_sd, 
         rse_bias = (ave_se - emp_sd) / emp_sd) %>% 
  ungroup()

#table
sim_sum %>% 
  select(model:N, 
         Bias = bias, `Monte Carlo Error` = bias_mcse, 
         `Standardized Bias` = std_bias, 
         `Relative SE Error` = rse_bias) %>% 
  knitr::kable(digits = 3L)

#graph
sim_outdata %>% 
  ggplot(aes(x = factor(N), y = est, color = model)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = alpha2)) + 
  facet_grid(alpha2_lab ~ phi22_lab, labeller = label_parsed) + 
  labs(x = "Sample Size (N)", y = "Estimates of Mean Slope")

#bias
sim_sum %>% 
  ggplot(aes(x = factor(N), y = std_bias, color = model)) + 
  geom_line(aes(group = model)) + 
  geom_hline(aes(yintercept = 0)) + 
  facet_grid(alpha2_lab ~ phi22_lab, labeller = label_parsed) + 
  labs(x = "Sample Size (N)", y = "Standardized Bias")

sim_sum %>% 
  ggplot(aes(x = factor(N), y = rse_bias, color = model)) + 
  geom_line(aes(group = model)) + 
  geom_hline(aes(yintercept = 0)) + 
  facet_grid(alpha2_lab ~ phi22_lab, labeller = label_parsed) + 
  labs(x = "Sample Size (N)", y = "Relative Standard Error Bias")

#sim design package
# Design factors:
DESIGNFACTOR <- expand.grid(
  N = c(50, 100, 200), 
  phi22 = c(0.1, 0.5), 
  alpha2 = c(0, 0.5)
)
# Add condition number:
DESIGNFACTOR <- rowid_to_column(DESIGNFACTOR, "cond")
DESIGNFACTOR

# Function for generating data:
gen_lgm_data <- function(condition, fixed_objects = NULL) {
  N <- condition$N
  phi22 <- condition$phi22
  alpha2 <- condition$alpha2
  alpha <- c(1, alpha2)
  Phi <- matrix(c(1, phi22 / 2,
                  phi22 / 2, phi22), nrow = 2)
  Lambda <- cbind(c(1, 1, 1, 1), 
                  c(0, 1, 2, 3))
  Theta <- diag(0.5, nrow = 4)
  # Generate latent factor scores
  eta <- rmnorm(N, mean = alpha, varcov = Phi)
  # Generate residuals:
  e <- rmnorm(N, varcov = Theta)
  # Compute outcome scores
  y <- tcrossprod(eta, Lambda) + e
  colnames(y) <- paste0("y", 1:4)
  # Make it a data frame
  as.data.frame(y)
}
# Test: generate data
# test_df <- gen_lgm_data(DESIGNFACTOR[1, ])

# lavaan syntax (to be passed through `fixed_objects`)
# True model
M1 <- 'i =~ 1*y1 + 1*y2 + 1*y3 + 1*y4
       s =~ 0*y1 + 1*y2 + 2*y3 + 3*y4
       i ~~ s'
# Model without random slopes
M2 <- 'i =~ 1*y1 + 1*y2 + 1*y3 + 1*y4
       s =~ 0*y1 + 1*y2 + 2*y3 + 3*y4
       s ~~ 0*i + 0*s'

# Analysis function
run_lgm <- function(condition, dat, fixed_objects = NULL) {
  # Run model 1
  m1_fit <- growth(fixed_objects$M1, data = dat)
  # Run model 2
  m2_fit <- growth(fixed_objects$M2, data = dat)
  # Extract parameter estimates and standard errors
  ret <- c(coef(m1_fit)["s~1"], 
           sqrt(vcov(m1_fit)["s~1", "s~1"]),
           coef(m2_fit)["s~1"], 
           sqrt(vcov(m2_fit)["s~1", "s~1"]))
  names(ret) <- c("m1_est", "m1_se", "m2_est", "m2_se")
  ret
}
# Test: run analyses for simulated data
# run_lgm(dat = test_df)

# Helper function for computing relative SE bias
rse_bias <- function(est_se, est) {
  est_se <- as.matrix(est_se)
  est <- as.matrix(est)
  est_se <- colMeans(est_se)
  emp_sd <- apply(est, 2L, sd)
  est_se / emp_sd - 1
}

# Evaluative function
evaluate_lgm <- function(condition, results, fixed_objects = NULL) {
  alpha2 <- condition$alpha2
  c(bias = bias(results[ , c(1, 3)], parameter = alpha2), 
    std_bias = bias(results[ , c(1, 3)], parameter = alpha2, 
                    type = "standardized"), 
    rmse = RMSE(results[ , c(1, 3)], parameter = alpha2), 
    rse_bias = rse_bias(results[ , c(2, 4)], results[ , c(1, 3)])
  )
}

# Put all together
sim_results <- runSimulation(DESIGNFACTOR, 
                             replications = 500, 
                             generate = gen_lgm_data, 
                             analyse = run_lgm, 
                             summarise = evaluate_lgm, 
                             fixed_objects = 
                               list(M1 = M1, M2 = M2))

# With parallel processing (4 cores in this example)
sim_results <- runSimulation(DESIGNFACTOR, 
                              replications = 500, 
                              generate = gen_lgm_data, 
                              analyse = run_lgm, 
                              summarise = evaluate_lgm, 
                              packages = c("mnormt", "lavaan"), 
                              parallel = TRUE, 
                              ncores = 4L, 
                              fixed_objects = 
                                list(M1 = M1, M2 = M2))
