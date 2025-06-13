suppressMessages(suppressWarnings(library(simulatr)))
suppressMessages(suppressWarnings(library(spacrt)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(MASS)))

#####################
# 1. Parameter grid
#####################
varying_params <- list(n = 5000,
                       gamma_0 = seq.int(from = -6, to = -2),
                       beta_0 = seq.int(from = -6, to = -2))

baseline_params <- list(n = 5e3, gamma_0 = -5, beta_0 = -5)

grid_ffd <- simulatr::create_param_grid_fractional_factorial(varying_params,
                                                             baseline_params)

get_ground_truth <- function(n){
  return(NULL)
}

rho <- c(seq.int(from = -4, to = 0), seq.int(from = 1, to = 4)/2)

parameter_grid <- grid_ffd %>% dplyr::select(1:3) %>%
  slice(rep(1:n(), each = length(rho))) %>%
  cbind(rho = rep(rho, times = nrow(grid_ffd))) %>%
  add_ground_truth(get_ground_truth)


#######################
# 2. Fixed parameters
#######################
fixed_parameters <- list(
  B = 50000,                      # number of data realizations
  seed = 10,                    # seed to set prior to generating data and running methods
  gamma_1 = 1,               # slope in X on Z model
  beta_1 = 1,                # slope in Y on Z model
  theta = 1                  # dispersion parameter
)


###################
# 3. Generate data
###################

# define data-generating model based on the Gaussian linear model
generate_data_nb <- function(n, gamma_0, gamma_1,
                             beta_0, beta_1, rho, theta){

  expit <- function(theta)(exp(theta)/(1 + exp(theta)))

  Z <- as.matrix(rnorm(n = n, mean = 0, sd = 1))
  X <- rbinom(n = n, size = 1, prob = expit(gamma_0 + gamma_1*Z))
  Y <- MASS::rnegbin(n = n, mu = exp(beta_0 + beta_1*Z + rho*X), theta = theta)

  return(list(X = X, Y = Y, Z = Z))
}


generate_data_function <- simulatr_function(
  f = generate_data_nb,
  arg_names = formalArgs(generate_data_nb),
  loop = TRUE
)

######################
# 4. Method functions
######################

GCM_spec_f <- simulatr_function(f = function(data) spacrt::GCM(data, "binomial","negative.binomial"),
                                arg_names = character(0), loop = T)
dCRT_spec_f <- simulatr_function(f = function(data) spacrt::dCRT(data, "binomial","negative.binomial",
                                                                 B = 10000,FALSE,FALSE),
                                 arg_names = character(0), loop = T)
spaCRT_spec_f <- simulatr_function(f = function(data) spacrt::spaCRT(data, "binomial","negative.binomial",
                                                                     FALSE,FALSE,5),
                                   arg_names = character(0), loop = T)
scoretest_glm_nb_spec_f <- simulatr_function(f = function(data) spacrt::score.test(data, "binomial",
                                                                                   "negative.binomial"),
                                             arg_names = character(0), loop = T)


run_method_functions <- list(GCM = GCM_spec_f,
                             dCRT = dCRT_spec_f,
                             spaCRT = spaCRT_spec_f,
                             score.test = scoretest_glm_nb_spec_f)

##########################
# 5. Evaluation functions
##########################

level.1.L <- function(output, ground_truth){
  return(suppressWarnings(mean(output$p.left, na.rm = TRUE) <= 5e-2))
}

level.1.R <- function(output, ground_truth){
  return(suppressWarnings(mean(output$p.right, na.rm = TRUE) <= 5e-2))
}

level.2.L <- function(output, ground_truth){
  return(suppressWarnings(mean(output$p.left, na.rm = TRUE) <= 1e-2))
}

level.2.R <- function(output, ground_truth){
  return(suppressWarnings(mean(output$p.right, na.rm = TRUE) <= 1e-2))
}

level.3.L <- function(output, ground_truth){
  return(suppressWarnings(mean(output$p.left, na.rm = TRUE) <= 5e-3))
}

level.3.R <- function(output, ground_truth){
  return(suppressWarnings(mean(output$p.right, na.rm = TRUE) <= 5e-3))
}

evaluation_functions <- list(level.1.L = level.1.L,
                             level.1.R = level.1.R,
                             level.2.L = level.2.L,
                             level.2.R = level.2.R,
                             level.3.L = level.3.L,
                             level.3.R = level.3.R)
