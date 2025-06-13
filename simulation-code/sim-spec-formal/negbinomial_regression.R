#####################
# 1. Parameter grid
#####################
get_ground_truth <- function(n){
  return(NULL)
}
varying_params <- list(n = 5000,
                       gamma_0 = seq.int(from = -6, to = -2),
                       beta_0 = seq.int(from = -6, to = -2))

baseline_params <- list(n = 5e3, gamma_0 = -5, beta_0 = -5)

# specify the parameter grid with the varying_params and basline_params
grid_ffd <- simulatr::create_param_grid_fractional_factorial(varying_params, baseline_params)

# specify the signal strength
rho <- c(seq.int(from = -4, to = 0), seq.int(from = 1, to = 4) / 2)

# the parameter grids
parameter_grid <- merge(expand.grid(theta = c(0.05, 1, 10), rho = rho), grid_ffd) |>
  dplyr::mutate(grid_id = 1:dplyr::n()) |>
  simulatr::add_ground_truth(get_ground_truth)

# save the parameter grid
saveRDS(parameter_grid, "simulation-code/sim-spec-formal/NB_parameter_grid.rds")

#######################
# 2. Fixed parameters
#######################
fixed_parameters <- list(
  B = 10000,                    # number of data realizations
  seed = 10,                    # seed to set prior to generating data and running methods
  gamma_1 = 1,                  # slope in X on Z model
  beta_1 = 1                    # slope in Y on Z model
)

###################
# 3. Generate data
###################
# define data-generating model based on the Gaussian linear model
generate_data_nb <- function(n, gamma_0, gamma_1, beta_0, beta_1, rho, theta){

  expit <- function(theta)(exp(theta)/(1 + exp(theta)))

  Z <- as.matrix(rnorm(n = n, mean = 0, sd = 1))
  X <- stats::rbinom(n = n, size = 1, prob = expit(gamma_0 + gamma_1*Z))
  Y <- MASS::rnegbin(n = n, mu = exp(beta_0 + beta_1 * Z + rho*X), theta = theta)

  return(list(X = X, Y = Y, Z = Z))
}

# compile the data generating funciton
generate_data_function <- simulatr::simulatr_function(
  f = generate_data_nb,
  arg_names = formalArgs(generate_data_nb),
  loop = TRUE
)

######################
# 4. Method functions
######################
# GCM method
GCM_spec_f <- simulatr::simulatr_function(f = function(data) {

  # store computation time
  computation_time <- system.time({
    GCM_results <- spacrtutils::GCM_internal(data,
                                             X_on_Z_fam = "binomial",
                                             Y_on_Z_fam = "negative.binomial",
                                             fitting_X_on_Z = "glm",
                                             fitting_Y_on_Z = "glm")
  })["elapsed"]

  # output p-values
  return(list(
    pvalue_list = GCM_results,
    computation_time = computation_time
  ))
}, arg_names = character(0), loop = T)

# dCRT method
dCRT_spec_f <- simulatr::simulatr_function(f = function(data) {

  # store computation time
  computation_time <- system.time({
    dCRT_results <- spacrtutils::dCRT_internal(data,
                                               X_on_Z_fam = "binomial",
                                               Y_on_Z_fam = "negative.binomial",
                                               fitting_X_on_Z = "glm",
                                               fitting_Y_on_Z = "glm",
                                               B = 10000)
  })["elapsed"]

  # output p-values
  return(list(
    pvalue_list = dCRT_results,
    computation_time = computation_time
  ))
}, arg_names = character(0), loop = T)

# spaCRT method
spaCRT_spec_f <- simulatr::simulatr_function(f = function(data) {

  # store computation time
  computation_time <- system.time({
    spaCRT_results <- spacrtutils::spaCRT_internal(data,
                                                   X_on_Z_fam = "binomial",
                                                   Y_on_Z_fam = "negative.binomial",
                                                   fitting_X_on_Z = "glm",
                                                   fitting_Y_on_Z = "glm")
  })["elapsed"]

  # output p-values
  return(list(
    pvalue_list = spaCRT_results,
    computation_time = computation_time
  ))
}, arg_names = character(0), loop = T)

# score test method
scoretest_glm_nb_spec_f <- simulatr::simulatr_function(f = function(data) {

  # store computation time
  computation_time <- system.time({
    score_results <- spacrtutils::score.test(data,
                                             X_on_Z_fam = "binomial",
                                             Y_on_Z_fam = "negative.binomial")
  })["elapsed"]

  # output p-values
  return(list(
    pvalue_list = score_results,
    computation_time = computation_time
  ))
}, arg_names = character(0), loop = T)

# compile the method funcitons
run_method_functions <- list(GCM = GCM_spec_f,
                             dCRT = dCRT_spec_f,
                             spaCRT = spaCRT_spec_f,
                             score.test = scoretest_glm_nb_spec_f)

# aggregate all the necessary components
simulatr_spec <- simulatr::simulatr_specifier(parameter_grid,
                                              fixed_parameters,
                                              generate_data_function,
                                              run_method_functions)

# ############################### check results ##################################
# check_results_nofit <- simulatr::check_simulatr_specifier_object(simulatr_spec, B_in = 1)
