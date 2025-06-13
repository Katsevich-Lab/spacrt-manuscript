library(simulatr)
library(spacrtutils)

#####################
# 1. Parameter grid
#####################
get_ground_truth <- function(n){
  return(NULL)
}

# create parameter grid
parameter_grid <- data.frame(eta = seq(0, 4, length.out = 5),
                             gamma_0 = -2.5) |>
  dplyr::bind_rows(data.frame(eta = seq(0, 2, length.out = 5),
                              gamma_0 = -2)) |>
  add_ground_truth(get_ground_truth) |>
  dplyr::mutate(grid_id = 1:dplyr::n())

#######################
# 2. Fixed parameters
#######################
fixed_parameters <- list(
  p = 3,                         # dimension
  n = 500,                       # sample size
  B = 50000,                     # number of data realizations
  seed = 2                       # seed to set prior to generating data and running methods
)

# save the parameter_grid
saveRDS(parameter_grid |>
          dplyr::mutate(
            p = fixed_parameters$p,
            alpha = fixed_parameters$alpha,
            B = fixed_parameters$B
          ), "simulation-code/sim-spec-formal/random_forest_parameter_grid.rds")

###################
# 3. Generate data
###################

# define data-generating model based gamma
generate_data <- function(n, p, gamma_0, eta){

  # sample covariate from uniform distribution
  Z <- matrix(runif(n * p, -1, 1), nrow = n, ncol = p)

  # prediction function
  nonlinear_pred <- function(z){
    # fill the details
    return(sin(pi * z[1]) + z[2]^2 * z[3])
  }

  # sample nonlinear model for X
  X <- rbinom(n = n, size = 1, prob = 1/(1 + exp(-gamma_0 + apply(Z, 1, nonlinear_pred))))

  # sample nonlinear model for Y
  Y <- rbinom(n = n, size = 1, prob = 1/(1 + exp(-gamma_0 + eta * X + apply(Z, 1, nonlinear_pred))))

  # return data
  return(list(X = X, Y = Y, Z = Z))
}

# need to call simulatr_function() to give simulatr a few more pieces of info
generate_data_function <- simulatr_function(
  f = generate_data,
  arg_names = formalArgs(generate_data),
  loop = TRUE
)

######################
# 4. Method functions
######################

GCM_rf_f <- function(data){

  # save the time
  computation_time <- system.time({
    GCM_outputs <- spacrtutils::GCM_internal(data = data,
                                             X_on_Z_fam = "binomial",
                                             Y_on_Z_fam = "binomial",
                                             fitting_X_on_Z = 'random_forest',
                                             fitting_Y_on_Z = 'random_forest')
  })["elapsed"]

  # output p-values
  return(list(
    pvalue_list = GCM_outputs,
    computation_time = computation_time
  ))
}

dCRT_rf_f <- function(data, B = 50000){

  # save the computation time
  computation_time <- system.time({
    dCRT_outputs <- spacrtutils::dCRT_internal(data = data,
                                               X_on_Z_fam = "binomial",
                                               Y_on_Z_fam = "binomial",
                                               fitting_X_on_Z = 'random_forest',
                                               fitting_Y_on_Z = 'random_forest',
                                               B = B)
  })["elapsed"]

  # output sorted p-values
  return(list(
    pvalue_list = dCRT_outputs,
    computation_time = computation_time
  ))
}

spaCRT_rf_f <- function(data){
  # save the computation time
  computation_time <- system.time({
    spaCRT_outputs <- spacrtutils::spaCRT_internal(data = data,
                                                   X_on_Z_fam = "binomial",
                                                   Y_on_Z_fam = "binomial",
                                                   fitting_X_on_Z = 'random_forest',
                                                   fitting_Y_on_Z = 'random_forest')
  })["elapsed"]

  # output sorted p-values
  return(list(
    pvalue_list = spaCRT_outputs,
    computation_time = computation_time
  ))
}


# assemble the method functions
GCM_rf_spec_f <- simulatr_function(f = GCM_rf_f, arg_names = character(0), loop = T)
dCRT_rf_spec_f <- simulatr_function(f = dCRT_rf_f, arg_names = character(0), loop = T)
spaCRT_rf_spec_f <- simulatr_function(f = spaCRT_rf_f, arg_names = character(0), loop = T)

# add all the functions
run_method_functions <- list(
  GCM = GCM_rf_spec_f,
  dCRT = dCRT_rf_spec_f,
  spaCRT = spaCRT_rf_spec_f
)

# aggregate all the necessary components
simulatr_spec <- simulatr_specifier(parameter_grid,
                                    fixed_parameters,
                                    generate_data_function,
                                    run_method_functions)

# ############################### check results ##################################
# check_results_nofit <- check_simulatr_specifier_object(simulatr_spec, B_in = 1)
