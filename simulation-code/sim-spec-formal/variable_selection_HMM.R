library(glmnet)
library(simulatr)
library(spacrtutils)
library(SNPknock)

#####################
# 1. Parameter grid
#####################

get_ground_truth <- function(n){
  return(NULL)
}

# create parameter grid
parameter_grid <- expand.grid(eta = seq(0, 1, length.out = 5),
                              gamma_0 = seq(-3, -2, length.out = 2),
                              beta_param = list(skew = stats::setNames(c(1, 3), c("alpha", "beta")),
                                                uniform = stats::setNames(c(1, 1), c("alpha", "beta")))) |>
  add_ground_truth(get_ground_truth) |>
  dplyr::mutate(grid_id = 1:dplyr::n())


#######################
# 2. Fixed parameters
#######################
fixed_parameters <- list(
  lasso_model = c("min", "1se"), # model selection parameter
  signal_prop = 0.1,             # signal proportion in
  alpha = 0.1,                   # FDR control required for knockoff
  K = 10,                        # number of latent states
  M = 2,                         # output distribution
  p = 500,                       # dimension
  n = 2000,                      # sample size
  B = 50,                        # number of data realizations
  seed = 2                       # seed to set prior to generating data and running methods
)

# save the parameter_grid
saveRDS(parameter_grid |>
          dplyr::mutate(
            p = fixed_parameters$p,
            signal_prop = fixed_parameters$signal_prop,
            alpha = fixed_parameters$alpha,
            B = fixed_parameters$B
          ), "simulation-code/sim-spec-formal/HMM_parameter_grid.rds")

###################
# 3. Generate data
###################

# define data-generating model based gamma
generate_data <- function(n, p, M, K, gamma_0, alpha,
                          eta, signal_prop, lasso_model, beta_param) {

  # generate pEmit, Q and pInit using spacrt package
  aux_info <- withr::with_seed(1, spacrtutils::help_dgp(p = p, K = K, M = M,
                                                        stay_prob = .9, beta_prior = TRUE,
                                                        alpha = beta_param["alpha"],
                                                        beta = beta_param["beta"]))

  # sample from HMM
  X <- SNPknock::sampleHMM(aux_info$pInit, aux_info$Q, aux_info$pEmit, n = n)

  # create effect size
  num_signal <- round(signal_prop * p)
  positive_signal <- withr::with_seed(1, runif(n = num_signal / 2, min = eta, max = eta))
  negative_signal <- withr::with_seed(1, runif(n = num_signal / 2, min = -eta, max = -eta))
  effect_size <- c(positive_signal, negative_signal)

  # sample high dimensional logistic model for Y
  Y <- rbinom(n = n, size = 1,
              prob = 1 / (1 + exp(-gamma_0 - colSums(t(X[, 1:num_signal]) * effect_size))))

  # transform X to inp file and obtain the file path
  hashing_id <- sprintf("eta%.3f-gamma%d-beta%d-%s",
                        eta, gamma_0, beta_param["beta"],
                        stringi::stri_rand_strings(n = 1, length = 20))
  X_file <- paste0(.get_config_path("LOCAL_CODE_DIR"),
                   sprintf("spacrt-project/simulation-code/sim-spec-formal/intermediate-files/%s-X.inp", hashing_id))
  SNPknock::writeXtoInp(X, phased=TRUE, out_file = X_file)

  # obtain the file path to the fastPhase
  fp_path <- paste0(.get_config_path("LOCAL_CODE_DIR"), "fastPhase/fastPHASE")

  # obtain the file path to the output path where the results will be saved
  out_path <- paste0(.get_config_path("LOCAL_CODE_DIR"),
                     "spacrt-project/simulation-code/sim-spec-formal/fastPhase-results")

  # return data
  return(list(X = X, Y = Y,
              X_file = X_file, fp_path = fp_path, out_path = out_path,
              hashing_id = hashing_id, alpha = alpha,
              pInit = aux_info$pInit, pEmit = aux_info$pEmit, Q = aux_info$Q,
              lasso_model = lasso_model))
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

GCM_HMM_f <- function(data){

  # save the time
  computation_time <- system.time({
    data$out_path <- paste0(data$out_path, sprintf("/GCM/%s", data$hashing_id))
    GCM_results <- spacrtutils::GCM_HMM(data)
  })["elapsed"]

  # output p-values
  return(list(
    pvalue_list = GCM_results,
    computation_time = computation_time
  ))
}

dCRT_HMM_f <- function(data){

  # save the computation time
  computation_time <- system.time({
    data$out_path <- paste0(data$out_path, sprintf("/dCRT/%s", data$hashing_id))
    data$B <- 5000
    dCRT_results <- spacrtutils::dCRT_HMM(data)
  })["elapsed"]

  # output sorted p-values
  return(list(
    pvalue_list = dCRT_results,
    computation_time = computation_time
  ))
}

spaCRT_HMM_f <- function(data){
  # save the computation time
  computation_time <- system.time({
    data$out_path <- paste0(data$out_path, sprintf("/spaCRT/%s", data$hashing_id))
    spaCRT_results <- spacrtutils::spaCRT_HMM(data)
  })["elapsed"]

  # output sorted p-values
  return(list(
    pvalue_list = spaCRT_results,
    computation_time = computation_time
  ))
}

knockoff_HMM_f <- function(data){

  # save the computation time
  computation_time <- system.time({
    data$out_path <- paste0(data$out_path, sprintf("/knockoff/%s", data$hashing_id))
    knockoff_results <- spacrtutils::knockoff_HMM(data)
  })["elapsed"]

  # return the selection set indices
  return(list(
    rejection_list = knockoff_results,
    computation_time = computation_time
  ))
}

# assemble the method functions
GCM_HMM_spec_f <- simulatr_function(f = GCM_HMM_f, arg_names = character(0), loop = T)
dCRT_HMM_spec_f <- simulatr_function(f = dCRT_HMM_f, arg_names = character(0), loop = T)
spaCRT_HMM_spec_f <- simulatr_function(f = spaCRT_HMM_f, arg_names = character(0), loop = T)
knockoff_HMM_spec_f <- simulatr_function(f = knockoff_HMM_f, arg_names = character(0), loop = T)

# add all the functions
run_method_functions <- list(
  GCM = GCM_HMM_spec_f,
  dCRT = dCRT_HMM_spec_f,
  spaCRT = spaCRT_HMM_spec_f,
  Knockoff = knockoff_HMM_spec_f
)

# aggregate all the necessary components
simulatr_spec <- simulatr_specifier(parameter_grid,
                                    fixed_parameters,
                                    generate_data_function,
                                    run_method_functions)

# ############################### check results ##################################
# check_results_nofit <- check_simulatr_specifier_object(simulatr_spec, B_in = 1)
