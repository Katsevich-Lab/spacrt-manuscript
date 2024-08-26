########################################################################
# Run a simulation locally
#
# This R script can be used to check a simulatr specifier object and/or
# to run a small-scale simulation on your local machine. You only need
# to update step 0 with the name of your simulation and whether you're
# checking or running the simulation.
#
# PARAMETERS:
#
# - sim_name: The name of the simulatr specifier file under code/sim_spec. Note
#             that this script will generate the simulatr specifier file from
#             the underlying .R file in code/sim_spec, if it does not already exist.
# - check:    Whether you are checking the simulatr specifier object (this is
#             the primary use case of this script), in which case to do only a
#             small number of replicates, or running the whole simulation, in
#             which case the number of replicates is read from the simulatr
#             specifier object.
#
# RUNNING THE SCRIPT:
#
# To run the script, add values for sim_name and check at the top of the script
# and then run the following command from your R console:
#
# source("check_simulation_local.R")
#
# OUTPUT:
#
# - The simulatr specifier object will be written to code/sim_spec/<sim_name>.rds.
# - If check = FALSE, the results of the simulation will be written to
#   results/<sim_name>/<sim_name>_results.rds.
########################################################################

# 0. Specify the name of the simulation and whether you're checking the simulatr
#    specifier object.
sim_name <- "ToE_grid_explore_FFD_bin_NB_normal_B_5000_5e3_n5_n5_disp_1_L_R_updated"
check <- TRUE

# 1. Set working directory to top level
# setwd(paste0(.get_config_path("LOCAL_CODE_DIR"), "spacrt-project"))

# 2. Set up R packages
# renv::activate()
# renv::restore()
library(simulatr)

# 2. Create and/or load simulatr specific object
sim_spec_script_fp <- paste0("simulation-code/sim-spec-trial/", sim_name, ".R")
sim_spec_fp <- paste0("simulation-code/sim-spec-trial/sim_spec_", sim_name, ".rds")
if(!file.exists(sim_spec_fp)){
  source(sim_spec_script_fp)
  simulatr_spec <- simulatr_specifier(
    parameter_grid,
    fixed_parameters,
    generate_data_function,
    run_method_functions,
    evaluation_functions
  )
  saveRDS(simulatr_spec, sim_spec_fp)
}
simulatr_spec <- readRDS(sim_spec_fp)

# 3. Run the simulation (or a subset of it)
if(check) B_in <- 5 else B_in <- NULL
results <- check_simulatr_specifier_object(simulatr_spec, B_in = B_in, return_data = FALSE, parallel = FALSE)

# 4. Save the results
if(!check){
  results_dir <- paste0("results/", sim_name)
  if(!dir.exists(results_dir)) dir.create(results_dir)
  saveRDS(results, file = paste0(results_dir, "/", sim_name, "_results.rds"))
}

