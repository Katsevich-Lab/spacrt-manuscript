library(sceptre)
library(dplyr)

# Read the command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign the arguments to variables in R
side <- args[1]
output_dir <- args[2]

# set the intermediate data directory
intermediate_data_dir <- paste0(.get_config_path("LOCAL_SPACRT_DATA_DIR"), 
                                "private/results/full_data/intermediate_data")


# null simulation
################################################################################
# load the sceptre object 
sceptre_object <- readRDS(sprintf("%s/sceptre_object_%s_grouped.rds",
                                  intermediate_data_dir, 
                                  side))

# run_discovery_analysis for sceptre_object
sceptre_object <- run_discovery_analysis(sceptre_object)

# select the n_nonzero_trt, n_nonzero_ctrl, p-value, response_id and grna_target
sceptre_discovery_result <- sceptre_object@discovery_result |>
  dplyr::select(response_id, 
                grna_target, 
                n_nonzero_trt, 
                n_nonzero_cntrl, 
                p_value)

# save discovery result
saveRDS(tibble::as_tibble(sceptre_discovery_result), 
        paste0(output_dir, 
               sprintf("/sceptre_result_%s_update_discovery.rds", side)))


