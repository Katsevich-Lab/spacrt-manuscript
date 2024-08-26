#######################################################################
#
# Run a simulation
#
# This script is a wrapper around the simulatr Nextflow pipeline, documented at
# https://katsevich-lab.github.io/simulatr/articles/example-remote.html.
#
# PARAMETERS:
#
# - sim_spec_dir: (Required) The directory containing either the simulatr
#                 specifier file or the simulatr specifier object.
# - results_dir:  (Required) The directory where the results of the simulation
# - sim_name:     (Required) The name of the simulatr specifier file.
# - profile:      (Optional) The Nextflow profile to use. The two options
#                 are "local" (interactive mode) and "standard" (batch mode).
# - B_check:      (Optional) See documentation linked above.
# - B:            (Optional) See documentation linked above.
# - max_gb:       (Optional) See documentation linked above.
# - max_hours:    (Optional) See documentation linked above.
#
# RUNNING THE SCRIPT:
#
# To run the script interactively, use the command
#
# bash simulation-code/run_simulation.sh --<param-name> <param-value> ...
#
# For example,
#
# bash simulation-code/run_simulation.sh --sim_spec_dir code/simulation/sim_spec \
#        --sim_name sample_sim --B_check 2
#
# To run the script in batch mode, use the command
#
# echo "bash simulation-code/run_simulation_pipeline.sh --<param-name> <param-value> ..." | qsub -N run_all
#
# OUTPUT:
#
# - The results of the simulation will be written to
#   results/<sim_name>/<sim_name>_results.rds.
#######################################################################

# 0. Read simulation parameters from command line

# set default arguments
profile="standard" # Nextflow profile (default standard)
B_check=3          # Number of replicates to use for benchmarking (default 3)
B=0                # Number of replicates for main simulation (default read from sim spec obj)
max_gb=8           # Maximum number of GB per process (default 7.5)
max_hours=4        # Maximum number of hours per process (default 4)
# mnemonic_name=""   # Nextflow mnemonic name of pipeline (useful for resuming a pipeline run earlier)

# read command line arguments
while [ $# -gt 0 ]; do
    if [[ $1 == "--"* ]]; then
        v="${1/--/}"
        declare "$v"="$2"
        shift
    fi
    shift
done

# check that required arguments are given
if [ -z "$sim_name" ]; then
  echo "Error: The simulation name needs to be given via the\
 command-line argument sim_name."
  exit
fi
if [ -z "$sim_spec_dir" ]; then
  echo "Error: The simulation specifier directory needs to be given via the\
 command-line argument sim_spec_dir."
  exit
fi
if [ -z "$results_dir" ]; then
  echo "Error: The results directory needs to be given via the\
 command-line argument results_dir."
  exit
fi
results_dir=$(realpath $results_dir)

# 1. Set up R packages

# update R packages
# Rscript -e 'renv::activate(); renv::restore()'
# extract the R library path from renv
Rscript -e '
renv::activate();
renv::restore();
write(paste0(.libPaths(), collapse = ":"), file = "temp.txt")'

RENV_R_LIBS_USER=$(cat temp.txt)
rm temp.txt

# export the R library path
export R_LIBS_USER=$RENV_R_LIBS_USER
# write the R library path into nextflow config
echo "env.R_LIBS_USER = \"$RENV_R_LIBS_USER\"" > nextflow.config

# 2. Create simulatr specifier object if necessary

sim_spec_script_fp=$(realpath "simulation-code/sim-spec-formal/${sim_name}.R")
sim_spec_fp=$(realpath "simulation-code/sim-spec-formal/sim_spec_${sim_name}.rds")

Rscript -e '
suppressMessages(suppressWarnings(library(simulatr)))
sim_spec_script_fp <- "'"$sim_spec_script_fp"'"
sim_spec_fp <- "'"$sim_spec_fp"'"
if (!file.exists(sim_spec_fp)) {
  if (!file.exists(sim_spec_script_fp)) {
    stop("Error: The simulatr specifier script was not found.")
  } else {
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
}'

# 3. Run the simulation

output_dir=$results_dir
mkdir -p $output_dir
output_filename=$sim_name"_results.rds"
if [ ! -f "$output_dir/$output_filename" ]; then
  echo "Running the "$sim_name" simulation..."
  nextflow pull katsevich-lab/simulatr-pipeline -r parallel-evaluation
  nextflow run katsevich-lab/simulatr-pipeline -r parallel-evaluation \
    --simulatr_specifier_fp $sim_spec_fp \
    --result_dir $output_dir \
    --result_file_name $output_filename \
    --B_check $B_check \
    --B $B \
    --max_gb $max_gb \
    --max_hours $max_hours \
    -profile $profile \
    -with-trace $output_dir/$sim_name"_trace.txt"
  if [ -f "$output_dir/$output_filename" ]; then
      nextflow clean
  fi
else
  echo $output_filename" already exists"
fi

echo "Done."
