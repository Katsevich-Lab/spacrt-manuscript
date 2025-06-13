################################################################################
################################################################################
############ Reproduce all analyses and figures from the manuscript ############
################################################################################
################################################################################

# Parse command line arguments
MODE="full"  # default mode
if [[ $1 == "plot-only" ]]; then
  MODE="plot-only"
elif [[ $# -gt 0 ]]; then
  echo "Usage: $0 [plot-only]"
  echo "  plot-only: Skip simulations, create figures from existing results"
  echo "  (no args): Run complete simulation pipeline"
  exit 1
fi

echo "Running in $MODE mode"

################ 0. Setup specific to computing environment ####################
# specify the directory that results files will be written to
RESULTS_DIR="$LOCAL_SPACRT_DATA_DIR"
# specify the preferred number of GB for each simulation process
MAX_GB=8
# specify the preferred running time for each simulation process
MAX_HOURS=4
# specify the profile for the simulation
PROFILE="standard"

################ 1. Install R packages #########################################

Rscript -e 'renv::activate(); renv::restore()'

if [[ "$MODE" == "plot-only" ]]; then
  echo "Plot-only mode: Skipping simulations, creating figures from existing results"
  
  # Check if simulation results exist
  if [[ ! -d "$RESULTS_DIR/formal-simulation-results" ]]; then
    echo "Error: Simulation results directory not found."
    echo "Please either:"
    echo "  1. Download simulation results from Dropbox and place in $RESULTS_DIR"
    echo "  2. Run with no arguments to generate results"
    exit 1
  fi
  
elif [[ "$MODE" == "full" ]]; then
  echo "Full mode: Running complete simulation pipeline"
  
  ################ 2. Run the numerical simulations ##############################
  # specify sim_spec_dir
  sim_spec_dir=simulation-code/sim-spec-formal
  
  # make two intermediate directories
  mkdir -p $sim_spec_dir/intermediate-files
  mkdir -p $sim_spec_dir/fastPhase-results
  
  # loop over different simulations
  for file in $sim_spec_dir/*.R; do
    sim_name=$(basename "$file" .R)
    output_dir="$RESULTS_DIR/formal-simulation-results/$sim_name"
    output_filename="${sim_name}_results.rds"
  
    echo "$output_dir/$output_filename"
  
    if [ ! -f "$output_dir/$output_filename" ]; then
      mkdir -p $output_dir
      bash simulation-code/run_new_simulation.sh \
        --sim_spec_dir $sim_spec_dir \
        --results_dir $output_dir \
        --sim_name $sim_name \
        --max_gb $MAX_GB \
        --max_hours $MAX_HOURS \
        --profile $PROFILE
      wait
      echo -e 'Task Complete.\nVery Good!\n....\n....\n'
    else
      echo "$output_filename already exists"
    fi
  done
  
  ################# 3. Delete the intermediate directories #######################
  rm -rf $sim_spec_dir/intermediate-files
  rm -rf $sim_spec_dir/fastPhase-results
  
fi

################ 4. Create the figures #########################################

echo "Creating simulation figures..."

Rscript -e 'source("simulation-code/plotting-code/NB_plotting.R")'
Rscript -e 'source("simulation-code/plotting-code/HMM_plotting.R")'
Rscript -e 'source("simulation-code/plotting-code/RF_plotting.R")'

echo "Simulation figures complete!"





