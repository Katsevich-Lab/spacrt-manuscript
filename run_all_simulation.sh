################################################################################
################################################################################
############ Reproduce all analyses and figures from the manuscript ############
################################################################################
################################################################################

################ 0. Setup specific to computing environment ####################
# specify the directory that results files will be written to
RESULTS_DIR="$LOCAL_SPACRT_DATA_DIR/private/results"
# specify the preferred number of GB for each simulation process
MAX_GB=16
# specify the preferred running time for each simulation process
MAX_HOURS=4
# # maybe specify other Nextflow options?
PROFILE="standard"


################ 1. Install R packages #########################################

Rscript -e 'renv::activate(); renv::restore()'

################ 2. Run the numerical simulations ##############################
sim_spec_dir=simulation-code/sim-spec-formal

for file in $sim_spec_dir/*.R; do
  sim_name=$(basename "$file" .R)
  output_dir="$RESULTS_DIR/simulation-results/formal-results/$sim_name"
  output_filename="${sim_name}_results.rds"

  echo "$output_dir/$output_filename"

  if [ ! -f "$output_dir/$output_filename" ]; then
    mkdir -p $output_dir
    bash simulation-code/run_simulation.sh \
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


################ 3. Create the figures #########################################

Rscript -e 'source("simulation-code/plotting-code/assemble-plots-NB-disp-5e-2.R")'
Rscript -e 'source("simulation-code/plotting-code/assemble-plots-NB-disp-1.R")'
Rscript -e 'source("simulation-code/plotting-code/assemble-plots-NB-disp-10.R")'















