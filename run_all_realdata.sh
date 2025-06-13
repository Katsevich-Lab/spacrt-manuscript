################################################################################
#
# Reproduce all analyses and figures from the manuscript
#
################################################################################

#!/bin/bash

#$ -N run_all_realdata         # Specify job name
#$ -l m_mem_free=24G           # to request 24 GB of memory per core.
#$ -cwd                        # to run the job in the current working directory.

# Parse command line arguments
MODE="full"  # default mode
if [[ $1 == "plot-only" ]]; then
  MODE="plot-only"
elif [[ $# -gt 0 ]]; then
  echo "Usage: $0 [plot-only]"
  echo "  plot-only: Skip analysis, create figures from existing results"
  echo "  (no args): Run complete analysis pipeline"
  exit 1
fi

echo "Running in $MODE mode"

################ 1. Install R packages #########################################

Rscript -e 'renv::activate(); renv::restore()'

output_dir=$LOCAL_SPACRT_DATA_DIR/full_data
min_gRNA_count=5
min_cutoff=7
max_cutoff=100
num_subsampled_gens=3000
num_subsampled_gens_timing=2

if [[ "$MODE" == "plot-only" ]]; then
  echo "Plot-only mode: Skipping analysis, creating figures from existing results"
  
  # Check if results exist
  if [[ ! -d "$output_dir" ]]; then
    echo "Error: Results directory $output_dir not found."
    echo "Please either:"
    echo "  1. Download results from Dropbox and place in $output_dir"
    echo "  2. Run with --mode full to generate results"
    exit 1
  fi
  
elif [[ "$MODE" == "full" ]]; then
  echo "Full mode: Running complete analysis pipeline"
  
  mkdir -p $output_dir
  
  ################ 2. Create sceptre object and intermediate data ################
  Rscript realdata-code/data_preprocessing.R $min_gRNA_count $num_subsampled_gens
  
  echo "Finish constructing sceptre objects and data extraction!"
  
  ################ 3. Run the sceptre simulation #################################
  
  echo "Begin SCEPTRE analysis!"
  
  # run sceptre analysis
  sceptre_directory=realdata-code/SCEPTRE
  qsub $sceptre_directory/sceptre.sh $sceptre_directory
  
  ################ 3. Run the three methods simulation ###########################
  
  echo "Start to run simulation for GCM, spaCRT and score!"
  
  # run the null simulation
  type_I_err_three_methods=realdata-code/GCM_SPACRT_SCORE/Type-I-error
  qsub $type_I_err_three_methods/null_three_methods.sh $type_I_err_three_methods $min_gRNA_count -N "Type_I_error_for_three_methods"
  
  echo "Finish the negative control analysis for GCM, spaCRT and score!"
  
  # run the power simulation
  echo "Rscript realdata-code/GCM_SPACRT_SCORE/power_three_methods.R "$min_gRNA_count""| qsub -N "power_three_methods"
  
  echo "Submit the positive control analysis for GCM, spaCRT and score!"
  
  ################ 4. Check if all the simulation finish #########################
  
  # check the latest submitted job finish or not; it does not complete, wait for it
  bash realdata-code/check.sh
  
  echo "Finish the analysis for all the methods!"
  
  ################ 5. Postprocess the results ####################################
  
  echo "Start to postprocess the results!"
  
  Rscript realdata-code/data_postprocessing.R $min_cutoff
  
  echo "Finish the postprocessing!"
  
  ################ 6. Computation time comparison ################################
  
  echo "Start to compare computation time!"
  
  Rscript realdata-code/time_comparison.R $min_gRNA_count $num_subsampled_gens_timing
  
  echo "Finish the computation time comparison!"
  
fi

################ 7. Create the figures #########################################

echo "Creating figures and tables..."

Rscript realdata-code/plotting-code.R $max_cutoff

# Create Table 3
Rscript realdata-code/sparsity_dataset.R

echo "Figures and tables complete!"
