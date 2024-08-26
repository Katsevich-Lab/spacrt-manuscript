#!/bin/bash

#$ -N run_all_sceptre          # Specify job name
#$ -l m_mem_free=24G           # to request 8 GB of memory per core.
#$ -cwd                        # to run the job in the current working directory.

source ~/.research_config

# Access the argument passed to the script
directory=$1

# make the output directory
output_dir=$LOCAL_SPACRT_DATA_DIR/private/results/full_data/sceptre
mkdir -p $output_dir

# 1. null simulation
sideness=("left" "right")
for side in "${sideness[@]}"
do
    echo "Rscript $directory/sceptre_null.R "$side" "$output_dir""| qsub -N "sceptre_null"_"$side" -l m_mem_free=12G
done

echo "Submit the job for analyzing real data for sceptre under null!"

# 2. power simulation
grouping=("grouped" "individual")
for grouping_type in "${grouping[@]}"
do
    Rscript -e '
    library(sceptre)
    # set the directory 
    intermediate_data_dir <- paste0(.get_config_path("LOCAL_SPACRT_DATA_DIR"), 
                                    "private/results/full_data/intermediate_data")
    # load the sceptre object 
    sceptre_object <- readRDS(sprintf("%s/sceptre_object_left_%s.rds",
                              intermediate_data_dir, 
                              "'"$grouping_type"'"))
  
    # run power check for sceptre_object
    sceptre_object <- run_power_check(sceptre_object)
  
    # save RDS files
    saveRDS(
      sceptre_object@power_result, 
      paste0("'"$output_dir"'", 
             sprintf("/sceptre_result_left_%s_update_power.rds", 
             "'"$grouping_type"'"))
    )'
done

echo "Finish the job for analyzing real data for sceptre under alternative!"
