################################################################################
#
# Reproduce all analyses and figures from the manuscript
#
################################################################################

#!/bin/bash

#$ -N run_all_realdata         # Specify job name
#$ -l m_mem_free=24G           # to request 24 GB of memory per core.
#$ -cwd                        # to run the job in the current working directory.

################ 1. Install R packages #########################################

Rscript -e 'renv::activate(); renv::restore()'

output_dir=$LOCAL_SPACRT_DATA_DIR/private/results/full_data
mkdir -p $output_dir

min_gRNA_count=5
min_cutoff=7
max_cutoff=100
num_subsampled_gens=3000
num_subsampled_gens_timing=2

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

################ 4. Postprocess the results ####################################

echo "Start to postprocess the results!"

Rscript realdata-code/data_postprocessing.R $min_cutoff

echo "Finish the postprocessing!"

################ 5. Create the figures #########################################

echo "Start to make the plots!"

Rscript realdata-code/plotting-code.R $max_cutoff

echo "Finish the plotting code!"

################ 6. Computation time comparison ################################

echo "Start to compare computation time!"

Rscript realdata-code/time_comparison.R $min_gRNA_count $num_subsampled_gens_timing

echo "Finish the computation time comparison!"
