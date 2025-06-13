source ~/.research_config

# Access the argument passed to the script
directory=$1
min_gRNA_count=$2

# 1. Extract the gRNA names from data

# write the gRNA list to gRNA.txt
Rscript -e '
intermediate_data_dir <- paste0(.get_config_path("LOCAL_SPACRT_DATA_DIR"), 
                                "full_data/intermediate_data")
negative_control_list <- readRDS(sprintf("%s/negative_control_pairs.rds", 
                                         intermediate_data_dir))
grna_group <- negative_control_list$grna$grna_group |> unique()
write(grna_group, file = "gRNA.txt", sep = "\n")
'

# load the gRNA.txt to the bash environment
gRNA=$(cat gRNA.txt)
rm gRNA.txt
IFS=$'\n' read -r -d '' -a RNA_list <<< "$gRNA"

# 2. For loop the gRNA names; save the output

# create the output folder for method
methods=("GCM" "spaCRT" "scoreglmnb")
for method in "${methods[@]}"
do
    output_dir=$LOCAL_SPACRT_DATA_DIR/full_data/$method
    mkdir -p $output_dir
    for RNA in "${RNA_list[@]}"
    do
    # create the output folder for gRNA
    output_gRNA_dir=$output_dir/$RNA
    mkdir -p $output_gRNA_dir
    
    # run the R script for each gRNA
    echo "Rscript $directory/run_all_gene.R "$method" "$RNA" "$output_gRNA_dir" "$min_gRNA_count""| qsub -N "$method"_"$RNA"
    done
    echo "Submit the null simulation for $method test!"
done

