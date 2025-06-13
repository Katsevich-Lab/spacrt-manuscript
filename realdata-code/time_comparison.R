library(dplyr)
library(spacrt)
library(MASS)
library(katlabutils)
library(ondisc)
library(statmod)

# Read the command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign the arguments to variables in R
min_gRNA_count <- as.numeric(args[1])
num_subsampled_gene <- as.numeric(args[2])

# Specify the output file path
output_dir <- sprintf("%s/full_data",
                      .get_config_path("LOCAL_SPACRT_DATA_DIR"))

# set the intermediate data directory
intermediate_data_dir <- paste0(.get_config_path("LOCAL_SPACRT_DATA_DIR"),
                                "full_data/intermediate_data")

# get RNA, gene and covariate dataframes
subsetted_ODM <- readRDS(sprintf("%s/subsetted_ODM.rds",
                                 intermediate_data_dir))
negative_control_list <- readRDS(sprintf("%s/negative_control_pairs.rds",
                                         intermediate_data_dir))
grna_group_list <- negative_control_list$grna$grna_group |> unique()

# get the covariate
transform_covariate <- readRDS(sprintf("%s/covariate_matrix.rds",
                                       intermediate_data_dir))


# run random 2 genes
set.seed(1)
gene_subsampled_list <- sample(negative_control_list$gene_id,
                               num_subsampled_gene, replace = FALSE)

# Create empty list
method_list <- c("GCM", "dCRT", "spaCRT", "score")
output <- array(data = NA,
                dim = c(length(gene_subsampled_list),
                        length(grna_group_list), length(method_list)),
                dimnames = list(gene = gene_subsampled_list,
                                gRNA = grna_group_list, method = method_list))

# for loop the grna_groups
for (grna_group in grna_group_list) {

  # obtain grna ids
  grna_ids <- negative_control_list$grna |>
    dplyr::filter(grna_group == !!grna_group) |>
    dplyr::select(grna_id) |>
    pull()

  # take the gRNA, take the threshold and extract the result
  grna_presences <- subsetted_ODM$grna[[grna_ids, ]] |>
    Matrix::t() |>
    Matrix::as.matrix() |>
    as.data.frame() |>
    stats::setNames(grna_ids) |>
    tibble::as_tibble() |>
    dplyr::transmute(
      dplyr::if_any(.cols = dplyr::everything(),
                    .fns = function(grna_count)(grna_count >= min_gRNA_count))) |>
    dplyr::pull()

  # loop over genes
  for(gene in gene_subsampled_list){

    # Get the gene data
    gene_expression <- subsetted_ODM$gene[[gene, ]] |>
      as.vector()

    # List the data
    data <- list(
      Y = gene_expression,
      X = as.numeric(grna_presences),
      Z = transform_covariate
    )

    # store the computation time
    output[gene, grna_group, "GCM"] <- system.time(
      spacrtutils::GCM_internal(data = data,
                                X_on_Z_fam = "binomial",
                                Y_on_Z_fam = "negative.binomial")
    )[["elapsed"]]

    output[gene, grna_group, "spaCRT"] <- system.time(
      spacrtutils::spaCRT_internal(data = data,
                                   X_on_Z_fam = "binomial",
                                   Y_on_Z_fam = "negative.binomial")
    )[["elapsed"]]

    output[gene, grna_group, "score"] <- system.time(
      spacrtutils::score.test(data = data,
                              X_on_Z_fam = "binomial",
                              Y_on_Z_fam = "negative.binomial")
    )[["elapsed"]]

    output[gene, grna_group, "dCRT"] <- system.time(
      spacrtutils::dCRT_internal(data,
                                 X_on_Z_fam = "binomial",
                                 Y_on_Z_fam = "negative.binomial",
                                 B = 1e5)
    )[["elapsed"]]
  }
}

# Convert array to a long format dataframe
output_df <- as.data.frame(as.table(output))

# Compute mean and standard deviation for each method
results_df <- output_df |>
  group_by(method) |>
  summarise(
    time_mean = mean(Freq, na.rm = TRUE),
    time_sd = sd(Freq, na.rm = TRUE)
  ) |>
  ungroup()

# save the output
saveRDS(results_df, paste0(output_dir, "/time_comparison.rds"))

