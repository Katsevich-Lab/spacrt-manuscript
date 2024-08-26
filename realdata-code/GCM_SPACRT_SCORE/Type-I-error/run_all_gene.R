# a second attempt for real data analysis with spacrt
library(dplyr)
library(spacrt)
library(MASS)
library(katlabutils)
library(statmod)
library(ondisc)

# Read the command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign the arguments to variables in R
infer_method <- args[1]
nt_grna_group <- args[2]
output_dir <- args[3]
min_gRNA_count <- as.numeric(args[4])

# set the intermediate data directory
intermediate_data_dir <- paste0(.get_config_path("LOCAL_SPACRT_DATA_DIR"),
                                "private/results/full_data/intermediate_data")

# get RNA, gene and covariate dataframes
subsetted_ODM <- readRDS(sprintf("%s/subsetted_ODM.rds",
                                 intermediate_data_dir))
negative_control_list <- readRDS(sprintf("%s/negative_control_pairs.rds",
                                         intermediate_data_dir))
# obtain grna ids
grna_ids <- negative_control_list$grna |>
  dplyr::filter(grna_group == !!nt_grna_group) |>
  dplyr::select(grna_id) |>
  pull()

# take the gRNA, take the threshold to be 5 and extract the result
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

# get the covariate
transform_covariate <- readRDS(sprintf("%s/covariate_matrix.rds",
                                       intermediate_data_dir))

gene_list <- negative_control_list$gene_id
num_gene <- length(gene_list)
# Create empty list
if(infer_method == "spaCRT"){
  output <- matrix(data = NA,
                   nrow = num_gene,
                   ncol = 5,
                   dimnames = list(gene = gene_list,
                                   output = c("left", "right", "both",
                                              "time", "GCM_default")))
}else{
  output <- matrix(data = NA,
                   nrow = num_gene,
                   ncol = 4,
                   dimnames = list(gene = gene_list,
                                   output = c("left", "right", "both",
                                              "time")))
}

# for loop the gene_subsampled_list
dispersion_list <- numeric(num_gene)
names(dispersion_list) <- gene_list
for(gene in gene_list){
  # Get the gene data
  gene_expression <- subsetted_ODM$gene[[gene,]] |>
    as.vector()

  # List the data
  data <- list(
    Y = gene_expression,
    X = grna_presences,
    Z = transform_covariate
  )

  ## Provide the auxiliary information for negative.binomial model
  aux_info_Y_on_Z <- nb_precomp(list(Y = data$Y, Z = data$Z))

  # use one of inference methods
  switch(infer_method,
         GCM = {

           # store the computation time
           output[gene, "time"] <- system.time(
             test_result <- GCM(data = data, X_on_Z_fam = "binomial",
                                Y_on_Z_fam = "negative.binomial")
           )[["elapsed"]]

           # compute and store the p-value
           output[gene, c("left", "right", "both")] <- c(test_result$p.left,
                                                      test_result$p.right,
                                                      test_result$p.both)

         },
         spaCRT = {

           # store the computation time
           output[gene, "time"] <- system.time(
             test_result <- spaCRT(data = data, X_on_Z_fam = "binomial",
                                   Y_on_Z_fam = "negative.binomial")
           )[["elapsed"]]

           # depending on the gcm is rejected or not, we use test_stat or p-value as output
           output[gene, c("left", "right", "both")] <- c(test_result$p.left,
                                                      test_result$p.right,
                                                      test_result$p.both)

           # store the if GCM is performed or not
           output[gene, ncol(output)] <- test_result$gcm.default
         },
         scoreglmnb = {

           # store the computation time
           output[gene, "time"] <- system.time(
             test_result <- score.test(data = data, X_on_Z_fam = "binomial",
                                       Y_on_Z_fam = "negative.binomial")
           )[["elapsed"]]

           # store the dispersion
           dispersion_list[gene] <- test_result$NB.disp.param

           # store the p-value
           output[gene, "left"] <- test_result$p.left
           output[gene, "right"] <- test_result$p.right
           output[gene, "both"] <- test_result$p.both
         }
  )
  print(gene)
  print(output[gene, ])
}

# save the dispersion if it is score test
if(infer_method == "scoreglmnb"){
  dispersion_path <- paste0(.get_config_path("LOCAL_SPACRT_DATA_DIR"),
                            "private/results/full_data/dispersion_list_glm_nb.rds")
  # check if there exists the rds or not; if not, create one
  if (!dir.exists(dispersion_path)) {
    # save the dispersion parameter
    saveRDS(dispersion_list, dispersion_path)
    cat("Result saved:", dispersion_path, "\n")
  } else {
    cat("Results already exists:", dispersion_path, "\n")
  }
}

# save the p-values
saveRDS(output, paste0(output_dir, "/output.rds"))

