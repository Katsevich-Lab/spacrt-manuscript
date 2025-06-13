library(ondisc)
library(sceptre)
library(dplyr)

# Read the command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign the arguments to variables in R
min_gRNA_count <- as.numeric(args[1])
num_gene <- as.numeric(args[2])

# load the data
gasperini_dir <- paste0(.get_config_path("LOCAL_GASPERINI_2019_V2_DATA_DIR"), "at-scale/processed/")

# 0. set the Gasperini directory load the group, gene, and grna data
# read the data in ondisc format
gasp_fp <- paste0(.get_config_path("LOCAL_GASPERINI_2019_V2_DATA_DIR"),
                  "at-scale/processed/")
gene_odm_fp <- paste0(gasp_fp, "gene/matrix.odm")
gene_metadata_fp <- paste0(gasp_fp, "gene/metadata.rds")
grna_odm_fp <- paste0(gasp_fp, "grna_expression/matrix.odm")
grna_metadata_fp <-  paste0(gasp_fp, "grna_expression/metadata.rds")

# read gene and grna ondisc matrices
gene_odm <- ondisc::read_odm(odm_fp = gene_odm_fp,
                             metadata_fp = gene_metadata_fp)
grna_odm <- ondisc::read_odm(odm_fp = grna_odm_fp,
                             metadata_fp = grna_metadata_fp)

# remove cells with 0 grna expression or gene expression
multimodal_odm <- ondisc::multimodal_ondisc_matrix(covariate_ondisc_matrix_list =
                                                     list(grna = grna_odm, gene = gene_odm))
ok_cells_grna <- (grna_odm |> ondisc::get_cell_covariates() |>
                    dplyr::pull(n_umis)) != 0L
ok_cells_gene <- (gene_odm |> ondisc::get_cell_covariates() |>
                    dplyr::pull(n_umis)) != 0L
multimodal_odm_sub <- multimodal_odm[, ok_cells_grna & ok_cells_gene]

# load the grna matrix
grna_mat <- multimodal_odm_sub@modalities$grna[[seq(1, nrow(multimodal_odm_sub@modalities$grna)),]]
rownames(grna_mat) <- multimodal_odm_sub@modalities$grna |> ondisc::get_feature_ids()

# get the global cell covariates
cell_covariates <- multimodal_odm_sub |>
  ondisc::get_cell_covariates() |>
  dplyr::select(p_mito = gene_p_mito, batch = gene_batch)
rownames(cell_covariates) <- NULL

# get the grna data frame
grna_feature_df <- grna_odm |>
  ondisc::get_feature_covariates() |>
  tibble::rownames_to_column(var = "grna_id") |>
  dplyr::select(grna_id, grna_group, target_type, target, target_gene) |>
  tibble::as_tibble()
target_type <- as.character(grna_feature_df$target_type)
target_gene <- as.character(grna_feature_df$target_gene)
target <- as.character(grna_feature_df$target)
x <- ifelse(target_type %in% c("candidate_enhancer", "known_enhancer"), target,
            ifelse(target_type == "gene_tss", target_gene, target_type))
x[is.na(x)] <- "not_mapped"

# From here we have to use the values in grna_group to replace the non-targeting
# element in grna_target
#################################################################################
grna_group_data_frame <- data.frame(grna_id = grna_feature_df$grna_id,
                                    grna_group = grna_feature_df$grna_group,
                                    grna_target = x)

# The following step is to avoid the invalidant message from sceptre
grna_group_data_frame <- grna_group_data_frame |>
  dplyr::mutate(
    grna_target = dplyr::if_else(grna_target == "non-targeting",
                                 grna_group,
                                 grna_target)
  )

# get the gene ids, grna ids and grna groups for positive control pairs
gene_full_list <- (multimodal_odm_sub@modalities$gene |> ondisc::get_feature_ids())
positive_control_grna_group_df <- grna_group_data_frame |>
  dplyr::filter(grna_target %in% gene_full_list)

gene_target <- positive_control_grna_group_df |>
  dplyr::select(grna_target) |>
  dplyr::pull() |>
  as.character() |>
  unique()

# get # num_gene random genes
set.seed(1)
genes <- gene_odm |> ondisc::get_feature_ids()
num_sample <- num_gene - length(gene_target)
sample_list <- sample(setdiff(genes, gene_target),
                      num_sample,
                      replace = FALSE)
gene_list <- c(
  sample_list,
  gene_target
  )

# load the gene matrix
gene_mat <- multimodal_odm_sub@modalities$gene[[gene_list,]]
rownames(gene_mat) <- gene_list

# create the sceptre_object
sceptre_object <- import_data(response_matrix = gene_mat,
                              grna_matrix = grna_mat,
                              grna_target_data_frame = grna_group_data_frame,
                              moi = "high",
                              extra_covariates = cell_covariates)

# get lists of NT gRNA groups
grna_nontarget_groups <- grna_feature_df |>
  dplyr::filter(target_type == "non-targeting") |>
  unique()

# specify the grna_target_data_frame
grna_nontarget_data_frame <- grna_nontarget_groups |>
  dplyr::rename(grna_target = grna_group) |>
  dplyr::select(-target_type)

# obtain the response_id versus grna_target pairs
discovery_pairs <- expand.grid(
  response_id = gene_list,
  grna_target = unique(grna_nontarget_groups$grna_group)
) |>
  dplyr::mutate(response_id = as.character(response_id),
                grna_target = as.character(grna_target))

# construct the positive control pairs
positive_control_pairs <- construct_positive_control_pairs(sceptre_object)

# set analysis parameters and run qc
sceptre_object_left_grouped <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  discovery_pairs = discovery_pairs,
  positive_control_pairs = positive_control_pairs,
  grna_integration_strategy = "union",
  side = "left") |>
  assign_grnas(method = "thresholding", threshold = min_gRNA_count) |>
  run_qc(p_mito_threshold = 1,
         response_n_umis_range = c(0, 1),
         response_n_nonzero_range = c(0,1))

sceptre_object_right_grouped <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  discovery_pairs = discovery_pairs,
  positive_control_pairs = positive_control_pairs,
  grna_integration_strategy = "union",
  side = "right") |>
  assign_grnas(method = "thresholding", threshold = min_gRNA_count) |>
  run_qc(p_mito_threshold = 1,
         response_n_umis_range = c(0, 1),
         response_n_nonzero_range = c(0,1))

sceptre_object_left_individual <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  discovery_pairs = discovery_pairs,
  positive_control_pairs = positive_control_pairs,
  grna_integration_strategy = "singleton",
  side = "left") |>
  assign_grnas(method = "thresholding", threshold = min_gRNA_count) |>
  run_qc(p_mito_threshold = 1,
         response_n_umis_range = c(0, 1),
         response_n_nonzero_range = c(0,1))

# positive control list
positive_control_list <- positive_control_grna_group_df |>
  dplyr::filter(grna_target %in% gene_list)

# negative control list
negative_control_list <- list(
  grna = grna_nontarget_groups |>
    dplyr::select(grna_id, grna_group),
  gene_id = gene_list
)

# create the directory
intermediate_data_dir <- paste0(.get_config_path("LOCAL_SPACRT_DATA_DIR"),
                                "full_data/intermediate_data")

if (!dir.exists(intermediate_data_dir)) {
  dir.create(intermediate_data_dir)
  cat("Directory created:", intermediate_data_dir, "\n")
} else {
  cat("Directory already exists:", intermediate_data_dir, "\n")
}

# save the sceptre objects
saveRDS(sceptre_object_left_grouped,
        sprintf("%s/sceptre_object_left_grouped.rds", intermediate_data_dir))
saveRDS(sceptre_object_right_grouped,
        sprintf("%s/sceptre_object_right_grouped.rds", intermediate_data_dir))
saveRDS(sceptre_object_left_individual,
        sprintf("%s/sceptre_object_left_individual.rds", intermediate_data_dir))

# save the positive/negative control lists
saveRDS(positive_control_list,
        sprintf("%s/positive_control_pairs.rds", intermediate_data_dir))
saveRDS(negative_control_list,
        sprintf("%s/negative_control_pairs.rds", intermediate_data_dir))

# save the subsetted ODM files as well as the covariate information
saveRDS(multimodal_odm_sub@modalities,
        sprintf("%s/subsetted_ODM.rds", intermediate_data_dir))

# remove the intercept term
covariate_matrix <- sceptre_object_left_grouped@covariate_matrix[, -1]
saveRDS(covariate_matrix,
        sprintf("%s/covariate_matrix.rds", intermediate_data_dir))
