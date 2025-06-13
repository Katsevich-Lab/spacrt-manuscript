# This is a Rscript postprocessing the simulation results
library(sceptre)
library(reshape2)
library(dplyr)
library(katlabutils)
library(tidyr)

# Read the command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign the arguments to variables in R
minimum_cutoff <- as.numeric(args[1])

# set directory
data_dir <- paste0(
  .get_config_path("LOCAL_SPACRT_DATA_DIR"),
  "full_data"
)

intermediate_data_dir <- paste0(.get_config_path("LOCAL_SPACRT_DATA_DIR"),
                                "full_data/intermediate_data")

# check the result directory exists or not; if not, create one.
result_dir <- paste0(
  .get_config_path("LOCAL_SPACRT_DATA_DIR"),
  "full_data/summary"
)
if (!dir.exists(result_dir)) {
  dir.create(result_dir)
  cat("Directory created:", result_dir, "\n")
} else {
  cat("Directory already exists:", result_dir, "\n")
}

# read the estimated dispersion parameter
dispersion_list <- readRDS(sprintf("%s/dispersion_list_glm_nb.rds", data_dir))

# data postprocessing for null simulation
################################################################################

# read the negative control list
negative_control_list <- readRDS(sprintf("%s/negative_control_pairs.rds",
                                         intermediate_data_dir))

# get lists of NT gRNA groups
grna_groups <- negative_control_list$grna$grna_group |> unique()

# obtain the gene lists
gene_list <- negative_control_list$gene_id

# compute the number of hypothesis
num_genes <- length(gene_list)
num_grna <- length(grna_groups)
methods <- c("GCM", "spaCRT", "scoreglmnb")
pvalue_matrix_left <- array(data = NA,
                            dim = c(length(methods), num_grna, num_genes))

dimnames(pvalue_matrix_left) <- list(method = methods,
                                     gRNA = grna_groups,
                                     gene = gene_list)

pvalue_matrix_right <- pvalue_matrix_left

# loop over grna groups and gene list
for (name_RNA in grna_groups) {
  GCM_result <- readRDS(sprintf("%s/GCM/%s/output.rds",
                                data_dir,
                                name_RNA))
  spaCRT_result <- readRDS(sprintf("%s/spaCRT/%s/output.rds",
                                   data_dir,
                                   name_RNA))
  score_result <- readRDS(sprintf("%s/scoreglmnb/%s/output.rds",
                                  data_dir,
                                  name_RNA))

  for (name_gene in gene_list) {
    # store the left-side p-value
    pvalue_matrix_left[methods, name_RNA, name_gene] <- c(unname(GCM_result[name_gene, "left"]),
                                                          unname(spaCRT_result[name_gene, "left"]),
                                                          unname(score_result[name_gene, "left"]))

    # store the right-side p-value
    pvalue_matrix_right[methods, name_RNA, name_gene] <- c(unname(GCM_result[name_gene, "right"]),
                                                           unname(spaCRT_result[name_gene, "right"]),
                                                           unname(score_result[name_gene, "right"]))
  }
}

# load the sceptre result
sceptre_result_left <- readRDS(sprintf("%s/sceptre/sceptre_result_left_update_discovery.rds",
                                       data_dir))
sceptre_result_right <- readRDS(sprintf("%s/sceptre/sceptre_result_right_update_discovery.rds",
                                        data_dir))

# melt the array
p_value_left <- reshape2::melt(pvalue_matrix_left)
p_value_right <- reshape2::melt(pvalue_matrix_right)

# cbind the left and right pvalues and impute the zero value in the p-value matrix using machine precision
p_value <- p_value_left |>
  left_join(p_value_right, by = c("method", "gRNA", "gene")) |>
  dplyr::rename(p_value_left = value.x,
                p_value_right = value.y) |>
  mutate(across(
    .cols = p_value_left:p_value_right,
    .fns = function(x) (pmax(x, .Machine$double.eps))
  )) |>
  pivot_wider(names_from = method,
              values_from = c(p_value_left, p_value_right))

# load the dispersion parameter and build the dispersion group
dispersion_mat <- cbind(round(dispersion_list, 2), gene_list)
colnames(dispersion_mat) <- c("dispersion", "response_id")

# filter sceptre result according to n_nonzero_cntrl_thresh
filter_sceptre <- sceptre_result_left |>
  dplyr::left_join(sceptre_result_right |>
                     dplyr::select("response_id", "grna_target", "p_value"),
                   by= c("response_id", "grna_target"),
                   suffix = c("_left", "_right")) |>
  left_join(data.frame(dispersion_mat), by = "response_id") |>
  filter(n_nonzero_trt >= minimum_cutoff) |>
  rename(p_value_left_sceptre = p_value_left,
         p_value_right_sceptre = p_value_right,
         gRNA = grna_target,
         gene = response_id)

# combined result
combined_result <- filter_sceptre |>
  left_join(p_value, by = c("gRNA", "gene")) |>
  dplyr::select(-c("n_nonzero_cntrl")) |>
  melt(value.name = "p-value",
       variable.name = "method",
       id.vars = c("gRNA", "gene","dispersion", "n_nonzero_trt")) |>
  group_by(method) |>
  mutate(side = dplyr::if_else(length(grep("left", method)) == 0,
                               "right-sided", "left-sided")) |>
  ungroup() |>
  dplyr::mutate(method = as.factor(sapply(strsplit(as.character(method), "_"),
                                          function(x) tail(x, 1)))) |>
  dplyr::mutate(method = dplyr::if_else(method == "scoreglmnb", "score", method)) |>
  dplyr::mutate(method = factor(method, levels = c("GCM", "score", "spaCRT", "sceptre")))

# save the result
saveRDS(combined_result,
        sprintf("%s/combined_result_null.rds", result_dir))

# data postprocessing for power simulation
################################################################################

# grouping_type list
grouping_list <- c("grouped", "individual")

# loop over grouping list
for (grouping_type in grouping_list) {

  # load the results
  pvalue_power_left <- readRDS(sprintf("%s/power_three_methods/pvalue_power_%s_left.rds",
                                       data_dir, grouping_type))
  sceptre_result_left_update_power <- readRDS(sprintf("%s/sceptre/sceptre_result_left_%s_update_power.rds",
                                                      data_dir, grouping_type))

  if(grouping_type == "grouped"){

    # combine the sceptre and the other three methods
    combined_result <- pvalue_power_left |> as.data.frame() |>
      tibble::rownames_to_column(var = "grna_target") |>
      as_tibble() |>
      dplyr::left_join(sceptre_result_left_update_power |>
                         dplyr::select(grna_target, n_nonzero_trt, sceptre = p_value),
                       by = "grna_target") |>
      tidyr::pivot_longer(cols = c("GCM", "spaCRT", "sceptre", "scoreglmnb"),
                          values_to = "p_value",
                          names_to = "method") |>
      dplyr::mutate(method = dplyr::if_else(method == "scoreglmnb", "score", method))
  }else{

    # combine the sceptre and the other three methods
    combined_result <- pvalue_power_left |> as.data.frame() |>
      tibble::rownames_to_column(var = "grna_gene") |>
      as_tibble() |>
      dplyr::left_join(sceptre_result_left_update_power |>
                         dplyr::mutate(grna_gene = sprintf("%s_%s", grna_id, grna_target)) |>
                         dplyr::select(grna_gene, n_nonzero_trt, sceptre = p_value),
                       by = "grna_gene") |>
      tidyr::pivot_longer(cols = c("GCM", "spaCRT", "sceptre", "scoreglmnb"),
                          values_to = "p_value",
                          names_to = "method") |>
      dplyr::mutate(method = dplyr::if_else(method == "scoreglmnb", "score", method))
  }

  # fill the 0 as machine precision value
  combined_result[combined_result == 0] <- .Machine$double.eps

  # rearrange the order for methods
  combined_result$method <- factor(combined_result$method,
                                   levels = c("GCM", "score",
                                              "spaCRT", "sceptre"))

  # save combined_result
  saveRDS(combined_result,
          sprintf("%s/combined_result_power_%s.rds", result_dir, grouping_type))
}

