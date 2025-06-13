# This is a Rscript handling the power analysis for Gasperini dataset
# get the sceptre object used in Tim's code
library(spacrt)
library(sceptre)
library(ondisc)

# Read the command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign the arguments to variables in R
min_gRNA_count <- as.numeric(args[1])

# set the intermediate data directory
intermediate_data_dir <- paste0(.get_config_path("LOCAL_SPACRT_DATA_DIR"),
                                "full_data/intermediate_data")

# set the result path
result_path <- paste0(
  .get_config_path("LOCAL_SPACRT_DATA_DIR"),
  "full_data/power_three_methods"
)

# create it if the directory does not exist
if (!dir.exists(result_path)) {
  dir.create(result_path)
  cat("Directory created:", result_path, "\n")
} else {
  cat("Directory already exists:", result_path, "\n")
}

# load the covariate matrix, positive control list and subsetted ODM
covariate_data <- readRDS(sprintf("%s/covariate_matrix.rds",
                                  intermediate_data_dir))
positive_control_list <- readRDS(sprintf("%s/positive_control_pairs.rds",
                                         intermediate_data_dir))
subsetted_ODM <- readRDS(sprintf("%s/subsetted_ODM.rds",
                                 intermediate_data_dir))

# create grouping list
grouping <- c("grouped", "individual")

# create methods list
method_list <- c("GCM", "spaCRT", "scoreglmnb")

# store the list of gene, grna group and grna id
gene_list <- positive_control_list$grna_target |> unique()
grna_group_list <- positive_control_list$grna_group |> unique()
grna_id_list <- positive_control_list$grna_id |> unique()

# loop over grouping strategy
for (grouping_type in grouping) {
  # compute how many pairs needed to be computed
  num_hyp <- dplyr::if_else(grouping_type == "grouped",
                            length(grna_group_list),
                            length(grna_id_list))

  # extract grna_gene pairs
  grna_gene <- apply(as.matrix(positive_control_list |>
                                 dplyr::select(-grna_group)), 1, function(x){
                                   sprintf("%s_%s", x[1], x[2])
                                 })

  # create empty set
  pvalue <- matrix(data = NA,
                   nrow = num_hyp,
                   ncol = length(method_list))

  if(grouping_type == "grouped"){
    dimnames(pvalue) <- list(gene = gene_list, methods = method_list)
  }else{
    dimnames(pvalue) <- list(grna_gene = grna_gene, methods = method_list)
  }

  # specify the time matrix
  timing <- pvalue

  # specify GCM default list
  default_GCM <- numeric(num_hyp)

  # specify the effective sample size list
  ess <- numeric(num_hyp)

  # specify the dispersion list
  dispersion <- numeric(num_hyp)

  # loop over number of hypothesis
  for (hyp in 1:num_hyp) {
    if(grouping_type == "grouped"){
      grna_ids <- positive_control_list |>
        dplyr::filter(grna_group == grna_group_list[hyp]) |>
        dplyr::select(grna_id) |>
        dplyr::pull()
    }else{
      grna_ids <- grna_id_list[hyp]
    }

    # obtain the grna presences
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

    # obtain the gene expression
    current_gene <- positive_control_list |>
      dplyr::filter(grna_id %in% grna_ids) |>
      dplyr::select(grna_target) |>
      dplyr::pull() |>
      unique()
    gene_expression <- subsetted_ODM$gene[[current_gene,]] |>
      as.vector()

    # create the data list
    data <- list(
      Y = gene_expression,
      X = as.numeric(grna_presences),
      Z = covariate_data
    )

    ## Provide the auxiliary information for negative.binomial model
    aux_info_Y_on_Z <- spacrt:::nb_precomp(V = data$Y, Z = data$Z)

    # run GCM method
    # store the computation time
    timing[hyp, "GCM"] <- system.time(
      test_result_GCM <- spacrtutils::GCM_internal(data = data,
                                                   X_on_Z_fam = "binomial",
                                                   Y_on_Z_fam = "negative.binomial")
    )[["elapsed"]]

    # compute and store the p-value
    pvalue[hyp, "GCM"] <- test_result_GCM$p.left


    # run spaCRT method
    # store the computation time
    timing[hyp, "spaCRT"] <- system.time(
      test_result_spaCRT <- spacrtutils::spaCRT_internal(data = data,
                                                         X_on_Z_fam = "binomial",
                                                         Y_on_Z_fam = "negative.binomial")
    )[["elapsed"]]

    # compute and store the p-value
    pvalue[hyp, "spaCRT"] <- test_result_spaCRT$p.left

    # store the default choice of GCM into dataframe
    default_GCM[hyp] <- 1 - test_result_spaCRT$spa.success

    # run score test
    # store the computation time
    timing[hyp, "scoreglmnb"] <- system.time(
      test_result_score <- spacrtutils::score.test(data = data,
                                                   X_on_Z_fam = "binomial",
                                                   Y_on_Z_fam = "negative.binomial")
    )[["elapsed"]]

    # compute and store the p-value
    pvalue[hyp, "scoreglmnb"] <- test_result_score$p.left

    # store the dispersion
    dispersion[hyp] <- test_result_score$NB.disp.param

    # print the result
    print(pvalue[hyp, ])
    print(default_GCM[hyp])
    print(timing[hyp, ])

    # compute the effective sample size
    ess[hyp] <- sum(grna_presences*gene_expression)
    print(ess[hyp])
  }

  # save the pvalue files
  saveRDS(pvalue,
          paste0(result_path, sprintf("/pvalue_power_%s_left.rds", grouping_type))
  )

  # save the timing files
  saveRDS(timing,
          paste0(result_path, sprintf("/timing_power_%s_left.rds.rds", grouping_type))
  )

  # save the gcm default
  saveRDS(default_GCM,
          paste0(result_path, sprintf("/default_GCM_power_%s.rds", grouping_type))
  )

  # save the effective sample size default
  saveRDS(ess,
          paste0(result_path, sprintf("/%s_effective_sample_size.rds",
                                      grouping_type))
  )

  # save the dispersion list
  saveRDS(dispersion,
          paste0(result_path, sprintf("/%s_dispersion.rds",
                                      grouping_type))
  )
}

