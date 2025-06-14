# load the sceptre result to obtain the n_nonzero_trt; this code is for Table 3
sceptre_result <- readRDS(sprintf("%s/full_data/sceptre/sceptre_result_left_update_discovery.rds",
                                  .get_config_path("LOCAL_SPACRT_DATA_DIR")))

# compute the sparsity level
sparsity_summary <- round(summary(sceptre_result$n_nonzero_trt), digits = 0)

# print the sparsity level
print(sparsity_summary)

# print the sparsity rate
quantile_list <- c(0.01, 0.99)
round(quantile(sceptre_result |>
                 dplyr::mutate(ratio = n_nonzero_trt / (n_nonzero_trt + n_nonzero_cntrl)) |>
                 dplyr::pull(),
               quantile_list),
      digits = 4)
