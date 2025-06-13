# THis is a Rscript checking the l2 error of fastPhase estimation
library(dplyr)
library(tidyr)
library(ggplot2)
library(katlabutils)
library(patchwork)
library(cowplot)

# specify the directory the plot will be saved
plots_to_dir <- "manuscript/figures-and-tables//simulation/HMM-variable-selection"

# make the directory if is not there
if(!dir.exists(plots_to_dir)){
  dir.create(plots_to_dir)
}

# define the whole parameter grid
parameters_results <- expand.grid(
  beta_parameter = c("uniform", "skew"),
  lambda_parameter = c("lambda.1se", "lambda.min")
)

############################# Preprocess the data ##############################
# specify the results dir
results_dir <- paste0(
  .get_config_path("LOCAL_SPACRT_DATA_DIR"),
  "formal-simulation-results/variable_selection_HMM/variable_selection_HMM_results.rds"
)

# load the results
results_hmm <- readRDS(results_dir)

# load parameter
parameter_grid <- readRDS("simulation-code/sim-spec-formal/HMM_parameter_grid.rds")

# joint two dfs
full_results <- results_hmm$results |>
  dplyr::left_join(parameter_grid, by = "grid_id")

# extract information
alpha <- unique(full_results$alpha)
p <- unique(full_results$p)
signal_prop <- unique(full_results$signal_prop)

########################## Plot the results ####################################
for (grid in 1:nrow(parameters_results)) {

  # define the results of interest
  results_to_plot_beta <- parameters_results$beta_parameter[grid]
  results_to_plot_lambda <- parameters_results$lambda_parameter[grid]

  # append the null and alternative set; define the prior
  joined_results <- full_results |>
    dplyr::group_by(eta, gamma_0, method, beta_param, run_id) |>
    dplyr::mutate(
      null_set = list(set = dplyr::if_else(eta == 0, 1, round(signal_prop * p) + 1): p),
      signal_set = dplyr::if_else(eta == 0, list(NULL), list(dplyr::setdiff(1:p, null_set[[1]]))),
      beta_param = names(beta_param)
    ) |>
    dplyr::filter(beta_param == results_to_plot_beta)

  ############################## summarize dCRT and GCM ########################
  # unnest the data frame
  unnested_result <- joined_results |>
    dplyr::filter(method != "Knockoff")|>
    tidyr::unnest_wider(output) |>
    tidyr::unnest_wider(pvalue_list) |>
    dplyr::group_by(run_id, method, eta, gamma_0, null_set, signal_set, computation_time) |>
    dplyr::summarise(
      min_rejection = min_pvalue[[1]]$both_pvalue,
      ose_rejection = ose_pvalue[[1]]$both_pvalue
    ) |>
    dplyr::ungroup() |>
    tidyr::pivot_longer(cols = c("min_rejection", "ose_rejection"),
                        names_to = "lasso_model",
                        values_to = "p_value")

  # compute rejection size (both-sided)
  dCRT_GCM_both <- unnested_result |>
    dplyr::group_by(method, run_id, eta, gamma_0, lasso_model, null_set, signal_set, computation_time) |>
    dplyr::summarise(
      rejection_size_BH = sum((stats::p.adjust(p_value, method = "BH") <= alpha)),
      power_BH = sum((stats::p.adjust(p_value, method = "BH") <= alpha)[signal_set[[1]]]),
      FDP = sum((stats::p.adjust(p_value, method = "BH") <= alpha)[null_set[[1]]]) / max(rejection_size_BH, 1)
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(method, eta, gamma_0, lasso_model, null_set, signal_set) |>
    dplyr::summarise(
      mean_rejection_size = mean(rejection_size_BH, na.rm = TRUE),
      mean_power = mean(power_BH, na.rm = TRUE),
      FDR = mean(FDP, na.rm = TRUE),
      se_FDP = sd(FDP, na.rm = TRUE) / sqrt(dplyr::n()),
      se_power = sd(power_BH, na.rm = TRUE) / sqrt(dplyr::n()),
      mean_computation = mean(computation_time),
      se_computation = sd(computation_time) / sqrt(dplyr::n())
    ) |>
    dplyr::ungroup()

  #################### knockoff results ########################################
  unnest_knockoff <- joined_results |>
    dplyr::filter(method == "Knockoff")|>
    tidyr::unnest_wider(output) |>
    tidyr::unnest_wider(rejection_list) |>
    dplyr::group_by(run_id, method, eta, gamma_0, null_set, signal_set, computation_time) |>
    dplyr::summarise(
      min_rejection = min_rejection[[1]],
      ose_rejection = ose_rejection[[1]]
    ) |>
    dplyr::ungroup() |>
    tidyr::pivot_longer(cols = c("min_rejection", "ose_rejection"),
                        names_to = "lasso_model",
                        values_to = "rejection_set")

  knockoff_results <- unnest_knockoff |>
    dplyr::group_by(method, run_id, eta, gamma_0, lasso_model, null_set, signal_set, computation_time) |>
    dplyr::summarise(
      rejection_size = length(rejection_set[[1]]),
      power_BH = length(dplyr::intersect(signal_set[[1]], rejection_set[[1]])),
      FDP = dplyr::if_else(length(rejection_set[[1]]) > 0,
                           length(dplyr::intersect(null_set[[1]], rejection_set[[1]])) / length(rejection_set[[1]]),
                           0)
    ) |>
    dplyr::ungroup() |>
    dplyr::group_by(method, eta, gamma_0, lasso_model, null_set, signal_set) |>
    dplyr::summarise(
      mean_rejection_size = mean(rejection_size),
      mean_power = mean(power_BH),
      FDR = mean(FDP),
      se_FDP = sd(FDP) / sqrt(dplyr::n()),
      se_power = sd(power_BH) / sqrt(dplyr::n()),
      mean_computation = mean(computation_time),
      se_computation = sd(computation_time) / sqrt(dplyr::n())
    ) |>
    dplyr::ungroup()

  ######################### append two results and plotting ####################
  merged_results <- rbind(dCRT_GCM_both, knockoff_results) |>
    dplyr::mutate(
      sparsity = sprintf("gamma[0] == %d", gamma_0)
    ) |>
    dplyr::mutate(
      sparsity = factor(sparsity, levels = c("gamma[0] == -2", "gamma[0] == -3")),
      lasso_model = dplyr::if_else(lasso_model == "min_rejection",
                                   "lambda.min", "lambda.1se")
    ) |>
    dplyr::mutate(method = factor(method, levels = c("GCM", "Knockoff", "dCRT", "spaCRT"))) |>
    dplyr::mutate(method = factor(method,
                                  levels = c("GCM", "Knockoff", "dCRT", "spaCRT"),
                                  labels = c("GCM test", "Knockoff", "dCRT", "spaCRT")))

  # plot power for lambda.1se
  avg_power_1se <- merged_results |>
    dplyr::mutate(
      power_upper = mean_power + 1.96 * se_power,
      power_lower = mean_power - 1.96 * se_power
    ) |>
    dplyr::filter(lasso_model == results_to_plot_lambda) |>
    ggplot(aes(x = eta, y = mean_power, color = method)) +
    ggh4x::facet_grid2(. ~ sparsity,  labeller = labeller(
      sparsity = label_parsed  # Parse LaTeX-like expressions
    )) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = power_lower, ymax = power_upper), width = 0.1) +
    labs(x = "Signal strength", y = "Averaged power", color = "Method") +
    theme_bw() +
    theme(legend.position = c(0.63, 0.7),
          legend.background = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 14))

  # plot FDR for lambda.1se
  FDR_plot_1se <- merged_results |>
    dplyr::mutate(
      FDR_upper = FDR + 1.96 * se_FDP,
      FDR_lower = FDR - 1.96 * se_FDP
    ) |>
    dplyr::filter(lasso_model == results_to_plot_lambda) |>
    ggplot(aes(x = eta, y = FDR, color = method)) +
    ggh4x::facet_grid2(. ~ sparsity,  labeller = labeller(
      sparsity = label_parsed  # Parse LaTeX-like expressions
    )) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = alpha, color = "red", linetype = "dashed") +
    geom_errorbar(aes(ymin = FDR_lower, ymax = FDR_upper), width = 0.1) +
    labs(x = "Signal strength", y = "FDR") +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 14))

  # plot the time comparison
  time_comparison <- merged_results |>
    dplyr::filter(lasso_model == results_to_plot_lambda) |>
    dplyr::mutate(metric = "Computation time") |>
    ggplot(aes(x = method, y = mean_computation, fill = method)) +
    geom_boxplot() +
    scale_y_log10() +
    facet_grid(. ~ metric) +
    labs(x = "Method", y = "Second per replication") +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 14))

  # Remove legend from FDR_plot
  FDR_plot_clean_1se <- FDR_plot_1se + theme(legend.position = "none")

  # Remove legend from avg_power (to prevent duplication)
  avg_power_clean_1se <- avg_power_1se

  # Create the left panel containing p1 and p2 stacked vertically
  joint_plot <- ((FDR_plot_clean_1se / avg_power_clean_1se) | time_comparison) +
    plot_layout(ncol = 2, widths = c(1.5, 1)) +
    plot_annotation(tag_levels = "a") +
    theme(legend.position = "none",
          legend.title = element_blank(),
          plot.margin = unit(c(0.05, 0.05, 0, 0.05), "cm"))

  # save the plots
  ggsave(sprintf("%s/HMM_simulation_%s_%s.pdf", plots_to_dir, results_to_plot_beta, results_to_plot_lambda),
         joint_plot, height = 6, width = 10)
}
