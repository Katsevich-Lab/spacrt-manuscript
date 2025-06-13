# THis is a Rscript checking the l2 error of fastPhase estimation
library(dplyr)
library(tidyr)
library(ggplot2)
library(katlabutils)
library(patchwork)
library(cowplot)
library(scales)

############################# Preprocess the data ##############################
# specify the results dir
results_dir <- paste0(
  .get_config_path("LOCAL_SPACRT_DATA_DIR"),
  "formal-simulation-results/random_forest/random_forest_results.rds"
)

# load the results
results_rf <- readRDS(results_dir)

# load parameter
parameter_grid <- readRDS('simulation-code/sim-spec-formal/random_forest_parameter_grid.rds')

# specify the directory to be saved
plots_dir <- "manuscript/figures-and-tables//simulation/RF-regression"
if(!dir.exists(plots_dir)){
  dir.create(plots_dir)
}

# set the significance level
level <- 0.005

# preprocess the results
joined_results <- results_rf$results %>%
  dplyr::left_join(parameter_grid, by = "grid_id") %>%
  dplyr::mutate(method = factor(method, levels = c('GCM', 'dCRT', 'spaCRT'))) %>%
  dplyr::arrange(grid_id,
                 run_id,
                 method) %>%
  mutate(p_left = purrr::map_dbl(output, ~ .x$pvalue_list$p.left)) %>%
  mutate(p_right = purrr::map_dbl(output, ~ .x$pvalue_list$p.right)) %>%
  mutate(comp_time = purrr::map_dbl(output, ~ .x$computation_time)) %>%
  select(-c(ground_truth, p, B, output))

# compute the rates of spaCRT using GCM default
results_rf$results |>
  dplyr::filter(method == "spaCRT") |>
  mutate(spa.success = purrr::map_dbl(output, ~ .x$pvalue_list$spa.success)) |>
  dplyr::group_by(grid_id) |>
  dplyr::summarise(GCM_default = 1 - mean(spa.success)) |>
  dplyr::ungroup() |>
  dplyr::summarise(max_GCM_default = max(GCM_default)) |>
  dplyr::pull()

# obtain the p-value dataframe
pvals_result_rf <- joined_results %>%
  mutate(p_left_marker = p_left <= level) %>%
  mutate(p_right_marker = p_right <= level) %>%
  mutate(gamma_val = dplyr::case_when(gamma_0 == -2.5 ~ "gamma[0]==-2.5",
                                      gamma_0 == -2 ~ "gamma[0]==-2")) |>
  mutate(gamma_val = factor(gamma_val,
                            levels = c("gamma[0]==-2.5",
                                       "gamma[0]==-2")))

###########################################################################################
## Power plot #############################################################################
###########################################################################################

# extract the default colors
default_colors <- hue_pal()(4)

# plot power
power_plots_rf <- pvals_result_rf %>%
  group_by(method, eta, gamma_val) %>%
  summarise(
    mean_left = mean(p_left_marker),
    mean_right = mean(p_right_marker),
    sd_left = sd(p_left_marker, na.rm = TRUE),
    sd_right = sd(p_right_marker, na.rm = TRUE),
    n = n(),
    se_left = sd_left / sqrt(n),
    se_right = sd_right / sqrt(n),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    toe_upper = mean_left + 1.96 * se_left,
    toe_lower = mean_left - 1.96 * se_left
  ) %>%
  ggplot(aes(x = eta, y = mean_left, colour = method)) +
  facet_grid(. ~ gamma_val, scales = "free_x", labeller = label_parsed) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = toe_lower, ymax = toe_upper), width = 0.02) +
  geom_hline(aes(yintercept = level), linetype = "dashed", color = "red") +
  scale_x_continuous() +
  theme_bw() +
  scale_color_manual(values = c("GCM" = default_colors[1],
                                "dCRT" = default_colors[3],
                                "spaCRT" = default_colors[4])) +
  labs(x = latex2exp::TeX("\\eta"), y = "Power") +
  theme(legend.position = c(0.65, 0.787),
        legend.key.size = unit(3, "mm"),
        legend.title = element_blank(),
        legend.box = "vertical",
        legend.box.background = element_rect(
          color = "black",    # Border color
          linewidth = 0.2,    # Border thickness
          linetype = "solid",
          fill = NA           # Transparent background
        ),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.box.spacing = unit(1, "mm"),  # Space between legend and plot
        plot.margin = unit(c(5, 5, 10, 5), "mm"),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))


###########################################################################################
## QQ plots ###############################################################################
###########################################################################################
# plot QQ-plot
qq_plots_rf <- pvals_result_rf %>%
  filter(eta == 0) %>%
  ggplot(mapping = aes(y = p_left, color = method)) +
  facet_grid(. ~ gamma_val, labeller = label_parsed) +
  stat_qq_points(ymin = 1e-8, size = 0.8) +
  stat_qq_band() +
  scale_x_continuous(trans = revlog_trans(10),
                     breaks = c(1e-1, 1e-2, 1e-3, 1e-4),
                     labels = c(expression(10^{-1}), expression(10^{-2}), expression(10^{-3}), expression(10^{-4}))) +
  scale_y_continuous(trans = revlog_trans(10),
                     breaks = c(1e-1, 1e-3, 1e-5, 1e-7),
                     labels = c(expression(10^{-1}), expression(10^{-3}), expression(10^{-5}), expression(10^{-7}))) +
  scale_color_manual(values = c("GCM" = default_colors[1],
                                "dCRT" = default_colors[3],
                                "spaCRT" = default_colors[4])) +
  geom_abline(col = "black") +
  labs(y = "Observed p-value",
       x = "Expected p-value") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))


###########################################################################################
## Time comparison ########################################################################
###########################################################################################
# plot time comparison
time_comparison_rf <- joined_results %>%
  dplyr::group_by(method, grid_id) %>%
  dplyr::summarise(mean_comp_time = mean(comp_time)) |>
  dplyr::ungroup() |>
  dplyr::mutate(plot_title = "Computation~time") |>
  ggplot(aes(x = method, y = mean_comp_time, fill = method)) +
  facet_grid(. ~ plot_title, labeller = label_parsed) +
  geom_boxplot() +
  scale_y_log10() +
  labs(x = "Method", y = "Second per replication") +
  theme_bw() +
  scale_fill_manual(values = c("GCM" = default_colors[1],
                               "dCRT" = default_colors[3],
                               "spaCRT" = default_colors[4])) +
  theme(legend.position = "none",
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))


# Combine three plots
(joint_plot <- ((qq_plots_rf / power_plots_rf) | time_comparison_rf) +
  plot_layout(ncol = 2, widths = c(1.8, 1)) +
  theme(#legend.position = "none",
        legend.title = element_blank(),
        plot.margin = unit(c(0.1, 0.15, 0.1, 0.1), "cm"))
)

# Save automatically
ggsave(sprintf("%s/simulation-summary.pdf", plots_dir),
       plot = joint_plot,
       height = 6, width = 8)

