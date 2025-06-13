library(tidyr)
library(dplyr)
library(cowplot)
library(ggplot2)
library(katlabutils)

disp.char <- '5e-2'
path_rds <- paste0(.get_config_path("LOCAL_SPACRT_DATA_DIR"),
                   'private/results/simulation-results/formal-results/',
                   'power_grid_bin_NB_normal_B_50000_disp_',
                   disp.char,'_L_R_formal/',
                   'power_grid_bin_NB_normal_B_50000_disp_',
                   disp.char,'_L_R_formal_results.rds')
results <- readRDS(path_rds)

# Type-I error and power -------------------------------------------------------
metrics <- results$metrics

metrics_filtered <- metrics |>
  filter(grepl("level", metric)) |>
  separate_wider_regex(
    cols = metric,
    patterns = c("level\\.", "level" = "[0-9]+", "\\.", "side" = "[A-Z]")
  ) |>
  filter(level == "2", beta_0 == -5) |>
  mutate(side = factor(side,
                       levels = c("L", "R"),
                       labels = c("left-sided", "right-sided"))) |>
  mutate(method = factor(method,
                         levels = c("GCM", "score.test", "dCRT", "spaCRT"),
                         labels = c("GCM test", "Score test", "dCRT", "spaCRT"))) |>
  select(-c(level, beta_0, se, n))

# Type-I error
p1 <- metrics_filtered |>
  filter(rho == 0) |>
  mutate(mean = ifelse(mean == 0, min(mean[mean > 0]), mean), .by = side) |>
  ggplot(aes(x = gamma_0, y = mean, color = method)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 1e-2, linetype = "dashed") +
  scale_y_log10() +
  facet_wrap(~side) +
  labs(x = expression(gamma[0]),
       y = "Type-I error") +
  theme_bw() +
  theme(legend.position = "none")

# Power
p2 <- metrics_filtered |>
  filter(gamma_0 == -3) |>
  filter((side == "left-sided" & rho <= 0) | (side == "right-sided" & rho >= 0)) |>
  ggplot(aes(x = rho, y = mean, color = method)) +
  geom_point() +
  geom_line() +
  facet_wrap(~side, scales = "free_x") +
  labs(x = expression(rho),
       y = "Power") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.box.spacing = unit(5, "pt"))

# QQ plots ---------------------------------------------------------------------

varying_params <- list(
  n = 5000,
  gamma_0 = seq.int(from = -6, to = -2),
  beta_0 = seq.int(from = -6, to = -2))

baseline_params <- list(n = 5e3, gamma_0 = -5, beta_0 = -5)

grid_ffd <- simulatr::create_param_grid_fractional_factorial(varying_params,
                                                             baseline_params)

rho <- c(seq.int(from = -4, to = 0), seq.int(from = 1, to = 4)/2)

parameter_grid <- grid_ffd |>
  dplyr::select(1:3) |>
  slice(rep(1:n(), each = length(rho))) |>
  cbind(rho = rep(rho, times = nrow(grid_ffd))) |>
  dplyr::mutate(grid_id = seq(1, length(rho), 1))

grid_id <- parameter_grid |>
  filter(beta_0 == -5, gamma_0 == -3, rho == 0) |>
  pull(grid_id)

results_gridrow <- results$results |>
  filter(grid_id == !!grid_id)

p_values <- results_gridrow |>
  rowwise() |>
  mutate(p_left = output$p.left, p_right = output$p.right) |>
  select(-c(output, grid_id)) |>
  ungroup() |>
  pivot_longer(cols = c("p_left", "p_right"),
               names_prefix = "p_",
               names_to = "side",
               values_to = "p_value")

p3 <- p_values |>
  mutate(side = factor(side,
                       levels = c("left", "right"),
                       labels = c("left-sided", "right-sided"))) |>
  mutate(method = factor(method,
                         levels = c("GCM", "score.test", "dCRT", "spaCRT"),
                         labels = c("GCM test", "Score test", "dCRT", "spaCRT"))) |>
  ggplot(mapping = aes(y = p_value, color = method)) +
  stat_qq_points(ymin = 1e-8, size = 0.8, max_pts_to_plot = 1000) +
  stat_qq_band() +
  geom_abline() +
  facet_wrap(~ side, ncol = 2) +
  labs(y = "Observed p-value",
       x = "Expected p-value") +
  scale_x_continuous(trans = revlog_trans(10),
                     breaks = c(1e-1, 1e-3, 1e-5),
                     labels = c(expression(10^{-1}), expression(10^{-3}), expression(10^{-5}))) +
  scale_y_continuous(trans = revlog_trans(10),
                     breaks = c(1e-1, 1e-3, 1e-5, 1e-7),
                     labels = c(expression(10^{-1}), expression(10^{-3}), expression(10^{-5}), expression(10^{-7}))) +
  theme_bw() +
  theme(legend.position = "none")

## Combine plots ----------------------------------------------------------------

# This code is for Figure 2
p <- plot_grid(p3, p1, p2, ncol = 1, rel_heights = c(1, 1, 1.2), align = "v", labels = "auto")

# save the plot
ggsave(filename = "manuscript/figures-and-tables/simulation-summary.pdf",
       plot = p,
       height = 7,
       width = 4.25)

## Computing time ---------------------------------------------------------------

# This code is for Figure 3
p4 <- metrics |>
  filter(metric == "hrs_per_rep", beta_0 == -5) |>
  mutate(method = factor(method,
                         levels = c("GCM", "score.test", "dCRT", "spaCRT"),
                         labels = c("GCM test", "Score test", "dCRT", "spaCRT"))) |>
  ggplot(aes(x = method, fill = method, y = mean * 60 * 60)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(x = "Method", y = "Second per replication") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_blank())

ggsave(filename = "manuscript/figures-and-tables/simulation-computing-times.pdf",
       plot = p4,
       height = 4,
       width = 3.5)
