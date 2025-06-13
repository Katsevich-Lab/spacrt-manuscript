suppressPackageStartupMessages(library(simulatr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(gtable))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(latex2exp))
suppressPackageStartupMessages(library(katlabutils))
suppressPackageStartupMessages(library(cowplot))

# specify the dispersion parameter
disp <- 10

if(disp == 1 | disp == 10) disp.char <- as.character(disp)
if(disp == 0.05) disp.char <- '5e-2'

path_rds <- paste0(.get_config_path("LOCAL_SPACRT_DATA_DIR"),
                   'private/results/simulation-results/formal-results/',
                   'power_grid_bin_NB_normal_B_50000_disp_',
                   disp.char,'_L_R_formal/',
                   'power_grid_bin_NB_normal_B_50000_disp_',
                   disp.char,'_L_R_formal_results.rds')

results_rds <- readRDS(path_rds)

print('Results extracted!')

B <- 50000

default_gamma <- -5
default_beta <- -5

plots_dir <- "manuscript/figures-and-tables/simulation"

# compute the default GCM proportion
max(results_rds$results |>
      filter(method == "spaCRT") |>
      group_by(grid_id) |>
      summarise(default = mean(sapply(output, function(a) a$gcm.default))) |>
      ungroup() |> dplyr::select(default) |> pull())

# 0.00596

###########################################################################################
## Type-I Error ###########################################################################
###########################################################################################

# The following code is for left-sided test
#####################

results_rds_metrics <- results_rds$metrics %>% filter(metric %in% c('level.1.L',
                                                                    'level.2.L',
                                                                    'level.3.L'))

# Type-I error plot
Type_I_result <- results_rds_metrics |>
  dplyr::arrange(method) %>%
  dplyr::filter(rho == 0) |>
  dplyr::mutate(metric = factor(metric)) |>
  dplyr::mutate(metric = forcats::fct_recode(metric,
                                             "alpha = 0.05" = "level.1.L",
                                             "alpha = 0.01" = "level.2.L",
                                             "alpha = 0.005" = "level.3.L")) |>
  dplyr::mutate(method = factor(method)) |>
  dplyr::mutate(method = forcats::fct_recode(method,
                                             "score" = "score.test")) |>
  dplyr::mutate(method = factor(method, levels = c("GCM",
                                                   "score",
                                                   "dCRT",
                                                   "spaCRT")))

######################
# plot the result for Type-I error; this code is for Figure 10
Type_I_result <- Type_I_result |>
  dplyr::mutate(
    level = case_when(
      metric == "alpha = 0.05" ~ 0.05,
      metric == "alpha = 0.01" ~ 0.01,
      metric == "alpha = 0.005" ~ 0.005,
      TRUE ~ NA_integer_  # Integer output for default case (NA as integer)
    )
  ) |>
  dplyr::mutate(mean = ifelse(mean == 0, 1/B, mean),
                sd = ifelse(se == 0, sqrt((1 - 1/B)) / B, se)) %>%
  dplyr::mutate(metric.num = as.numeric(unlist(strsplit(as.character(metric),
                                                        split='= ', fixed=TRUE))[2]))

gamma_plot <- Type_I_result |>
  dplyr::filter(gamma_0 == -5) |>
  dplyr::mutate(gamma_0 = sprintf("r = %d", gamma_0)) |>
  ggplot(aes(x = beta_0,
             y = mean,
             color = method)) +
  ggh4x::facet_grid2(gamma_0 ~ metric, scales = "free_y") +
  geom_point(size = 0.8) +
  geom_line(linewidth = 0.5) +
  geom_hline(aes(yintercept = level), linetype = "dashed", color = "red") +
  scale_x_continuous() +
  scale_y_log10() +
  theme_bw() +
  labs(x = "Varying b") +
  theme(strip.text.x = element_blank(),
        strip.text.y = element_text(margin = margin(0,0.08,0,0.08, "cm")),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none"
  )


beta_plot <- Type_I_result |>
  dplyr::filter(beta_0 == -5) |>
  dplyr::mutate(beta_0 = sprintf("b = %d", beta_0)) |>
  ggplot(aes(x = gamma_0,
             y = mean,
             color = method)) +
  scale_x_continuous() +
  scale_y_log10() +
  ggh4x::facet_grid2(beta_0 ~ metric, scales = "free_y") +
  geom_point(size = 0.8) +
  geom_line(linewidth = 0.5) +
  geom_hline(aes(yintercept = level), linetype = "dashed", color = "red") +
  theme_bw() +
  labs(x = "Varying r") +
  theme(strip.text.x = element_blank(),
        strip.text.y = element_text(margin = margin(0,0.08,0,0.08, "cm")),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none"
  )


plot_aux <- Type_I_result |>
  dplyr::filter(beta_0 == -5) |>
  dplyr::mutate(beta_0 = sprintf("b = %d", beta_0)) |>
  ggplot(aes(x = gamma_0,
             y = mean,
             color = method)) +
  scale_x_continuous() +
  scale_y_log10() +
  ggh4x::facet_grid2(metric ~ beta_0, scales = "free_y") +
  geom_point(size = 0.5) +
  geom_line() +
  geom_hline(aes(yintercept = level), linetype = "dashed", color = "red") +
  theme_bw() +
  theme(strip.text.x = element_text(margin = margin(0.1,0,0.1,0,"cm")),
        strip.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="bottom"
  )

# get legend from plot_aux
legend <- cowplot::get_plot_component(plot_aux, 'guide-box-bottom', return_all = TRUE)

# combine the plot
plot <- cowplot::plot_grid(gamma_plot, beta_plot,
                           nrow = 2, align = "v")

#create common x and y labels
y.grob <- textGrob("Type I error",
                   gp=gpar(col="black"), rot=90)


# add to plot
final_plot <- grid.arrange(arrangeGrob(plot, left = y.grob, nrow=1),
                           cowplot::plot_grid(legend, nrow = 1),
                           top = textGrob(TeX(paste0("Left-sided test with \\theta = ",
                                                     disp)),
                                          gp=grid::gpar(fontface=2)),
                           nrow=2, heights=c(9, 1))

ggsave(paste0(sprintf('%s/Type-I-error/plot-bin-NB-normal-B-50000-n-5000-5e3-n5-n5-disp-',
                      plots_dir), disp.char,'-Type-I-error-LEFT.pdf'),
       plot = final_plot,
       width = 6, height = 4.5)

rm(Type_I_result, gamma_plot, beta_plot, plot_aux, plot, final_plot)

################################################################################
# The following code is for right-sided test; this code is for Figure 11
################################################################################

results_rds_metrics <- results_rds$metrics %>% filter(metric %in% c('level.1.R',
                                                                    'level.2.R',
                                                                    'level.3.R'))

####### Type-I error plot
Type_I_result <- results_rds_metrics |>
  dplyr::arrange(method) %>%
  dplyr::filter(rho == 0) |>
  dplyr::mutate(metric = factor(metric)) |>
  dplyr::mutate(metric = forcats::fct_recode(metric,
                                             "alpha = 0.05" = "level.1.R",
                                             "alpha = 0.01" = "level.2.R",
                                             "alpha = 0.005" = "level.3.R")) |>
  dplyr::mutate(method = factor(method)) |>
  dplyr::mutate(method = forcats::fct_recode(method,
                                             "score" = "score.test")) |>
  dplyr::mutate(method = factor(method, levels = c("GCM",
                                                   "score",
                                                   "dCRT",
                                                   "spaCRT")))

########

# Type-I error plot
Type_I_result <- results_rds_metrics |>
  dplyr::arrange(method) %>%
  dplyr::filter(rho == 0) |>
  dplyr::mutate(metric = factor(metric)) |>
  dplyr::mutate(metric = forcats::fct_recode(metric,
                                             "alpha = 0.05" = "level.1.R",
                                             "alpha = 0.01" = "level.2.R",
                                             "alpha = 0.005" = "level.3.R")) |>
  dplyr::mutate(method = factor(method)) |>
  dplyr::mutate(method = forcats::fct_recode(method,
                                             "score" = "score.test")) |>
  dplyr::mutate(method = factor(method, levels = c("GCM",
                                                   "score",
                                                   "dCRT",
                                                   "spaCRT")))

######################
# plot the result for Type-I error
Type_I_result <- Type_I_result |>
  dplyr::mutate(
    level = case_when(
      metric == "alpha = 0.05" ~ 0.05,
      metric == "alpha = 0.01" ~ 0.01,
      metric == "alpha = 0.005" ~ 0.005,
      TRUE ~ NA_integer_  # Integer output for default case (NA as integer)
    )
  ) |>
  dplyr::mutate(mean = ifelse(mean == 0, 1/B, mean),
                sd = ifelse(se == 0, sqrt((1 - 1/B)) / B, se)) %>%
  dplyr::mutate(metric.num = as.numeric(unlist(strsplit(as.character(metric),
                                                        split='= ', fixed=TRUE))[2]))

gamma_plot <- Type_I_result |>
  dplyr::filter(gamma_0 == -5) |>
  dplyr::mutate(gamma_0 = sprintf("r = %d", gamma_0)) |>
  ggplot(aes(x = beta_0,
             y = mean,
             color = method)) +
  ggh4x::facet_grid2(gamma_0 ~ metric, scales = "free_y") +
  geom_point(size = 0.8) +
  geom_line(linewidth = 0.5) +
  geom_hline(aes(yintercept = level), linetype = "dashed", color = "red") +
  scale_x_continuous() +
  scale_y_log10() +
  theme_bw() +
  labs(x = "Varying b") +
  theme(strip.text.x = element_blank(),
        strip.text.y = element_text(margin = margin(0,0.08,0,0.08, "cm")),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none"
  )


beta_plot <- Type_I_result |>
  dplyr::filter(beta_0 == -5) |>
  dplyr::mutate(beta_0 = sprintf("b = %d", beta_0)) |>
  ggplot(aes(x = gamma_0,
             y = mean,
             color = method)) +
  scale_x_continuous() +
  scale_y_log10() +
  ggh4x::facet_grid2(beta_0 ~ metric, scales = "free_y") +
  geom_point(size = 0.8) +
  geom_line(linewidth = 0.5) +
  geom_hline(aes(yintercept = level), linetype = "dashed", color = "red") +
  theme_bw() +
  labs(x = "Varying r") +
  theme(strip.text.x = element_blank(),
        strip.text.y = element_text(margin = margin(0,0.08,0,0.08, "cm")),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none"
  )


plot_aux <- Type_I_result |>
  dplyr::filter(beta_0 == -5) |>
  dplyr::mutate(beta_0 = sprintf("b = %d", beta_0)) |>
  ggplot(aes(x = gamma_0,
             y = mean,
             color = method)) +
  scale_x_continuous() +
  scale_y_log10() +
  ggh4x::facet_grid2(metric ~ beta_0, scales = "free_y") +
  geom_point(size = 0.5) +
  geom_line() +
  geom_hline(aes(yintercept = level), linetype = "dashed", color = "red") +
  theme_bw() +
  theme(strip.text.x = element_text(margin = margin(0.08,0,0.08,0,"cm")),
        strip.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="bottom"
  )

# get legend from plot_aux
legend <- cowplot::get_plot_component(plot_aux, 'guide-box-bottom', return_all = TRUE)

# combine the plot
plot <- cowplot::plot_grid(gamma_plot, beta_plot,
                           nrow = 2, align = "v")

#create common x and y labels
y.grob <- textGrob("Type I error",
                   gp=gpar(col="black"), rot=90)

#add to plot
final_plot <- grid.arrange(arrangeGrob(plot, left = y.grob, nrow=1),
                           cowplot::plot_grid(legend, nrow = 1),
                           top = textGrob(TeX(paste0("Right-sided test with \\theta = ",
                                                     disp)),
                                          gp=grid::gpar(fontface=2)),
                           nrow=2, heights=c(9, 1))

ggsave(paste0(sprintf('%s/Type-I-error/plot-bin-NB-normal-B-50000-n-5000-5e3-n5-n5-disp-',
                      plots_dir),
              disp.char,'-Type-I-error-RIGHT.pdf'),
       plot = final_plot,
       width = 6, height = 4.5)

print("Type-I error plots generated!")

rm(Type_I_result, gamma_plot, beta_plot, plot_aux, plot, final_plot)

###########################################################################################
## Power ##################################################################################
###########################################################################################


# The following chunk of code is to produce left-sided plot; this code is for Figure 12

results_rds_metrics <- results_rds$metrics %>% filter(metric %in% c('level.1.L',
                                                                    'level.2.L',
                                                                    'level.3.L'))

power_result <- results_rds_metrics |>
  dplyr::filter(rho <= 0) %>%
  dplyr::mutate(metric = factor(metric)) |>
  dplyr::mutate(metric = forcats::fct_recode(metric,
                                             "alpha = 0.05" = "level.1.L",
                                             "alpha = 0.01" = "level.2.L",
                                             "alpha = 0.005" = "level.3.L")) |>
  dplyr::mutate(method = factor(method)) |>
  dplyr::mutate(method = forcats::fct_recode(method,
                                             "score" = "score.test")) |>
  dplyr::mutate(method = factor(method, levels = c("GCM",
                                                   "score",
                                                   "dCRT",
                                                   "spaCRT")))

######################
power_result <- power_result |>
  dplyr::mutate(
    level = case_when(
      metric == "alpha = 0.05" ~ 0.05,
      metric == "alpha = 0.01" ~ 0.01,
      metric == "alpha = 0.005" ~ 0.005,
      TRUE ~ NA_integer_  # Integer output for default case (NA as integer)
    )
  ) |>
  dplyr::mutate(mean = ifelse(mean == 0, 1/B, mean),
                sd = ifelse(se == 0, sqrt((1 - 1/B)) / B, se))


# varying beta
gamma_rho_plot <- power_result |>
  dplyr::filter(gamma_0 == -5 & beta_0 >= -4) |>
  dplyr::mutate(beta_0 = sprintf("b = %d", beta_0)) |>
  dplyr::mutate(metric.num = dplyr::recode(metric,
                                           "alpha = 0.05" = 0.05,
                                           "alpha = 0.01" = 0.01,
                                           "alpha = 0.005" = 0.005)) %>%
  dplyr::mutate(metric = sprintf("alpha == %.3f", level)) |>
  ggplot(aes(x = rho,
             y = mean,
             color = method)) +
  scale_x_continuous() +
  scale_y_log10() +
  ggh4x::facet_grid2(metric ~ beta_0, scales = "free_y", labeller = labeller(
    metric = label_parsed  # Parse LaTeX-like expressions for 'eps'
  )) +
  geom_point(size = 0.5) +
  geom_line() +
  geom_hline(aes(yintercept = metric.num), linetype = "dashed", color = "red") +
  theme_bw() +
  labs(x = TeX("varying \\rho"), y = "Rejection rate",
       title = TeX(paste0("Left-sided test with \\theta = ",
                          disp))) +
  theme(strip.text.x = element_text(margin = margin(0.08, 0, 0.08, 0, "cm")),
        strip.text.y = element_text(margin = margin(0,0.08,0,0.08, "cm")),
        legend.position="none",
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

# varying gamma
beta_rho_plot <- power_result |>
  dplyr::filter(beta_0 == -5 & gamma_0 >= -4) |>
  dplyr::mutate(gamma_0 = sprintf("r = %d", gamma_0)) |>
  dplyr::mutate(metric.num = dplyr::recode(metric,
                                           "alpha = 0.05" = 0.05,
                                           "alpha = 0.01" = 0.01,
                                           "alpha = 0.005" = 0.005)) %>%
  dplyr::mutate(metric = sprintf("alpha == %.3f", level)) |>
  ggplot(aes(x = rho,
             y = mean,
             color = method)) +
  scale_x_continuous() +
  scale_y_log10() +
  ggh4x::facet_grid2(metric ~ gamma_0, scales = "free_y", labeller = labeller(
    metric = label_parsed  # Parse LaTeX-like expressions for 'eps'
  )) +
  geom_point(size = 0.5) +
  geom_line() +
  geom_hline(aes(yintercept = metric.num), linetype = "dashed", color = "red") +
  theme_bw() +
  labs(x = TeX("varying \\rho"), y = "Rejection rate",
       title = TeX(paste0("Left-sided test with b = -5, \\theta = ",
                          disp))) +
  theme(strip.text.x = element_text(margin = margin(0.08,0,0.08,0, "cm")),
        strip.text.y = element_text(margin = margin(0,0.08,0,0.08, "cm")),
        legend.position="bottom",
        panel.grid.minor = element_blank(),
        plot.title = element_blank())

# same plots can be done for right-sided plot.

ggsave(paste0(sprintf('%s/power/plot-bin-NB-normal-B-50000-n-5000-5e3-n5-n5-disp-',
                      plots_dir),
              disp.char,'-power-fixed-gamma-LEFT.pdf'),
       plot = gamma_rho_plot,
       width = 4.5, height = 4.5)

ggsave(paste0(sprintf('%s/power/plot-bin-NB-normal-B-50000-n-5000-5e3-n5-n5-disp-',
                      plots_dir),
              disp.char,'-power-fixed-beta-LEFT.pdf'),
       plot = beta_rho_plot,
       width = 4.5, height = 4.5)

rm(power_result, gamma_rho_plot, beta_rho_plot)

# The following code is used to produce the right-sided test; this code is for Figure 13

results_rds_metrics <- results_rds$metrics %>% filter(metric %in% c('level.1.R',
                                                                    'level.2.R',
                                                                    'level.3.R'))

power_result <- results_rds_metrics |>
  dplyr::filter(rho >= 0) %>%
  dplyr::mutate(metric = factor(metric)) |>
  dplyr::mutate(metric = forcats::fct_recode(metric,
                                             "alpha = 0.05" = "level.1.R",
                                             "alpha = 0.01" = "level.2.R",
                                             "alpha = 0.005" = "level.3.R")) |>
  dplyr::mutate(method = factor(method)) |>
  dplyr::mutate(method = forcats::fct_recode(method,
                                             "score" = "score.test")) |>
  dplyr::mutate(method = factor(method, levels = c("GCM",
                                                   "score",
                                                   "dCRT",
                                                   "spaCRT")))

######################
power_result <- power_result |>
  dplyr::mutate(
    level = case_when(
      metric == "alpha = 0.05" ~ 0.05,
      metric == "alpha = 0.01" ~ 0.01,
      metric == "alpha = 0.005" ~ 0.005,
      TRUE ~ NA_integer_  # Integer output for default case (NA as integer)
    )
  ) |>
  dplyr::mutate(mean = ifelse(mean == 0, 1/B, mean),
                sd = ifelse(se == 0, sqrt((1 - 1/B)) / B, se))

# varying beta
gamma_rho_plot <- power_result |>
  dplyr::filter(gamma_0 == -5 & beta_0 >= -4) |>
  dplyr::mutate(beta_0 = sprintf("b = %d", beta_0)) |>
  dplyr::mutate(metric.num = dplyr::recode(metric,
                                           "alpha = 0.05" = 0.05,
                                           "alpha = 0.01" = 0.01,
                                           "alpha = 0.005" = 0.005)) %>%
  dplyr::mutate(metric = sprintf("alpha == %.3f", level)) |>
  ggplot(aes(x = rho,
             y = mean,
             color = method)) +
  scale_x_continuous() +
  scale_y_log10() +
  ggh4x::facet_grid2(metric ~ beta_0, scales = "free_y", labeller = labeller(
    metric = label_parsed  # Parse LaTeX-like expressions for 'eps'
  )) +
  geom_point(size = 0.5) +
  geom_line() +
  geom_hline(aes(yintercept = metric.num), linetype = "dashed", color = "red") +
  theme_bw() +
  labs(x = TeX("varying \\rho"), y = "Rejection rate",
       title = TeX(paste0("Right-sided test with \\theta = ",
                          disp))) +
  theme(strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")),
        strip.text.y = element_text(margin = margin(0,0.1,0,0.1, "cm")),
        legend.position="none",
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

# varying gamma
beta_rho_plot <- power_result |>
  dplyr::filter(beta_0 == -5 & gamma_0 >= -4) |>
  dplyr::mutate(gamma_0 = sprintf("r = %d", gamma_0)) |>
  dplyr::mutate(metric.num = dplyr::recode(metric,
                                           "alpha = 0.05" = 0.05,
                                           "alpha = 0.01" = 0.01,
                                           "alpha = 0.005" = 0.005)) %>%
  dplyr::mutate(metric = sprintf("alpha == %.3f", level)) |>
  ggplot(aes(x = rho,
             y = mean,
             color = method)) +
  scale_x_continuous() +
  scale_y_log10() +
  ggh4x::facet_grid2(metric ~ gamma_0, scales = "free_y", labeller = labeller(
    metric = label_parsed  # Parse LaTeX-like expressions for 'eps'
  )) +
  geom_point(size = 0.5) +
  geom_line() +
  geom_hline(aes(yintercept = metric.num), linetype = "dashed", color = "red") +
  theme_bw() +
  labs(x = TeX("varying \\rho"), y = "Rejection rate",
       title = TeX(paste0("Right-sided test with b = -5, \\theta = ",
                          disp))) +
  theme(strip.text.x = element_text(margin = margin(0.1,0,0.1,0,"cm")),
        strip.text.y = element_text(margin = margin(0,0.1,0,0.1, "cm")),
        legend.position="bottom",
        panel.grid.minor = element_blank(),
        plot.title = element_blank())

# same plots can be done for right-sided plot.

ggsave(paste0(sprintf('%s/power/plot-bin-NB-normal-B-50000-n-5000-5e3-n5-n5-disp-',
                      plots_dir),
              disp.char,'-power-fixed-gamma-RIGHT.pdf'),
       plot = gamma_rho_plot,
       width = 4.5, height = 4.5)

ggsave(paste0(sprintf('%s/power/plot-bin-NB-normal-B-50000-n-5000-5e3-n5-n5-disp-',
                      plots_dir),
              disp.char,'-power-fixed-beta-RIGHT.pdf'),
       plot = beta_rho_plot,
       width = 4.5, height = 4.5)

print("Power plots generated!")

rm(power_result, gamma_rho_plot, beta_rho_plot)

###########################################################################################
## QQ-plots ###############################################################################
###########################################################################################


varying_params <- list(
  n = 5000,
  gamma_0 = seq.int(from = -6, to = -2),
  beta_0 = seq.int(from = -6, to = -2))

baseline_params <- list(n = 5e3, gamma_0 = -5, beta_0 = -5)

grid_ffd <- simulatr::create_param_grid_fractional_factorial(varying_params,
                                                             baseline_params)


rho <- c(seq.int(from = -4, to = 0), seq.int(from = 1, to = 4)/2)


parameter_grid <- grid_ffd %>% dplyr::select(1:3) %>%
  slice(rep(1:n(), each = length(rho))) %>%
  cbind(rho = rep(rho, times = nrow(grid_ffd))) %>%
  dplyr::mutate(grid_id = seq(1, length(rho), 1))

parameter_grid_null <- parameter_grid %>% filter(rho == 0)

results_rds_null <- results_rds$results |> filter(grid_id %in%
                                                    as.integer(parameter_grid_null$grid_id))

p_value_df_null <- results_rds_null |>
  dplyr::arrange(method) %>%
  tidyr::unnest(output) %>%
  dplyr::mutate(row_id = 1:nrow(.))

row_id_GCM <- p_value_df_null %>% filter(method == 'GCM') %>% pull(row_id)
row_id_dCRT <- p_value_df_null %>% filter(method == 'dCRT') %>% pull(row_id)
row_id_spaCRT <- p_value_df_null %>% filter(method == 'spaCRT') %>% pull(row_id)
row_id_score <- p_value_df_null %>% filter(method == 'score.test') %>% pull(row_id)

row_id_GCM_pvalue <- row_id_GCM[sort(c(seq(2,length(row_id_GCM),
                                           by = length(row_id_GCM)/(B*nrow(grid_ffd))),
                                       seq(3,length(row_id_GCM),
                                           by = length(row_id_GCM)/(B*nrow(grid_ffd)))))]
row_id_dCRT_pvalue <- row_id_dCRT[sort(c(seq(2,length(row_id_dCRT),
                                             by = length(row_id_dCRT)/(B*nrow(grid_ffd))),
                                         seq(3,length(row_id_dCRT),
                                             by = length(row_id_dCRT)/(B*nrow(grid_ffd)))))]
row_id_spaCRT_pvalue <- row_id_spaCRT[sort(c(seq(2,length(row_id_spaCRT),
                                                 by = length(row_id_spaCRT)/(B*nrow(grid_ffd))),
                                             seq(3,length(row_id_spaCRT),
                                                 by = length(row_id_spaCRT)/(B*nrow(grid_ffd)))))]
row_id_score_pvalue <- row_id_score[sort(c(seq(2,length(row_id_score),
                                               by = length(row_id_score)/(B*nrow(grid_ffd))),
                                           seq(3,length(row_id_score),
                                               by = length(row_id_score)/(B*nrow(grid_ffd)))))]

rm(row_id_GCM, row_id_dCRT, row_id_spaCRT, row_id_score)

p_value_df_null <- p_value_df_null %>%
  filter(row_id %in% c(row_id_GCM_pvalue,
                       row_id_dCRT_pvalue,
                       row_id_spaCRT_pvalue,
                       row_id_score_pvalue)) %>%
  arrange(row_id) |>
  dplyr::mutate(sideness = rep(c("left", "right"),
                               B * length(unique(p_value_df_null$method)) *
                                 nrow(parameter_grid_null)),
                p_value = unlist(output)) |>
  dplyr::select(-output) %>%
  dplyr::left_join(parameter_grid_null, by = "grid_id") |>
  dplyr::mutate(method = forcats::fct_recode(method,
                                             "score" = "score.test")) |>
  dplyr::mutate(method = factor(method, levels = c("GCM",
                                                   "score",
                                                   "dCRT",
                                                   "spaCRT")))

my_theme <- theme_bw() + theme(axis.line = element_line(color = "black"),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.border = element_blank(),
                               panel.background = element_blank(),
                               plot.title = element_text(hjust = 0.5, size=11))

################################################################################
# The following code is to plot the left-sided test; this code is for Figure 10
################################################################################

plot_gamma <- p_value_df_null |>
  dplyr::filter(gamma_0 == default_gamma &
                  sideness == "left" &
                  beta_0 >= -4) |>
  dplyr::mutate(beta_0 = sprintf("b = %d", beta_0)) |>
  ggplot(mapping = aes(y = p_value, color = method)) +
  stat_qq_points(ymin = 1e-8, size = 0.8) +
  stat_qq_band() +
  ggh4x::facet_grid2(. ~ beta_0, scales = "free", independent = "x") +
  scale_x_continuous(trans = revlog_trans(10),
                     breaks = c(1e-1, 1e-3, 1e-5)) +
  scale_y_continuous(trans = revlog_trans(10),
                     breaks = c(1e-1, 1e-3, 1e-5, 1e-7)) +
  geom_abline(col = "black") +
  my_theme +
  theme(strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")),
        strip.text.y = element_text(margin = margin(0,0.1,0,0.1, "cm")),
        legend.position = "none",
        legend.margin = margin(t = -0.2, r = -0.2, unit = 'cm'),
        legend.text = element_text(size = 8),  # Change legend text size
        legend.title = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm")) +
  guides(color = guide_legend(
    keywidth = 0.0,
    keyheight = 0.15,
    default.unit = "inch",
    override.aes = list(size = 2.5)))

plot_beta <- p_value_df_null |>
  dplyr::filter(beta_0 == default_beta &
                  sideness == "left" &
                  gamma_0 >= -4) |>
  dplyr::mutate(gamma_0 = sprintf("r = %d", gamma_0)) |>
  ggplot(mapping = aes(y = p_value, color = method)) +
  stat_qq_points(ymin = 1e-8, size = 0.8) +
  stat_qq_band() +
  ggh4x::facet_grid2(. ~ gamma_0, scales = "free", independent = "x") +
  scale_x_continuous(trans = revlog_trans(10),
                     breaks = c(1e-1, 1e-3, 1e-5)) +
  scale_y_continuous(trans = revlog_trans(10),
                     breaks = c(1e-1, 1e-3, 1e-5, 1e-7)) +
  geom_abline(col = "black") +
  my_theme +
  theme(strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")),
        strip.text.y = element_text(margin = margin(0,0.1,0,0.1, "cm")),
        legend.position = "none",
        legend.margin = margin(t = -0.2, r = -0.2, unit = 'cm'),
        legend.text = element_text(size = 8),  # Change legend text size
        legend.title = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm")) +
  guides(color = guide_legend(
    keywidth = 0.0,
    keyheight = 0.15,
    default.unit = "inch",
    override.aes = list(size = 2.5)))

plot_aux <- p_value_df_null |>
  dplyr::filter(beta_0 == -5 &
                  sideness == "left" &
                  gamma_0 >= -4) |>
  ggplot(mapping = aes(y = p_value, color = method)) +
  stat_qq_points(ymin = 1e-8, size = 0.8) +
  stat_qq_band() +
  ggh4x::facet_grid2(. ~ gamma_0, scales = "free", independent = "x") +
  scale_x_continuous(trans = revlog_trans(10),
                     breaks = c(1e-1, 1e-3, 1e-5)) +
  scale_y_continuous(trans = revlog_trans(10),
                     breaks = c(1e-1, 1e-3, 1e-5, 1e-7)) +
  geom_abline(col = "black") +
  my_theme +
  theme(strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")),
        strip.text.y = element_text(margin = margin(0,0.1,0,0.1, "cm")),
        legend.position = "bottom",
        legend.margin=margin(t = -0.2, r = -0.2, unit = 'cm'),
        legend.text = element_text(size = 8),  # Change legend text size
        legend.title = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm")) +
  guides(color = guide_legend(
    keywidth = 0.0,
    keyheight = 0.15,
    default.unit = "inch",
    override.aes = list(size = 2.5)))

# combine the plot
plot <- cowplot::plot_grid(plot_gamma, plot_beta,
                           nrow = 2, align = "h")

# get legend from plot_aux
legend <- cowplot::get_plot_component(plot_aux, 'guide-box-bottom', return_all = TRUE)

#create common x and y labels
y.grob <- textGrob("Observed p-value",
                   gp=gpar(col="black"), rot=90)

x.grob <- textGrob("Expected null p-value",
                   gp=gpar(col="black"))


#add to plot
final_plot <- grid.arrange(arrangeGrob(plot, left = y.grob, bottom = x.grob, nrow=1),
                           plot_grid(legend, nrow = 1),
                           nrow=2, heights=c(9, 1))


ggsave(paste0(sprintf('%s/QQ/plot-bin-NB-normal-B-50000-n-5000-5e3-n5-n5-disp-',
                      plots_dir),
              disp.char,'-QQ-LEFT.pdf'),
       plot = final_plot,
       width = 6, height = 4.5)

rm(plot_gamma, plot_beta, plot_aux, plot, final_plot)

################################################################################
# The following code is to plot the right-sided test; this code is for Figure 11
################################################################################
plot_gamma <- p_value_df_null |>
  dplyr::filter(gamma_0 == default_gamma &
                  sideness == "right" &
                  beta_0 >= -4) |>
  dplyr::mutate(beta_0 = sprintf("b = %d", beta_0)) |>
  ggplot(mapping = aes(y = p_value, color = method)) +
  stat_qq_points(ymin = 1e-8, size = 0.8) +
  stat_qq_band() +
  ggh4x::facet_grid2(. ~ beta_0, scales = "free", independent = "x") +
  scale_x_continuous(trans = revlog_trans(10),
                     breaks = c(1e-1, 1e-3, 1e-5)) +
  scale_y_continuous(trans = revlog_trans(10),
                     breaks = c(1e-1, 1e-3, 1e-5, 1e-7)) +
  geom_abline(col = "black") +
  my_theme +
  theme(strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")),
        strip.text.y = element_text(margin = margin(0,0.1,0,0.1, "cm")),
        legend.position = "none",
        legend.margin = margin(t = -0.2, r = -0.2, unit = 'cm'),
        legend.text = element_text(size = 8),  # Change legend text size
        legend.title = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm")) +
  guides(color = guide_legend(
    keywidth = 0.0,
    keyheight = 0.15,
    default.unit = "inch",
    override.aes = list(size = 2.5)))

plot_beta <- p_value_df_null |>
  dplyr::filter(beta_0 == default_beta &
                  sideness == "right" &
                  gamma_0 >= -4) |>
  dplyr::mutate(gamma_0 = sprintf("r = %d", gamma_0)) |>
  ggplot(mapping = aes(y = p_value, color = method)) +
  stat_qq_points(ymin = 1e-8, size = 0.8) +
  stat_qq_band() +
  ggh4x::facet_grid2(. ~ gamma_0, scales = "free", independent = "x") +
  scale_x_continuous(trans = revlog_trans(10),
                     breaks = c(1e-1, 1e-3, 1e-5)) +
  scale_y_continuous(trans = revlog_trans(10),
                     breaks = c(1e-1, 1e-3, 1e-5, 1e-7)) +
  geom_abline(col = "black") +
  my_theme +
  theme(strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")),
        strip.text.y = element_text(margin = margin(0,0.1,0,0.1, "cm")),
        legend.position = "none",
        legend.margin = margin(t = -0.2, r = -0.2, unit = 'cm'),
        legend.text = element_text(size = 8),  # Change legend text size
        legend.title = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm")) +
  guides(color = guide_legend(
    keywidth = 0.0,
    keyheight = 0.15,
    default.unit = "inch",
    override.aes = list(size = 2.5)))

plot_aux <- p_value_df_null |>
  dplyr::filter(beta_0 == default_beta &
                  sideness == "right" &
                  gamma_0 >= -4) |>
  ggplot(mapping = aes(y = p_value, color = method)) +
  stat_qq_points(ymin = 1e-8, size = 0.8) +
  stat_qq_band() +
  ggh4x::facet_grid2(. ~ gamma_0, scales = "free", independent = "x") +
  scale_x_continuous(trans = revlog_trans(10),
                     breaks = c(1e-1, 1e-3, 1e-5)) +
  scale_y_continuous(trans = revlog_trans(10),
                     breaks = c(1e-1, 1e-3, 1e-5, 1e-7)) +
  geom_abline(col = "black") +
  my_theme +
  theme(strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")),
        strip.text.y = element_text(margin = margin(0,0.1,0,0.1, "cm")),
        legend.position = "bottom",
        legend.margin=margin(t = -0.2, r = -0.2, unit = 'cm'),
        legend.text = element_text(size = 8),  # Change legend text size
        legend.title = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm")) +
  guides(color = guide_legend(
    keywidth = 0.0,
    keyheight = 0.15,
    default.unit = "inch",
    override.aes = list(size = 2.5)))

# combine the plot
plot <- cowplot::plot_grid(plot_gamma, plot_beta,
                           nrow = 2, align = "h")

# get legend from plot_aux
legend <- cowplot::get_plot_component(plot_aux, 'guide-box-bottom', return_all = TRUE)

#create common x and y labels
y.grob <- textGrob("Observed p-value",
                   gp=gpar(col="black"), rot=90)

x.grob <- textGrob("Expected null p-value",
                   gp=gpar(col="black"))


#add to plot
final_plot <- grid.arrange(arrangeGrob(plot, left = y.grob, bottom = x.grob, nrow=1),
                           plot_grid(legend, nrow = 1),
                           nrow=2, heights=c(9, 1))


ggsave(paste0(sprintf('%s/QQ/plot-bin-NB-normal-B-50000-n-5000-5e3-n5-n5-disp-',
                      plots_dir),
              disp.char,'-QQ-RIGHT.pdf'),
       plot = final_plot,
       width = 6, height = 4.5)

print("QQ plots generated!")

rm(plot_gamma, plot_beta, plot_aux, plot, final_plot)



################################################################################
# The following code is to generate plot for comparing dCRT and spaCRT p-values (Figure 22, 23)
################################################################################

# rename the left and right sided tests
p_value_df_null <- p_value_df_null |>
  mutate(sideness = if_else(sideness == "left", "left-sided", "right-sided"))

# set the seed
set.seed(1)

# plot the line and dot figure
partial_line_plot <- p_value_df_null |>
  dplyr::filter(gamma_0 == -5 & beta_0 == -5) |>
  tidyr::pivot_wider(id_cols =  c("sideness", "run_id"),
                     names_from = "method",
                     values_from = "p_value") |>
  slice_sample(n = 20000) |>
  ggplot(aes(x = dCRT, y = spaCRT)) +
  scale_x_continuous(trans = revlog_trans(10),
                     breaks = c(1e-1, 1e-2, 1e-3)) +
  scale_y_continuous(trans = revlog_trans(10),
                     breaks = c(1e-1, 1e-2, 1e-3)) +
  geom_abline(col = "black") +
  facet_wrap(.~sideness) +
  geom_point(size = 0.8) + geom_abline() + my_theme

# save the plot
ggsave(paste0(sprintf('%s/QQ/plot-bin-NB-normal-B-50000-n-5000-5e3-n5-n5-disp-',
                      plots_dir),
              disp.char,'partial-dCRT-spaCRT.pdf'),
       plot = partial_line_plot,
       width = 6, height = 3)

# full line plot for varying gamma
full_line_plot_1 <- p_value_df_null |>
  dplyr::filter(beta_0 == -5) |>
  tidyr::pivot_wider(id_cols =  c("gamma_0", "sideness", "run_id"),
                     names_from = "method",
                     values_from = "p_value") |>
  dplyr::mutate(gamma_0 = sprintf("gamma[0] == %d", gamma_0)) |>
  slice_sample(n = 20000) |>
  ggplot(aes(x = dCRT, y = spaCRT)) +
  scale_x_continuous(trans = revlog_trans(10),
                     breaks = c(1e-1, 1e-2, 1e-3)) +
  scale_y_continuous(trans = revlog_trans(10),
                     breaks = c(1e-1, 1e-2, 1e-3)) +
  geom_abline(col = "black") +
  ggh4x::facet_grid2(gamma_0~sideness,  labeller = labeller(
    gamma_0 = label_parsed  # Parse LaTeX-like expressions for 'eps'
  )) +
  geom_point(size = 0.8) + my_theme

# save the plot
ggsave(paste0(sprintf('%s/QQ/plot-bin-NB-normal-B-50000-n-5000-5e3-n5-n5-disp-',
                      plots_dir),
              disp.char,'full-dCRT-spaCRT-varying-gamma.pdf'),
       plot = full_line_plot_1,
       width = 4, height = 7)

# full line plot for varying beta
full_line_plot_2 <- p_value_df_null |>
  dplyr::filter(gamma_0 == -5) |>
  tidyr::pivot_wider(id_cols =  c("beta_0", "sideness", "run_id"),
                     names_from = "method",
                     values_from = "p_value") |>
  dplyr::mutate(beta_0 = sprintf("beta[0] == %d", beta_0)) |>
  slice_sample(n = 20000) |>
  ggplot(aes(x = dCRT, y = spaCRT)) +
  scale_x_continuous(trans = revlog_trans(10),
                     breaks = c(1e-1, 1e-2, 1e-3)) +
  scale_y_continuous(trans = revlog_trans(10),
                     breaks = c(1e-1, 1e-2, 1e-3)) +
  geom_abline(col = "black") +
  facet_grid(beta_0~sideness,  labeller = labeller(
    beta_0 = label_parsed  # Parse LaTeX-like expressions for 'eps'
  )) +
  geom_point(size = 0.8) + my_theme

# save the plot
ggsave(paste0(sprintf('%s/QQ/plot-bin-NB-normal-B-50000-n-5000-5e3-n5-n5-disp-',
                      plots_dir),
              disp.char,'full-dCRT-spaCRT-varying-beta.pdf'),
       plot = full_line_plot_2,
       width = 4, height = 7)
