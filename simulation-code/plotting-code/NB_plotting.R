library(dplyr)
library(ggplot2)
library(katlabutils)
library(cowplot)
library(grid)
library(kableExtra)
library(latex2exp)
library(gridExtra)
library(patchwork)

# read simulation results
source(".Rprofile")
path_rds <- paste0(.get_config_path("LOCAL_SPACRT_DATA_DIR"),
                   'formal-simulation-results/',
                   "negbinomial_regression/negbinomial_regression_results.rds")
results_rds <- readRDS(path_rds)

# load the parameter grid
parameter_grid <- readRDS("simulation-code/sim-spec-formal/NB_parameter_grid.rds")

# specify the directory to be saved
plots_dir <- "manuscript/figures-and-tables//simulation/NB-regression"

# specify plotting theme
my_theme <- theme_bw() + theme(axis.line = element_line(color = "black"),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.border = element_blank(),
                               panel.background = element_blank(),
                               plot.title = element_text(hjust = 0.5, size=11))

# obtain the default gamma_0 and beta_0 parameters
default_gamma <- parameter_grid |> dplyr::filter(!arm_gamma_0) |>
  dplyr::select(gamma_0) |>
  dplyr::pull() |> unique()
default_beta <- parameter_grid |> dplyr::filter(!arm_beta_0) |>
  dplyr::select(beta_0) |>
  dplyr::pull() |> unique()

########################## 0. Create the directory #############################
QQ_plot_dir <- sprintf("%s/QQ", plots_dir)
power_plot_dir <- sprintf("%s/power", plots_dir)
FDR_plot_dir <- sprintf("%s/FDR", plots_dir)
summary_plot_dir <- sprintf("%s/summary", plots_dir)
for (plot_dir in c(QQ_plot_dir, power_plot_dir, FDR_plot_dir, summary_plot_dir)) {
  if(dir.exists(plot_dir)){
    next
  }else{
    dir.create(plot_dir)
  }
}

########################## 1. Summarize the results ############################
joined_results <- results_rds$results |>
  dplyr::left_join(parameter_grid, by = "grid_id") |>
  dplyr::mutate(method = forcats::fct_recode(method, "score" = "score.test")) |>
  dplyr::mutate(method = factor(method, levels = c("GCM", "score", "dCRT", "spaCRT")))

# specify the significance level
alpha <- 0.005
FDR <- 0.1
B <- max(joined_results$run_id)

# compute the rates of spaCRT using GCM default
joined_results |>
  dplyr::filter(method == "spaCRT") |>
  dplyr::select(theta, rho, gamma_0, beta_0, output) |>
  tidyr::unnest_wider(output) |>
  dplyr::select(-computation_time) |>
  tidyr::unnest_wider(pvalue_list) |>
  dplyr::group_by(theta, rho, gamma_0, beta_0) |>
  dplyr::summarise(GCM_default = 1 - mean(spa.success)) |>
  dplyr::ungroup() |>
  dplyr::summarise(max_GCM_default = max(GCM_default)) |>
  dplyr::pull()

########################## 2. Loop over dispersion parameter ###################
for(size_parameter in unique(parameter_grid$theta)){

  # filter the results
  current_results <- joined_results |>
    dplyr::filter(theta == size_parameter) |>
    tidyr::unnest_wider(output)

  ############################ 3. QQ-plots #####################################
  # extract the results for varying gamma
  varying_beta_null <- current_results |>
    dplyr::filter(gamma_0 == default_gamma & beta_0 >= -4 & rho == 0) |>
    dplyr::mutate(beta_0 = sprintf("beta[0] == %d", beta_0)) |>
    tidyr::unnest_wider(pvalue_list)

  # extract the results for varying beta
  varying_gamma_null <- current_results |>
    dplyr::filter(beta_0 == default_beta & gamma_0 >= -4 & rho == 0) |>
    dplyr::mutate(gamma_0 = sprintf("gamma[0] == %d", gamma_0)) |>
    tidyr::unnest_wider(pvalue_list)

  ############################ 3.1 Left-sided test #############################
  # varying beta_0 parameter
  plot_beta <- varying_beta_null |>
    ggplot(mapping = aes(y = p.left, color = method)) +
    stat_qq_points(ymin = 1e-8, size = 0.8) +
    stat_qq_band() +
    ggh4x::facet_grid2(. ~ beta_0, scales = "free", independent = "x",
                       labeller = labeller(
                         beta_0 = label_parsed  # Parse LaTeX-like expressions
                       )) +
    scale_x_continuous(trans = revlog_trans(10),
                       breaks = c(1e-1, 1e-3, 1e-5)) +
    scale_y_continuous(trans = revlog_trans(10),
                       breaks = c(1e-1, 1e-3, 1e-5, 1e-7)) +
    geom_abline(col = "black") +
    my_theme +
    theme(strip.text.x = element_text(margin = margin(0.05, 0, 0.05, 0, "cm"), size = 12),
          strip.text.y = element_text(margin = margin(0, 0.05, 0, 0.05, "cm"), size = 12),
          legend.position = "none",
          axis.text.y = element_text(size = 12),
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

  # varying gamma_0 plot
  plot_gamma <- varying_gamma_null |>
    ggplot(mapping = aes(y = p.left, color = method)) +
    stat_qq_points(ymin = 1e-8, size = 0.8) +
    stat_qq_band() +
    ggh4x::facet_grid2(. ~ gamma_0, scales = "free", independent = "x",
                       labeller = labeller(
                         gamma_0 = label_parsed  # Parse LaTeX-like expressions
                       )) +
    scale_x_continuous(trans = revlog_trans(10), breaks = c(1e-1, 1e-3, 1e-5)) +
    scale_y_continuous(trans = revlog_trans(10), breaks = c(1e-1, 1e-3, 1e-5, 1e-7)) +
    geom_abline(col = "black") +
    my_theme +
    theme(strip.text.x = element_text(margin = margin(0.05, 0, 0.05, 0, "cm"), size = 12),
          strip.text.y = element_text(margin = margin(0, 0.05, 0, 0.05, "cm"), size = 12),
          legend.position = "none",
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm")) +
    guides(color = guide_legend(
      keywidth = 0.0,
      keyheight = 0.15,
      default.unit = "inch",
      override.aes = list(size = 2.5)
    ))

  # combine the plot
  plot <- cowplot::plot_grid(plot_gamma, plot_beta, nrow = 2, align = "h")

  # get legend from plot_aux
  legend <- cowplot::get_plot_component(plot_gamma +
                                          theme(legend.position = "bottom",
                                                legend.text = element_text(size = 12)),
                                        'guide-box-bottom', return_all = TRUE)

  #create common x and y labels
  y.grob <- textGrob("Observed p-value", gp = gpar(col = "black", fontsize = 15), rot = 90)
  x.grob <- textGrob("Expected null p-value", gp = gpar(col = "black", fontsize = 15))

  # compile final plot
  final_plot <- grid.arrange(arrangeGrob(plot, left = y.grob, bottom = x.grob, nrow=1),
                             plot_grid(legend, nrow = 1), nrow = 2, heights=c(9, 1))

  # save the plot
  ggsave(paste0(sprintf('%s/disp-', QQ_plot_dir), size_parameter,'-QQ-LEFT.pdf'),
         plot = final_plot, width = 6, height = 4.5)

  ########################## 3.2 Right-sided test ##############################
  # varying beta
  plot_beta <- varying_beta_null |>
    ggplot(mapping = aes(y = p.right, color = method)) +
    stat_qq_points(ymin = 1e-8, size = 0.8) +
    stat_qq_band() +
    ggh4x::facet_grid2(. ~ beta_0, scales = "free", independent = "x",
                       labeller = labeller(
                         beta_0 = label_parsed  # Parse LaTeX-like expressions
                       )) +
    scale_x_continuous(trans = revlog_trans(10),
                       breaks = c(1e-1, 1e-3, 1e-5)) +
    scale_y_continuous(trans = revlog_trans(10),
                       breaks = c(1e-1, 1e-3, 1e-5, 1e-7)) +
    geom_abline(col = "black") +
    my_theme +
    theme(strip.text.x = element_text(margin = margin(0.05, 0, 0.05, 0, "cm"), size = 12),
          strip.text.y = element_text(margin = margin(0, 0.05, 0, 0.05, "cm"), size = 12),
          legend.position = "none",
          axis.text.y = element_text(size = 12),
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

  # varying gamma_0 parameter
  plot_gamma <- varying_gamma_null |>
    ggplot(mapping = aes(y = p.right, color = method)) +
    stat_qq_points(ymin = 1e-8, size = 0.8) +
    stat_qq_band() +
    ggh4x::facet_grid2(. ~ gamma_0, scales = "free", independent = "x",
                       labeller = labeller(
                         gamma_0 = label_parsed  # Parse LaTeX-like expressions
                       )) +
    scale_x_continuous(trans = revlog_trans(10), breaks = c(1e-1, 1e-3, 1e-5)) +
    scale_y_continuous(trans = revlog_trans(10), breaks = c(1e-1, 1e-3, 1e-5, 1e-7)) +
    geom_abline(col = "black") +
    my_theme +
    theme(strip.text.x = element_text(margin = margin(0.05, 0, 0.05, 0, "cm"), size = 12),
          strip.text.y = element_text(margin = margin(0, 0.05, 0, 0.05, "cm"), size = 12),
          legend.position = "none",
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm")) +
    guides(color = guide_legend(
      keywidth = 0.0,
      keyheight = 0.15,
      default.unit = "inch",
      override.aes = list(size = 2.5)
    ))

  # combine the plot
  plot <- cowplot::plot_grid(plot_gamma, plot_beta, nrow = 2, align = "h")

  # get legend from plot_aux
  legend <- cowplot::get_plot_component(plot_gamma +
                                          theme(legend.position = "bottom",
                                                legend.text = element_text(size = 12)),
                                        'guide-box-bottom', return_all = TRUE)

  #create common x and y labels
  y.grob <- textGrob("Observed p-value", gp = gpar(col = "black", fontsize = 15), rot = 90)
  x.grob <- textGrob("Expected null p-value", gp = gpar(col = "black", fontsize = 15))

  # compile final plot
  final_plot <- grid.arrange(arrangeGrob(plot, left = y.grob, bottom = x.grob, nrow=1),
                             plot_grid(legend, nrow = 1), nrow = 2, heights=c(9, 1))

  # save the plot
  ggsave(paste0(sprintf('%s/disp-', QQ_plot_dir), size_parameter,'-QQ-RIGHT.pdf'),
         plot = final_plot, width = 6, height = 4.5)

  ########################### 3. Power-plots ###################################
  # extract the results for varying gamma
  varying_gamma_power <- current_results |>
    dplyr::filter(beta_0 == default_beta & gamma_0 >= -4) |>
    dplyr::mutate(gamma_0 = sprintf("gamma[0] == %d", gamma_0)) |>
    tidyr::unnest_wider(pvalue_list) |>
    dplyr::group_by(gamma_0, rho, method) |>
    dplyr::summarise(
      rejection_left = max(1 / B, mean(p.left <= alpha)),
      rejection_right = max(1 / B, mean(p.right <= alpha))
    ) |>
    dplyr::ungroup()

  # extract the results for varying beta
  varying_beta_power <- current_results |>
    dplyr::filter(gamma_0 == default_gamma & beta_0 >= -4) |>
    dplyr::mutate(beta_0 = sprintf("beta[0] == %d", beta_0)) |>
    tidyr::unnest_wider(pvalue_list) |>
    dplyr::group_by(beta_0, rho, method) |>
    dplyr::summarise(
      rejection_left = max(1 / B, mean(p.left <= alpha)),
      rejection_right = max(1 / B, mean(p.right <= alpha))
    ) |>
    dplyr::ungroup()

  ######################### 3.1 Left-sided test ################################
  # varying gamma
  plot_gamma <- varying_gamma_power |>
    dplyr::filter(rho <= 0) |>
    ggplot(aes(x = rho, y = rejection_left, color = method)) +
    scale_x_continuous() +
    scale_y_log10(
      labels = scales::label_number(accuracy = 0.001)
    ) +
    ggh4x::facet_grid2(. ~ gamma_0, scales = "free_y",
                       labeller = labeller(
                         gamma_0 = label_parsed  # Parse LaTeX-like expressions
                       )) +
    geom_point(size = 0.5) +
    geom_line() +
    geom_hline(aes(yintercept = alpha), linetype = "dashed", color = "red") +
    theme_bw() +
    theme(strip.text.x = element_text(margin = margin(0.05, 0, 0.05, 0, "cm"), size = 12),
          strip.text.y = element_text(margin = margin(0, 0.05, 0, 0.05, "cm"), size = 12),
          legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank())

  # varying beta
  plot_beta <- varying_beta_power |>
    dplyr::filter(rho <= 0) |>
    ggplot(aes(x = rho, y = rejection_left, color = method)) +
    scale_x_continuous() +
    scale_y_log10(
      labels = scales::label_number(accuracy = 0.001)
    ) +
    ggh4x::facet_grid2(. ~ beta_0, scales = "free_y",
                       labeller = labeller(
                         beta_0 = label_parsed  # Parse LaTeX-like expressions
                       )) +
    geom_point(size = 0.5) +
    geom_line() +
    geom_hline(aes(yintercept = alpha), linetype = "dashed", color = "red") +
    theme_bw() +
    theme(strip.text.x = element_text(margin = margin(0.05, 0, 0.05, 0, "cm"), size = 12),
          strip.text.y = element_text(margin = margin(0, 0.05, 0, 0.05, "cm"), size = 12),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())

  # combine the plot
  plot <- cowplot::plot_grid(plot_gamma, plot_beta, nrow = 2, align = "h")

  # get legend from plot_aux
  legend <- cowplot::get_plot_component(plot_gamma +
                                          theme(legend.position = "bottom",
                                                legend.text = element_text(size = 12)),
                                        'guide-box-bottom', return_all = TRUE)

  #create common x and y labels
  y.grob <- textGrob("Rejection rate", gp = gpar(col = "black", fontsize = 15), rot = 90)
  x.grob <- textGrob(TeX("varying \\rho"), gp = gpar(col = "black", fontsize = 15))

  # compile final plot
  final_plot <- grid.arrange(arrangeGrob(plot, left = y.grob, bottom = x.grob, nrow=1),
                             plot_grid(legend, nrow = 1), nrow = 2, heights=c(9, 1))

  # save the plot
  ggsave(paste0(sprintf('%s/disp-', power_plot_dir), size_parameter,'-power-LEFT.pdf'),
         plot = final_plot, width = 6, height = 5)

  ############################ 3.2 Right-sided test ############################
  # varying gamma
  plot_gamma <- varying_gamma_power |>
    dplyr::filter(rho >= 0) |>
    ggplot(aes(x = rho, y = rejection_right, color = method)) +
    scale_x_continuous() +
    scale_y_log10(
      labels = scales::label_number(accuracy = 0.001)
    ) +
    ggh4x::facet_grid2(. ~ gamma_0, scales = "free_y",
                       labeller = labeller(
                         gamma_0 = label_parsed  # Parse LaTeX-like expressions
                       )) +
    geom_point(size = 0.5) +
    geom_line() +
    geom_hline(aes(yintercept = alpha), linetype = "dashed", color = "red") +
    theme_bw() +
    theme(strip.text.x = element_text(margin = margin(0.05, 0, 0.05, 0, "cm"), size = 12),
          strip.text.y = element_text(margin = margin(0, 0.05, 0, 0.05, "cm"), size = 12),
          legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank())

  # varying beta
  plot_beta <- varying_beta_power |>
    dplyr::filter(rho >= 0) |>
    ggplot(aes(x = rho, y = rejection_right, color = method)) +
    scale_x_continuous() +
    scale_y_log10(
      labels = scales::label_number(accuracy = 0.001)
    ) +
    ggh4x::facet_grid2(. ~ beta_0, scales = "free_y",
                       labeller = labeller(
                         beta_0 = label_parsed  # Parse LaTeX-like expressions
                       )) +
    geom_point(size = 0.5) +
    geom_line() +
    geom_hline(aes(yintercept = alpha), linetype = "dashed", color = "red") +
    theme_bw() +
    theme(strip.text.x = element_text(margin = margin(0.05, 0, 0.05, 0, "cm"), size = 12),
          strip.text.y = element_text(margin = margin(0, 0.05, 0, 0.05, "cm"), size = 12),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())

  # combine the plot
  plot <- cowplot::plot_grid(plot_gamma, plot_beta, nrow = 2, align = "h")

  # get legend from plot_aux
  legend <- cowplot::get_plot_component(plot_gamma +
                                          theme(legend.position = "bottom",
                                                legend.text = element_text(size = 12)),
                                        'guide-box-bottom', return_all = TRUE)

  #create common x and y labels
  y.grob <- textGrob("Rejection rate", gp = gpar(col = "black", fontsize = 15), rot = 90)
  x.grob <- textGrob(TeX("varying \\rho"), gp = gpar(col = "black", fontsize = 15))

  # compile final plot
  final_plot <- grid.arrange(arrangeGrob(plot, left = y.grob, bottom = x.grob, nrow=1),
                             plot_grid(legend, nrow = 1), nrow = 2, heights=c(9, 1))

  # save the plot
  ggsave(paste0(sprintf('%s/disp-', power_plot_dir), size_parameter,'-power-RIGHT.pdf'),
         plot = final_plot, width = 6, height = 5)

  ######################### 4. FDR computation #################################
  # extract all null p-values
  null_df <- current_results |>
    dplyr::filter(rho == 0) |>
    dplyr::select(method, pvalue_list, gamma_0, beta_0, run_id) |>
    tidyr::unnest_longer(pvalue_list) |>
    dplyr::filter(pvalue_list_id %in% c("p.left", "p.right")) |>
    dplyr::mutate(sideness = dplyr::case_when(
      pvalue_list_id == "p.left" ~ "Left-sided",
      pvalue_list_id == "p.right" ~ "Right-sided"
    )) |>
    dplyr::rename(p_value = pvalue_list)

  # plot the rejection rate
  M <- 50
  partition_size <- B / M

  # compute the number of rejections
  reject_df <- null_df |>
    dplyr::group_by(method, gamma_0, beta_0, sideness) |>
    # Calculate the partition group within each group
    dplyr::mutate(partition = ceiling(run_id / partition_size)) |>
    dplyr::group_by(method, sideness, gamma_0, beta_0, partition) |>
    dplyr::summarise(
      n_reject_bonf = sum(p.adjust(p_value, method = "bonferroni") <= FDR, na.rm = TRUE),
      n_reject_bh = sum(p.adjust(p_value, method = "BH") <= FDR, na.rm = TRUE)
    ) |>
    dplyr::ungroup() |>
    # Now average over the 50 partitions within each group
    dplyr::group_by(method, sideness, gamma_0, beta_0) |>
    dplyr::summarise(
      n_reject_bonf = mean(n_reject_bonf, na.rm = TRUE),
      n_reject_bh = mean(n_reject_bh, na.rm = TRUE)
    ) |>
    dplyr::ungroup()

  # Varying gamma_0
  varying_gamma_rejection <- reject_df |>
    rename(BH = n_reject_bh, Bonferroni = n_reject_bonf) |>
    tidyr::pivot_longer(
      cols = c("BH", "Bonferroni"), names_to = "correction", values_to = "rejection"
    ) |>
    dplyr::group_by(method, sideness, correction, beta_0, gamma_0) |>
    dplyr::mutate(rejection = rejection + runif(1, 1e-2, 2e-2)) |>
    dplyr::ungroup() |>
    dplyr::filter(gamma_0 == default_gamma) |>
    ggplot(aes(x = beta_0, y = rejection, color = method)) +
    scale_x_continuous() +
    scale_y_log10() +
    ggh4x::facet_grid2(correction~sideness, scales = "free_y") +
    geom_point(size = 1.5) +
    geom_line() +
    theme_bw() +
    labs(x = TeX("varying $\\gamma_0$ "), y = "Rejection rate",
         title = TeX("Varying $\\gamma_0$")) +
    theme(strip.text.x = element_text(margin = margin(0.05, 0, 0.05, 0, "cm"), size = 12),
          strip.text.y = element_text(margin = margin(0.05, 0.05, 0, 0.05, "cm"), size = 12),
          legend.position = "none",
          axis.title.x = element_text(size = 12),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          panel.grid.minor = element_blank(),
          plot.title = element_blank())

  # Varying beta_0
  varying_beta_rejection <- reject_df |>
    rename(BH = n_reject_bh, Bonferroni = n_reject_bonf) |>
    tidyr::pivot_longer(
      cols = c("BH", "Bonferroni"), names_to = "correction", values_to = "rejection"
    ) |>
    dplyr::group_by(method, sideness, correction, beta_0, gamma_0) |>
    dplyr::mutate(rejection = rejection + runif(1, 1e-2,2e-2)) |>
    dplyr::ungroup() |>
    dplyr::filter(beta_0 == default_beta) |>
    ggplot(aes(x = gamma_0, y = rejection, color = method)) +
    scale_x_continuous() +
    scale_y_log10() +
    ggh4x::facet_grid2(correction~sideness, scales = "free_y") +
    geom_point(size = 1.5) +
    geom_line() +
    theme_bw() +
    labs(x = TeX("varying $\\beta_0$ "), y = "Rejection rate",
         title = TeX("Varying $\\beta_0$")) +
    theme(strip.text.x = element_text(margin = margin(0.05, 0, 0.05, 0, "cm"), size = 12),
          strip.text.y = element_text(margin = margin(0.05, 0.05, 0, 0.05, "cm"), size = 12),
          legend.position = "bottom",
          axis.title.x = element_text(size = 12),
          axis.title.y = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          plot.title = element_blank())

  # Extract the legends from the plots without background
  legend <- get_plot_component(varying_beta_rejection,'guide-box', return_all = TRUE)

  # Create a shared Y-axis title
  y_title <- textGrob("Averaged number of rejections", rot = 90, gp = gpar(fontsize = 15))

  # Combine the two legends into one
  combined_legend <- plot_grid(legend[[3]], ncol = 1)

  # Add the X and Y axis titles to the combined plots
  plots_with_titles <- grid.arrange(
    arrangeGrob(
      plot_grid(
        varying_gamma_rejection,
        varying_beta_rejection + theme(legend.position = "none"),
        ncol = 1,
        rel_heights = c(2, 2) # Adjust plot proportions
      ),
      left = y_title)
  )

  # Combine with the legend at the bottom
  final_plot <- plot_grid(
    plots_with_titles,
    combined_legend,
    ncol = 1,
    rel_heights = c(8, 1) # Adjust as necessary
  )

  # save the plot
  ggsave(paste0(sprintf("%s/disp-", FDR_plot_dir),
                size_parameter, "-multiple_testing_rejection_plot.pdf"),
         final_plot, height = 8, width = 5)

  ################# 5. Scatter plot for dCRT-spaCRT ############################
  # plot the line and dot figure
  dCRT_spaCRT_comparison <- null_df |>
    dplyr::filter(gamma_0 == default_gamma & beta_0 == default_beta) |>
    tidyr::pivot_wider(id_cols =  c("sideness", "run_id"),
                       names_from = "method",
                       values_from = "p_value") |>
    ggplot(aes(x = dCRT, y = spaCRT)) +
    scale_x_continuous(trans = revlog_trans(10), breaks = c(1e-1, 1e-2, 1e-3)) +
    scale_y_continuous(trans = revlog_trans(10), breaks = c(1e-1, 1e-2, 1e-3)) +
    geom_abline(col = "black") +
    facet_wrap(.~sideness) +
    geom_point(size = 0.8) + geom_abline() + my_theme +
    theme(strip.text.x = element_text(margin = margin(0.05, 0, 0.05, 0, "cm"), size = 12),
          strip.text.y = element_text(margin = margin(0.05, 0.05, 0, 0.05, "cm"), size = 12))

  # save the plot
  ggsave(paste0(sprintf('%s/disp-', QQ_plot_dir), size_parameter,'-dCRT-spaCRT.pdf'),
         plot = dCRT_spaCRT_comparison, width = 6, height = 3)

  # full line plot for varying gamma
  dCRT_spaCRT_varying_gamma <- null_df |>
    dplyr::filter(beta_0 == default_beta) |>
    tidyr::pivot_wider(id_cols =  c("gamma_0", "sideness", "run_id"),
                       names_from = "method",
                       values_from = "p_value") |>
    dplyr::mutate(gamma_0 = sprintf("gamma[0] == %d", gamma_0)) |>
    ggplot(aes(x = dCRT, y = spaCRT)) +
    scale_x_continuous(trans = revlog_trans(10), breaks = c(1e-1, 1e-2, 1e-3)) +
    scale_y_continuous(trans = revlog_trans(10), breaks = c(1e-1, 1e-2, 1e-3)) +
    geom_abline(col = "black") +
    ggh4x::facet_grid2(gamma_0~sideness,  labeller = labeller(
      gamma_0 = label_parsed  # Parse LaTeX-like expressions
    )) +
    geom_point(size = 0.8) + my_theme +
    theme(strip.text.x = element_text(margin = margin(0.05, 0, 0.05, 0, "cm"), size = 12),
          strip.text.y = element_text(margin = margin(0.05, 0.05, 0, 0.05, "cm"), size = 12))

  # save the plot
  ggsave(paste0(sprintf('%s/disp-', QQ_plot_dir), size_parameter,'-dCRT-spaCRT-varying-gamma.pdf'),
         plot = dCRT_spaCRT_varying_gamma, width = 4, height = 7)

  # full line plot for varying beta
  dCRT_spaCRT_varying_beta <- null_df |>
    dplyr::filter(gamma_0 == default_gamma) |>
    tidyr::pivot_wider(id_cols =  c("beta_0", "sideness", "run_id"),
                       names_from = "method",
                       values_from = "p_value") |>
    dplyr::mutate(beta_0 = sprintf("beta[0] == %d", beta_0)) |>
    ggplot(aes(x = dCRT, y = spaCRT)) +
    scale_x_continuous(trans = revlog_trans(10), breaks = c(1e-1, 1e-2, 1e-3)) +
    scale_y_continuous(trans = revlog_trans(10), breaks = c(1e-1, 1e-2, 1e-3)) +
    geom_abline(col = "black") +
    ggh4x::facet_grid2(beta_0~sideness,  labeller = labeller(
      beta_0 = label_parsed  # Parse LaTeX-like expressions
    )) +
    geom_point(size = 0.8) + my_theme +
    theme(strip.text.x = element_text(margin = margin(0.05, 0, 0.05, 0, "cm"), size = 12),
          strip.text.y = element_text(margin = margin(0.05, 0.05, 0, 0.05, "cm"), size = 12))

  # save the plot
  ggsave(paste0(sprintf('%s/disp-', QQ_plot_dir), size_parameter,'-dCRT-spaCRT-varying-beta.pdf'),
         plot = dCRT_spaCRT_varying_beta, width = 4, height = 7)
}

############################# 6. Do summary plots ##############################
filtered_results <- joined_results |>
  dplyr::filter(theta == 0.05 & beta_0 == default_beta) |>
  tidyr::unnest_wider(output) |>
  dplyr::select(method, pvalue_list, gamma_0, run_id, rho, computation_time) |>
  tidyr::unnest_longer(pvalue_list) |>
  dplyr::filter(pvalue_list_id %in% c("p.left", "p.right")) |>
  dplyr::mutate(sideness = dplyr::case_when(
    pvalue_list_id == "p.left" ~ "Left-sided",
    pvalue_list_id == "p.right" ~ "Right-sided"
  )) |>
  dplyr::rename(p_value = pvalue_list)

############################## 6.1 Type-I error ################################
p1 <- filtered_results |>
  dplyr::filter(rho == 0) |>
  dplyr::group_by(gamma_0, method, sideness) |>
  dplyr::summarise(
    rejection = ifelse(mean(p_value <= alpha) == 0, 1 / B, mean(p_value <= alpha))
  ) |>
  dplyr::ungroup() |>
  ggplot(aes(x = gamma_0, y = rejection, color = method)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = alpha, linetype = "dashed") +
  scale_y_log10(breaks = c(1e-1, 1e-2, 1e-3, 1e-4),
                labels = c(expression(10^{-1}), expression(10^{-2}), expression(10^{-3}), expression(10^{-4}))) +
  facet_wrap( ~ sideness) +
  labs(x = expression(gamma[0]), y = "Type-I error") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))

############################### 6.2 Power ######################################
p2 <- filtered_results |>
  dplyr::filter(gamma_0 == -3) |>
  dplyr::group_by(gamma_0, method, sideness, rho) |>
  dplyr::summarise(
    rejection = ifelse(mean(p_value <= alpha) == 0, 1 / B, mean(p_value <= alpha))
  ) |>
  dplyr::ungroup() |>
  filter((sideness == "Left-sided" & rho <= 0) | (sideness == "Right-sided" & rho >= 0)) |>
  ggplot(aes(x = rho, y = rejection, color = method)) +
  geom_point() +
  geom_line() +
  facet_wrap( ~ sideness, scales = "free_x") +
  labs(x = expression(rho), y = "Power") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))

############################## 6.3 QQ plots ####################################
p3 <- filtered_results |>
  dplyr::filter(rho == 0 & gamma_0 == -3) |>
  dplyr::mutate(method = factor(method, levels = c("GCM", "score", "dCRT", "spaCRT"),
                                labels = c("GCM test", "Score test", "dCRT", "spaCRT"))) |>
  ggplot(mapping = aes(y = p_value, color = method)) +
  stat_qq_points(ymin = 1e-8, size = 1.5, max_pts_to_plot = 1000) +
  stat_qq_band() +
  geom_abline() +
  facet_wrap(~ sideness, ncol = 2) +
  labs(y = "Observed p-value", x = "Expected p-value") +
  scale_x_continuous(trans = revlog_trans(10), breaks = c(1e-1, 1e-3, 1e-5),
                     labels = c(expression(10^{-1}), expression(10^{-3}), expression(10^{-5}))) +
  scale_y_continuous(trans = revlog_trans(10), breaks = c(1e-1, 1e-3, 1e-5, 1e-7),
                     labels = c(expression(10^{-1}), expression(10^{-3}), expression(10^{-5}), expression(10^{-7}))) +
  theme_bw() +
  theme(legend.position = c(0.65, 0.65),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))

############################ 6.4 Computing time ################################
p4 <- joined_results |>
  dplyr::filter(theta == 0.05) |>
  tidyr::unnest_wider(output) |>
  mutate(method = factor(method,
                         levels = c("GCM", "score", "dCRT", "spaCRT"),
                         labels = c("GCM test", "Score test", "dCRT", "spaCRT"))) |>
  dplyr::group_by(method, gamma_0, rho, beta_0) |>
  dplyr::summarise(mean_computation_time = mean(computation_time)) |>
  dplyr::ungroup() |>
  dplyr::mutate(metric = "Computation time") |>
  ggplot(aes(x = method, fill = method, y = mean_computation_time)) +
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
        plot.title = element_blank(),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))

################################ 6.5 scatter plot ##############################
gamma <- -3
beta <- -5
rho <- 1
theta <- 0.05

# filter the results
results_partial <- joined_results |>
  dplyr::filter(gamma_0 == gamma, beta_0 == beta, rho == !!rho, theta == !!theta) |>
  dplyr::filter(method %in% c("dCRT", "spaCRT")) |>
  dplyr::rowwise() |>
  dplyr::mutate(p.left = output$pvalue_list$p.left, p.right = output$pvalue_list$p.right) |>
  dplyr::select(method, run_id, p.left, p.right)

# give the scatter plot
p5 <- results_partial |>
  tidyr::pivot_longer(cols = c(p.left, p.right),
                      names_to = "side",
                      names_prefix = "p.",
                      values_to = "p") |>
  tidyr::pivot_wider(names_from = method, values_from = c(p)) |>
  filter(side == "right") |>
  filter(run_id %in% sample.int(n = length(unique(run_id)), size = 10000)) |>
  dplyr::mutate(metric = "P-value approximation accuracy") |>
  ggplot(aes(x = dCRT, y = spaCRT)) +
  geom_point(alpha = 0.25) +
  facet_grid(. ~ metric) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  scale_x_continuous(trans = katlabutils::revlog_trans(base = 10),
                     breaks = c(1, 1e-1, 1e-2, 1e-3, 1e-4),
                     labels = c(1, expression(10^{-1}), expression(10^{-2}), expression(10^{-3}), expression(10^{-4}))) +
  scale_y_continuous(trans = katlabutils::revlog_trans(base = 10),
                     breaks = c(1e-1, 1e-3, 1e-5, 1e-7),
                     labels = c(expression(10^{-1}), expression(10^{-3}), expression(10^{-5}), expression(10^{-7}))) +
  labs(x = "p-value (dCRT)", y = "p-value (spaCRT)") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_blank(),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14))

######################## 6.6 Combine plots #####################################
joint_plot <- ((p3 / p1 / p2) | (p5/p4)) +
  plot_layout(ncol = 2, widths = c(1.5, 1)) +
  plot_annotation(tag_levels = "a") +
  theme(legend.position = "none", legend.title = element_blank())

# save the plot
ggsave(filename = sprintf("%s/simulation-summary.pdf", summary_plot_dir),
       plot = joint_plot, height = 8, width = 10)
