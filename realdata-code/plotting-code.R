# This is Rscript aggregating all the results and plotting
library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)
library(dplyr)
library(reshape2)
library(OneR)
library(katlabutils)
library(patchwork)
library(tidyr)
library(tibble)
library(kableExtra)

# Read the command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign the arguments to variables in R
max_cutoff <- as.numeric(args[1])

# plots directory
plots_dir <- "manuscript/figures-and-tables//real-data"

if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
  cat("Directory created:", plots_dir, "\n")
} else {
  cat("Directory already exists:", plots_dir, "\n")
}

# set the data directory
data_dir <- paste0(
  .get_config_path("LOCAL_SPACRT_DATA_DIR"),
  "full_data/summary"
)


################################################################################
# QQ-plot for null simulation
################################################################################
# load the result
combined_result_null <- readRDS(sprintf("%s/combined_result_null.rds",
                                        data_dir))

# compute the rate of GCM default values
combined_result_null |> dplyr::filter(method %in% c("GCM", "spaCRT")) |>
  tidyr::pivot_wider(names_from = method, values_from = `p-value`, id_cols = c("gRNA", "gene", "side")) |>
  dplyr::mutate(GCM_default = ifelse(spaCRT == GCM, TRUE, FALSE)) |>
  dplyr::summarise(GCM_default_rate = mean(GCM_default)) |>
  dplyr::pull()

# overall calibration plot; this code is for plotting Figure 4
condensed_no_stratification <- combined_result_null |>
  ggplot(mapping = aes(y = `p-value`)) +
  stat_qq_points(aes(color = side), ymin = 1e-8, size = 0.8) +
  stat_qq_band() +
  geom_abline(col = "black") +
  scale_x_continuous(trans = revlog_trans(10),
                     breaks = c(1, 1e-2, 1e-4, 1e-6),
                     labels = c(expression(10^{0}), expression(10^{-2}), expression(10^{-4}), expression(10^{-6}))) +
  scale_y_continuous(trans = revlog_trans(10),
                     breaks = c(1, 1e-2, 1e-4, 1e-6, 1e-8),
                     labels = c(expression(10^{0}), expression(10^{-2}), expression(10^{-4}), expression(10^{-6}), expression(10^{-8}))) +
  facet_wrap(~method) +
  guides(color = guide_legend(position = "inside",
                              override.aes = list(size = 1))) +
  labs(y = "Observed p-value",
       x = "Expected p-value") +
  theme_bw() +
  theme(legend.position.inside = c(0.15,0.9),
        legend.background = element_blank(),
        legend.title = element_blank(),
        panel.spacing = unit(0.75, "lines"),
        plot.margin = margin(0.1, 0.5, 0.1, 0.1, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0.5, 0.5, 0.5, 0.5),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text.x = element_text(margin = margin(0.05, 0, 0.05, 0, "cm"), size = 12),
        strip.text.y = element_text(margin = margin(0.05, 0.05, 0, 0.05, "cm"), size = 12))

# non-stratified power plot; this plot is for Figure 5
combined_result_power <- readRDS(sprintf("%s/combined_result_power_individual.rds",
                                         data_dir))

# compute the rate of failure of spaCRT
combined_result_power |> dplyr::filter(method %in% c("GCM", "spaCRT")) |>
  tidyr::pivot_wider(names_from = method, values_from = p_value, id_cols = c("grna_gene")) |>
  dplyr::mutate(GCM_default = ifelse(spaCRT == GCM, TRUE, FALSE)) |>
  dplyr::summarise(
    GCM_default_rate = mean(GCM_default)
  ) |>
  dplyr::pull()

# plot the p-values
p <- combined_result_power |>
  mutate(method = factor(method,
                         levels = c("GCM", "score", "sceptre", "spaCRT"),
                         labels = c("GCM test", "Score test", "sceptre", "spaCRT"))) |>
  ggplot(aes(x = method,
             y = pmax(p_value, 1e-100),
             fill = method)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  labs(x = "Method",
       y = "p-value") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_blank(),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))

######################## combine two plots #####################################
combined_plot <- (condensed_no_stratification | p) +
  plot_layout(ncol = 2, widths = c(1.5, 1)) +
  plot_annotation(tag_levels = "a") +
  theme(legend.position = "none", legend.title = element_blank())

ggsave(sprintf("%s/combined_summary.pdf", plots_dir),
       combined_plot,
       width = 10,
       height = 5)

################################################################################
# Table for rejection with BH and Bonferroni correction
################################################################################
# compute the number of discoveries
alpha <- 0.1
reject_df <- combined_result_null |>
  mutate(p_value = `p-value`,
         side = if_else(side == "left-sided", "left", "right")) |>
  dplyr::select(method, side, p_value) |>
  group_by(method, side) |>
  summarise(n_reject_bonf = sum(p.adjust(p_value, method = "bonferroni") <= alpha, na.rm = TRUE),
            n_reject_bh = sum(p.adjust(p_value, method = "BH") <= alpha, na.rm = TRUE)) |>
  ungroup()

reject_df |>
  tidyr::pivot_wider(names_from = "side",
                     values_from = c("n_reject_bonf", "n_reject_bh")) |>
  mutate(method = as.character(method),
         method = ifelse(method == "score", "Score test", method),
         method = ifelse(method == "GCM", "GCM test", method),
         method = ifelse(method == "sceptre", "dCRT", method),
         method = factor(method, levels = c("GCM test", "Score test", "dCRT", "spaCRT"))) |>
  dplyr::select(method, n_reject_bonf_left, n_reject_bh_left,
                n_reject_bonf_right, n_reject_bh_right) |>
  knitr::kable(
    format = "latex",
    row.names = NA,
    booktabs = TRUE,
    digits = 1,
    caption = "Number of rejections for negative control pairs on the Gasperini data.",
    label = "real_data_rejection",
    col.names = NULL
  ) |>
  kableExtra::kable_styling(latex_options = "hold_position") |>
  add_header_above(c("Method" = 1, "Bonferroni" = 1, "BH" = 1, "Bonferroni" = 1, "BH" = 1)) |>
  add_header_above(c(" " = 1, "Left-sided test" = 2, "Right-sided test" = 2)) |>
  add_header_above(c(" " = 1, "Number of rejections" = 4)) |>
  kableExtra::save_kable(sprintf("%s/real-data-num-rejection.tex", plots_dir))



################################################################################
# Table for timing results
################################################################################

output_dir <- paste0(
  .get_config_path("LOCAL_SPACRT_DATA_DIR"),
  "full_data/"
)
timing_results <- readRDS(paste0(output_dir, "/time_comparison.rds"))

# This code is for creating table 4
timing_results |>
  mutate(method = as.character(method),
         method = ifelse(method == "score", "Score test", method),
         method = ifelse(method == "GCM", "GCM test", method),
         method = factor(method, levels = c("GCM test", "Score test", "dCRT", "spaCRT"))) |>
  rename("Method" = method,
         "Mean" = time_mean,
         "Std dev" = time_sd) |>
  knitr::kable(format = "latex",
               row.names = NA,
               booktabs = TRUE,
               digits = 1,
               caption = "Computation times (in seconds) per test on the Gasperini data.",
               label = "time_comparison") |>
  kableExtra::kable_styling(latex_options = "hold_position") |>
  kableExtra::save_kable(sprintf("%s/gasperini-computing-times.tex", plots_dir))

################################################################################
# Motivating example for introduction
################################################################################

p1 <- combined_result_null |>
  filter(method != "score", side == "left-sided") |>
  mutate(method = factor(method,
                         levels = c("spaCRT", "sceptre", "GCM"),
                         labels = c("spaCRT", "dCRT", "GCM test"))) |>
  ggplot(mapping = aes(y = `p-value`)) +
  stat_qq_points(aes(color = method), ymin = 1e-8, size = 0.8, show.legend = FALSE) +
  stat_qq_band() +
  geom_abline(col = "black") +
  scale_color_manual(values = c("#C77CFF", "#00BFC4", "#F8766D")) +
  scale_x_continuous(trans = revlog_trans(10),
                     breaks = c(1, 1e-2, 1e-4, 1e-6),
                     labels = c(expression(10^{0}), expression(10^{-2}), expression(10^{-4}), expression(10^{-6}))) +
  scale_y_continuous(trans = revlog_trans(10),
                     breaks = c(1, 1e-2, 1e-4, 1e-6, 1e-8),
                     labels = c(expression(10^{0}), expression(10^{-2}), expression(10^{-4}), expression(10^{-6}), expression(10^{-8}))) +
  labs(y = "Observed p-value",
       x = "Expected p-value") +
  theme_bw()

p2 <- timing_results |>
  filter(method != "score") |>
  mutate(method = factor(method,
                         levels = c("GCM", "dCRT", "spaCRT"),
                         labels = c("GCM test", "dCRT", "spaCRT"))) |>
  ggplot(aes(x = method, y = time_mean, fill = method)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("#F8766D", "#00BFC4", "#C77CFF")) +
  labs(x = "Method",
       y = "Mean seconds per test") +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# This code is for plotting Figure 1
p <- p1 + plot_spacer() + p2 + plot_layout(ncol = 3, guides = "collect", widths = c(2.5, 0.25, 1)) &
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.margin = unit(c(0.05, 0.05, 0, 0.05), "cm"))

ggsave(sprintf("%s/motivating_example.pdf", plots_dir),
       plot = p,
       height = 3,
       width = 4.25)

################################################################################
# The following code is used to create Figure 18 and 19 in Appendix
################################################################################

# specify the theme
my_theme <- theme_bw() + theme(axis.line = element_line(color = "black"),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.border = element_blank(),
                               panel.background = element_blank(),
                               plot.title = element_text(hjust = 0.5, size=11))

# specify colors
pannel_cols <- c("darkorchid4", "dodgerblue4", "dodgerblue2", "dodgerblue1", "deepskyblue")

# plot against effective sample size under truncated effective sample size
filtered_result <- combined_result_null |>
  dplyr::filter(n_nonzero_trt <= max_cutoff) |>
  dplyr::mutate(n_nonzero_trt_bin = cut_number(n_nonzero_trt, n = 5))

# plot for each method
n_to_sample <- filtered_result |>
  group_by(n_nonzero_trt_bin, method, side) |>
  summarize(count = n()) |>
  pull(count) |> min()

set.seed(1)
to_plot <- filtered_result |>
  group_by(n_nonzero_trt_bin, method, side) |>
  sample_n(n_to_sample)

# plot
text_added_data <- to_plot |>
  dplyr::arrange(n_nonzero_trt_bin) |>
  dplyr::mutate(text_label = dplyr::if_else(method == "GCM" & side == "left-sided",
                                            "Effective sample size", NA))

plot_GCM <- text_added_data |>
  dplyr::filter(method == "GCM") |>
  ggplot(mapping = aes(y = `p-value` ,
                       col = n_nonzero_trt_bin
  )) +
  stat_qq_points(ymin = 1e-8, size = 0.8) +
  stat_qq_band() +
  ggh4x::facet_grid2(method ~ side, scales = "free", independent = "x") +
  scale_x_continuous(trans = revlog_trans(10)) +
  scale_y_continuous(trans = revlog_trans(10), breaks = c(1, 1e-2, 1e-4, 1e-6)) +
  geom_abline(col = "black") +
  my_theme +
  theme(strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")),
        strip.text.y = element_text(margin = margin(0,0.1,0,0.1, "cm")),
        legend.position = "none",
        legend.margin=margin(t = -0.2, r = -0.2, unit = 'cm'),
        legend.text = element_text(size = 8),  # Change legend text size
        legend.title = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) +
  guides(color = guide_legend(
    keywidth = 0.0,
    keyheight = 0.15,
    default.unit = "inch",
    override.aes = list(size = 2.5))) +
  scale_color_manual(values = pannel_cols, name = "Dispersion")

plot_score_glm_nb <- text_added_data |>
  dplyr::filter(method == "score") |>
  ggplot(mapping = aes(y = `p-value` ,
                       col = n_nonzero_trt_bin
  )) +
  stat_qq_points(ymin = 1e-8, size = 0.8) +
  stat_qq_band() +
  ggh4x::facet_grid2(method ~ side, scales = "free", independent = "x") +
  scale_x_continuous(trans = revlog_trans(10)) +
  scale_y_continuous(trans = revlog_trans(10), breaks = c(1, 1e-2, 1e-4, 1e-6)) +
  geom_abline(col = "black") +
  my_theme +
  theme(strip.text.y = element_text(margin = margin(0,0.1,0,0.1, "cm")),
        legend.position = "none",
        legend.margin=margin(t = -0.2, r = -0.2, unit = 'cm'),
        legend.text = element_text(size = 8),  # Change legend text size
        legend.title = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text.x = element_blank()) +
  guides(color = guide_legend(
    keywidth = 0.0,
    keyheight = 0.15,
    default.unit = "inch",
    override.aes = list(size = 2.5))) +
  scale_color_manual(values = pannel_cols, name = "Dispersion")

plot_spaCRT <- text_added_data |>
  dplyr::filter(method == "spaCRT") |>
  ggplot(mapping = aes(y = `p-value` ,
                       col = n_nonzero_trt_bin
  )) +
  stat_qq_points(ymin = 1e-8, size = 0.8) +
  stat_qq_band() +
  ggh4x::facet_grid2(method ~ side, scales = "free", independent = "x") +
  scale_x_continuous(trans = revlog_trans(10)) +
  scale_y_continuous(trans = revlog_trans(10), breaks = c(1, 1e-2, 1e-4, 1e-6)) +
  geom_abline(col = "black") +
  my_theme +
  theme(strip.text.y = element_text(margin = margin(0,0.1,0,0.1, "cm")),
        legend.position = "none",
        legend.margin=margin(t = -0.2, r = -0.2, unit = 'cm'),
        legend.text = element_text(size = 8),  # Change legend text size
        legend.title = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text.x = element_blank()) +
  guides(color = guide_legend(
    keywidth = 0.0,
    keyheight = 0.15,
    default.unit = "inch",
    override.aes = list(size = 2.5))) +
  scale_color_manual(values = pannel_cols, name = "Dispersion")

plot_sceptre <- text_added_data |>
  dplyr::filter(method == "sceptre") |>
  ggplot(mapping = aes(y = `p-value` ,
                       col = n_nonzero_trt_bin
  )) +
  stat_qq_points(ymin = 1e-8, size = 0.8) +
  stat_qq_band() +
  ggh4x::facet_grid2(method ~ side, scales = "free", independent = "x") +
  scale_x_continuous(trans = revlog_trans(10), breaks = c(1, 1e-2, 1e-4)) +
  scale_y_continuous(trans = revlog_trans(10), breaks = c(1, 1e-2, 1e-4, 1e-6)) +
  geom_abline(col = "black") +
  my_theme +
  theme(strip.text.y = element_text(margin = margin(0,0.1,0,0.1, "cm")),
        legend.position = "none",
        legend.margin=margin(t = -0.2, r = -0.2, unit = 'cm'),
        legend.text = element_text(size = 8),  # Change legend text size
        legend.title = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text.x = element_blank(),
        strip.background.x = element_blank(),
        strip.placement.x = "outside") +
  guides(color = guide_legend(
    keywidth = 0.0,
    keyheight = 0.15,
    default.unit = "inch",
    override.aes = list(size = 2.5))) +
  scale_color_manual(values = pannel_cols, name = "Dispersion")

plot_aux <- text_added_data |>
  dplyr::filter(method == "sceptre") |>
  ggplot(mapping = aes(y = `p-value` ,
                       col = n_nonzero_trt_bin
  )) +
  stat_qq_points(ymin = 1e-8, size = 0.8) +
  stat_qq_band() +
  ggh4x::facet_grid2(method ~ side, scales = "free", independent = "x") +
  scale_x_continuous(trans = revlog_trans(10)) +
  scale_y_continuous(trans = revlog_trans(10)) +
  geom_abline(col = "black") +
  my_theme +
  theme(legend.position = "bottom",
        legend.margin=margin(t = -0.5, r = -0.5, unit = 'cm'),
        legend.text = element_text(size = 10),  # Change legend text size
        legend.title = element_text(size = 10)) +
  guides(color = guide_legend(
    keywidth = 0.0,
    keyheight = 0.15,
    default.unit = "inch",
    override.aes = list(size = 2.5))) +
  scale_color_manual(values = pannel_cols, name = "Effective sample size")

# get legend from plot_aux
legend <- get_plot_component(plot_aux, 'guide-box-bottom', return_all = TRUE)

# combine the plot
plot <- plot_grid(plot_GCM, plot_score_glm_nb, plot_spaCRT, plot_sceptre,
                  nrow = 4, align = "h")

#create common x and y labels
y.grob <- textGrob("Observed p-value",
                   gp=gpar(col="black"), rot=90)

x.grob <- textGrob("Expected null p-value",
                   gp=gpar(col="black"))


#add to plot; this plot is for creating Figure 18
final_plot <- grid.arrange(arrangeGrob(plot, left = y.grob, bottom = x.grob, nrow=1),
                           plot_grid(legend, nrow = 1), nrow=2, heights=c(10, 1))

# save the plot
ggsave(sprintf("%s/facet_plot_different_withglmnb_%d.pdf",
               plots_dir,
               max_cutoff),
       plot = final_plot,
       height = 9,
       width = 6)

# plot against dispersion parameter
combined_result_null <- combined_result_null |>
  dplyr::mutate(dispersion_group = cut_number(round(as.numeric(dispersion), 2),
                                              n = 5))

n_to_sample <- combined_result_null |>
  group_by(dispersion_group, method, side) |>
  summarize(count = n()) |>
  pull(count) |> min()

set.seed(1)
to_plot <- combined_result_null |>
  group_by(dispersion_group, method, side) |>
  sample_n(n_to_sample)

# add effective sample size text label
text_added_data <- to_plot |>
  dplyr::arrange(dispersion_group) |>
  dplyr::mutate(text_label = dplyr::if_else(method == "GCM" & side == "left-sided",
                                            "Estimated dispersion", NA))

plot_GCM <- text_added_data |>
  dplyr::filter(method == "GCM") |>
  ggplot(mapping = aes(y = `p-value` ,
                       col = dispersion_group
  )) +
  stat_qq_points(ymin = 1e-8, size = 0.8) +
  stat_qq_band() +
  ggh4x::facet_grid2(method ~ side, scales = "free", independent = "x") +
  scale_x_continuous(trans = revlog_trans(10)) +
  scale_y_continuous(trans = revlog_trans(10), breaks = c(1, 1e-2, 1e-4, 1e-6)) +
  geom_abline(col = "black") +
  my_theme +
  theme(strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")),
        strip.text.y = element_text(margin = margin(0,0.1,0,0.1, "cm")),
        legend.position = "none",
        legend.margin=margin(t = -0.2, r = -0.2, unit = 'cm'),
        legend.text = element_text(size = 8),  # Change legend text size
        legend.title = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank()) +
  guides(color = guide_legend(
    keywidth = 0.0,
    keyheight = 0.15,
    default.unit = "inch",
    override.aes = list(size = 2.5))) +
  scale_color_manual(values = pannel_cols, name = "Size parameter")

plot_score_glm_nb <- text_added_data |>
  dplyr::filter(method == "score") |>
  ggplot(mapping = aes(y = `p-value` ,
                       col = dispersion_group
  )) +
  stat_qq_points(ymin = 1e-8, size = 0.8) +
  stat_qq_band() +
  ggh4x::facet_grid2(method ~ side, scales = "free", independent = "x") +
  scale_x_continuous(trans = revlog_trans(10)) +
  scale_y_continuous(trans = revlog_trans(10), breaks = c(1, 1e-2, 1e-4, 1e-6)) +
  geom_abline(col = "black") +
  my_theme +
  theme(strip.text.y = element_text(margin = margin(0,0.1,0,0.1, "cm")),
        legend.position = "none",
        legend.margin=margin(t = -0.2, r = -0.2, unit = 'cm'),
        legend.text = element_text(size = 8),  # Change legend text size
        legend.title = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text.x = element_blank()) +
  guides(color = guide_legend(
    keywidth = 0.0,
    keyheight = 0.15,
    default.unit = "inch",
    override.aes = list(size = 2.5))) +
  scale_color_manual(values = pannel_cols, name = "Size parameter")

plot_spaCRT <- text_added_data |>
  dplyr::filter(method == "spaCRT") |>
  ggplot(mapping = aes(y = `p-value` ,
                       col = dispersion_group
  )) +
  stat_qq_points(ymin = 1e-8, size = 0.8) +
  stat_qq_band() +
  ggh4x::facet_grid2(method ~ side, scales = "free", independent = "x") +
  scale_x_continuous(trans = revlog_trans(10)) +
  scale_y_continuous(trans = revlog_trans(10), breaks = c(1, 1e-2, 1e-4, 1e-6)) +
  geom_abline(col = "black") +
  my_theme +
  theme(strip.text.y = element_text(margin = margin(0,0.1,0,0.1, "cm")),
        legend.position = "none",
        legend.margin=margin(t = -0.2, r = -0.2, unit = 'cm'),
        legend.text = element_text(size = 8),  # Change legend text size
        legend.title = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text.x = element_blank()) +
  guides(color = guide_legend(
    keywidth = 0.0,
    keyheight = 0.15,
    default.unit = "inch",
    override.aes = list(size = 2.5))) +
  scale_color_manual(values = pannel_cols, name = "Size parameter")

plot_sceptre <- text_added_data |>
  dplyr::filter(method == "sceptre") |>
  ggplot(mapping = aes(y = `p-value` ,
                       col = dispersion_group
  )) +
  stat_qq_points(ymin = 1e-8, size = 0.8) +
  stat_qq_band() +
  ggh4x::facet_grid2(method ~ side, scales = "free", independent = "x") +
  scale_x_continuous(trans = revlog_trans(10), breaks = c(1, 1e-2, 1e-4)) +
  scale_y_continuous(trans = revlog_trans(10), breaks = c(1, 1e-2, 1e-4, 1e-6)) +
  geom_abline(col = "black") +
  my_theme +
  theme(strip.text.y = element_text(margin = margin(0,0.1,0,0.1, "cm")),
        legend.position = "none",
        legend.margin=margin(t = -0.2, r = -0.2, unit = 'cm'),
        legend.text = element_text(size = 8),  # Change legend text size
        legend.title = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text.x = element_blank(),
        strip.background.x = element_blank(),
        strip.placement.x = "outside") +
  guides(color = guide_legend(
    keywidth = 0.0,
    keyheight = 0.15,
    default.unit = "inch",
    override.aes = list(size = 2.5))) +
  scale_color_manual(values = pannel_cols, name = "Size parameter")

plot_aux <- text_added_data |>
  dplyr::filter(method == "sceptre") |>
  ggplot(mapping = aes(y = `p-value` ,
                       col = dispersion_group
  )) +
  stat_qq_points(ymin = 1e-8, size = 0.8) +
  stat_qq_band() +
  ggh4x::facet_grid2(method ~ side, scales = "free", independent = "x") +
  scale_x_continuous(trans = revlog_trans(10)) +
  scale_y_continuous(trans = revlog_trans(10)) +
  geom_abline(col = "black") +
  my_theme +
  theme(legend.position = "bottom",
        legend.margin=margin(t = -0.5, r = -0.5, unit = 'cm'),
        legend.text = element_text(size = 10),  # Change legend text size
        legend.title = element_text(size = 10)) +
  guides(color = guide_legend(
    keywidth = 0.0,
    keyheight = 0.15,
    default.unit = "inch",
    override.aes = list(size = 2.5))) +
  scale_color_manual(values = pannel_cols, name = "Size parameter")

# get legend from plot_aux
legend <- get_plot_component(plot_aux, 'guide-box-bottom', return_all = TRUE)

# combine the plot; this code is for plotting Figure 19
plot <- plot_grid(plot_GCM, plot_score_glm_nb, plot_spaCRT, plot_sceptre,
                  nrow = 4, align = "h")

#create common x and y labels
y.grob <- textGrob("Observed p-value",
                   gp=gpar(col="black"), rot=90)

x.grob <- textGrob("Expected null p-value",
                   gp=gpar(col="black"))


# add to plot
final_plot <- grid.arrange(arrangeGrob(plot, left = y.grob, bottom = x.grob, nrow=1),
                           plot_grid(legend, nrow = 1), nrow=2, heights=c(10, 1))



# save the plot
ggsave(sprintf("%s/facet_plot_different_withglmnb_dispersion.pdf", plots_dir),
       plot = final_plot,
       height = 9,
       width = 6)

