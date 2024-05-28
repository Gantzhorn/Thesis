# Title: Main file for the thesis
# Author: Anders Gantzhorn Kristensen (University of Copenhagen, andersgantzhorn@gmail.com)
# Date: 2024-01-31 (Last Updated: 2024-05-28)
#-----------------------------------------------------------------------------------------------------------------------------#
# Project: Tipping Point Estimation in Ecological Systems using Stochastic Differential Equations
# Description: This script builds an runs all the experiments and results shown in the thesis.
#-----------------------------------------------------------------------------------------------------------------------------#
# License: MIT License (for more information, see LICENSE file in the repository).
# Dependencies: tidyverse and requirements from sourced files.
#-----------------------------------------------------------------------------------------------------------------------------#

library(tidyverse)
source("source_code/tipping_simulations.R")
source("source_code/model_fitting.R")
source("source_code/model_diagnostics.R")

thesis_theme <- ggthemes::theme_base() + 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.background = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank()
  )

ggplot2::theme_set(thesis_theme)

# Colorblind friendly plot palette
thesis_palette <- c("#E69F00", "#56B4E9", "#009E73", "#CCB400", "#0072B2",
                    "#D55E00", "#CC79A7", "#F0E442", "#D55E87", "#6E016B")

# Auxilliary methods used in the main-file
colwise_median <- function(estimates) {
  apply(estimates, 2, function(x) quantile(x, probs = 0.5))
}

colwise_quantile <- function(estimates, probs = 0.975) {
  apply(estimates, 2, function(x) quantile(x, probs = probs))
}


###------------------------------------------------------------------------###
###----------------------------Method section------------------------------###
###------------------------------------------------------------------------###
# Create plot of path for figures in "Tipping in stochastic systems"
# Choose parameters appropriate for any of the diffusions
true_param <- c(1.5, 0.4, -0.2, 0.1)
actual_dt <- 1/12 * 1 / 10
tau <- 50
t_0 <- 25


simulate_multiple_times <- function(simulate_function, n_times, ...) {
  simulate_res <- do.call(rbind,
          lapply(1:n_times, function(i) {
            simulate_function(..., seed = methodPlotseed + i) |>
              mutate(sample_id = i)
          })
  )
  simulate_res$sample_id <- factor(simulate_res$sample_id)
  simulate_res
}

numSim <- 3
methodPlotseed <- 210424 
xs_all <- bind_rows(
  simulate_multiple_times(simulate_additive_noise_tipping_model, numSim, actual_dt, true_param, tau, t_0)  |> 
    mutate(Model = "Additive"),
  simulate_multiple_times(simulate_squareroot_noise_tipping_model, numSim, actual_dt, true_param, tau, t_0) |>
    mutate(Model = "Square root"),
  simulate_multiple_times(simulate_linear_noise_tipping_model, numSim, actual_dt, true_param, tau, t_0) |>
    mutate(Model = "Linear"),
  simulate_multiple_times(simulate_t_distribution_tipping_model, numSim, actual_dt, true_param, tau, t_0) |>
    mutate(Model = "t-distribution"),
  simulate_multiple_times(simulate_F_distribution_tipping_model, numSim, actual_dt, true_param, tau, t_0) |>
    mutate(Model = "F-distribution"),
  simulate_multiple_times(simulate_jacobi_diffusion_tipping_model, numSim, actual_dt, true_param, tau, t_0) |>
    mutate(Model = "Jacobi diffusion")
) |> 
  mutate(Model = factor(Model,
                        levels = c("Additive", "Square root", "Linear",
                                   "t-distribution", "F-distribution", "Jacobi diffusion")))

fixed_point_lines <- tibble(
  t = seq(min(xs_all$t), max(xs_all$t), length.out = 100),
  upper = true_param[2] + sqrt(-true_param[3] * (1 - pmax(0, t - t_0) / tau) / true_param[1]),
  lower = true_param[2] - sqrt(-true_param[3] * (1 - pmax(0, t - t_0) / tau) / true_param[1]),
)

sample_paths_plot_small_scale <- ggplot(xs_all, aes(x = t, y = X_t, color = Model)) +
  geom_line(data = fixed_point_lines, aes(x = t, y = upper), linewidth = 1, linetype = "solid", color = "grey25") +
  geom_line(data = fixed_point_lines, aes(x = t, y = lower), linewidth = 0.85, linetype = "dashed", color = "grey25") +
  geom_step(linewidth = 0.5) + 
  geom_hline(yintercept = true_param[2], linetype = "dashed") +
  facet_grid(sample_id ~ Model) +
  ylim(0, 1) +
  scale_color_manual(values = thesis_palette) +
  labs(y = expression(X[t])) + 
  guides(color = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) + 
  theme(legend.text=element_text(size=20),
        legend.title = element_text(size = 22),
        axis.title = element_text(size = 22),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        legend.key.width = unit(2, "lines"),
        legend.key.height = unit(1, "lines"),
        legend.position = "bottom")

# ggsave(sample_paths_plot_small_scale, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "sample_paths_plot_small_scale.jpeg",
#        height = 8, width = 13, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)

# For the appendix create the same plot but reversed
true_param_reversed <- c(-1.5, 0.4, 0.2, 0.1)

methodPlotseed <- 210424 
xs_all_reversed <- bind_rows(
  simulate_multiple_times(simulate_additive_noise_tipping_model, numSim, actual_dt, true_param_reversed, tau, t_0)  |> 
    mutate(Model = "Additive"),
  simulate_multiple_times(simulate_squareroot_noise_tipping_model, numSim, actual_dt, true_param_reversed, tau, t_0) |>
    mutate(Model = "Square root"),
  simulate_multiple_times(simulate_linear_noise_tipping_model, numSim, actual_dt, true_param_reversed, tau, t_0) |>
    mutate(Model = "Linear"),
  simulate_multiple_times(simulate_t_distribution_tipping_model, numSim, actual_dt, true_param_reversed, tau, t_0) |>
    mutate(Model = "t-distribution"),
  simulate_multiple_times(simulate_F_distribution_tipping_model, numSim, actual_dt, true_param_reversed, tau, t_0) |>
    mutate(Model = "F-distribution"),
  simulate_multiple_times(simulate_jacobi_diffusion_tipping_model, numSim, actual_dt, true_param_reversed, tau, t_0) |>
    mutate(Model = "Jacobi diffusion")
) |> 
  mutate(Model = factor(Model,
                        levels = c("Additive", "Square root", "Linear",
                                   "t-distribution", "F-distribution", "Jacobi diffusion")))

sample_paths_plot_small_scale_reversed <- ggplot(xs_all_reversed, aes(x = t, y = X_t, color = Model)) +
  geom_line(data = fixed_point_lines, aes(x = t, y = upper), linewidth = 1, linetype = "dashed", color = "grey25") +
  geom_line(data = fixed_point_lines, aes(x = t, y = lower), linewidth = 0.85, linetype = "solid", color = "grey25") +
  geom_step(linewidth = 0.5) + 
  geom_hline(yintercept = true_param[2], linetype = "dashed") +
  facet_grid(sample_id ~ Model) +
  ylim(-0.2, 0.8) +
  scale_color_manual(values = thesis_palette) +
  labs(y = expression(X[t])) + 
  guides(color = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) + 
  theme(legend.text=element_text(size=20),
        legend.title = element_text(size = 22),
        axis.title = element_text(size = 22),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        legend.key.width = unit(2, "lines"),
        legend.key.height = unit(1, "lines"),
        legend.position = "bottom")

# ggsave(sample_paths_plot_small_scale_reversed, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "sample_paths_plot_small_scale_reversed.jpeg",
#        height = 8, width = 13, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)

# Create plot of path for figures in "Tipping in stochastic systems"
# Choose parameters appropriate for all but jacobi diffusion
true_param <- c(1.5, 3, -1.75, 0.1)

actual_dt <- 1/12 * 1 / 10
tau <- 50
t_0 <- 25

numSim <- 2
methodPlotseed <- 210424 
xs_all <- bind_rows(
  simulate_multiple_times(simulate_additive_noise_tipping_model, numSim, actual_dt, true_param, tau, t_0)  |> 
    mutate(Model = "Additive"),
  simulate_multiple_times(simulate_squareroot_noise_tipping_model, numSim, actual_dt, true_param, tau, t_0) |>
    mutate(Model = "Square root"),
  simulate_multiple_times(simulate_linear_noise_tipping_model, numSim, actual_dt, true_param, tau, t_0) |>
    mutate(Model = "Linear"),
  simulate_multiple_times(simulate_t_distribution_tipping_model, numSim, actual_dt, true_param, tau, t_0) |>
    mutate(Model = "t-distribution"),
  simulate_multiple_times(simulate_F_distribution_tipping_model, numSim, actual_dt, true_param, tau, t_0) |>
    mutate(Model = "F-distribution")
)|> 
  mutate(Model = factor(Model,
                        levels = c("Additive", "Square root", "Linear",
                                   "t-distribution", "F-distribution")))

fixed_point_lines_big <- tibble(
  t = seq(min(xs_all$t), max(xs_all$t), length.out = 100),
  upper = true_param[2] + sqrt(-true_param[3] * (1 - pmax(0, t - t_0) / tau) / true_param[1]),
  lower = true_param[2] - sqrt(-true_param[3] * (1 - pmax(0, t - t_0) / tau) / true_param[1]),
)


sample_paths_plot_big_scale <- ggplot(xs_all, aes(x = t, y = X_t, color = Model)) +
  geom_line(data = fixed_point_lines_big, aes(x = t, y = upper), linewidth = 1, linetype = "solid", color = "grey25") +
  geom_line(data = fixed_point_lines_big, aes(x = t, y = lower), linewidth = 0.75, linetype = "dashed", color = "grey25") +
  geom_step(linewidth = 0.5) + 
  geom_hline(yintercept = true_param[2], linetype = "dashed") +
  facet_grid(sample_id ~ Model) +
  scale_color_manual(values = thesis_palette) +
  labs(y = expression(X[t])) + 
  theme(legend.text=element_text(size=20),
        legend.title = element_text(size = 22),
        axis.title = element_text(size = 22),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        legend.key.width = unit(2, "lines"),
        legend.key.height = unit(1, "lines"),
        legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) + 
  ylim(min(fixed_point_lines_big$lower)-0.1, max(fixed_point_lines_big$upper)+1)

# ggsave(sample_paths_plot_big_scale, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "sample_paths_plot_big_scale.jpeg",
#        height = 8, width = 13, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)

# For the appendix create the same plot but reversed
true_param_reversed <- c(-1.5, 3, 1.75, 0.1)

methodPlotseed <- 210424 
xs_all_reversed <- bind_rows(
  simulate_multiple_times(simulate_additive_noise_tipping_model, numSim, actual_dt, true_param_reversed, tau, t_0)  |> 
    mutate(Model = "Additive"),
  simulate_multiple_times(simulate_squareroot_noise_tipping_model, numSim, actual_dt, true_param_reversed, tau, t_0) |>
    mutate(Model = "Square root"),
  simulate_multiple_times(simulate_linear_noise_tipping_model, numSim, actual_dt, true_param_reversed, tau, t_0) |>
    mutate(Model = "Linear"),
  simulate_multiple_times(simulate_t_distribution_tipping_model, numSim, actual_dt, true_param_reversed, tau, t_0) |>
    mutate(Model = "t-distribution"),
  simulate_multiple_times(simulate_F_distribution_tipping_model, numSim, actual_dt, true_param_reversed, tau, t_0) |>
    mutate(Model = "F-distribution")
) |> 
  mutate(Model = factor(Model,
                        levels = c("Additive", "Square root", "Linear",
                                   "t-distribution", "F-distribution")))

sample_paths_plot_big_scale_reversed <- ggplot(xs_all_reversed, aes(x = t, y = X_t, color = Model)) +
  geom_line(data = fixed_point_lines_big, aes(x = t, y = upper), linewidth = 1, linetype = "dashed", color = "grey25") +
  geom_line(data = fixed_point_lines_big, aes(x = t, y = lower), linewidth = 0.85, linetype = "solid", color = "grey25") +
  geom_step(linewidth = 0.5) + 
  geom_hline(yintercept = true_param[2], linetype = "dashed") +
  facet_grid(sample_id ~ Model) +
  ylim(1.5, 5) +
  scale_color_manual(values = thesis_palette) +
  labs(y = expression(X[t])) + 
  guides(color = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) + 
  theme(legend.text=element_text(size=20),
        legend.title = element_text(size = 22),
        axis.title = element_text(size = 22),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        legend.key.width = unit(2, "lines"),
        legend.key.height = unit(1, "lines"),
        legend.position = "bottom")

# ggsave(sample_paths_plot_big_scale_reversed, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "sample_paths_plot_big_scale_reversed.jpeg",
#        height = 8, width = 13, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)


# Potential, flow and bifurcation plots for "Deterministic Dynamical Systems and bifurcations
lambda_double_well_vec <- c(0, 0.25, 2 / (3 * sqrt(3)))

double_well <- function(x, lambda){
  x^4 / 4 - x^2 / 2 + lambda * x
}

double_well_deriv <- function(x, lambda){
  -x^3 + x - lambda
}

double_well_doublederiv <- function(x, lambda){
  -3 * x^2 + 1
}

# Generate x values
x_values_vec <- seq(from = -1.2, to = 1.2, by = 0.05)

# Create an empty data frame to store results
double_well_plot_data <- tibble(x_values = rep(x_values_vec, length(lambda_double_well_vec)),
                                x_prime = NA,
                                x = NA,
                                force = NA,
                                lambda = factor(rep(lambda_double_well_vec, each = length(x_values_vec))))



# Loop over lambda parameters and calculate force values
double_well_plot_data <- lapply(lambda_double_well_vec, function(lambda) {
  force_values <- double_well_deriv(x_values_vec, lambda)
  tibble(
    x_values = x_values_vec,
    x_prime = x_values_vec,
    force = force_values,
    lambda = factor(lambda, levels = lambda_double_well_vec)
  )
}) |> bind_rows()

# Find roots for each lambda parameter
roots_neglambda <- purrr::map2(lambda_double_well_vec, lambda_double_well_vec, ~ {
  roots <- polyroot(c(-.y, 1, 0, -1))
  real_roots <- roots[dplyr::near(Im(roots), 0)]
  data.frame(lambda = .x, root = sort(Re(real_roots)))  # Enrich each set of roots with lambda information
})

# Combine the list into a single data frame
roots_data_neglambda <- do.call(rbind, roots_neglambda) 

roots_data_neglambda <- roots_data_neglambda |>
  mutate(
    root_color = case_when(
      double_well_doublederiv(root) < 0 ~ "black",
      TRUE ~ "white"
    ),
    lambda = factor(lambda)
  ) |>
  as_tibble()

# Make half circle
circleFun <- function(center=c(0,0), diameter=1, npoints=100, start=0, end=2)
{
  tt <- seq(start*pi, end*pi, length.out=npoints)
  tibble(x = center[1] + diameter / 2 * cos(tt),
         y = center[2] + diameter / 2 * sin(tt))
}
dat_half_circle <- circleFun(c(0.577, 0), 0.1, start = 1/2, end = 3/2)
dat_half_circle <- rbind(dat_half_circle, dat_half_circle[1,])
dat_half_circle$lambda <- factor(lambda_double_well_vec[3]) 


double_well_plot_neg <- ggplot(double_well_plot_data, aes(x = x_prime, y = force)) +
  geom_line(aes(col = lambda), linewidth = 1.25) +
  geom_hline(yintercept = 0) +
  geom_point(data = roots_data_neglambda, aes(x = root, y = 0, fill = root_color), size = 3, pch = 21) +
  scale_fill_manual(values = c("black", "white"), guide = "none") + 
  geom_polygon(data = dat_half_circle, aes(x, y), inherit.aes = FALSE, fill = "white", color = "black") +
  labs(x = "x", y = expression(dot(x))) +
  facet_grid(. ~ lambda) + 
  scale_color_manual(values = thesis_palette, labels = c("0", "0.25", expression(lambda * phantom()[" c"]))) +
  theme(strip.text = element_blank(),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 18)) +
  guides(color = guide_legend(title = expression(lambda), override.aes = list(linewidth = 5))) + 
  coord_fixed()



# ggsave(double_well_plot_neg, path = paste0(getwd(), "/tex_files/figures"), filename = "double_well_plot_neg.jpeg",
#        height = 8, width = 10, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)

# Plot the potential

potential_plot_data <- tibble(
  xs = rep(seq(-1.5, 1.5, length.out = 100), length(lambda_double_well_vec)),
  lambda = rep(lambda_double_well_vec, each = 100),
  double_well = double_well(xs, lambda)
)

# Plotting code
potential_plot_data |>
  mutate(lambda = factor(lambda)) |> 
  ggplot(aes(x = xs, y = double_well, col = lambda)) +
  geom_line(linewidth = 1.15) +
  scale_color_manual(values = thesis_palette, labels = c("0", "0.25", expression(lambda[c]))) +
  guides(color = guide_legend(
    title = expression(lambda)
  )) +
  labs(
    y = expression(V(x, lambda)),
    col = expression(lambda)
  ) + theme(ggthemes::theme_base())

# Plot the normal form
normal_form <- function(x, A, lambda, m){
  -(A*(x-m)^2 + lambda)
}

normal_form_deriv <- function(x, A, lambda, m){
  2 * A * (x - m)
}

A <- 0.8
lambda <- c(-2.5, 0, 2.5)
m <- -1.5
xs <- seq(-2.5 + m, m + 2.5, length.out = 100)



normal_form_plot_data <- tibble(x = rep(xs, times = 3), lambdas = rep(lambda, 100),
                                ys = normal_form(x, A, lambdas, m)) |> 
  mutate(lambda = factor(lambdas))

roots_normal_form <- tibble(roots = c(m + c(1,-1) * sqrt(-lambda[1] / A), m + sqrt(-lambda[2] / A)), lambda = 
                              c(rep("-2.5",2), "0"), sign =  normal_form_deriv(roots, A, lambda = lambda, m)) |> 
  mutate(lambda = factor(lambda), color = ifelse(sign<0, "black", "white")) 

dat_half_circle1 <- circleFun(c(m, 0), 0.225, start = 1/2, end = 3/2)
dat_half_circle1 <- rbind(dat_half_circle1, dat_half_circle1[1,])
dat_half_circle1$lambda <- factor(0) 


normal_form_plot_data |> 
  ggplot(aes(x = x, y = ys, col = lambda)) + 
  geom_line(linewidth = 1.25) + 
  geom_hline(yintercept = 0) + 
  scale_color_manual(values = thesis_palette) +
  geom_point(data = roots_normal_form, aes(x = roots, y = 0, fill = color), size = 3, pch = 21, color = "black") +
  scale_fill_manual(values = c("black", "white"), guide = "none") + 
  geom_polygon(data = dat_half_circle1, aes(x, y), inherit.aes = FALSE, fill = "black", color = "black") +
  facet_wrap(~lambda) + 
  labs(x = "x", y = expression(dot(x))) + 
  theme(strip.text = element_blank(),
        legend.text = element_text(size = 14),
        panel.spacing = unit(2, "lines"),
        axis.text.y = element_text(size = 12)) +
  guides(color = guide_legend(title = expression(lambda))) + 
  coord_fixed()

# Bifurcation diagram
bifurcation_diagram <- tibble(lambda = seq(-1.5, 0, length.out = 5000), stable = m + sqrt(A*abs(lambda)), 
                              unstable = m - sqrt(A*abs(lambda))) |>
  pivot_longer(-lambda, names_to = "stability", values_to = "fixed_point") |> 
  ggplot(aes(x = lambda, y = fixed_point, linetype = stability)) + 
  geom_line(linewidth = 1.25) +
  geom_hline(yintercept = m, linetype = "dashed") +
  geom_text(data = tibble(stability = "stable"), aes(x = -1.45, y = m, label = "italic(m)"),
            parse = TRUE, vjust = -0.5, hjust = 0, size = 4.5) +
  labs(linetype = "Type of fixed point") + 
  ylab(expression(mu)) + xlab(expression(lambda)) + 
  theme(axis.title.y = element_text(size = 16, angle = 90),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(margin = margin(t = 20), size = 20))

# ggsave(bifurcation_diagram, path = paste0(getwd(), "/tex_files/figures"), filename = "bifurcation_diagram.jpeg",
#        height = 6, width = 10, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)

# Plot of different nu-parameters
# Parameters used in the simulation
lambda_0 <- -2
m_nu_plot <- -1.5
A_nu_plot <- 0.8
tau <- 100
t_0 <- 10
log_steps <- seq(-1.5, 1.5, length.out = 7)

# Calculating nu_values using exponential function to mirror around 1
nu_values <- round(sort(unique(c(exp(log_steps), 1/exp(log_steps)))), 2)

t_values <- seq(0, tau + t_0, length.out = 500)

# Initialize a data frame to hold the values
nu_plot_data <- data.frame()

# Calculate the function values for each nu
for (nu in nu_values) {
  lambda_t <- lambda_0 * (1 - ifelse(t_values > t_0, (t_values - t_0) / tau, 0))^nu
  fixed_point_stable <- m_nu_plot + sqrt(-lambda_t / A_nu_plot)
  fixed_point_unstable <- m_nu_plot - sqrt(-lambda_t / A_nu_plot)
  # Append to the data frame
  nu_plot_data <- rbind(nu_plot_data, data.frame(t = t_values, lambda = lambda_t,
                                                 nu = as.factor(nu), mu_stable = fixed_point_stable,
                                                 mu_unstable = fixed_point_unstable))
}
nu_plot_data <- as_tibble(nu_plot_data) |> pivot_longer(-c(t, lambda, nu), values_to = "Value", names_to = "Fixed_type")

nu_plot <- ggplot(nu_plot_data, aes(x = t, y = Value, color = nu, linetype = Fixed_type)) +
  geom_line(linewidth = 1.5) +
  geom_hline(yintercept = m_nu_plot, linetype = "dashed") +
  labs(x = "Time (t)",
       y = expression(mu),
       color = expression(nu)) + 
  scale_color_manual(values = thesis_palette) +
  scale_linetype(guide = "none") + 
  guides(color = guide_legend(override.aes = list(linewidth = 5))) +
  theme(plot.margin = unit(c(0, 0, 0, 0.5), "cm"),  panel.border = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 18), legend.title = element_text(size = 20),
        legend.text = element_text(size = 17))
  

# ggsave(nu_plot, path = paste0(getwd(), "/tex_files/figures"), filename = "nu_plot.jpeg",
#        height = 7, width = 11, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)


###------------------------------------------------------------------------###
###------------------------Simulation studies------------------------------###
###------------------------------------------------------------------------###

#Overview of the estimation methods

# General study of performance
# Parameters where all models work
true_param <- c(1.5, 0.4, -0.2, 0.1)
tau <- 30
t_0 <- 30
mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) * sqrt(abs(true_param[3] / true_param[1]))
alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
stationary_part_true_param <- c(alpha0, mu0, true_param[4])
actual_dts <- c(1/5, 1/10, 1/25, 1/100, 1/250) 
dynamic_part_true_param <- c(tau, true_param[1])
numSim <- 50

col_wise_median_OU_like <- matrix(data = NA, nrow = length(actual_dts), ncol = 4)
col_wise_median_OU_score <- matrix(data = NA, nrow = length(actual_dts), ncol = 4)
col_wise_median_sqrt_eq <- matrix(data = NA, nrow = length(actual_dts), ncol = 4)
col_wise_median_sqrt_strang <- matrix(data = NA, nrow = length(actual_dts), ncol = 4)
col_wise_median_linear_eq <- matrix(data = NA, nrow = length(actual_dts), ncol = 4)
col_wise_median_linear_alt_strang <- matrix(data = NA, nrow = length(actual_dts), ncol = 4)
col_wise_median_linear_strang <- matrix(data = NA, nrow = length(actual_dts), ncol = 4)
col_wise_median_t_strang <- matrix(data = NA, nrow = length(actual_dts), ncol = 4)
col_wise_median_F_strang <- matrix(data = NA, nrow = length(actual_dts), ncol = 4)
col_wise_median_Jacobi_Strang <- matrix(data = NA, nrow = length(actual_dts), ncol = 4)

col_wise_median_OU_dynamic <- matrix(data = NA, nrow = length(actual_dts), ncol = 3)
col_wise_median_sqrt_dynamic <- matrix(data = NA, nrow = length(actual_dts), ncol = 3)
col_wise_median_sqrt_alt_dynamic <- matrix(data = NA, nrow = length(actual_dts), ncol = 3)
col_wise_median_linear_dynamic <- matrix(data = NA, nrow = length(actual_dts), ncol = 3)
col_wise_median_linear_alt_dynamic <- matrix(data = NA, nrow = length(actual_dts), ncol = 3)
col_wise_median_t_dynamic <- matrix(data = NA, nrow = length(actual_dts), ncol = 3)
col_wise_median_F_dynamic <- matrix(data = NA, nrow = length(actual_dts), ncol = 3)
col_wise_median_Jacobi_dynamic <- matrix(data = NA, nrow = length(actual_dts), ncol = 3)


set.seed(110524)
for (j in seq_along(actual_dts)){
  
  OU_likelihood_estim <- matrix(NA, nrow = numSim, ncol = 4)
  OU_score_estim <- matrix(NA, nrow = numSim, ncol = 4)
  sqrt_martingale_estim <- matrix(NA, nrow = numSim, ncol = 4)
  sqrt_strang_estim <- matrix(NA, nrow = numSim, ncol = 4)
  linear_martingale_estim <- matrix(NA, nrow = numSim, ncol = 4)
  linear_alt_strang_estim <- matrix(NA, nrow = numSim, ncol = 4)
  linear_strang_estim <- matrix(NA, nrow = numSim, ncol = 4)
  t_strang_estim <- matrix(NA, nrow = numSim, ncol = 4)
  F_strang_estim <- matrix(NA, nrow = numSim, ncol = 4)
  Jacobi_strang_estim <- matrix(NA, nrow = numSim, ncol = 4)
  
  OU_dynamic_estim <- matrix(data = NA, nrow = numSim, ncol = 3)
  sqrt_dynamic_estim <- matrix(data = NA, nrow = numSim, ncol = 3)
  sqrt_alt_dynamic_estim <- matrix(data = NA, nrow = numSim, ncol = 3)
  Linear_dynamic_estim <- matrix(data = NA, nrow = numSim, ncol = 3)
  Linear_alt_dynamic_estim <- matrix(data = NA, nrow = numSim, ncol = 3)
  t_dynamic_estim <- matrix(data = NA, nrow = numSim, ncol = 3)
  F_dynamic_estim <- matrix(data = NA, nrow = numSim, ncol = 3)
  Jacobi_dynamic_estim <- matrix(data = NA, nrow = numSim, ncol = 3)
  
  
  for (i in 1:numSim){
    cat("Grinding iteration", i, "\n")
    success <- FALSE
    while(!success){
    tryCatch({
      
    seed <- sample.int(100000, 1)
    random_noise_start_value <- runif(3, min = .9, 1.1)
    random_noise_start_value_dynamic <- runif(2, min = .9, 1.1)
    
    # Jacobi model
    sim_res_jacobi <- simulate_jacobi_diffusion_tipping_model(step_length = actual_dts[j],
                                                              par = true_param, tau = tau, t_0 = t_0,
                                                              beyond_tipping = 0, seed = seed)
    
    jacobi_result <- optimize_stationary_likelihood(
      likelihood_fun = jacobi_diffusion_strang_splitting,
      data = 2 * asin(sqrt(sim_res_jacobi$X_t[sim_res_jacobi$t<t_0])),
      init_par = stationary_part_true_param * random_noise_start_value,
      delta = actual_dts[j],
      exp_sigma = TRUE, 
      control = list(reltol = sqrt(.Machine$double.eps) / 10))
    if(jacobi_result$objective >= 0){stop("Jacobi failed")}
    
    Jacobi_strang_estim[i, 1:3] <- abs(jacobi_result$par - stationary_part_true_param) /
      stationary_part_true_param
    
    Jacobi_strang_estim[i, 4] <- microbenchmark::microbenchmark(optimize_stationary_likelihood(
      likelihood_fun = jacobi_diffusion_strang_splitting,
      data = 2 * asin(sqrt(sim_res_jacobi$X_t[sim_res_jacobi$t<t_0])),
      init_par = stationary_part_true_param * random_noise_start_value,
      delta = actual_dts[j],
      exp_sigma = TRUE, 
      control = list(reltol = sqrt(.Machine$double.eps) / 10)), times = 1, unit = "us")$time / 1e9
    
    Jacobi_dynamic_estim[i, 1:2] <-  abs(optimize_dynamic_likelihood(
                      likelihood_fun = jacobi_diffusion_transform_dynamic_likelihood,
                      data = 2 * asin(sqrt(sim_res_jacobi$X_t[sim_res_jacobi$t>t_0])),
                      init_par = dynamic_part_true_param * random_noise_start_value_dynamic,
                      delta = actual_dts[j],
                      alpha0 = jacobi_result$par[1],
                      mu0 = jacobi_result$par[2],
                      sigma = jacobi_result$par[3],
                      control = list(reltol = sqrt(.Machine$double.eps) / 10))$par - dynamic_part_true_param) /
      dynamic_part_true_param
    
    Jacobi_dynamic_estim[i, 3] <-  microbenchmark::microbenchmark(optimize_dynamic_likelihood(
      likelihood_fun = jacobi_diffusion_transform_dynamic_likelihood,
      data = 2 * asin(sqrt(sim_res_jacobi$X_t[sim_res_jacobi$t>t_0])),
      init_par = dynamic_part_true_param * random_noise_start_value_dynamic,
      delta = actual_dts[j],
      alpha0 = jacobi_result$par[1],
      mu0 = jacobi_result$par[2],
      sigma = jacobi_result$par[3],
      control = list(reltol = sqrt(.Machine$double.eps) / 10)), times = 1, unit = "us")$time / 1e9


    # OU model
    sim_res_add <- simulate_additive_noise_tipping_model(step_length = actual_dts[j],
                                          par = true_param,
                                          tau = tau, t_0 = t_0,
                                          beyond_tipping = 0, seed = seed)

    OU_result <- optimize_stationary_likelihood(likelihood_fun = OU_likelihood,
                                                data = sim_res_add$X_t[sim_res_add$t < t_0],
                                                init_par = stationary_part_true_param * random_noise_start_value,
                                                delta = actual_dts[j], exp_sigma = FALSE,
                                                control = list(reltol = sqrt(.Machine$double.eps)/10))
    
    OU_likelihood_estim[i, 1:3] <- abs(OU_result$par - stationary_part_true_param) /stationary_part_true_param
    
    OU_likelihood_estim[i, 4] <- microbenchmark::microbenchmark(
      optimize_stationary_likelihood(likelihood_fun = OU_likelihood,
      data = sim_res_add$X_t[sim_res_add$t < t_0],
      init_par = stationary_part_true_param * random_noise_start_value,
      delta = actual_dts[j], exp_sigma = FALSE,
      control = list(reltol = sqrt(.Machine$double.eps) / 10)),
      times = 1, unit = "us")$time / 1e9
    
    OU_score_estim[i, 1:3] <- abs(nleqslv::nleqslv(x = stationary_part_true_param * random_noise_start_value,
                     fn = OU_score,
                     data = sim_res_add$X_t[sim_res_add$t < t_0],
                     delta = actual_dts[j])$x -
                       stationary_part_true_param) / stationary_part_true_param
    
    OU_score_estim[i, 4] <- microbenchmark::microbenchmark(
      nleqslv::nleqslv(x = stationary_part_true_param * random_noise_start_value,
                       fn = OU_score,
                       data = sim_res_add$X_t[sim_res_add$t < t_0],
                       delta = actual_dts[j])$x,
      times = 1, unit = "us")$time / 1e9
    
    OU_dynamic_estim[i, 1:2] <-  abs(optimize_dynamic_likelihood(
      likelihood_fun = OU_dynamic_likelihood,
      data = sim_res_add$X_t[sim_res_add$t > t_0],
      init_par = dynamic_part_true_param,
      delta = actual_dts[j],
      alpha0 = OU_result$par[1],
      mu0 = OU_result$par[2],
      sigma = OU_result$par[3],
      control = list(reltol = sqrt(.Machine$double.eps) / 10))$par - dynamic_part_true_param) / 
      dynamic_part_true_param
    
    OU_dynamic_estim[i, 3] <-  microbenchmark::microbenchmark(optimize_dynamic_likelihood(
      likelihood_fun = OU_dynamic_likelihood,
      data = sim_res_add$X_t[sim_res_add$t > t_0],
      init_par = dynamic_part_true_param * random_noise_start_value_dynamic,
      delta = actual_dts[j],
      alpha0 = OU_result$par[1],
      mu0 = OU_result$par[2],
      sigma = OU_result$par[3],
      control = list(reltol = sqrt(.Machine$double.eps) / 10)), times = 1, unit = "us")$time / 1e9

    
    # Sqrt-model
    
    sim_res_sqrt <- simulate_squareroot_noise_tipping_model(step_length = actual_dts[j],
                                            par = true_param, tau = tau, t_0 = t_0,
                                            beyond_tipping = 0, seed = seed)

    sqrt_martingale_estim[i, 1:3] <- abs(CIR_quadratic_martingale(sim_res_sqrt$X_t[sim_res_sqrt$t < t_0],
                                                                  actual_dts[j]) -
                                        stationary_part_true_param) /stationary_part_true_param
    
    sqrt_martingale_estim[i, 4] <- median(microbenchmark::microbenchmark(
      CIR_quadratic_martingale(sim_res_sqrt$X_t[sim_res_sqrt$t < t_0], actual_dts[j]),
      times = 1, unit = "us")$time / 1e9)
    
    sqrt_strang_result <- optimize_stationary_likelihood(likelihood_fun = CIR_strang_splitting,
                                                         data = 2 * sqrt(sim_res_sqrt$X_t[sim_res_sqrt$t < t_0]),
                                                         init_par = stationary_part_true_param * random_noise_start_value,
                                                         delta = actual_dts[j], exp_sigma = TRUE,
                                                         control = list(reltol = sqrt(.Machine$double.eps)/10))
    
    sqrt_strang_estim[i, 1:3] <- abs(sqrt_strang_result$par -
                                    stationary_part_true_param) / stationary_part_true_param
    
    sqrt_strang_estim[i, 4] <- microbenchmark::microbenchmark(
      optimize_stationary_likelihood(likelihood_fun = CIR_strang_splitting,
                                     data = 2 * sqrt(sim_res_sqrt$X_t[sim_res_sqrt$t < t_0]),
                                     init_par = stationary_part_true_param * random_noise_start_value,
                                     delta = actual_dts[j], exp_sigma = TRUE,
                                     control = list(reltol = sqrt(.Machine$double.eps)/10)),
                                     times = 1, unit = "us")$time / 1e9
    
    sqrt_dynamic_estim[i, 1:2] <- abs(optimize_dynamic_likelihood(
                likelihood_fun = CIR_transform_dynamic_likelihood,
                data = 2 * sqrt(sim_res_sqrt$X_t[sim_res_sqrt$t > t_0]),
                init_par = dynamic_part_true_param * random_noise_start_value_dynamic,
                delta = actual_dts[j],
                alpha0 = sqrt_strang_result$par[1],
                mu0 = sqrt_strang_result$par[2],
                sigma = sqrt_strang_result$par[3],
                control = list(reltol = sqrt(.Machine$double.eps) / 10))$par - dynamic_part_true_param) /
      dynamic_part_true_param

    sqrt_dynamic_estim[i, 3] <-  microbenchmark::microbenchmark(optimize_dynamic_likelihood(
      likelihood_fun = CIR_transform_dynamic_likelihood,
      data = 2 * sqrt(sim_res_sqrt$X_t[sim_res_sqrt$t > t_0]),
      init_par = dynamic_part_true_param * random_noise_start_value_dynamic,
      delta = actual_dts[j],
      alpha0 = sqrt_strang_result$par[1],
      mu0 = sqrt_strang_result$par[2],
      sigma = sqrt_strang_result$par[3],
      control = list(reltol = sqrt(.Machine$double.eps) / 10)),
      times = 1, unit = "us")$time / 1e9
    
    
    sqrt_alt_dynamic_estim[i, 1:2] <-  abs(optimize_dynamic_likelihood(
      likelihood_fun = CIR_dynamic_likelihood,
      data = sim_res_sqrt$X_t[sim_res_sqrt$t > t_0],
      init_par = dynamic_part_true_param * random_noise_start_value_dynamic,
      delta = actual_dts[j],
      alpha0 = sqrt_strang_result$par[1],
      mu0 = sqrt_strang_result$par[2],
      sigma = sqrt_strang_result$par[3],
      control = list(reltol = sqrt(.Machine$double.eps) / 10))$par - dynamic_part_true_param) /
      dynamic_part_true_param
    
    sqrt_alt_dynamic_estim[i, 3] <-  microbenchmark::microbenchmark(optimize_dynamic_likelihood(
      likelihood_fun = CIR_dynamic_likelihood,
      data = sim_res_sqrt$X_t[sim_res_sqrt$t > t_0],
      init_par = dynamic_part_true_param * random_noise_start_value_dynamic,
      delta = actual_dts[j],
      alpha0 = sqrt_strang_result$par[1],
      mu0 = sqrt_strang_result$par[2],
      sigma = sqrt_strang_result$par[3],
      control = list(reltol = sqrt(.Machine$double.eps) / 10)),
      times = 1, unit = "us")$time / 1e9
    
    # Linear noise model
    sim_res_linear <- simulate_linear_noise_tipping_model(step_length = actual_dts[j],
                                        par = true_param, tau = tau, t_0 = t_0,
                                        beyond_tipping = 0, seed = seed)
    

    linear_martingale_estim[i, 1:3] <- abs(nleqslv::nleqslv(x = stationary_part_true_param * random_noise_start_value,
                     fn = mean_reverting_GBM_martingale,
                     data = sim_res_linear$X_t[sim_res_linear$t < t_0],
                     delta = actual_dts[j])$x -
                       stationary_part_true_param) /stationary_part_true_param

    linear_martingale_estim[i, 4] <- microbenchmark::microbenchmark(
      nleqslv::nleqslv(x = stationary_part_true_param * random_noise_start_value,
                       fn = mean_reverting_GBM_martingale,
                       data = sim_res_linear$X_t[sim_res_linear$t < t_0],
                       delta = actual_dts[j])$x, times = 1, unit = "us")$time / 1e9
    
    linear_strang_result <- optimize_stationary_likelihood(mean_reverting_GBM_strang, 
                          log(sim_res_linear$X_t[sim_res_linear$t<t_0]),
                          init_par = stationary_part_true_param * random_noise_start_value,
                          delta = actual_dts[j],
                          exp_sigma = FALSE,
                          control = list(reltol = sqrt(.Machine$double.eps)/10))

    linear_strang_estim[i, 1:3] <- abs(linear_strang_result$par - stationary_part_true_param) /
      stationary_part_true_param
    
    linear_strang_estim[i, 4] <- microbenchmark::microbenchmark(
      optimize_stationary_likelihood(mean_reverting_GBM_strang, 
                       log(sim_res_linear$X_t[sim_res_linear$t<t_0]),
                       init_par = stationary_part_true_param * random_noise_start_value,
                       delta = actual_dts[j],
                       exp_sigma = FALSE,
                       control = list(reltol = sqrt(.Machine$double.eps)/10)),
      times = 1, unit = "us")$time / 1e9

    linear_alt_strang_result <- optimize_stationary_likelihood(mean_reverting_GBM_alt_strang,
                                   sim_res_linear$X_t[sim_res_linear$t<t_0],
                                   init_par = stationary_part_true_param * random_noise_start_value,
                                   delta = actual_dts[j],
                                   exp_sigma = FALSE,
                                   control = list(reltol = sqrt(.Machine$double.eps)/10))

    linear_alt_strang_estim[i, 1:3] <- abs(linear_alt_strang_result$par - stationary_part_true_param) /
      stationary_part_true_param
    
    linear_alt_strang_estim[i, 4] <- microbenchmark::microbenchmark(
      optimize_stationary_likelihood(mean_reverting_GBM_alt_strang,
                       sim_res_linear$X_t[sim_res_linear$t<t_0],
                       init_par = stationary_part_true_param * random_noise_start_value,
                       delta = actual_dts[j],
                       exp_sigma = FALSE,
                       control = list(reltol = sqrt(.Machine$double.eps)/10)), times = 1, unit = "us")$time / 1e9

    Linear_alt_dynamic_estim[i, 1:2] <- abs(optimize_dynamic_likelihood(
    likelihood_fun = mean_reverting_GBM_dynamic_likelihood,
    data = sim_res_linear$X_t[sim_res_linear$t > t_0],
    init_par = dynamic_part_true_param * random_noise_start_value_dynamic,
    delta = actual_dts[j],
    alpha0 = linear_alt_strang_result$par[1],
    mu0 =  linear_alt_strang_result$par[2],
    sigma =  linear_alt_strang_result$par[3],
    control = list(reltol = sqrt(.Machine$double.eps) / 10))$par - dynamic_part_true_param) /
      dynamic_part_true_param
    
    Linear_alt_dynamic_estim[i, 3] <- microbenchmark::microbenchmark(
      optimize_dynamic_likelihood(
        likelihood_fun = mean_reverting_GBM_dynamic_likelihood,
        data = sim_res_linear$X_t[sim_res_linear$t > t_0],
        init_par = dynamic_part_true_param * random_noise_start_value_dynamic,
        delta = actual_dts[j],
        alpha0 = linear_alt_strang_result$par[1],
        mu0 =  linear_alt_strang_result$par[2],
        sigma =  linear_alt_strang_result$par[3],
        control = list(reltol = sqrt(.Machine$double.eps) / 10)),
      times = 1, unit = "us")$time / 1e9

    Linear_dynamic_estim[i, 1:2] <- abs(optimize_dynamic_likelihood(
          likelihood_fun = mean_reverting_GBM_transform_dynamic_likelihood,
          data = log(sim_res_linear$X_t[sim_res_linear$t > t_0]),
          init_par = dynamic_part_true_param * random_noise_start_value_dynamic,
          delta = actual_dts[j],
          alpha0 = linear_strang_result$par[1],
          mu0 = linear_strang_result$par[2],
          sigma = linear_strang_result$par[3],
          control = list(reltol = sqrt(.Machine$double.eps)/10))$par - dynamic_part_true_param) /
      dynamic_part_true_param
    
    Linear_dynamic_estim[i, 3] <- microbenchmark::microbenchmark(
      optimize_dynamic_likelihood(
        likelihood_fun = mean_reverting_GBM_transform_dynamic_likelihood,
        data = log(sim_res_linear$X_t[sim_res_linear$t > t_0]),
        init_par = dynamic_part_true_param * random_noise_start_value_dynamic,
        delta = actual_dts[j],
        alpha0 = linear_strang_result$par[1],
        mu0 = linear_strang_result$par[2],
        sigma = linear_strang_result$par[3],
        control = list(reltol = sqrt(.Machine$double.eps)/10)),
      times = 1, unit = "us")$time / 1e9
    
    
    # t-diffusion model
    sim_res_t_distribution <- simulate_t_distribution_tipping_model(step_length = actual_dts[j],
                                          par = true_param, tau = tau, t_0 = t_0,
                                          beyond_tipping = 0, seed = seed)
    
    t_strang_result <- optimize_stationary_likelihood(
      likelihood_fun = t_diffusion_strang_splitting,
      data = asinh(sim_res_t_distribution$X_t[sim_res_t_distribution$t < t_0]),
      init_par = stationary_part_true_param * random_noise_start_value,
      delta = actual_dts[j],
      exp_sigma = TRUE,
      control = list(reltol = sqrt(.Machine$double.eps)/10))
    
    t_strang_estim[i, 1:3] <- abs(t_strang_result$par - stationary_part_true_param) / stationary_part_true_param
    
    t_strang_estim[i, 4] <- microbenchmark::microbenchmark(optimize_stationary_likelihood(
      likelihood_fun = t_diffusion_strang_splitting,
      data = asinh(sim_res_t_distribution$X_t[sim_res_t_distribution$t < t_0]),
      init_par = stationary_part_true_param * random_noise_start_value,
      delta = actual_dts[j],
      exp_sigma = TRUE,
      control = list(reltol = sqrt(.Machine$double.eps)/10)),
      times = 1, unit = "us")$time / 1e9
    
    t_dynamic_estim[i, 1:2] <-  abs(optimize_dynamic_likelihood(
                      likelihood_fun = t_transform_dynamic_likelihood,
                      data = asinh(sim_res_t_distribution$X_t[sim_res_t_distribution$t > t_0]),
                      init_par = dynamic_part_true_param * random_noise_start_value_dynamic,
                      delta = actual_dts[j],
                      alpha0 = t_strang_result$par[1],
                      mu0 = t_strang_result$par[2],
                      sigma = t_strang_result$par[3],
                      control = list(reltol = sqrt(.Machine$double.eps) / 10))$par - dynamic_part_true_param) / 
      dynamic_part_true_param
    
    t_dynamic_estim[i, 3] <- microbenchmark::microbenchmark(
      optimize_dynamic_likelihood(
        likelihood_fun = t_transform_dynamic_likelihood,
        data = asinh(sim_res_t_distribution$X_t[sim_res_t_distribution$t > t_0]),
        init_par = dynamic_part_true_param * random_noise_start_value_dynamic,
        delta = actual_dts[j],
        alpha0 = t_strang_result$par[1],
        mu0 = t_strang_result$par[2],
        sigma = t_strang_result$par[3],
        control = list(reltol = sqrt(.Machine$double.eps) / 10)),
      times = 1, unit = "us")$time / 1e9
    

    # F-diffusion model
    F_sim_dynamic <- simulate_F_distribution_tipping_model(step_length = actual_dts[j],
                                          par = true_param, tau = tau, t_0 = t_0,
                                          beyond_tipping = 0, seed = seed)
    
    F_strang_result <- optimize_stationary_likelihood(
      likelihood_fun = F_diffusion_strang_splitting,
      data = 2 * asinh(sqrt(F_sim_dynamic$X_t[F_sim_dynamic$t<t_0])),
      init_par = stationary_part_true_param * random_noise_start_value,
      delta = actual_dts[j],
      exp_sigma = TRUE,
      control = list(reltol = sqrt(.Machine$double.eps)/10))
    
   F_strang_estim[i, 1:3] <-  abs(F_strang_result$par - stationary_part_true_param) / stationary_part_true_param
    
    F_strang_estim[i, 4] <- microbenchmark::microbenchmark(optimize_stationary_likelihood(
     likelihood_fun = F_diffusion_strang_splitting,
     data = 2 * asinh(sqrt(F_sim_dynamic$X_t[F_sim_dynamic$t<t_0])),
     init_par = stationary_part_true_param * random_noise_start_value,
     delta = actual_dts[j],
     exp_sigma = TRUE,
     control = list(reltol = sqrt(.Machine$double.eps)/10)), times = 1, unit = "us")$time / 1e9
    
    F_dynamic_estim[i, 1:2] <-  abs(optimize_dynamic_likelihood(
                    likelihood_fun = F_transform_dynamic_likelihood,
                    data = 2 * asinh(sqrt(F_sim_dynamic$X_t[F_sim_dynamic$t>t_0])),
                    init_par = dynamic_part_true_param * random_noise_start_value_dynamic,
                    delta = actual_dts[j],
                    alpha0 = F_strang_result$par[1],
                    mu0 = F_strang_result$par[2],
                    sigma = F_strang_result$par[3],
                    control = list(reltol = sqrt(.Machine$double.eps) / 10))$par - dynamic_part_true_param) / 
      dynamic_part_true_param
    
    F_dynamic_estim[i, 3] <- microbenchmark::microbenchmark(
      optimize_dynamic_likelihood(
        likelihood_fun = F_transform_dynamic_likelihood,
        data = 2 * asinh(sqrt(F_sim_dynamic$X_t[F_sim_dynamic$t>t_0])),
        init_par = dynamic_part_true_param * random_noise_start_value_dynamic,
        delta = actual_dts[j],
        alpha0 = F_strang_result$par[1],
        mu0 = F_strang_result$par[2],
        sigma = F_strang_result$par[3],
        control = list(reltol = sqrt(.Machine$double.eps) / 10)),
      times = 1, unit = "us")$time / 1e9

    success <- TRUE
    }, error = function(e) {
      cat("Error occurred at iteration", i, ":", conditionMessage(e), "\n")
      
      Sys.sleep(1)
    })
    }
  }
  col_wise_median_OU_like[j, ] <- colwise_median(OU_likelihood_estim)
  col_wise_median_OU_score[j, ] <- colwise_median(OU_score_estim)
  col_wise_median_sqrt_eq[j, ] <- colwise_median(sqrt_martingale_estim)
  col_wise_median_sqrt_strang[j, ] <- colwise_median(sqrt_strang_estim)
  col_wise_median_linear_eq[j, ] <- colwise_median(linear_martingale_estim)
  col_wise_median_linear_alt_strang[j, ] <- colwise_median(linear_alt_strang_estim)
  col_wise_median_linear_strang[j, ] <- colwise_median(linear_strang_estim)
  col_wise_median_t_strang[j, ] <- colwise_median(t_strang_estim)
  col_wise_median_F_strang[j, ] <- colwise_median(F_strang_estim)
  col_wise_median_Jacobi_Strang[j, ] <- colwise_median(Jacobi_strang_estim)
  
  col_wise_median_OU_dynamic[j, ] <- colwise_median(OU_dynamic_estim)
  col_wise_median_sqrt_dynamic[j, ] <- colwise_median(sqrt_dynamic_estim)
  col_wise_median_sqrt_alt_dynamic[j, ] <-colwise_median(sqrt_alt_dynamic_estim)
  col_wise_median_linear_dynamic[j, ] <- colwise_median(Linear_dynamic_estim)
  col_wise_median_linear_alt_dynamic[j, ] <- colwise_median(Linear_alt_dynamic_estim)
  col_wise_median_t_dynamic[j, ] <- colwise_median(t_dynamic_estim)
  col_wise_median_F_dynamic[j, ] <- colwise_median(F_dynamic_estim)
  col_wise_median_Jacobi_dynamic[j, ] <- colwise_median(Jacobi_dynamic_estim)
  
}

## Stationary analysis
OU_likelihood_tibble <- as_tibble(col_wise_median_OU_like) |>
  mutate(delta = actual_dts, Model = "Additive", Type = "Likelihood", Method = "MLE")

OU_score_tibble <- as_tibble(col_wise_median_OU_score) |>
  mutate(delta = actual_dts, Model = "Additive", Type = "Estimation eq.", Method = "Score")

sqrt_estimation_equation_tibble <- as_tibble(col_wise_median_sqrt_eq) |>
  mutate(delta = actual_dts, Model = "Square-root", Type = "Estimation eq.", Method = "Martingale")

sqrt_Strang_tibble <- as_tibble(col_wise_median_sqrt_strang) |>
  mutate(delta = actual_dts, Model = "Square-root", Type = "Likelihood", Method = "Strang")

Linear_estimation_equation_tibble <- as_tibble(col_wise_median_linear_eq) |>
  mutate(delta = actual_dts, Model = "Linear", Type = "Estimation eq.", Method = "Martingale")

Linear_alt_Strang_tibble <- as_tibble(col_wise_median_linear_alt_strang) |>
  mutate(delta = actual_dts, Model = "Linear", Type = "Likelihood (Alt.)", Method = "Strang (Alternative)")

Linear_Strang_tibble <- as_tibble(col_wise_median_linear_strang)  |>
  mutate(delta = actual_dts, Model = "Linear", Type = "Likelihood", Method = "Strang")

t_diffusion_Strang_tibble <- as_tibble(col_wise_median_t_strang)  |>
  mutate(delta = actual_dts, Model = "t-diffusion", Type = "Likelihood", Method = "Strang")

F_diffusion_Strang_tibble <- as_tibble(col_wise_median_F_strang)  |>
  mutate(delta = actual_dts, Model = "F-diffusion", Type = "Likelihood", Method = "Strang")

Jacobi_diffusion_Strang_tibble <- as_tibble(col_wise_median_Jacobi_Strang) |>
  mutate(delta = actual_dts, Model = "Jacobi-diffusion", Type = "Likelihood", Method = "Strang")

Stationary_estimation_all <- bind_rows(OU_likelihood_tibble, OU_score_tibble, sqrt_estimation_equation_tibble,
          sqrt_Strang_tibble, Linear_estimation_equation_tibble,
          Linear_alt_Strang_tibble, Linear_Strang_tibble, t_diffusion_Strang_tibble,
          F_diffusion_Strang_tibble, Jacobi_diffusion_Strang_tibble) |>
  rename(alpha0 = V1, mu0 = V2, sigma = V3, running_time = V4) |> 
  mutate(Model = factor(Model), Type = factor(Type), Method = factor(Method)) |> 
  pivot_longer(-c(delta, Model, Type, Method), names_to = "Parameter", values_to = "ARE") |> 
  mutate(Model = factor(Model, levels = c("Additive", "Square-root", "Linear", 
                                   "t-diffusion", "F-diffusion", "Jacobi-diffusion")),
         Type = factor(Type, levels = c("Likelihood", "Estimation eq.", "Likelihood (Alt.)")))

# if(!file.exists("data/Stationary_estimation_all.csv")){
#   utils::write.table(Stationary_estimation_all, file="data/Stationary_estimation_all.csv",
#                      sep = ",", row.names = FALSE)
# } else{
#   Stationary_estimation_all <- read_csv("data/Stationary_estimation_all.csv")
# }

Stationary_estimation_all <- Stationary_estimation_all |> 
  mutate(Model = factor(Model, levels = c("Additive", "Square-root", "Linear", 
       "t-diffusion", "F-diffusion", "Jacobi-diffusion")),
       Type = factor(Type, levels = c("Likelihood", "Estimation eq.", "Likelihood (Alt.)")))

parameter_precision_stationary <- Stationary_estimation_all |> filter(Parameter != "running_time") |> 
  ggplot(aes(x = t_0 * 1 / delta, y = ARE, color = Parameter, linetype = Type)) +
  geom_line(linewidth = 1.5) +  
  geom_point(size = 3) +
  facet_wrap(~Model) + 
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_log10(breaks = t_0 * 1 / actual_dts) +
  scale_color_manual(values = thesis_palette[4:6],
                     labels = expression(alpha*phantom(.)[0], mu*phantom(.)[0], sigma)) +
  labs(x = "N", y = "Absolute Relative Error") +  
  theme(strip.text = element_text(face = "bold", size = 18), panel.spacing = unit(2, "lines"),
        axis.title.x = element_text(face = "bold", size = 18),
        axis.title.y = element_text(face = "bold", size = 18),
        legend.text  = element_text(size = 18),
        legend.title = element_text(face = "bold", size = 18),
        axis.text = element_text(face = "bold", size = 14)) +
  guides(color = guide_legend(override.aes = list(shape = NA, linewidth = 5)), 
  linetype = guide_legend(override.aes = list(linewidth = 0.75))) +
  scale_linetype_manual(values=c("solid", "dashed", "dotted"))

# ggsave(parameter_precision_stationary, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "parameter_precision_stationary.jpeg",
#        height = 6, width = 14, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)

estimation_duration_stationary <- Stationary_estimation_all |>
  filter(Parameter == "running_time") |>
  ggplot(aes(x = t_0 * 1 / delta, y = ARE, linetype = Type, color = Model)) +
  geom_line(linewidth = 2.25) +  
  geom_point(size = 3.5) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_log10(breaks = t_0 * 1 / actual_dts) +
  labs(x = "N", y = "Running time (s)") +  
  theme(axis.title.x = element_text(face = "bold", size = 18), 
        axis.title.y = element_text(face = "bold", size = 18),
  legend.text  = element_text(size = 17),
  legend.title = element_text(face = "bold", size = 18),
  axis.text = element_text(face = "bold", size = 16),
  panel.border = element_blank(),
  ) +
  scale_color_manual(values = thesis_palette[1:6]) +
  guides(color = guide_legend(override.aes = list(shape = NA, linewidth = 5)), 
         linetype = guide_legend(override.aes = list(linewidth = 0.75))) +
  scale_linetype_manual(values=c("solid", "dashed", "dotted"))
  
# ggsave(estimation_duration_stationary, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "estimation_duration_stationary.jpeg",
#        height = 6, width = 13, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)

## Dynamic analysis
OU_dynamic_tibble <- as_tibble(col_wise_median_OU_dynamic) |>
  mutate(delta = actual_dts, Model = "Additive", Type = "Strang")

sqrt_dynamic_tibble <- as_tibble(col_wise_median_sqrt_dynamic) |>
  mutate(delta = actual_dts, Model = "Square-root", Type = "Strang")

sqrt_alt_dynamic_tibble <- as_tibble(col_wise_median_sqrt_alt_dynamic) |>
  mutate(delta = actual_dts, Model = "Square-root", Type = "Strang (Alt.)")

linear_dynamic_tibble <- as_tibble(col_wise_median_linear_dynamic) |>
  mutate(delta = actual_dts, Model = "Linear", Type = "Strang")

linear_alt_dynamic_tibble <- as_tibble(col_wise_median_linear_alt_dynamic) |>
  mutate(delta = actual_dts, Model = "Linear", Type = "Strang (Alt.)")

t_dynamic_tibble <- as_tibble(col_wise_median_t_dynamic) |>
  mutate(delta = actual_dts, Model = "t-diffusion", Type = "Strang")

F_dynamic_tibble <- as_tibble(col_wise_median_F_dynamic) |>
  mutate(delta = actual_dts, Model = "F-diffusion", Type = "Strang")

Jacobi_dynamic_tibble <- as_tibble(col_wise_median_Jacobi_dynamic) |>
  mutate(delta = actual_dts, Model = "Jacobi-diffusion", Type = "Strang")
 
Dynamic_estimation_all <- bind_rows(OU_dynamic_tibble, sqrt_dynamic_tibble, sqrt_alt_dynamic_tibble,
                                       linear_dynamic_tibble, linear_alt_dynamic_tibble,
                                       t_dynamic_tibble, F_dynamic_tibble, Jacobi_dynamic_tibble) |>
                rename(tau = V1, A = V2, running_time = V3) |> 
                mutate(Model = factor(Model), Type = factor(Type)) |> 
                pivot_longer(-c(delta, Model, Type), names_to = "Parameter", values_to = "ARE") |> 
                mutate(Model = factor(Model, levels = c("Additive", "Square-root", "Linear", 
                    "t-diffusion", "F-diffusion", "Jacobi-diffusion")),
Type = factor(Type, levels = c("Strang", "Strang (Alt.)")))

# if(!file.exists("data/Dynamic_estimation_all.csv")){
#   utils::write.table(Dynamic_estimation_all, file="data/Dynamic_estimation_all.csv",
#                      sep = ",", row.names = FALSE)
# } else{
#   Dynamic_estimation_all <- read_csv("data/Dynamic_estimation_all.csv")
# }

Dynamic_estimation_all <- Dynamic_estimation_all |> 
mutate(Model = factor(Model, levels = c("Additive", "Square-root", "Linear", 
                                        "t-diffusion", "F-diffusion", "Jacobi-diffusion")),
       Type = factor(Type, levels = c("Strang", "Strang (Alt.)")))


parameter_precision_dynamic <- Dynamic_estimation_all |> filter(Parameter != "running_time") |> 
  ggplot(aes(x = tau / delta, y = ARE, color = Parameter, linetype = Type)) +
  geom_line(linewidth = 2.25) +  
  geom_point(size = 3) +
  facet_wrap(~Model) + 
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_log10(breaks = tau * 1 / actual_dts) +
  scale_color_manual(values = thesis_palette[4:6],
                     labels = expression(A, tau*phantom(.)[c])) +
  labs(x = "N", y = "Absolute Relative Error") +  
  theme(strip.text = element_text(face = "bold", size = 18), panel.spacing = unit(2, "lines"),
        axis.title.x = element_text(face = "bold", size = 18),
        axis.title.y = element_text(face = "bold", size = 18),
        legend.text  = element_text(size = 18),
        legend.title = element_text(face = "bold", size = 18),
        axis.text = element_text(face = "bold", size = 14)) + 
  guides(color = guide_legend(override.aes = list(shape = NA, linewidth = 5)), 
         linetype = guide_legend(override.aes = list(linewidth = 0.75))) +
  scale_linetype_manual(values=c("solid", "dashed", "dotted"))

# ggsave(parameter_precision_dynamic, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "parameter_precision_dynamic.jpeg",
#        height = 6, width = 13, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)

estimation_duration_dynamic <-  Dynamic_estimation_all |> filter(Parameter == "running_time") |>
  ggplot(aes(x = tau * 1 / delta, y = ARE, linetype = Type, color = Model)) +
  geom_line(linewidth = 1.95) +  
  geom_point(size = 3.5) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_log10(breaks = tau * 1 / actual_dts) +
  labs(x = "N", y = "Running time (s)") + 
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold"),
        legend.text  = element_text(size = 14),
        axis.text = element_text(face = "bold", size = 14),
        panel.border = element_blank()) +
  scale_color_manual(values = thesis_palette[1:6]) +
  guides(color = guide_legend(override.aes = list(shape = NA, linewidth = 5)), 
         linetype = guide_legend(override.aes = list(linewidth = 0.75))) +
  scale_linetype_manual(values=c("solid", "dashed", "dotted"))

# ggsave(estimation_duration_dynamic, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "estimation_duration_dynamic.jpeg",
#        height = 6, width = 13, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)


##### Estimation of the nu-parameter
true_param_base <- c(0.87, -1.51, -3, sqrt(0.15))
mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) * sqrt(abs(true_param[3] / true_param[1]))
alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
stationary_part_true_param <- c(alpha0, mu0, true_param[4])
actual_dt <- c(1/3, 1/5, 1/25, 1/50, 1/100)
tau <- 132
t_0 <- 54

nus <- c(0.75, 0.9, 1, 1/0.9, 1/0.75)
numSim <- 50


t_error_counts <- matrix(data = NA, nrow = length(actual_dt), ncol = length(nus))
both_error_count <- matrix(data = NA, nrow = length(actual_dt), ncol = length(nus))
add_error_counts <- matrix(data = NA, nrow = length(actual_dt), ncol = length(nus))
median_add_dynamic_nus <- matrix(data = NA, nrow = length(actual_dt), ncol = length(nus))
median_add_dynamic_A <- matrix(data = NA, nrow = length(actual_dt), ncol = length(nus))
median_add_dynamic_tau <- matrix(data = NA, nrow = length(actual_dt), ncol = length(nus))
median_t_dynamic_nus <- matrix(data = NA, nrow = length(actual_dt), ncol = length(nus))
median_t_dynamic_A <- matrix(data = NA, nrow = length(actual_dt), ncol = length(nus))
median_t_dynamic_tau <- matrix(data = NA, nrow = length(actual_dt), ncol = length(nus))
set.seed(130524)

for(n in seq_along(actual_dt)){
t_dynamic_nus <- matrix(data = NA, nrow = numSim, ncol = length(nus))
add_dynamic_nus <- matrix(data = NA, nrow = numSim, ncol = length(nus)) 
t_dynamic_A <- matrix(data = NA, nrow = numSim, ncol = length(nus))
add_dynamic_A <- matrix(data = NA, nrow = numSim, ncol = length(nus))
t_dynamic_tau <- matrix(data = NA, nrow = numSim, ncol = length(nus))
add_dynamic_tau <- matrix(data = NA, nrow = numSim, ncol = length(nus))
t_error_count_vec <- numeric(length = length(nus))
add_error_count_vec <- numeric(length = length(nus))
both_error_count_vec <- numeric(length = length(nus))
for (i in seq_along(nus)){
cat("Grinding nu number: ", i, "\n")
true_param <- c(true_param_base, nus[i])
dynamic_part_true_param <- c(tau, true_param[1], true_param[5])
  for(j in 1:numSim){
    cat("Currently at: ", j, " in the simulation\n")
    success <- FALSE
    while(!success){
    t_failed <- FALSE
    OU_failed <- FALSE
    seed <- sample.int(100000, size = 1)
    random_noise_start <- runif(2, min = 0.9, max = 1.1)
    sim_res_t_distribution <- simulate_t_distribution_tipping_model(
    step_length = actual_dt[n],
    par = true_param,
    tau = tau, t_0 = t_0,
    beyond_tipping = 0, seed = seed)
  
    sim_res_add <- simulate_additive_noise_tipping_model(
    step_length = actual_dt[n],
     par = true_param,
     tau = tau, t_0 = t_0,
     beyond_tipping = 0, seed = seed)
    

    t_dist_estim_param <- optimize_stationary_likelihood(
    likelihood_fun = t_diffusion_strang_splitting,
    data = asinh(sim_res_t_distribution$X_t[sim_res_t_distribution$t < t_0]),
    init_par = stationary_part_true_param,
    delta = actual_dt[n],
    exp_sigma = TRUE)$par
    
    additive_estim_param <- optimize_stationary_likelihood(
    likelihood_fun = OU_likelihood,
    data = sim_res_add$X_t[sim_res_add$t < t_0],
    init_par = stationary_part_true_param,
    delta = actual_dt[n],
    exp_sigma = TRUE)$par
    
    tryCatch({
      t_result <- optimize_dynamic_likelihood(
        likelihood_fun = t_transform_dynamic_likelihood,
        data = asinh(sim_res_t_distribution$X_t[sim_res_t_distribution$t > t_0]),
        init_par = c(c(tau, true_param[1]) * random_noise_start, 1),
        delta = actual_dt[n],
        alpha0 = t_dist_estim_param[1],
        mu0 = t_dist_estim_param[2],
        sigma = t_dist_estim_param[3],
        control = list(reltol = sqrt(.Machine$double.eps) / 1000))
      if(t_result$objective >=0){stop("t failed")}
    t_dynamic_tau[j, i] <- abs(t_result$par[1] - tau) / tau
    t_dynamic_A[j, i] <-   abs(t_result$par[2] - true_param[1]) / true_param[1]
    t_dynamic_nus[j, i] <- abs(t_result$par[3]  - nus[i]) / nus[i]
    }, error = function(e){
      cat("t failed \n")
      t_error_count_vec[i] <<- t_error_count_vec[i] + 1
      t_failed <<- TRUE
    })
    
    tryCatch({
      add_result <- optimize_dynamic_likelihood(
        likelihood_fun = OU_dynamic_likelihood,
        data  = sim_res_add$X_t[sim_res_add$t > t_0],
        init_par = c(c(tau, true_param[1]) * random_noise_start, 1),
        delta = actual_dt[n],
        alpha0 = additive_estim_param[1],
        mu0 = additive_estim_param[2],
        sigma = additive_estim_param[3],
        control = list(reltol = sqrt(.Machine$double.eps) / 1000))
      add_dynamic_tau[j, i] <- abs(add_result$par[1] - tau) / tau
      add_dynamic_A[j, i] <-   abs(add_result$par[2] - true_param[1]) / true_param[1]
      add_dynamic_nus[j, i] <- abs(add_result$par[3]  - nus[i]) / nus[i]
    }, error = function(e){
      cat("OU failed \n")
      add_error_count_vec[i] <<- add_error_count_vec[i] + 1
      OU_failed <<- TRUE
    })

    if(OU_failed && t_failed){
      print("Both failed")
      both_error_count_vec[i] <- both_error_count_vec[i] + 1
    }
    if(!t_failed && !OU_failed){
    success <- TRUE
    }
    }
  }
}
t_error_counts[n, ] <- t_error_count_vec
add_error_counts[n, ] <- add_error_count_vec
both_error_count[n, ]<- both_error_count_vec
median_add_dynamic_nus[n, ] <- colwise_median(add_dynamic_nus)
median_t_dynamic_nus[n, ] <-  colwise_median(t_dynamic_nus)
median_add_dynamic_tau[n, ] <- colwise_median(add_dynamic_tau)
median_t_dynamic_tau[n, ] <-  colwise_median(t_dynamic_tau)
median_add_dynamic_A[n, ] <-  colwise_median(add_dynamic_A)
median_t_dynamic_A[n, ] <- colwise_median(t_dynamic_A)
}

add_dynamic_nus_tibble <- as_tibble(median_add_dynamic_nus) |> 
  mutate(Model = "Additive", N = tau / actual_dt, Parameter = "nu")

t_dynamic_nus_tibble <- as_tibble(median_t_dynamic_nus) |> 
  mutate(Model = "t-diffusion", N = tau / actual_dt, Parameter = "nu")

add_dynamic_tau_tibble <- as_tibble(median_add_dynamic_tau) |> 
  mutate(Model = "Additive", N = tau / actual_dt, Parameter = "tau")

t_dynamic_tau_tibble <- as_tibble(median_t_dynamic_tau) |> 
  mutate(Model = "t-diffusion", N = tau / actual_dt, Parameter = "tau")

add_dynamic_A_tibble <- as_tibble(median_add_dynamic_A) |> 
  mutate(Model = "Additive", N = tau / actual_dt, Parameter = "A")

t_dynamic_A_tibble <- as_tibble(median_t_dynamic_A) |> 
  mutate(Model = "t-diffusion", N = tau / actual_dt, Parameter = "A")

add_error_counts_tibble <- as_tibble(add_error_counts) |> 
  mutate(Model = "Additive", N = tau / actual_dt)

t_error_counts_tibble <- as_tibble(t_error_counts) |> 
  mutate(Model = "t-diffusion", N = tau / actual_dt)

both_error_counts_tibble <- as_tibble(both_error_count) |> 
  mutate(Model = "Both", N = tau / actual_dt)


error_count_tibble <- bind_rows(add_error_counts_tibble, t_error_counts_tibble, both_error_counts_tibble)

names(error_count_tibble) <- c(nus, "Model", "N")

error_count_tibble_long <- error_count_tibble |> 
  pivot_longer(-c(Model, N), names_to = "nu", values_to = "Error_count") |> 
  mutate(nu = factor(round(as.numeric(nu), 2)),
         Model = factor(Model, levels = c("Additive", "t-diffusion", "Both")),
  error_proportion = Error_count / (Error_count + numSim))

# if(!file.exists("data/error_count_nu_data.csv")){
#   utils::write.table(error_count_tibble_long, file="data/error_count_nu_data.csv",
#                      sep = ",", row.names = FALSE)
# } else{
#   error_count_tibble_long <- read_csv("data/error_count_nu_data.csv")
# }
error_count_tibble_long <- error_count_tibble_long |> 
mutate(nu = factor(round(as.numeric(nu), 2)),
       Model = factor(Model, levels = c("Additive", "t-diffusion", "Both")))

error_count_plot <- error_count_tibble_long |>
  ggplot(aes(x = nu, y = error_proportion, fill = Model)) +
  geom_col(position = "dodge2", col = "black", linewidth = 0.5) +
  facet_wrap(~factor(N)) +
  scale_fill_manual(values = thesis_palette) +
  scale_y_continuous(labels = scales::percent) + 
  labs(x = expression(nu*phantom(.)[sim]), y = "") +  
  theme(strip.text = element_text(face = "bold", size = 18), panel.spacing = unit(2, "lines"),
        axis.title.x = element_text(size = 22), axis.title.y = element_text(face = "bold"),
        legend.title = element_text(size = 20),
        legend.text  = element_text(size = 20), axis.text = element_text(size = 16)) +
  guides(fill = guide_legend(override.aes = list(size = 5, shape = NA, col = NA)))

# ggsave(error_count_plot, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "error_count_plot.jpeg",
#        height = 6, width = 13, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)

combined_nus_tibble <- bind_rows(add_dynamic_nus_tibble, t_dynamic_nus_tibble,
                                 add_dynamic_tau_tibble, t_dynamic_tau_tibble,
                                 add_dynamic_A_tibble, t_dynamic_A_tibble) |> 
  mutate(Model = factor(Model), Parameter = factor(Parameter))

names(combined_nus_tibble) <- c(nus, "Model", "N", "Parameter")

combined_nus_tibble_long <- combined_nus_tibble |> 
  pivot_longer(-c(Model, N, Parameter), names_to = "nu", values_to = "ARE") |> 
  mutate(nu = factor(round(as.numeric(nu), 2)))

# if(!file.exists("data/nu_estimation_ARE.csv")){
#   utils::write.table(combined_nus_tibble_long, file="data/nu_estimation_ARE.csv",
#                      sep = ",", row.names = FALSE)
# } else{
#   combined_nus_tibble_long <- read_csv("data/nu_estimation_ARE.csv")
# }


combined_nus_tibble_long <- combined_nus_tibble_long |> 
  mutate(nu = factor(round(as.numeric(nu), 2)))

combined_nus_plot <- combined_nus_tibble_long |> ggplot(aes(x = N, y = ARE, col = nu)) +
  geom_line(linewidth = 2.25) + geom_point(size = 3) + 
  facet_grid(Model ~ Parameter, 
             scales = "free_y", labeller = label_parsed) +
  scale_x_log10(breaks = tau / actual_dt) +
  scale_color_manual(values = thesis_palette) +
  labs(x = "N", y = "Absolute Relative Error", color = expression(nu*phantom(.)[sim])) +  
  theme(strip.text = element_text(face = "bold", size = 22),
        panel.spacing = unit(2.25, "lines"),
        axis.title.x = element_text(face = "bold", size = 20),
        axis.title.y = element_text(face = "bold", size = 20),
        legend.text  = element_text(size = 20),
        legend.title = element_text(size = 22),
        axis.text = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold", angle = 315)) +
  guides(color = guide_legend(override.aes = list(shape = NA, linewidth = 7.5)))
  
# ggsave(combined_nus_plot, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "combined_nus_plot.jpeg",
#        height = 6, width = 13, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)


# Early Warning signals
# Try t-distribution
true_param <- c(0.8, -2, -3, 0.15)
actual_dt <- 1 / 12
tau <- 175
t_0 <- 25
mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) * sqrt(abs(true_param[3] / true_param[1]))
alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
stationary_part_true_param <- c(alpha0, mu0, true_param[4])
dynamic_part_true_param <- c(tau, true_param[1])

log_steps <- seq(-1, 1, length.out = 9)

nu_values <- unique(round(sort(c(exp(log_steps), 1/exp(log_steps))), 2))
numSim <- 100

time_to_tipping <- seq(-125, -25, by = 25)

results_t <- list()

results_add <- list()

set.seed(110524)
for (nu in nu_values){
  cat("Currently at: ", nu)
  true_param_inner <- c(true_param, nu)
  mu0 <- true_param[2] + ifelse(true_param_inner[1] >= 0, 1, -1) * sqrt(abs(true_param_inner[3] / true_param_inner[1]))
  alpha0 <- 2 * sqrt(abs(true_param_inner[1] * true_param_inner[3]))
  stationary_part_true_param <- c(alpha0, mu0, true_param_inner[4])
  dynamic_part_true_param <- c(tau, true_param_inner[1], nu)
  
  accuracy_tau_t <- matrix(data = NA, nrow = numSim, ncol = length(time_to_tipping))
  accuracy_tau_add <- matrix(data = NA, nrow = numSim, ncol = length(time_to_tipping))
  
  for (i in seq_along(time_to_tipping)){
    print(time_to_tipping[i])
    for (j in 1:numSim){
      if(j%%5 == 0 | j == 1){print(j)}
     
    success <- FALSE
    while (!success) {
      tryCatch({
        random_seed <- sample.int(100000, size = 1)
        sim_res_t_distribution <- simulate_t_distribution_tipping_model(
          actual_dt, true_param_inner, tau,
          t_0, beyond_tipping = time_to_tipping[i],
          seed = random_seed)
        
        sim_res_additive <- simulate_additive_noise_tipping_model(
          actual_dt, true_param_inner, tau,
          t_0, beyond_tipping = time_to_tipping[i],
          seed = random_seed)
        
        # Stationary part
        
        t_dist_estim_param <- optimize_stationary_likelihood(
          likelihood_fun = t_diffusion_strang_splitting,
          data = asinh(sim_res_t_distribution$X_t[sim_res_t_distribution$t < t_0]),
          init_par = stationary_part_true_param,
          delta = actual_dt,
          exp_sigma = TRUE)$par
        
        additive_estim_param <- optimize_stationary_likelihood(
          likelihood_fun = OU_likelihood,
          data = sim_res_additive$X_t[sim_res_additive$t < t_0],
          init_par = stationary_part_true_param,
          delta = actual_dt,
          exp_sigma = TRUE)$par
        
        
        # Dynamic part
        noise_starting_values <- c(runif(n = 2, min = 0.9, 1.1), nu)
        accuracy_tau_t[j, i] <- (optimize_dynamic_likelihood(likelihood_fun = t_transform_dynamic_likelihood,
                       data = asinh(sim_res_t_distribution$X_t[sim_res_t_distribution$t > t_0]),
                       init_par = dynamic_part_true_param + noise_starting_values,
                       delta = actual_dt,
                       alpha0 = t_dist_estim_param[1],
                       mu0 = t_dist_estim_param[2],
                       sigma = t_dist_estim_param[3])$par[1] - tau) / tau

        accuracy_tau_add[j, i] <- (optimize_dynamic_likelihood(likelihood_fun = OU_dynamic_likelihood,
                       data  = sim_res_additive$X_t[sim_res_additive$t > t_0],
                       init_par = dynamic_part_true_param + noise_starting_values,
                       delta = actual_dt,
                       alpha0 = additive_estim_param[1],
                       mu0 = additive_estim_param[2],
                       sigma = additive_estim_param[3])$par[1] - tau) / tau
        
        
        success <- TRUE
      }, error = function(e) {
        cat("Error occurred at iteration", i, ":", conditionMessage(e), "\n")
      })
    }
    
    
    }
  }
  # Calculate means and quantiles
  column_wise_means_t <- apply(accuracy_tau_t, 2, function(x) quantile(x, probs = 0.5))
  column_wise_lower_quantiles_t <- apply(accuracy_tau_t, 2, function(x) quantile(x, probs = 0.1))
  column_wise_upper_quantiles_t <- apply(accuracy_tau_t, 2, function(x) quantile(x, probs = 0.9))
  
  column_wise_means_add <-apply(accuracy_tau_add, 2, function(x) quantile(x, probs = 0.5))
  column_wise_lower_quantiles_add <- apply(accuracy_tau_add, 2, function(x) quantile(x, probs = 0.1))
  column_wise_upper_quantiles_add <- apply(accuracy_tau_add, 2, function(x) quantile(x, probs = 0.9))
  
  # Store results
  results_t[[as.character(nu)]] <- tibble(
    x = time_to_tipping,
    median = column_wise_means_t,
    lower = column_wise_lower_quantiles_t,
    upper = column_wise_upper_quantiles_t,
    nu = nu
  )
  
  results_add[[as.character(nu)]] <- tibble(
    x = time_to_tipping,
    median = column_wise_means_add,
    lower = column_wise_lower_quantiles_add,
    upper = column_wise_upper_quantiles_add,
    nu = nu
  )
  
}

combined_data_t <- bind_rows(results_t, .id = "nu_label")

combined_data_t$nu <- as.factor(combined_data_t$nu)

combined_data_t <- combined_data_t |> mutate(model = "t-diffusion")

combined_data_add <- bind_rows(results_add, .id = "nu_label")

combined_data_add$nu <- as.factor(combined_data_add$nu)

combined_data_add <- combined_data_add |> mutate(model = "Additive")

combined_data <- bind_rows(combined_data_add, combined_data_t)

ggplot(combined_data, aes(x = x, y = median, fill = model)) +
  geom_line(linewidth = 1.25) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  facet_wrap(~nu, scales = "free_y") + 
  xlab("Observations until time") + ylab("Relative deviation from tipping time")


### Numerical Strang splitting
## F-diffusion based model
F_lamperti_drift <- function(y, par){
  beta  <- par[1]
  mu    <- par[2]
  sigma <- par[3]
  
  - 1 / sinh(y) * ((beta + sigma^2 / 2) * cosh(y) - beta * (2 * mu + 1))
}

linear_lamperti_drift <- function(y, par){
  beta <- par[1]
  mu <- par[2]
  sigma <- par[3]
  
  -(beta * (1 - mu * exp(-y)) + sigma^2 / 2)
}

# Choose parameters appropriate for any of the diffusions
true_param <- c(1.5, 0.4, -0.2, 0.15)
actual_dts <- c(1/5, 1/50, 1/150, 1/500, 1/1000)
tau <- 1
t_0 <- 30

# ## Stationary part
# # Parameters for stationary part
mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) * sqrt(abs(true_param[3] / true_param[1]))
alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))

stationary_part_true_param <- c(alpha0, mu0, true_param[4])

numSim <- 50
# F-diffusion setup

ARE_F_closed_alpha <- matrix(data = NA, nrow = numSim, ncol = length(actual_dts))
ARE_F_numeric_alpha <- matrix(data = NA, nrow = numSim, ncol = length(actual_dts))

ARE_F_closed_mu <- matrix(data = NA, nrow = numSim, ncol = length(actual_dts))
ARE_F_numeric_mu <- matrix(data = NA, nrow = numSim, ncol = length(actual_dts))

ARE_F_closed_sigma <- matrix(data = NA, nrow = numSim, ncol = length(actual_dts))
ARE_F_numeric_sigma <- matrix(data = NA, nrow = numSim, ncol = length(actual_dts))

col_wise_median_closed_ARE_F_closed_alpha <- numeric(length(actual_dts))
col_wise_median_closed_ARE_F_numeric_alpha <- numeric(length(actual_dts))
col_wise_median_closed_ARE_F_closed_mu <- numeric(length(actual_dts))
col_wise_median_closed_ARE_F_numeric_mu <- numeric(length(actual_dts))
col_wise_median_closed_ARE_F_closed_sigma <- numeric(length(actual_dts))
col_wise_median_closed_ARE_F_numeric_sigma <- numeric(length(actual_dts))

col_wise_median_closed_F <- numeric(length(actual_dts))
col_wise_median_numeric_F <- numeric(length(actual_dts))

col_wise_upper_closed_F <- numeric(length(actual_dts))
col_wise_upper_numeric_F <- numeric(length(actual_dts))

col_wise_lower_closed_F <- numeric(length(actual_dts))
col_wise_lower_numeric_F <- numeric(length(actual_dts))

closed_form_dist_F <- matrix(data = NA, nrow = numSim, ncol = length(actual_dts))
numeric_dist_F <- matrix(data = NA, nrow = numSim, ncol = length(actual_dts))

# Linear setup
ARE_Linear_closed_alpha <- matrix(data = NA, nrow = numSim, ncol = length(actual_dts))
ARE_Linear_numeric_alpha <- matrix(data = NA, nrow = numSim, ncol = length(actual_dts))

ARE_Linear_closed_mu <- matrix(data = NA, nrow = numSim, ncol = length(actual_dts))
ARE_Linear_numeric_mu <- matrix(data = NA, nrow = numSim, ncol = length(actual_dts))

ARE_Linear_closed_sigma <- matrix(data = NA, nrow = numSim, ncol = length(actual_dts))
ARE_Linear_numeric_sigma <- matrix(data = NA, nrow = numSim, ncol = length(actual_dts))

col_wise_median_closed_ARE_Linear_closed_alpha <- numeric(length(actual_dts))
col_wise_median_closed_ARE_Linear_numeric_alpha <- numeric(length(actual_dts))
col_wise_median_closed_ARE_Linear_closed_mu <- numeric(length(actual_dts))
col_wise_median_closed_ARE_Linear_numeric_mu <- numeric(length(actual_dts))
col_wise_median_closed_ARE_Linear_closed_sigma <- numeric(length(actual_dts))
col_wise_median_closed_ARE_Linear_numeric_sigma <- numeric(length(actual_dts))

col_wise_median_closed_Linear <- numeric(length(actual_dts))
col_wise_median_numeric_Linear <- numeric(length(actual_dts))

col_wise_upper_closed_Linear <- numeric(length(actual_dts))
col_wise_upper_numeric_Linear <- numeric(length(actual_dts))

col_wise_lower_closed_Linear <- numeric(length(actual_dts))
col_wise_lower_numeric_Linear <- numeric(length(actual_dts))

closed_form_dist_Linear <- matrix(data = NA, nrow = numSim, ncol = length(actual_dts))
numeric_dist_Linear <- matrix(data = NA, nrow = numSim, ncol = length(actual_dts))

# Count number of fails
F_closed_fails <-  matrix(data = 0, nrow = numSim, ncol = length(actual_dts))
F_numeric_fails <-matrix(data = 0, nrow = numSim, ncol = length(actual_dts))
Linear_closed_fails <-matrix(data = 0, nrow = numSim, ncol = length(actual_dts))
Linear_numeric_fails <-matrix(data = 0, nrow = numSim, ncol = length(actual_dts))

# Count number of total function and gradient evaluations
F_closed_function_count <-  matrix(data = 0, nrow = numSim, ncol = length(actual_dts))
F_closed_gradient_count <-  matrix(data = 0, nrow = numSim, ncol = length(actual_dts))
F_numeric_function_count <-  matrix(data = 0, nrow = numSim, ncol = length(actual_dts))
F_numeric_gradient_count <-  matrix(data = 0, nrow = numSim, ncol = length(actual_dts))
Linear_closed_function_count <-  matrix(data = 0, nrow = numSim, ncol = length(actual_dts))
Linear_closed_gradient_count <-  matrix(data = 0, nrow = numSim, ncol = length(actual_dts))
Linear_numeric_function_count <-  matrix(data = 0, nrow = numSim, ncol = length(actual_dts))
Linear_numeric_gradient_count <-  matrix(data = 0, nrow = numSim, ncol = length(actual_dts))


set.seed(14052024)
for (j in seq_along(actual_dts)){
cat("Grinding delta number: ", j, "\n")
for (i in 1:numSim){
  if(i %% 5 == 1){
  cat("At simulation number: ", i, "\n")
  }
  success <- FALSE
  while(!success){
  F_closed_fail <- FALSE
  F_numeric_fail <- FALSE
  Linear_closed_fail <- FALSE
  Linear_numeric_fail <- FALSE
  simSeed <- sample.int(100000, 1)
  F_sim_dynamic <- simulate_F_distribution_tipping_model(actual_dts[j], true_param, t_0 = t_0, 
                                                         tau = tau, seed = simSeed)
  
  sim_res_linear <- simulate_linear_noise_tipping_model(actual_dts[j], true_param,
                                                        tau, t_0, seed = simSeed)
  
  random_noise_numeric_test <- runif(3, min = 0.75, max =1.25)
  tryCatch({
  F_closed_res <- optimize_stationary_likelihood(
                   likelihood_fun = F_diffusion_strang_splitting, 
                   init_par = stationary_part_true_param * random_noise_numeric_test, 
                   data = 2 * asinh(sqrt(F_sim_dynamic$X_t[F_sim_dynamic$t<t_0])),
                   delta = actual_dts[j], exp_sigma = TRUE, 
                   return_all = TRUE)
  }, error = function(e){
    F_closed_fails[i, j] <<- F_closed_fails[i, j] + 1
    F_closed_fail <<- TRUE
  }
  )
  ARE_F_closed_alpha[i, j] <- abs(F_closed_res$par[1] - stationary_part_true_param[1]) /
    stationary_part_true_param[1]
  ARE_F_closed_mu[i, j] <- abs(F_closed_res$par[2] - stationary_part_true_param[2]) /
    stationary_part_true_param[2]
  ARE_F_closed_sigma[i, j] <- abs(F_closed_res$par[3] - stationary_part_true_param[3]) /
    stationary_part_true_param[3]
  
  F_closed_function_count[i, j] <-  unname(F_closed_res$counts[1])
  F_closed_gradient_count[i, j] <-  unname(F_closed_res$counts[2])

  
  tryCatch({
  F_numeric_res <- optimize_stationary_likelihood(
                   likelihood_fun = numerical_strang_splitting,
                   init_par = stationary_part_true_param * random_noise_numeric_test, 
                   data = 2 * asinh(sqrt(F_sim_dynamic$X_t[F_sim_dynamic$t<t_0])),
                   delta = actual_dts[j],
                   drift_lamperti_sde = F_lamperti_drift,
                   exp_sigma = FALSE, return_all = TRUE)
  }, error = function(e){
    F_numeric_fails[i, j] <<- F_numeric_fails[i, j] + 1
    F_numeric_fail <<- TRUE
  }
  )
  
  ARE_F_numeric_alpha[i, j] <- abs(F_numeric_res$par[1] - stationary_part_true_param[1]) /
    stationary_part_true_param[1]
  ARE_F_numeric_mu[i, j] <- abs(F_numeric_res$par[2] - stationary_part_true_param[2]) /
    stationary_part_true_param[2]
  ARE_F_numeric_sigma[i, j] <- abs(F_numeric_res$par[3] - stationary_part_true_param[3]) /
    stationary_part_true_param[3]
  
  F_numeric_function_count[i, j] <-  unname(F_numeric_res$counts[1])
  F_numeric_gradient_count[i, j] <-  unname(F_numeric_res$counts[2])
  
  tryCatch({
  Linear_closed_res <-  optimize_stationary_likelihood(
                                   likelihood_fun =  mean_reverting_GBM_strang,
                                   data = log(sim_res_linear$X_t[sim_res_linear$t<t_0]),
                                   init_par = stationary_part_true_param * random_noise_numeric_test,
                                   delta = actual_dts[j],
                                   exp_sigma = FALSE, return_all = TRUE)
  }, error = function(e){
    Linear_closed_fails[i, j] <<- Linear_closed_fails[i, j] + 1
    Linear_closed_fail <<- TRUE
  }
  )
  
  ARE_Linear_closed_alpha[i, j] <- abs(Linear_closed_res$par[1] - stationary_part_true_param[1]) /
    stationary_part_true_param[1]
  ARE_Linear_closed_mu[i, j] <- abs(Linear_closed_res$par[2] - stationary_part_true_param[2]) /
    stationary_part_true_param[2]
  ARE_Linear_closed_sigma[i, j] <- abs(Linear_closed_res$par[3] - stationary_part_true_param[3]) /
    stationary_part_true_param[3]
  
  Linear_closed_function_count[i, j] <-  unname(Linear_closed_res$counts[1])
  Linear_closed_gradient_count[i, j] <-  unname(Linear_closed_res$counts[2])
  
  tryCatch({
  Linear_numeric_res <- optimize_stationary_likelihood(
                       likelihood_fun = numerical_strang_splitting,
                       data = log(sim_res_linear$X_t[sim_res_linear$t<t_0]),
                       init_par = stationary_part_true_param * random_noise_numeric_test,
                       delta = actual_dts[j],
                       exp_sigma = FALSE, return_all = TRUE,
                       drift_lamperti_sde = linear_lamperti_drift)
  }, error = function(e){
    Linear_numeric_fails[i, j] <<- Linear_numeric_fails[i, j] + 1
    Linear_numeric_fail <<- TRUE
  }
  )
  
  ARE_Linear_numeric_alpha[i, j] <- abs(Linear_numeric_res$par[1] - stationary_part_true_param[1]) /
    stationary_part_true_param[1]
  ARE_Linear_numeric_mu[i, j] <- abs(Linear_numeric_res$par[2] - stationary_part_true_param[2]) /
    stationary_part_true_param[2]
  ARE_Linear_numeric_sigma[i, j] <- abs(Linear_numeric_res$par[3] - stationary_part_true_param[3]) /
    stationary_part_true_param[3]
  
  Linear_numeric_function_count[i, j] <-  unname(Linear_numeric_res$counts[1])
  Linear_numeric_gradient_count[i, j] <-  unname(Linear_numeric_res$counts[2])
  
  if(!F_closed_fail && !F_numeric_fail && !Linear_closed_fail && !Linear_numeric_fail){
    success <- TRUE
  }
  }
  # We have ensured everthing goes right before benchmarking the times
  
  closed_form_dist_F[i, j] <- microbenchmark::microbenchmark(
         optimize_stationary_likelihood(likelihood_fun = F_diffusion_strang_splitting, 
         init_par = stationary_part_true_param * random_noise_numeric_test, 
         data = 2 * asinh(sqrt(F_sim_dynamic$X_t[F_sim_dynamic$t<t_0])),
         delta = actual_dts[j], exp_sigma = TRUE), times = 1, unit = "us")$time / 1e9
  
  numeric_dist_F[i, j] <- microbenchmark::microbenchmark(
    optimize_stationary_likelihood(likelihood_fun = numerical_strang_splitting,
     init_par = stationary_part_true_param * random_noise_numeric_test, 
     data = 2 * asinh(sqrt(F_sim_dynamic$X_t[F_sim_dynamic$t<t_0])),
     delta = actual_dts[j],
     drift_lamperti_sde = F_lamperti_drift,
     exp_sigma = FALSE), times = 1, unit = "us")$time / 1e9
  
  closed_form_dist_Linear[i, j] <- microbenchmark::microbenchmark(
    optimize_stationary_likelihood(
      likelihood_fun =  mean_reverting_GBM_strang,
      data = log(sim_res_linear$X_t[sim_res_linear$t<t_0]),
      init_par = stationary_part_true_param * random_noise_numeric_test,
      delta = actual_dts[j],
      exp_sigma = FALSE), times = 1, unit = "us")$time / 1e9
  
  numeric_dist_Linear[i, j] <- microbenchmark::microbenchmark(
    optimize_stationary_likelihood(
      likelihood_fun = numerical_strang_splitting,
      data = log(sim_res_linear$X_t[sim_res_linear$t<t_0]),
      init_par = stationary_part_true_param * random_noise_numeric_test,
      delta = actual_dts[j],
      exp_sigma = FALSE,
      drift_lamperti_sde = linear_lamperti_drift), times = 1, unit = "us")$time / 1e9
  
}
}
# Fail analysis
numeric_fail_tibble <- tibble(F_closed_fails = colSums(F_closed_fails),
       F_numeric_fails = colSums(F_numeric_fails), 
       Linear_closed_fails = colSums(Linear_closed_fails),
       Linear_numeric_fails = colSums(Linear_numeric_fails),
       N = t_0 / actual_dts)

# if(!file.exists("data/count_numeric_fails.csv")){
#   utils::write.table(numeric_fail_tibble, file="data/count_numeric_fails.csv",
#                      sep = ",", row.names = FALSE)
# } else{
#   numeric_fail_tibble <- read_csv("data/count_numeric_fails.csv")
# }



# Evaluate the ARE
col_wise_median_closed_ARE_F_closed_alpha <-  colwise_median(ARE_F_closed_alpha)
col_wise_median_closed_ARE_F_numeric_alpha <-  colwise_median(ARE_F_numeric_alpha)
col_wise_median_closed_ARE_F_closed_mu <-  colwise_median(ARE_F_closed_mu)
col_wise_median_closed_ARE_F_numeric_mu <-  colwise_median(ARE_F_numeric_mu)
col_wise_median_closed_ARE_F_closed_sigma <-  colwise_median(ARE_F_closed_sigma)
col_wise_median_closed_ARE_F_numeric_sigma <-  colwise_median(ARE_F_numeric_sigma)
col_wise_median_closed_ARE_Linear_closed_alpha <-  colwise_median(ARE_Linear_closed_alpha)
col_wise_median_closed_ARE_Linear_numeric_alpha <-  colwise_median(ARE_Linear_numeric_alpha)
col_wise_median_closed_ARE_Linear_closed_mu <-  colwise_median(ARE_Linear_closed_mu)
col_wise_median_closed_ARE_Linear_numeric_mu <-  colwise_median(ARE_Linear_numeric_mu)
col_wise_median_closed_ARE_Linear_closed_sigma <-  colwise_median(ARE_Linear_closed_sigma)
col_wise_median_closed_ARE_Linear_numeric_sigma <-  colwise_median(ARE_Linear_numeric_sigma)

ARE_median_result <- bind_rows(
  as_tibble(col_wise_median_closed_ARE_F_closed_alpha) |>
    mutate(Model = "F-diffusion", Parameter = "alpha", Type = "Closed", N = t_0 / actual_dts),
  
  as_tibble(col_wise_median_closed_ARE_F_numeric_alpha) |>
    mutate(Model = "F-diffusion", Parameter = "alpha", Type = "Numeric", N = t_0 / actual_dts),
  
  as_tibble(col_wise_median_closed_ARE_F_closed_mu) |>
    mutate(Model = "F-diffusion", Parameter = "mu", Type = "Closed", N = t_0 / actual_dts),
  
  as_tibble(col_wise_median_closed_ARE_F_numeric_mu) |>
    mutate(Model = "F-diffusion", Parameter = "mu", Type = "Numeric", N = t_0 / actual_dts),
  
  as_tibble(col_wise_median_closed_ARE_F_closed_sigma) |>
    mutate(Model = "F-diffusion", Parameter = "sigma", Type = "Closed", N = t_0 / actual_dts),
  
  as_tibble(col_wise_median_closed_ARE_F_numeric_sigma) |>
    mutate(Model = "F-diffusion", Parameter = "sigma", Type = "Numeric", N = t_0 / actual_dts),
  
  as_tibble(col_wise_median_closed_ARE_Linear_closed_alpha) |>
    mutate(Model = "Linear", Parameter = "alpha", Type = "Closed", N = t_0 / actual_dts),
  
  as_tibble(col_wise_median_closed_ARE_Linear_numeric_alpha) |>
    mutate(Model = "Linear", Parameter = "alpha", Type = "Numeric", N = t_0 / actual_dts),
  
  as_tibble(col_wise_median_closed_ARE_Linear_closed_mu) |>
    mutate(Model = "Linear", Parameter = "mu", Type = "Closed", N = t_0 / actual_dts),
  
  as_tibble(col_wise_median_closed_ARE_Linear_numeric_mu) |>
    mutate(Model = "Linear", Parameter = "mu", Type = "Numeric", N = t_0 / actual_dts),
  
  as_tibble(col_wise_median_closed_ARE_Linear_closed_sigma) |>
    mutate(Model = "Linear", Parameter = "sigma", Type = "Closed", N = t_0 / actual_dts),
  
  as_tibble(col_wise_median_closed_ARE_Linear_numeric_sigma) |>
    mutate(Model = "Linear", Parameter = "sigma", Type = "Numeric", N = t_0 / actual_dts)
)|> mutate(Model = factor(Model), Parameter = factor(Parameter), Type = factor(Type))

# if(!file.exists("data/ARE_numerical_vs_closed_form_result_median.csv")){
#   utils::write.table(ARE_median_result, file="data/ARE_numerical_vs_closed_form_result_median.csv",
#                      sep = ",", row.names = FALSE)
# } else{
#   ARE_median_result <- read_csv("data/ARE_numerical_vs_closed_form_result_median.csv")
# }


ARE_median_result |> ggplot(aes(x = N, y = value, col = Type)) +
  geom_point(size = 3) + geom_line(linewidth = 1.25) + 
  scale_x_log10(breaks = t_0 / actual_dts) + scale_y_log10() + 
  facet_grid(Model~Parameter, scales = "free_y") + scale_color_manual(values = thesis_palette)

# Distribution

ARE_F_closed_alpha_tibble <- as_tibble(ARE_F_closed_alpha) |>
  mutate(Model = "F-diffusion", Parameter = "alpha", Type = "Closed")
names(ARE_F_closed_alpha_tibble) <- c(t_0 / actual_dts, "Model", "Parameter", "Type")

ARE_F_numeric_alpha_tibble <- as_tibble(ARE_F_numeric_alpha) |>
  mutate(Model = "F-diffusion", Parameter = "alpha", Type = "Numeric")
names(ARE_F_numeric_alpha_tibble) <- c(t_0 / actual_dts, "Model", "Parameter", "Type")

ARE_F_closed_mu_tibble <- as_tibble(ARE_F_closed_mu) |>
  mutate(Model = "F-diffusion", Parameter = "mu", Type = "Closed")
names(ARE_F_closed_mu_tibble) <- c(t_0 / actual_dts, "Model", "Parameter", "Type")

ARE_F_numeric_mu_tibble <- as_tibble(ARE_F_numeric_mu) |>
  mutate(Model = "F-diffusion", Parameter = "mu", Type = "Numeric")
names(ARE_F_numeric_mu_tibble) <- c(t_0 / actual_dts, "Model", "Parameter", "Type")

ARE_F_closed_sigma_tibble <- as_tibble(ARE_F_closed_sigma) |>
  mutate(Model = "F-diffusion", Parameter = "sigma", Type = "Closed")
names(ARE_F_closed_sigma_tibble) <- c(t_0 / actual_dts, "Model", "Parameter", "Type")

ARE_F_numeric_sigma_tibble <- as_tibble(ARE_F_numeric_sigma) |>
  mutate(Model = "F-diffusion", Parameter = "sigma", Type = "Numeric")
names(ARE_F_numeric_sigma_tibble) <- c(t_0 / actual_dts, "Model", "Parameter", "Type")

ARE_Linear_closed_alpha_tibble <- as_tibble(ARE_Linear_closed_alpha) |>
  mutate(Model = "Linear", Parameter = "alpha", Type = "Closed")
names(ARE_Linear_closed_alpha_tibble) <- c(t_0 / actual_dts, "Model", "Parameter", "Type")

ARE_Linear_numeric_alpha_tibble <-  as_tibble(ARE_Linear_numeric_alpha) |>
  mutate(Model = "Linear", Parameter = "alpha", Type = "Numeric")
names(ARE_Linear_numeric_alpha_tibble) <- c(t_0 / actual_dts, "Model", "Parameter", "Type")

ARE_Linear_closed_mu_tibble <- as_tibble(ARE_Linear_closed_mu) |>
  mutate(Model = "Linear", Parameter = "mu", Type = "Closed")
names(ARE_Linear_closed_mu_tibble) <- c(t_0 / actual_dts, "Model", "Parameter", "Type")

ARE_Linear_numeric_mu_tibble <- as_tibble(ARE_Linear_numeric_mu) |>
  mutate(Model = "Linear", Parameter = "mu", Type = "Numeric")
names(ARE_Linear_numeric_mu_tibble) <- c(t_0 / actual_dts, "Model", "Parameter", "Type")

ARE_Linear_closed_sigma_tibble <- as_tibble(ARE_Linear_closed_sigma) |>
  mutate(Model = "Linear", Parameter = "sigma", Type = "Closed")
names(ARE_Linear_closed_sigma_tibble) <- c(t_0 / actual_dts, "Model", "Parameter", "Type")

ARE_Linear_numeric_sigma_tibble <- as_tibble(ARE_Linear_numeric_sigma) |>
  mutate(Model = "Linear", Parameter = "sigma", Type = "Numeric")
names(ARE_Linear_numeric_sigma_tibble) <- c(t_0 / actual_dts, "Model", "Parameter", "Type")


ARE_dist_result <- bind_rows(ARE_F_closed_alpha_tibble, ARE_F_numeric_alpha_tibble,
          ARE_F_closed_mu_tibble, ARE_F_numeric_mu_tibble, 
          ARE_F_closed_sigma_tibble, ARE_F_numeric_sigma_tibble, 
          ARE_Linear_closed_alpha_tibble, ARE_Linear_numeric_alpha_tibble,
          ARE_Linear_closed_mu_tibble, ARE_Linear_numeric_mu_tibble, 
          ARE_Linear_closed_sigma_tibble, ARE_Linear_numeric_sigma_tibble) |> 
  pivot_longer(cols = -c(Model, Parameter, Type), names_to = "N", values_to = "ARE") |> 
  mutate(Model = factor(Model), Parameter = factor(Parameter), Type = factor(Type), N = as.numeric(N))

# if(!file.exists("data/ARE_numerical_vs_closed_form_result.csv")){
#   utils::write.table(ARE_dist_result, file="data/ARE_numerical_vs_closed_form_result.csv",
#                      sep = ",", row.names = FALSE)
# } else{
#   ARE_dist_result <- read_csv("data/ARE_numerical_vs_closed_form_result.csv")
# }

ARE_dist_result_plot_Linear <- ARE_dist_result |> filter(Model == "Linear") |> 
  mutate(Type = ifelse(Type == "Numeric", "Numerical", Type)) |> 
  ggplot(aes(x = Type, y = ARE, fill = Parameter)) +
  geom_violin(scale = "width") + 
  facet_wrap(~factor(N), ncol = 5) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold", size = 18),
        legend.text  = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.title = element_text(size = 18),
        axis.text = element_text(face = "bold", size = 13)) +
  scale_fill_manual(values = thesis_palette,
                    labels = expression(alpha*phantom(.)[0], mu*phantom(.)[0], sigma)) +
  xlab("")

# ggsave(ARE_dist_result_plot_Linear, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "ARE_dist_result_plot_Linear.jpeg",
#        height = 6, width = 13, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)

ARE_dist_result_plot_F <- ARE_dist_result |> filter(Model == "F-diffusion") |> 
  mutate(Type = ifelse(Type == "Numeric", "Numerical", Type)) |> 
  ggplot(aes(x = Type, y = ARE, fill = Parameter)) +
  geom_violin(scale = "width") + 
  facet_wrap(~factor(N), ncol = 5) + 
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold", size = 18),
        legend.text  = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.title = element_text(size = 18),
        axis.text = element_text(face = "bold", size = 13)) +
  scale_fill_manual(values = thesis_palette,
                    labels = expression(alpha*phantom(.)[0], mu*phantom(.)[0], sigma)) +
  xlab("")

# ggsave(ARE_dist_result_plot_F, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "ARE_dist_result_plot_F.jpeg",
#        height = 6, width = 13, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)


# Evaluate the running times
col_wise_median_closed_F <- colwise_median(closed_form_dist_F)
col_wise_median_numeric_F <- colwise_median(numeric_dist_F)
col_wise_median_closed_Linear <- colwise_median(closed_form_dist_Linear)
col_wise_median_numeric_Linear <- colwise_median(numeric_dist_Linear)
# Median tibbles
median_closed_F_tibble <- as_tibble(col_wise_median_closed_F) |>
  mutate(Type = "Closed form", Model = "F-diffusion", N = t_0 / actual_dts, quantile = "Median")
median_numeric_F_tibble <- as_tibble(col_wise_median_numeric_F) |>
  mutate(Type = "Numerical", Model = "F-diffusion", N = t_0 / actual_dts, quantile = "Median")

median_closed_Linear_tibble <- as_tibble(col_wise_median_closed_Linear) |>
  mutate(Type = "Closed form", Model = "Linear", N = t_0 / actual_dts, quantile = "Median")
median_numeric_Linear_tibble <- as_tibble(col_wise_median_numeric_Linear) |>
  mutate(Type = "Numerical", Model = "Linear", N = t_0 / actual_dts, quantile = "Median")

# Upper quantiles
col_wise_upper_closed_F <- colwise_quantile(closed_form_dist_F, probs = 0.9)
col_wise_upper_numeric_F <- colwise_quantile(numeric_dist_F, probs = 0.9)
col_wise_upper_closed_Linear <- colwise_quantile(closed_form_dist_Linear, probs = 0.9)
col_wise_upper_numeric_Linear <- colwise_quantile(numeric_dist_Linear, probs = 0.9)

upper_closed_F_tibble <- as_tibble(col_wise_upper_closed_F) |>
  mutate(Type = "Closed form", Model = "F-diffusion", N = t_0 / actual_dts, quantile = "Upper")
upper_numeric_F_tibble <- as_tibble(col_wise_upper_numeric_F) |>
  mutate(Type = "Numerical", Model = "F-diffusion", N = t_0 / actual_dts, quantile = "Upper")

upper_closed_Linear_tibble <- as_tibble(col_wise_upper_closed_Linear) |>
  mutate(Type = "Closed form", Model = "Linear", N = t_0 / actual_dts, quantile = "Upper")
upper_numeric_Linear_tibble <- as_tibble(col_wise_upper_numeric_Linear) |>
  mutate(Type = "Numerical", Model = "Linear", N = t_0 / actual_dts, quantile = "Upper")

# Lower quantiles
col_wise_lower_closed_F <- colwise_quantile(closed_form_dist_F, probs = 0.1)
col_wise_lower_numeric_F <- colwise_quantile(numeric_dist_F, probs = 0.1)
col_wise_lower_closed_Linear <- colwise_quantile(closed_form_dist_Linear, probs = 0.1)
col_wise_lower_numeric_Linear <- colwise_quantile(numeric_dist_Linear, probs = 0.1)

lower_closed_F_tibble <- as_tibble(col_wise_lower_closed_F) |>
  mutate(Type = "Closed form", Model = "F-diffusion", N = t_0 / actual_dts, quantile = "Lower")
lower_numeric_F_tibble <- as_tibble(col_wise_lower_numeric_F) |>
  mutate(Type = "Numerical", Model = "F-diffusion", N = t_0 / actual_dts, quantile = "Lower")

lower_closed_Linear_tibble <- as_tibble(col_wise_lower_closed_Linear) |>
  mutate(Type = "Closed form", Model = "Linear", N = t_0 / actual_dts, quantile = "Lower")
lower_numeric_Linear_tibble <- as_tibble(col_wise_lower_numeric_Linear) |>
  mutate(Type = "Numerical", Model = "Linear", N = t_0 / actual_dts, quantile = "Lower")

# Combine all results into a single tibble
Running_result <- bind_rows(
  median_closed_F_tibble,
  median_numeric_F_tibble,
  median_closed_Linear_tibble,
  median_numeric_Linear_tibble,
  upper_closed_F_tibble,
  upper_numeric_F_tibble,
  upper_closed_Linear_tibble,
  upper_numeric_Linear_tibble,
  lower_closed_F_tibble,
  lower_numeric_F_tibble,
  lower_closed_Linear_tibble,
  lower_numeric_Linear_tibble
) |> pivot_wider(values_from = value, names_from = quantile)

# if(!file.exists("data/Running_result_numeric.csv")){
#   utils::write.table(Running_result, file="data/Running_result_numeric.csv",
#                      sep = ",", row.names = FALSE)
# } else{
#   Running_result <- read_csv("data/Running_result_numeric.csv")
# }

Running_result_numeric_plot <- Running_result |> 
  ggplot(aes(x = N, y = Median, col = Type, fill = Type)) +
  geom_line(linewidth = 2) + geom_point(size = 4) + 
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15) + 
  facet_wrap(~Model) +
  scale_x_log10(breaks = t_0 / actual_dts) + scale_y_log10() + 
  scale_fill_manual(values = thesis_palette) + scale_color_manual(values = thesis_palette) +
  theme(strip.text = element_text(face = "bold", size = 20),
        panel.spacing = unit(3, "lines"),
        axis.title.x = element_text(face = "bold", size = 18),
        axis.title.y = element_text(face = "bold", size = 18),
        legend.title = element_text(face = "bold", size = 20),
        legend.text  = element_text(size = 22), axis.text = element_text(size = 15),
        axis.text.y = element_text(size = 18)) + 
  ylab("Running time (s)")

# ggsave(Running_result_numeric_plot, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "Running_result_numeric.jpeg",
#        height = 6, width = 16, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)


# Number of iterations
function_gradient_count_result <-  bind_rows(
as_tibble(F_closed_function_count) |> 
  mutate(Model = "F-diffusion", Type = "Function", Method =  "Closed form"),
as_tibble(F_closed_gradient_count) |> 
  mutate(Model = "F-diffusion", Type = "Gradient", Method =  "Closed form"),
as_tibble(F_numeric_function_count) |> 
  mutate(Model = "F-diffusion", Type = "Function", Method =  "Numerical"),
as_tibble(F_numeric_gradient_count) |> 
  mutate(Model = "F-diffusion", Type = "Gradient", Method =  "Numerical"),
as_tibble(Linear_closed_function_count) |> 
  mutate(Model = "Linear", Type = "Function", Method =  "Closed form"),
as_tibble(Linear_closed_gradient_count) |> 
  mutate(Model = "Linear", Type = "Gradient", Method =  "Closed form"),
as_tibble(Linear_numeric_function_count) |> 
  mutate(Model = "Linear", Type = "Function", Method =  "Numerical"),
as_tibble(Linear_numeric_gradient_count) |> 
  mutate(Model = "Linear", Type = "Gradient", Method =  "Numerical")
) |> mutate(Model = factor(Model), Type = factor(Type), Method = factor(Method))
names(function_gradient_count_result) <- c(t_0 / actual_dts, "Model", "Type", "Method")

function_gradient_count_result_long <- function_gradient_count_result |>
  pivot_longer(-c(Model, Type, Method), values_to = "Count", names_to = "N") |> 
  mutate(N = as.numeric(N))


# if(!file.exists("data/gradientAndFunctionCount.csv")){
#   utils::write.table(function_gradient_count_result_long, file="data/gradientAndFunctionCount.csv",
#                      sep = ",", row.names = FALSE)
# } else{
#   function_gradient_count_result_long <- read_csv("data/gradientAndFunctionCount.csv")
# }

function_gradient_count_Linear_plot <- function_gradient_count_result_long |>
  filter(Model == "Linear") |> 
  mutate(Method = ifelse(Method == "Numeric", "Numerical", Method)) |> 
  ggplot(aes(x = Method, y = Count, fill = Type)) +
  geom_boxplot(alpha = 0.85) + scale_y_log10() +
  facet_wrap(~factor(N), ncol = 5) +
  scale_fill_manual(values = thesis_palette) +
  xlab("") +
  ylab("Number of Evaluations") +
  theme(
    strip.text = element_text(size = 22, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 22),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.key.size = unit(2, 'lines')
  )

# ggsave(function_gradient_count_Linear_plot, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "function_gradient_count_Linear_plot.jpeg",
#        height = 8, width = 16, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)

function_gradient_count_F_plot <- function_gradient_count_result_long |>
  filter(Model == "F-diffusion") |>
  mutate(Method = ifelse(Method == "Numeric", "Numerical", Method)) |>  
  ggplot(aes(x = Method, y = Count, fill = Type)) +
  geom_boxplot(alpha = 0.85) + scale_y_log10() +
  facet_wrap(~factor(N), ncol = 5) +
  scale_fill_manual(values = thesis_palette) +
  xlab("") +
  ylab("Number of Evaluations") +
  theme(
    strip.text = element_text(size = 22, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 22),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.key.size = unit(2, 'lines')
  )

# ggsave(function_gradient_count_F_plot, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "function_gradient_count_F_plot.jpeg",
#        height = 6, width = 16, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)

summarized_function_gradient_count_result <- function_gradient_count_result_long |>
  group_by(Model, Type, Method, N) |> 
  summarise(Median = median(Count), upper = quantile(Count, 0.9), lower = quantile(Count, 0.1)) |> 
  ungroup()

summarized_function_gradient_count_result |>
  ggplot(aes(x = N, y = Median, col = Method, fill = Method)) +
  geom_line(linewidth = 1.25) + geom_point(size = 2.5) + 
  scale_x_log10(breaks = t_0 / actual_dts) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .2) + 
  facet_grid(Type~Model, scales = "free_y") + 
  scale_fill_manual(values = thesis_palette) + 
  scale_color_manual(values = thesis_palette) +
  ylab("Evaluations") + xlab("N")
  

# Model misspecification
tau_values <- c(35, 60, 90, 132, 175, 200, 250)

# Initialize lists to store results
results_list <- list(
  quantiles_OU = list(),
  quantiles_OU_dynamic = list(),
  quantiles_t = list(),
  quantiles_t_dynamic = list(),
  stationary_estim_param = list(),
  dynamic_estim_param = list()
)

# Define the true parameters
true_param <- c(0.87, -1.51, -2.69, sqrt(0.30))
mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) * sqrt(abs(true_param[3] / true_param[1]))
alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
stationary_part_true_param <- c(alpha0, mu0 , true_param[4])
actual_dt <- 1/12
t_0 <- 54

# Parameters for the simulation
numSim <- 200
considered_quantiles <- c(0.01, 0.025, 0.05, 0.1,
                          0.175, 0.25, 0.33, 
                          0.5, 0.67, 0.75, 0.825,
                          0.9, 0.95, 0.975, 0.99)

set.seed(150524)

for (tau in tau_values) {
  cat("Processing tau =", tau, "\n")
  
  # Initialize matrices for each value of tau
  quantiles_OU <- matrix(NA, nrow = numSim, ncol = length(considered_quantiles))
  quantiles_OU_dynamic <- matrix(NA, nrow = numSim, ncol = length(considered_quantiles))
  quantiles_t <- matrix(NA, nrow = numSim, ncol = length(considered_quantiles))
  quantiles_t_dynamic <- matrix(NA, nrow = numSim, ncol = length(considered_quantiles))
  
  # Initialize matrices for parameter estimates
  stationary_estim_param <- matrix(NA, nrow = numSim, ncol = 6)
  dynamic_estim_param <- matrix(NA, nrow = numSim, ncol = 4)
  
  for (i in 1:numSim) {
    if (i == 1 || i %% 5 == 0) {
      cat("Grinding ", i, "\n")
    }
    success <- FALSE
    while (!success) {
      seedOU <- sample.int(100000, size = 1)
      sim_res_add <- simulate_additive_noise_tipping_model(actual_dt, true_param, tau, t_0, seed = seedOU)
      
      mu_init <- mean(sim_res_add$X_t[sim_res_add$t <= t_0])
      
      rho_init <- sum((sim_res_add$X_t[sim_res_add$t <= t_0 & sim_res_add$t > 0] - mu_init) *
                        (sim_res_add$X_t[sim_res_add$t < t_0] - mu_init)) / sum((sim_res_add$X_t[sim_res_add$t <= t_0] - mu_init)^2)
      
      alpha_init <- -log(rho_init) / actual_dt
      
      gamma_sq_init <- mean((sim_res_add$X_t[sim_res_add$t <= t_0 & sim_res_add$t > 0] - 
                               sim_res_add$X_t[sim_res_add$t < t_0] * rho_init - mu_init * (1 - rho_init))^2) / 
        (1 - rho_init^2)
      
      sigma_init <- sqrt(gamma_sq_init * 2 * alpha_init)
      
      initital_param_stationary <- c(alpha_init, mu_init, sigma_init)
      
      noise_dynamic_part <- runif(2, min = .85, max = 1.15)
      failOU <- FALSE
      failOU_dynamic <- FALSE
      failt <- FALSE
      failt_dynamic <- FALSE
      fail_quantiles <- FALSE
      
      tryCatch({
        stationary_estim_param[i, 1:3] <- optimize_stationary_likelihood(
          likelihood_fun = OU_likelihood,
          data = sim_res_add$X_t[sim_res_add$t <= t_0],
          init_par = initital_param_stationary,
          delta = actual_dt, exp_sigma = TRUE)$par
      }, error = function(e){
        failOU <<- TRUE
      })
      
      tryCatch({
        dynamic_estim_param[i, 1:2] <- optimize_dynamic_likelihood(
          likelihood_fun = OU_dynamic_likelihood,
          data = sim_res_add$X_t[sim_res_add$t > t_0],
          init_par = dynamic_part_true_param * noise_dynamic_part,
          delta = actual_dt,
          alpha0 = stationary_estim_param[i, 1],
          mu0 = stationary_estim_param[i, 2],
          sigma = stationary_estim_param[i, 3],
          control = list(reltol = sqrt(.Machine$double.eps) / 100))$par
      }, error = function(e){
        failOU_dynamic <<- TRUE
      })
      
      tryCatch({
        stationary_estim_param[i, 4:6] <- optimize_stationary_likelihood(
          likelihood_fun = t_diffusion_strang_splitting,
          data = asinh(sim_res_add$X_t[sim_res_add$t < t_0]),
          init_par = initital_param_stationary,
          delta = actual_dt,
          exp_sigma = TRUE)$par
      }, error = function(e){
        failt <<- TRUE
      })
      
      tryCatch({
        dynamic_estim_param[i, 3:4] <- optimize_dynamic_likelihood(
          likelihood_fun = t_transform_dynamic_likelihood,
          data = asinh(sim_res_add$X_t[sim_res_add$t > t_0]),
          init_par = dynamic_part_true_param * noise_dynamic_part,
          delta = actual_dt,
          alpha0 = stationary_estim_param[i, 4],
          mu0 = stationary_estim_param[i, 5],
          sigma = stationary_estim_param[i, 6],
          control = list(reltol = sqrt(.Machine$double.eps) / 100))$par
      }, error = function(e){
        failt_dynamic <<- TRUE
      })
      
    tryCatch({
    
    quantiles_OU[i, ] <- unname(quantile(
      OU_likelihood_resid(
        par = stationary_estim_param[i, 1:3],
        data = sim_res_add$X_t[sim_res_add$t <= t_0],
        delta = actual_dt),
      probs = considered_quantiles))
    
    quantiles_OU_dynamic[i, ] <- unname(quantile(
      OU_dynamic_likelihood_resid(
        par = dynamic_estim_param[i, 1:2],
        data = sim_res_add$X_t[sim_res_add$t > t_0],
        delta = actual_dt,
        alpha0 = stationary_estim_param[i, 1],
        mu0 = stationary_estim_param[i, 2],
        sigma = stationary_estim_param[i, 3]
      ),
      probs = considered_quantiles))
    
    t_quantiles <-  t_diffusion_strang_splitting_resid(
        par = stationary_estim_param[i, 4:6],
        data = asinh(sim_res_add$X_t[sim_res_add$t <= t_0]),
        delta = actual_dt)
    
    quantiles_t[i, ] <- unname(quantile(t_quantiles[!is.na(t_quantiles)],
                                        probs = considered_quantiles))
    
    t_dynamic_quantiles <- t_transform_dynamic_likelihood_resid(
      par = dynamic_estim_param[i, 3:4],
      data = asinh(sim_res_add$X_t[sim_res_add$t > t_0]),
      delta = actual_dt,
      alpha0 = stationary_estim_param[i, 4],
      mu0 = stationary_estim_param[i, 5],
      sigma = stationary_estim_param[i, 6])
    
    quantiles_t_dynamic[i, ] <- unname(quantile(t_dynamic_quantiles[!is.na(t_dynamic_quantiles)],
                                                probs = considered_quantiles))

    }, error = function(e){
      fail_quantiles <- TRUE
    })
    if (any(is.na(quantiles_t_dynamic[i, ]) | is.nan(quantiles_t_dynamic[i, ]))) {

      fail_quantiles <- TRUE
    }
    if (!failOU && !failOU_dynamic && !failt_dynamic && !failt && !fail_quantiles) {
      success <- TRUE
    }
  }
}
  
  # Store the results for each tau value in the list
  results_list$quantiles_OU[[as.character(tau)]] <- quantiles_OU
  results_list$quantiles_OU_dynamic[[as.character(tau)]] <- quantiles_OU_dynamic
  results_list$quantiles_t[[as.character(tau)]] <- quantiles_t
  results_list$quantiles_t_dynamic[[as.character(tau)]] <- quantiles_t_dynamic
  results_list$stationary_estim_param[[as.character(tau)]] <- stationary_estim_param
  results_list$dynamic_estim_param[[as.character(tau)]] <- dynamic_estim_param
}


combined_data_quantiles <- list()

for (tau in tau_values) {

  df_OU <- as.data.frame(results_list$quantiles_OU[[as.character(tau)]]) |>
    mutate(tau = tau, Model = "OU")
  df_OU_dynamic <- as.data.frame(results_list$quantiles_OU_dynamic[[as.character(tau)]]) |>
    mutate(tau = tau, Model = "OU_dynamic")
  df_t <- as.data.frame(results_list$quantiles_t[[as.character(tau)]]) |>
    mutate(tau = tau, Model = "t")
  df_t_dynamic <- as.data.frame(results_list$quantiles_t_dynamic[[as.character(tau)]]) |>
    mutate(tau = tau, Model = "t_dynamic")
  

  combined_tau_df <- bind_rows(df_OU, df_OU_dynamic, df_t, df_t_dynamic)
  

  combined_data_quantiles[[as.character(tau)]] <- combined_tau_df
}

# Bind all the data frames in the list into one big data frame
final_combined_tau_quantiles <- as_tibble(bind_rows(combined_data_quantiles))

names(final_combined_tau_quantiles) <- c(considered_quantiles, "tau", "Model")

# if(!file.exists("data/tau_quantiles.csv")){
#   utils::write.table(final_combined_tau_quantiles,
#                      file="data/tau_quantiles.csv", sep = ",", row.names = FALSE)
# } else{
#   final_combined_tau_quantiles <- read_csv("data/tau_quantiles.csv")
# }


final_combined_tau_quantiles_long <- final_combined_tau_quantiles |> pivot_longer(-c(tau, Model), names_to = "Quantile", values_to = "Empiric_quantile") |> 
  mutate(Model = factor(Model), tau = factor(tau), Quantile = as.numeric(Quantile),
         theoretical_quantile = qnorm(Quantile))

quantiles_plot_tau <- final_combined_tau_quantiles_long |> 
  filter(str_detect(Model, "dynamic")) |> 
  ggplot(aes(x = theoretical_quantile, y = Empiric_quantile, col = tau)) +
  geom_smooth(method = "lm", se = F, linewidth = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 1) +
  facet_wrap(~Model, labeller = labeller(Model = c("OU_dynamic" = "Additive",
                                                   "t_dynamic" = "t-diffusion"))) +
  scale_color_manual(values = thesis_palette) +
  guides(col = guide_legend(override.aes = list(linewidth = 8))) + 
  labs(color = expression(tau)) + xlab("Theoretical quantiles") + 
  ylab("Empirical quantiles") + 
  theme(
    strip.text = element_text(size = 24, face = "bold"),
    axis.title = element_text(face = "bold", size = 24),
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 28),
    axis.text = element_text(size = 24)
  )

# ggsave(quantiles_plot_tau, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "quantiles_plot_tau.jpeg",
#        height = 7, width = 16, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)

stationary_parameters_misspecification <- map2_dfr(results_list$stationary_estim_param, tau_values, ~{
  as_tibble(.x) |> mutate(tau = .y)
})

names(stationary_parameters_misspecification) <- c("alpha_OU", "mu_OU", "sigma_OU",
                                                   "alpha_t", "mu_t", "sigma_t", "tau")

dynamic_parameters_misspecification <- map2_dfr(results_list$dynamic_estim_param, tau_values, ~{
  as_tibble(.x) |> mutate(tau = .y)
})

names(dynamic_parameters_misspecification) <- c("tau_OU", "A_OU", "tau_t", "A_t", "tau")


# if(!file.exists("data/dynamic_parameters_misspecification.csv")){
#   utils::write.table(dynamic_parameters_misspecification,
#                      file="data/dynamic_parameters_misspecification.csv", sep = ",", row.names = FALSE)
# } else{
#   dynamic_parameters_misspecification <- read_csv("data/dynamic_parameters_misspecification.csv")
# }

dynamic_parameters_misspecification_long <- dynamic_parameters_misspecification |> 
  mutate(tau_OU = (tau_OU - tau) / tau,
  A_OU = (A_OU - true_param[1]) / true_param[1],
  tau_t = (tau_t - tau) / tau,
  A_t = (A_t - true_param[1]) / true_param[1])  |> 
  pivot_longer(cols = -tau, 
               names_to = c("Metric", "Model"), 
               names_sep = "_", 
               values_to = "RE")



# Overview of extreme cases
dynamic_parameters_misspecification_long |>
  filter(abs(RE) >= 1.5, Metric != "A") |> group_by(tau, Model) |>
  summarise(count = n()) |> pivot_wider(values_from = count, names_from = Model) |>
  mutate(t = ifelse(is.na(t), 0, t)) |> xtable::xtable()

RE_dist_tau <- dynamic_parameters_misspecification_long |> 
  filter(abs(RE) < 1.5, Metric != "A") |> 
  ggplot(aes(x = factor(tau), y = RE, fill = Model)) +
  geom_violin() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) + 
  scale_fill_manual(values = thesis_palette, labels = c("Additive", "t-diffusion")) +
  ylab("Relative error") +
  xlab(expression(tau)) +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 20)),
    axis.title.y = element_text(size = 25, face = "bold"),
    axis.text = element_text(size = 24, face = "bold"),
    legend.title = element_text(size = 26, face = "bold"),
    legend.text = element_text(size = 24)
  )

# ggsave(RE_dist_tau, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "RE_dist_tau.jpeg",
#        height = 7, width = 16, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)

###------------------------------------------------------------------------###
###-----------------------------Result-------------------------------------###
###------------------------------------------------------------------------###
### Estimation on the AMOC

# Map of the AMOC
bbox <- c(xmin = -70, xmax = 5, ymin = 45, ymax = 70)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

NorthAtlanticOcean <- ggplot(data = world) +
  geom_sf(fill = "white", color = "black") + 
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE) +
  labs(title = "", x = "Longitude", y = "Latitude") + 
  theme(
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 16)
  )

ggsave(NorthAtlanticOcean, path = paste0(getwd(), "/tex_files/figures"),
       filename = "NorthAtlanticOcean.jpeg",
       height = 8, width = 14, dpi = 300, units = "in", device = "jpeg",
       limitsize = FALSE, scale = 1)

AMOC_data <- readr::read_table("data/AMOCdata.txt")

AMOC_data <- dplyr::rename_all(AMOC_data, ~ gsub('"', '',.))

AMOC_data <- mutate(AMOC_data, AMOC3 = AMOC0 - 3 * GM)



AMOC_data_longer <- AMOC_data |>
  pivot_longer(cols = -time, names_to = "Type", values_to = "Value") 

AMOC_data_plot <- AMOC_data_longer  |> 
  filter(!(Type %in% c("AMOC1", "AMOC3")))  |> 
  mutate(Type = case_when(
    Type == "AMOC0" ~ "SG SST anomaly",
    Type == "GM" ~ "GM SST anomaly",
    Type == "AMOC2" ~ "AMOC fingerprint"
  ))  |> 
  mutate(Type = factor(Type, levels = c("SG SST anomaly", "GM SST anomaly", "AMOC fingerprint")))  |> 
  ggplot(aes(x = time, y = Value, col = Type)) +
  geom_line(linewidth = 1) + 
  facet_wrap(~Type, ncol = 1, scales = "free_y") + 
  scale_x_continuous(breaks = scales::pretty_breaks(10)) +
  scale_color_manual(values = c(thesis_palette[1],thesis_palette[6], thesis_palette[3])) +
  theme(
    legend.title = element_blank(),
    strip.text = element_blank(),
    legend.key.size = unit(1, 'lines'),
    legend.position = "bottom",
    legend.text = element_text(face = "bold", size = 22),
    axis.text = element_text(face = "bold", size = 24),
    axis.title = element_text(face = "bold", size = 26),
  ) + 
  xlab("Year") + ylab("[K]") + 
  guides(col = guide_legend(override.aes = list(linewidth = 10)))

# ggsave(AMOC_data_plot, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "AMOC_data_plot.jpeg",
#        height = 12, width = 22, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)

AMOC_alt_plot <- AMOC_data_longer  |> 
  filter((Type %in% c("AMOC1", "AMOC2", "AMOC3")))  |> 
  mutate(Type = factor(Type, levels = c("AMOC1", "AMOC2", "AMOC3")))  |> 
  ggplot(aes(x = time, y = Value, col = Type)) +
  geom_line(linewidth = 0.6) + 
  facet_wrap(~Type) + 
  scale_x_continuous(breaks = scales::pretty_breaks(5)) +
  scale_color_manual(values = c(thesis_palette[1],thesis_palette[6], thesis_palette[3])) +
  theme(
    legend.title = element_blank(),
    strip.text = element_blank(),
    legend.key.size = unit(2, 'lines'),
    legend.text = element_text(face = "bold", size = 22),
    axis.text = element_text(face = "bold", size = 22),
    axis.title = element_text(face = "bold", size = 22),
  ) + 
  xlab("Year") + ylab("[K]") + 
  guides(col = guide_legend(override.aes = list(linewidth = 4)))

# ggsave(AMOC_alt_plot, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "AMOC_alt_plot.jpeg",
#        height = 7, width = 14, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)

# Temporal resolution 
actual_dt <- 1 / 12 

# Sweeping over starting years for ramping:

end_year <- 1950
start_year <- 1910
t_0s <- seq(start_year, end_year, length.out = (end_year + 1 - start_year)) 

# # MSE Method 
# MSE <- numeric(length(length(t_0s)))
# for (i in seq_along(t_0s)) {
#     segment_data <- filter(AMOC_data, time > t_0s[i])
#     lambda_fit <- lm(AMOC2 ~ time, data = segment_data)
#     MSE[i] <- mean(resid(lambda_fit)^2)
# }
# 
# MSE_tibble <- tibble(MSE = MSE) |> mutate(Start_Year = t_0s)
# 
# MSE_tibble |> 
#   ggplot(aes(x = Start_Year, y = MSE)) + geom_line()
# 
# t_0_v1 <- MSE_tibble$Start_Year[which.min(MSE_tibble$MSE)]


# Set starting point of tau to be the time between initialization of ramping and max year in data set
dynamic_part_starting_param <- matrix(NA, nrow = length(t_0s), ncol = 2)
dynamic_part_starting_param[, 1] <- max(AMOC_data$time) - t_0s 
dynamic_part_starting_param[, 2] <- 1

stationary_params <- matrix(nrow = length(t_0s), ncol = 3)
dynamic_params <- matrix(nrow = length(t_0s), ncol = 2)
stationary_part_likelihood <- numeric(length(t_0s))
dynamic_part_likelihood <- numeric(length(t_0s))

for (i in seq_along(t_0s)){
# Stationary part
print(i)
  
stationary_estim <- optimize_stationary_likelihood(
    likelihood_fun = t_diffusion_strang_splitting,
    data = asinh(AMOC_data$AMOC2[AMOC_data$time < t_0s[i]]),
    init_par = OU_init_params(data = AMOC_data$AMOC2[AMOC_data$time < t_0s[i]], delta = actual_dt),
    delta = actual_dt,
    exp_sigma = TRUE,
    control = list(reltol = sqrt(.Machine$double.eps) / 1000))
stationary_params[i, ] <- stationary_estim$par
stationary_part_likelihood[i] <- stationary_estim$objective


# Dynamic part

dynamic_estim <- optimize_dynamic_likelihood(
  likelihood_fun = t_transform_dynamic_likelihood,
  data = asinh(AMOC_data$AMOC2[AMOC_data$time >= t_0s[i]]),
  init_par = dynamic_part_starting_param[i , ],
  delta = actual_dt,
  alpha0 = stationary_estim$par[1],
  mu0 = stationary_estim$par[2],
  sigma = stationary_estim$par[3],
  method = "BFGS",
  control = list(reltol = sqrt(.Machine$double.eps) / 1000))
print(dynamic_estim$objective)
dynamic_params[i, ] <- dynamic_estim$par

dynamic_part_likelihood[i] <- dynamic_estim$objective
}

likelihood_tibble <- tibble(stationary_likelihood = stationary_part_likelihood,
       dynamic_likelihood = dynamic_part_likelihood,
       index = t_0s,
       n_stationary = index - min(AMOC_data$time),
       n_dynamic = max(AMOC_data$time) - index,
       alpha0 = stationary_params[, 1],
       mu0 = stationary_params[, 2],
       sigma = stationary_params[, 3],
       tau = dynamic_params[, 1],
       tipping_point = tau + t_0s,
       A = dynamic_params[, 2]) |> 
  filter(dynamic_part_likelihood != 50000) |>
  mutate(dynamic_likelihood = dynamic_part_likelihood / n_dynamic,
         stationary_likelihood = stationary_likelihood / n_stationary) |> 
  select(-c(n_stationary, n_dynamic))

stationary_min_row <- likelihood_tibble[which.min(likelihood_tibble$stationary_likelihood),]
dynamic_min_row <- likelihood_tibble[which.min(likelihood_tibble$dynamic_likelihood),]

ramping_year_likelihood_plot <- likelihood_tibble |> 
  select(index, stationary_likelihood, dynamic_likelihood) |> 
  pivot_longer(cols = -index, names_to = c("Part", "Parameter"), names_sep = "_", values_to = "Objective") |> 
  select(Part, Objective, index) |>
  ggplot(aes(x = index, y = Objective, col = Part)) + 
  geom_step(linewidth = 1.25) + 
  geom_hline(
    data = tibble(index = c(dynamic_min_row$index, stationary_min_row$index),
                  Objective = c(dynamic_min_row$dynamic_likelihood, stationary_min_row$stationary_likelihood),
                  Part = c("dynamic", "stationary")),
    aes(yintercept = Objective, col = Part), linetype = "dashed", linewidth = 0.8) + 
  facet_wrap(~Part) +
  geom_vline(
    data = tibble(index = c(dynamic_min_row$index, stationary_min_row$index),
                    Objective = c(dynamic_min_row$dynamic_likelihood, stationary_min_row$stationary_likelihood),
                    Part = c("stationary", "dynamic")),
             aes(xintercept = index, col = Part), linetype = "dashed", linewidth = 0.8, show.legend = FALSE) + 
  geom_point(
    data = tibble(index = c(dynamic_min_row$index, stationary_min_row$index),
                  Objective = c(dynamic_min_row$dynamic_likelihood, stationary_min_row$stationary_likelihood),
                  Part = c("dynamic", "stationary")), size = 3, show.legend = FALSE) +
  geom_text(
    data = tibble(index = c(stationary_min_row$index),
                  Objective = c(stationary_min_row$stationary_likelihood),
                  Part = c("stationary")),
    aes(label = index), nudge_y = -0.015 , show.legend = FALSE, size=5) +
  geom_text(
    data = tibble(index = c(dynamic_min_row$index),
                  Objective = c(dynamic_min_row$dynamic_likelihood),
                  Part = c("dynamic")),
    aes(label = index), nudge_y = -0.015, show.legend = FALSE, size=5) +
  scale_color_manual(labels = c("Dynamic", "Stationary"), values = thesis_palette[5:6]) + 
  guides(col = guide_legend(override.aes = list(linewidth = 4))) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +  
  xlab(expression(t[0])) +
  theme(
    strip.text = element_blank(),
    legend.text = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 18),
    axis.text = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 18),
    theme(panel.spacing = unit(2, "lines"))
  )

# ggsave(ramping_year_likelihood_plot, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "ramping_year_likelihood_plot.jpeg",
#        height = 5, width = 15, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)

# Add min likelihood points to estimate plot 
min_dynamic_value_wide <- likelihood_tibble |>
  slice_min(dynamic_likelihood) |> 
  select(-c(stationary_likelihood, dynamic_likelihood))

min_dynamic_value <- min_dynamic_value_wide |> 
  pivot_longer(cols = -index, names_to = "Parameter", values_to = "Estimate") |> 
  mutate(Part = "Dynamic")

min_stationary_value_wide <- likelihood_tibble |>
  slice_min(stationary_likelihood) |> 
  select(-c(stationary_likelihood, dynamic_likelihood))

min_stationary_value <- min_stationary_value_wide |> 
  pivot_longer(cols = -index, names_to = "Parameter", values_to = "Estimate") |> 
  mutate(Part = "Stationary")

# Table
min_dynamic_value |> inner_join(min_stationary_value, by = "Parameter") |> 
  select(Parameter, Estimate.x, Estimate.y) |> 
  mutate(`Relative difference` = 100 * (Estimate.x - Estimate.y) / Estimate.y) |> 
  rename(`Estimates 1924` = Estimate.x, `Estimate_1912` = Estimate.y) |> 
  xtable::xtable()
                                
# Combine the extracted values
highlight_points <- bind_rows(min_dynamic_value, min_stationary_value)



estimators_base_plot <- likelihood_tibble |> 
  select(-c(stationary_likelihood, dynamic_likelihood)) |> 
  pivot_longer(cols = -index, names_to = "Parameter", values_to = "Estimate") |>
  ggplot(aes(x = index, y = Estimate)) + 
  geom_step(linewidth = 1.25, aes(col = Parameter)) +
  facet_wrap(Parameter~., scales = "free_y") + 
  xlab(expression(t[0])) +
  theme(
    strip.text = element_blank(),
    legend.text = element_text(size = 18),
    legend.title = element_text(face = "bold", size = 18),
    axis.text = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 18),
    theme(panel.spacing = unit(2, "lines"))
  ) + 
  scale_color_manual(labels = c("A", expression(alpha),  expression(mu), expression(sigma), expression(tau),
                               "Tipping year"), 
                               values = thesis_palette[-c(5,6,7, 8)]) + 
  guides(col = guide_legend(override.aes = list(linewidth = 4))) +
  scale_y_continuous(breaks = scales::pretty_breaks())

estimators_full_plot <- estimators_base_plot +  ggnewscale::new_scale_colour() + 
geom_point(data = highlight_points, aes(x = index, y = Estimate, col = Part), size = 3) +
  geom_hline(data = highlight_points, aes(yintercept = Estimate, col = Part),
             linewidth = 0.8, linetype = "dashed", show.legend = FALSE) +
  scale_color_manual(labels = c("1924", "1912"), values = thesis_palette[5:6]) +
  guides(col = guide_legend(override.aes = list(linewidth = 4), title = expression(t[0])))

# ggsave(estimators_full_plot, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "estimators_full_plot.jpeg",
#        height = 8, width = 15, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)


## QQ-analysis
QQ_data_parts <- bind_rows(
  # Part based on t_0 = 1924
  tibble(sample = t_diffusion_strang_splitting_resid(
  par = c(min_dynamic_value_wide$alpha0, min_dynamic_value_wide$mu0, min_dynamic_value_wide$sigma), 
  data = asinh(AMOC_data$AMOC2[AMOC_data$time < min_dynamic_value_wide$index]),
  delta = actual_dt
)) |> mutate(Optimized = "Dynamic", Part = "Stationary part"),

  tibble(sample = t_transform_dynamic_likelihood_resid(
   par =  c(min_dynamic_value_wide$tau, min_dynamic_value_wide$A),
   data = asinh(AMOC_data$AMOC2[AMOC_data$time >= min_dynamic_value_wide$index]),
   delta = actual_dt,
   alpha0 = min_dynamic_value_wide$alpha0,
   mu0 = min_dynamic_value_wide$mu0,
   sigma = min_dynamic_value_wide$sigma)
   )|> mutate(Optimized = "Dynamic", Part = "Dynamic part"),

  # Part based on t_0 = 1912
  tibble(sample = t_diffusion_strang_splitting_resid(
  par = c(min_stationary_value_wide$alpha0, min_stationary_value_wide$mu0, min_stationary_value_wide$sigma), 
  data = asinh(AMOC_data$AMOC2[AMOC_data$time < min_stationary_value_wide$index]),
  delta = actual_dt
)) |> mutate(Optimized = "Stationary", Part = "Stationary part"),

tibble(sample = t_transform_dynamic_likelihood_resid(
  par =  c(min_stationary_value_wide$tau, min_stationary_value_wide$A),
  data = asinh(AMOC_data$AMOC2[AMOC_data$time >= min_stationary_value_wide$index]),
  delta = actual_dt,
  alpha0 = min_stationary_value_wide$alpha0,
  mu0 = min_stationary_value_wide$mu0,
  sigma = min_stationary_value_wide$sigma)
)|> mutate(Optimized = "Stationary", Part = "Dynamic part")
) |> 
  mutate(Part = factor(Part, levels = c("Stationary part", "Dynamic part")),
         Optimized = factor(Optimized, levels = c("Stationary", "Dynamic")))

QQ_plot_parts <- QQ_data_parts |> 
  ggplot(aes(sample = sample)) + geom_qq(aes(col = Optimized), alpha = .8, size = 1.3) + geom_qq_line() + 
  facet_wrap(~Part) + 
  scale_color_manual(labels = c("1912", "1924"), values = thesis_palette[5:6]) +
  theme(
    strip.text = element_text(face = "bold", size = 22),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22),
    axis.title = element_text(face = "bold", size = 20)
  ) + xlab("Theoretical Quantiles") + ylab("Empirical Quantiles") +
  guides(col = guide_legend(override.aes = list(alpha = 1, size = 5), title = expression(t[0])))

quantiles_t0 <- c(0.01, 0.025, 0.167, 0.33, 0.5, 0.67, 0.833, 0.975, 0.99)

QQ_table_data <- QQ_data_parts |> group_by(Part, Optimized) |> 
  reframe(empirical_quantiles = 
              quantile(sample, probs = quantiles_t0),
          probs = factor(quantiles_t0))

# Table for quantiles
bind_rows(
  tibble(Part = "Theoretical", Optimized = "Theoretical", probs = factor(quantiles_t0),
         quantiles = qnorm(p = quantiles_t0)) |> 
    pivot_wider(id_cols = c(Part, Optimized), values_from = quantiles, names_from = probs), 
QQ_table_data |> 
  pivot_wider(id_cols = c(Part, Optimized), values_from = empirical_quantiles, names_from = probs)
) |> xtable::xtable()


# ggsave(QQ_plot_parts, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "QQ_plot_parts.jpeg",
#        height = 8, width = 15, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)

# Estimation in the AMOC

t_0 <- 1924
# Simulate from the model to construct parametric bootstrap confidence intervals
numSim <- 2000

# Define parameters - note that we shift the first year, 1870, to be year 0.
T_0 <- t_0 - 1870
tau_estim <- min_dynamic_value_wide$tau
A_estim <- min_dynamic_value_wide$A
alpha_0_estim <- min_dynamic_value_wide$alpha0
mu0_estim <- min_dynamic_value_wide$mu0
sigma_estim <- min_dynamic_value_wide$sigma

lambda_0_estim <- -alpha_0_estim^2 / (4 * A_estim)
m_estim <- mu0_estim - alpha_0_estim / (2 * A_estim)

tc_estim <- tau_estim + t_0
# Parameters needed for sim method along with time to tipping from largest year in data set (we stop simulating at that point)
sim_param <- c(A_estim, m_estim, lambda_0_estim, sigma_estim)
time_to_tipping <- tc_estim - max(AMOC_data$time) 
 


original_estim <- tibble(true_value = c(A_estim, 
                                        alpha_0_estim, lambda_0_estim,  m_estim,
                                        mu0_estim, sigma_estim, tau_estim, tc_estim))

estim_matrix <- matrix(data = NA, nrow = numSim, ncol = 8)
#tc_vector <- numeric(numSim)
set.seed(07052024)
for (i in 1:numSim) {
  if (i %% 5 == 1) {
    cat("Currently grinding", i, "....\n")
  }
  
  success <- FALSE
  
  while (!success) {
    tryCatch({
      random_seed <- sample(100000, size = 1)
      sim_t <- simulate_t_distribution_tipping_model(step_length = actual_dt, par = sim_param, tau = tau_estim,
                                             t_0 = T_0, beyond_tipping = -time_to_tipping, seed = random_seed)
      
      # Stationary part
      sim_t_stationary_estim <- optimize_stationary_likelihood(
        likelihood_fun = t_diffusion_strang_splitting,
        data = asinh(sim_t$X_t[sim_t$t < T_0]),
        init_par = c(alpha_0_estim, mu0_estim, sigma_estim),
        delta = actual_dt,
        exp_sigma = TRUE,
        method = "BFGS"
      )$par
      
      # Dynamic part
      sim_t_dynamic_estim <- optimize_dynamic_likelihood(
        likelihood_fun = t_transform_dynamic_likelihood,
        data = asinh(sim_t$X_t[sim_t$t >= T_0]),
        init_par = c(tau_estim, A_estim),
        delta = actual_dt,
        alpha0 = sim_t_stationary_estim[1],
        mu0 = sim_t_stationary_estim[2],
        sigma = sim_t_stationary_estim[3],
        method = "BFGS",
        control = list(reltol = sqrt(.Machine$double.eps) / 1000)
      )$par
      
      estim_matrix[i, 1:8] <- c(
        sim_t_dynamic_estim[2], 
        sim_t_stationary_estim[1],
        -sim_t_stationary_estim[1]^2 / (4 * sim_t_dynamic_estim[2]),
        sim_t_stationary_estim[2] - sim_t_stationary_estim[1] / (2 * sim_t_dynamic_estim[2]),
        sim_t_stationary_estim[2],
        sim_t_stationary_estim[3],
        sim_t_dynamic_estim[1],
        T_0 + sim_t_dynamic_estim[1] + 1870
      )
      
      success <- TRUE
    }, error = function(e) {
      cat("Error occurred at iteration", i, ":", conditionMessage(e), "\n")
    })
  }
  print(estim_matrix[i, 8])
}


estim_tibble <- as_tibble(estim_matrix)

names(estim_tibble) <- c("A", "alpha_0", "lambda_0", "m", "mu0", "sigma", "tau", "tc")

if(!file.exists("data/original_estim.csv")){
  utils::write.table(original_estim, file="data/original_estim.csv", sep = ",", row.names = FALSE)
} else{
  original_estim <- read_csv("data/original_estim.csv")
}

if(!file.exists("data/estim_tibble.csv")){
  utils::write.table(estim_tibble, file="data/estim_tibble.csv", sep = ",", row.names = FALSE)
} else{
  estim_tibble <- read_csv("data/estim_tibble.csv")
}

original_estim <- original_estim |> mutate(Parameter = names(estim_tibble))

# Remove a couple of outliers
estim_tibble_trimmed <- estim_tibble |> filter(A < 10 & lambda_0 > -7.5)

estim_tibble_long <- estim_tibble_trimmed |> 
  pivot_longer(cols = everything(), names_to = "Parameter", values_to = "value")

estim_tibble_plot <- estim_tibble_long |> ggplot(aes(x = value, y = after_stat(density), fill = Parameter)) +
  geom_histogram(bins = 30, alpha = 2, col = "black", linewidth = .2) + 
  geom_vline(data = original_estim, mapping = aes(xintercept = true_value), linewidth = 1, linetype = "dashed") +
  facet_wrap(~Parameter, scales = "free", ncol = 4) +
  labs(x = "Estimate", y = "") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) + 
  theme(strip.text.x = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        axis.title.x = element_text(face = "bold", size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)
        ) + 
  scale_fill_manual(values = thesis_palette,
                    labels = c("A", expression(alpha*phantom(.)[0]), expression(lambda*phantom(.)[0]), "m",
                               expression(mu*phantom(.)[0]), expression(sigma), expression(tau*phantom(.)[c], "Tipping year")))

# ggplot2::ggsave(filename = "tex_files/figures/estim_tibble_plot.jpeg",
#                 plot = estim_tibble_plot,
#                 width = 15,
#                 height = 7.5,
#                 dpi = 300)

# Extreme observations
estim_tibble |> filter(!(A < 10 & lambda_0 > -7.5)) |> arrange(lambda_0) |> 
  xtable::xtable()

# Quantiles for tipping time
tibble::tibble(
Quantile = c("2.5%", "16.5%", "50%", "83.%", "97.5%"),
value = quantile(estim_tibble_trimmed$tc, prob = c(0.025, 0.165, 0.5, 0.835, 0.975))
) |> pivot_wider(names_from = Quantile, values_from = value) |> 
  xtable::xtable()
  
  
qqplot_before_t0 <- tibble::tibble(obsSample =
                 t_diffusion_strang_splitting_resid(
                 par = c(original_estim$true_value[2], original_estim$true_value[5], original_estim$true_value[6]),
                 data = asinh(AMOC_data$AMOC2[AMOC_data$time < t_0]),
                 delta = actual_dt)) |>
                 ggplot2::ggplot(ggplot2::aes(sample = obsSample)) +
                 ggplot2::geom_qq() + ggplot2::geom_qq_line() + 
                 xlab("Theoretical Quantiles") + ylab("Sample Quantiles")  


qqplot_after_t0 <-  tibble::tibble(obsSample =
                 t_transform_dynamic_likelihood_resid(
                 par = c(original_estim$true_value[7], original_estim$true_value[1]),
                 data = asinh(AMOC_data$AMOC2[AMOC_data$time >= t_0]),
                 delta = actual_dt,
                 alpha0 = original_estim$true_value[2],
                 mu0 =  original_estim$true_value[5], 
                 sigma = original_estim$true_value[6])) |>
                 ggplot2::ggplot(ggplot2::aes(sample = obsSample)) +
                 ggplot2::geom_qq() + ggplot2::geom_qq_line() + 
                 xlab("Theoretical Quantiles") + ylab("Sample Quantiles")

qqplot_combined <- gridExtra::grid.arrange(qqplot_before_t0, qqplot_after_t0, ncol = 2)


# ggplot2::ggsave(filename = "tex_files/figures/qqplot_combined.jpeg",
#                 plot = qqplot_combined,
#                 width = 15,
#                 height = 7.5,
#                 dpi = 300)

# Using additive model but without penalization 

# Set starting point of tau to be the time between initialization of ramping and max year in data set
dynamic_part_starting_param_additive <- matrix(NA, nrow = length(t_0s), ncol = 2)
dynamic_part_starting_param_additive[, 1] <- max(AMOC_data$time) - t_0s 
dynamic_part_starting_param_additive[, 2] <- 1

stationary_params_additive <- matrix(nrow = length(t_0s), ncol = 3)
dynamic_params_additive <- matrix(nrow = length(t_0s), ncol = 2)
stationary_part_likelihood_additive <- numeric(length(t_0s))
dynamic_part_likelihood_additive <- numeric(length(t_0s))

for (i in seq_along(t_0s)){
  # Stationary part
  print(i)
  
  stationary_estim <- optimize_stationary_likelihood(
    likelihood_fun = OU_likelihood,
    data = AMOC_data$AMOC2[AMOC_data$time < t_0s[i]],
    init_par = OU_init_params(data = AMOC_data$AMOC2[AMOC_data$time < t_0s[i]], delta = actual_dt),
    delta = actual_dt,
    exp_sigma = TRUE,
    control = list(reltol = sqrt(.Machine$double.eps) / 1000))
  stationary_params_additive[i, ] <- stationary_estim$par
  stationary_part_likelihood_additive[i] <- stationary_estim$objective
  

  # Dynamic part
  
  dynamic_estim <- optimize_dynamic_likelihood(
    likelihood_fun = OU_dynamic_likelihood,
    data = AMOC_data$AMOC2[AMOC_data$time >= t_0s[i]],
    init_par = dynamic_part_starting_param_additive[i , ],
    delta = actual_dt,
    alpha0 = stationary_estim$par[1],
    mu0 = stationary_estim$par[2],
    sigma = stationary_estim$par[3],
    method = "BFGS",
    control = list(reltol = sqrt(.Machine$double.eps) / 1000))
  print(dynamic_estim$objective)
  dynamic_params_additive[i, ] <- dynamic_estim$par
  print(dynamic_estim$par)
  dynamic_part_likelihood_additive[i] <- dynamic_estim$objective
}


likelihood_tibble_additive <- tibble(stationary_likelihood = stationary_part_likelihood_additive,
                            dynamic_likelihood = dynamic_part_likelihood_additive,
                            index = t_0s,
                            n_stationary_add = (index - 1870),
                            n_dynamic_add = (max(AMOC_data$time) - index),
                            alpha0 = stationary_params_additive[, 1],
                            mu0 = stationary_params_additive[, 2],
                            sigma = stationary_params_additive[, 3],
                            tau = dynamic_params_additive[, 1],
                            tipping_point = tau + t_0s,
                            A = dynamic_params_additive[, 2]#,
                            #nu = dynamic_params_additive[, 3]
                            ) |> 
  mutate(dynamic_likelihood = dynamic_likelihood / n_dynamic_add,
         stationary_likelihood = stationary_likelihood / n_stationary_add) |> 
  select(-c(n_stationary_add, n_dynamic_add))

stationary_min_row_additive <- likelihood_tibble_additive[which.min(likelihood_tibble_additive$stationary_likelihood),]
dynamic_min_row_additive <- likelihood_tibble_additive[which.min(likelihood_tibble_additive$dynamic_likelihood),]


ramping_year_likelihood_plot_additive <- likelihood_tibble_additive |> 
  select(index, stationary_likelihood, dynamic_likelihood) |> 
  pivot_longer(cols = -index, names_to = c("Part", "Parameter"), names_sep = "_", values_to = "Objective") |> 
  select(Part, Objective, index) |>
  ggplot(aes(x = index, y = Objective, col = Part)) + 
  geom_step(linewidth = 1.25) + 
  geom_hline(
    data = tibble(index = c(dynamic_min_row_additive$index, stationary_min_row_additive$index),
                  Objective = c(dynamic_min_row_additive$dynamic_likelihood, stationary_min_row_additive$stationary_likelihood),
                  Part = c("dynamic", "stationary")),
    aes(yintercept = Objective, col = Part), linetype = "dashed", linewidth = 0.8) + 
  facet_wrap(~Part) +
  geom_vline(
    data = tibble(index = c(dynamic_min_row_additive$index, stationary_min_row_additive$index),
                  Objective = c(dynamic_min_row_additive$dynamic_likelihood, stationary_min_row_additive$stationary_likelihood),
                  Part = c("stationary", "dynamic")),
    aes(xintercept = index, col = Part), linetype = "dashed", linewidth = 0.8, show.legend = FALSE) + 
  geom_point(
    data = tibble(index = c(dynamic_min_row_additive$index, stationary_min_row_additive$index),
                  Objective = c(dynamic_min_row_additive$dynamic_likelihood, stationary_min_row_additive$stationary_likelihood),
                  Part = c("dynamic", "stationary")), size = 3, show.legend = FALSE) +
  geom_text(
    data = tibble(index = c(stationary_min_row_additive$index),
                  Objective = c(stationary_min_row_additive$stationary_likelihood),
                  Part = c("stationary")),
    aes(label = index), nudge_y = -0.015 , show.legend = FALSE, size=5) +
  geom_text(
    data = tibble(index = c(dynamic_min_row_additive$index),
                  Objective = c(dynamic_min_row_additive$dynamic_likelihood),
                  Part = c("dynamic")),
    aes(label = index), nudge_y = -0.015, show.legend = FALSE, size=5) +
  scale_color_manual(labels = c("Dynamic", "Stationary"), values = thesis_palette[5:6]) + 
  guides(col = guide_legend(override.aes = list(linewidth = 4))) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +  
  xlab(expression(t[0])) +
  theme(
    strip.text = element_blank(),
    legend.text = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 18),
    axis.text = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 18),
    theme(panel.spacing = unit(2, "lines"))
  )

# Add min likelihood points to estimate plot 
min_dynamic_value_additive_wide <- likelihood_tibble_additive |>
  slice_min(dynamic_likelihood) |> 
  select(-c(stationary_likelihood, dynamic_likelihood))

min_dynamic_additive_value <- min_dynamic_value_additive_wide |> 
  pivot_longer(cols = -index, names_to = "Parameter", values_to = "Estimate") |> 
  mutate(Part = "Dynamic")

min_stationary_additive_value_wide <- likelihood_tibble_additive |>
  slice_min(stationary_likelihood) |> 
  select(-c(stationary_likelihood, dynamic_likelihood))

min_stationary_additive_value <- min_stationary_additive_value_wide |> 
  pivot_longer(cols = -index, names_to = "Parameter", values_to = "Estimate") |> 
  mutate(Part = "Stationary")

# Table
min_dynamic_additive_value |> inner_join(min_stationary_additive_value, by = "Parameter") |> 
  select(Parameter, Estimate.x, Estimate.y) |> 
  mutate(`Relative difference` = 100 * (Estimate.x - Estimate.y) / Estimate.y) |> 
  rename(`Estimates 1924` = Estimate.x, `Estimate_1912` = Estimate.y) |> 
  xtable::xtable()

# Combine the extracted values
highlight_points_additive <- bind_rows(min_dynamic_additive_value, min_stationary_additive_value)



estimators_base_plot_additive <- likelihood_tibble_additive |> 
  select(-c(stationary_likelihood, dynamic_likelihood)) |> 
  pivot_longer(cols = -index, names_to = "Parameter", values_to = "Estimate") |>
  ggplot(aes(x = index, y = Estimate)) + 
  geom_step(linewidth = 1.25, aes(col = Parameter)) +
  facet_wrap(Parameter~., scales = "free_y") + 
  xlab(expression(t[0])) +
  theme(
    strip.text = element_blank(),
    legend.text = element_text(size = 18),
    legend.title = element_text(face = "bold", size = 18),
    axis.text = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 18),
    theme(panel.spacing = unit(2, "lines"))
  ) + 
  scale_color_manual(labels = c("A", expression(alpha),  expression(mu), expression(sigma), expression(tau),
                                "Tipping year"), 
                     values = thesis_palette[-c(5,6,7, 8)]) + 
  guides(col = guide_legend(override.aes = list(linewidth = 4))) +
  scale_y_continuous(breaks = scales::pretty_breaks())

estimators_full_plot_additive <- estimators_base_plot_additive +  ggnewscale::new_scale_colour() + 
  geom_point(data = highlight_points_additive, aes(x = index, y = Estimate, col = Part), size = 3) +
  geom_hline(data = highlight_points_additive, aes(yintercept = Estimate, col = Part),
             linewidth = 0.8, linetype = "dashed", show.legend = FALSE) +
  scale_color_manual(labels = c("1915", "1943"), values = thesis_palette[5:6]) +
  guides(col = guide_legend(override.aes = list(linewidth = 4), title = expression(t[0])))

## QQ plot
## QQ-analysis
QQ_data_parts_additive <- bind_rows(
  # Part based on t_0 = 1943
  tibble(sample = OU_likelihood_resid(
    par = c(min_dynamic_value_additive_wide$alpha0, min_dynamic_value_additive_wide$mu0, min_dynamic_value_additive_wide$sigma), 
    data = AMOC_data$AMOC2[AMOC_data$time < min_dynamic_value_additive_wide$index],
    delta = actual_dt
  )) |> mutate(Optimized = "Dynamic", Part = "Stationary part"),
  
  tibble(sample = OU_dynamic_likelihood_resid(
    par =  c(min_dynamic_value_additive_wide$tau, min_dynamic_value_additive_wide$A),
    data = AMOC_data$AMOC2[AMOC_data$time >= min_dynamic_value_additive_wide$index],
    delta = actual_dt,
    alpha0 = min_dynamic_value_additive_wide$alpha0,
    mu0 = min_dynamic_value_additive_wide$mu0,
    sigma = min_dynamic_value_additive_wide$sigma)
  )|> mutate(Optimized = "Dynamic", Part = "Dynamic part"),
  
  # Part based on t_0 = 1912
  tibble(sample = OU_likelihood_resid(
    par = c(min_stationary_additive_value_wide$alpha0, min_stationary_additive_value_wide$mu0, min_stationary_additive_value_wide$sigma), 
    data = AMOC_data$AMOC2[AMOC_data$time < min_stationary_additive_value_wide$index],
    delta = actual_dt
  )) |> mutate(Optimized = "Stationary", Part = "Stationary part"),
  
  tibble(sample = OU_dynamic_likelihood_resid(
    par =  c(min_stationary_additive_value_wide$tau, min_stationary_additive_value_wide$A),
    data = AMOC_data$AMOC2[AMOC_data$time >= min_stationary_additive_value_wide$index],
    delta = actual_dt,
    alpha0 = min_stationary_additive_value_wide$alpha0,
    mu0 = min_stationary_additive_value_wide$mu0,
    sigma = min_stationary_additive_value_wide$sigma)
  )|> mutate(Optimized = "Stationary", Part = "Dynamic part")
) |> 
  mutate(Part = factor(Part, levels = c("Stationary part", "Dynamic part")),
         Optimized = factor(Optimized, levels = c("Stationary", "Dynamic")))

QQ_plot_parts_additve <- QQ_data_parts_additive |> 
  ggplot(aes(sample = sample)) + geom_qq(aes(col = Optimized), alpha = .8, size = 1.3) + geom_qq_line() + 
  facet_wrap(~Part) + 
  scale_color_manual(labels = c("1915", "1943"), values = thesis_palette[5:6]) +
  theme(
    strip.text = element_text(face = "bold", size = 22),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22),
    axis.title = element_text(face = "bold", size = 20)
  ) + xlab("Theoretical Quantiles") + ylab("Empirical Quantiles") +
  guides(col = guide_legend(override.aes = list(alpha = 1, size = 5), title = expression(t[0])))

QQ_plot_parts

t_0_additive <- 1915
actual_dt <- 1 / 12
# Simulate from the model to construct parametric bootstrap confidence intervals
numSim <- 2000

# Define parameters - note that we shift the first year, 1870, to be year 0.
T_0_additive <- t_0_additive - 1870
tau_estim_additive <- min_dynamic_value_additive_wide$tau
A_estim_additive <- min_dynamic_value_additive_wide$A
alpha_0_estim_additive <- min_dynamic_value_additive_wide$alpha0
mu0_estim_additive <- min_dynamic_value_additive_wide$mu0
sigma_estim_additive <- min_dynamic_value_additive_wide$sigma

lambda_0_estim_additive <- -alpha_0_estim_additive^2 / (4 * A_estim_additive)
m_estim_additive <- mu0_estim_additive - alpha_0_estim_additive / (2 * A_estim_additive)

tc_estim_additive <- tau_estim_additive + t_0_additive
# Parameters needed for sim method along with time to tipping from largest year in data set (we stop simulating at that point)
sim_param_additive <- c(A_estim_additive, m_estim_additive, lambda_0_estim_additive, sigma_estim_additive)
time_to_tipping_additive <- tc_estim_additive - max(AMOC_data$time) 



original_estim_additive <- tibble(true_value = c(A_estim_additive, 
                                        alpha_0_estim_additive, lambda_0_estim_additive,  m_estim_additive,
                                        mu0_estim_additive, sigma_estim_additive, tau_estim_additive, tc_estim_additive))

#dynamic_part_init <- c(max(AMOC_data$time) - t_0, 1)
estim_matrix_additive <- matrix(data = NA, nrow = numSim, ncol = 8)
#tc_vector <- numeric(numSim)
set.seed(07052024)
for (i in 1:numSim) {
  if (i %% 5 == 1) {
    cat("Currently grinding", i, "....\n")
  }
  
  success <- FALSE
  
  while (!success) {
    tryCatch({
      random_seed <- sample(100000, size = 1)
      OU_sim <- simulate_additive_noise_tipping_model(step_length = actual_dt,
                                                     par = sim_param_additive, tau = tau_estim_additive,
                                                     t_0 = T_0_additive,
                                                     beyond_tipping = -time_to_tipping_additive, seed = random_seed)
      
      # Stationary part
      sim_OU_stationary_estim <- optimize_stationary_likelihood(
        likelihood_fun = OU_likelihood,
        data = OU_sim$X_t[OU_sim$t < T_0_additive],
        init_par = c(alpha_0_estim_additive, mu0_estim_additive, sigma_estim_additive),
        delta = actual_dt,
        exp_sigma = TRUE,
        method = "BFGS"
      )$par
      
      # Dynamic part
      sim_OU_dynamic_estim <- optimize_dynamic_likelihood(
        likelihood_fun = OU_dynamic_likelihood,
        data = OU_sim$X_t[OU_sim$t >= T_0_additive],
        init_par = c(tau_estim_additive, A_estim_additive),
        delta = actual_dt,
        alpha0 = sim_OU_stationary_estim[1],
        mu0 = sim_OU_stationary_estim[2],
        sigma = sim_OU_stationary_estim[3],
        method = "BFGS",
        control = list(reltol = sqrt(.Machine$double.eps) / 1000)
      )$par
      
      estim_matrix_additive[i, 1:8] <- c(
        sim_OU_dynamic_estim[2], 
        sim_OU_stationary_estim[1],
        -sim_OU_stationary_estim[1]^2 / (4 * sim_OU_dynamic_estim[2]),
        sim_OU_stationary_estim[2] - sim_OU_stationary_estim[1] / (2 * sim_OU_dynamic_estim[2]),
        sim_OU_stationary_estim[2],
        sim_OU_stationary_estim[3],
        sim_OU_dynamic_estim[1],
        T_0_additive + sim_OU_dynamic_estim[1] + 1870
      )
      
      success <- TRUE
    }, error = function(e) {
      cat("Error occurred at iteration", i, ":", conditionMessage(e), "\n")
    })
  }
  print(estim_matrix_additive[i, 8])
}


estim_tibble_additive <- as_tibble(estim_matrix_additive)

names(estim_tibble_additive) <- c("A", "alpha_0", "lambda_0", "m", "mu0", "sigma", "tau", "tc")

if(!file.exists("data/original_estim_additive.csv")){
  utils::write.table(original_estim_additive, file="data/original_estim_additive.csv", sep = ",", row.names = FALSE)
} else{
  original_estim_additive <- read_csv("data/original_estim_additive.csv")
}

if(!file.exists("data/estim_tibble_additive.csv")){
  utils::write.table(estim_tibble_additive, file="data/estim_tibble_additive.csv", sep = ",", row.names = FALSE)
} else{
  estim_tibble_additive <- read_csv("data/estim_tibble_additive.csv")
}

original_estim_additive <- original_estim_additive |> mutate(Parameter = names(estim_tibble_additive))

estim_tibble_additive_long <- estim_tibble_additive |> 
  pivot_longer(cols = everything(), names_to = "Parameter", values_to = "value")

estim_tibble_additive_plot <- estim_tibble_additive_long |> ggplot(aes(x = value, y = after_stat(density), fill = Parameter)) +
  geom_histogram(bins = 30, alpha = 2, col = "black", linewidth = .2) + 
  geom_vline(data = original_estim_additive, mapping = aes(xintercept = true_value), linewidth = 1, linetype = "dashed") +
  facet_wrap(~Parameter, scales = "free", ncol = 4) +
  labs(x = "Estimate", y = "") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) + 
  theme(strip.text.x = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        axis.title.x = element_text(face = "bold", size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)
  ) + 
  scale_fill_manual(values = thesis_palette,
                    labels = c("A", expression(alpha*phantom(.)[0]), expression(lambda*phantom(.)[0]), "m",
                               expression(mu*phantom(.)[0]), expression(sigma), expression(tau*phantom(.)[c], "Tipping year")))


# ggsave(estim_tibble_additive_plot, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "estim_tibble_additive_plot.jpeg",
#        height = 7, width = 15, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)

tibble::tibble(
  Quantile = c("2.5%", "16.5%", "50%", "83.5%", "97.5%"),
  value = quantile(estim_tibble_additive$tc, prob = c(0.025, 0.165, 0.5, 0.835, 0.975))
) |> pivot_wider(names_from = Quantile, values_from = value) |> 
  xtable::xtable()

# OU vs t-diffusion

OU_vs_t_diffusion_QQ_plot <- bind_rows(QQ_data_parts_additive |> mutate(model = "Additive"),
          QQ_data_parts |> mutate(model = "t-diffusion")) |> 
  filter(Optimized == "Dynamic") |> 
  ggplot(aes(sample = sample, col = model)) +
  geom_qq() + 
  geom_qq_line(aes(group = 1), linewidth = 0.85, color = "black", linetype = "dashed") +
  facet_wrap(~Part) + scale_color_manual(values = thesis_palette) + 
  theme(
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 20),
    strip.text = element_text(face = "bold", size = 20),
    axis.title = element_text(face = "bold", size = 20)
  ) + 
  guides(col = guide_legend(override.aes = list(size = 4))) + 
  xlab("Theoretical Quantiles") + ylab("Empirical Quantiles") 

# ggsave(OU_vs_t_diffusion_QQ_plot, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "OU_vs_t_diffusion_QQ_plot.jpeg",
#        height = 6, width = 14, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)

OU_vs_t_diffusion_table <- bind_rows(QQ_data_parts_additive |> mutate(model = "Additive"),
          QQ_data_parts |> mutate(model = "t-diffusion")) |> 
  filter(Optimized == "Dynamic") |> 
  group_by(Part, model) |> 
  reframe(quantile = quantile(sample, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
          probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)) |> 
  pivot_wider(id_cols = c(Part, model), values_from = quantile, names_from = probs) |> 
  bind_rows(tibble(Part = "Theoretical",
                   model = "-",
                   quantile = qnorm(c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)),
                   probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)) |>
                   pivot_wider(id_cols = c(Part, model), values_from = quantile, names_from = probs))

# Robustness analysis

# Set starting point of tau to be the time between initialization of ramping and max year in data set

dynamic_part_starting_param_AMOC1 <- matrix(NA, nrow = length(t_0s), ncol = 2)
dynamic_part_starting_param_AMOC1[, 1] <- max(AMOC_data$time) - t_0s 
dynamic_part_starting_param_AMOC1[, 2] <- 1

stationary_params_AMOC1 <- matrix(nrow = length(t_0s), ncol = 3)
dynamic_params_AMOC1 <- matrix(nrow = length(t_0s), ncol = 2)
stationary_part_likelihood_AMOC1 <- numeric(length(t_0s))
dynamic_part_likelihood_AMOC1 <- numeric(length(t_0s))

for (i in seq_along(t_0s)){
  # Stationary part
  print(i)
  
  stationary_estim_AMOC1 <- optimize_stationary_likelihood(
    likelihood_fun = t_diffusion_strang_splitting,
    data = asinh(AMOC_data$AMOC1[AMOC_data$time < t_0s[i]]),
    init_par = OU_init_params(data = AMOC_data$AMOC1[AMOC_data$time < t_0s[i]], delta = actual_dt),
    delta = actual_dt,
    exp_sigma = TRUE,
    control = list(reltol = sqrt(.Machine$double.eps) / 1000))
  stationary_params_AMOC1[i, ] <- stationary_estim_AMOC1$par
  stationary_part_likelihood_AMOC1[i] <- stationary_estim_AMOC1$objective
  
  
  # Dynamic part
  
  dynamic_estim_AMOC1 <- optimize_dynamic_likelihood(
    likelihood_fun = t_transform_dynamic_likelihood,
    data = asinh(AMOC_data$AMOC1[AMOC_data$time >= t_0s[i]]),
    init_par = dynamic_part_starting_param[i , ],
    delta = actual_dt,
    alpha0 = stationary_estim$par[1],
    mu0 = stationary_estim$par[2],
    sigma = stationary_estim$par[3],
    method = "BFGS",
    control = list(reltol = sqrt(.Machine$double.eps) / 1000))
  print(dynamic_estim_AMOC1$objective)
  dynamic_params_AMOC1[i, ] <- dynamic_estim_AMOC1$par
  
  dynamic_part_likelihood_AMOC1[i] <- dynamic_estim_AMOC1$objective
}

likelihood_tibble_AMOC1 <- tibble(stationary_likelihood = stationary_part_likelihood,
                            dynamic_likelihood = dynamic_part_likelihood,
                            index = t_0s,
                            n_stationary = index - min(AMOC_data$time),
                            n_dynamic = max(AMOC_data$time) - index,
                            alpha0 = stationary_params_AMOC1[, 1],
                            mu0 = stationary_params_AMOC1[, 2],
                            sigma = stationary_params_AMOC1[, 3],
                            tau = dynamic_params_AMOC1[, 1],
                            tipping_point = tau + t_0s,
                            A = dynamic_params_AMOC1[, 2]) |> 
  filter(dynamic_part_likelihood != 50000) |>
  mutate(dynamic_likelihood = dynamic_part_likelihood / n_dynamic,
         stationary_likelihood = stationary_likelihood / n_stationary) |> 
  select(-c(n_stationary, n_dynamic))

stationary_min_row_AMOC1 <- likelihood_tibble_AMOC1[which.min(likelihood_tibble_AMOC1$stationary_likelihood),]
dynamic_min_row_AMOC1 <- likelihood_tibble_AMOC1[which.min(likelihood_tibble_AMOC1$dynamic_likelihood),]



# Define parameters - note that we shift the first year, 1870, to be year 0.
T_0 <- t_0 - 1870
tau_estim_AMOC1 <- dynamic_min_row_AMOC1$tau
A_estim_AMOC1 <- dynamic_min_row_AMOC1$A
alpha_0_estim_AMOC1 <- dynamic_min_row_AMOC1$alpha0
mu0_estim_AMOC1 <- dynamic_min_row_AMOC1$mu0
sigma_estim_AMOC1 <- dynamic_min_row_AMOC1$sigma

lambda_0_estim_AMOC1 <- -alpha_0_estim_AMOC1^2 / (4 * A_estim_AMOC1)
m_estim_AMOC1 <- mu0_estim_AMOC1 - alpha_0_estim_AMOC1 / (2 * A_estim_AMOC1)

tc_estim_AMOC1 <- tau_estim_AMOC1 + t_0
# Parameters needed for sim method along with time to tipping from largest year in data set (we stop simulating at that point)
sim_param_AMOC1 <- c(A_estim_AMOC1, m_estim_AMOC1, lambda_0_estim_AMOC1, sigma_estim_AMOC1)
time_to_tipping_AMOC1 <- tc_estim_AMOC1 - max(AMOC_data$time) 



original_estim_AMOC1 <- tibble(true_value = c(A_estim_AMOC1, 
                                        alpha_0_estim_AMOC1, lambda_0_estim_AMOC1,  m_estim_AMOC1,
                                        mu0_estim_AMOC1, sigma_estim_AMOC1, tau_estim_AMOC1, tc_estim_AMOC1))

numSim <- 2000
estim_matrix_AMOC1 <- matrix(data = NA, nrow = numSim, ncol = 8)

set.seed(07052024)
for (i in 1:numSim) {
  if (i %% 5 == 1) {
    cat("Currently grinding", i, "....\n")
  }
  
  success <- FALSE
  
  while (!success) {
    tryCatch({
      random_seed <- sample(100000, size = 1)
      sim_t <- simulate_t_distribution_tipping_model(step_length = actual_dt,
                                                     par = sim_param_AMOC1, tau = tau_estim_AMOC1,
                                                     t_0 = T_0,
                                                     beyond_tipping = -time_to_tipping_AMOC1, seed = random_seed)
      
      # Stationary part
      sim_t_stationary_estim <- optimize_stationary_likelihood(
        likelihood_fun = t_diffusion_strang_splitting,
        data = asinh(sim_t$X_t[sim_t$t < T_0]),
        init_par = c(alpha_0_estim_AMOC1, mu0_estim_AMOC1, sigma_estim_AMOC1),
        delta = actual_dt,
        exp_sigma = TRUE,
        method = "BFGS"
      )$par
      
      # Dynamic part
      sim_t_dynamic_estim <- optimize_dynamic_likelihood(
        likelihood_fun = t_transform_dynamic_likelihood,
        data = asinh(sim_t$X_t[sim_t$t >= T_0]),
        init_par = c(tau_estim_AMOC1, A_estim_AMOC1),
        delta = actual_dt,
        alpha0 = sim_t_stationary_estim[1],
        mu0 = sim_t_stationary_estim[2],
        sigma = sim_t_stationary_estim[3],
        method = "BFGS",
        control = list(reltol = sqrt(.Machine$double.eps) / 1000)
      )$par
      
      estim_matrix_AMOC1[i, 1:8] <- c(
        sim_t_dynamic_estim[2], 
        sim_t_stationary_estim[1],
        -sim_t_stationary_estim[1]^2 / (4 * sim_t_dynamic_estim[2]),
        sim_t_stationary_estim[2] - sim_t_stationary_estim[1] / (2 * sim_t_dynamic_estim[2]),
        sim_t_stationary_estim[2],
        sim_t_stationary_estim[3],
        sim_t_dynamic_estim[1],
        T_0 + sim_t_dynamic_estim[1] + 1870
      )
      
      success <- TRUE
    }, error = function(e) {
      cat("Error occurred at iteration", i, ":", conditionMessage(e), "\n")
    })
  }
  print(estim_matrix_AMOC1[i, 8])
}


estim_tibble_AMOC1 <- as_tibble(estim_matrix_AMOC1)

names(estim_tibble_AMOC1) <- c("A", "alpha_0", "lambda_0", "m", "mu0", "sigma", "tau", "tc")

if(!file.exists("data/original_estim_AMOC1.csv")){
  utils::write.table(original_estim_AMOC1, file="data/original_estim_AMOC1.csv", sep = ",", row.names = FALSE)
} else{
  original_estim_AMOC1 <- read_csv("data/original_estim_AMOC1.csv")
}

if(!file.exists("data/estim_tibble_AMOC1.csv")){
  utils::write.table(estim_tibble_AMOC1, file="data/estim_tibble_AMOC1.csv", sep = ",", row.names = FALSE)
} else{
  estim_tibble_AMOC1 <- read_csv("data/estim_tibble_AMOC1.csv")
}

original_estim_AMOC1 <- original_estim_AMOC1 |> mutate(Parameter = names(estim_tibble_AMOC1))

estim_tibble_AMOC1_long <- estim_tibble_AMOC1 |> filter(A < 10, A > 0, lambda_0 > -7.5) |> 
  pivot_longer(cols = everything(), names_to = "Parameter", values_to = "value")

estim_tibble_AMOC1_plot <- estim_tibble_AMOC1_long |>
  ggplot(aes(x = value, y = after_stat(density), fill = Parameter)) +
  geom_histogram(bins = 30, alpha = 2, col = "black", linewidth = .2) + 
  geom_vline(data = original_estim_AMOC1, mapping = aes(xintercept = true_value),
             linewidth = 1, linetype = "dashed") +
  facet_wrap(~Parameter, scales = "free", ncol = 4) +
  labs(x = "Estimate", y = "") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) + 
  theme(strip.text.x = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        axis.title.x = element_text(face = "bold", size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)
  ) + 
  scale_fill_manual(values = thesis_palette,
                    labels = c("A", expression(alpha*phantom(.)[0]), expression(lambda*phantom(.)[0]), "m",
                               expression(mu*phantom(.)[0]), expression(sigma), expression(tau*phantom(.)[c], "Tipping year")))


# ggsave(estim_tibble_AMOC1_plot, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "estim_tibble_AMOC1_plot.jpeg",
#        height = 7, width = 15, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)

# Set starting point of tau to be the time between initialization of ramping and max year in data set

dynamic_part_starting_param_AMOC3 <- matrix(NA, nrow = length(t_0s), ncol = 2)
dynamic_part_starting_param_AMOC3[, 1] <- max(AMOC_data$time) - t_0s 
dynamic_part_starting_param_AMOC3[, 2] <- 1

stationary_params_AMOC3 <- matrix(nrow = length(t_0s), ncol = 3)
dynamic_params_AMOC3 <- matrix(nrow = length(t_0s), ncol = 2)
stationary_part_likelihood_AMOC3 <- numeric(length(t_0s))
dynamic_part_likelihood_AMOC3 <- numeric(length(t_0s))

for (i in seq_along(t_0s)){
  # Stationary part
  print(i)
  
  stationary_estim_AMOC3 <- optimize_stationary_likelihood(
    likelihood_fun = t_diffusion_strang_splitting,
    data = asinh(AMOC_data$AMOC3[AMOC_data$time < t_0s[i]]),
    init_par = OU_init_params(data = AMOC_data$AMOC3[AMOC_data$time < t_0s[i]], delta = actual_dt),
    delta = actual_dt,
    exp_sigma = TRUE,
    control = list(reltol = sqrt(.Machine$double.eps) / 1000))
  stationary_params_AMOC3[i, ] <- stationary_estim_AMOC3$par
  stationary_part_likelihood_AMOC3[i] <- stationary_estim_AMOC3$objective
  
  
  # Dynamic part
  
  dynamic_estim_AMOC3 <- optimize_dynamic_likelihood(
    likelihood_fun = t_transform_dynamic_likelihood,
    data = asinh(AMOC_data$AMOC3[AMOC_data$time >= t_0s[i]]),
    init_par = dynamic_part_starting_param[i , ],
    delta = actual_dt,
    alpha0 = stationary_estim_AMOC3$par[1],
    mu0 = stationary_estim_AMOC3$par[2],
    sigma = stationary_estim_AMOC3$par[3],
    method = "BFGS",
    control = list(reltol = sqrt(.Machine$double.eps) / 1000))
  print(dynamic_estim_AMOC3$objective)
  dynamic_params_AMOC3[i, ] <- dynamic_estim_AMOC3$par
  
  dynamic_part_likelihood_AMOC3[i] <- dynamic_estim_AMOC3$objective
}

likelihood_tibble_AMOC3 <- tibble(stationary_likelihood = stationary_part_likelihood_AMOC3,
                                  dynamic_likelihood = dynamic_part_likelihood_AMOC3,
                                  index = t_0s,
                                  n_stationary = index - min(AMOC_data$time),
                                  n_dynamic = max(AMOC_data$time) - index,
                                  alpha0 = stationary_params_AMOC3[, 1],
                                  mu0 = stationary_params_AMOC3[, 2],
                                  sigma = stationary_params_AMOC3[, 3],
                                  tau = dynamic_params_AMOC3[, 1],
                                  tipping_point = tau + t_0s,
                                  A = dynamic_params_AMOC3[, 2]) |> 
  filter(dynamic_likelihood != 50000) |>
  mutate(dynamic_likelihood = dynamic_likelihood / n_dynamic,
         stationary_likelihood = stationary_likelihood / n_stationary) |> 
  select(-c(n_stationary, n_dynamic))

stationary_min_row_AMOC3 <- likelihood_tibble_AMOC3[which.min(likelihood_tibble_AMOC3$stationary_likelihood),]
dynamic_min_row_AMOC3 <- likelihood_tibble_AMOC3[which.min(likelihood_tibble_AMOC3$dynamic_likelihood),]

# Define parameters - note that we shift the first year, 1870, to be year 0.
T_0 <- dynamic_min_row_AMOC3$index - 1870
tau_estim_AMOC3 <- dynamic_min_row_AMOC3$tau
A_estim_AMOC3 <- dynamic_min_row_AMOC3$A
alpha_0_estim_AMOC3 <- dynamic_min_row_AMOC3$alpha0
mu0_estim_AMOC3 <- dynamic_min_row_AMOC3$mu0
sigma_estim_AMOC3 <- dynamic_min_row_AMOC3$sigma

lambda_0_estim_AMOC3 <- -alpha_0_estim_AMOC3^2 / (4 * A_estim_AMOC3)
m_estim_AMOC3 <- mu0_estim_AMOC3 - alpha_0_estim_AMOC3 / (2 * A_estim_AMOC3)

tc_estim_AMOC3 <- tau_estim_AMOC3 + dynamic_min_row_AMOC3$index
# Parameters needed for sim method along with time to tipping from largest year in data set (we stop simulating at that point)
sim_param_AMOC3 <- c(A_estim_AMOC3, m_estim_AMOC3, lambda_0_estim_AMOC3, sigma_estim_AMOC3)
time_to_tipping_AMOC3 <- tc_estim_AMOC3 - max(AMOC_data$time) 



original_estim_AMOC3 <- tibble(true_value = c(A_estim_AMOC3, 
                                              alpha_0_estim_AMOC3, lambda_0_estim_AMOC3,  m_estim_AMOC3,
                                              mu0_estim_AMOC3, sigma_estim_AMOC3, tau_estim_AMOC3, tc_estim_AMOC3))

numSim <- 2000
estim_matrix_AMOC3 <- matrix(data = NA, nrow = numSim, ncol = 8)

set.seed(07052024)
for (i in 1:numSim) {
  if (i %% 5 == 1) {
    cat("Currently grinding", i, "....\n")
  }
  
  success <- FALSE
  
  while (!success) {
    tryCatch({
      random_seed <- sample(100000, size = 1)
      sim_t <- simulate_t_distribution_tipping_model(step_length = actual_dt,
                                                     par = sim_param_AMOC3, tau = tau_estim_AMOC3,
                                                     t_0 = T_0,
                                                     beyond_tipping = -time_to_tipping_AMOC3, seed = random_seed)
      
      # Stationary part
      sim_t_stationary_estim <- optimize_stationary_likelihood(
        likelihood_fun = t_diffusion_strang_splitting,
        data = asinh(sim_t$X_t[sim_t$t < T_0]),
        init_par = c(alpha_0_estim_AMOC3, mu0_estim_AMOC3, sigma_estim_AMOC3),
        delta = actual_dt,
        exp_sigma = TRUE,
        method = "BFGS"
      )$par
      
      # Dynamic part
      sim_t_dynamic_estim <- optimize_dynamic_likelihood(
        likelihood_fun = t_transform_dynamic_likelihood,
        data = asinh(sim_t$X_t[sim_t$t >= T_0]),
        init_par = c(tau_estim_AMOC3, A_estim_AMOC3),
        delta = actual_dt,
        alpha0 = sim_t_stationary_estim[1],
        mu0 = sim_t_stationary_estim[2],
        sigma = sim_t_stationary_estim[3],
        method = "BFGS",
        control = list(reltol = sqrt(.Machine$double.eps) / 1000)
      )$par
      
      estim_matrix_AMOC3[i, 1:8] <- c(
        sim_t_dynamic_estim[2], 
        sim_t_stationary_estim[1],
        -sim_t_stationary_estim[1]^2 / (4 * sim_t_dynamic_estim[2]),
        sim_t_stationary_estim[2] - sim_t_stationary_estim[1] / (2 * sim_t_dynamic_estim[2]),
        sim_t_stationary_estim[2],
        sim_t_stationary_estim[3],
        sim_t_dynamic_estim[1],
        T_0 + sim_t_dynamic_estim[1] + 1870
      )
      
      success <- TRUE
    }, error = function(e) {
      cat("Error occurred at iteration", i, ":", conditionMessage(e), "\n")
    })
  }
  print(estim_matrix_AMOC3[i, 8])
}


estim_tibble_AMOC3 <- as_tibble(estim_matrix_AMOC3)

names(estim_tibble_AMOC3) <- c("A", "alpha_0", "lambda_0", "m", "mu0", "sigma", "tau", "tc")

if(!file.exists("data/original_estim_AMOC3.csv")){
  utils::write.table(original_estim_AMOC3, file="data/original_estim_AMOC3.csv", sep = ",", row.names = FALSE)
} else{
  original_estim_AMOC3 <- read_csv("data/original_estim_AMOC3.csv")
}

if(!file.exists("data/estim_tibble_AMOC3.csv")){
  utils::write.table(estim_tibble_AMOC3, file="data/estim_tibble_AMOC3.csv", sep = ",", row.names = FALSE)
} else{
  estim_tibble_AMOC3 <- read_csv("data/estim_tibble_AMOC3.csv")
}

original_estim_AMOC3 <- original_estim_AMOC3 |> mutate(Parameter = names(estim_tibble_AMOC3))

estim_tibble_AMOC3_long <- estim_tibble_AMOC3 |> filter(A < 10, A > 0) |> 
  pivot_longer(cols = everything(), names_to = "Parameter", values_to = "value")

estim_tibble_AMOC3_plot <- estim_tibble_AMOC3_long |> ggplot(aes(x = value, y = after_stat(density), fill = Parameter)) +
  geom_histogram(bins = 30, alpha = 2, col = "black", linewidth = .2) + 
  geom_vline(data = original_estim_AMOC3, mapping = aes(xintercept = true_value), linewidth = 1, linetype = "dashed") +
  facet_wrap(~Parameter, scales = "free", ncol = 4) +
  labs(x = "Estimate", y = "") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) + 
  theme(strip.text.x = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        axis.title.x = element_text(face = "bold", size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)
  ) + 
  scale_fill_manual(values = thesis_palette,
                    labels = c("A", expression(alpha*phantom(.)[0]), expression(lambda*phantom(.)[0]), "m",
                               expression(mu*phantom(.)[0]), expression(sigma), expression(tau*phantom(.)[c], "Tipping year")))

# ggsave(estim_tibble_AMOC3_plot, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "estim_tibble_AMOC3_plot.jpeg",
#        height = 7, width = 15, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)


# AMOC1 for additive model
dynamic_part_starting_param_additive_AMOC1 <- matrix(NA, nrow = length(t_0s), ncol = 2)
dynamic_part_starting_param_additive_AMOC1[, 1] <- max(AMOC_data$time) - t_0s 
dynamic_part_starting_param_additive_AMOC1[, 2] <- 1

stationary_params_additive_AMOC1 <- matrix(nrow = length(t_0s), ncol = 3)
dynamic_params_additive_AMOC1 <- matrix(nrow = length(t_0s), ncol = 2)
stationary_part_likelihood_additive_AMOC1 <- numeric(length(t_0s))

dynamic_part_likelihood_additive_AMOC1 <- numeric(length(t_0s))

for (i in seq_along(t_0s)){
  # Stationary part
  print(i)
  
  stationary_estim <- optimize_stationary_likelihood(
    likelihood_fun = OU_likelihood,
    data = AMOC_data$AMOC1[AMOC_data$time < t_0s[i]],
    init_par = OU_init_params(data = AMOC_data$AMOC1[AMOC_data$time < t_0s[i]], delta = actual_dt),
    delta = actual_dt,
    exp_sigma = TRUE,
    control = list(reltol = sqrt(.Machine$double.eps) / 1000))
  stationary_params_additive_AMOC1[i, ] <- stationary_estim$par
  stationary_part_likelihood_additive_AMOC1[i] <- stationary_estim$objective
  
  
  # Dynamic part
  
  dynamic_estim <- optimize_dynamic_likelihood(
    likelihood_fun = OU_dynamic_likelihood,
    data = AMOC_data$AMOC1[AMOC_data$time >= t_0s[i]],
    init_par = dynamic_part_starting_param_additive[i , ],
    delta = actual_dt,
    alpha0 = stationary_estim$par[1],
    mu0 = stationary_estim$par[2],
    sigma = stationary_estim$par[3],
    method = "BFGS",
    control = list(reltol = sqrt(.Machine$double.eps) / 1000))
  print(dynamic_estim$objective)
  dynamic_params_additive_AMOC1[i, ] <- dynamic_estim$par
  print(dynamic_estim$par)
  dynamic_part_likelihood_additive_AMOC1[i] <- dynamic_estim$objective
}


likelihood_tibble_additive_AMOC1 <- tibble(stationary_likelihood = stationary_part_likelihood_additive_AMOC1,
                                     dynamic_likelihood = dynamic_part_likelihood_additive_AMOC1,
                                     index = t_0s,
                                     n_stationary_add = (index - 1870),
                                     n_dynamic_add = (max(AMOC_data$time) - index),
                                     alpha0 = stationary_params_additive_AMOC1[, 1],
                                     mu0 = stationary_params_additive_AMOC1[, 2],
                                     sigma = stationary_params_additive_AMOC1[, 3],
                                     tau = dynamic_params_additive_AMOC1[, 1],
                                     tipping_point = tau + t_0s,
                                     A = dynamic_params_additive_AMOC1[, 2]
) |> 
  mutate(dynamic_likelihood = dynamic_likelihood / n_dynamic_add,
         stationary_likelihood = stationary_likelihood / n_stationary_add) |> 
  select(-c(n_stationary_add, n_dynamic_add))

stationary_min_row_additive_AMOC1 <- likelihood_tibble_additive_AMOC1[which.min(likelihood_tibble_additive_AMOC1$stationary_likelihood),]
dynamic_min_row_additive_AMOC1 <- likelihood_tibble_additive_AMOC1[which.min(likelihood_tibble_additive_AMOC1$dynamic_likelihood),]


t_0_additive_AMOC1 <- dynamic_min_row_additive_AMOC1$index
actual_dt <- 1 / 12
# Simulate from the model to construct parametric bootstrap confidence intervals
numSim <- 2000

# Define parameters - note that we shift the first year, 1870, to be year 0.
T_0_additive_AMOC1 <- t_0_additive_AMOC1 - 1870
tau_estim_additive_AMOC1 <- dynamic_min_row_additive_AMOC1$tau
A_estim_additive_AMOC1 <- dynamic_min_row_additive_AMOC1$A
alpha_0_estim_additive_AMOC1 <- dynamic_min_row_additive_AMOC1$alpha0
mu0_estim_additive_AMOC1 <- dynamic_min_row_additive_AMOC1$mu0
sigma_estim_additive_AMOC1 <- dynamic_min_row_additive_AMOC1$sigma

lambda_0_estim_additive_AMOC1 <- -alpha_0_estim_additive_AMOC1^2 / (4 * A_estim_additive_AMOC1)
m_estim_additive_AMOC1 <- mu0_estim_additive_AMOC1 - alpha_0_estim_additive_AMOC1 / (2 * A_estim_additive_AMOC1)

tc_estim_additive_AMOC1 <- tau_estim_additive_AMOC1 + t_0_additive_AMOC1
# Parameters needed for sim method along with time to tipping from largest year 
# in data set (we stop simulating at that point)
sim_param_additive_AMOC1 <- c(A_estim_additive_AMOC1, m_estim_additive_AMOC1, 
                              lambda_0_estim_additive_AMOC1, sigma_estim_additive_AMOC1)
time_to_tipping_additive_AMOC1 <- tc_estim_additive_AMOC1 - max(AMOC_data$time) 



original_estim_additive_AMOC1 <- tibble(true_value = c(A_estim_additive_AMOC1, 
                           alpha_0_estim_additive_AMOC1, lambda_0_estim_additive_AMOC1,  m_estim_additive_AMOC1,
                           mu0_estim_additive_AMOC1, sigma_estim_additive_AMOC1, tau_estim_additive_AMOC1, tc_estim_additive_AMOC1))


estim_matrix_additive_AMOC1 <- matrix(data = NA, nrow = numSim, ncol = 8)
set.seed(07052024)
for (i in 1:numSim) {
  if (i %% 5 == 1) {
    cat("Currently grinding", i, "....\n")
  }
  
  success <- FALSE
  
  while (!success) {
    tryCatch({
      random_seed <- sample(100000, size = 1)
      OU_sim <- simulate_additive_noise_tipping_model(step_length = actual_dt,
                                                      par = sim_param_additive_AMOC1, tau = tau_estim_additive_AMOC1,
                                                      t_0 = T_0_additive_AMOC1,
                                                      beyond_tipping = -time_to_tipping_additive_AMOC1, seed = random_seed)
      
      # Stationary part
      sim_OU_stationary_estim <- optimize_stationary_likelihood(
        likelihood_fun = OU_likelihood,
        data = OU_sim$X_t[OU_sim$t < T_0_additive_AMOC1],
        init_par = c(alpha_0_estim_additive_AMOC1, mu0_estim_additive_AMOC1, sigma_estim_additive_AMOC1),
        delta = actual_dt,
        exp_sigma = TRUE,
        method = "BFGS"
      )$par
      
      # Dynamic part
      sim_OU_dynamic_estim <- optimize_dynamic_likelihood(
        likelihood_fun = OU_dynamic_likelihood,
        data = OU_sim$X_t[OU_sim$t >= T_0_additive_AMOC1],
        init_par = c(tau_estim_additive_AMOC1, A_estim_additive_AMOC1),
        delta = actual_dt,
        alpha0 = sim_OU_stationary_estim[1],
        mu0 = sim_OU_stationary_estim[2],
        sigma = sim_OU_stationary_estim[3],
        method = "BFGS",
        control = list(reltol = sqrt(.Machine$double.eps) / 1000)
      )$par
      
      estim_matrix_additive_AMOC1[i, 1:8] <- c(
        sim_OU_dynamic_estim[2], 
        sim_OU_stationary_estim[1],
        -sim_OU_stationary_estim[1]^2 / (4 * sim_OU_dynamic_estim[2]),
        sim_OU_stationary_estim[2] - sim_OU_stationary_estim[1] / (2 * sim_OU_dynamic_estim[2]),
        sim_OU_stationary_estim[2],
        sim_OU_stationary_estim[3],
        sim_OU_dynamic_estim[1],
        T_0_additive_AMOC1 + sim_OU_dynamic_estim[1] + 1870
      )
      
      success <- TRUE
    }, error = function(e) {
      cat("Error occurred at iteration", i, ":", conditionMessage(e), "\n")
    })
  }
  print(estim_matrix_additive_AMOC1[i, 8])
}


estim_tibble_additive_AMOC1 <- as_tibble(estim_matrix_additive_AMOC1)

names(estim_tibble_additive_AMOC1) <- c("A", "alpha_0", "lambda_0", "m", "mu0", "sigma", "tau", "tc")

if(!file.exists("data/original_estim_additive_AMOC1.csv")){
  utils::write.table(original_estim_additive_AMOC1, file="data/original_estim_additive_AMOC1.csv",
                     sep = ",", row.names = FALSE)
} else{
  original_estim_additive_AMOC1 <- read_csv("data/original_estim_additive_AMOC1.csv")
}

if(!file.exists("data/estim_tibble_additive_AMOC1.csv")){
  utils::write.table(estim_tibble_additive_AMOC1, file="data/estim_tibble_additive_AMOC1.csv",
                     sep = ",", row.names = FALSE)
} else{
  estim_tibble_additive_AMOC1 <- read_csv("data/estim_tibble_additive_AMOC1.csv")
}

original_estim_additive_AMOC1 <- original_estim_additive_AMOC1 |> mutate(Parameter = names(estim_tibble_additive_AMOC1))

estim_tibble_additive_AMOC1_long <- estim_tibble_additive_AMOC1 |> 
  pivot_longer(cols = everything(), names_to = "Parameter", values_to = "value")

estim_tibble_additive_AMOC1_plot <- estim_tibble_additive_AMOC1_long |> ggplot(aes(x = value, y = after_stat(density), fill = Parameter)) +
  geom_histogram(bins = 30, alpha = 2, col = "black", linewidth = .2) + 
  geom_vline(data = original_estim_additive_AMOC1, mapping = aes(xintercept = true_value), linewidth = 1, linetype = "dashed") +
  facet_wrap(~Parameter, scales = "free", ncol = 4) +
  labs(x = "Estimate", y = "") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) + 
  theme(strip.text.x = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        axis.title.x = element_text(face = "bold", size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)
  ) + 
  scale_fill_manual(values = thesis_palette,
                    labels = c("A", expression(alpha*phantom(.)[0]), expression(lambda*phantom(.)[0]), "m",
                               expression(mu*phantom(.)[0]), expression(sigma), expression(tau*phantom(.)[c], "Tipping year")))

# ggsave(estim_tibble_additive_AMOC1_plot, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "estim_tibble_additive_AMOC1_plot.jpeg",
#        height = 7, width = 15, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)



# AMOC1 for additive model
dynamic_part_starting_param_additive_AMOC3 <- matrix(NA, nrow = length(t_0s), ncol = 2)
dynamic_part_starting_param_additive_AMOC3[, 1] <- max(AMOC_data$time) - t_0s 
dynamic_part_starting_param_additive_AMOC3[, 2] <- 1

stationary_params_additive_AMOC3 <- matrix(nrow = length(t_0s), ncol = 3)
dynamic_params_additive_AMOC3 <- matrix(nrow = length(t_0s), ncol = 2)
stationary_part_likelihood_additive_AMOC3 <- numeric(length(t_0s))

dynamic_part_likelihood_additive_AMOC3 <- numeric(length(t_0s))

for (i in seq_along(t_0s)){
  # Stationary part
  print(i)
  
  stationary_estim <- optimize_stationary_likelihood(
    likelihood_fun = OU_likelihood,
    data = AMOC_data$AMOC3[AMOC_data$time < t_0s[i]],
    init_par = OU_init_params(data = AMOC_data$AMOC3[AMOC_data$time < t_0s[i]], delta = actual_dt),
    delta = actual_dt,
    exp_sigma = TRUE,
    control = list(reltol = sqrt(.Machine$double.eps) / 1000))
  stationary_params_additive_AMOC3[i, ] <- stationary_estim$par
  stationary_part_likelihood_additive_AMOC3[i] <- stationary_estim$objective
  
  
  # Dynamic part
  
  dynamic_estim <- optimize_dynamic_likelihood(
    likelihood_fun = OU_dynamic_likelihood,
    data = AMOC_data$AMOC3[AMOC_data$time >= t_0s[i]],
    init_par = dynamic_part_starting_param_additive[i , ],
    delta = actual_dt,
    alpha0 = stationary_estim$par[1],
    mu0 = stationary_estim$par[2],
    sigma = stationary_estim$par[3],
    method = "BFGS",
    control = list(reltol = sqrt(.Machine$double.eps) / 1000))
  print(dynamic_estim$objective)
  dynamic_params_additive_AMOC3[i, ] <- dynamic_estim$par
  print(dynamic_estim$par)
  dynamic_part_likelihood_additive_AMOC3[i] <- dynamic_estim$objective
}


likelihood_tibble_additive_AMOC3 <- tibble(stationary_likelihood = stationary_part_likelihood_additive_AMOC3,
                                           dynamic_likelihood = dynamic_part_likelihood_additive_AMOC3,
                                           index = t_0s,
                                           n_stationary_add = (index - 1870),
                                           n_dynamic_add = (max(AMOC_data$time) - index),
                                           alpha0 = stationary_params_additive_AMOC3[, 1],
                                           mu0 = stationary_params_additive_AMOC3[, 2],
                                           sigma = stationary_params_additive_AMOC3[, 3],
                                           tau = dynamic_params_additive_AMOC3[, 1],
                                           tipping_point = tau + t_0s,
                                           A = dynamic_params_additive_AMOC3[, 2]
) |> 
  mutate(dynamic_likelihood = dynamic_likelihood / n_dynamic_add,
         stationary_likelihood = stationary_likelihood / n_stationary_add) |> 
  select(-c(n_stationary_add, n_dynamic_add))

stationary_min_row_additive_AMOC3 <- likelihood_tibble_additive_AMOC3[which.min(likelihood_tibble_additive_AMOC3$stationary_likelihood),]
dynamic_min_row_additive_AMOC3 <- likelihood_tibble_additive_AMOC3[which.min(likelihood_tibble_additive_AMOC3$dynamic_likelihood),]


t_0_additive_AMOC3 <- dynamic_min_row_additive_AMOC3$index
actual_dt <- 1 / 12
# Simulate from the model to construct parametric bootstrap confidence intervals
numSim <- 2000

# Define parameters - note that we shift the first year, 1870, to be year 0.
T_0_additive_AMOC3 <- t_0_additive_AMOC3 - 1870
tau_estim_additive_AMOC3 <- dynamic_min_row_additive_AMOC3$tau
A_estim_additive_AMOC3 <- dynamic_min_row_additive_AMOC3$A
alpha_0_estim_additive_AMOC3 <- dynamic_min_row_additive_AMOC3$alpha0
mu0_estim_additive_AMOC3 <- dynamic_min_row_additive_AMOC3$mu0
sigma_estim_additive_AMOC3 <- dynamic_min_row_additive_AMOC3$sigma

lambda_0_estim_additive_AMOC3 <- -alpha_0_estim_additive_AMOC3^2 / (4 * A_estim_additive_AMOC3)
m_estim_additive_AMOC3 <- mu0_estim_additive_AMOC3 - alpha_0_estim_additive_AMOC3 / (2 * A_estim_additive_AMOC3)

tc_estim_additive_AMOC3 <- tau_estim_additive_AMOC3 + t_0_additive_AMOC3
# Parameters needed for sim method along with time to tipping from largest year in data set (we stop simulating at that point)
sim_param_additive_AMOC3 <- c(A_estim_additive_AMOC3, m_estim_additive_AMOC3, 
                              lambda_0_estim_additive_AMOC3, sigma_estim_additive_AMOC3)
time_to_tipping_additive_AMOC3 <- tc_estim_additive_AMOC3 - max(AMOC_data$time) 



original_estim_additive_AMOC3 <- tibble(true_value = c(A_estim_additive_AMOC3, 
                                                       alpha_0_estim_additive_AMOC3, lambda_0_estim_additive_AMOC3,  m_estim_additive_AMOC3,
                                                       mu0_estim_additive_AMOC3, sigma_estim_additive_AMOC3, tau_estim_additive_AMOC3, tc_estim_additive_AMOC3))



#dynamic_part_init <- c(max(AMOC_data$time) - t_0, 1)
estim_matrix_additive_AMOC3 <- matrix(data = NA, nrow = numSim, ncol = 8)
#tc_vector <- numeric(numSim)
set.seed(07052024)
for (i in 1:numSim) {
  if (i %% 5 == 1) {
    cat("Currently grinding", i, "....\n")
  }
  
  success <- FALSE
  
  while (!success) {
    tryCatch({
      random_seed <- sample(100000, size = 1)
      OU_sim <- simulate_additive_noise_tipping_model(step_length = actual_dt,
                                                      par = sim_param_additive_AMOC3, tau = tau_estim_additive_AMOC3,
                                                      t_0 = T_0_additive_AMOC3,
                                                      beyond_tipping = -time_to_tipping_additive_AMOC3, seed = random_seed)
      
      # Stationary part
      sim_OU_stationary_estim <- optimize_stationary_likelihood(
        likelihood_fun = OU_likelihood,
        data = OU_sim$X_t[OU_sim$t < T_0_additive_AMOC3],
        init_par = c(alpha_0_estim_additive_AMOC3, mu0_estim_additive_AMOC3, sigma_estim_additive_AMOC3),
        delta = actual_dt,
        exp_sigma = TRUE,
        method = "BFGS"
      )$par
      
      # Dynamic part
      sim_OU_dynamic_estim <- optimize_dynamic_likelihood(
        likelihood_fun = OU_dynamic_likelihood,
        data = OU_sim$X_t[OU_sim$t >= T_0_additive_AMOC3],
        init_par = c(tau_estim_additive_AMOC3, A_estim_additive_AMOC3),
        delta = actual_dt,
        alpha0 = sim_OU_stationary_estim[1],
        mu0 = sim_OU_stationary_estim[2],
        sigma = sim_OU_stationary_estim[3],
        method = "BFGS",
        control = list(reltol = sqrt(.Machine$double.eps) / 1000)
      )$par
      
      estim_matrix_additive_AMOC3[i, 1:8] <- c(
        sim_OU_dynamic_estim[2], 
        sim_OU_stationary_estim[1],
        -sim_OU_stationary_estim[1]^2 / (4 * sim_OU_dynamic_estim[2]),
        sim_OU_stationary_estim[2] - sim_OU_stationary_estim[1] / (2 * sim_OU_dynamic_estim[2]),
        sim_OU_stationary_estim[2],
        sim_OU_stationary_estim[3],
        sim_OU_dynamic_estim[1],
        T_0_additive_AMOC3 + sim_OU_dynamic_estim[1] + 1870
      )
      
      success <- TRUE
    }, error = function(e) {
      cat("Error occurred at iteration", i, ":", conditionMessage(e), "\n")
    })
  }
  print(estim_matrix_additive_AMOC3[i, 8])
}


estim_tibble_additive_AMOC3 <- as_tibble(estim_matrix_additive_AMOC3)

names(estim_tibble_additive_AMOC3) <- c("A", "alpha_0", "lambda_0", "m", "mu0", "sigma", "tau", "tc")

if(!file.exists("data/original_estim_additive_AMOC3.csv")){
  utils::write.table(original_estim_additive_AMOC3, file="data/original_estim_additive_AMOC3.csv",
                     sep = ",", row.names = FALSE)
} else{
  original_estim_additive_AMOC3 <- read_csv("data/original_estim_additive_AMOC3.csv")
}

if(!file.exists("data/estim_tibble_additive_AMOC3.csv")){
  utils::write.table(estim_tibble_additive_AMOC3, file="data/estim_tibble_additive_AMOC3.csv",
                     sep = ",", row.names = FALSE)
} else{
  estim_tibble_additive_AMOC3 <- read_csv("data/estim_tibble_additive_AMOC3.csv")
}

original_estim_additive_AMOC3 <- original_estim_additive_AMOC3 |> mutate(Parameter = names(estim_tibble_additive_AMOC3))

estim_tibble_additive_AMOC3_long <- estim_tibble_additive_AMOC3 |> filter(A > 0, A < 5) |> 
  pivot_longer(cols = everything(), names_to = "Parameter", values_to = "value")

estim_tibble_additive_AMOC3_plot <- estim_tibble_additive_AMOC3_long |> 
  ggplot(aes(x = value, y = after_stat(density), fill = Parameter)) +
  geom_histogram(bins = 30, alpha = 2, col = "black", linewidth = .2) + 
  geom_vline(data = original_estim_additive_AMOC3, mapping = aes(xintercept = true_value), linewidth = 1, linetype = "dashed") +
  facet_wrap(~Parameter, scales = "free", ncol = 4) +
  labs(x = "Estimate", y = "") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) + 
  theme(strip.text.x = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        axis.title.x = element_text(face = "bold", size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)
  ) + 
  scale_fill_manual(values = thesis_palette,
                    labels = c("A", expression(alpha*phantom(.)[0]), expression(lambda*phantom(.)[0]), "m",
                               expression(mu*phantom(.)[0]), expression(sigma), expression(tau*phantom(.)[c], "Tipping year")))


# ggsave(estim_tibble_additive_AMOC3_plot, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "estim_tibble_additive_AMOC3_plot.jpeg",
#        height = 7, width = 15, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)

Estimates_original_paper <- readxl::read_excel("data/estim_matrix_AMOC_paper.xlsx")

# Compare fingerprint distribution of tipping year for both models
tau_distribution_all <- bind_rows(
  tibble(A = estim_tibble$A, tipping_year = estim_tibble$tc, Fingerprint = "AMOC2", Model = "t-diffusion"),
  tibble(A = filter(estim_tibble_AMOC1, A < 10, A > 0)$A,
         tipping_year = filter(estim_tibble_AMOC1, A < 10, A > 0)$tc,
         Fingerprint = "AMOC1", Model = "t-diffusion"),
  tibble(A = filter(estim_tibble_AMOC3, A < 10, A > 0)$A,
         tipping_year = filter(estim_tibble_AMOC3, A < 10, A > 0)$tc,
         Fingerprint = "AMOC3", Model = "t-diffusion"),
  tibble(A = estim_tibble_additive$A, tipping_year = estim_tibble_additive$tc,
         Fingerprint = "AMOC2", Model = "Additive"),
  tibble(A =  filter(estim_tibble_additive_AMOC1, A < 10, A > 0)$A,
         tipping_year = filter(estim_tibble_additive_AMOC1, A < 10, A > 0)$tc,
         Fingerprint = "AMOC1", Model = "Additive"),
  tibble(A = filter(estim_tibble_additive_AMOC3, A < 10, A > 0)$A, tipping_year = filter(estim_tibble_additive_AMOC3, A < 10, A > 0)$tc,
         Fingerprint = "AMOC3", Model = "Additive"),
  tibble(A = Estimates_original_paper$a, tipping_year = Estimates_original_paper$tc,
         Fingerprint = "AMOC2", Model = "Penalized additive")
)


distribution_tau_models <- tau_distribution_all |> ggplot(aes(x = tipping_year, y = Model, fill = Fingerprint)) + geom_violin() +
  scale_fill_manual(values = c(thesis_palette[1],thesis_palette[6], thesis_palette[3]))  +
  ylab("") + xlab("Tipping year") + 
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  theme(
    legend.title = element_blank(),
    strip.text = element_blank(),
    legend.key.size = unit(2, 'lines'),
    legend.text = element_text(face = "bold", size = 18),
    axis.text = element_text(face = "bold", size = 20),
    axis.title = element_text(face = "bold", size = 22),
  ) + 
  guides(col = guide_legend(override.aes = list(linewidth = 4)))

# ggsave(distribution_tau_models, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "distribution_tau_models.jpeg",
#        height = 7, width = 15, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)

tau_quantiles_table <- tau_distribution_all |> group_by(Fingerprint, Model) |> 
  reframe(quantiles = quantile(tipping_year, prob = c(0.025, 0.165, 0.5, 0.835, 0.975)),
          prob = c(0.025, 0.165, 0.5, 0.835, 0.975)) |> 
  pivot_wider(id_cols = c(Fingerprint, Model), values_from = quantiles, names_from = prob) 

xtable::xtable(tau_quantiles_table)

tau_survival_curve <- tau_distribution_all |> group_by(Fingerprint, Model) |> 
  reframe(quantiles = quantile(tipping_year, prob = seq(0, 1, by = 0.01)),
          prob = seq(0, 1, by = 0.01)) 

tau_estimates_models <- bind_rows(
  original_estim |> mutate(Fingerprint = "AMOC2", Model = "t-diffusion"),
  original_estim_AMOC1 |> mutate(Fingerprint = "AMOC1", Model = "t-diffusion"),
  original_estim_AMOC3 |> mutate(Fingerprint = "AMOC3", Model = "t-diffusion"),
  original_estim_additive |> mutate(Fingerprint = "AMOC2", Model = "Additive"),
  original_estim_additive_AMOC1 |> mutate(Fingerprint = "AMOC1", Model = "Additive"),
  original_estim_additive_AMOC3 |> mutate(Fingerprint = "AMOC3", Model = "Additive"),
  tibble(true_value = 2057, Parameter = "tc") |> mutate(Fingerprint = "AMOC2", Model = "Penalized additive")
) |> filter(Parameter == "tc")


surival_curve_first97.5 <- tau_survival_curve |>
  bind_rows(
    tau_survival_curve |> filter(prob == 1, Model != "Penalized additive") |> group_by(Model) |> slice_max(quantiles) |> 
      select(Model, quantiles, prob) |> crossing(tibble(Fingerprint = c("AMOC3", "AMOC2")))
  ) |> filter(prob < 0.975) |> 
  ggplot(aes(x = quantiles, y = 1 - prob, col = Fingerprint)) + 
  geom_step(linewidth = 1.5) + 
  geom_vline(data = tau_estimates_models, mapping = aes(xintercept = true_value, col = Fingerprint), linetype = "dashed",
             linewidth = 1.5) + 
  geom_vline(xintercept = 2024, linetype = "dotted", linewidth = 1.5) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
  facet_wrap(~Model, scales = "free_x") +
  scale_color_manual(values = c(thesis_palette[1],thesis_palette[6], thesis_palette[3])) +
  scale_y_continuous(labels = scales::percent) +
  theme(
    legend.title = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text =  element_text(face = "bold", size = 24),
    legend.text = element_text(face = "bold", size = 22),
    axis.text = element_text(face = "bold", size = 20),
    axis.title = element_text(face = "bold", size = 24),
    legend.position = "bottom"
  ) + 
  guides(col = guide_legend(override.aes = list(linewidth = 8))) + 
  xlab("Year") + ylab("Survival prob.")

# ggsave(surival_curve_first97.5, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "surival_curve_first97.5.jpeg",
#          height = 7, width = 15, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)


###------------------------------------------------------------------------###
###-----------------------------Discussion---------------------------------###
###------------------------------------------------------------------------###

true_param <- c(0.8, -1.5, -2, 0.3)
actual_dt <- 1/12 * 1/2
tau <- 100
t_0 <- 10
nus <- round(sort(unique(c(exp(seq(-1.5, 1.5, length.out = 7)), 1/exp(seq(-1.5, 1.5, length.out = 7))))), 2)
beyond_tippings <- -tau * c(0.9, 0.7, 0.5, 0.1)

simulate_and_label_nu_disussion <- function(nu, beyond_tipping, seed) {
  sim_result <- simulate_additive_noise_tipping_model(
    step_length = actual_dt, par = c(true_param,nu), tau = tau, t_0 = t_0, beyond_tipping = beyond_tipping, 
    seed = seed
  )

  as_tibble(sim_result) |>
    mutate(nu = nu, beyond_tipping = beyond_tipping)
}


combinations <- crossing(nu = nus, beyond_tipping = beyond_tippings)

# Apply the simulate_and_label function to each combination
all_simulations <- combinations |> 
  pmap_dfr(~simulate_and_label_nu_disussion(..1, ..2, seed = 200524))

beyond_tipping_shifted <- tibble(
  beyond_tipping = beyond_tippings,
  prev_beyond_tipping = lag(beyond_tippings) + tau + t_0
) |>
  filter(!is.na(prev_beyond_tipping))

# Merge the shifted data frame with the simulation data
all_simulations <- all_simulations  |> 
  left_join(beyond_tipping_shifted, by = "beyond_tipping")

# Plot
mu_simulations_discussion_plot <- all_simulations |> 
  ggplot(mapping = aes(x = t, y = X_t, col = factor(nu))) + 
  geom_step(linewidth = 0.75) + 
  geom_vline(aes(xintercept = prev_beyond_tipping)
             , col = "black ", linetype = "dashed", linewidth = 0.8) +
  #geom_line(aes(y = mu_t), linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = mu_t - 2 * sqrt(-lambda_t/true_param[1])), linetype = "dashed", linewidth = 1) + 
  facet_wrap(~factor(beyond_tipping), scales = "free_x") +
  ylim(c(min(all_simulations$mu_t - 2 * sqrt(-all_simulations$lambda_t/true_param[1])) - 0.1,
         max(all_simulations$X_t) + 0.01)) + 
  scale_color_manual(values = thesis_palette) +
  xlab("time (t)") + ylab(expression(X[t])) +
  theme(
    strip.text = element_blank(),
    legend.text = element_text(face = "bold", size = 16),
    axis.text = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 18)
  ) + 
  scale_x_continuous(breaks = scales::pretty_breaks()) + 
  guides(color = guide_legend(override.aes = list(linewidth = 5), title = expression(nu)))


# ggsave(mu_simulations_discussion_plot, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "mu_simulations_discussion_plot.jpeg",
#        height = 7, width = 11, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)

# Correlation between A and the tipping year
correlation_between_A_and_tau_plot <- tau_distribution_all |> filter(A < 100, tipping_year < 2300) |>
  mutate(Model = ifelse(Model == "Penalized additive", "Pen. additive", Model)) |> 
  ggplot(aes(x = A, y = tipping_year, col = Fingerprint)) + 
  geom_smooth(method = "lm", se = F, linewidth = 2) + 
  geom_point(alpha = .4) + 
  geom_hline(yintercept = 2024, linetype = "dashed", linewidth = 2) +
  facet_wrap(~Model, scales = "free") + 
  scale_color_manual(values = c(thesis_palette[1],thesis_palette[6], thesis_palette[3])) + 
  theme(
    legend.title = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text =  element_text(face = "bold", size = 24),
    legend.text = element_text(face = "bold", size = 22),
    axis.text = element_text(face = "bold", size = 20),
    axis.title = element_text(face = "bold", size = 24)
  ) + 
  guides(col = guide_legend(override.aes = list(size = 5, linewidth = 0, alpha = 1))) +
  ylab("Year") + xlab("A")

# ggsave(correlation_between_A_and_tau_plot, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "correlation_between_A_and_tau_plot.jpeg",
#        height = 7, width = 15, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)

# Estimate with the penalization in the model like in the paper
# AMOC2
paper_stationary <- optimize_stationary_likelihood(
  likelihood_fun = OU_likelihood,
  data = AMOC_data$AMOC2[AMOC_data$time <= 1924],
  init_par = OU_init_params(AMOC_data$AMOC2[AMOC_data$time <= 1924], delta = 1/12),
  delta = 1 / 12, 
  exp_sigma = TRUE
)

# Dynamic part 

dynamic_estim_penalty <- optimize_dynamic_likelihood(
  likelihood_fun = OU_dynamic_likelihood,
  data = AMOC_data$AMOC2[AMOC_data$time > 1924],
  init_par = c(100, 1),
  delta = 1 / 12,
  alpha0 = paper_stationary$par[1],
  mu0 = paper_stationary$par[2],
  sigma = paper_stationary$par[3],
  pen = 0.004)

print(dynamic_estim_penalty$par)

# Trace plots of optimization

optimization_path_for_plot <- list(parameters = list(), objectives = numeric())

OU_dynamic_likelihood_trace <- function(par, ...) {
  obj_value <- OU_dynamic_likelihood(par, ...)

  optimization_path_for_plot$parameters <<- append(optimization_path_for_plot$parameters, list(par))
  optimization_path_for_plot$objectives <<- c(optimization_path_for_plot$objectives, obj_value)
  
  obj_value
}

dynamic_estim_penalty <- optimize_dynamic_likelihood(
  likelihood_fun = OU_dynamic_likelihood_trace,
  data = AMOC_data$AMOC2[AMOC_data$time > 1924],
  init_par = c(100, 1),
  delta = 1 / 12,
  alpha0 = paper_stationary$par[1],
  mu0 = paper_stationary$par[2],
  sigma = paper_stationary$par[3],
  pen = 0.004,
  method = "Nelder-Mead")

trace_optimization_tibble <- do.call(rbind, optimization_path_for_plot$parameters) |> as_tibble() |> 
  mutate(n = row_number(), objective = optimization_path_for_plot$objectives) |> 
   rename(tau = V1, A = V2) |>
  mutate(tau = (tau - 100) / 100,
         A = (A - 1) / 1)

trace_optimization_long <- trace_optimization_tibble |> 
  pivot_longer(-c(objective, n), names_to = "Variable", values_to = "Estimate")

trace_plot <- trace_optimization_long |> 
  ggplot(aes(x = n, y = Estimate, col = Variable)) + geom_step(linewidth = 1.5) +
  facet_wrap(~Variable, scales = "free") +
  scale_y_continuous(labels = scales::percent) +
  scale_color_manual(labels = c("A", expression(tau*phantom(.)[c])), values = thesis_palette) + 
  theme(strip.text = element_blank(),
        legend.title = element_blank(),
        legend.key.height = unit(3, "lines"),
        legend.text = element_text(size = 28),
        axis.title = element_text(face = "bold", size = 20),
        axis.text = element_text(size = 18)) + 
  guides(color = guide_legend(override.aes = list(linewidth = 7.5))) +
  xlab("N") + ylab("Relative deviation from init. value")

# ggsave(trace_plot, path = paste0(getwd(), "/tex_files/figures"),
#        filename = "trace_plot.jpeg",
#        height = 7, width = 15, dpi = 300, units = "in", device = "jpeg",
#        limitsize = FALSE, scale = 1)

## showcase of score along with likelihood
true_param <- c(0.87, -1.51, -2.69, 1)
tau <- 150
t_0 <- 20

mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) *  sqrt(abs(true_param[3] / true_param[1]))
alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
stationary_part_true_param <- c(alpha0, mu0, true_param[4])

# Potentially different values of actual_dt
actual_dts <- c(0.05)

# Number of repetitions
n_reps <- 100

results <- list()
set.seed(220524)
for (dt in actual_dts) {
  print(dt)
  for (rep in 1:n_reps) {
    if((rep %% 5) == 0){print(rep)}
    seed = sample.int(1000000, 1)
    sim_res_add <- simulate_additive_noise_tipping_model(dt, true_param, tau, t_0)
    
    data <- sim_res_add$X_t[sim_res_add$t < t_0]
  
    mle_optim <- optimize_stationary_likelihood(OU_likelihood, data, stationary_part_true_param, dt, exp_sigma = FALSE)
    
    score_nleqslv <- nleqslv::nleqslv(x = stationary_part_true_param, fn = OU_score, data = data, delta = dt)$x
    
    mle_optim_score <- optimize_stationary_likelihood(OU_likelihood, data, stationary_part_true_param, dt, exp_sigma = FALSE, gr = OU_score)
    
    results[[paste(dt, rep)]] <- list(
      actual_dt = dt,
      rep = rep,
      mle_optim = mle_optim$par,
      score_nleqslv = score_nleqslv,
      mle_optim_score = mle_optim_score$par
    )
  }
}

# Convert results to data frame
results_df_OU <- do.call(rbind, lapply(names(results), function(key) {
  res <- results[[key]]
  data.frame(
    actual_dt = res$actual_dt,
    rep = res$rep,
    Method = rep(c("Likelihood", "Score", "Likelihood & score"), each = length(res$mle_optim)),
    Parameter = rep(c("alpha0", "mu0", "sigma"), times = 3),
    Estimate = c(res$mle_optim, res$score_nleqslv, res$mle_optim_score)
  )
}))

# Plotting the distribution of estimates
results_df_OU |> group_by(Method, Parameter) |> 
  reframe(quantile = quantile(Estimate, probs = c(0.025, 0.165, 0.5, 0.835, 0.975)),
          probs = c(0.025, 0.165, 0.5, 0.835, 0.975)) |> 
  pivot_wider(id_cols = c(Method, Parameter), values_from = quantile, names_from = probs) |> 
  xtable::xtable()