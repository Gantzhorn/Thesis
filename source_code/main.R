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


###------------------------------------------------------------------------###
###----------------------------Method section------------------------------###
###------------------------------------------------------------------------###
# Create plot of path for figures in "Saddle-node bifurcation and Tipping Point Estimation
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
  geom_line(data = fixed_point_lines, aes(x = t, y = lower), linewidth = 0.75, linetype = "dashed", color = "grey25") +
  geom_step(linewidth = 0.5) + 
  geom_hline(yintercept = true_param[2], linetype = "dashed") +
  facet_grid(sample_id ~ Model) +
  ylim(0, 1) +
  scale_color_manual(values = thesis_palette) +
  labs(y = expression(X[t])) + 
  guides(color = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) + 
  theme(legend.text=element_text(size=18),
        legend.title = element_text(size = 22),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        legend.key.width = unit(2, "lines"),
        legend.key.height = unit(1, "lines"))

ggsave(sample_paths_plot_small_scale, path = paste0(getwd(), "/tex_files/figures"),
       filename = "sample_paths_plot_small_scale.jpeg",
       height = 8, width = 13, dpi = 300, units = "in", device = "jpeg",
       limitsize = FALSE, scale = 1)

# Create plot of path for figures in "Saddle-node bifurcation and Tipping Point Estimation
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
  theme(
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    legend.key.width = unit(2, "lines"),
    legend.key.height = unit(1, "lines")
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) + 
  ylim(min(fixed_point_lines_big$lower)-0.1, max(fixed_point_lines_big$upper)+1)

ggsave(sample_paths_plot_big_scale, path = paste0(getwd(), "/tex_files/figures"),
       filename = "sample_paths_plot_big_scale.jpeg",
       height = 8, width = 13, dpi = 300, units = "in", device = "jpeg",
       limitsize = FALSE, scale = 1)


#nu plot 
# Parameters
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
  geom_line(linewidth = 1) +
  geom_hline(yintercept = m_nu_plot, linetype = "dashed") +
  labs(x = "Time (t)",
       y = expression(mu),
       color = expression(nu)) + 
  scale_color_manual(values = thesis_palette) +
  scale_linetype(guide = "none") + 
  guides(color = guide_legend(override.aes = list(linewidth = 5))) +
  theme(plot.margin = unit(c(0, 0, 0, 0.5), "cm"))

ggsave(nu_plot, path = paste0(getwd(), "/tex_files/figures"), filename = "nu_plot.jpeg",
       height = 6, width = 10, dpi = 300, units = "in", device = "jpeg",
       limitsize = FALSE, scale = 1)

# Potential, flow and bifurcation plots
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
}) %>% bind_rows()

# Find roots for each lambda parameter
roots_neglambda <- purrr::map2(lambda_double_well_vec, lambda_double_well_vec, ~ {
  roots <- polyroot(c(-.y, 1, 0, -1))
  real_roots <- roots[dplyr::near(Im(roots), 0)]
  data.frame(lambda = .x, root = sort(Re(real_roots)))  # Enrich each set of roots with lambda information
})

# Combine the list into a single data frame
roots_data_neglambda <- do.call(rbind, roots_neglambda) 

roots_data_neglambda <- roots_data_neglambda %>%
  mutate(
    root_color = case_when(
      #abs(double_well_doublederiv(root)) < 0.1 ~ "grey",
      double_well_doublederiv(root) < 0 ~ "black",
      TRUE ~ "white"
    ),
    lambda = factor(lambda)
  ) %>%
  as_tibble()

circleFun <- function(center=c(0,0), diameter=1, npoints=100, start=0, end=2)
{
  tt <- seq(start*pi, end*pi, length.out=npoints)
  tibble(x = center[1] + diameter / 2 * cos(tt),
         y = center[2] + diameter / 2 * sin(tt))
}
dat_half_circle <- circleFun(c(0.577, 0), 0.1, start = 1/2, end = 3/2)
dat_half_circle <- rbind(dat_half_circle, dat_half_circle[1,])  # Ensure closure by repeating the first point at the end
dat_half_circle$lambda <- factor(lambda_double_well_vec[3])  # Adjust as per your specific lambda settings


double_well_plot_neg <- ggplot(double_well_plot_data, aes(x = x_prime, y = force)) +
  geom_line(aes(col = lambda), linewidth = 1) +
  geom_hline(yintercept = 0) +
  geom_point(data = roots_data_neglambda, aes(x = root, y = 0, fill = root_color), size = 3, pch = 21) +
  scale_fill_manual(values = c("black", "white"), guide = "none") + 
  geom_polygon(data = dat_half_circle, aes(x, y), inherit.aes = FALSE, fill = "white", color = "black") +
  labs(x = "x", y = expression(dot(x))) +
  facet_grid(. ~ lambda) + 
  scale_color_manual(values = thesis_palette, labels = c("0", "0.25", expression(lambda * phantom()[" c"]))) +
  theme(strip.text = element_blank(),
        legend.text = element_text(size = 14)) +
  guides(color = guide_legend(title = expression(lambda))) + 
  coord_fixed()



ggsave(double_well_plot_neg, path = paste0(getwd(), "/tex_files/figures"), filename = "double_well_plot_neg.jpeg",
       height = 8, width = 10, dpi = 300, units = "in", device = "jpeg",
       limitsize = FALSE, scale = 1)

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
  theme(axis.title.y = element_text(size = 18),
        axis.title.x = element_text(margin = margin(t = 20), size = 18))



ggsave(bifurcation_diagram, path = paste0(getwd(), "/tex_files/figures"), filename = "bifurcation_diagram.jpeg",
       height = 6, width = 10, dpi = 300, units = "in", device = "jpeg",
       limitsize = FALSE, scale = 1)



###------------------------------------------------------------------------###
###------------------------Simulation studies------------------------------###
###------------------------------------------------------------------------###

# General study of performance
# Parameters where all work
true_param <- c(1.5, 0.4, -0.2, 0.1)
mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) * sqrt(abs(true_param[3] / true_param[1]))
alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
stationary_part_true_param <- c(alpha0, mu0, true_param[4])
actual_dt <- 1/12
tau <- 10
t_0 <- 30
numSim <- 100

OU_likelihood_estim <- matrix(NA, nrow = numSim, ncol = 3)

OU_score_estim <- matrix(NA, nrow = numSim, ncol = 3)

sqrt_martingale_estim <- matrix(NA, nrow = numSim, ncol = 3)

sqrt_strang_estim <- matrix(NA, nrow = numSim, ncol = 3)

linear_martingale_estim <- matrix(NA, nrow = numSim, ncol = 3)

linear_alt_strang_estim <- matrix(NA, nrow = numSim, ncol = 3)

linear_strang_estim <- matrix(NA, nrow = numSim, ncol = 3)

t_strang_estim <- matrix(NA, nrow = numSim, ncol = 3)

F_strang_estim <- matrix(NA, nrow = numSim, ncol = 3)

Jacobi_strang_estim <- matrix(NA, nrow = numSim, ncol = 3)


set.seed(110524)
for (i in 1:numSim){
  cat("Grinding iteration", i, "\n")
  success <- FALSE
  while(!success){
  tryCatch({
    
  seed <- sample.int(100000, 1)
  random_noise_start_value <- runif(3, min = .8, 1.2)
  
  # OU model
  sim_res_add <- simulate_additive_noise_tipping_model(step_length = actual_dt,
                                        par = true_param,
                                        tau = tau, t_0 = t_0,
                                        beyond_tipping = 0, seed = seed)
  
  OU_likelihood_estim[i, ] <- abs(optimize_stationary_likelihood(likelihood_fun = OU_likelihood,
                                 data = sim_res_add$X_t[sim_res_add$t < t_0],
                                 init_par = stationary_part_true_param * random_noise_start_value,
                                 delta = actual_dt, exp_sigma = FALSE,
                                 control = list(reltol = sqrt(.Machine$double.eps)/10))$par -
                                   stationary_part_true_param) /stationary_part_true_param
  
  OU_score_estim[i, ] <- abs(nleqslv::nleqslv(x = stationary_part_true_param * random_noise_start_value,
                   fn = OU_score,
                   data = sim_res_add$X_t[sim_res_add$t < t_0],
                   delta = actual_dt)$x -
                     stationary_part_true_param) / stationary_part_true_param
  
  # Sqrt-model
  sim_res_sqrt <- simulate_squareroot_noise_tipping_model(step_length = actual_dt,
                                          par = true_param, tau = tau, t_0 = t_0,
                                          beyond_tipping = 0, seed = seed)
  
  sqrt_martingale_estim[i, ] <- abs(CIR_quadratic_martingale(sim_res_sqrt$X_t[sim_res_sqrt$t < t_0], actual_dt) -
                                      stationary_part_true_param) /stationary_part_true_param

  sqrt_strang_estim[i, ] <- abs(optimize_stationary_likelihood(likelihood_fun = CIR_strang_splitting,
                                 data = 2 * sqrt(sim_res_sqrt$X_t[sim_res_sqrt$t < t_0]),
                                 init_par = stationary_part_true_param * random_noise_start_value,
                                 delta = actual_dt, exp_sigma = TRUE,
                                 control = list(reltol = sqrt(.Machine$double.eps)/10))$par -
                                  stationary_part_true_param) / stationary_part_true_param
  
  # Linear noise model
  sim_res_linear <- simulate_linear_noise_tipping_model(step_length = actual_dt,
                                      par = true_param, tau = tau, t_0 = t_0,
                                      beyond_tipping = 0, seed = seed)
  
  linear_martingale_estim[i, ] <- abs(nleqslv::nleqslv(x = stationary_part_true_param * random_noise_start_value,
                   fn = mean_reverting_GBM_martingale,
                   data = sim_res_linear$X_t[sim_res_linear$t < t_0],
                   delta = actual_dt)$x -
                     stationary_part_true_param) /stationary_part_true_param

  linear_strang_estim[i, ] <- abs(optimize_stationary_likelihood(mean_reverting_GBM_strang, 
                                 log(sim_res_linear$X_t[sim_res_linear$t<t_0]),
                                 init_par = stationary_part_true_param * random_noise_start_value,
                                 delta = actual_dt,
                                 exp_sigma = FALSE,
                                 control = list(reltol = sqrt(.Machine$double.eps)/10))$par -
    stationary_part_true_param) /stationary_part_true_param

  linear_alt_strang_estim[i, ] <- abs(optimize_stationary_likelihood(mean_reverting_GBM_alt_strang,
                                sim_res_linear$X_t[sim_res_linear$t<t_0],
                                init_par = stationary_part_true_param * random_noise_start_value,
                                delta = actual_dt,
                                exp_sigma = FALSE,
                                control = list(reltol = sqrt(.Machine$double.eps)/10))$par -
                                  stationary_part_true_param) / stationary_part_true_param
  
  
  # t-diffusion model
  sim_res_t_distribution <- simulate_t_distribution_tipping_model(step_length = actual_dt,
                                        par = true_param, tau = tau, t_0 = t_0,
                                        beyond_tipping = 0, seed = seed)
  
  t_strang_estim[i, ] <- abs(optimize_stationary_likelihood(
                likelihood_fun = t_diffusion_strang_splitting,
                data = asinh(sim_res_t_distribution$X_t[sim_res_t_distribution$t < t_0]),
                init_par = stationary_part_true_param * random_noise_start_value,
                delta = actual_dt,
                exp_sigma = TRUE,
                control = list(reltol = sqrt(.Machine$double.eps)/10))$par -
                  stationary_part_true_param) / stationary_part_true_param
  
  # F-diffusion model
  F_sim_dynamic <- simulate_F_distribution_tipping_model(step_length = actual_dt,
                                        par = true_param, tau = tau, t_0 = t_0,
                                        beyond_tipping = 0, seed = seed)
  
 F_strang_estim[i, ] <-  abs(optimize_stationary_likelihood(
    likelihood_fun = F_diffusion_strang_splitting,
    data = 2 * asinh(sqrt(F_sim_dynamic$X_t[F_sim_dynamic$t<t_0])),
    init_par = stationary_part_true_param * random_noise_start_value,
    delta = actual_dt,
    exp_sigma = TRUE,
    control = list(reltol = sqrt(.Machine$double.eps)/10))$par -
      stationary_part_true_param) / stationary_part_true_param
  
  # Jacobi model
  sim_res_jacobi <- simulate_jacobi_diffusion_tipping_model(step_length = actual_dt,
                                          par = true_param, tau = tau, t_0 = t_0,
                                          beyond_tipping = 0, seed = seed)

  Jacobi_strang_estim[i, ] <- abs(optimize_stationary_likelihood(
    likelihood_fun = jacobi_diffusion_strang_splitting,
    data = 2 * asin(sqrt(sim_res_jacobi$X_t[sim_res_jacobi$t<t_0])),
    init_par = stationary_part_true_param * random_noise_start_value,
    delta = actual_dt,
    exp_sigma = TRUE, 
    control = list(reltol = sqrt(.Machine$double.eps)/10))$par - 
      stationary_part_true_param) / stationary_part_true_param
  
  success <- TRUE
  }, error = function(e) {
    cat("Error occurred at iteration", i, ":", conditionMessage(e), "\n")
    
    Sys.sleep(1)
  })
  }
}

OU_likelihood_tibble <- as_tibble(OU_likelihood_estim) |>
  mutate(Model = "Additive", Type = "Likelihood", Method = "MLE")

OU_score_tibble <- as_tibble(OU_score_estim) |>
  mutate(Model = "Additive", Type = "Estimation equation", Method = "Score")

sqrt_estimation_equation_tibble <- as_tibble(sqrt_martingale_estim) |>
  mutate(Model = "Square-root", Type = "Estimation equation", Method = "Martingale")

sqrt_Strang_tibble <- as_tibble(sqrt_strang_estim) |>
  mutate(Model = "Square-root", Type = "Likelihood", Method = "Strang")

Linear_estimation_equation_tibble <- as_tibble(linear_martingale_estim) |>
  mutate(Model = "Linear", Type = "Estimation equation", Method = "Martingale")

Linear_alt_Strang_tibble <- as_tibble(linear_alt_strang_estim) |>
  mutate(Model = "Linear", Type = "Likelihood", Method = "Strang (Alternative)")

Linear_Strang_tibble <- as_tibble(linear_strang_estim)  |>
  mutate(Model = "Linear", Type = "Likelihood", Method = "Strang")

t_diffusion_Strang_tibble <- as_tibble(t_strang_estim)  |>
  mutate(Model = "t-diffusion", Type = "Likelihood", Method = "Strang")

F_diffusion_Strang_tibble <- as_tibble(F_strang_estim)  |>
  mutate(Model = "F-diffusion", Type = "Likelihood", Method = "Strang")

Jacobi_diffusion_Strang_tibble <- as_tibble(Jacobi_strang_estim) |>
  mutate(Model = "Jacobi-diffusion", Type = "Likelihood", Method = "Strang")

Stationary_estimation_all <- bind_rows(OU_likelihood_tibble, OU_score_tibble, sqrt_estimation_equation_tibble,
          sqrt_Strang_tibble, Linear_estimation_equation_tibble, 
          Linear_alt_Strang_tibble, Linear_Strang_tibble, t_diffusion_Strang_tibble,
          F_diffusion_Strang_tibble, Jacobi_diffusion_Strang_tibble) |> 
  rename(alpha0 = V1, mu0 = V2, sigma = V3) |>  
  mutate(Model = factor(Model), Type = factor(Type), Method = factor(Method)) |> 
  pivot_longer(-c(Model, Type, Method), names_to = "Parameter", values_to = "Estimate")

Stationary_estimation_all |> filter(Estimate < 5, Model %in% c("Additive", "Linear", "Square-root")) |> 
  ggplot(aes(x = Parameter, y = Estimate, fill = Method)) + geom_violin() +
  facet_wrap(~Model) + scale_y_log10() + scale_fill_manual(values = thesis_palette[3:8])

Stationary_estimation_all |> filter(Estimate < 5,!( Model%in% c("Additive", "Linear", "Square-root"))) |> 
  ggplot(aes(x = Parameter, y = Estimate, fill= Model)) + geom_violin() +
  scale_y_log10() + scale_fill_manual(values = thesis_palette[3:8])



# Precision of estimator as a function of how early we have information up to.
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
          likelihood_fun = t_diffusion_strang_splitting,
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
                       data  =sim_res_additive$X_t[sim_res_additive$t > t_0],
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


###------------------------------------------------------------------------###
###-----------------------------Result-------------------------------------###
###------------------------------------------------------------------------###

### Estimation on the AMOC

AMOC_data <- readr::read_table("data/AMOCdata.txt")

AMOC_data <- dplyr::rename_all(AMOC_data, ~ gsub('"', '',.))

AMOC_data <- mutate(AMOC_data, AMOC3 = AMOC0 - 3 * GM)

AMOC_data |> select(-GM) |>  pivot_longer(cols = -time, names_to = "AMOC_type", values_to = "Value") |> 
  ggplot(aes(x = time, y = Value, col = AMOC_type)) + geom_step() + facet_wrap(~AMOC_type) +
  scale_color_manual(values = thesis_palette)

actual_dt <- 1 / 12 # Observations every month
t_0 <- 1924 # Use same baseline data as paper.

# Stationary part
stationary_part_starting_param <- c(1, 1, 1)


stationary_part_estim_param <-  optimize_stationary_likelihood(
              likelihood_fun = t_diffusion_strang_splitting,
              data = asinh(AMOC_data$AMOC2[AMOC_data$time < t_0]),
              init_par = stationary_part_starting_param,
              delta = actual_dt,
              exp_sigma = TRUE)$par

# Dynamic part
dynamic_part_starting_param <- c(100, 1) 
dynamic_part_estim_param <- optimize_dynamic_likelihood(
                  likelihood_fun = t_transform_dynamic_likelihood,
                  data = asinh(AMOC_data$AMOC2[AMOC_data$time >= t_0]),
                  init_par = dynamic_part_starting_param,
                  delta = actual_dt,
                  alpha0 = stationary_part_estim_param[1],
                  mu0 = stationary_part_estim_param[2],
                  sigma = stationary_part_estim_param[3],
                  method = "BFGS",
                  control = list(reltol = sqrt(.Machine$double.eps) / 100000))$par

# Simulate from the model to construct parametric bootstrap confidence intervals
numSim <- 1000

# Define parameters - note that we shift the year to index 0.
T_0 <- t_0 - 1870
tau_estim <- dynamic_part_estim_param[1]
A_estim <- dynamic_part_estim_param[2]
alpha_0_estim <- stationary_part_estim_param[1]
mu0_estim <- stationary_part_estim_param[2]
sigma_estim <- stationary_part_estim_param[3]

lambda_0_estim <- -alpha_0_estim^2 / (4 * A_estim)
m_estim <- mu0_estim - alpha_0_estim / (2 * A_estim)

tc_estim <- tau_estim + t_0
sim_param <- c(A_estim, m_estim, lambda_0_estim, sigma_estim)
time_to_tipping <- tc_estim - max(AMOC_data$time)

original_estim <- tibble(true_value = c(A_estim, alpha_0_estim, lambda_0_estim,  m_estim, mu0_estim, sigma_estim, tau_estim, tc_estim))

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
        init_par = stationary_part_starting_param,
        delta = actual_dt,
        exp_sigma = TRUE,
        method = "BFGS"
      )$par
      
      # Dynamic part
      sim_t_dynamic_estim <- optimize_dynamic_likelihood(
        likelihood_fun = t_transform_dynamic_likelihood,
        data = asinh(sim_t$X_t[sim_t$t >= T_0]),
        init_par = dynamic_part_starting_param,
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
      
      attempts <- attempts + 1
      
      Sys.sleep(1)
    })
  }
  print(estim_matrix[i, 8])
}


estim_tibble <- as_tibble(estim_matrix)

names(estim_tibble) <- c("A", "alpha_0", "lambda_0", "m", "mu0", "sigma", "tau", "tc")

if(!file.exists("data/estim_tibble.csv")){
  utils::write.table(estim_tibble, file="data/estim_tibble.csv", sep = ",", row.names = FALSE)
} else{
  estim_tibble <- read_csv("data/estim_tibble.csv")
}

original_estim <- original_estim |> mutate(parameter = names(estim_tibble))

estim_tibble_long <- estim_tibble %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")

estim_tibble_plot <- estim_tibble_long |> ggplot(aes(x = value)) +
  geom_histogram(bins = 30, fill = thesis_palette[5], alpha = 0.7, col = "black") + 
  geom_vline(data = original_estim, mapping = aes(xintercept = true_value), linewidth = 1) +
  facet_wrap(~parameter, scales = "free_x", ncol = 4) +
  labs(x = "Estimate", y = "Count") +
  theme(strip.text.x = element_text(size = 12, face = "bold"),
        panel.spacing = unit(1.5, "lines"))

ggplot2::ggsave(filename = "tex_files/figures/estim_tibble_plot.jpeg",
                plot = estim_tibble_plot,
                width = 15,
                height = 7.5,
                dpi = 300)

  
qqplot_before_t0 <- tibble::tibble(obsSample =
                 t_diffusion_strang_splitting_resid(stationary_part_estim_param,
                 data = asinh(AMOC_data$AMOC2[AMOC_data$time < t_0]),
                 delta = actual_dt)) |>
                 ggplot2::ggplot(ggplot2::aes(sample = obsSample)) +
                 ggplot2::geom_qq() + ggplot2::geom_qq_line() + 
                 xlab("Theoretical Quantiles") + ylab("Sample Quantiles")  


qqplot_after_t0 <-  tibble::tibble(obsSample =
                 t_transform_dynamic_likelihood_resid(dynamic_part_estim_param,
                 data = asinh(AMOC_data$AMOC2[AMOC_data$time >= t_0]),
                 delta = actual_dt,
                 alpha0 = stationary_part_estim_param[1],
                 mu0 = stationary_part_estim_param[2],
                 sigma = stationary_part_estim_param[3])) |>
                 ggplot2::ggplot(ggplot2::aes(sample = obsSample)) +
                 ggplot2::geom_qq() + ggplot2::geom_qq_line() + 
                 xlab("Theoretical Quantiles") + ylab("Sample Quantiles")

qqplot_combined <- gridExtra::grid.arrange(qqplot_before_t0, qqplot_after_t0, ncol = 2)


ggplot2::ggsave(filename = "tex_files/figures/qqplot_combined.jpeg",
                plot = qqplot_combined,
                width = 15,
                height = 7.5,
                dpi = 300)
