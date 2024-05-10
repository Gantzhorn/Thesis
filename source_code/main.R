library(tidyverse)


source("source_code/tipping_simulations.R")
source("source_code/model_fitting.R")
source("source_code/model_diagnostics.R")

thesis_theme <- ggthemes::theme_base() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),  # Adds border around each facet
    plot.background = element_blank(),       # No background for the plot
    panel.background = element_blank(),      # No background for the panels (facets)
    strip.background = element_blank()       # No background for the facet labels area
  )

ggplot2::theme_set(thesis_theme)

# Colorblind friendly plot palette
thesis_palette <- c("#E69F00", "#56B4E9", "#009E73", "#CCB400", "#0072B2",
                    "#D55E00", "#CC79A7", "#F0E442", "#D55E87", "#6E016B")



# Create plot of path for figures in "Saddle-node bifurcation and Tipping Point Estimation
# Choose parameters appropriate for any of the diffusions
true_param <- c(1.5, 0.4, -0.2, 0.1)
actual_dt <- 1/12 * 1 / 10
tau <- 50
t_0 <- 25


simulate_multiple_times <- function(simulate_function, n_times, ...) {
  simulate_res <- do.call(rbind,
          lapply(1:n_times, function(i) {
            #set.seed(methodPlotseed + i)
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


 # Experiment with different estimation methods on simulated data

# Additive noise model
true_param <- c(0.87, -1.51, -2.69, 0.2)
actual_dt <- 1/12
tau <- 50
t_0 <- 25

# For numerical likelihood
M <- 1 / actual_dt
N <- 10

sim_res_add <- simulate_additive_noise_tipping_model(actual_dt, true_param, tau, t_0, beyond_tipping = -5)
#sim_res_add <- filter(sim_res_add, t < t_0 + tau / 2)

## Parameters for stationary part
mu0 = true_param[2] + sqrt(abs(true_param[3]) / true_param[1])
alpha0 <- 2 * sqrt(true_param[1] * abs(true_param[3]))
stationary_part_true_param <- c(alpha0, mu0, true_param[4])
## Parameters for dynamic part
dynamic_part_true_param <- c(tau, true_param[1])

# Simulation setup
numSim <- 100
stationaryPart_Likelihood_Optimization <- matrix(NA, nrow = numSim, ncol = 3)
stationaryPart_Score_root <- matrix(NA, nrow = numSim, ncol = 3)

dynamicPart_Likelihood_Strang <- matrix(NA, nrow = numSim, ncol = 2)
dynamicPart_Likelihood_Sim <- matrix(NA, nrow = numSim, ncol = 2)



for (i in 1:numSim){
  cat("Grinding", i, "\n")
  sim_res_add <- simulate_additive_noise_tipping_model(actual_dt, true_param, tau, t_0, beyond_tipping = -15)
  # Stationary part
  stationaryPart_Likelihood_Optimization[i, ] <- (optimize_stationary_likelihood(likelihood_fun = OU_likelihood,
                                                  data = sim_res_add$X_t[sim_res_add$t < t_0],
                                                  init_par = stationary_part_true_param,
                                                  delta = actual_dt, exp_sigma = FALSE) -
                                                    stationary_part_true_param) / stationary_part_true_param

  stationaryPart_Score_root[i, ] <- (nleqslv::nleqslv(x = stationary_part_true_param,
                 fn = OU_Score,
                 data = sim_res_add$X_t[sim_res_add$t < t_0],
                 delta = actual_dt)$x - stationary_part_true_param) / stationary_part_true_param

  ## Dynamic part

  dynamicPart_Likelihood_Strang[i, ] <- (optimize_dynamic_likelihood(likelihood_fun = OU_dynamic_likelihood,
                            data = sim_res_add$X_t[sim_res_add$t > t_0],
                            init_par = dynamic_part_true_param,
                            delta = actual_dt,
                            alpha0 = stationary_part_true_param[1],
                            mu0 = stationary_part_true_param[2],
                            sigma = stationary_part_true_param[3], exp_sigma = TRUE) -
                            dynamic_part_true_param) / dynamic_part_true_param

  dynamicPart_Likelihood_Sim[i, ] <- (optimize_dynamic_simulation_likelihood(likelihood_fun = OU_dynamic_simulation_likelihood,
                                       data = sim_res_add$X_t[sim_res_add$t > t_0],
                                       times = sim_res_add$t[sim_res_add$t > t_0],
                                       M = M, N = N, init_par = dynamic_part_true_param,
                                       alpha0 = stationary_part_true_param[1],
                                       mu0 = stationary_part_true_param[2],
                                       sigma = stationary_part_true_param[3], t_0 = t_0) - 
                                       dynamic_part_true_param) / dynamic_part_true_param
}

# stationary plot
stationaryPart_Likelihood_tibble <-  as_tibble(stationaryPart_Likelihood_Optimization) |> mutate(type = "Likelihood")
names(stationaryPart_Likelihood_tibble) <- c("beta", "mu", "sigma", "type")
stationaryPart_Score_tibble <-  as_tibble(stationaryPart_Score_root) |> mutate(type = "Score")
names(stationaryPart_Score_tibble) <- c("beta", "mu", "sigma", "type")


bind_rows(stationaryPart_Likelihood_tibble, stationaryPart_Score_tibble) |>
  pivot_longer(-type, names_to = "Estimate", values_to = "Value") |> 
  ggplot(aes(y = Value, x = type, fill = type)) + 
  geom_boxplot() + 
  facet_wrap(~Estimate, scales = "free_y") + theme(
    axis.title.x = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.text.x = element_blank()
  )

# Dynamic plot
dynamicPart_Strang_tibble <-  as_tibble(dynamicPart_Likelihood_Strang) |> mutate(type = "Strang")
names(dynamicPart_Strang_tibble) <- c("tau", "A", "type")
dynamicPart_Simulation_tibble <-  as_tibble(dynamicPart_Likelihood_Sim) |> mutate(type = "Simulation")
names(dynamicPart_Simulation_tibble) <- c("tau", "A", "type")


bind_rows(dynamicPart_Strang_tibble, dynamicPart_Simulation_tibble) |> 
  pivot_longer(-type, names_to = "Estimate", values_to = "Value") |> 
  filter(Value < 600) |> 
  ggplot(aes(y = Value, x = type, fill = type)) + 
  geom_boxplot() + 
  facet_wrap(~Estimate, scales = "free_y") + theme(
    axis.title.x = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.text.x = element_blank()
  )

# Precision of estimator as a function of how early we have information up to.
# Try t-distribution
true_param <- c(0.3, -2, -3, 0.1)
actual_dt <- 0.1
tau <- 100
t_0 <- 50
mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) * sqrt(abs(true_param[3] / true_param[1]))
alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
stationary_part_true_param <- c(alpha0, mu0, true_param[4])
dynamic_part_true_param <- c(tau, true_param[1])

numSim <- 10

time_to_tipping <- seq(-50, 0, by = 5)

results <- list()
nu_param_values <- c(0.6, 0.8, 1, 1.2, 1.4)

for (nu in nu_param_values){
  true_param <- c(0.3, -2, -3, 0.1, nu)
  mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) * sqrt(abs(true_param[3] / true_param[1]))
  alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
  stationary_part_true_param <- c(alpha0, mu0, true_param[4])
  dynamic_part_true_param <- c(tau, true_param[1], nu)
  
  
  
  for (i in seq_along(time_to_tipping)){
    print(time_to_tipping[i])
    for (j in 1:numSim){
      if(j%%5 == 0 | j == 1){print(j)}
      sim_res_t_distribution <- simulate_t_distribution_tipping_model(
        actual_dt, true_param, tau,
        t_0, beyond_tipping = time_to_tipping[i])
      
    # Stationary part
    
      t_dist_estim_param <- optimize_stationary_likelihood(
        likelihood_fun = t_diffusion_strang_splitting,
        data = asinh(sim_res_t_distribution$X_t[sim_res_t_distribution$t < t_0]),
        init_par = stationary_part_true_param,
        delta = actual_dt,
        exp_sigma = TRUE)
    
    # Dynamic part
    
    accuracy_tau[j, i] <- (optimize_dynamic_likelihood(likelihood_fun = t_transform_dynamic_likelihood,
                                data = asinh(sim_res_t_distribution$X_t[sim_res_t_distribution$t > t_0]),
                                init_par = dynamic_part_true_param,
                                delta = actual_dt,
                                alpha0 = t_dist_estim_param[1],
                                mu0 = t_dist_estim_param[2],
                                sigma = t_dist_estim_param[3])[1] - tau) / tau
    }
  }
  # Calculate means and quantiles
  column_wise_means <- colMeans(accuracy_tau)
  column_wise_lower_quantiles <- apply(accuracy_tau, 2, function(x) quantile(x, probs = 0.1))
  column_wise_upper_quantiles <- apply(accuracy_tau, 2, function(x) quantile(x, probs = 0.9))
  
  # Store results
  results[[as.character(nu)]] <- tibble(
    x = time_to_tipping,
    mean = column_wise_means,
    lower = column_wise_lower_quantiles,
    upper = column_wise_upper_quantiles,
    nu = nu
  )
}

combined_data <- bind_rows(results, .id = "nu_label")

combined_data$nu <- as.factor(combined_data$nu)

ggplot(combined_data, aes(x = x, y = mean)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  facet_wrap(~nu, scales = "free_y") + 
  xlab("Observations until time") + ylab("Relative deviation from tipping time") +
  ggtitle("Relative deviation from tipping time by fifth parameter value")

### Plots for sample paths
true_param_t <- c(0.7, -2, -3, 0.15)
actual_dt <- 0.1
tau <- 100
t_0 <- 50
sim_res_t_distribution <- simulate_t_distribution_tipping_model(actual_dt, true_param_t, tau, t_0)
sim_res_t_distribution |> ggplot2::ggplot(ggplot2::aes(x = t, y = X_t)) +
  ggplot2::geom_step() + ggplot2::geom_hline(yintercept = true_param_t[2], linetype = "dashed") +
  ggplot2::geom_vline(xintercept = t_0)

true_param_sqrt <- c(0.7, 2, -3, 0.15)

sim_res_sqrt <- simulate_squareroot_noise_tipping_model(actual_dt, true_param_sqrt, tau, t_0)
sim_res_sqrt |> ggplot2::ggplot(ggplot2::aes(x = t, y = X_t)) +
  ggplot2::geom_step() + ggplot2::geom_hline(yintercept = true_param_sqrt[2], linetype = "dashed") +
  ggplot2::geom_vline(xintercept = t_0)

#nu plot 
# Parameters
lambda_0 <- -2 
tau <- 75
t_0 <- 10
nu_values <- c(0.05, 0.1, 0.5, 0.8, 1, 1.2, 1.7, 2.2, 3)

# Create a sequence of t values from 0 to 100
t_values <- seq(0, tau + t_0, length.out = 500)

# Initialize a data frame to hold the values
plot_data <- data.frame()

# Calculate the function values for each nu
for (nu in nu_values) {
  # Apply the function to each value of t
  lambda_t <- lambda_0 * (1 - ifelse(t_values > t_0, (t_values - t_0) / tau, 0))^nu
  
  # Append to the data frame
  plot_data <- rbind(plot_data, data.frame(t = t_values, lambda = lambda_t, nu = as.factor(nu)))
}

# Plotting using ggplot2
ggplot(plot_data, aes(x = t, y = lambda, color = nu, group = nu)) +
  geom_line() +
  labs(title = "Plot of the function for different values of nu",
       x = "Time (t)",
       y = "Function value",
       color = "Nu")




# For negative lambda
#A <- 1
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
  
  # Initialize a flag to track whether the computation was successful
  success <- FALSE
  
  # Retry the computation until successful or a maximum number of attempts is reached
  attempts <- 0
  max_attempts <- 3
  
  while (!success && attempts < max_attempts) {
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
