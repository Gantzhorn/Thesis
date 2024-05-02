library(tidyverse)

thesis_theme <- ggthemes::theme_base() +
  theme(panel.border = element_blank(),       
        plot.background = element_blank(),     
        panel.background = element_blank())

ggplot2::theme_set(thesis_theme)
source("source_code/tipping_simulations.R")
source("source_code/model_fitting.R")
source("source_code/model_diagnostics.R")

# Colorblind friendly plot palette
thesis_palette <- c("#E69F00", "#56B4E9", "#009E73", "#CCB400", "#0072B2", "#D55E00", "#CC79A7")

# Create plot of path for figures in "Saddle-node bifurcation and Tipping Point Estimation
# Choose parameters appropriate for any of the diffusions
true_param <- c(1.5, 0.2 , -0.4, 0.15)
actual_dt <- 1/12 * 1 / 10
tau <- 50
t_0 <- 25
methodPlotseed <- 210424 

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
)


ggplot(xs_all, aes(x = t, y = X_t, color = Model)) +
  geom_step(linewidth = 0.5) + 
  geom_hline(yintercept = true_param[2], linetype = "dashed") +
  facet_grid(sample_id ~ Model) +
  ylim(0, 1) +
  scale_color_manual(values = thesis_palette) +
  labs(y = expression(X[t])) + 
  theme(
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    legend.key.width = unit(2, "lines"),
    legend.key.height = unit(1, "lines")
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 5))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))

# Create plot of path for figures in "Saddle-node bifurcation and Tipping Point Estimation
# Choose parameters appropriate for all but jacobi diffusion
true_param <- c(0.2, 1.75, -2, 0.15)

actual_dt <- 0.005
tau <- 100
t_0 <- 50
methodPlotseed <- 210424 

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
)


ggplot(xs_all, aes(x = t, y = X_t, color = Model)) +
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
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  ylim(true_param[2] - 1, true_param[2] + sqrt(abs(true_param[3] / true_param[1])) + 1)

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
A <- 1
C <- -0.5

double_well <- function(x){
  -A * x^3 - C * x
}

double_well_deriv <- function(x){
  - 3 * A * x^2 - C
}

roots_neglambda <- polyroot(c(0, -C, 0, -A))

real_roots_neglambda <- roots_neglambda[dplyr::near(Im(roots_neglambda), 0)]
real_roots_neglambda <- sort(Re(real_roots_neglambda))

double_well_approx <- function(x){
  -3 * real_roots_neglambda[1] * (x - real_roots_neglambda[1])^2 - real_roots_neglambda[1] * C 
}

# Generate data
x_values <- seq(real_roots_neglambda[1] - 0.2, real_roots_neglambda[3] + 0.2, by = 0.01)

colors <- ifelse(double_well_deriv(real_roots_neglambda) < 0, "black", "white")

data <- data.frame(x = x_values, force = double_well(x_values), approx = double_well_approx(x_values))

data_pivot <- data |> 
  pivot_longer(cols = c(force, approx), names_to = "Type", values_to = "Value")

roots_data_neglambda <- data.frame(x = real_roots_neglambda, y = rep(0, length(real_roots_neglambda)), col = colors)

double_well_plot_neg <- ggplot(data, aes(x = x, y = force)) +
  geom_line(linewidth = 1) +
  labs(x = "x",
       y = expression(dot(x))) + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(data = roots_data_neglambda, aes(x = x, y = y, fill = col), size = 3, pch = 21) +
  scale_fill_manual(values = c("black", "white")) + 
  theme(legend.position = "none",
        axis.title.y = element_text(angle = 0, vjust = 0.5))

ggsave(path = paste0(getwd(), "/tex_files/figures"), filename = "double_well_plot_neg.jpeg")

# Double well potential plot positive lambda
# Parameters
A <- 1
C <- 0.5

roots_poslambda <- polyroot(c(0, -C, 0, -A))

roots_poslambda <- roots_poslambda[dplyr::near(Im(roots_poslambda), 0)]
roots_poslambda <- sort(Re(roots_poslambda))

double_well_approx <- function(x){
  -3 * roots_poslambda[1] * (x - roots_poslambda[1])^2 - roots_poslambda[1] * C 
}

# Generate data

colors <- ifelse(double_well_deriv(roots_poslambda) < 0, "black", "white")

data <- data.frame(x = x_values, force = double_well(x_values), approx = double_well_approx(x_values))

data_pivot <- data |> 
  pivot_longer(cols = c(force, approx), names_to = "Type", values_to = "Value")

roots_data_poslambda <- data.frame(x = roots_poslambda, y = rep(0, length(roots_poslambda)), col = colors)

double_well_plot_pos <- ggplot(data, aes(x = x, y = force)) +
  geom_line(linewidth = 1) +
  labs(x = "x",
       y = expression(dot(x))) + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(data = roots_data_poslambda, aes(x = x, y = y, fill = col), size = 3, pch = 21) +
  scale_fill_manual(values = c("black", "white")) + 
  theme(legend.position = "none",
        axis.title.y = element_text(angle = 0, vjust = 0.5))

ggsave(plot = double_well_plot_pos, 
       path = paste0(getwd(), "/tex_files/figures"),
       filename = "double_well_plot_pos.jpeg")

double_well_plot_pos <- double_well_plot_pos + ylab("") 

double_well_plot_combined <- gridExtra::grid.arrange(double_well_plot_neg, double_well_plot_pos, ncol = 2)

ggsave(filename = "double_well_plot_combined.jpeg", plot = double_well_plot_combined,
       path = paste0(getwd(), "/tex_files/figures"),
       width = 10,
       height = 5,
       dpi = 300)


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
              exp_sigma = TRUE)



# Dynamic part
dynamic_part_starting_param <- c(100, 1) + runif(2, min = -1) * c(30, 0.3)
#dynamic_part_starting_param <- c(83.619117, 1.202179)
dynamic_part_estim_param <- optimize_dynamic_likelihood(
                  likelihood_fun = t_transform_dynamic_likelihood,
                  data = asinh(AMOC_data$AMOC2[AMOC_data$time >= t_0]),
                  init_par = dynamic_part_starting_param,
                  delta = actual_dt,
                  alpha0 = stationary_part_estim_param[1],
                  mu0 = stationary_part_estim_param[2],
                  sigma = stationary_part_estim_param[3],
                  method = "BFGS",
                  control = list(reltol = sqrt(.Machine$double.eps) / 100000))

t_transform_dynamic_likelihood(par = c(100,1),
                               data = asinh(AMOC_data$AMOC2[AMOC_data$time >= t_0]),
                               delta = actual_dt,
                               alpha0 = stationary_part_estim_param[1],
                               mu0 = stationary_part_estim_param[2],
                               sigma = stationary_part_estim_param[3])

library(dplyr)
library(ggplot2)

# Define the parameter ranges
par1_values <- seq(95, 200, length.out = 100)  # from 90 to 200
par2_values <- seq(0.5, 1.3, length.out = 100)  # from 0.7 to 1.3

# Create a grid of parameter combinations
param_grid <- expand.grid(par1 = par1_values, par2 = par2_values)

# Function to apply
evaluate_function <- function(par1, par2) {
  t_transform_dynamic_likelihood(
    par = c(par1, par2),
    data = asinh(AMOC_data$AMOC2[AMOC_data$time >= t_0]),
    delta = actual_dt,
    alpha0 = stationary_part_estim_param[1],
    mu0 = stationary_part_estim_param[2],
    sigma = stationary_part_estim_param[3]
  )
}

# Apply the function to each row in the parameter grid
param_grid$result <- mapply(evaluate_function, param_grid$par1, param_grid$par2)

ggplot(param_grid, aes(x = par1, y = par2, fill = result)) +
  geom_tile() +
  labs(title = "Heatmap of t_transform_dynamic_likelihood",
       x = "Parameter 1 (par1)",
       y = "Parameter 2 (par2)",
       fill = "Function Value") +
  scale_fill_gradient(low = "blue", high = "red")

dynamic_part_estim_param
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

for (i in 1:numSim){
if(i %% 5 == 1){cat("Currently grinding", i, "....\n" )}
sim_t <- simulate_t_distribution_tipping_model(step_length = actual_dt, par = sim_param, tau = tau_estim,
                                      t_0 = T_0, beyond_tipping = -time_to_tipping, seed = i)


# Stationary part
sim_t_stationary_estim <-  optimize_stationary_likelihood(
  likelihood_fun = t_diffusion_strang_splitting,
  data = asinh(sim_t$X_t[sim_t$t < T_0]),
  init_par = stationary_part_starting_param,
  delta = actual_dt,
  exp_sigma = TRUE,
  method = "BFGS"
  )

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
  control = list(reltol = sqrt(.Machine$double.eps) / 1000))



estim_matrix[i, 1] <- sim_t_dynamic_estim[2]

estim_matrix[i, 2] <- sim_t_stationary_estim[1]
estim_matrix[i, 3] <- -sim_t_stationary_estim[1]^2 / (4 * sim_t_dynamic_estim[2])
estim_matrix[i, 4] <- sim_t_stationary_estim[2] - sim_t_stationary_estim[1] / (2 * sim_t_dynamic_estim[2])
estim_matrix[i, 5] <- sim_t_stationary_estim[2]
estim_matrix[i, 6] <- sim_t_stationary_estim[3]
estim_matrix[i, 7] <- sim_t_dynamic_estim[1]
estim_matrix[i, 8] <- T_0 + sim_t_dynamic_estim[1] + 1870 # Shift time back
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
