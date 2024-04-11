library(tidyverse)
ggplot2::theme_set(ggthemes::theme_base())
source("source_code/tipping_simulations.R")
source("source_code/model_fitting.R")
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
stationaryPart_Likelihood_tibble <-  as_tibble(stationaryPart_Likelihood_Optimization) %>% mutate(type = "Likelihood")
names(stationaryPart_Likelihood_tibble) <- c("beta", "mu", "sigma", "type")
stationaryPart_Score_tibble <-  as_tibble(stationaryPart_Score_root) %>% mutate(type = "Score")
names(stationaryPart_Score_tibble) <- c("beta", "mu", "sigma", "type")


bind_rows(stationaryPart_Likelihood_tibble, stationaryPart_Score_tibble) %>%
  pivot_longer(-type, names_to = "Estimate", values_to = "Value") %>% 
  ggplot(aes(y = Value, x = type, fill = type)) + 
  geom_boxplot() + 
  facet_wrap(~Estimate, scales = "free_y") + theme(
    axis.title.x = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.text.x = element_blank()
  )

# Dynamic plot
dynamicPart_Strang_tibble <-  as_tibble(dynamicPart_Likelihood_Strang) %>% mutate(type = "Strang")
names(dynamicPart_Strang_tibble) <- c("tau", "A", "type")
dynamicPart_Simulation_tibble <-  as_tibble(dynamicPart_Likelihood_Sim) %>% mutate(type = "Simulation")
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


