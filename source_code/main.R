library(tidyverse)
ggplot2::theme_set(ggthemes::theme_base())
source("source_code/tipping_simulations.R")
source("source_code/model_fitting.R")
# Experiment with different estimation methods on simulated data

# Additive noise model
true_param <- c(0.87, -1.51, -2.69, 0.2)
actual_dt <- 1/12
tau <- 100
t_0 <- 25

# For numerical likelihood
M <- 3
N <- 5

sim_res_add <- simulate_additive_noise_tipping_model(actual_dt, true_param, tau, t_0)

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
  sim_res_add <- simulate_additive_noise_tipping_model(actual_dt, true_param, tau, t_0)
  # Stationary part
  stationaryPart_Likelihood_Optimization[i, ] <- (optimize_stationary_likelihood(likelihood_fun = OU_likelihood,
                                                  data = sim_res_add$X_weak_2.0[sim_res_add$t < t_0],
                                                  init_par = stationary_part_true_param,
                                                  delta = actual_dt, exp_sigma = FALSE) -
                                                    stationary_part_true_param) / stationary_part_true_param

  stationaryPart_Score_root[i, ] <- (nleqslv::nleqslv(x = stationary_part_true_param,
                 fn = OU_Score,
                 data = sim_res_add$X_weak_2.0[sim_res_add$t < t_0],
                 delta = actual_dt)$x - stationary_part_true_param) / stationary_part_true_param

  ## Dynamic part

  dynamicPart_Likelihood_Strang[i, ] <- (optimize_dynamic_likelihood(likelihood_fun = OU_dynamic_likelihood,
                            data = sim_res_add$X_weak_2.0[sim_res_add$t > t_0],
                            init_par = dynamic_part_true_param,
                            delta = actual_dt,
                            alpha0 = stationary_part_true_param[1],
                            mu0 = stationary_part_true_param[2],
                            sigma = stationary_part_true_param[3], exp_sigma = TRUE) -
                            dynamic_part_true_param) / dynamic_part_true_param

  dynamicPart_Likelihood_Sim[i, ] <- (optimize_dynamic_simulation_likelihood(likelihood_fun = OU_dynamic_simulation_likelihood,
                                       data = sim_res_add$X_weak_2.0[sim_res_add$t > t_0],
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
