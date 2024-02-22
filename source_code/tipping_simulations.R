# Title: Simulation of Stochastic Differential Equations for Tipping Points Analysis
# Author: Anders Gantzhorn Kristensen (University of Copenhagen, andersgantzhorn@gmail.com)
# Date: 2024-01-31 (Last Updated: 2024-02-22)
#-----------------------------------------------------------------------------------------------------------------------------#
# Project: Tipping Point Estimation in Ecological Systems using Stochastic Differential Equations
# Description: This script implements methods for simulating different types of Stochastic differential equation to model
# tipping points. By means of the scalar weak order 2.0 Itô-taylor method it aims to accurately simulate.
# sample paths from the rather complex bifurcation tipping point model with various noise structures. 
#-----------------------------------------------------------------------------------------------------------------------------#
# License: MIT License (for more information, see LICENSE file in the repository).
# Dependencies: ggplot2, ggthemes, tibble - All part of the tidyverse.
# Usage: See example usage in script below.
# Acknowledgements: For algorithm see: Applied Stochastic Differential Equations Särkkä, Simo and Solin, Arno - Algorithm 8.5.
# Keywords: Stochastic Differential Equations, Tipping Points, Simulation, Itô-taylor weak order 2.0
#-----------------------------------------------------------------------------------------------------------------------------#
ggplot2::theme_set(ggthemes::theme_base())

simulate_additive_noise_tipping_model <- function(step_length, param,
                                                  tau = 100, t_0 = 10,
                                                  X_0 = NA, beyond_tipping = 0){
  A <- param[1]
  m <- param[2]
  lambda_0 <- param[3]
  sigma <- param[4]
  
  if(is.na(X_0)){
    X_0 <- m + sqrt(abs(lambda_0) / A)
  }
  
  #Tipping point with added time if specified
  total_time <- tau + t_0 + beyond_tipping 
  
  N <- as.integer((total_time) / step_length) 
  
  # Initialize the process
  X_weak_2.0 <- numeric(N + 1) # N + 1 to include the initial point
  X_weak_2.0[1] <- X_0
  
  dW <- rnorm(N, mean = 0, sd = sqrt(step_length))
  time <- step_length * 0:N
  lambda_t <- numeric(N + 1)
  lambda_t <- lambda_0 * (1 - (time >= t_0) * (time - t_0) / tau)
  alpha_t <- 2 * sqrt(A * abs(lambda_t))
  mu_t <- m + sqrt(abs(lambda_t)/A)
  
  for(i in (2:(N+1))){
    X_weak_2.0[i] <- X_weak_2.0[i - 1] - (A*(X_weak_2.0[i - 1] - m)^2 + lambda_t[i - 1]) * step_length +
      sigma * dW[i - 1] - 2 * A *(X_weak_2.0[i - 1] - m) * sigma * (1 / 2 * dW[i - 1] * step_length) + 
      1 / 2 * (2 * A * (X_weak_2.0[i - 1] - m) * (A * (X_weak_2.0[i - 1] - m)^2 + lambda_t[i - 1]) - A * sigma^2) * step_length^2
    
  }
  tibble::tibble(t = time,
                 X_weak_2.0,
                 lambda_t, alpha_t, mu_t)
}

## Additive noise test simulations
# true_param <- c(0.015, 300, -2, 1)
# actual_dt <- 0.001
# tau <- 100
# t_0 <- 50
# sim_res_add <- simulate_additive_noise_tipping_model(actual_dt, true_param, tau, t_0)
# sim_res_add |> ggplot2::ggplot(ggplot2::aes(x = t, y = X_weak_2.0)) +
#   ggplot2::geom_step() + ggplot2::geom_hline(yintercept = true_param[2], linetype = "dashed") +
#   ggplot2::geom_vline(xintercept = t_0)
# 
# sim_res_add |>  ggplot2::ggplot(ggplot2::aes(x = t, y = lambda_t)) +
#   ggplot2::geom_step() + ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
#   ggplot2::geom_vline(xintercept = tau + t_0, linetype = "dashed")  +
#   ggplot2::geom_vline(xintercept = t_0, linetype = "dashed")

simulate_squareroot_noise_tipping_model<- function(step_length, param,
                                                   tau = 100, t_0 = 10,
                                                   X_0 = NA, beyond_tipping = 0){
  A <- param[1]
  m <- param[2]
  lambda_0 <- param[3]
  sigma <- param[4]

  if(is.na(X_0)){
    X_0 <- m + sqrt(abs(lambda_0) / A)
  }
  
  #Tipping point with added time if specified
  total_time <- tau + t_0 + beyond_tipping 
  
  N <- as.integer((total_time) / step_length) 
  
  # Initialize the process
  X_weak_2.0 <- numeric(N + 1) # N + 1 to include the initial point
  X_weak_2.0[1] <- X_0
  
  dW <- rnorm(N, mean = 0, sd = sqrt(step_length))
  time <- step_length * 0:N
  lambda_t <- numeric(N + 1)
  lambda_t <- lambda_0 * (1 - (time >= t_0) * (time - t_0) / tau)
  alpha_t <- 2 * sqrt(A * abs(lambda_t))
  mu_t <- m + sqrt(abs(lambda_t)/A)
  
  for(i in (2:(N+1))){
    X_weak_2.0[i] <- X_weak_2.0[i - 1] - (A*(X_weak_2.0[i - 1] - m)^2 + lambda_t[i - 1]) * step_length +
      sigma * sqrt(X_weak_2.0[i - 1]) * dW[i - 1] +
      1/4 * sigma^2 * (dW[i - 1]^2 - step_length) -
      2 * (A * (X_weak_2.0[i - 1] - m)) * sigma * sqrt(X_weak_2.0[i - 1]) * (1/2 * dW[i - 1] * step_length) +
      ((A * (X_weak_2.0[i - 1] - m)^2 + lambda_t[i - 1]) * (A * (X_weak_2.0[i - 1] - m)) +
         1/2 * A * sigma^2 * X_weak_2.0[i - 1]) * (step_length)^2 -
      1/2 * ((A * (X_weak_2.0[i - 1] - m)^2 + lambda_t[i - 1]) * 1/sqrt(X_weak_2.0[i - 1]) +
               sigma^2/(4 * sqrt(X_weak_2.0[i - 1]))) * (1/2 * dW[i - 1] * step_length)
    
  }
  tibble::tibble(t = time,
         X_weak_2.0,
         lambda_t, alpha_t, mu_t)
}

## Square-root noise test simulations
# true_param <- c(0.015, 300, -2, 0.1)
# actual_dt <- 0.001
# tau <- 100
# t_0 <- 50
# sim_res_sqrt <- simulate_squareroot_noise_tipping_model(actual_dt, true_param, tau, t_0)
# sim_res_sqrt |> ggplot2::ggplot(ggplot2::aes(x = t, y = X_weak_2.0)) +
#   ggplot2::geom_step() + ggplot2::geom_hline(yintercept = true_param[2], linetype = "dashed") +
#   ggplot2::geom_vline(xintercept = t_0)
# 
# sim_res_sqrt |>  ggplot2::ggplot(ggplot2::aes(x = t, y = lambda_t)) +
#   ggplot2::geom_step() + ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
#   ggplot2::geom_vline(xintercept = tau + t_0, linetype = "dashed")  +
#   ggplot2::geom_vline(xintercept = t_0, linetype = "dashed")

simulate_linear_noise_tipping_model<- function(step_length, param,
                                               tau = 100, t_0 = 10,
                                               X_0 = NA, beyond_tipping = 0){
  A <- param[1]
  m <- param[2]
  lambda_0 <- param[3]
  sigma <- param[4]
  
  if(is.na(X_0)){
    X_0 <- m + sqrt(abs(lambda_0) / A)
  }
  
  total_time <- tau + t_0 + beyond_tipping #Tipping point
  
  N <- as.integer((total_time) / step_length) 
  
  # Initialize the process
  X_weak_2.0 <- numeric(N + 1) # N + 1 to include the initial point
  X_weak_2.0[1] <- X_0
  
  dW <- rnorm(N, mean = 0, sd = sqrt(step_length))
  time <- step_length * 0:N
  lambda_t <- numeric(N + 1)
  lambda_t <- lambda_0 * (1 - (time >= t_0) * (time - t_0) / tau)
  alpha_t <- 2 * sqrt(A * abs(lambda_t))
  mu_t <- m + sqrt(abs(lambda_t)/A)
  
  for(i in (2:(N+1))){
    # Weak-order 2.0 Itô taylor based sampler
    X_weak_2.0[i] <- X_weak_2.0[i - 1] -(A*(X_weak_2.0[i - 1] - m)^2 + lambda_t[i - 1]) * step_length +
      sigma * X_weak_2.0[i - 1] * dW[i - 1] +
      1/2 * sigma^2 * X_weak_2.0[i - 1] * (dW[i - 1]^2 - step_length) -
      2 * A * (X_weak_2.0[i - 1] - m) * sigma * X_weak_2.0[i - 1] * (1/2 * dW[i - 1] * step_length) +
      1/2 * (2 * A * (X_weak_2.0[i - 1] - m) * (A * (X_weak_2.0[i - 1] - m)^2 + lambda_t[i - 1]) - A * sigma^2 * X_weak_2.0[i - 1]) * step_length^2 +
      (-(A * (X_weak_2.0[i - 1] - m)^2 + lambda_t[i - 1]) * sigma) * (dW[i - 1] * step_length - (1/2 * dW[i - 1] * step_length))
    
    
  }
  tibble::tibble(t = time,
         X_weak_2.0,
         lambda_t, alpha_t, mu_t)
}

## Linear noise test simulations
# true_param <- c(0.005, 50, -0.05, 0.005)
# actual_dt <- 0.001
# tau <- 100
# t_0 <- 50
# sim_res_linear <- simulate_linear_noise_tipping_model(actual_dt, true_param, tau, t_0)
# sim_res_linear |> ggplot2::ggplot(ggplot2::aes(x = t, y = X_weak_2.0)) +
#   ggplot2::geom_step() + ggplot2::geom_hline(yintercept = true_param[2], linetype = "dashed") +
#   ggplot2::geom_vline(xintercept = t_0)
# 
# sim_res_linear |>  ggplot2::ggplot(ggplot2::aes(x = t, y = lambda_t)) +
#   ggplot2::geom_step() + ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
#   ggplot2::geom_vline(xintercept = tau + t_0, linetype = "dashed")  +
#   ggplot2::geom_vline(xintercept = t_0, linetype = "dashed")
