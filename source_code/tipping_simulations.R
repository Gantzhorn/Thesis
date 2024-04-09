# Title: Simulation of Stochastic Differential Equations for Tipping Points Analysis
# Author: Anders Gantzhorn Kristensen (University of Copenhagen, andersgantzhorn@gmail.com)
# Date: 2024-01-31 (Last Updated: 2024-04-09)
#-----------------------------------------------------------------------------------------------------------------------------#
# Project: Tipping Point Estimation in Ecological Systems using Stochastic Differential Equations
# Description: This script implements methods for simulating different types of Stochastic differential equation to model
# tipping points. By means of the scalar weak order 2.0 Itô-taylor method it aims to accurately simulate.
# sample paths from the rather complex bifurcation tipping point model with various noise structures. 
#-----------------------------------------------------------------------------------------------------------------------------#
# License: MIT License (for more information, see LICENSE file in the repository).
# Dependencies: ggplot2, ggthemes, tibble - All part of the tidyverse.
# Usage: See usage in the bottom of the script
# Acknowledgements: For algorithm see: Applied Stochastic Differential Equations Särkkä, Simo and Solin, Arno - Algorithm 8.5.
# Keywords: Stochastic Differential Equations, Tipping Points, Simulation, Itô-taylor weak order 2.0
#-----------------------------------------------------------------------------------------------------------------------------#
# Implementation for model with additive noise structure
simulate_additive_noise_tipping_model <- function(step_length, par,
                                                  tau = 100, t_0 = 10,
                                                  X_0 = NA, beyond_tipping = 0){
  A        <- par[1]
  m        <- par[2]
  lambda_0 <- par[3]
  sigma    <- par[4]
  nu       <- if(length(par) == 5) par[5] else 1
  
  if(is.na(X_0)){
    X_0 <- m + ifelse(A >=0, 1, -1) * sqrt(abs(lambda_0 / A))
  }
  
  #Tipping point with added time if specified
  total_time <- tau + t_0 + beyond_tipping 
  N          <- as.integer((total_time) / step_length) 
  
  # Initialize the process
  X_t    <- numeric(N + 1) # N + 1 to include the initial point
  X_t[1] <- X_0
  
  dW            <- rnorm(N, mean = 0, sd = sqrt(step_length))
  time          <- step_length * 0:N
  lambda_t      <- numeric(N + 1)
  lambda_t      <- lambda_0 * (1 - (time >= t_0) * (time - t_0) / tau)^nu
  alpha_t       <- 2 * sqrt(abs(A * lambda_t))
  mu_t          <- m + sqrt(abs(lambda_t/A))
  
  for(i in (2:(N+1))){
    X_t[i] <- X_t[i - 1] - (A*(X_t[i - 1] - m)^2 + lambda_t[i - 1]) * step_length +
      sigma * dW[i - 1] - 2 * A *(X_t[i - 1] - m) * sigma * (1 / 2 * dW[i - 1] * step_length) + 
      1 / 2 * (2 * A * (X_t[i - 1] - m) * (A * (X_t[i - 1] - m)^2 + lambda_t[i - 1]) - A * sigma^2) * step_length^2
    
  }
  
  tibble::tibble(t = time,
                 X_t,
                 lambda_t,
                 alpha_t,
                 mu_t)
}

#-----------------------------------------------------------------------------------------------------------------------------#
# Implementation for model with square-root noise structure
simulate_squareroot_noise_tipping_model<- function(step_length, par,
                                                   tau = 100, t_0 = 10,
                                                   X_0 = NA, beyond_tipping = 0){
  A        <- par[1]
  m        <- par[2]
  lambda_0 <- par[3]
  sigma    <- par[4]
  nu       <- if(length(par) == 5) par[5] else 1
  
  if(is.na(X_0)){
    X_0 <- m + ifelse(A >=0, 1, -1) * sqrt(abs(lambda_0 / A))
  }
  
  #Tipping point with added time if specified
  total_time <- tau + t_0 + beyond_tipping
  N          <- as.integer((total_time) / step_length) 
  
  # Initialize the process
  X_t    <- numeric(N + 1) # N + 1 to include the initial point
  X_t[1] <- X_0
  
  dW            <- rnorm(N, mean = 0, sd = sqrt(step_length))
  time          <- step_length * 0:N
  lambda_t      <- numeric(N + 1)
  lambda_t      <- lambda_0 * (1 - (time >= t_0) * (time - t_0) / tau)^nu
  alpha_t       <- 2 * sqrt(abs(A * lambda_t))
  mu_t          <- m + sqrt(abs(lambda_t/A))
  
  for(i in (2:(N+1))){
    X_t[i] <- X_t[i - 1] - (A*(X_t[i - 1] - m)^2 + lambda_t[i - 1]) * step_length +
      sigma * sqrt(X_t[i - 1]) * dW[i - 1] +
      1/4 * sigma^2 * (dW[i - 1]^2 - step_length) -
      2 * A * (X_t[i - 1] - m) * sigma * sqrt(X_t[i - 1]) * (1/2 * dW[i - 1] * step_length) +
      ((A * (X_t[i - 1] - m)^2 + lambda_t[i - 1]) * (A * (X_t[i - 1] - m)) -
         1/2 * A * sigma^2 * X_t[i - 1]) * (step_length)^2 -
      1 / (4 * sqrt(X_t[i - 1])) * ((A * (X_t[i - 1] - m)^2 + lambda_t[i - 1]) + sigma^3 / 4) *
      (dW[i - 1] * step_length)
  }
  tibble::tibble(t = time,
         X_t,
         lambda_t,
         alpha_t,
         mu_t)
}

#-----------------------------------------------------------------------------------------------------------------------------#
# Implementation for model with linear noise structure
simulate_linear_noise_tipping_model<- function(step_length, par,
                                               tau = 100, t_0 = 10,
                                               X_0 = NA, beyond_tipping = 0){
  A        <- par[1]
  m        <- par[2]
  lambda_0 <- par[3]
  sigma    <- par[4]
  nu       <- if(length(par) == 5) par[5] else 1
  
  if(is.na(X_0)){
    X_0 <- m + ifelse(A >= 0, 1, -1) * sqrt(abs(lambda_0 / A))
  }
  
  #Tipping point with added time if specified
  total_time <- tau + t_0 + beyond_tipping 
  N          <- as.integer((total_time) / step_length) 
  
  # Initialize the process
  X_t    <- numeric(N + 1) # N + 1 to include the initial point
  X_t[1] <- X_0
  
  dW            <- rnorm(N, mean = 0, sd = sqrt(step_length))
  time          <- step_length * 0:N
  lambda_t      <- numeric(N + 1)
  lambda_t      <- lambda_0 * (1 - (time >= t_0) * (time - t_0) / tau)^nu
  alpha_t       <- 2 * sqrt(abs(A * lambda_t))
  mu_t          <- m + sqrt(abs(lambda_t / A))
  
  for(i in (2:(N+1))){
    # Weak-order 2.0 Itô taylor based sampler
    X_t[i] <- X_t[i - 1] -(A*(X_t[i - 1] - m)^2 + lambda_t[i - 1]) * step_length +
      sigma * X_t[i - 1] * dW[i - 1] +
      1/2 * sigma^2 * X_t[i - 1] * (dW[i - 1]^2 - step_length) -
      2 * A * (X_t[i - 1] - m) * sigma * X_t[i - 1] * (1/2 * dW[i - 1] * step_length) +
      1/2 * (2 * A * (X_t[i - 1] - m) * (A * (X_t[i - 1] - m)^2 + lambda_t[i - 1]) -
               A * sigma^2 * X_t[i - 1]^2) * step_length^2 -
      ((A * (X_t[i - 1] - m)^2 + lambda_t[i - 1]) * sigma) *
      (1/2 * dW[i - 1] * step_length)
    
    
  }
  tibble::tibble(t = time,
         X_t,
         lambda_t,
         alpha_t,
         mu_t)
}

#-----------------------------------------------------------------------------------------------------------------------------#
# Implementation for model with t-distributed stationary distribution
simulate_t_distribution_tipping_model <- function(step_length, par,
                             tau = 100, t_0 = 10,
                             X_0 = NA, beyond_tipping = 0){
  A        <- par[1]
  m        <- par[2]
  lambda_0 <- par[3]
  sigma    <- par[4]
  nu       <- if(length(par) == 5) par[5] else 1
  
  if(is.na(X_0)){
    X_0 <- m + ifelse(A >=0, 1, -1) * sqrt(abs(lambda_0 / A))
  }
  
  #Tipping point with added time if specified
  total_time <- tau + t_0 + beyond_tipping 
  N          <- as.integer((total_time) / step_length) 
  
  # Initialize the process
  X_t    <- numeric(N + 1) # N + 1 to include the initial point
  X_t[1] <- X_0
  
  dW            <- rnorm(N, mean = 0, sd = sqrt(step_length))
  time          <- step_length * 0:N
  lambda_t      <- numeric(N + 1)
  lambda_t      <- lambda_0 * (1 - (time >= t_0) * (time - t_0) / tau)^nu
  alpha_t       <- 2 * sqrt(abs(A * lambda_t))
  mu_t          <- m + sqrt(abs(lambda_t/A))
  
  for(i in (2:(N+1))){
    X_t[i] <- X_t[i - 1] - (A * (X_t[i - 1] - m)^2 + lambda_t[i - 1]) * step_length + 
      sigma * sqrt(X_t[i - 1]^2 + 1) * dW[i - 1] + 
      1 / 2 * sigma^2 * X_t[i - 1] * (dW[i - 1]^2 - step_length) - 
      A * (X_t[i - 1] - m) * sigma * sqrt(X_t[i - 1]^2 + 1) * dW[i - 1] * step_length + 
      (((A * (X_t[i - 1] - m)^2) + lambda_t[i - 1]) * (A * (X_t[i - 1] - m)) -
                 1 / 2 * A * sigma^2 * (X_t[i - 1]^2 + 1)) * step_length^2 + 
      (- (A * (X_t[i - 1] - m)^2 + lambda_t[i - 1]) * sigma * X_t[i - 1] + 
         sigma^3 / 2) * (dW[i - 1] * step_length) / (2 * sqrt(X_t[i - 1]^2 + 1))
  }
  
  tibble::tibble(t = time,
                 X_t,
                 lambda_t,
                 alpha_t,
                 mu_t)
}

#-----------------------------------------------------------------------------------------------------------------------------#
# Implementation for model with F-distributed stationary distribution

simulate_F_distribution_tipping_model <- function(step_length, par, tau = 100, t_0 = 10, X_0 = NA, beyond_tipping = 0){
  A        <- par[1]
  m        <- par[2]
  lambda_0 <- par[3]
  sigma    <- par[4]
  nu       <- if(length(par) == 5) par[5] else 1
  
  if(is.na(X_0)){
    X_0 <- m + ifelse(A >=0, 1, -1) * sqrt(abs(lambda_0 / A))
  }
  
  #Tipping point with added time if specified
  total_time <- tau + t_0 + beyond_tipping 
  N          <- as.integer((total_time) / step_length) 
  
  # Initialize the process
  X_t    <- numeric(N + 1) # N + 1 to include the initial point
  X_t[1] <- X_0
  
  dW            <- rnorm(N, mean = 0, sd = sqrt(step_length))
  time          <- step_length * 0:N
  lambda_t      <- numeric(N + 1)
  lambda_t      <- lambda_0 * (1 - (time >= t_0) * (time - t_0) / tau)^nu
  alpha_t       <- 2 * sqrt(abs(A * lambda_t))
  mu_t          <- m + sqrt(abs(lambda_t/A))
  
  for(i in (2:(N+1))){
    X_t[i] <- X_t[i - 1] - (A * (X_t[i - 1] - m)^2 + lambda_t[i - 1]) * step_length + 
      sigma * sqrt(X_t[i - 1]^2 + X_t[i - 1]) * dW[i - 1] + 
      1 / 4 * sigma^2 * (2 * X_t[i - 1] + 1) * (dW[i - 1]^2 - step_length) -
      A * (X_t[i - 1] - m) * sigma * sqrt(X_t[i - 1] * (X_t[i - 1] + 1)) * dW[i - 1] * step_length +
      ((A * (X_t[i - 1] - m)^2 + lambda_t[i - 1]) * A * (X_t[i - 1] - m) - 
         sigma^2 * A / 2 * (X_t[i - 1] * (X_t[i - 1] + 1))) * step_length^2 -
      1 / (4 * sqrt(X_t[i - 1] * (X_t[i - 1] + 1))) *
      ((A * (X_t[i - 1] - m)^2 + lambda_t[i - 1]) * sigma * (2 * X_t[i - 1] + 1) +
      sigma^3 / 4) * step_length * dW[i - 1]
  }
  
  tibble::tibble(t = time,
                 X_t,
                 lambda_t,
                 alpha_t,
                 mu_t)
}

#-----------------------------------------------------------------------------------------------------------------------------#
# Implementation for model based on the jacobi diffusion

simulate_jacobi_diffusion_tipping_model <- function(step_length, par, tau = 100, t_0 = 10, X_0 = NA, beyond_tipping = 0){
  A        <- par[1]
  m        <- par[2]
  lambda_0 <- par[3]
  sigma    <- par[4]
  nu       <- if(length(par) == 5) par[5] else 1
  
  if(is.na(X_0)){
    X_0 <- m + ifelse(A >=0, 1, -1) * sqrt(abs(lambda_0 / A))
  }
  
  #Tipping point with added time if specified
  total_time <- tau + t_0 + beyond_tipping 
  N          <- as.integer((total_time) / step_length) 
  
  # Initialize the process
  X_t    <- numeric(N + 1) # N + 1 to include the initial point
  X_t[1] <- X_0
  
  dW            <- rnorm(N, mean = 0, sd = sqrt(step_length))
  time          <- step_length * 0:N
  lambda_t      <- numeric(N + 1)
  lambda_t      <- lambda_0 * (1 - (time >= t_0) * (time - t_0) / tau)^nu
  alpha_t       <- 2 * sqrt(abs(A * lambda_t))
  mu_t          <- m + sqrt(abs(lambda_t/A))
  
  for(i in (2:(N+1))){
    X_t[i] <- X_t[i - 1] - (A * (X_t[i - 1] - m)^2 + lambda_t[i - 1]) * step_length + 
      sigma * sqrt(X_t[i - 1] * (1 - X_t[i - 1])) * dW[i - 1] +
      sigma^2 / 4 * (1 - 2 * X_t[i - 1]) * (dW[i - 1] - step_length) -
      A * (X_t[i - 1] - m) * sigma * sqrt(X_t[i - 1] * (1 - X_t[i - 1])) * dW[i - 1] * step_length +
      ((A * (X_t[i - 1] - m)^2 + lambda_t[i - 1]) * (A * (X_t[i - 1] - m)) -
         sigma^2 * A / 2 * (X_t[i - 1] * (1 - X_t[i - 1]))) * step_length^2 - 
      ((A * (X_t[i - 1] - m) + lambda_t[i - 1]) * sigma * (1 - 2 * X_t[i - 1]) + sigma^3 / 2) *
      dW[i - 1] * step_length / (4 * sqrt(X_t[i - 1] * (1 - X_t[i - 1])))
  
    
  }
  tibble::tibble(t = time,
                 X_t,
                 lambda_t,
                 alpha_t,
                 mu_t)
}

#-----------------------------------------------------------------------------------------------------------------------------#
## Example usage of each of the methods
# ggplot2::theme_set(ggthemes::theme_base())

# # Additive noise example
# true_param <- c(0.87, -1.51, -2.69, 0.3, 1)
# actual_dt <- 0.001
# tau <- 100
# t_0 <- 50
# sim_res_add <- simulate_additive_noise_tipping_model(actual_dt, true_param, tau, t_0)
# sim_res_add |> ggplot2::ggplot(ggplot2::aes(x = t, y = X_t)) +
#   ggplot2::geom_step() + ggplot2::geom_hline(yintercept = true_param[2], linetype = "dashed") +
#   ggplot2::geom_vline(xintercept = t_0)
# 
# sim_res_add |> ggplot2::ggplot(ggplot2::aes(x = t, y = lambda_t)) +
#   ggplot2::geom_step() + ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
#   ggplot2::geom_vline(xintercept = c(t_0, tau + t_0))


## Square-root noise example
# true_param <- c(-1, 1.5, 1, 0.1, 0.85)
# 
# actual_dt <- 0.005
# tau <- 100
# t_0 <- 50
# sim_res_sqrt <- simulate_squareroot_noise_tipping_model(actual_dt, true_param, tau, t_0)
# sim_res_sqrt |> ggplot2::ggplot(ggplot2::aes(x = t, y = X_t)) +
#   ggplot2::geom_step() + ggplot2::geom_hline(yintercept = true_param[2], linetype = "dashed") +
#   ggplot2::geom_vline(xintercept = t_0)
# 
# sim_res_sqrt |> ggplot2::ggplot(ggplot2::aes(x = t, y = lambda_t)) +
#   ggplot2::geom_step() + ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
#   ggplot2::geom_vline(xintercept = c(t_0, tau + t_0))


## Linear noise example
# true_param <- c(0.15, 4.5, -1, 0.05)
# 
# actual_dt <- 0.005
# tau <- 100
# t_0 <- 50
# sim_res_linear <- simulate_linear_noise_tipping_model(actual_dt, true_param, tau, t_0)
# 
# sim_res_linear |> ggplot2::ggplot(ggplot2::aes(x = t, y = X_t)) +
#   ggplot2::geom_step() + ggplot2::geom_hline(yintercept = true_param[2], linetype = "dashed") +
#   ggplot2::geom_vline(xintercept = t_0)
# 
# sim_res_linear |> ggplot2::ggplot(ggplot2::aes(x = t, y = lambda_t)) +
#   ggplot2::geom_step() + ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
#   ggplot2::geom_vline(xintercept = c(t_0, tau + t_0))

## t-distributed example
# true_param <- c(0.5, -2, -3, 0.125)
# actual_dt <- 0.005
# tau <- 100
# t_0 <- 50
# sim_res_t_distribution <- simulate_t_distribution_tipping_model(actual_dt, true_param, tau, t_0)
# sim_res_t_distribution |> ggplot2::ggplot(ggplot2::aes(x = t, y = X_t)) +
#   ggplot2::geom_step() + ggplot2::geom_hline(yintercept = true_param[2], linetype = "dashed") +
#   ggplot2::geom_vline(xintercept = t_0)
# 
# sim_res_t_distribution |> ggplot2::ggplot(ggplot2::aes(x = t, y = lambda_t)) +
#   ggplot2::geom_step() + ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
#   ggplot2::geom_vline(xintercept = c(t_0, tau + t_0))


## F-distributed example
# true_param <- c(0.5, 3, -3, 0.1)
# actual_dt <- 0.005
# tau <- 100
# t_0 <- 50
# sim_res_F_distribution <- simulate_F_distribution_tipping_model(actual_dt, true_param, tau, t_0)
# sim_res_F_distribution |> ggplot2::ggplot(ggplot2::aes(x = t, y = X_t)) +
#   ggplot2::geom_step() + ggplot2::geom_hline(yintercept = true_param[2], linetype = "dashed") +
#   ggplot2::geom_vline(xintercept = t_0)
# 
# sim_res_F_distribution |> ggplot2::ggplot(ggplot2::aes(x = t, y = lambda_t)) +
#   ggplot2::geom_step() + ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
#   ggplot2::geom_vline(xintercept = c(t_0, tau + t_0))


## Jacobi diffusion 
# true_param <- c(3, 0.2, -1 , 0.1)
# actual_dt <- 0.01
# tau <- 100
# t_0 <- 50
# 
# sim_res_jacobi <- simulate_jacobi_diffusion_tipping_model(actual_dt, true_param, tau, t_0)
# 
# sim_res_jacobi |> ggplot2::ggplot(ggplot2::aes(x = t, y = X_t)) +
#   ggplot2::geom_step() + ggplot2::geom_hline(yintercept = true_param[2], linetype = "dashed") +
#   ggplot2::geom_vline(xintercept = t_0)
# 
# sim_res_jacobi |> ggplot2::ggplot(ggplot2::aes(x = t, y = lambda_t)) +
#   ggplot2::geom_step() + ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
#   ggplot2::geom_vline(xintercept = c(t_0, tau + t_0))
