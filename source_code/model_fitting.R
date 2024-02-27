# Title: Model fitting of Stochastic Differential Equations with tipping points.
# Author: Anders Gantzhorn Kristensen (University of Copenhagen, andersgantzhorn@gmail.com)
# Date: 2024-02-22 (Last Updated: 2024-02-22)
#-----------------------------------------------------------------------------------------------------------------------------#
# Project: Tipping Point Estimation in Ecological Systems using Stochastic Differential Equations
# Description: This script implements optimizers, estimation methods and negative log-likelihood functions for estimation
# of tipping points in stochastic differential equations with additive-, square-root- and linear noise structure
#-----------------------------------------------------------------------------------------------------------------------------#
# License: MIT License (for more information, see LICENSE file in the repository).
# Dependencies: None
# Usage: See example usage in the bottom of the script
# Acknowledgements: A few of the methods was made by me and Lucas Støjko for the mandatory project during the course:
# "NMAK23005U Inference for Stochastic Differential Equations" at the University of Copenhagen in 2023-2024.
# Additionally, some of the code is inspired by the code from Susanne - and Peter Ditlevsens paper: 
# "Warning of a forthcoming collapse of the Atlantic meridional overturning circulation"
# Keywords: Stochastic Differential Equations, Likelihood methods, Estimation, Numerical Optimization

#-----------------------------------------------------------------------------------------------------------------------------#
# Implementation of methods that are based on the additive noise term process
OU_likelihood <- function(par, data, delta){
  alpha0 <- par[1]
  mu0 <- par[2]
  sigma <- par[3]
  N <- length(data)
  
  Xupp <- data[2:N]
  Xlow <- data[1:(N - 1)]
  
  gamma2 <- sigma^2 / (2 * alpha0)
  rho0   <- exp(-alpha0 * delta)
  m_part <- Xupp - Xlow * rho0 - mu0 * (1 - rho0)
  v_part <- gamma2 * (1 - rho0^2)
  # Negative log-likelihood
  N * (log(v_part)) + sum(m_part^2 / v_part) 
}

score_OU <- function(par, data, dt){
  x0 <- data[1:(length(data) - 1)]
  x1 <- data[2:length(data)]
  n <- length(x0)
  beta <- par[1]
  mu <-  par[2]
  sigma <- par[3]
  
  gamma_square <- sigma^2 / (2 * beta)
  rho <- exp(- beta * dt)
  
  eq1 <- (1 - rho) / (gamma_square * (1 - rho^2)) * sum(x1 - x0 * rho - mu * (1 - rho))
  
  eq2 <- n * rho / (1 - rho^2) + 
    sum((x1 - x0 * rho - mu * (1 - rho)) * (x0 - mu)) / (gamma_square * (1 - rho^2)) - 
    rho * sum((x1 - x0 * rho - mu * (1 - rho))^2) / (gamma_square * (1 - rho^2)^2)
  
  eq3 <- - n / (2 * gamma_square) + sum((x1 - x0 * rho - mu * (1 - rho))^2) / (2 * gamma_square^2 * (1 - rho^2))
  
  c(eq1, eq2, eq3)
}


OU_dynamic_likelihood <-  function(par, data, delta, alpha0, mu0, sigma, pen = 0){
  tau     <- par[1]
  A       <- exp(par[2])

  N       <- length(data)
  Xupp    <- data[2:N]
  Xlow    <- data[1:(N-1)]
  time    <- delta * (1:(N-1))
  
  m       <- mu0 - alpha0/(2 * A)
  lambda0 <- -alpha0^2 / (4 * A)
  lam_seq <- lambda0 * (1 - time / tau)
  alpha_seq  <- 2 * sqrt(A * abs(lam_seq))
  gamma2_seq <- sigma^2 / (2 * alpha_seq)
  rho_seq    <- exp(-alpha_seq * delta)
  mu_seq     <- m + sqrt(abs(lam_seq) / A)
  
  ## Calculating the Strang splitting scheme pseudo likelihood
  fh_half_tmp_low <- A * delta * (Xlow - mu_seq) / 2
  fh_half_tmp_upp <- A * delta * (Xupp - mu_seq) / 2
  
  fh_half     <- (mu_seq * fh_half_tmp_low + Xlow) / (fh_half_tmp_low + 1)
  
  fh_half_inv <- (mu_seq * fh_half_tmp_upp - Xupp) / (fh_half_tmp_upp - 1)
  
  mu_h        <- fh_half * rho_seq + mu_seq * (1 - rho_seq)
  m_part      <- fh_half_inv - mu_h
  v_part    <- gamma2_seq * (1 - rho_seq^2)
  det_Dfh_half_inv <- 1 / (fh_half_tmp_upp-1)^2
  # Negative log-likelihood
  sum(log(v_part)) + sum(m_part^2 / v_part) -
    2 * sum(log(det_Dfh_half_inv)) + pen * N * (1 / A - 1) * (A < exp(1))
}

#-----------------------------------------------------------------------------------------------------------------------------#
# Implementation of methods that are based on the square-root noise term process

CIR_quadratic_martingale <- function(data, delta) {
  N <- length(data)
  M <- N - 1
  Xlow <- data[1:M]
  Xupp <- data[2:N]
  S1 <- sum(Xupp)
  S2 <- sum(1 / Xlow)
  S3 <- sum(Xupp / Xlow)
  
  ebdel <- (M * S3 - S1 * S2) / (M^2 - S1 * S2)
  beta <- -log(ebdel) / dt
  
  mu <- (S3 - M * ebdel) / ((1 - ebdel) * S2)
  
  S4 <- sum((Xupp - mu + ebdel * (mu - Xlow))^2 / Xlow)
  S5 <- sum((Xupp * (ebdel - ebdel^2) / beta + mu * (1 - ebdel)^2 / (2 * beta)) / Xupp)
  
  sigma <- S4 / S5
  
  c(beta = unname(beta), mu = unname(mu), sigma = unname(sqrt(sigma)))
}

CIR_dynamic_likelihood <- function(par, data, delta, alpha0, mu0, sigma, pen = 0){
  tau     <-  par[1]
  A       <-  exp(par[2])

  N       <- length(data)
  Xupp    <- data[2:N]
  Xlow    <- data[1:(N-1)]
  time    <- delta * (1:(N-1))
  
  m       <-  mu0 - alpha0 / (2 * A)
  lambda0 <-  -alpha0^2 / (4 * A)
  lam_seq    <- lambda0 * (1 - time / tau)
  alpha_seq  <- 2 * sqrt(A * abs(lam_seq))
  mu_seq     <- m + sqrt(abs(lam_seq) / A)
  
  # ## Calculating the Strang splitting scheme pseudo likelihood
  fh_half_tmp_low = A * delta * (Xlow - mu_seq) / 2
  fh_half_tmp_upp = A * delta * (Xupp - mu_seq) / 2
  
  fh_half     = (mu_seq * fh_half_tmp_low + Xlow) / (fh_half_tmp_low+1)
  
  fh_half_inv = (mu_seq * fh_half_tmp_upp - Xupp) / (fh_half_tmp_upp-1)
  
  det_Dfh_half_inv = 1 / (fh_half_tmp_upp-1)^2
  
  sd.part <- sqrt(delta) * sigma * sqrt(fh_half)
  
  mu.part <- fh_half - alpha_seq * (fh_half - mu_seq) * delta
  
  negloglik <- -sum(stats::dnorm(fh_half_inv, mu.part, sd.part, log = TRUE)) - sum(log(abs(det_Dfh_half_inv)))
  negloglik
}

#-----------------------------------------------------------------------------------------------------------------------------#
# Implementation of methods that are based on the linear noise term process
mean_reverting_GMB_EM_likelihood <- function(par, data, delta){
  beta <- par[1]
  mu <- par[2]
  sigma <- exp(par[3])
  
  N <- length(data)
  
  Xupp <- data[2:N]
  Xlow <- data[1:(N - 1)]
  
  sd.part <- sqrt(delta) * sigma * Xlow
  
  mu.part <- Xlow - beta * (Xlow - mu) * delta
  -sum(stats::dnorm(Xupp, mu.part, sd.part, log = TRUE))
}

mean_reverting_GMB_Kessler_likelihood <- function(par, data, delta){
  beta <- par[1]
  mu <- par[2]
  sigma <- exp(par[3])
  
  N <- length(data)
  
  Xupp <- data[2:N]
  Xlow <- data[1:(N - 1)]
  
  mu_ks <- exp(-beta * delta) * (Xlow - mu) + mu
  
  
  var_ks <- ((Xlow-mu)^2  * (1 - exp(-sigma^2 * delta)) -
               2 * mu * sigma^2 / (beta - sigma^2) * (Xlow - mu) * (1 - exp(-(sigma^2 - beta) * delta)) - 
               mu^2 * sigma^2 / (2 * beta - sigma^2) * (1 - exp(-(sigma^2 - 2 * beta) * delta))) * exp(-(2 * beta - sigma^2) * delta)
  
  sd_ks <- sqrt(var_ks)
  
  -sum(dnorm(Xupp, mean = mu_ks, sd = sd_ks, log = TRUE)) 
}

mean_reverting_GMB_martingale <- function(par, data, delta){
  x0 <- data[1:(length(data) - 1)]
  x1 <- data[2:length(data)]
  
  beta <- par[1]
  mu <-  par[2]
  sigma <- par[3]
  
  F_i <- exp(-beta * delta) * (x0 - mu) + mu
  
  phi_i <- exp(-(2 * beta - sigma^2) * delta) * (
    (x0 - mu)^2 * (1 - exp(-sigma^2 * delta)) - 
      2 * mu * sigma^2 / (beta - sigma^2) * (x0 - mu) * (1 - exp(-(sigma^2 - beta) * delta)) -
      mu^2 * sigma^2 / (2 * beta - sigma^2) * (1 - exp(-(sigma^2 - 2*beta) * delta))
  )
  
  eq1 <- beta / sigma^2 * sum(1 / x0 * (x1 - F_i))
  eq2 <- - 1 / sigma^2 * sum(1 / x0 * (x1 - F_i))
  eq3 <- 1 / sigma^3 * sum(1 / x0^2 * ((x1 - F_i)^2 - phi_i))
  
  c(eq1, eq2, eq3)
}

mean_reverting_GMB_dynamic_likelihood <- function(par, data, delta, alpha0, mu0, sigma, pen = 0){
  tau     <-  par[1]
  A       <-  exp(par[2])
  
  N       <- length(data)
  Xupp    <- data[2:N]
  Xlow    <- data[1:(N-1)]
  time    <- delta * (1:(N-1))
  
  m       <-  mu0 - alpha0 / (2 * A)
  lambda0 <-  -alpha0^2 / (4 * A)
  lam_seq    <- lambda0 * (1 - time / tau)
  alpha_seq  <- 2 * sqrt(A * abs(lam_seq))
  mu_seq     <- m + sqrt(abs(lam_seq) / A)
  
  # ## Calculating the Strang splitting scheme pseudo likelihood
  fh_half_tmp_low = A * delta * (Xlow - mu_seq) / 2
  fh_half_tmp_upp = A * delta * (Xupp - mu_seq) / 2
  
  fh_half     = (mu_seq * fh_half_tmp_low + Xlow)/(fh_half_tmp_low+1)
  
  fh_half_inv = (mu_seq * fh_half_tmp_upp - Xupp)/(fh_half_tmp_upp-1)
  
  det_Dfh_half_inv = 1/(fh_half_tmp_upp-1)^2
  
  sd.part <- sqrt(delta) * sigma * fh_half
  
  mu.part <- fh_half - alpha_seq * (fh_half - mu_seq) * delta
  
  negloglik <- -sum(stats::dnorm(fh_half_inv, mu.part, sd.part, log = TRUE)) - sum(log(abs(det_Dfh_half_inv)))
  negloglik
}

#-----------------------------------------------------------------------------------------------------------------------------#
## General optimizers for stationary - and dynamic part of process
optimize_OU <- function(data, init_par, delta){
  res_optim <- optim(init_par, fn = OU_likelihood,
                     method = "BFGS",
                     data = data, delta = delta)
  res_optim$par
}

optimize_stationary_likelihood <- function(likelihood_fun, data, init_par, delta){
  res_optim <- optim(init_par, fn = likelihood_fun,
                     method = "BFGS",
                     data = data, delta = delta)
  
  res_optim$par[3] <- exp(res_optim$par[3])
  res_optim$par
}

optimize_dynamic_likelihood <- function(likelihood_fun, data,
                                        init_par, delta,
                                        alpha0, mu0, sigma, pen = 0){

  res_optim <- stats::optim(init_par, fn = likelihood_fun,
                            method = "Nelder-Mead",
                            data = data, delta = delta,
                            alpha0 = alpha0 , mu0 = mu0, sigma = sigma, pen = pen)
  res_optim$par[2] <- exp(res_optim$par[2])
  res_optim$par
}

#-----------------------------------------------------------------------------------------------------------------------------#
# Example usage
# source("tipping_simulations.R")
# 
# ## Additive noise model
# 
# true_param <- c(0.87, -1.51, -2.69, 0.3)
# actual_dt <- 0.001
# tau <- 100
# t_0 <- 50
# sim_res_add <- simulate_additive_noise_tipping_model(actual_dt, true_param, tau, t_0)
# 
# ## Stationary part
# ## Parameters for stationary part
# mu0 = true_param[2] + sqrt(abs(true_param[3]) / true_param[1])
# alpha0 <- 2 * sqrt(true_param[1] * abs(true_param[3]))
# stationary_part_true_param <- c(alpha0, mu0, true_param[4])
# 
# optimize_OU(data = sim_res_add$X_weak_2.0[sim_res_add$t < t_0],
#                                init_par =  stationary_part_true_param,
#                                delta = actual_dt) - stationary_part_true_param
# 
# nleqslv::nleqslv(x = stationary_part_true_param,
#                  fn = score_OU,
#                  data = sim_res_add$X_weak_2.0[sim_res_add$t < t_0],
#                  dt = actual_dt)$x - stationary_part_true_param
# 
# ## Dynamic part
# dynamic_part_true_param <- c(tau, true_param[1])
# optimize_dynamic_likelihood(likelihood_fun = OU_dynamic_likelihood,
#                             data = sim_res_add$X_weak_2.0[sim_res_add$t > t_0],
#                             init_par = dynamic_part_true_param,
#                             delta = actual_dt,
#                             alpha0 = stationary_part_true_param[1],
#                             mu0 = stationary_part_true_param[2],
#                             sigma = stationary_part_true_param[3]) - dynamic_part_true_param


## Square-root noise model
# true_param <- c(0.5, 1.5, -2.69, 0.3)
# actual_dt <- 0.001
# tau <- 100
# t_0 <- 50
# sim_res_sqrt <- simulate_squareroot_noise_tipping_model(actual_dt, true_param, tau, t_0)


## Linear noise model
# true_param <- c(0.1, 1.5, -3, 0.15)
# actual_dt <- 0.001
# tau <- 100
# t_0 <- 50
# sim_res_linear <- simulate_linear_noise_tipping_model(actual_dt, true_param, tau, t_0)

