# Title: Model diagnostics of Stochastic Differential Equations with tipping points.
# Author: Anders Gantzhorn Kristensen (University of Copenhagen, andersgantzhorn@gmail.com)
# Date: 2024-02-28 (Last Updated: 2024-02-28)
#-----------------------------------------------------------------------------------------------------------------------------#
# Project: Tipping Point Estimation in Ecological Systems using Stochastic Differential Equations
# Description: This script implements diagnostics for the models considered in "model_fitting.R".
#-----------------------------------------------------------------------------------------------------------------------------#
# License: MIT License (for more information, see LICENSE file in the repository).
# Dependencies: None
# Usage: See example usage in the bottom of the script
# Acknowledgements: A few of the methods was made by me and Lucas Støjko for the mandatory project during the course:
# "NMAK23005U Inference for Stochastic Differential Equations" at the University of Copenhagen in 2023-2024.
# Additionally, some of the code is inspired by the code from Susanne - and Peter Ditlevsens paper: 
# "Warning of a forthcoming collapse of the Atlantic meridional overturning circulation"
# Keywords: Stochastic Differential Equations, Likelihood methods, Estimation, diagnostics.

#-----------------------------------------------------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------------------------------------------------#
# Implementation of methods that are based on the square-root noise term process
CIR_quadratic_martingale_resid <- function(par, data, delta) {
  
  data <- - 4 * par[1] / (par[3]^2 * expm1(-par[1] * delta)) * data
  
  Xlow <- data[1:(length(data) - 1)]
  Xupp <- data[2:length(data)]
  
  beta <- par[1]
  mu <- par[2]
  sigma <- par[3]

  qnorm(pchisq(Xupp, df = 4 * mu * beta / sigma^2, ncp = Xlow * exp(-beta * delta)))
}

CIR_strang_splitting_resid <- function(par, data, delta) {
  x0 <- data[1:(length(data) - 1)]
  x1 <- data[2:length(data)]
  
  beta <- par[1]
  mu <- par[2]
  sigma <- par[3]
  
  A <- - (sigma^2 + 2 * beta) / 2
  b <- - log((2 * beta + sigma^2) / (2 * beta * mu))
  
  diff_f <- function(t, y){-beta + beta * mu * exp(-y) - sigma^2 / 2 + (sigma^2 + 2 * beta) / 2 * y +
      (sigma^2 + 2 * beta) / 2 * log((sigma^2 + 2 * beta) / (2 * beta * mu))}
  
  f_h <- runge_kutta(x0, delta / 2, diff_f, n = 1)

  inv_f <- runge_kutta(x1, -delta / 2, diff_f, n = 1)
  
  mu_f <- exp(A * delta) * (f_h - b) + b
  sd_f <- sigma * sqrt(expm1(2 * A * delta) / (2 * A))

  qnorm(pnorm(inv_f, mean = mu_f, sd = sd_f))
}

CIR_dynamic_likelihood_resid <- function(par, data, delta, alpha0, mu0, sigma, pen = 0){
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
  fh_half_tmp_low <-  A * delta * (Xlow - mu_seq) / 2
  fh_half_tmp_upp <-  A * delta * (Xupp - mu_seq) / 2
  
  fh_half     <-  (mu_seq * fh_half_tmp_low + Xlow) / (fh_half_tmp_low+1)
  
  fh_half_inv <-  (mu_seq * fh_half_tmp_upp - Xupp) / (fh_half_tmp_upp-1)
  
  det_Dfh_half_inv <- 1 / (fh_half_tmp_upp-1)^2
  phi_dot <- fh_half / alpha_seq * (exp(-alpha_seq * delta) - exp(-2 * alpha_seq * delta)) +
    mu_seq / (2 * alpha_seq) * (-expm1(-alpha_seq * delta))^2
  
  if(any(phi_dot < 0) | any(is.na(phi_dot))){return(100000)}
  
  sd.part <- sigma * sqrt(phi_dot)
  
  mu.part <- exp(-alpha_seq * delta) * (fh_half - mu_seq) + mu_seq
  qnorm(pnorm(fh_half_inv, mean = mu.part, sd = sd.part))
}

#-----------------------------------------------------------------------------------------------------------------------------#
# Example usage
ggplot2::theme_set(ggthemes::theme_base())
source("source_code/tipping_simulations.R")
source("source_code/model_fitting.R")
# Square-root noise model
true_param <- c(1, 3, -3, 0.3)
actual_dt <- 0.01
tau <- 100
t_0 <- 50
sim_res_sqrt <- simulate_squareroot_noise_tipping_model(actual_dt, true_param, tau, t_0)
# 
# # Stationary part
# # Parameters for stationary part
mu0 = true_param[2] + sqrt(abs(true_param[3]) / true_param[1])
alpha0 <- 2 * sqrt(true_param[1] * abs(true_param[3]))
stationary_part_true_param <- c(alpha0, mu0, true_param[4])
CIR_stationary_part_estimated_param_martingale <- CIR_quadratic_martingale(sim_res_sqrt$X_weak_2.0[sim_res_sqrt$t < t_0],
                                                                           actual_dt)
# tibble::tibble(obsSample = CIR_quadratic_martingale_resid(CIR_stationary_part_estimated_param_martingale,
#                                data = sim_res_sqrt$X_weak_2.0[sim_res_sqrt$t < t_0],
#                                delta = actual_dt)) |> ggplot2::ggplot(ggplot2::aes(sample = obsSample)) +
#   ggplot2::geom_qq() + ggplot2::geom_qq_line()
# 
# 
#   
# CIR_stationary_part_estimated_param_strang <- optimize_stationary_likelihood(likelihood_fun = CIR_strang_splitting,
#                                data = 2*sqrt(sim_res_sqrt$X_weak_2.0[sim_res_sqrt$t < t_0]),
#                                init_par = stationary_part_true_param,
#                                delta = actual_dt, exp_sigma = FALSE)
# 
# tibble::tibble(obsSample = CIR_strang_splitting_resid(CIR_stationary_part_estimated_param_strang,
#                            2*sqrt(sim_res_sqrt$X_weak_2.0[sim_res_sqrt$t < t_0]),
#                            actual_dt)) |> ggplot2::ggplot(ggplot2::aes(sample = obsSample)) +
#   ggplot2::geom_qq() + ggplot2::geom_qq_line()

# dynamic_part_true_param <- c(tau, true_param[1])
# CIR_dynamic_part_estimated_param_strang <- optimize_dynamic_likelihood(likelihood_fun = CIR_dynamic_likelihood,
#                             data = sim_res_sqrt$X_weak_2.0[sim_res_sqrt$t > t_0],
#                             init_par = dynamic_part_true_param,
#                             delta = actual_dt,
#                             alpha0 = stationary_part_true_param[1],
#                             mu0 = stationary_part_true_param[2],
#                             sigma = stationary_part_true_param[3])
# 
# tibble::tibble(obsSample = CIR_dynamic_likelihood_resid(CIR_dynamic_part_estimated_param_strang,
#                              data = sim_res_sqrt$X_weak_2.0[sim_res_sqrt$t > t_0],
#                              delta = actual_dt,
#                              alpha0 = stationary_part_estimated_param[1],
#                              mu0 = stationary_part_estimated_param[2],
#                              sigma = stationary_part_estimated_param[3])) |> ggplot2::ggplot(ggplot2::aes(sample = obsSample)) +
#     ggplot2::geom_qq() + ggplot2::geom_qq_line()
  
