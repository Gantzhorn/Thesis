# Title: Model diagnostics of Stochastic Differential Equations with tipping points.
# Author: Anders Gantzhorn Kristensen (University of Copenhagen, andersgantzhorn@gmail.com)
# Date: 2024-02-28 (Last Updated: 2024-05-28)
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
# Implementation of methods that are based on the additive noise process

OU_likelihood_resid <- function(par, data, delta){
  alpha0 <- par[1]
  mu0    <- par[2]
  sigma  <- par[3]
  
  N <- length(data)

  Xupp <- data[2:N]
  Xlow <- data[1:(N - 1)]
  
  gamma2 <- sigma^2 / (2 * alpha0)
  rho0   <- exp(-alpha0 * delta)
  
  v_part <- gamma2 * (1 - rho0^2)
  m_part <- Xlow * rho0 + mu0 * (1 - rho0)

  qnorm(p = pnorm(Xupp, mean = m_part, sd = sqrt(v_part)))
}

OU_dynamic_likelihood_resid <-  function(par, data, delta,
                                         alpha0, mu0, sigma){
  tau     <- par[1]
  A       <- par[2]
  nu      <- if(length(par) == 3) par[3] else 1
  
  
  N       <- length(data)
  Xupp    <- data[2:N]
  Xlow    <- data[1:(N-1)]
  time    <- delta * (1:(N-1))
  
  m          <- mu0 - alpha0 / (2 * A)
  lambda0    <- -alpha0^2 / (4 * A)
  lam_seq    <- lambda0 * (1 - time / tau)^nu
  alpha_seq  <- 2 * sqrt(abs(A * lam_seq))
  gamma2_seq <- sigma^2 / (2 * alpha_seq)
  rho_seq    <- exp(-alpha_seq * delta)
  mu_seq     <- m + ifelse(A >= 0, 1, -1) * sqrt(abs(lam_seq / A))
  
  ## Calculating the Strang splitting scheme pseudo likelihood
  fh_half_tmp_low <- A * delta * (Xlow - mu_seq) / 2
  fh_half_tmp_upp <- A * delta * (Xupp - mu_seq) / 2
  
  fh_half     <- (mu_seq * fh_half_tmp_low + Xlow) / (fh_half_tmp_low + 1)
  
  fh_half_inv <- (mu_seq * fh_half_tmp_upp - Xupp) / (fh_half_tmp_upp - 1)
  
  mu_h        <- fh_half * rho_seq + mu_seq * (1 - rho_seq)
  
  v_part    <- gamma2_seq * (1 - rho_seq^2)
  det_Dfh_half_inv <- 1 / (fh_half_tmp_upp-1)^2
  
  qnorm(pnorm(fh_half_inv, mean = mu_h, sd = sqrt(v_part)))
}

#-----------------------------------------------------------------------------------------------------------------------------#
# Implementation of methods that are based on the square-root noise term process
CIR_quadratic_martingale_resid <- function(par, data, delta) {
  beta  <- par[1]
  mu    <- par[2]
  sigma <- par[3]
  
  data <- - 4 * beta / (sigma^2 * expm1(-beta * delta)) * data
  
  Xlow <- data[1:(length(data) - 1)]
  Xupp <- data[2:length(data)]

  qnorm(pchisq(Xupp, df = 4 * mu * beta / sigma^2,
               ncp = Xlow * exp(-beta * delta)))
}

CIR_strang_splitting_resid <- function(par, data, delta) {
  Xlow <- data[1:(length(data) - 1)]
  Xupp <- data[2:length(data)]
  
  beta  <- par[1]
  mu    <- par[2]
  sigma <- par[3]
  
  A <- - beta
  b <- 2 * sqrt((beta * mu - sigma^2 / 2) / beta)
  
  diff_f <- function(t, y){beta / 2 * y + (2 * beta * mu - sigma^2) / y -
      2 * sqrt(beta * (beta * mu - sigma^2 / 2))}
  
  f_h   <- runge_kutta(Xlow, delta / 2, diff_f, n = 1)
  inv_f <- runge_kutta(Xupp, -delta / 2, diff_f, n = 1)
  
  mu_f <- exp(A * delta) * (f_h - b) + b
  sd_f <- sigma * sqrt(expm1(2 * A * delta) / (2 * A))

  qnorm(pnorm(inv_f, mean = mu_f, sd = sd_f))
}

CIR_dynamic_likelihood_resid <- function(par, data, delta,
                                         alpha0, mu0, sigma){
  tau     <-  par[1]
  A       <-  par[2]
  nu      <- if(length(par) == 3) par[3] else 1
  
  N       <- length(data)
  Xupp    <- data[2:N]
  Xlow    <- data[1:(N-1)]
  time    <- delta * (1:(N-1))
  
  m          <- mu0 - alpha0 / (2 * A)
  lambda0    <- -alpha0^2 / (4 * A)
  lam_seq    <- lambda0 * (1 - time / tau)^nu
  alpha_seq  <- 2 * sqrt(abs(A * lam_seq))
  mu_seq     <- m + ifelse(A >= 0, 1, -1) * sqrt(abs(lam_seq / A))
  
  # ## Calculating the Strang splitting scheme pseudo likelihood
  fh_half_tmp_low <-  A * delta * (Xlow - mu_seq) / 2
  fh_half_tmp_upp <-  A * delta * (Xupp - mu_seq) / 2
  
  fh_half     <-  (mu_seq * fh_half_tmp_low + Xlow) / (fh_half_tmp_low+1)
  
  fh_half_inv <-  (mu_seq * fh_half_tmp_upp - Xupp) / (fh_half_tmp_upp-1)
  
  det_Dfh_half_inv <- 1 / (fh_half_tmp_upp-1)^2
  phi_dot <- fh_half / alpha_seq * (exp(-alpha_seq * delta) - exp(-2 * alpha_seq * delta)) +
    mu_seq / (2 * alpha_seq) * (-expm1(-alpha_seq * delta))^2
  
  sd.part <- sigma * sqrt(phi_dot)
  
  mu.part <- exp(-alpha_seq * delta) * (fh_half - mu_seq) + mu_seq
  qnorm(pnorm(fh_half_inv, mean = mu.part, sd = sd.part))
}

CIR_transform_dynamic_likelihood_resid <- function(par, data, delta,
                                                   alpha0, mu0, sigma){
  tau     <-  par[1]
  A       <-  par[2]
  nu      <- if(length(par) == 3) par[3] else 1
  
  N       <- length(data)
  Xupp    <- data[2:N]
  Xlow    <- data[1:(N-1)]
  time    <- delta * (1:(N-1))
  
  m          <- mu0 - alpha0 / (2 * A)
  lambda0    <- -alpha0^2 / (4 * A)
  lam_seq    <- lambda0 * (1 - time / tau)^nu
  
  inner_sqrt_argument <- - A * (lam_seq + sigma^2 / 4)
  
  A_linear_part <- - sqrt(4 * inner_sqrt_argument)
  
  b <- 2 * sqrt(m + sqrt(inner_sqrt_argument) / A)
  
  diff_f <- function(t, y){-2 / y * (A * (y^2 / 4 - m)^2 + lam_seq + sigma^2 / 4) - 
      A_linear_part * (y - b)}
  
  # Solution to ODE
  f_h <- runge_kutta(Xlow, delta / 2, diff_f, n = 1)
  
  mu_f <- exp(A_linear_part * delta) * (f_h - b) + b
  sd_f <- sigma * sqrt(expm1(2 * A_linear_part * delta) / (2 * A_linear_part))
  
  # Inverse of non-linear ODE
  inv_f <- runge_kutta(Xupp, -delta / 2, diff_f, n = 1)
  
  # Strang likelihood
  qnorm(pnorm(inv_f, mean = mu_f, sd = sd_f))
}


#-----------------------------------------------------------------------------------------------------------------------------#
# Methods that are based on the linear noise process
# For martingale estimation equation 
mean_reverting_GBM_Kessler_likelihood_resid <- function(par, data, delta){
  beta <- par[1]
  mu <- par[2]
  sigma <- par[3]
  
  N <- length(data)
  
  Xupp <- data[2:N]
  Xlow <- data[1:(N - 1)]
  
  mu_ks <- exp(-beta * delta) * (Xlow - mu) + mu
  
  
  var_ks <- ((Xlow-mu)^2  * (1 - exp(-sigma^2 * delta)) -
               2 * mu * sigma^2 / (beta - sigma^2) * (Xlow - mu) *
               (1 - exp(-(sigma^2 - beta) * delta)) - 
               mu^2 * sigma^2 / (2 * beta - sigma^2) *
               (1 - exp(-(sigma^2 - 2 * beta) * delta))) *
    exp(-(2 * beta - sigma^2) * delta)
  
  sd_ks <- sqrt(var_ks)
  
  qnorm(pnorm(Xupp, mean = mu_ks, sd = sd_ks))
}

mean_reverting_GBM_alt_strang_resid <- function(par, data, delta,
                                                exp_sigma = FALSE) {
  Xlow <- data[1:(length(data) - 1)]
  Xupp <- data[2:length(data)]
  
  beta  <- par[1]
  mu    <- par[2]
  if(exp_sigma){sigma <- exp(par[3])} else{sigma <- par[3]}
  
  # Solution to non-linear ODE
  f_h    <- (Xlow - mu) * exp(-beta * delta / 2) + mu
  
  mu_log <- log(f_h) - sigma^2 / 2 * delta
  
  sd_log <- sigma * sqrt(delta)
  # Inverse to non-linear ODE
  inv_f  <- (Xupp - mu) * exp(beta * delta / 2) + mu
  
  # Strang likelihood
  qnorm(plnorm(inv_f, meanlog = mu_log, sdlog = sd_log))
}

mean_reverting_GBM_strang_resid <- function(par, data, delta,
                                            exp_sigma = FALSE) {
  Xlow <- data[1:(length(data) - 1)]
  Xupp <- data[2:length(data)]
  
  beta  <- par[1]
  mu    <- par[2]
  sigma <- par[3]
  
  A <- - (sigma^2 + 2 * beta) / 2
  b <- - log((2 * beta + sigma^2) / (2 * beta * mu))
  
  diff_f <- function(t, y){-beta + beta * mu * exp(-y) -
      sigma^2 / 2 - A * (y - b)}
  
  
  # Solution to non-linear ODE
  f_h    <- runge_kutta(Xlow, delta / 2, diff_f, n = 1)
  mu_f   <- exp(A * delta) * (f_h - b) + b
  sd_f   <- sigma * sqrt(expm1(2 * A * delta) / (2 * A))
  
  # Inverse to non-linear ODE
  inv_f  <- runge_kutta(Xupp, -delta / 2, diff_f, n = 1)
  
  # Derivative of inverse using Richardson Extrapolation
  inv_f2 <- runge_kutta(Xupp + 0.01, -delta / 2, diff_f, n = 1)
  inv_f3 <- runge_kutta(Xupp - 0.01, -delta / 2, diff_f, n = 1)
  df     <- (inv_f2 - inv_f3) / (2 * 0.01) 
  
  # Strang likelihood
  qnorm(pnorm(inv_f, mean = mu_f, sd = sd_f))
}

mean_reverting_GBM_dynamic_likelihood_resid <- function(par, data, delta,
                                                        alpha0, mu0, sigma){
  tau     <-  par[1]
  A       <-  par[2]
  nu      <- if(length(par) == 3) par[3] else 1
  
  
  N       <- length(data)
  Xupp    <- data[2:N]
  Xlow    <- data[1:(N-1)]
  time    <- delta * (1:(N-1))
  
  m          <- mu0 - alpha0 / (2 * A)
  lambda0    <- -alpha0^2 / (4 * A)
  lam_seq    <- lambda0 * (1 - time / tau)^nu
  alpha_seq  <- 2 * sqrt(abs(A * lam_seq))
  mu_seq     <- m + ifelse(A >= 0, 1, -1) * sqrt(abs(lam_seq / A))
  
  # Calculating the Strang splitting scheme pseudo likelihood
  fh_half_tmp_low <-  A * delta * (Xlow - mu_seq) / 2
  fh_half_tmp_upp <-  A * delta * (Xupp - mu_seq) / 2
  
  fh_half     <-  (mu_seq * fh_half_tmp_low + Xlow)/(fh_half_tmp_low+1)
  
  fh_half_inv <-  (mu_seq * fh_half_tmp_upp - Xupp)/(fh_half_tmp_upp-1)
  
  det_Dfh_half_inv <-  1/(fh_half_tmp_upp-1)^2
  
  sd.part <-  exp(-(alpha_seq - sigma^2) * delta) * sqrt(
    (fh_half - mu_seq)^2 * (1 - exp(-sigma^2 * delta)) - 
      2 * mu_seq * sigma^2 / (alpha_seq - sigma^2) * 
          (fh_half - mu_seq) * (1 - exp(-(sigma^2 - alpha_seq) * delta)) -
      mu_seq^2 * sigma^2 / (2 * alpha_seq - sigma^2) * 
          (1 - exp(-(sigma^2 - 2*alpha_seq) * delta))
  )
  
  mu.part <- exp(-alpha_seq * delta) * (fh_half - mu_seq) + mu_seq
  
  qnorm(pnorm(fh_half_inv, mean = mu.part, sd = sd.part))
}

mean_reverting_GBM_transform_dynamic_likelihood_resid <- function(par, data, delta,
                                                                  alpha0, mu0, sigma){
  tau     <-  par[1]
  A       <-  par[2]
  nu      <- if(length(par) == 3) par[3] else 1
  
  N       <- length(data)
  Xupp    <- data[2:N]
  Xlow    <- data[1:(N-1)]
  time    <- delta * (1:(N-1))
  
  m          <- mu0 - alpha0 / (2 * A)
  lambda0    <- -alpha0^2 / (4 * A)
  lam_seq    <- lambda0 * (1 - time / tau)^nu

  sqrt_argument <- sigma^4 / 4 - A * (2 * m * sigma^2 + 4 * lam_seq)
  
  exp_b <- (2 * m * A - 1 / 2 * sigma^2 + sqrt(sqrt_argument)) / (2 * A)
  b <- log(exp_b)
  
  A_linear_part <- - 1 / exp_b * (A * exp_b^2 - lam_seq - A * m^2)
  
  diff_f <- function(t, y){-((A * (exp(y) - m)^2 + lam_seq) / exp(y) + 1 / 2 * sigma^2) - 
      A_linear_part * (y - b)}
  
  # Solution to ODE
  f_h <- runge_kutta(Xlow, delta / 2, diff_f, n = 1)
  
  mu_f <- exp(A_linear_part * delta) * (f_h - b) + b
  sd_f <- sigma * sqrt(expm1(2 * A_linear_part * delta) / (2 * A_linear_part))
  
  # Inverse of non-linear ODE
  inv_f <- runge_kutta(Xupp, -delta / 2, diff_f, n = 1)
  
  qnorm(pnorm(inv_f, mean = mu_f, sd = sd_f))
}

#-----------------------------------------------------------------------------------------------------------------------------#

# Methods that are based on the t-distributed stationary process
t_diffusion_strang_splitting_resid <- function(par, data, delta) {
  Xlow <- data[1:(length(data) - 1)]
  Xupp <- data[2:length(data)]
  
  beta  <- par[1]
  mu    <- par[2]
  sigma <- par[3]
  
  asinh_term <- asinh(2 * mu * beta / (2 * beta + sigma^2))
  
  diff_f <- function(t, y){(y - tanh(y) - asinh_term) * (beta + 1 / 2 * sigma^2) +
      mu * beta / cosh(y)}
  
  # Solution to ODE
  f_h <- runge_kutta(Xlow, delta / 2, diff_f, n = 1)
  
  A <- -(beta + 1 / 2 * sigma^2)
  b <- asinh_term
  
  mu_f <- exp(A * delta) * (f_h - b) + b
  sd_f <- sigma * sqrt(expm1(2 * A * delta) / (2 * A))
  
  # Inverse of non-linear ODE
  inv_f <- runge_kutta(Xupp, -delta / 2, diff_f, n = 1)
  
  qnorm(pnorm(inv_f, mean = mu_f, sd = sd_f))
}

t_transform_dynamic_likelihood_resid <- function(par, data, delta,
                                                 alpha0, mu0, sigma){
  tau     <-  par[1]
  A       <-  par[2]
  nu      <- if(length(par) == 3) par[3] else 1
  
  N       <- length(data)
  Xupp    <- data[2:N]
  Xlow    <- data[1:(N-1)]
  time    <- delta * (1:(N-1))
  
  m          <- mu0 - alpha0 / (2 * A)
  lambda0    <- -alpha0^2 / (4 * A)
  lam_seq    <- lambda0 * (1 - time / tau)^nu
  
  A_linear_part <- -sqrt(sigma^4 / 4 - 2 * A * (2 * lam_seq + m * sigma^2))
  
  b <- asinh((sqrt(sigma^4 - 8 * A * (2 * lam_seq + m * sigma^2)) + 4 * A * m - sigma^2) /
               (4 * A))
  
  diff_f <- function(t, y){-1 / cosh(y) * (A * (sinh(y) - m)^2 + lam_seq + sigma^2 / 2 * sinh(y)) -
      A_linear_part * (y - b)}
  
  # Solution to ODE
  f_h <- runge_kutta(Xlow, delta / 2, diff_f, n = 1)
  
  mu_f <- exp(A_linear_part * delta) * (f_h - b) + b
  sd_f <- sigma * sqrt(expm1(2 * A_linear_part * delta) / (2 * A_linear_part))
  
  # Inverse of non-linear ODE
  inv_f <- runge_kutta(Xupp, -delta / 2, diff_f, n = 1)

  qnorm(pnorm(inv_f, mean = mu_f, sd = sd_f))
}

#-----------------------------------------------------------------------------------------------------------------------------#
# Implementation of methods that are based on the process with stationary F-distribution

F_diffusion_strang_splitting_resid <- function(par, data, delta) {
  Xlow <- data[1:(length(data) - 1)]
  Xupp <- data[2:length(data)]
  
  beta  <- par[1]
  mu    <- par[2]
  sigma <- par[3]
  
  A <- - (beta + sigma^2 / 2)
  b <- acosh(2 * beta * (2 * mu + 1) / (2 * beta + sigma^2))
  
  diff_f <- function(t,y){-A * (y - b - cosh(y) / sinh(y)) + beta * (2 * mu + 1) / sinh(y)}
  
  # Solution to ODE
  f_h <- runge_kutta(Xlow, delta / 2, diff_f, n = 1)
  
  mu_f <- exp(A * delta) * (f_h - b) + b
  sd_f <- sigma * sqrt(expm1(2 * A * delta) / (2 * A))
  
  # Inverse of non-linear ODE
  inv_f <- runge_kutta(Xupp, -delta / 2, diff_f, n = 1)
  
  qnorm(pnorm(inv_f, mean = mu_f, sd = sd_f))
}

F_transform_dynamic_likelihood_resid <- function(par, data, delta,
                                                 alpha0, mu0, sigma){
  tau     <-  par[1]
  A       <-  par[2]
  nu      <- if(length(par) == 3) par[3] else 1
  
  N       <- length(data)
  Xupp    <- data[2:N]
  Xlow    <- data[1:(N-1)]
  time    <- delta * (1:(N-1))
  
  m          <- mu0 - alpha0 / (2 * A)
  lambda0    <- -alpha0^2 / (4 * A)
  lam_seq    <- lambda0 * (1 - time / tau)^nu
  
  sqrt_argument <- sigma^4 / 4 - A * (sigma^2 * (2 * m + 1) + 4 * lam_seq)
  
  b <- acosh(((A * (2 * m + 1) - sigma^2 / 2) + sqrt(sqrt_argument)) / A)
  A_linear_part <- -sqrt(sqrt_argument)
  
  diff_f <- function(t, y){-1 / sinh(y) * (
    A / 2 * cosh(y)^2 + (sigma^2 / 2 - A * (2 * m + 1)) * cosh(y) +
      2 * lam_seq + 2 * A * m^2 + A / 2 + 2 * A * m) -
      A_linear_part * (y - b)}
  
  # Solution to ODE
  f_h <- runge_kutta(Xlow, delta / 2, diff_f, n = 1)
  
  mu_f <- exp(A_linear_part * delta) * (f_h - b) + b
  sd_f <- sigma * sqrt(expm1(2 * A_linear_part * delta) / (2 * A_linear_part))
  
  # Inverse of non-linear ODE
  inv_f <- runge_kutta(Xupp, -delta / 2, diff_f, n = 1)
  
  qnorm(pnorm(inv_f, mean = mu_f, sd = sd_f))
}

##-----------------------------------------------------------------------------------------------------------------------------#
# Implementation of methods that are based on the jacobi diffusion

jacobi_diffusion_strang_splitting_resid <- function(par, data, delta) {
  Xlow <- data[1:(length(data) - 1)]
  Xupp <- data[2:length(data)]
  
  beta  <- par[1]
  mu    <- par[2]
  sigma <- par[3]
  
  A <- 1/2 * sigma^2 - beta
  
  b <- acos(beta * (2 * mu - 1) / (sigma^2 / 2 - beta))
  
  diff_f <- function(t, y){-1 / sin(y) * (beta * (1 - cos(y) - 2 * mu) +
                                            sigma^2 / 2 * cos(y)) - A * (y - b)}
  
  # Solution to ODE
  f_h <- runge_kutta(Xlow, delta / 2, diff_f, n = 1)
  
  mu_f <- exp(A * delta) * (f_h - b) + b
  sd_f <- sigma * sqrt(expm1(2 * A * delta) / (2 * A))
  
  # Inverse of non-linear ODE
  inv_f <- runge_kutta(Xupp, -delta / 2, diff_f, n = 1)
  
  qnorm(pnorm(inv_f, mean = mu_f, sd = sd_f))
}

jacobi_diffusion_transform_dynamic_likelihood_resid <- function(par, data, delta,
                                                                alpha0, mu0, sigma){
  tau     <-  par[1]
  A       <-  par[2]
  nu      <- if(length(par) == 3) par[3] else 1
  
  N       <- length(data)
  Xupp    <- data[2:N]
  Xlow    <- data[1:(N-1)]
  time    <- delta * (1:(N-1))
  
  m          <- mu0 - alpha0 / (2 * A)
  lambda0    <- -alpha0^2 / (4 * A)
  lam_seq    <- lambda0 * (1 - time / tau)^nu
  
  sqrt_argument <- sigma^4 / 4 + sigma^2 * A * (2 * m - 1) - 4 * A * lam_seq
  
  acos_argument <- (-A * (2 * m - 1) - sigma^2 / 2 - sqrt(sqrt_argument)) / A

  A_linear_part <- -sqrt(sqrt_argument)
  b <- acos(acos_argument)
  
  diff_f <- function(t, y){-1 / sin(y) * 
      (A / 2 * cos(y)^2 + (sigma^2 / 2 + 2 * A * m - A) * cos(y) +
         A * (1 / 2 + 2 * m^2 - 2 * m) + 2 * lam_seq) -
      A_linear_part * (y - b)}
  
  # Solution to ODE
  f_h <- runge_kutta(Xlow, delta / 2, diff_f, n = 1)
  
  mu_f <- exp(A_linear_part * delta) * (f_h - b) + b
  sd_f <- sigma * sqrt(expm1(2 * A_linear_part * delta) / (2 * A_linear_part))
  
  # Inverse of non-linear ODE
  inv_f <- runge_kutta(Xupp, -delta / 2, diff_f, n = 1)
  
  qnorm(pnorm(inv_f, mean = mu_f, sd = sd_f))
}


#-----------------------------------------------------------------------------------------------------------------------------#
# Example usage
# ggplot2::theme_set(ggthemes::theme_base())
# source("source_code/tipping_simulations.R")
# source("source_code/model_fitting.R")

# Additive noise model

# true_param <- c(0.87, -1.51, -2.69, 0.3)
# actual_dt <- 1/12
# tau <- 100
# t_0 <- 50
# sim_res_add <- simulate_additive_noise_tipping_model(actual_dt, true_param, tau, t_0)
# sample_n(sim_res_add, min(nrow(sim_res_add), 10000)) |> ggplot(aes(x = t, y = X_t)) + geom_step()
## Stationary part
## Parameters for stationary part
# mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) * sqrt(abs(true_param[3] / true_param[1]))
# alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
# stationary_part_true_param <- c(alpha0, mu0, true_param[4])
# 
# OU_par <- optimize_stationary_likelihood(likelihood_fun = OU_likelihood, data = sim_res_add$X_t[sim_res_add$t < t_0],
#                                init_par = stationary_part_true_param,
#                                delta = actual_dt, exp_sigma = FALSE)$par
# 
# tibble::tibble(obsSample = OU_likelihood_resid(par = OU_par,
#          data = sim_res_add$X_t[sim_res_add$t < t_0],
#          delta = actual_dt)) |> ggplot2::ggplot(ggplot2::aes(sample = obsSample)) +
#     ggplot2::geom_qq() + ggplot2::geom_qq_line()
# 
# # # Dynamic part
# dynamic_part_true_param <- c(tau, true_param[1])
# OU_dynamic_par <- optimize_dynamic_likelihood(likelihood_fun = OU_dynamic_likelihood,
#                             data = sim_res_add$X_t[sim_res_add$t > t_0],
#                             init_par = dynamic_part_true_param,
#                             delta = actual_dt,
#                             alpha0 = OU_par[1],
#                             mu0 = OU_par[2],
#                             sigma = OU_par[3])$par
# #
# tibble::tibble(obsSample = OU_dynamic_likelihood_resid(OU_dynamic_par,
#                             data = sim_res_add$X_t[sim_res_add$t > t_0],
#                             delta = actual_dt,
#                             alpha0 =  OU_par[1],
#                             mu0 = OU_par[2],
#                             sigma = OU_par[3])) |>
#   ggplot2::ggplot(ggplot2::aes(sample = obsSample)) +
#   ggplot2::geom_qq() + ggplot2::geom_qq_line()

#-----------------------------------------------------------------------------------------------------------------------------#

# # Square-root noise model
# true_param <- c(1, 3, -2, 0.1)
# actual_dt <- 0.01
# tau <- 100
# t_0 <- 50
# sim_res_sqrt <- simulate_squareroot_noise_tipping_model(actual_dt, true_param, tau, t_0)
# sample_n(sim_res_sqrt, min(nrow(sim_res_sqrt), 10000)) |> ggplot(aes(x = t, y = X_t)) + geom_step()
#
# # Stationary part
# # Parameters for stationary part
# mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) * sqrt(abs(true_param[3] / true_param[1]))
# alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
# stationary_part_true_param <- c(alpha0, mu0, true_param[4])
# 
# CIR_stationary_part_estimated_param_martingale <-
#   CIR_quadratic_martingale(sim_res_sqrt$X_t[sim_res_sqrt$t < t_0],actual_dt)
# 
# tibble::tibble(obsSample = CIR_quadratic_martingale_resid(CIR_stationary_part_estimated_param_martingale,
#                                data = sim_res_sqrt$X_t[sim_res_sqrt$t < t_0],
#                                delta = actual_dt)) |> ggplot2::ggplot(ggplot2::aes(sample = obsSample)) +
#   ggplot2::geom_qq() + ggplot2::geom_qq_line()
# 
# CIR_stationary_part_estimated_param_strang <- optimize_stationary_likelihood(likelihood_fun = CIR_strang_splitting,
#                                data = 2*sqrt(sim_res_sqrt$X_t[sim_res_sqrt$t < t_0]),
#                                init_par = stationary_part_true_param,
#                                delta = actual_dt, exp_sigma = FALSE)$par
# 
# tibble::tibble(obsSample = CIR_strang_splitting_resid(CIR_stationary_part_estimated_param_strang,
#                            2*sqrt(sim_res_sqrt$X_t[sim_res_sqrt$t < t_0]),
#                            actual_dt)) |> ggplot2::ggplot(ggplot2::aes(sample = obsSample)) +
#   ggplot2::geom_qq() + ggplot2::geom_qq_line()
#
# Dynamic part
# 
# dynamic_part_true_param <- c(tau, true_param[1])
# CIR_dynamic_part_estimated_param_strang <- optimize_dynamic_likelihood(likelihood_fun = CIR_dynamic_likelihood,
#                             data = sim_res_sqrt$X_t[sim_res_sqrt$t > t_0],
#                             init_par = dynamic_part_true_param,
#                             delta = actual_dt,
#                             alpha0 = CIR_stationary_part_estimated_param_strang[1],
#                             mu0 = CIR_stationary_part_estimated_param_strang[2],
#                             sigma = CIR_stationary_part_estimated_param_strang[3])$par
# 
# tibble::tibble(obsSample = CIR_dynamic_likelihood_resid(CIR_dynamic_part_estimated_param_strang,
#                 data = sim_res_sqrt$X_t[sim_res_sqrt$t > t_0],
#                 delta = actual_dt,
#                 alpha0 = CIR_stationary_part_estimated_param_strang[1],
#                 mu0 = CIR_stationary_part_estimated_param_strang[2],
#                 sigma = CIR_stationary_part_estimated_param_strang[3])) |>
#   ggplot2::ggplot(ggplot2::aes(sample = obsSample)) +
#   ggplot2::geom_qq() + ggplot2::geom_qq_line()

# CIR_dynamic_part_estimated_param_strang_alt <-
#   optimize_dynamic_likelihood(likelihood_fun = CIR_transform_dynamic_likelihood,
#                               data = 2 * sqrt(sim_res_sqrt$X_t[sim_res_sqrt$t > t_0]),
#                               init_par = dynamic_part_true_param,
#                               delta = actual_dt,
#                               alpha0 = CIR_stationary_part_estimated_param_strang[1],
#                               mu0 = CIR_stationary_part_estimated_param_strang[2],
#                               sigma = CIR_stationary_part_estimated_param_strang[3])$par
#
# tibble::tibble(obsSample = CIR_transform_dynamic_likelihood_resid(CIR_dynamic_part_estimated_param_strang_alt,
#                              data = 2 * sqrt(sim_res_sqrt$X_t[sim_res_sqrt$t > t_0]),
#                              delta = actual_dt,
#                              alpha0 = CIR_stationary_part_estimated_param_strang[1],
#                              mu0 = CIR_stationary_part_estimated_param_strang[2],
#                              sigma = CIR_stationary_part_estimated_param_strang[3])) |>
#     ggplot2::ggplot(ggplot2::aes(sample = obsSample)) +
#     ggplot2::geom_qq() + ggplot2::geom_qq_line()

#-----------------------------------------------------------------------------------------------------------------------------#

# # Linear noise model
# true_param <- c(0.5, 5, -2, 0.025)
# actual_dt <- 1 / 12
# tau <- 150
# t_0 <- 50
# sim_res_linear <- simulate_linear_noise_tipping_model(actual_dt, true_param, tau, t_0)
# sample_n(sim_res_linear, min(nrow(sim_res_linear), 10000)) |>
#   ggplot2::ggplot(ggplot2::aes(x = t, y = X_t)) + ggplot2::geom_step()
# # 
# # Stationary part
# # Parameters for stationary part
# mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) * sqrt(abs(true_param[3] / true_param[1]))
# alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
# stationary_part_true_param <- c(alpha0, mu0, true_param[4])

# GBM_stationary_part_estimated_param_martingale <- nleqslv::nleqslv(x = stationary_part_true_param,
#                                                                    fn = mean_reverting_GBM_martingale,
#                  data = sim_res_linear$X_t[sim_res_linear$t < t_0],
#                  delta = actual_dt)$x
#
# tibble::tibble(obsSample =
#                  mean_reverting_GBM_Kessler_likelihood_resid(GBM_stationary_part_estimated_param_martingale,
#                  data = sim_res_linear$X_t[sim_res_linear$t < t_0],
#                  delta = actual_dt)) |>
#                  ggplot2::ggplot(ggplot2::aes(sample = obsSample)) +
#                  ggplot2::geom_qq() + ggplot2::geom_qq_line()

# GBM_stationary_part_estimated_param_strang <-
#   optimize_stationary_likelihood(mean_reverting_GBM_strang, log(sim_res_linear$X_t[sim_res_linear$t<t_0]),
#                                init_par = stationary_part_true_param, delta = actual_dt,
#                                exp_sigma = FALSE)$par
# 
# tibble::tibble(obsSample = mean_reverting_GBM_strang_resid(GBM_stationary_part_estimated_param_martingale,
#                     data = log(sim_res_linear$X_t[sim_res_linear$t < t_0]),
#                    delta = actual_dt)) |>
#   ggplot2::ggplot(ggplot2::aes(sample = obsSample)) +
#   ggplot2::geom_qq() + ggplot2::geom_qq_line()



# Dynamic part
# dynamic_part_true_param <- c(tau, true_param[1])

# GBM_dynamic_part_estimated_param_strang <- optimize_dynamic_likelihood(
#                             likelihood_fun = mean_reverting_GBM_dynamic_likelihood,
#                             data = sim_res_linear$X_t[sim_res_linear$t > t_0],
#                             init_par = dynamic_part_true_param,
#                             delta = actual_dt,
#                             alpha0 = GBM_stationary_part_estimated_param_martingale[1],
#                             mu0 = GBM_stationary_part_estimated_param_martingale[2],
#                             sigma = GBM_stationary_part_estimated_param_martingale[3])$par


# tibble::tibble(obsSample =
#                  mean_reverting_GBM_dynamic_likelihood_resid(GBM_dynamic_part_estimated_param_strang,
#                  data= sim_res_linear$X_t[sim_res_linear$t > t_0],
#                  delta = actual_dt,
#                  alpha0 = GBM_stationary_part_estimated_param_martingale[1],
#                  mu0 = GBM_stationary_part_estimated_param_martingale[2],
#                  sigma = GBM_stationary_part_estimated_param_martingale[3])) |>
#                  ggplot2::ggplot(ggplot2::aes(sample = obsSample)) +
#                  ggplot2::geom_qq() + ggplot2::geom_qq_line()

# GBM_transform_dynamic_part_estimated_param_strang <-  optimize_dynamic_likelihood(
#                             likelihood_fun = mean_reverting_GBM_transform_dynamic_likelihood,
#                             data = log(sim_res_linear$X_t[sim_res_linear$t > t_0]),
#                             init_par = dynamic_part_true_param,
#                             delta = actual_dt,
#                             alpha0 = GBM_stationary_part_estimated_param_martingale[1],
#                             mu0 = GBM_stationary_part_estimated_param_martingale[2],
#                             sigma = GBM_stationary_part_estimated_param_martingale[3])$par

# tibble::tibble(obsSample =
#                  mean_reverting_GBM_transform_dynamic_likelihood_resid(
#                  par = GBM_transform_dynamic_part_estimated_param_strang,
#                  data = log(sim_res_linear$X_t[sim_res_linear$t > t_0]),
#                  delta = actual_dt,
#                  alpha0 = GBM_stationary_part_estimated_param_martingale[1],
#                  mu0 = GBM_stationary_part_estimated_param_martingale[2],
#                  sigma = GBM_stationary_part_estimated_param_martingale[3])) |>
#                  ggplot2::ggplot(ggplot2::aes(sample = obsSample)) +
#                  ggplot2::geom_qq() + ggplot2::geom_qq_line()

#-----------------------------------------------------------------------------------------------------------------------------#

## t-distributed stationary distribution model

# true_param <- c(0.8, -1.51, -2.67, 0.15)
# actual_dt <- 1 / 12
# tau <- 100
# t_0 <- 50
# sim_res_t_distribution <- simulate_t_distribution_tipping_model(actual_dt, true_param, tau, t_0)
# sim_res_t_distribution |> ggplot2::ggplot(ggplot2::aes(x = t, y = X_t)) +
#   ggplot2::geom_step() + ggplot2::geom_hline(yintercept = true_param[2], linetype = "dashed") +
#   ggplot2::geom_vline(xintercept = t_0)

# Stationary part
# Parameters for stationary part
# mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) * sqrt(abs(true_param[3] / true_param[1]))
# alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
# stationary_part_true_param <- c(alpha0, mu0, true_param[4])

# stationary_part_estim_param <- optimize_stationary_likelihood(
#   likelihood_fun = t_diffusion_strang_splitting,
#   data = asinh(sim_res_t_distribution$X_t[sim_res_t_distribution$t < t_0]),
#   init_par = stationary_part_true_param,
#   delta = actual_dt,
#   exp_sigma = TRUE)$par

# tibble::tibble(obsSample =
#                  t_diffusion_strang_splitting_resid(stationary_part_estim_param,
#                  data = asinh(sim_res_t_distribution$X_t[sim_res_t_distribution$t < t_0]),
#                  delta = actual_dt)) |>
#                  ggplot2::ggplot(ggplot2::aes(sample = obsSample)) +
#                  ggplot2::geom_qq() + ggplot2::geom_qq_line()

# Dynamic part
# dynamic_part_true_param <- c(tau, true_param[1])
# dynamic_part_estim_param <- optimize_dynamic_likelihood(likelihood_fun = t_transform_dynamic_likelihood,
#                             data = asinh(sim_res_t_distribution$X_t[sim_res_t_distribution$t > t_0]),
#                             init_par = dynamic_part_true_param,
#                             delta = actual_dt,
#                             alpha0 = stationary_part_estim_param[1],
#                             mu0 = stationary_part_estim_param[2],
#                             sigma = stationary_part_estim_param[3])$par

# tibble::tibble(obsSample =
#                  t_transform_dynamic_likelihood_resid(dynamic_part_estim_param,
#                  data= asinh(sim_res_t_distribution$X_t[sim_res_t_distribution$t > t_0]),
#                  delta = actual_dt,
#                  alpha0 = stationary_part_estim_param[1],
#                  mu0 = stationary_part_estim_param[2],
#                  sigma = stationary_part_estim_param[3])) |>
#                  ggplot2::ggplot(ggplot2::aes(sample = obsSample)) +
#                  ggplot2::geom_qq() + ggplot2::geom_qq_line()

#-----------------------------------------------------------------------------------------------------------------------------#

## F-distributed stationary distribution model
# 
# actual_dt <- 1 / 12
# t_0 <- 50
# tau <- 100
# true_param <- c(0.75, 1.5, -2.5, 0.1)
# 
# F_sim_dynamic <- simulate_F_distribution_tipping_model(actual_dt, true_param, t_0 = t_0, tau = tau)
# 
# F_sim_dynamic |> ggplot(aes(x = t, y = X_t)) + geom_step()

# ## Stationary part
# # Parameters for stationary part
# mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) * sqrt(abs(true_param[3] / true_param[1]))
# alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
# 
# stationary_part_true_param <- c(alpha0, mu0, true_param[4])

# stationary_part_estim_param <- optimize_stationary_likelihood(
#   likelihood_fun = F_diffusion_strang_splitting,
#   data = 2 * asinh(sqrt(F_sim_dynamic$X_t[F_sim_dynamic$t<t_0])),
#   init_par = stationary_part_true_param,
#   delta = actual_dt,
#   exp_sigma = TRUE)$par

# tibble::tibble(obsSample =
#                  F_diffusion_strang_splitting_resid(stationary_part_estim_param,
#                  data = 2 * asinh(sqrt(F_sim_dynamic$X_t[F_sim_dynamic$t<t_0])),
#                  delta = actual_dt)) |>
#                  ggplot2::ggplot(ggplot2::aes(sample = obsSample)) +
#                  ggplot2::geom_qq() + ggplot2::geom_qq_line()

# ## Dynamic part
# dynamic_part_true_param <- c(tau, true_param[1])
# dynamic_part_estim_param <- optimize_dynamic_likelihood(likelihood_fun = F_transform_dynamic_likelihood,
#                             data = 2 * asinh(sqrt(F_sim_dynamic$X_t[F_sim_dynamic$t>t_0])),
#                             init_par = dynamic_part_true_param,
#                             delta = actual_dt,
#                             alpha0 = stationary_part_estim_param[1],
#                             mu0 = stationary_part_estim_param[2],
#                             sigma = stationary_part_estim_param[3])$par

# tibble::tibble(obsSample =
#                  F_transform_dynamic_likelihood_resid(dynamic_part_estim_param,
#                  data= 2 * asinh(sqrt(F_sim_dynamic$X_t[F_sim_dynamic$t>t_0])),
#                  delta = actual_dt,
#                  alpha0 = stationary_part_estim_param[1],
#                  mu0 = stationary_part_estim_param[2],
#                  sigma = stationary_part_estim_param[3])) |>
#                  ggplot2::ggplot(ggplot2::aes(sample = obsSample)) +
#                  ggplot2::geom_qq() + ggplot2::geom_qq_line()

#-----------------------------------------------------------------------------------------------------------------------------#

## Jacobi diffusion 
# true_param <- c(5, 0.1, -1, 0.15)
# actual_dt <- 1 / 12
# tau <- 100
# t_0 <- 50
# 
# sim_res_jacobi <- simulate_jacobi_diffusion_tipping_model(actual_dt, true_param, tau, t_0)
# sim_res_jacobi |> ggplot2::ggplot(ggplot2::aes(x = t, y = X_t)) +
#   ggplot2::geom_step() + ggplot2::geom_hline(yintercept = true_param[2], linetype = "dashed") +
#   ggplot2::geom_vline(xintercept = t_0)

# # Stationary part
# # Parameters for stationary part
# mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) * sqrt(abs(true_param[3] / true_param[1]))
# alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
# stationary_part_true_param <- c(alpha0, mu0, true_param[4])

# stationary_part_estim_param <- optimize_stationary_likelihood(
#   likelihood_fun = jacobi_diffusion_strang_splitting,
#   data = 2 * asin(sqrt(sim_res_jacobi$X_t[sim_res_jacobi$t<t_0])),
#   init_par = stationary_part_true_param,
#   delta = actual_dt,
#   exp_sigma = TRUE)$par

# tibble::tibble(obsSample =
#                  jacobi_diffusion_strang_splitting_resid(stationary_part_estim_param,
#                  data = 2 * asin(sqrt(sim_res_jacobi$X_t[sim_res_jacobi$t<t_0])),
#                  delta = actual_dt)) |>
#                  ggplot2::ggplot(ggplot2::aes(sample = obsSample)) +
#                  ggplot2::geom_qq() + ggplot2::geom_qq_line()

# dynamic_part_true_param <- c(tau, true_param[1])
# dynamic_part_estim_param <- optimize_dynamic_likelihood(likelihood_fun = jacobi_diffusion_transform_dynamic_likelihood,
#                             data = 2 * asin(sqrt(sim_res_jacobi$X_t[sim_res_jacobi$t>t_0])),
#                             init_par = dynamic_part_true_param,
#                             delta = actual_dt,
#                             alpha0 = stationary_part_estim_param[1],
#                             mu0 = stationary_part_estim_param[2],
#                             sigma = stationary_part_estim_param[3])$par

# tibble::tibble(obsSample =
#                  jacobi_diffusion_transform_dynamic_likelihood_resid(dynamic_part_estim_param,
#                  data= 2 * asin(sqrt(sim_res_jacobi$X_t[sim_res_jacobi$t>t_0])),
#                  delta = actual_dt,
#                  alpha0 = stationary_part_estim_param[1],
#                  mu0 = stationary_part_estim_param[2],
#                  sigma = stationary_part_estim_param[3])) |>
#                  ggplot2::ggplot(ggplot2::aes(sample = obsSample)) +
#                  ggplot2::geom_qq() + ggplot2::geom_qq_line()

