# Title: Model fitting of Stochastic Differential Equations with tipping points.
# Author: Anders Gantzhorn Kristensen (University of Copenhagen, andersgantzhorn@gmail.com)
# Date: 2024-02-22 (Last Updated: 2024-05-06)
#-----------------------------------------------------------------------------------------------------------------------------#
# Project: Tipping Point Estimation in Ecological Systems using Stochastic Differential Equations
# Description: This script implements optimizers, estimation methods and negative log-likelihood functions for estimation
# of tipping points in stochastic differential equations with additive-, square-root- and linear noise structure
#-----------------------------------------------------------------------------------------------------------------------------#
# License: MIT License (for more information, see LICENSE file in the repository).
# Dependencies: nleqslv, numDeriv.
# Usage: See example usage in the bottom of the script
# Acknowledgements: A few of the methods was made by me and Lucas Støjko for the mandatory project during the course:
# "NMAK23005U Inference for Stochastic Differential Equations" at the University of Copenhagen in 2023-2024.
# Additionally, some of the code is inspired by the code from Susanne - and Peter Ditlevsen's paper: 
# "Warning of a forthcoming collapse of the Atlantic meridional overturning circulation"
# Keywords: Stochastic Differential Equations, Likelihood methods, Estimation, Numerical Optimization

#-----------------------------------------------------------------------------------------------------------------------------#
# Helper functions
runge_kutta <- function(y0, h, func, n = 1) {
  h <- h / n
  for(i in 1:n) {
    k1 <- func(0, y0)
    k2 <- func(0.5 * h, y0 + 0.5 * h * k1)
    k3 <- func(0.5 * h, y0 + 0.5 * h * k2)
    k4 <- func(h, y0 + h * k3)
    
    y0 <- y0 + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
  }
  y0
}
#-----------------------------------------------------------------------------------------------------------------------------#
# General method that allows estimation allows almost completely numerical
# computation of the Strang method where the splitting is chosen such that
# the linear part is the linearization around the fixed point of the drift

numerical_strang_splitting <- function(par,
                                      data,
                                      delta,
                                      drift_lamperti_sde,
                                      exp_sigma = FALSE) {
  Xlow <- data[1:(length(data) - 1)]
  Xupp <- data[2:length(data)]
  
  beta  <- par[1]
  mu    <- par[2]
  if(exp_sigma){sigma <- exp(par[3])} else{sigma <- par[3]}
  
  # Median is chosen for the initial point.
  # Is calculated at every function call - optimization might be possible
  b <- nleqslv::nleqslv(x = median(data), fn = drift_lamperti_sde, par = par)$x
  
  A <- numDeriv::grad(func = drift_lamperti_sde, x = b, par = par)
  
  diff_f <- function(t, y){drift_lamperti_sde(y, par) - A * (y - b)}

  # Solution to ODE
  f_h <- runge_kutta(Xlow, delta / 2, diff_f, n = 1)
  
  mu_f <- exp(A * delta) * (f_h - b) + b
  sd_f <- sigma * sqrt(expm1(2 * A * delta) / (2 * A))
  
  # Inverse of non-linear ODE
  inv_f <- runge_kutta(Xupp, -delta / 2, diff_f, n = 1)
  
  # Derivative of inverse using Richardson Extrapolation.
  inv_f2 <- runge_kutta(Xupp + 0.01, -delta / 2, diff_f, n = 1)
  inv_f3 <- runge_kutta(Xupp - 0.01, -delta / 2, diff_f, n = 1)
  df     <- (inv_f2 - inv_f3) / (2 * 0.01)

  # Strang likelihood
  -mean(stats::dnorm(inv_f, mean = mu_f, sd = sd_f, log = TRUE), na.rm = TRUE) - 
    mean(log(abs(df)), na.rm = TRUE)
}

# Implementation of methods that are based on the additive noise term process
OU_init_params <- function(data, delta){
  N <- length(data)
  Xlow <- data[1:(length(data) - 1)]
  Xupp <- data[2:length(data)]
  
  mu_hat <- mean(data)
  
  rho_hat <- sum((Xupp- mu_hat) * (Xlow - mu_hat)) / 
    sum((Xlow - mu_hat)^2)
  
  alpha_hat <- - log(rho_hat) / delta
  
  gamma_sq_hat <- mean((Xupp - Xlow * rho_hat - mu_hat * (1 - rho_hat))^2) /
    (1 - rho_hat^2)
  
  sigma_hat <- sqrt(gamma_sq_hat * 2 * alpha_hat)
  
  c(alpha_hat, mu_hat, sigma_hat)
}


OU_likelihood <- function(par, data, delta, exp_sigma = FALSE){
  alpha0 <- par[1]
  mu0    <- par[2]
  if(exp_sigma){sigma <- exp(par[3])} else{sigma <- par[3]}

  N      <- length(data)
  
  Xupp <- data[2:N]
  Xlow <- data[1:(N - 1)]
  
  gamma2 <- sigma^2 / (2 * alpha0)
  rho0   <- exp(-alpha0 * delta)
  
  v_part <- gamma2 * (1 - rho0^2)
  m_part <- Xlow * rho0 + mu0 * (1 - rho0)
  
  # Negative log-likelihood
  -sum(dnorm(x = Xupp, mean = m_part, sd = sqrt(v_part), log = TRUE), na.rm = TRUE)
}

OU_score <- function(par, data, delta, exp_sigma = FALSE){
  Xlow <- data[1:(length(data) - 1)]
  Xupp <- data[2:length(data)]
  N  <- length(Xlow)
  
  alpha0 <- par[1]
  mu0    <- par[2]
  if(exp_sigma){sigma <- exp(par[3])} else{sigma <- par[3]}
  
  gamma2 <- sigma^2 / (2 * alpha0)
  rho          <- exp(-alpha0 * delta)
  
  eq1 <- (1 - rho) / (gamma2 * (1 - rho^2)) * sum(Xupp - Xlow * rho - mu0 * (1 - rho))
  
  eq2 <- N * rho / (1 - rho^2) + 
    sum((Xupp - Xlow * rho - mu0 * (1 - rho)) * (Xlow - mu0)) / (gamma2 * (1 - rho^2)) - 
    rho * sum((Xupp - Xlow * rho - mu0 * (1 - rho))^2) / (gamma2 * (1 - rho^2)^2)
  
  eq3 <- - N / (2 * gamma2) + sum((Xupp - Xlow * rho - mu0 * (1 - rho))^2) /
    (2 * gamma2^2 * (1 - rho^2))
  
  c(eq1, eq2, eq3)
}


OU_dynamic_likelihood <-  function(par, data, delta, 
                                   alpha0, mu0, sigma,
                                   pen = 0){
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
  
  fh_half     <- (mu_seq * fh_half_tmp_low + Xlow) / 
                    (fh_half_tmp_low + 1)
  
  fh_half_inv <- (mu_seq * fh_half_tmp_upp - Xupp) /
                    (fh_half_tmp_upp - 1)
  
  mu_h        <- fh_half * rho_seq + mu_seq * (1 - rho_seq)
  
  v_part           <- gamma2_seq * (1 - rho_seq^2)
  det_Dfh_half_inv <- 1 / (fh_half_tmp_upp-1)^2
  
  # Negative log-likelihood
  -sum(stats::dnorm(fh_half_inv, mu_h, sqrt(v_part), log = TRUE), na.rm = TRUE) -
    sum(log(abs(det_Dfh_half_inv)), na.rm = TRUE) +
    pen * N * (1/A - 1) * (abs(A) < 1) # penalty for discussion part
}

#-----------------------------------------------------------------------------------------------------------------------------#
# Implementation of methods that are based on the square-root noise term process

CIR_quadratic_martingale <- function(data, delta) {
  N <- length(data)
  M <- N - 1
  
  Xlow <- data[1:M]
  Xupp <- data[2:N]
  
  S1    <- sum(Xupp)
  S2    <- sum(1 / Xlow)
  S3    <- sum(Xupp / Xlow)
  ebdel <- (M * S3 - S1 * S2) / (M^2 - S1 * S2)
  beta  <- -log(ebdel) / delta
  
  mu <- (S3 - M * ebdel) / ((1 - ebdel) * S2)
  
  S4 <- sum((Xupp - mu + ebdel * (mu - Xlow))^2 / Xlow)
  S5 <- sum((Xupp * (ebdel - ebdel^2) / beta + mu * (1 - ebdel)^2 / (2 * beta)) / Xupp)
  
  sigma <- S4 / S5
  
  c(beta, mu, sqrt(sigma))
}

CIR_strang_splitting <- function(par, data, delta,
                                 exp_sigma = FALSE) {
  Xlow <- data[1:(length(data) - 1)]
  Xupp <- data[2:length(data)]
  
  beta  <- par[1]
  mu    <- par[2]
  if(exp_sigma){sigma <- exp(par[3])} else{sigma <- par[3]}
  
  if((4 * beta * mu - sigma^2) < 0 | sigma < 0){return(100000)}
  
  A <- - beta
  b <- 2 * sqrt((beta * mu - sigma^2 / 2) / beta)
  
  diff_f <- function(t, y){beta / 2 * y + (2 * beta * mu - sigma^2) / y -
      2 * sqrt(beta * (beta * mu - sigma^2 / 2))}
  
  # Solution to ODE
  f_h <- runge_kutta(Xlow, delta / 2, diff_f, n = 1)


  mu_f <- exp(A * delta) * (f_h - b) + b
  sd_f <- sigma * sqrt(expm1(2 * A * delta) / (2 * A))
  
  # Inverse of non-linear ODE
  inv_f <- runge_kutta(Xupp, -delta / 2, diff_f, n = 1)
  
  # Derivative of inverse using Richardson Extrapolation.
  inv_f2 <- runge_kutta(Xupp + 0.01, -delta / 2, diff_f, n = 1)
  inv_f3 <- runge_kutta(Xupp - 0.01, -delta / 2, diff_f, n = 1)
  df     <- (inv_f2 - inv_f3) / (2 * 0.01)
  
  # Strang likelihood
  -sum(stats::dnorm(inv_f, mean = mu_f, sd = sd_f, log = TRUE), na.rm = TRUE) -
    sum(log(abs(df)), na.rm = TRUE)
}


CIR_dynamic_likelihood <- function(par, data, delta,
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
  
  fh_half     <- (mu_seq * fh_half_tmp_low + Xlow) / (fh_half_tmp_low+1)
  
  fh_half_inv <- (mu_seq * fh_half_tmp_upp - Xupp) / (fh_half_tmp_upp-1)
  
  det_Dfh_half_inv <- 1 / (fh_half_tmp_upp-1)^2
  phi_dot <- fh_half / alpha_seq * (exp(-alpha_seq * delta) - exp(-2 * alpha_seq * delta)) +
    mu_seq / (2 * alpha_seq) * (-expm1(-alpha_seq * delta))^2

  if(any(phi_dot < 0) | any(is.na(phi_dot))){return(100000 + runif(1, min = -5, max = 1))}
  
  sd.part  <- sigma * sqrt(phi_dot)
  
  mu.part  <- exp(-alpha_seq * delta) * (fh_half - mu_seq) + mu_seq
  
  -sum(stats::dnorm(fh_half_inv, mu.part, sd.part, log = TRUE), na.rm = TRUE) -
    sum(log(abs(det_Dfh_half_inv)), na.rm = TRUE)
}

CIR_transform_dynamic_likelihood <- function(par, data, delta,
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
  
  # Derivative of inverse using Richardson Extrapolation.
  inv_f2 <- runge_kutta(Xupp + 0.01, -delta / 2, diff_f, n = 1)
  inv_f3 <- runge_kutta(Xupp - 0.01, -delta / 2, diff_f, n = 1)
  df     <- (inv_f2 - inv_f3) / (2 * 0.01)
  # Strang likelihood
  -sum(stats::dnorm(inv_f, mean = mu_f, sd = sd_f, log = TRUE), na.rm = TRUE) -
    sum(log(abs(df)), na.rm = TRUE)
}

#-----------------------------------------------------------------------------------------------------------------------------#
# Implementation of methods that are based on the linear noise term process
mean_reverting_GBM_martingale <- function(par, data, delta){
  Xlow <- data[1:(length(data) - 1)]
  Xupp <- data[2:length(data)]
  
  beta  <- par[1]
  mu    <- par[2]
  sigma <- par[3]
  
  F_i   <- exp(-beta * delta) * (Xlow - mu) + mu

  phi_i <- exp(-(2 * beta - sigma^2) * delta) *
    (Xlow^2 - (2 * beta * mu)/(beta - sigma^2) * Xlow + 2 * beta^2 * mu^2/ ((beta - sigma^2) * (2 * beta - sigma^2))) + 
    2 * beta * mu / (beta - sigma^2) * (exp(-beta * delta) * (Xlow - mu) + mu) - 
    2 * beta^2 * mu^2 / ((beta - sigma^2) * (2 * beta - sigma^2)) - 
    exp(-2 * beta * delta) * (Xlow - mu)^2 - mu^2 - 2 * mu * exp(-beta * delta) * (Xlow - mu)
  
  
  eq1 <- beta / sigma^2 * sum(1 / Xlow * (Xupp - F_i))
  eq2 <- - 1 / sigma^2 * sum(1 / Xlow * (Xupp - F_i))
  eq3 <- 1 / sigma^3 * sum(1 / Xlow^2 * ((Xupp - F_i)^2 - phi_i))
  
  c(eq1, eq2, eq3)
}

mean_reverting_GBM_alt_strang <- function(par, data, delta,
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
  -sum(stats::dlnorm(inv_f, meanlog = mu_log, sdlog = sd_log, log = TRUE), na.rm = TRUE) -
    beta * delta * ((length(data) - 1)) / 2
}

mean_reverting_GBM_strang <- function(par, data, delta,
                                      exp_sigma = FALSE) {
  Xlow <- data[1:(length(data) - 1)]
  Xupp <- data[2:length(data)]
  
  beta  <- par[1]
  mu    <- par[2]
  if(exp_sigma){sigma <- exp(par[3])} else{sigma <- par[3]}
  
  if((sigma^2 + 2 * beta) / (2 * beta * mu) < 0 |
     (2 * beta + sigma^2) / (2 * beta * mu) < 0){return(50000)}
  
  A <- - (sigma^2 + 2 * beta) / 2
  b <- - log((2 * beta + sigma^2) / (2 * beta * mu))
  
  diff_f <- function(t, y){-beta + beta * mu * exp(-y) - sigma^2 / 2 - A * (y - b)}

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
  -sum(stats::dnorm(inv_f, mean = mu_f, sd = sd_f, log = TRUE), na.rm = TRUE) -
    sum(log(abs(df)), na.rm = TRUE)
}

mean_reverting_GBM_dynamic_likelihood <- function(par, data, delta,
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
  
  fh_half          <- (mu_seq * fh_half_tmp_low + Xlow) / (fh_half_tmp_low + 1)
  
  fh_half_inv      <- (mu_seq * fh_half_tmp_upp - Xupp) / (fh_half_tmp_upp - 1)
  
  det_Dfh_half_inv <- 1 / (fh_half_tmp_upp - 1)^2
  
  sd.part          <- exp(-(alpha_seq - sigma^2) * delta) * sqrt(
    (fh_half - mu_seq)^2 * (1 - exp(-sigma^2 * delta)) - 
      2 * mu_seq * sigma^2 / (alpha_seq - sigma^2) *
          (fh_half - mu_seq) * (1 - exp(-(sigma^2 - alpha_seq) * delta)) -
      mu_seq^2 * sigma^2 / (2 * alpha_seq - sigma^2) *
          (1 - exp(-(sigma^2 - 2*alpha_seq) * delta))
  )
  
  mu.part <- exp(-alpha_seq * delta) * (fh_half - mu_seq) + mu_seq
  
  -sum(stats::dnorm(fh_half_inv, mu.part, sd.part, log = TRUE), na.rm = TRUE) - 
    sum(log(abs(det_Dfh_half_inv)), na.rm = TRUE)
}

mean_reverting_GBM_transform_dynamic_likelihood <- function(par, data, delta,
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
  
  sqrt_argument <- pmax(sigma^4 / 4 - A * (2 * m * sigma^2 + 4 * lam_seq), 0.001)
  
  exp_b <- max((2 * m * A - 1 / 2 * sigma^2 +
                  sqrt(sqrt_argument)) / (2 * A), 0.001)
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
  
  # Derivative of inverse using Richardson Extrapolation.
  inv_f2 <- runge_kutta(Xupp + 0.01, -delta / 2, diff_f, n = 1)
  inv_f3 <- runge_kutta(Xupp - 0.01, -delta / 2, diff_f, n = 1)
  df     <- (inv_f2 - inv_f3) / (2 * 0.01)
  
  # Strang likelihood
 -sum(stats::dnorm(inv_f, mean = mu_f, sd = sd_f, log = TRUE), na.rm = TRUE) - 
    sum(log(abs(df)), na.rm = TRUE)
}

#-----------------------------------------------------------------------------------------------------------------------------#
# Implementation of methods that are based on the process with stationary t-distribution
t_diffusion_strang_splitting <- function(par, data, delta,
                                         exp_sigma = FALSE) {
  Xlow <- data[1:(length(data) - 1)]
  Xupp <- data[2:length(data)]
  
  beta  <- par[1]
  mu    <- par[2]
  if(exp_sigma){sigma <- exp(par[3])} else{sigma <- par[3]}
  
  asinh_term <- asinh(2 * mu * beta / (2 * beta + sigma^2))
  
  diff_f <- function(t, y){(y - tanh(y) - asinh_term) *
      (beta + 1 / 2 * sigma^2) + mu * beta / cosh(y)}
  
  # Solution to ODE
  f_h <- runge_kutta(Xlow, delta / 2, diff_f, n = 1)
  
  A <- - (beta + 1 / 2 * sigma^2)
  b <- asinh_term
  
  mu_f <- exp(A * delta) * (f_h - b) + b
  sd_f <- sigma * sqrt(expm1(2 * A * delta) / (2 * A))
  
  # Inverse of non-linear ODE
  inv_f <- runge_kutta(Xupp, -delta / 2, diff_f, n = 1)
  
  # Derivative of inverse using Richardson Extrapolation.
  inv_f2 <- runge_kutta(Xupp + 0.01, -delta / 2, diff_f, n = 1)
  inv_f3 <- runge_kutta(Xupp - 0.01, -delta / 2, diff_f, n = 1)
  df     <- (inv_f2 - inv_f3) / (2 * 0.01)
  
  # Strang likelihood
  -sum(stats::dnorm(inv_f, mean = mu_f, sd = sd_f, log = TRUE), na.rm = TRUE) -
    sum(log(abs(df)), na.rm = TRUE)
}

t_transform_dynamic_likelihood <- function(par, data, delta, alpha0, mu0, sigma){
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
  
  if(any((sigma^4 - 8 * A * (2 * lam_seq + m * sigma^2)) < 0) |
     any(is.na(sigma^4 - 8 * A * (2 * lam_seq + m * sigma^2)))){
    return(50000)
    }
  
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
  
  # Derivative of inverse using Richardson Extrapolation.
  inv_f2 <- runge_kutta(Xupp + 0.01, -delta / 2, diff_f, n = 1)
  inv_f3 <- runge_kutta(Xupp - 0.01, -delta / 2, diff_f, n = 1)
  df     <- (inv_f2 - inv_f3) / (2 * 0.01)
  
  # Strang likelihood
  -sum(stats::dnorm(inv_f, mean = mu_f, sd = sd_f, log = TRUE), na.rm = TRUE) -
    sum(log(abs(df)), na.rm = TRUE)
}

#-----------------------------------------------------------------------------------------------------------------------------#
# Implementation of methods that are based on the process with stationary F-distribution

F_diffusion_strang_splitting <- function(par, data, delta, exp_sigma = FALSE) {
  Xlow <- data[1:(length(data) - 1)]
  Xupp <- data[2:length(data)]
  
  beta  <- par[1]
  mu    <- par[2]
  if(exp_sigma){sigma <- exp(par[3])} else{sigma <- par[3]}
  
  A <- - (beta + sigma^2 / 2)
  b <- acosh(2 * beta * (2 * mu + 1) / (2 * beta + sigma^2))
  
  diff_f <- function(t,y){-A * (y - b - cosh(y) / sinh(y)) + beta * (2 * mu + 1) / sinh(y)}
  # Solution to ODE
  f_h <- runge_kutta(Xlow, delta / 2, diff_f, n = 1)
  
  mu_f <- exp(A * delta) * (f_h - b) + b
  sd_f <- sigma * sqrt(expm1(2 * A * delta) / (2 * A))
  
  # Inverse of non-linear ODE
  inv_f <- runge_kutta(Xupp, -delta / 2, diff_f, n = 1)
  
  # Derivative of inverse using Richardson Extrapolation.
  inv_f2 <- runge_kutta(Xupp + 0.01, -delta / 2, diff_f, n = 1)
  inv_f3 <- runge_kutta(Xupp - 0.01, -delta / 2, diff_f, n = 1)
  df     <- (inv_f2 - inv_f3) / (2 * 0.01)

  # Strang likelihood
  -mean(stats::dnorm(inv_f, mean = mu_f, sd = sd_f, log = TRUE), na.rm = TRUE) -
    mean(log(abs(df)), na.rm = TRUE)
}

F_transform_dynamic_likelihood <- function(par, data, delta, alpha0, mu0, sigma){
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
  
  # Derivative of inverse using Richardson Extrapolation.
  inv_f2 <- runge_kutta(Xupp + 0.01, -delta / 2, diff_f, n = 1)
  inv_f3 <- runge_kutta(Xupp - 0.01, -delta / 2, diff_f, n = 1)
  df     <- (inv_f2 - inv_f3) / (2 * 0.01)
  
  # Strang likelihood
  -sum(stats::dnorm(inv_f, mean = mu_f, sd = sd_f, log = TRUE), na.rm = TRUE) -
    sum(log(abs(df)), na.rm = TRUE)
}

#-----------------------------------------------------------------------------------------------------------------------------#
# Implementation of methods that are based on the jacobi diffusion

jacobi_diffusion_strang_splitting <- function(par, data, delta, exp_sigma = FALSE) {
  Xlow <- data[1:(length(data) - 1)]
  Xupp <- data[2:length(data)]
  
  beta  <- par[1]
  mu    <- par[2]
  if(exp_sigma){sigma <- exp(par[3])} else{sigma <- par[3]}
  
  if(abs(beta * (2 * mu - 1) / (sigma^2 / 2 - beta))>1)
    {return(50000+runif(1, min = -1, max = 1))}
  
  A <- 1/2 * sigma^2 - beta
  
  b <- acos(beta * (2 * mu - 1) / (sigma^2 / 2 - beta))
  
  diff_f <- function(t, y){-1 / sin(y) * (beta * (1 - cos(y) - 2 * mu) +
                                            sigma^2 / 2 * cos(y)) - A * (y - b)}
  
  # Solution to ODE
  f_h <- runge_kutta(Xlow, delta / 2, diff_f, n = 1)
  
  #if(any(is.na(f_h))){return(50000)}
  
  mu_f <- exp(A * delta) * (f_h - b) + b
  sd_f <- sigma * sqrt(expm1(2 * A * delta) / (2 * A))
  
  # Inverse of non-linear ODE
  inv_f <- runge_kutta(Xupp, -delta / 2, diff_f, n = 1)
  
  # Derivative of inverse using Richardson Extrapolation.
  inv_f2 <- runge_kutta(Xupp + 0.01, -delta / 2, diff_f, n = 1)
  inv_f3 <- runge_kutta(Xupp - 0.01, -delta / 2, diff_f, n = 1)
  df     <- (inv_f2 - inv_f3) / (2 * 0.01)
  
  # Strang likelihood
  -sum(stats::dnorm(inv_f, mean = mu_f, sd = sd_f, log = TRUE), na.rm = TRUE) -
    sum(log(abs(df)), na.rm = TRUE)
}

jacobi_diffusion_transform_dynamic_likelihood <- function(par, data, delta, alpha0, mu0, sigma){
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
  
  # Derivative of inverse using Richardson Extrapolation.
  inv_f2 <- runge_kutta(Xupp + 0.01, -delta / 2, diff_f, n = 1)
  inv_f3 <- runge_kutta(Xupp - 0.01, -delta / 2, diff_f, n = 1)
  df     <- (inv_f2 - inv_f3) / (2 * 0.01)
  
  # Strang likelihood
  -sum(stats::dnorm(inv_f, mean = mu_f, sd = sd_f, log = TRUE), na.rm = TRUE) -
    sum(log(abs(df)), na.rm = TRUE)
}

#-----------------------------------------------------------------------------------------------------------------------------#
## General optimizers for stationary - and dynamic part of process
optimize_stationary_likelihood <- function(likelihood_fun,
                                           data,
                                           init_par, 
                                           delta,
                                           exp_sigma = TRUE,
                                           method = "BFGS",
                                           return_all = FALSE,
                                           ...){
  res_optim <- stats::optim(init_par, fn = likelihood_fun,
                            method = method,
                            data = data, delta = delta, exp_sigma = exp_sigma, ...)
  
  if(exp_sigma){
    res_optim$par[3] <- exp(res_optim$par[3])
  }
  
  if(return_all) {
    return(res_optim)
  } else {
    return(list(par = res_optim$par, objective = res_optim$value))
  }
}

optimize_dynamic_likelihood <- function(likelihood_fun, data,
                                        init_par, delta,
                                        alpha0, mu0, sigma, 
                                        method = "BFGS", 
                                        return_all = FALSE,...){

  res_optim <- stats::optim(init_par, fn = likelihood_fun,
                            method = method,
                            data = data, delta = delta,
                            alpha0 = alpha0 , mu0 = mu0, sigma = sigma, ...)
  
  if(return_all) {
    return(res_optim)
  } else {
    return(list(par = res_optim$par, objective = res_optim$value))
  }
}

optimize_dynamic_simulation_likelihood <- function(likelihood_fun, data, times, M, N,
                                                    init_par, alpha0, mu0, sigma, t_0){
  res_optim <- optim(par = init_par, fn = likelihood_fun,
        method = "L-BFGS-B", data = data,
        times = times, M = M, N = N,
        alpha0 = alpha0, mu0 = mu0,
        sigma = sigma, t_0 = t_0,
        lower = rep(0, length(init_par)),
        upper = rep(Inf, length(init_par)))$par
  res_optim
}

#-----------------------------------------------------------------------------------------------------------------------------#
# Example usage
# source("source_code/tipping_simulations.R")

# Additive noise model
# true_param <- c(0.87, -1.51, -2.69, 0.2, 0.5)
# actual_dt <- 0.005
# tau <- 150
# t_0 <- 50
# sim_res_add <- simulate_additive_noise_tipping_model(actual_dt, true_param, tau, t_0)
# sample_n(sim_res_add, min(nrow(sim_res_add), 10000)) |> ggplot(aes(x = t, y = X_t)) + geom_step() +
#   geom_line(aes(y = mu_t))
# Stationary part
# Parameters for stationary part
# mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) *  sqrt(abs(true_param[3] / true_param[1]))
# alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
# stationary_part_true_param <- c(alpha0, mu0, true_param[4])
# 
# optimize_stationary_likelihood(likelihood_fun = OU_likelihood,
#                                data = sim_res_add$X_t[sim_res_add$t < t_0],
#                                init_par = stationary_part_true_param,
#                                delta = actual_dt, exp_sigma = FALSE)

# 
# nleqslv::nleqslv(x = stationary_part_true_param,
#                  fn = OU_score,
#                  data = sim_res_add$X_t[sim_res_add$t < t_0],
#                  delta = actual_dt)$x


## Dynamic part
# dynamic_part_true_param <- c(tau, true_param[1], 1)
# 
# 
# optimize_dynamic_likelihood(likelihood_fun = OU_dynamic_likelihood,
#                             data = sim_res_add$X_t[sim_res_add$t > t_0],
#                             init_par = dynamic_part_true_param,
#                             delta = actual_dt,
#                             alpha0 = stationary_part_true_param[1],
#                             mu0 = stationary_part_true_param[2],
#                             sigma = stationary_part_true_param[3])



#-----------------------------------------------------------------------------------------------------------------------------#

# # Square-root noise model
# true_param <- c(0.8, 2, -1, 0.05)
# actual_dt <- 0.005
# tau <- 120
# t_0 <- 50
# sim_res_sqrt <- simulate_squareroot_noise_tipping_model(actual_dt, true_param, tau, t_0)
# sample_n(sim_res_sqrt, min(nrow(sim_res_sqrt), 10000)) |> ggplot(aes(x = t, y = X_t)) + geom_step()

## Stationary part
## Parameters for stationary part
# mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) * sqrt(abs(true_param[3] / true_param[1]))
# alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
# stationary_part_true_param <- c(alpha0, mu0, true_param[4])
# 
# CIR_quadratic_martingale(sim_res_sqrt$X_t[sim_res_sqrt$t < t_0], actual_dt)
# 
# 
# optimize_stationary_likelihood(likelihood_fun = CIR_strang_splitting,
#                                data = 2 * sqrt(sim_res_sqrt$X_t[sim_res_sqrt$t < t_0]),
#                                init_par = stationary_part_true_param,
#                                delta = actual_dt, exp_sigma = FALSE)
#
## Dynamic part
# 
# dynamic_part_true_param <- c(tau, true_param[1])
# optimize_dynamic_likelihood(likelihood_fun = CIR_dynamic_likelihood,
#                             data = sim_res_sqrt$X_t[sim_res_sqrt$t > t_0],
#                             init_par = dynamic_part_true_param,
#                             delta = actual_dt,
#                             alpha0 = stationary_part_true_param[1],
#                             mu0 = stationary_part_true_param[2],
#                             sigma = stationary_part_true_param[3])
# 
# optimize_dynamic_likelihood(likelihood_fun = CIR_transform_dynamic_likelihood,
#                             data = 2 * sqrt(sim_res_sqrt$X_t[sim_res_sqrt$t > t_0]),
#                             init_par = dynamic_part_true_param,
#                             delta = actual_dt,
#                             alpha0 = stationary_part_true_param[1],
#                             mu0 = stationary_part_true_param[2],
#                             sigma = stationary_part_true_param[3],
#                             control = list(reltol = sqrt(.Machine$double.eps) / 1000))
#-----------------------------------------------------------------------------------------------------------------------------#

## Linear noise model
# true_param <- c(0.2, 3, -1, 0.2)
# actual_dt <- 0.1
# tau <- 100
# t_0 <- 50
# sim_res_linear <- simulate_linear_noise_tipping_model(actual_dt, true_param, tau, t_0)
# sample_n(sim_res_linear, min(nrow(sim_res_linear), 10000)) |> ggplot(aes(x = t, y = X_t)) + geom_step()

## Stationary part
## Parameters for stationary part
# mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) * sqrt(abs(true_param[3] / true_param[1]))
# alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
# stationary_part_true_param <- c(alpha0, mu0, true_param[4])
# 
# nleqslv::nleqslv(x = stationary_part_true_param, fn = mean_reverting_GBM_martingale,
#                  data = sim_res_linear$X_t[sim_res_linear$t < t_0],
#                  delta = actual_dt)$x
#
# optimize_stationary_likelihood(mean_reverting_GBM_strang, log(sim_res_linear$X_t[sim_res_linear$t<t_0]),
#                                init_par = stationary_part_true_param, delta = actual_dt,
#                                exp_sigma = TRUE)
# 
# optimize_stationary_likelihood(mean_reverting_GBM_alt_strang,
#                                sim_res_linear$X_t[sim_res_linear$t<t_0],
#                               init_par = stationary_part_true_param, delta = actual_dt,
#                               exp_sigma = TRUE)

## Dynamic part
# dynamic_part_true_param <- c(tau, true_param[1])
# optimize_dynamic_likelihood(likelihood_fun = mean_reverting_GBM_dynamic_likelihood,
#                             data = sim_res_linear$X_t[sim_res_linear$t > t_0],
#                             init_par = dynamic_part_true_param,
#                             delta = actual_dt,
#                             alpha0 = stationary_part_true_param[1],
#                             mu0 = stationary_part_true_param[2],
#                             sigma = stationary_part_true_param[3],
#                             control = list(reltol = sqrt(.Machine$double.eps) / 100))
# 
# optimize_dynamic_likelihood(likelihood_fun = mean_reverting_GBM_transform_dynamic_likelihood,
#                             data = log(sim_res_linear$X_t[sim_res_linear$t > t_0]),
#                             init_par = dynamic_part_true_param,
#                             delta = actual_dt,
#                             alpha0 = stationary_part_true_param[1],
#                             mu0 = stationary_part_true_param[2],
#                             sigma = stationary_part_true_param[3],
#                             control = list(reltol = sqrt(.Machine$double.eps)))

#-----------------------------------------------------------------------------------------------------------------------------#

## t-distributed stationary distribution model
#
# true_param <- c(0.5, -2, -3, 0.35)
# actual_dt <- 0.01
# tau <- 80
# t_0 <- 50
# sim_res_t_distribution <- simulate_t_distribution_tipping_model(actual_dt, true_param, tau, t_0)
# sim_res_t_distribution |> ggplot2::ggplot(ggplot2::aes(x = t, y = X_t)) +
#   ggplot2::geom_step() + ggplot2::geom_hline(yintercept = true_param[2], linetype = "dashed") +
#   ggplot2::geom_vline(xintercept = t_0)

# # Stationary part
# # Parameters for stationary part
# mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) * sqrt(abs(true_param[3] / true_param[1]))
# alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
# stationary_part_true_param <- c(alpha0, mu0, true_param[4])
# 
# optimize_stationary_likelihood(
#               likelihood_fun = t_diffusion_strang_splitting,
#               data = asinh(sim_res_t_distribution$X_t[sim_res_t_distribution$t < t_0]),
#               init_par = stationary_part_true_param,
#               delta = actual_dt,
#               exp_sigma = TRUE,
#               control = list(reltol = sqrt(.Machine$double.eps) / 1000))

# Dynamic part
# dynamic_part_true_param <- c(tau, true_param[1])
#
# optimize_dynamic_likelihood(likelihood_fun = t_transform_dynamic_likelihood,
#                   data = asinh(sim_res_t_distribution$X_t[sim_res_t_distribution$t > t_0]),
#                   init_par = dynamic_part_true_param,
#                   delta = actual_dt,
#                   alpha0 = stationary_part_true_param[1],
#                   mu0 = stationary_part_true_param[2],
#                   sigma = stationary_part_true_param[3],
#                   control = list(reltol = sqrt(.Machine$double.eps) / 10))

#-----------------------------------------------------------------------------------------------------------------------------#

## F-distributed stationary distribution model
# 
# true_param <- c(0.5, 1, -2, 0.1)
# actual_dt <- 0.01
# t_0 <- 100
# tau <- 100
# F_sim_dynamic <- simulate_F_distribution_tipping_model(actual_dt, true_param, t_0 = t_0, tau = tau)
# F_sim_dynamic |> ggplot(aes(x = t, y = X_t)) + geom_step()
# ## Stationary part
# # Parameters for stationary part
# mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) * sqrt(abs(true_param[3] / true_param[1]))
# alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
# 
# stationary_part_true_param <- c(alpha0, mu0, true_param[4])
# 
# optimize_stationary_likelihood(
#   likelihood_fun = F_diffusion_strang_splitting,
#   data = 2 * asinh(sqrt(F_sim_dynamic$X_t[F_sim_dynamic$t<t_0])),
#   init_par = stationary_part_true_param,
#   delta = actual_dt,
#   exp_sigma = TRUE)

# ## Dynamic part
# dynamic_part_true_param <- c(tau, true_param[1])
# optimize_dynamic_likelihood(likelihood_fun = F_transform_dynamic_likelihood,
#                             data = 2 * asinh(sqrt(F_sim_dynamic$X_t[F_sim_dynamic$t>t_0])),
#                             init_par = dynamic_part_true_param,
#                             delta = actual_dt,
#                             alpha0 = stationary_part_true_param[1],
#                             mu0 = stationary_part_true_param[2],
#                             sigma = stationary_part_true_param[3],
#                             control = list(reltol = sqrt(.Machine$double.eps) / 100))
# 
#-----------------------------------------------------------------------------------------------------------------------------#

## Jacobi diffusion 
# true_param <- c(4, 0.1, -2, 0.3)
# actual_dt <- 0.01
# tau <- 100
# t_0 <- 50
# 
# sim_res_jacobi <- simulate_jacobi_diffusion_tipping_model(actual_dt, true_param, tau, t_0)
# sim_res_jacobi |> ggplot2::ggplot(ggplot2::aes(x = t, y = X_t)) +
#   ggplot2::geom_step() + ggplot2::geom_hline(yintercept = true_param[2], linetype = "dashed") +
#   ggplot2::geom_vline(xintercept = t_0)

## Stationary part
## Parameters for stationary part
# mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) * sqrt(abs(true_param[3] / true_param[1]))
# alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
# stationary_part_true_param <- c(alpha0, mu0, true_param[4])

# optimize_stationary_likelihood(
#   likelihood_fun = jacobi_diffusion_strang_splitting,
#   data = 2 * asin(sqrt(sim_res_jacobi$X_t[sim_res_jacobi$t<t_0])),
#   init_par = stationary_part_true_param,
#   delta = actual_dt,
#   exp_sigma = TRUE)

# dynamic_part_true_param <- c(tau, true_param[1])
# optimize_dynamic_likelihood(likelihood_fun = jacobi_diffusion_transform_dynamic_likelihood,
#                             data = 2 * asin(sqrt(sim_res_jacobi$X_t[sim_res_jacobi$t>t_0])),
#                             init_par = dynamic_part_true_param,
#                             delta = actual_dt,
#                             alpha0 = stationary_part_true_param[1],
#                             mu0 = stationary_part_true_param[2],
#                             sigma = stationary_part_true_param[3])

