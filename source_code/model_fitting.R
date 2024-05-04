# Title: Model fitting of Stochastic Differential Equations with tipping points.
# Author: Anders Gantzhorn Kristensen (University of Copenhagen, andersgantzhorn@gmail.com)
# Date: 2024-02-22 (Last Updated: 2024-04-09)
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
runge_kutta <- function(y0, h, f, n = 1) {
  h <- h / n
  
  for(i in 1:n) {
    k1 <- f(0, y0)
    k2 <- f(0.5 * h, y0 + 0.5 * h * k1)
    k3 <- f(0.5 * h, y0 + 0.5 * h * k2)
    k4 <- f(h, y0 + h * k3)
    
    y0 <- y0 + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
  }
  y0
}
#-----------------------------------------------------------------------------------------------------------------------------#

# General method that allows estimation allows almost completely numeric computation of the Strang method
# where the splitting is chosen such that the linear part is the linearization around the fixed point of the drift
numeric_strang_splitting <- function(par, data, delta, drift_lamperti_sde, exp_sigma = FALSE) {
  x0 <- data[1:(length(data) - 1)]
  x1 <- data[2:length(data)]
  
  beta  <- par[1]
  mu    <- par[2]
  if(exp_sigma){sigma <- exp(par[3])} else{sigma <- par[3]}
  
  
  b <- nleqslv::nleqslv(x = mu, fn = drift_lamperti_sde, par = par)$x
  
  A <- numDeriv::grad(func = drift_lamperti_sde, x = b, par = par)
  
  diff_f <- function(t, y){drift_lamperti_sde(y, par) - A * (y - b)}
  
  # Solution to ODE
  f_h <- runge_kutta(x0, delta / 2, diff_f, n = 1)
  
  
  mu_f <- exp(A * delta) * (f_h - b) + b
  sd_f <- sigma * sqrt(expm1(2 * A * delta) / (2 * A))
  
  # Inverse of non-linear ODE
  inv_f <- runge_kutta(x1, -delta / 2, diff_f, n = 1)
  
  # Derivative of inverse using Richardson Extrapolation.
  inv_f2 <- runge_kutta(x1 + 0.01, -delta / 2, diff_f, n = 1)
  inv_f3 <- runge_kutta(x1 - 0.01, -delta / 2, diff_f, n = 1)
  df     <- (inv_f2 - inv_f3) / (2 * 0.01)
  
  # Strang likelihood
  -sum(stats::dnorm(inv_f, mean = mu_f, sd = sd_f, log = TRUE), na.rm = TRUE) - 
    sum(log(abs(df)), na.rm = TRUE)
}

# Implementation of methods that are based on the additive noise term process
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
  x0 <- data[1:(length(data) - 1)]
  x1 <- data[2:length(data)]
  N  <- length(x0)
  
  alpha0 <- par[1]
  mu0    <-  par[2]
  if(exp_sigma){sigma <- exp(par[3])} else{sigma <- par[3]}
  
  gamma_square <- sigma^2 / (2 * alpha0)
  rho          <- exp(-alpha0 * delta)
  
  eq1 <- (1 - rho) / (gamma_square * (1 - rho^2)) * sum(x1 - x0 * rho - mu0 * (1 - rho))
  
  eq2 <- N * rho / (1 - rho^2) + 
    sum((x1 - x0 * rho - mu0 * (1 - rho)) * (x0 - mu0)) / (gamma_square * (1 - rho^2)) - 
    rho * sum((x1 - x0 * rho - mu0 * (1 - rho))^2) / (gamma_square * (1 - rho^2)^2)
  
  eq3 <- - N / (2 * gamma_square) + sum((x1 - x0 * rho - mu0 * (1 - rho))^2) / (2 * gamma_square^2 * (1 - rho^2))
  
  c(eq1, eq2, eq3)
}


OU_dynamic_likelihood <-  function(par, data, delta, alpha0, mu0, sigma){
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
    sum(log(abs(det_Dfh_half_inv)), na.rm = TRUE)
}

OU_dynamic_numeric_score <- function(par, data, delta, alpha0, mu0, sigma){
  numDeriv::grad(OU_dynamic_likelihood,
                 x = par,
                 delta = delta,
                 data = data,
                 alpha0 = alpha0,
                 mu0 = mu0,
                 sigma = sigma
  )
}

OU_dynamic_score <- function(par, data, delta, alpha0, mu0, sigma){
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
  
  #fh_half_tmp_low <- A * delta * (Xlow - mu_seq)
  #fh_half_tmp_upp <- A * delta * (Xupp - mu_seq)


  # eq1 <- sum(
  #     (1 - rho_seq) / (gamma2_seq * (1 - rho_seq^2)) * 
  #       (fh_half_tmp_upp - fh_half_tmp_low * rho_seq - mu_seq * (1 - rho_seq))
  # )

  # eq2 <- (sum(
  #   rho_seq / (1 - rho_seq^2)
  # ) + sum(
  #     (fh_half_tmp_upp - fh_half_tmp_low * rho_seq - mu_seq * (1 - rho_seq)) *
  #       (fh_half_tmp_low - mu_seq) / (gamma2_seq * (1 - rho_seq^2))
  #   ) - sum(
  #     rho_seq * (fh_half_tmp_upp - fh_half_tmp_low * rho_seq - mu_seq * 
  #                  (1 - rho_seq))^2 / (gamma2_seq * (1 - rho_seq^2)^2)
  #   )
  # ) / N

  # eq3 <- (-sum(
  #   1 / (2 * gamma2_seq)
  #   ) + sum(
  #   (fh_half_tmp_upp - fh_half_tmp_low * rho_seq - mu_seq * (1 - rho_seq))^2 /
  #     (2 * gamma2_seq^2 * (1 - rho_seq^2)))
  #   )
  # 

  ## Calculate the Strang-based score - divided by N for numeric stability
  eq1 <- sum(
      (1 - rho_seq) / (gamma2_seq * (1 - rho_seq^2)) *
        (Xupp - Xlow * rho_seq - mu_seq * (1 - rho_seq))
  ) / N
  
  eq2 <- (sum(
    rho_seq / (1 - rho_seq^2)
   ) + sum(
      (Xupp - Xlow * rho_seq - mu_seq * (1 - rho_seq)) * (Xlow - mu_seq) /
        (gamma2_seq * (1 - rho_seq^2))
    ) - sum(
      rho_seq * (Xupp - Xlow * rho_seq - mu_seq *(1 - rho_seq))^2 /
        (gamma2_seq * (1 - rho_seq^2)^2)
    )
  ) / N
  
  eq3 <- (-sum(
    1 / (2 * gamma2_seq)
    ) + sum(
    (Xupp - Xlow * rho_seq - mu_seq * (1 - rho_seq))^2 /
      (2 * gamma2_seq^2 * (1 - rho_seq^2)))
    ) / N

  
  if(length(par) == 3){return(c(eq1, eq2, eq3))} # If we want to estimate nu, use all three equations.
  c(eq1, eq3)
}


OU_dynamic_simulation_likelihood <- function(par, data, times, M, N, alpha0, mu0, sigma, t_0){
  tau <- par[1]
  A   <- par[2]
  nu  <- if(length(par) == 3) par[3] else 1
  
  m_tip   <- mu0 - alpha0/(2 * A)
  lambda0 <- -alpha0^2 / (4 * A)
  delta   <- 1 / M
  #time   <- delta * (1:(N-1))
  lam_seq <- lambda0 * (1 - (times[-length(times)]-t_0) / tau)^nu
  lam_mat <- t(matrix(rep(lam_seq, length.out = N * length(lam_seq)), nrow = length(lam_seq), ncol = N))
  # Number of data points
  numData <- length(data) - 1
  
  
  dW <- array(dqrng::dqrnorm(N * M * numData, mean = 0, sd = sqrt(delta)), dim = c(M, N, numData))
  # Initialize the process array
  X  <- array(NA, dim = c(M, N, numData))
  
  X[1, , ] <- rep(data[-length(data)], each = N)
  for (m in 2:M){
    X[m, , ] <- X[m - 1, , ] - (A * (X[m - 1, , ] - m_tip)^2 + lam_mat) * delta + sigma * dW[m - 1, , ]
  }
  
  second_to_last_X <- t(X[M - 1, , ])
  
  sigma_t   <- sigma * sqrt(delta)
  
  mu_t      <- data[1:numData] - (A * (second_to_last_X - m_tip)^2 + lam_seq) * delta
  
  densities <- dnorm(data[-1], mean = mu_t, sd = sigma_t, log = TRUE)

  loglik    <- -mean(colMeans(densities))
  if(is.na(loglik) | is.infinite(loglik)){
    print("NA Likelihood")
    return(1000000)
  }
  loglik
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

CIR_alt_strang_splitting <- function(par, data, delta, exp_sigma = FALSE) {
  x0 <- data[1:(length(data) - 1)]
  x1 <- data[2:length(data)]
  
  beta  <- par[1]
  mu    <- par[2]
  if(exp_sigma){sigma <- exp(par[3])} else{sigma <- par[3]}
  
  if(4 * beta * mu - sigma^2 < 0){return(100000)}
  
  diff_f <- function(t, y){y/2 * ((4 * beta * mu - sigma^2)/log(y) - beta * log(y))}
  
  f <- runge_kutta(x0, delta / 2, diff_f)

  mu_log <- log(f)

  sd_log <- sigma * sqrt(delta)
  
  inv_f  <- runge_kutta(x1, -delta / 2, diff_f)
  
  inv_f2 <- runge_kutta(x1 + 0.01, -delta / 2, diff_f)
  inv_f3 <- runge_kutta(x1 - 0.01, -delta / 2, diff_f)
  
  df <- (inv_f2 - inv_f3) / (2 * 0.01)
  
  -sum(stats::dlnorm(inv_f, meanlog = mu_log, sdlog = sd_log, log = TRUE), na.rm = TRUE) -
    sum(log(abs(df)), na.rm = TRUE)
}

CIR_strang_splitting <- function(par, data, delta, exp_sigma = FALSE) {
  x0 <- data[1:(length(data) - 1)]
  x1 <- data[2:length(data)]
  
  beta  <- par[1]
  mu    <- par[2]
  if(exp_sigma){sigma <- exp(par[3])} else{sigma <- par[3]}
  
  if((4 * beta * mu - sigma^2) < 0 | sigma < 0){return(100000)}
  diff_f <- function(t, y){beta / 2 * y + (2 * beta * mu - sigma^2) / y - 2 * sqrt(beta * (beta * mu - sigma^2 / 2))}
  
  # Solution to ODE
  f_h <- runge_kutta(x0, delta / 2, diff_f, n = 1)
  A <- - beta
  b <- 2 * sqrt((beta * mu - sigma^2 / 2) / beta)

  mu_f <- exp(A * delta) * (f_h - b) + b
  sd_f <- sigma * sqrt(expm1(2 * A * delta) / (2 * A))
  
  # Inverse of non-linear ODE
  inv_f <- runge_kutta(x1, -delta / 2, diff_f, n = 1)
  
  # Derivative of inverse using Richardson Extrapolation.
  inv_f2 <- runge_kutta(x1 + 0.01, -delta / 2, diff_f, n = 1)
  inv_f3 <- runge_kutta(x1 - 0.01, -delta / 2, diff_f, n = 1)
  df     <- (inv_f2 - inv_f3) / (2 * 0.01)
  
  # Strang likelihood
  -sum(stats::dnorm(inv_f, mean = mu_f, sd = sd_f, log = TRUE), na.rm = TRUE) -
    sum(log(abs(df)), na.rm = TRUE)
}


CIR_dynamic_likelihood <- function(par, data, delta, alpha0, mu0, sigma){
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

CIR_transform_dynamic_likelihood <- function(par, data, delta, alpha0, mu0, sigma){
  tau     <-  par[1]
  A       <-  par[2]
  nu      <- if(length(par) == 3) par[3] else 1
  
  N       <- length(data)
  x1    <- data[2:N]
  x0    <- data[1:(N-1)]
  time    <- delta * (1:(N-1))
  
  m          <- mu0 - alpha0 / (2 * A)
  lambda0    <- -alpha0^2 / (4 * A)
  lam_seq    <- lambda0 * (1 - time / tau)^nu
  
  #if(any((-A*(sigma^2 + 4 * lam_seq))<0)){return(50000)}
  
  A_linear_part <- 1 / (2 * m + sqrt(-A*(sigma^2 + 4 * lam_seq))) * 
    (lam_seq * (1 - 4 * A^2) + sigma^2 * (1 / 4 - A^2)) - A * sqrt(-A*(sigma^2 + 4 * lam_seq))
  
  b <- sqrt(4 * m + 2 * sqrt(-A*(sigma^2 + 4 * lam_seq)))
  
  diff_f <- function(t, y){-2 / y * (A * (y^2 / 4 - m)^2 + lam_seq + sigma^2 / 4) - 
      A_linear_part * (y - b)}
  
  # Solution to ODE
  f_h <- runge_kutta(x0, delta / 2, diff_f, n = 1)
  
  mu_f <- exp(A_linear_part * delta) * (f_h - b) + b
  sd_f <- sigma * sqrt(expm1(2 * A_linear_part * delta) / (2 * A_linear_part))
  
  # Inverse of non-linear ODE
  inv_f <- runge_kutta(x1, -delta / 2, diff_f, n = 1)
  
  # Derivative of inverse using Richardson Extrapolation.
  inv_f2 <- runge_kutta(x1 + 0.01, -delta / 2, diff_f, n = 1)
  inv_f3 <- runge_kutta(x1 - 0.01, -delta / 2, diff_f, n = 1)
  df     <- (inv_f2 - inv_f3) / (2 * 0.01)
  # Strang likelihood
  -sum(stats::dnorm(inv_f, mean = mu_f, sd = sd_f, log = TRUE), na.rm = TRUE) -
    sum(log(abs(df)), na.rm = TRUE)
}

CIR_dynamic_simulation_likelihood <- function(par, data, times, M, N, alpha0, mu0, sigma, t_0){
  tau <- par[1]
  A   <- par[2]
  nu  <- if(length(par) == 3) par[3] else 1
  
  m_tip   <- mu0 - alpha0/(2 * A)
  lambda0 <- -alpha0^2 / (4 * A)
  
  delta   <- 1 / M
  lam_seq <- lambda0 * (1 - (times[-length(times)]-t_0) / tau)^nu
  lam_mat <- t(matrix(rep(lam_seq, length.out = N * length(lam_seq)), nrow = length(lam_seq), ncol = N))
  numData <- length(data) - 1
  
  
  dW <- array(dqrng::dqrnorm(N * M * numData, mean = 0, sd = sqrt(delta)), dim = c(M, N, numData))
  # Initialize the process array
  X <- array(NA, dim = c(M, N, numData))
  
  X[1, , ] <- rep(data[-length(data)], each = N)
  for (m in 2:M){
    if(any(X[m - 1, , ] < 0 | any(is.na(sqrt(X[m - 1, , ]))))){ 
      return(50000 + runif(1))
    }
    
    X[m, , ] <- X[m - 1, , ] - (A * (abs(X[m - 1, , ]) - m_tip)^2 + lam_mat) * delta + sigma * sqrt(abs(X[m - 1, , ])) * dW[m - 1, , ]
  }
  
  second_to_last_X <- t(X[M - 1, , ])
  
  sigma_t   <- sigma * sqrt(delta) * sqrt(second_to_last_X)
  
  mu_t      <- data[1:numData] - (A * (second_to_last_X - m_tip)^2 + lam_seq) * delta
  
  densities <- dnorm(data[-1], mean = mu_t, sd = sigma_t, log = TRUE)
  
  loglik    <- -mean(colMeans(densities))
  
  if(is.na(loglik) | is.infinite(loglik)){
    print("NA Likelihood")
    return(1000000)
  }
  
  loglik
}

#-----------------------------------------------------------------------------------------------------------------------------#
# Implementation of methods that are based on the linear noise term process
mean_reverting_GBM_EM_likelihood <- function(par, data, delta){
  beta  <- par[1]
  mu    <- par[2]
  sigma <- exp(par[3])
  
  N <- length(data)
  
  Xupp <- data[2:N]
  Xlow <- data[1:(N - 1)]
  
  sd.part <- sqrt(delta) * sigma * Xlow
  
  mu.part <- Xlow - beta * (Xlow - mu) * delta
  -sum(stats::dnorm(Xupp, mu.part, sd.part, log = TRUE), na.rm = TRUE)
}

mean_reverting_GBM_Kessler_likelihood <- function(par, data, delta){
  beta  <- par[1]
  mu    <- par[2]
  sigma <- exp(par[3])
  
  N <- length(data)
  
  Xupp <- data[2:N]
  Xlow <- data[1:(N - 1)]
  
  
  
  mu_ks  <- exp(-beta * delta) * (Xlow - mu) + mu
  
  var_ks <- ((Xlow-mu)^2  * (1 - exp(-sigma^2 * delta)) -
               2 * mu * sigma^2 / (beta - sigma^2) * (Xlow - mu) * (1 - exp(-(sigma^2 - beta) * delta)) - 
               mu^2 * sigma^2 / (2 * beta - sigma^2) * (1 - exp(-(sigma^2 - 2 * beta) * delta))) * exp(-(2 * beta - sigma^2) * delta)
  
  sd_ks  <- sqrt(var_ks)
  
  -sum(dnorm(Xupp, mean = mu_ks, sd = sd_ks, log = TRUE), na.rm = TRUE) 
}

mean_reverting_GBM_martingale <- function(par, data, delta){
  x0 <- data[1:(length(data) - 1)]
  x1 <- data[2:length(data)]
  
  beta  <- par[1]
  mu    <-  par[2]
  sigma <- par[3]
  
  F_i   <- exp(-beta * delta) * (x0 - mu) + mu

  phi_i <- exp(-(2 * beta - sigma^2) * delta) *
    (x0^2 - (2 * beta * mu)/(beta - sigma^2) * x0 + 2 * beta^2 * mu^2/ ((beta - sigma^2) * (2 * beta - sigma^2))) + 
    2 * beta * mu / (beta - sigma^2) * (exp(-beta * delta) * (x0 - mu) + mu) - 
    2 * beta^2 * mu^2 / ((beta - sigma^2) * (2 * beta - sigma^2)) - 
    exp(-2 * beta * delta) * (x0 - mu)^2 - mu^2 - 2 * mu * exp(-beta * delta) * (x0 - mu)
  
  
  eq1 <- beta / sigma^2 * sum(1 / x0 * (x1 - F_i))
  eq2 <- - 1 / sigma^2 * sum(1 / x0 * (x1 - F_i))
  eq3 <- 1 / sigma^3 * sum(1 / x0^2 * ((x1 - F_i)^2 - phi_i))
  
  c(eq1, eq2, eq3)
}

mean_reverting_GBM_strang <- function(par, data, delta, exp_sigma = FALSE) {
  x0 <- data[1:(length(data) - 1)]
  x1 <- data[2:length(data)]
  
  beta  <- par[1]
  mu    <- par[2]
  if(exp_sigma){sigma <- exp(par[3])} else{sigma <- par[3]}
  
  if((sigma^2 + 2 * beta) / (2 * beta * mu) < 0 | (2 * beta + sigma^2) / (2 * beta * mu) < 0){return(50000)}
  
  diff_f <- function(t, y){-beta + beta * mu * exp(-y) - sigma^2 / 2 + (sigma^2 + 2 * beta) / 2 * y +
      (sigma^2 + 2 * beta) / 2 * log((sigma^2 + 2 * beta) / (2 * beta * mu))}
  
  A <- - (sigma^2 + 2 * beta) / 2
  b <- - log((2 * beta + sigma^2) / (2 * beta * mu))
  
  # Solution to non-linear ODE
  f_h    <- runge_kutta(x0, delta / 2, diff_f, n = 1)
  mu_f   <- exp(A * delta) * (f_h - b) + b
  sd_f   <- sigma * sqrt(expm1(2 * A * delta) / (2 * A))
  
  # Inverse to non-linear ODE
  inv_f  <- runge_kutta(x1, -delta / 2, diff_f, n = 1)
  
  # Derivative of inverse using Richardson Extrapolation
  inv_f2 <- runge_kutta(x1 + 0.01, -delta / 2, diff_f, n = 1)
  inv_f3 <- runge_kutta(x1 - 0.01, -delta / 2, diff_f, n = 1)
  df     <- (inv_f2 - inv_f3) / (2 * 0.01) 
  
  # Strang likelihood
  -sum(stats::dnorm(inv_f, mean = mu_f, sd = sd_f, log = TRUE), na.rm = TRUE) -
    sum(log(abs(df)), na.rm = TRUE)
}

mean_reverting_GBM_dynamic_likelihood <- function(par, data, delta, alpha0, mu0, sigma){
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
      2 * mu_seq * sigma^2 / (alpha_seq - sigma^2) * (fh_half - mu_seq) * (1 - exp(-(sigma^2 - alpha_seq) * delta)) -
      mu_seq^2 * sigma^2 / (2 * alpha_seq - sigma^2) * (1 - exp(-(sigma^2 - 2*alpha_seq) * delta))
  )
  
  mu.part          <- exp(-alpha_seq * delta) * (fh_half - mu_seq) + mu_seq
  
  -sum(stats::dnorm(fh_half_inv, mu.part, sd.part, log = TRUE), na.rm = TRUE) - 
    sum(log(abs(det_Dfh_half_inv)), na.rm = TRUE)
}

mean_reverting_GBM_transform_dynamic_likelihood <- function(par, data, delta, alpha0, mu0, sigma){
  tau     <-  par[1]
  A       <-  par[2]
  nu      <- if(length(par) == 3) par[3] else 1
  
  N       <- length(data)
  x1      <- data[2:N]
  x0      <- data[1:(N-1)]
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
  f_h <- runge_kutta(x0, delta / 2, diff_f, n = 1)
  
  mu_f <- exp(A_linear_part * delta) * (f_h - b) + b
  sd_f <- sigma * sqrt(expm1(2 * A_linear_part * delta) / (2 * A_linear_part))
  
  # Inverse of non-linear ODE
  inv_f <- runge_kutta(x1, -delta / 2, diff_f, n = 1)
  
  # Derivative of inverse using Richardson Extrapolation.
  inv_f2 <- runge_kutta(x1 + 0.01, -delta / 2, diff_f, n = 1)
  inv_f3 <- runge_kutta(x1 - 0.01, -delta / 2, diff_f, n = 1)
  df     <- (inv_f2 - inv_f3) / (2 * 0.01)
  
  # Strang likelihood
 -sum(stats::dnorm(inv_f, mean = mu_f, sd = sd_f, log = TRUE), na.rm = TRUE) - 
    sum(log(abs(df)), na.rm = TRUE)
}

mean_reverting_GBM_dynamic_simulation_likelihood <- function(par, data, times, M, N, alpha0, mu0, sigma, t_0){
  tau <- par[1]
  A  <- par[2]
  nu <- if(length(par) == 3) par[3] else 1
  
  m_tip   <- mu0 - alpha0/(2 * A)
  lambda0 <- -alpha0^2 / (4 * A)
  delta   <- 1 / M

  lam_seq <- lambda0 * (1 - (times[-length(times)]-t_0) / tau)^nu
  lam_mat <- t(matrix(rep(lam_seq, length.out = N * length(lam_seq)), nrow = length(lam_seq), ncol = N))
  # Number of data points
  numData <- length(data) - 1
  
  
  dW <- array(dqrng::dqrnorm(N * M * numData, mean = 0, sd = sqrt(delta)), dim = c(M, N, numData))
  # Initialize the process array
  X <- array(NA, dim = c(M, N, numData))
  
  X[1, , ] <- rep(data[-length(data)], each = N)
  for (m in 2:M){
    if(any(X[m - 1, , ] < 0 | any(is.na(X[m - 1, , ])))){ 
      return(50000 + runif(1))
    }
    X[m, , ] <- X[m - 1, , ] - (A * (X[m - 1, , ] - m_tip)^2 + lam_mat) * delta +
      sigma * X[m - 1, , ] * dW[m - 1, , ]
  }
  
  second_to_last_X <- t(X[M - 1, , ])
  
  sigma_t   <- sigma * sqrt(delta) * second_to_last_X
  
  mu_t      <- data[1:numData] - (A * (second_to_last_X - m_tip)^2 + lam_seq) * delta
  
  densities <- dnorm(data[-1], mean = mu_t, sd = sigma_t, log = TRUE)
  
  # Compute the mean log likelihood for each time point
  loglik <- -mean(colMeans(densities))
  print(loglik)
  if(is.na(loglik) | is.infinite(loglik)){
    print("NA Likelihood")
    return(1000000)
  }
  loglik
}

#-----------------------------------------------------------------------------------------------------------------------------#
# Implementation of methods that are based on the process with stationary t-distribution
t_diffusion_strang_splitting <- function(par, data, delta, exp_sigma = FALSE) {
  x0 <- data[1:(length(data) - 1)]
  x1 <- data[2:length(data)]
  
  beta  <- par[1]
  mu    <- par[2]
  if(exp_sigma){sigma <- exp(par[3])} else{sigma <- par[3]}
  
  asinh_argument <- (2 * mu * beta / (2 * beta + sigma^2))
  
  asinh_term <- asinh(asinh_argument)
  
  diff_f <- function(t, y){(y - tanh(y) - asinh_term) * (beta + 1 / 2 * sigma^2) + mu * beta / cosh(y)}
  
  # Solution to ODE
  f_h <- runge_kutta(x0, delta / 2, diff_f, n = 1)
  
  A <- - (beta + 1 / 2 * sigma^2)
  b <- asinh_term
  
  mu_f <- exp(A * delta) * (f_h - b) + b
  sd_f <- sigma * sqrt(expm1(2 * A * delta) / (2 * A))
  
  # Inverse of non-linear ODE
  inv_f <- runge_kutta(x1, -delta / 2, diff_f, n = 1)
  
  # Derivative of inverse using Richardson Extrapolation.
  inv_f2 <- runge_kutta(x1 + 0.01, -delta / 2, diff_f, n = 1)
  inv_f3 <- runge_kutta(x1 - 0.01, -delta / 2, diff_f, n = 1)
  df     <- (inv_f2 - inv_f3) / (2 * 0.01)
  
  # Strang likelihood
  -sum(stats::dnorm(inv_f, mean = mu_f, sd = sd_f, log = TRUE), na.rm = TRUE) -
    sum(log(abs(df)), na.rm = TRUE)
}

t_diffusion_martingale <- function(par, data, delta){
  x0 <- data[1:(length(data) - 1)]
  x1 <- data[2:length(data)]
  
  beta  <- par[1]
  mu    <-  par[2]
  sigma <- par[3]
  
  F_i   <- exp(-beta * delta) * (x0 - mu) + mu
  
  phi_i <- exp(-2 * beta * delta) * expm1(sigma^2 * delta) * x0^2 +
    (2 * beta * mu / (sigma^2 - beta) * expm1(-(beta - sigma^2) * delta) +
       2 * mu * expm1(-beta * delta) ) * x0 * exp(-beta * delta) + 
    (sigma^2 * (beta - sigma^2) - 2 * beta^2 * mu^2) / ( (2 * beta - sigma^2) * (beta - sigma^2) ) *
    expm1(-(2 * beta - sigma^2) * delta) -
    (2 * beta / (sigma^2 - beta) * expm1(-beta * delta) - 
       1 + exp(-beta * delta) * (1 - expm1(-beta * delta))) * mu^2
  
  
  
  eq1 <-  sum(1 / (x0^2 + 1) * (x1 - F_i))
  eq2 <-  -sum(x0 / (x0^2 + 1) * (x1 - F_i))
  eq3 <-  sum(1 / ((x0^2 + 1)) * ((x1 - F_i)^2 - phi_i))
  
  c(eq1, eq2, eq3)
}

# t_dynamic_likelihood <- function(par, data, delta, alpha0, mu0, sigma){
#   tau     <-  par[1]
#   A       <-  par[2]
#   nu      <- if(length(par) == 3) par[3] else 1
#   
#   N       <- length(data)
#   Xupp    <- data[2:N]
#   Xlow    <- data[1:(N-1)]
#   time    <- delta * (1:(N-1))
#   
#   m          <- mu0 - alpha0 / (2 * A)
#   lambda0    <- -alpha0^2 / (4 * A)
#   lam_seq    <- lambda0 * (1 - time / tau)^nu
#   alpha_seq  <- 2 * sqrt(abs(A * lam_seq))
#   mu_seq     <- m + ifelse(A >= 0, 1, -1) * sqrt(abs(lam_seq / A))
#   
#   # ## Calculating the Strang splitting scheme pseudo likelihood
#   fh_half_tmp_low <-  A * delta * (Xlow - mu_seq) / 2
#   fh_half_tmp_upp <-  A * delta * (Xupp - mu_seq) / 2
#   
#   fh_half     <- (mu_seq * fh_half_tmp_low + Xlow) / (fh_half_tmp_low+1)
#   
#   fh_half_inv <- (mu_seq * fh_half_tmp_upp - Xupp) / (fh_half_tmp_upp-1)
#   
#   det_Dfh_half_inv <- 1 / (fh_half_tmp_upp-1)^2
#   
#   phi_dot <- exp(-2 * alpha_seq * delta) * expm1(sigma^2 * delta) * fh_half^2 +
#     (2 * alpha_seq * mu_seq / (sigma^2 - alpha_seq) * expm1(-(alpha_seq - sigma^2) * delta) +
#        2 * mu_seq * expm1(-alpha_seq * delta) ) * fh_half * exp(-alpha_seq * delta) + 
#     (sigma^2 * (alpha_seq - sigma^2) - 2 * alpha_seq^2 * mu_seq^2) / ( (2 * alpha_seq - sigma^2) * (alpha_seq - sigma^2) ) *
#     expm1(-(2 * alpha_seq - sigma^2) * delta) -
#     (2 * alpha_seq / (sigma^2 - alpha_seq) * expm1(-alpha_seq * delta) - 
#        1 + exp(-alpha_seq * delta) * (1 - expm1(-alpha_seq * delta))) * mu_seq^2
#   
#   if(any(phi_dot < 0) | any(is.na(phi_dot))){return(100000 + runif(1, min = -5, max = 1))}
#   
#   sd.part  <- sigma * sqrt(phi_dot)
#   
#   mu.part  <- exp(-alpha_seq * delta) * (fh_half - mu_seq) + mu_seq
#   
#   -sum(stats::dnorm(fh_half_inv, mu.part, sd.part, log = TRUE)) - sum(log(abs(det_Dfh_half_inv)))
# }

t_transform_dynamic_likelihood <- function(par, data, delta, alpha0, mu0, sigma){
  tau     <-  par[1]
  A       <-  par[2]
  nu      <- if(length(par) == 3) par[3] else 1
  
  N       <- length(data)
  x1    <- data[2:N]
  x0    <- data[1:(N-1)]
  time    <- delta * (1:(N-1))
  
  m          <- mu0 - alpha0 / (2 * A)
  lambda0    <- -alpha0^2 / (4 * A)
  lam_seq    <- lambda0 * (1 - time / tau)^nu
  
  if(any((sigma^4 - 8 * A * (2 * lam_seq + m * sigma^2)) < 0) |
     any(is.na(sigma^4 - 8 * A * (2 * lam_seq + m * sigma^2)))){
    par <- par + runif(2, min = -1) * c(10, 0.2)
    return(50000 + runif(1, min = -5, max = 5))
    }
  
  fix_points <- asinh((sqrt(sigma^4 - 8 * A * (2 * lam_seq + m * sigma^2)) + 4 * A * m - sigma^2) /
                          (4 * A))
  
  A_linear_part <- - sigma^2 + sqrt(sigma^4 / 4 - 2 * A * (2 * lam_seq + m * sigma^2))
  b <- fix_points
  
  diff_f <- function(t, y){-1 / cosh(y) * (A * (sinh(y) - m)^2 + lam_seq + sigma^2 / 2 * sinh(y)) -
      A_linear_part * (y - b)}
  
  # Solution to ODE
  f_h <- runge_kutta(x0, delta / 2, diff_f, n = 1)
  
  mu_f <- exp(A_linear_part * delta) * (f_h - b) + b
  sd_f <- sigma * sqrt(expm1(2 * A_linear_part * delta) / (2 * A_linear_part))
  
  # Inverse of non-linear ODE
  inv_f <- runge_kutta(x1, -delta / 2, diff_f, n = 1)
  
  # Derivative of inverse using Richardson Extrapolation.
  inv_f2 <- runge_kutta(x1 + 0.01, -delta / 2, diff_f, n = 1)
  inv_f3 <- runge_kutta(x1 - 0.01, -delta / 2, diff_f, n = 1)
  df     <- (inv_f2 - inv_f3) / (2 * 0.01)
  
  # Strang likelihood
  -sum(stats::dnorm(inv_f, mean = mu_f, sd = sd_f, log = TRUE), na.rm = TRUE) -
    sum(log(abs(df)), na.rm = TRUE)
}

#-----------------------------------------------------------------------------------------------------------------------------#
# Implementation of methods that are based on the process with stationary F-distribution

F_diffusion_strang_splitting <- function(par, data, delta, exp_sigma = FALSE) {
  x0 <- data[1:(length(data) - 1)]
  x1 <- data[2:length(data)]
  
  beta  <- par[1]
  mu    <- par[2]
  if(exp_sigma){sigma <- exp(par[3])} else{sigma <- par[3]}
  
  A <- - (beta + sigma^2 / 2)
  b <- acosh(2 * beta * (2 * mu + 1) / (2 * beta + sigma^2))
  
  diff_f <- function(t,y){-A * (y - b - cosh(y) / sinh(y)) + beta * (2 * mu + 1) / sinh(y)}
  # Solution to ODE
  f_h <- runge_kutta(x0, delta / 2, diff_f, n = 1)
  
  mu_f <- exp(A * delta) * (f_h - b) + b
  sd_f <- sigma * sqrt(expm1(2 * A * delta) / (2 * A))
  
  # Inverse of non-linear ODE
  inv_f <- runge_kutta(x1, -delta / 2, diff_f, n = 1)
  
  # Derivative of inverse using Richardson Extrapolation.
  inv_f2 <- runge_kutta(x1 + 0.01, -delta / 2, diff_f, n = 1)
  inv_f3 <- runge_kutta(x1 - 0.01, -delta / 2, diff_f, n = 1)
  df     <- (inv_f2 - inv_f3) / (2 * 0.01)

  # Strang likelihood
  -sum(stats::dnorm(inv_f, mean = mu_f, sd = sd_f, log = TRUE), na.rm = TRUE) -
    sum(log(abs(df)), na.rm = TRUE)
}

F_transform_dynamic_likelihood <- function(par, data, delta, alpha0, mu0, sigma){
  tau     <-  par[1]
  A       <-  par[2]
  nu      <- if(length(par) == 3) par[3] else 1
  
  N       <- length(data)
  x1    <- data[2:N]
  x0    <- data[1:(N-1)]
  time    <- delta * (1:(N-1))
  
  m          <- mu0 - alpha0 / (2 * A)
  lambda0    <- -alpha0^2 / (4 * A)
  lam_seq    <- lambda0 * (1 - time / tau)^nu
  
  sqrt_arg <-  pmax(sigma^4 / 4 + m^2 - sigma^2 * m - A * (sigma^2 + 4 * lam_seq + 4 * A * m^2), 0.01)
  
  zeta <- (m + A - sigma^2 / 2) - sqrt(sqrt_arg)
  fix_points <- acosh(max(zeta / A, 1.01))
  
  A_linear_part <- - zeta + (m + A - sigma^2 / 2) 
  b <- fix_points
  diff_f <- function(t, y){-1 / sinh(y) * (A / 2 * cosh(y)^2 + (sigma^2 / 2 - m - A) * cosh(y) + 
                                             2 * lam_seq + m + (2 * m^2 + 1 / 2) * A) -  
      A_linear_part * (y - b)}
  
  # Solution to ODE
  f_h <- runge_kutta(x0, delta / 2, diff_f, n = 1)
  
  mu_f <- exp(A_linear_part * delta) * (f_h - b) + b
  sd_f <- sigma * sqrt(expm1(2 * A_linear_part * delta) / (2 * A_linear_part))
  
  # Inverse of non-linear ODE
  inv_f <- runge_kutta(x1, -delta / 2, diff_f, n = 1)
  
  # Derivative of inverse using Richardson Extrapolation.
  inv_f2 <- runge_kutta(x1 + 0.01, -delta / 2, diff_f, n = 1)
  inv_f3 <- runge_kutta(x1 - 0.01, -delta / 2, diff_f, n = 1)
  df     <- (inv_f2 - inv_f3) / (2 * 0.01)
  
  # Strang likelihood
  -sum(stats::dnorm(inv_f, mean = mu_f, sd = sd_f, log = TRUE), na.rm = TRUE) -
    sum(log(abs(df)), na.rm = TRUE)
}

#-----------------------------------------------------------------------------------------------------------------------------#
# Implementation of methods that are based on the jacobi diffusion

jacobi_diffusion_strang_splitting <- function(par, data, delta, exp_sigma = FALSE) {
  x0 <- data[1:(length(data) - 1)]
  x1 <- data[2:length(data)]
  
  beta  <- par[1]
  mu    <- par[2]
  if(exp_sigma){sigma <- exp(par[3])} else{sigma <- par[3]}
  
  if(abs(beta * (2 * mu - 1) / (sigma^2 / 2 - beta))>1){return(50000)}
  
  A <- 1/2 * sigma^2 - beta
  
  b <- acos(beta * (2 * mu - 1) / (sigma^2 / 2 - beta))
  
  diff_f <- function(t, y){-1 / sin(y) * (beta * (1 - cos(y) - 2 * mu) + sigma^2 / 2 * cos(y)) - A * (y - b)}
  
  # Solution to ODE
  f_h <- runge_kutta(x0, delta / 2, diff_f, n = 1)
  
  if(any(is.na(f_h))){return(50000)}
  
  mu_f <- exp(A * delta) * (f_h - b) + b
  sd_f <- sigma * sqrt(expm1(2 * A * delta) / (2 * A))
  
  # Inverse of non-linear ODE
  inv_f <- runge_kutta(x1, -delta / 2, diff_f, n = 1)
  
  # Derivative of inverse using Richardson Extrapolation.
  inv_f2 <- runge_kutta(x1 + 0.01, -delta / 2, diff_f, n = 1)
  inv_f3 <- runge_kutta(x1 - 0.01, -delta / 2, diff_f, n = 1)
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
  x1    <- data[2:N]
  x0    <- data[1:(N-1)]
  time    <- delta * (1:(N-1))
  
  m          <- mu0 - alpha0 / (2 * A)
  lambda0    <- -alpha0^2 / (4 * A)
  lam_seq    <- lambda0 * (1 - time / tau)^nu
  
  sqrt_arg <- sigma^4 / 4 - sigma^2 * A - 2 * sigma^2 * A * m - 2 * A * lam_seq
  zeta <- (A - sigma^2 / 2 - 2 * A) + sqrt(max(sqrt_arg, 0.01))
  

  if(any(sqrt_arg < 0) || any(is.na(sqrt_arg))){return(50000)}
  
  fix_points <- acos(zeta / A)
  
  A_linear_part <- zeta + sigma^2 / 2 - 2 * A * m - A
  b <- fix_points
  
  diff_f <- function(t, y){-1 / sin(y) * 
      (A / 2 * cos(y)^2 + (sigma^2 / 2 + 2 * A * m - A) * cos(y) +
       A * (1 / 2 + 2 * m^2 - 2 * m) + 2 * lam_seq) -
      A_linear_part * (y - b)}
  
  # Solution to ODE
  f_h <- runge_kutta(x0, delta / 2, diff_f, n = 1)
  
  mu_f <- exp(A_linear_part * delta) * (f_h - b) + b
  sd_f <- sigma * sqrt(expm1(2 * A_linear_part * delta) / (2 * A_linear_part))
  
  # Inverse of non-linear ODE
  inv_f <- runge_kutta(x1, -delta / 2, diff_f, n = 1)
  
  # Derivative of inverse using Richardson Extrapolation.
  inv_f2 <- runge_kutta(x1 + 0.01, -delta / 2, diff_f, n = 1)
  inv_f3 <- runge_kutta(x1 - 0.01, -delta / 2, diff_f, n = 1)
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
                                           ...){
  res_optim <- stats::optim(init_par, fn = likelihood_fun,
                     method = method,
                     data = data, delta = delta, exp_sigma = exp_sigma, ...)
  if(exp_sigma){
  res_optim$par[3] <- exp(res_optim$par[3])
  }
  res_optim$par
}

optimize_dynamic_likelihood <- function(likelihood_fun, data,
                                        init_par, delta,
                                        alpha0, mu0, sigma, 
                                        method = "BFGS", ...){

  res_optim <- stats::optim(init_par, fn = likelihood_fun,
                            method = method,
                            data = data, delta = delta,
                            alpha0 = alpha0 , mu0 = mu0, sigma = sigma, ...)
  list(par = res_optim$par, objective = res_optim$value)
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
# true_param <- c(0.87, -1.51, -2.69, 0.2)
# actual_dt <- 0.005
# tau <- 150
# t_0 <- 50
# sim_res_add <- simulate_additive_noise_tipping_model(actual_dt, true_param, tau, t_0)
# sample_n(sim_res_add, min(nrow(sim_res_add), 10000)) |> ggplot(aes(x = t, y = X_t)) + geom_step()
# # Stationary part
# # Parameters for stationary part
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
#
## Dynamic part
# dynamic_part_true_param <- c(tau, true_param[1])
# 
# 
# optimize_dynamic_likelihood(likelihood_fun = OU_dynamic_likelihood,
#                             data = sim_res_add$X_t[sim_res_add$t > t_0],
#                             init_par = dynamic_part_true_param,
#                             delta = actual_dt,
#                             alpha0 = stationary_part_true_param[1],
#                             mu0 = stationary_part_true_param[2],
#                             sigma = stationary_part_true_param[3])

# nleqslv::nleqslv(OU_dynamic_score,
#                  x = dynamic_part_true_param,
#                  data = sim_res_add$X_t[sim_res_add$t > t_0],
#                  delta = actual_dt,
#                  alpha0 = stationary_part_true_param[1],
#                  mu0 = stationary_part_true_param[2],
#                  sigma = stationary_part_true_param[3])$x

#-----------------------------------------------------------------------------------------------------------------------------#

# # Square-root noise model
# true_param <- c(1, 2, -1, 0.05)
# actual_dt <- 0.005
# tau <- 150
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
# optimize_stationary_likelihood(CIR_alt_strang_splitting,
#                                exp(2*sqrt(sim_res_sqrt$X_t[sim_res_sqrt$t < t_0])),
#                                stationary_part_true_param, actual_dt,
#                                exp_sigma = FALSE)
# 
# optimize_stationary_likelihood(likelihood_fun = CIR_strang_splitting,
#                                data = 2 * sqrt(sim_res_sqrt$X_t[sim_res_sqrt$t < t_0]),
#                                init_par = stationary_part_true_param,
#                                delta = actual_dt, exp_sigma = TRUE)
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
#                             init_par =dynamic_part_true_param,
#                             delta = actual_dt,
#                             alpha0 = stationary_part_true_param[1],
#                             mu0 = stationary_part_true_param[2],
#                             sigma = stationary_part_true_param[3],
#                             control = list(reltol = sqrt(.Machine$double.eps) / 1000))
#-----------------------------------------------------------------------------------------------------------------------------#

## Linear noise model
# true_param <- c(-0.05, 100, 2, 0.01)
# actual_dt <- 0.001
# tau <- 150
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
# true_param <- c(0.1, -2, -3, 0.1)
# actual_dt <- 0.005
# tau <- 100
# t_0 <- 50
# sim_res_t_distribution <- simulate_t_distribution_tipping_model(actual_dt, true_param, tau, t_0)
# sim_res_t_distribution |> ggplot2::ggplot(ggplot2::aes(x = t, y = X_t)) +
#   ggplot2::geom_step() + ggplot2::geom_hline(yintercept = true_param[2], linetype = "dashed") +
#   ggplot2::geom_vline(xintercept = t_0)
# 
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
# optimize_dynamic_likelihood(likelihood_fun = t_dynamic_likelihood,
#                             data = sim_res_t_distribution$X_t[sim_res_t_distribution$t > t_0],
#                             init_par = dynamic_part_true_param,
#                             delta = actual_dt,
#                             alpha0 = stationary_part_true_param[1],
#                             mu0 = stationary_part_true_param[2],
#                             sigma = stationary_part_true_param[3],
#                             method = "BFGS")
# 
optimize_dynamic_likelihood(likelihood_fun = t_transform_dynamic_likelihood,
                  data = asinh(sim_res_t_distribution$X_t[sim_res_t_distribution$t > t_0]),
                  init_par = dynamic_part_true_param,
                  delta = actual_dt,
                  alpha0 = stationary_part_true_param[1],
                  mu0 = stationary_part_true_param[2],
                  sigma = stationary_part_true_param[3],
                  control = list(reltol = sqrt(.Machine$double.eps) / 1000))

#-----------------------------------------------------------------------------------------------------------------------------#

## F-distributed stationary distribution model

# actual_dt <- 0.01
# t_0 <- 50
# tau <- 50
# true_param <- c(0.3, 2, -4, 0.15)
# 
# F_sim_dynamic <- simulate_F_distribution_tipping_model(actual_dt, true_param, t_0 = t_0, tau = tau)
# 
# sample_n(F_sim_dynamic, 10000) |> ggplot(aes(x = t, y = X_t)) + geom_step()
# 
## Stationary part
## Parameters for stationary part
# mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) * sqrt(abs(true_param[3] / true_param[1]))
# alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
# 
# stationary_part_true_param <- c(alpha0, mu0, true_param[4])

# optimize_stationary_likelihood(
#   likelihood_fun = F_diffusion_strang_splitting,
#   data = 2 * asinh(sqrt(F_sim_dynamic$X_t[F_sim_dynamic$t<t_0])),
#   init_par = stationary_part_true_param,
#   delta = actual_dt,
#   exp_sigma = FALSE,
#   control = list(reltol = sqrt(.Machine$double.eps) / 1000))

# ## Dynamic part
# dynamic_part_true_param <- c(tau, true_param[1])
# optimize_dynamic_likelihood(likelihood_fun = F_transform_dynamic_likelihood,
#                             data = 2 * asinh(sqrt(F_sim_dynamic$X_t[F_sim_dynamic$t>t_0])),
#                             init_par = dynamic_part_true_param,
#                             delta = actual_dt,
#                             alpha0 = stationary_part_true_param[1],
#                             mu0 = stationary_part_true_param[2],
#                             sigma = stationary_part_true_param[3],
#                             control = list(reltol = sqrt(.Machine$double.eps) / 1000))

#-----------------------------------------------------------------------------------------------------------------------------#

## Jacobi diffusion 
# true_param <- c(5, 0.1, -0.8, 0.15)
# actual_dt <- 0.005
# tau <- 100
# t_0 <- 50
# 
# sim_res_jacobi <- simulate_jacobi_diffusion_tipping_model(actual_dt, true_param, tau, t_0,
#                                                           beyond_tipping = - 35)
# sim_res_jacobi |> ggplot2::ggplot(ggplot2::aes(x = t, y = X_t)) +
#   ggplot2::geom_step() + ggplot2::geom_hline(yintercept = true_param[2], linetype = "dashed") +
#   ggplot2::geom_vline(xintercept = t_0)

## Stationary part
## Parameters for stationary part
# mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) * sqrt(abs(true_param[3] / true_param[1]))
# alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
# stationary_part_true_param <- c(alpha0, mu0, true_param[4])
# 
# optimize_stationary_likelihood(
#   likelihood_fun = jacobi_diffusion_strang_splitting,
#   data = 2 * asin(sqrt(sim_res_jacobi$X_t[sim_res_jacobi$t<t_0])),
#   init_par = stationary_part_true_param,
#   delta = actual_dt,
#   exp_sigma = TRUE)
# # 
# dynamic_part_true_param <- c(tau, true_param[1])
# optimize_dynamic_likelihood(likelihood_fun = jacobi_diffusion_transform_dynamic_likelihood,
#                             data = 2 * asin(sqrt(sim_res_jacobi$X_t[sim_res_jacobi$t>t_0])),
#                             init_par = dynamic_part_true_param,
#                             delta = actual_dt,
#                             alpha0 = stationary_part_true_param[1],
#                             mu0 = stationary_part_true_param[2],
#                             sigma = stationary_part_true_param[3],
#                             control = list(reltol = sqrt(.Machine$double.eps) / 1000))

