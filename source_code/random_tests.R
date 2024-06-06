# Title: Random test and
# Author: Anders Gantzhorn Kristensen (University of Copenhagen, andersgantzhorn@gmail.com)
# Date: 2024-04-20 (Last Updated: 2024-05-28)
#-----------------------------------------------------------------------------------------------------------------------------#
# Project: Tipping Point Estimation in Ecological Systems using Stochastic Differential Equations
# Description: This script holds legacy code that did not end up in the final thesis and various
# test conducted on the methods.
#-----------------------------------------------------------------------------------------------------------------------------#
# License: MIT License (for more information, see LICENSE file in the repository).
# Dependencies: Requirements from sourced files.
#-----------------------------------------------------------------------------------------------------------------------------#
source("source_code/tipping_simulations.R")
source("source_code/model_fitting.R")
source("source_code/model_diagnostics.R")
#-----------------------------------------------------------------------------------------------------------------------------#
# Tests
#-----------------------------------------------------------------------------------------------------------------------------#
CIR_lamperti_drift <- function(y, par){
  beta  <- par[1]
  mu    <- par[2]
  sigma <- par[3]
  -2 / y * (beta * (y^2 / 4 - mu) + sigma^2 / 2)
}

true_param <- c(1, 2, -1, 0.05)
actual_dt <- 0.005
tau <- 100
t_0 <- 50
sim_res_sqrt <- simulate_squareroot_noise_tipping_model(actual_dt, true_param, tau, t_0)
sample_n(sim_res_sqrt, min(nrow(sim_res_sqrt), 10000)) |> ggplot(aes(x = t, y = X_t)) + geom_step()
# #
# ## Stationary part
# ## Parameters for stationary part
mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) * sqrt(abs(true_param[3] / true_param[1]))
alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
stationary_part_true_param <- c(alpha0, mu0, true_param[4])

CIR_quadratic_martingale(sim_res_sqrt$X_t[sim_res_sqrt$t < t_0], actual_dt)# - stationary_part_true_param
optimize_stationary_likelihood(CIR_alt_strang_splitting, exp(2*sqrt(sim_res_sqrt$X_t[sim_res_sqrt$t < t_0])),
                               stationary_part_true_param, actual_dt)# - stationary_part_true_param

optimize_stationary_likelihood(likelihood_fun = CIR_strang_splitting,
                               data = 2 * sqrt(sim_res_sqrt$X_t[sim_res_sqrt$t < t_0]),
                               init_par = stationary_part_true_param,
                               delta = actual_dt, exp_sigma = FALSE)
                     
optimize_stationary_likelihood(likelihood_fun = numerical_strang_splitting,
                               data = 2 * sqrt(sim_res_sqrt$X_t[sim_res_sqrt$t < t_0]),
                               init_par = stationary_part_true_param,
                               delta = actual_dt,
                               exp_sigma = FALSE,
                               drift_lamperti_sde = CIR_lamperti_drift)

true_param <- c(-0.05, 100, 2, 0.01, 0.65)
actual_dt <- 0.001
tau <- 150
t_0 <- 50
sim_res_linear <- simulate_linear_noise_tipping_model(actual_dt, true_param, tau, t_0)
sample_n(sim_res_linear, min(nrow(sim_res_linear), 10000)) |> ggplot(aes(x = t, y = X_t)) + geom_step()

# ## Stationary part
## Parameters for stationary part
mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) * sqrt(abs(true_param[3] / true_param[1]))
alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
stationary_part_true_param <- c(alpha0, mu0, true_param[4])

nleqslv::nleqslv(x = stationary_part_true_param, fn = mean_reverting_GMB_martingale,
                 data = sim_res_linear$X_t[sim_res_linear$t < t_0],
                 delta = actual_dt)$x

optimize_stationary_likelihood(mean_reverting_GMB_strang, log(sim_res_linear$X_t[sim_res_linear$t<t_0]),
                               init_par = stationary_part_true_param, delta = actual_dt,
                               exp_sigma = TRUE)

GBM_lamperti_drift <- function(y, par){
  beta  <- par[1]
  mu    <- par[2]
  sigma <- par[3]
  
  -beta * (1 - mu * exp(-y) + sigma^2 / 2)
}

optimize_stationary_likelihood(numerical_strang_splitting,
                               data = log(sim_res_linear$X_t[sim_res_linear$t<t_0]),
                               init_par = stationary_part_true_param, delta = actual_dt,
                               exp_sigma = FALSE,
                               drift_lamperti_sde = GBM_lamperti_drift)


t_lamperti_drift <- function(y, par){
  beta  <- par[1]
  mu    <- par[2]
  sigma <- par[3]
  
  - 1 / cosh(y) * (sinh(y) * (beta + 1 / 2 * sigma^2) - beta * mu)
}

true_param <- c(0.1, -2, -3, 0.1)
actual_dt <- 0.005
tau <- 100
t_0 <- 50
sim_res_t_distribution <- simulate_t_distribution_tipping_model(actual_dt, true_param, tau, t_0)
sim_res_t_distribution |> ggplot2::ggplot(ggplot2::aes(x = t, y = X_t)) +
  ggplot2::geom_step() + ggplot2::geom_hline(yintercept = true_param[2], linetype = "dashed") +
  ggplot2::geom_vline(xintercept = t_0)

# Stationary part
# Parameters for stationary part
mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) * sqrt(abs(true_param[3] / true_param[1]))
alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
stationary_part_true_param <- c(alpha0, mu0, true_param[4])

optimize_stationary_likelihood(
  likelihood_fun = t_diffusion_strang_splitting,
  data = asinh(sim_res_t_distribution$X_t[sim_res_t_distribution$t < t_0]),
  init_par = stationary_part_true_param,
  delta = actual_dt,
  exp_sigma = TRUE)

optimize_stationary_likelihood(
  likelihood_fun = numerical_strang_splitting,
  data = asinh(sim_res_t_distribution$X_t[sim_res_t_distribution$t < t_0]),
  init_par = stationary_part_true_param,
  delta = actual_dt,
  exp_sigma = FALSE,
  drift_lamperti_sde = t_lamperti_drift)


## F-distributed stationary distribution model
F_lamperti_drift <- function(y, par){
  beta  <- par[1]
  mu    <- par[2]
  sigma <- par[3]
  
  - 1 / sinh(y) * ((beta + sigma^2 / 2) * cosh(y) - beta * (2 * mu + 1))
}

actual_dt <- 0.0005
t_0 <- 50
tau <- 100
true_param <- c(0.3, 2, -4, 0.15)

F_sim_dynamic <- simulate_F_distribution_tipping_model(actual_dt, true_param, t_0 = t_0, tau = tau)

sample_n(F_sim_dynamic, 10000) |> ggplot(aes(x = t, y = X_t)) + geom_step()
# 
# ## Stationary part
# # Parameters for stationary part
mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) * sqrt(abs(true_param[3] / true_param[1]))
alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))

stationary_part_true_param <- c(alpha0, mu0, true_param[4])

optimize_stationary_likelihood(
  likelihood_fun = F_diffusion_strang_splitting,
  data = 2 * asinh(sqrt(F_sim_dynamic$X_t[F_sim_dynamic$t<t_0])),
  init_par = stationary_part_true_param,
  delta = actual_dt,
  exp_sigma = FALSE)

optimize_stationary_likelihood(
  likelihood_fun = numerical_strang_splitting,
  data = 2 * asinh(sqrt(F_sim_dynamic$X_t[F_sim_dynamic$t<t_0])),
  init_par = stationary_part_true_param,
  delta = actual_dt,
  exp_sigma = FALSE,
  drift_lamperti_sde = F_lamperti_drift)

numerical_strang_splitting(par = stationary_part_true_param,
                         data = 2 * asinh(sqrt(F_sim_dynamic$X_t[F_sim_dynamic$t<t_0])),
                         delta = actual_dt,
                         drift_lamperti_sde = F_lamperti_drift)


dummy_func <- function(y, par){
  beta  <- par[1]
mu    <- par[2]
sigma <- par[3]
(beta + sigma^2 / 2) * cosh(y) - beta * (2 * mu + 1)
}

# Jacobi
true_param <- c(5, 0.1, -0.8, 0.05)
actual_dt <- 0.1
tau <- 100
t_0 <- 50

sim_res_jacobi <- simulate_jacobi_diffusion_tipping_model(actual_dt, true_param, tau, t_0)
sim_res_jacobi |> ggplot2::ggplot(ggplot2::aes(x = t, y = X_t)) +
  ggplot2::geom_step() + ggplot2::geom_hline(yintercept = true_param[2], linetype = "dashed") +
  ggplot2::geom_vline(xintercept = t_0)
# 
# # Stationary part
# # Parameters for stationary part
mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) * sqrt(abs(true_param[3] / true_param[1]))
alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
stationary_part_true_param <- c(alpha0, mu0, true_param[4])



jacobi_lamperti_drift <- function(y, par){
  beta <- par[1]
  mu <- par[2]
  sigma <- par[3]
  
  - 1 / sin(y) * (sigma^2 / 2 - beta) * cos(y) + beta - 2 * beta * mu
}

optimize_stationary_likelihood(
  likelihood_fun = jacobi_diffusion_strang_splitting,
  data = 2 * asin(sqrt(sim_res_jacobi$X_t[sim_res_jacobi$t<t_0])),
  init_par = stationary_part_true_param,
  delta = actual_dt,
  exp_sigma = TRUE)

optimize_stationary_likelihood(
  likelihood_fun = numerical_strang_splitting,
  data = 2 * asin(sqrt(sim_res_jacobi$X_t[sim_res_jacobi$t<t_0])),
  init_par = stationary_part_true_param,
  delta = actual_dt,
  exp_sigma = FALSE,
  drift_lamperti_sde = jacobi_lamperti_drift)

pure_numerical_strang_splitting <- function(par,
                                     data,
                                     delta,
                                     diffusion_term,
                                     xi,
                                     exp_sigma = FALSE) {
  x0 <- data[1:(length(data) - 1)]
  x1 <- data[2:length(data)]
  
  beta  <- par[1]
  mu    <- par[2]
  if(exp_sigma){sigma <- exp(par[3])} else{sigma <- par[3]}
  
  lamperti_transform_derivative <- function(y){
    1 / diffusion_term(y)
  }
  
  # lamperti_transform <- function(y){
  #   sapply(y, function(x) integrate(f = lamperti_transform_derivative, lower = xi, upper = x,
  #                                   rel.tol = .Machine$double.eps^0.25)$value )
  # }
  
  lamperti_transform_new <- function(y) {
    # Define the number of subdivisions
    n <- 20
    
    h <- (y - xi) / n
    
    yj <- seq.int(xi, y, length.out = n + 1)
    yj <- yj[-1]
    yj <- yj[-length(yj)]
    (h / 3) * (lamperti_transform_derivative(xi) + 
                           2 * sum(lamperti_transform_derivative(yj[seq.int(2, length(yj), 2)])) +
                           4 * sum(lamperti_transform_derivative(yj[seq.int(1, length(yj), 2)])) +
                           lamperti_transform_derivative(y))
  }
  
  inverse_equation <- function(z, y){lamperti_transform_new(y) - z}

  inverse_lamperti_transform <- function(z) {
    rootSolve::multiroot(f = inverse_equation, start = (z - xi)/2 , z = z)$root
    #uniroot(inverse_equation, interval = interval, z = z)$root
  }
  
  
  lamperti_transform_second_derivative <- function(y){
  
    numDeriv::grad(func = lamperti_transform_derivative, x = y)
  }
  
  drift_lamperti_sde <- function(y){
    (- par[1] * (inverse_lamperti_transform(y) - par[2])) *
      lamperti_transform_derivative(inverse_lamperti_transform(y)) +
      par[3]^2 / 2 * lamperti_transform_second_derivative(inverse_lamperti_transform(y))
  }

  b <- nleqslv::nleqslv(x = drift_lamperti_sde(), fn = drift_lamperti_sde)$x
  
  A <- numDeriv::grad(func = drift_lamperti_sde, x = b)

  diff_f <- function(t, y){drift_lamperti_sde(y) - A * (y - b)}
  
  # Solution to ODE
  f_h <- runge_kutta(x0, delta / 2, diff_f, n = 1)
  
  
  mu_f <- exp(A * delta) * (f_h - b) + b
  sd_f <- sigma * sqrt(expm1(2 * A * delta) / (2 * A))
  browser()
  # Inverse of non-linear ODE
  inv_f <- runge_kutta(x1, -delta / 2, diff_f, n = 1)
  
  # Derivative of inverse using Richardson Extrapolation.
  inv_f2 <- runge_kutta(x1 + 0.01, -delta / 2, diff_f, n = 1)
  inv_f3 <- runge_kutta(x1 - 0.01, -delta / 2, diff_f, n = 1)
  df     <- (inv_f2 - inv_f3) / (2 * 0.01)
  
  # Strang likelihood
  loglike <- -mean(stats::dnorm(inv_f, mean = mu_f, sd = sd_f, log = TRUE), na.rm = TRUE) - 
    mean(log(abs(df)), na.rm = TRUE)
  print(loglike)
  print(par)
  loglike
}

jacobi_diffusion_term <- function(y){
  sqrt(y * (1 - y))
}

optimize_stationary_likelihood(
  likelihood_fun = pure_numerical_strang_splitting,
  data = sim_res_jacobi$X_t[sim_res_jacobi$t<t_0],
  init_par = stationary_part_true_param,
  delta = actual_dt,
  exp_sigma = FALSE,
  diffusion_term = jacobi_diffusion_term,
  xi = 0.01)

#-----------------------------------------------------------------------------------------------------------------------------#
# Legacy code
#-----------------------------------------------------------------------------------------------------------------------------#

# Old Likelihood implementations

# Simulated likelihood methods
OU_dynamic_simulation_likelihood <- function(par, data, times,
                                             M, N, alpha0, mu0,
                                             sigma, t_0){
  tau <- par[1]
  A   <- par[2]
  nu  <- if(length(par) == 3) par[3] else 1
  
  m_tip   <- mu0 - alpha0/(2 * A)
  lambda0 <- -alpha0^2 / (4 * A)
  delta   <- 1 / M
  
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

# Other likelihoods

CIR_alt_strang_splitting <- function(par, data, delta, exp_sigma = FALSE) {
  Xlow <- data[1:(length(data) - 1)]
  Xupp <- data[2:length(data)]
  
  beta  <- par[1]
  mu    <- par[2]
  if(exp_sigma){sigma <- exp(par[3])} else{sigma <- par[3]}
  
  # Discourage method to pick undefined values
  if(4 * beta * mu - sigma^2 < 0){return(100000)}
  
  diff_f <- function(t, y){y/2 * ((4 * beta * mu - sigma^2)/log(y) - beta * log(y))}
  
  f <- runge_kutta(Xlow, delta / 2, diff_f)
  
  mu_log <- log(f)
  
  sd_log <- sigma * sqrt(delta)
  
  inv_f  <- runge_kutta(Xupp, -delta / 2, diff_f)
  
  inv_f2 <- runge_kutta(Xupp + 0.01, -delta / 2, diff_f)
  inv_f3 <- runge_kutta(Xupp - 0.01, -delta / 2, diff_f)
  
  df <- (inv_f2 - inv_f3) / (2 * 0.01)
  
  -sum(stats::dlnorm(inv_f, meanlog = mu_log, sdlog = sd_log, log = TRUE), na.rm = TRUE) -
    sum(log(abs(df)), na.rm = TRUE)
}

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

t_diffusion_martingale <- function(par, data, delta){
  Xlow <- data[1:(length(data) - 1)]
  Xupp <- data[2:length(data)]
  
  beta  <- par[1]
  mu    <-  par[2]
  sigma <- par[3]
  
  F_i   <- exp(-beta * delta) * (Xlow - mu) + mu
  
  phi_i <- exp(-2 * beta * delta) * expm1(sigma^2 * delta) * Xlow^2 +
    (2 * beta * mu / (sigma^2 - beta) * expm1(-(beta - sigma^2) * delta) +
       2 * mu * expm1(-beta * delta) ) * Xlow * exp(-beta * delta) + 
    (sigma^2 * (beta - sigma^2) - 2 * beta^2 * mu^2) / ( (2 * beta - sigma^2) * (beta - sigma^2) ) *
    expm1(-(2 * beta - sigma^2) * delta) -
    (2 * beta / (sigma^2 - beta) * expm1(-beta * delta) - 
       1 + exp(-beta * delta) * (1 - expm1(-beta * delta))) * mu^2
  
  
  
  eq1 <-  sum(1 / (Xlow^2 + 1) * (Xupp - F_i))
  eq2 <-  -sum(Xlow / (Xlow^2 + 1) * (Xupp - F_i))
  eq3 <-  sum(1 / ((Xlow^2 + 1)) * ((Xupp - F_i)^2 - phi_i))
  
  c(eq1, eq2, eq3)
}

t_dynamic_likelihood <- function(par, data, delta, alpha0, mu0, sigma){
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
  
  phi_dot <- exp(-2 * alpha_seq * delta) * expm1(sigma^2 * delta) * fh_half^2 +
    (2 * alpha_seq * mu_seq / (sigma^2 - alpha_seq) * expm1(-(alpha_seq - sigma^2) * delta) +
       2 * mu_seq * expm1(-alpha_seq * delta) ) * fh_half * exp(-alpha_seq * delta) +
    (sigma^2 * (alpha_seq - sigma^2) - 2 * alpha_seq^2 * mu_seq^2) / ( (2 * alpha_seq - sigma^2) * (alpha_seq - sigma^2) ) *
    expm1(-(2 * alpha_seq - sigma^2) * delta) -
    (2 * alpha_seq / (sigma^2 - alpha_seq) * expm1(-alpha_seq * delta) -
       1 + exp(-alpha_seq * delta) * (1 - expm1(-alpha_seq * delta))) * mu_seq^2
  
  if(any(phi_dot < 0) | any(is.na(phi_dot))){return(100000 + runif(1, min = -5, max = 1))}
  
  sd.part  <- sigma * sqrt(phi_dot)
  
  mu.part  <- exp(-alpha_seq * delta) * (fh_half - mu_seq) + mu_seq
  
  -sum(stats::dnorm(fh_half_inv, mu.part, sd.part, log = TRUE)) - sum(log(abs(det_Dfh_half_inv)))
}


# New idea
OU_dynamic_likelihood_new <- function(par, data, delta,
                                             alpha0, mu0, sigma){
  tau     <-  par[1]
  A       <-  par[2]
  
  nu      <- if(length(par) == 3) exp(par[3]) else 1
  
  N       <- length(data)
  Xupp    <- data[2:N]
  Xlow    <- data[1:(N-1)]
  time    <- delta * (1:(N-1))
  
  m          <- mu0 - alpha0 / (2 * A)
  lambda0    <- -alpha0^2 / (4 * A)
  lam_seq    <- lambda0 * (1 - time / tau)^nu
  
  sd.part <- (1 - A * (Xlow - m) * delta) * sigma * sqrt(delta)
  
  mu.part <- Xlow - (A * (Xlow - m)^2 + lam_seq) * delta * (1 - A * (Xlow - m) * delta) - 
    A / 2 * sigma^2 * (delta)^2
  
  -sum(stats::dnorm(Xupp, mu.part, sd.part, log = TRUE), na.rm = TRUE)
}

OU_dynamic_likelihood_new_numeric_grad <- function(par, data, delta,
                                                   alpha0, mu0, sigma){#, grMethod = "Richardson"){
  numDeriv::grad(func = OU_dynamic_likelihood_new, x = par,
                 method = "Richardson",  
                 delta = delta,
                 data = data,
                 alpha0 = alpha0,
                 mu0 = mu0,
                 sigma = sigma)
}

# Additive noise model
true_param <- c(0.87, -1.51, -2.69, 0.1, 1.5)
actual_dt <- 0.1
tau <- 150
t_0 <- 50
sim_res_add <- simulate_additive_noise_tipping_model(actual_dt, true_param, tau, t_0)
sample_n(sim_res_add, min(nrow(sim_res_add), 10000)) |> ggplot(aes(x = t, y = X_t)) + geom_step() +
  geom_line(aes(y = mu_t))
# Stationary part
# Parameters for stationary part
mu0 <- true_param[2] + ifelse(true_param[1] >= 0, 1, -1) *  sqrt(abs(true_param[3] / true_param[1]))
alpha0 <- 2 * sqrt(abs(true_param[1] * true_param[3]))
stationary_part_true_param <- c(alpha0, mu0, true_param[4])

stationary_param_estim <- optimize_stationary_likelihood(likelihood_fun = OU_likelihood,
                               data = sim_res_add$X_t[sim_res_add$t < t_0],
                               init_par = stationary_part_true_param,
                               delta = actual_dt, exp_sigma = FALSE)$par

dynamic_part_true_param <- c(tau, true_param[1], true_param[5])


optimize_dynamic_likelihood(likelihood_fun = OU_dynamic_likelihood,
                            data = sim_res_add$X_t[sim_res_add$t > t_0],
                            init_par = c(150, 1, 1),
                            delta = actual_dt,
                            alpha0 = stationary_param_estim[1],
                            mu0 = stationary_param_estim[2],
                            sigma = stationary_param_estim[3],
                            control = list(reltol = sqrt(.Machine$double.eps) / 1000)
                            )

optimize_dynamic_likelihood(likelihood_fun = OU_dynamic_likelihood_new,
                            data = sim_res_add$X_t[sim_res_add$t > t_0],
                            init_par = c(150, 1, 0),
                            delta = actual_dt,
                            alpha0 = stationary_param_estim[1],
                            mu0 = stationary_param_estim[2],
                            sigma = stationary_param_estim[3],
                            control = list(reltol = sqrt(.Machine$double.eps) / 1000))

nleqslv::nleqslv(x = c(160, 1, 1), fn = OU_dynamic_likelihood_new_numeric_grad,
                 data = sim_res_add$X_t[sim_res_add$t > t_0],
                 delta = actual_dt,
                 alpha0 = stationary_param_estim[1],
                 mu0 = stationary_param_estim[2],
                 sigma = stationary_param_estim[3])$x

sum_gradients <- function(par){
  sum(OU_dynamic_likelihood_new_numeric_grad(par,
                                             data = sim_res_add$X_t[sim_res_add$t > t_0],
                                             delta = actual_dt,
                                             alpha0 = stationary_param_estim[1],
                                             mu0 = stationary_param_estim[2],
                                             sigma = stationary_param_estim[3]))^2
}

optimx::optimx(par = c(150, 1, 1), fn = sum_gradients, method = "L-BFGS-B")

sim_res_add$t[sim_res_add$t > t_0] - t_0

stats::nls(X_t ~ stationary_param_estim[2] -
             stationary_param_estim[1] / (2 * A) * (1 + sqrt((1 - (t - t_0) / tau)^nu)),
           data = filter(sim_res_add, t > t_0), start = list(A = dynamic_part_true_param[1],
                                            tau = dynamic_part_true_param[2],
                                            nu = dynamic_part_true_param[3]))

tibble(t = sim_res_add$t[sim_res_add$t > t_0],
       mu_estim_real = true_param[2] -
         sqrt(-true_param[3] *(1 - (sim_res_add$t[sim_res_add$t > t_0] - t_0) /
                   dynamic_part_true_param[1])^dynamic_part_true_param[3] / true_param[1]),
       
       mu_estim = stationary_part_true_param[2] +
         stationary_part_true_param[1] / (2 * dynamic_part_true_param[2]) * (1 +
  (sqrt((1 - (sim_res_add$t[sim_res_add$t > t_0] - t_0) /
             dynamic_part_true_param[1])^dynamic_part_true_param[3]))),
  mu_real = sim_res_add$mu_t[sim_res_add$t > t_0])

