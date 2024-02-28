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

resid_strang_CIR <- function(par, data, delta) {
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
