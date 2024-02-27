# Title: Running Statistics
# Author: Anders Gantzhorn Kristensen (University of Copenhagen, andersgantzhorn@gmail.com)
# Date: 2024-02-22 (Last Updated: 2024-02-22)
#-----------------------------------------------------------------------------------------------------------------------------#
# Project: Tipping Point Estimation in Ecological Systems using Stochastic Differential Equations
# Description: This script implement various methods for analyzing sample paths of stochastic differential equations
#-----------------------------------------------------------------------------------------------------------------------------#
# License: MIT License (for more information, see LICENSE file in the repository).
# Dependencies: zoo
# Usage: See example usage in the bottom of the script
# Acknowledgements: "zoo: S3 Infrastructure for Regular and Irregular Time Series" - Zeileis A, Grothendieck G (2005).
# Keywords: Time series, Running statistics
#-----------------------------------------------------------------------------------------------------------------------------#

running_variance <- function(data, width = NA){
  # Make default size based on dataset
  if(is.na(width)){
    width <- round(length(data)/10)
  }
  zoo::rollapply(data, width = width, FUN = var, fill = NA)
}

running_mean <- function(data, k = NA){
  # Make default size based on dataset
  if(is.na(k)){
    k <- round(length(data)/10)
  }
  zoo::rollmean(data, k = k, fill = NA)
}

running_autocorrelation <- function(data, width = NA){
  # Make default size based on dataset
  if(is.na(width)){
    width <- round(length(data)/10)
  }
  autocorr_lag1_fun <- function(x){acf(x,  lag.max = 1, plot = FALSE)$acf[2]}
  zoo::rollapply(data, width = width, FUN = autocorr_lag1_fun, fill = NA)
}

#-----------------------------------------------------------------------------------------------------------------------------#

## Example usage
## Take methods from tipping_simulations.R
# source("tipping_simulations.R")
# true_param <- c(0.87, -1.51, -2.69, sqrt(0.3))
# actual_dt <- 1/12
# tau <- 100
# t_0 <- 50
# sim_res_add <- simulate_additive_noise_tipping_model(actual_dt, true_param, tau, t_0)
# sim_res_add |> ggplot2::ggplot(ggplot2::aes(x = t, y = X_weak_2.0)) +
#   ggplot2::geom_step() + ggplot2::geom_hline(yintercept = true_param[2], linetype = "dashed") +
#   ggplot2::geom_vline(xintercept = t_0)
# run_var <- running_variance(sim_res_add$X_weak_2.0)
# 
# tibble::tibble(running_variance = run_var, t = sim_res_add$t) |> 
#   ggplot2::ggplot(ggplot2::aes(x = t, y = running_variance)) + ggplot2::geom_line() + 
#   ggplot2::geom_smooth(se = F, method = "lm")
# 
# 
# run_acf <- running_autocorrelation(sim_res_add$X_weak_2.0)
# tibble::tibble(running_autocorrelation = run_acf, t = sim_res_add$t) |> 
#   ggplot2::ggplot(ggplot2::aes(x = t, y = running_autocorrelation)) + ggplot2::geom_line() + 
#   ggplot2::geom_smooth(se = F, method = "lm")
