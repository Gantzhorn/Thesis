numeric_strang_splitting <- function(par, data, delta, drift_lamperti_sde) {
  x0 <- data[1:(length(data) - 1)]
  x1 <- data[2:length(data)]
  
  beta  <- par[1]
  mu    <- par[2]
  sigma <- par[3]

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
  -sum(stats::dnorm(inv_f, mean = mu_f, sd = sd_f, log = TRUE)) - sum(log(abs(df)))
}


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
                     
optimize_stationary_likelihood(likelihood_fun = numeric_strang_splitting,
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

optimize_stationary_likelihood(numeric_strang_splitting,
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
  likelihood_fun = numeric_strang_splitting,
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
  likelihood_fun = numeric_strang_splitting,
  data = 2 * asinh(sqrt(F_sim_dynamic$X_t[F_sim_dynamic$t<t_0])),
  init_par = stationary_part_true_param,
  delta = actual_dt,
  exp_sigma = FALSE,
  drift_lamperti_sde = F_lamperti_drift)

numeric_strang_splitting(par = stationary_part_true_param,
                         data = 2 * asinh(sqrt(F_sim_dynamic$X_t[F_sim_dynamic$t<t_0])),
                         delta = actual_dt,
                         drift_lamperti_sde = F_lamperti_drift)


dummy_func <- function(y, par){
  beta  <- par[1]
mu    <- par[2]
sigma <- par[3]
(beta + sigma^2 / 2) * cosh(y) - beta * (2 * mu + 1)
}


(beta + sigma^2 / 2) * cosh(nleqslv::nleqslv(x = mu, fn = dummy_func, par = par)$x) - beta * (2 * mu + 1)


nleqslv::nleqslv(x = mu, fn = dummy_func, par = par)$x
numDeriv::grad(func = F_lamperti_drift, x = b, par = par)

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
  likelihood_fun = numeric_strang_splitting,
  data = 2 * asin(sqrt(sim_res_jacobi$X_t[sim_res_jacobi$t<t_0])),
  init_par = stationary_part_true_param,
  delta = actual_dt,
  exp_sigma = FALSE,
  drift_lamperti_sde = jacobi_lamperti_drift)

pure_numeric_strang_splitting <- function(par,
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
  likelihood_fun = pure_numeric_strang_splitting,
  data = sim_res_jacobi$X_t[sim_res_jacobi$t<t_0],
  init_par = stationary_part_true_param,
  delta = actual_dt,
  exp_sigma = FALSE,
  diffusion_term = jacobi_diffusion_term,
  xi = 0.01)


library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
# Define the bounding box for the map (North Atlantic region)
bbox <- c(xmin = -70, xmax = 5, ymin = 45, ymax = 70)

# Get world data and transform to simple features (sf)
world <- ne_countries(scale = "medium", returnclass = "sf")

# Create a polygon for the subpolar gyre (approximation)
# subpolar_gyre <- st_polygon(list(rbind(
#   c(-45, 60),  # Near southern tip of Greenland
#   c(-30, 64),  # Southern tip of Iceland
#   c(-10, 60),  # Waters northwest of the UK
#   c(-20, 55),  # North of Newfoundland
#   c(-40, 55),  # South of Greenland
#   c(-45, 60)   # Closing the loop near southern tip of Greenland
# ))) %>% 
#   st_sfc(crs = st_crs(world))

# # Convert to sf object
# subpolar_gyre_sf <- st_sf(geometry = subpolar_gyre)

# Plot the map
NorthAtlanticOcean <- ggplot(data = world) +
  geom_sf(fill = "white", color = "black") +  # landmasses
  #geom_sf(data = subpolar_gyre_sf, fill = "blue", alpha = 0.5) +  # subpolar gyre
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE) +
  labs(title = "", x = "Longitude", y = "Latitude") + 
  theme(
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 16)
  )

ggsave(NorthAtlanticOcean, path = paste0(getwd(), "/tex_files/figures"),
       filename = "NorthAtlanticOcean.jpeg",
       height = 8, width = 14, dpi = 300, units = "in", device = "jpeg",
       limitsize = FALSE, scale = 1)


